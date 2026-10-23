#' Control parameters for nbreg()
#'
#' @param tol Outer convergence tolerance: stop when the log-likelihood
#'   changes by less than `tol` between two outer iterations.
#' @param maxit Maximum number of outer iterations (one IRLS fit for
#'   \eqn{\beta} followed by one Newton-Raphson fit for \eqn{\theta}).
#' @param inner_tol Tolerance for the inner IRLS and Newton-Raphson loops
#'   (largest absolute change in \eqn{\beta}, or in \eqn{\log\theta}).
#' @param inner_maxit Maximum number of iterations of each inner loop.
#' @param theta_max Upper limit for \eqn{\theta}. If the estimate reaches it,
#'   the data show no overdispersion (the Poisson limit) and a warning is given.
#' @param trace Logical; print the log-likelihood and \eqn{\theta} at each
#'   outer iteration?
#' @return A list of control parameters.
#' @examples
#' nbreg_control(tol = 1e-10, trace = TRUE)
#' @export
nbreg_control <- function(tol = 1e-8, maxit = 100, inner_tol = 1e-10, inner_maxit = 100,
                          theta_max = 1e5, trace = FALSE) {
  check_positive(tol, "tol")
  check_positive(inner_tol, "inner_tol")
  check_positive(maxit, "maxit", min = 0.5)        # i.e. maxit >= 1
  check_positive(inner_maxit, "inner_maxit", min = 0.5)
  check_positive(theta_max, "theta_max")
  list(tol = tol, maxit = as.integer(maxit), inner_tol = inner_tol,
       inner_maxit = as.integer(inner_maxit), theta_max = theta_max, trace = isTRUE(trace))
}

#' Negative binomial regression by IRLS and Newton-Raphson
#'
#' Fits the NB2 regression model
#' \deqn{\log\mu_i = x_i^\top\beta + o_i, \qquad \mathrm{Var}(Y_i) = \mu_i + \mu_i^2/\theta}
#' by maximum likelihood, alternating
#' * IRLS for \eqn{\beta} with \eqn{\theta} fixed: working weights
#'   \eqn{w_i = \theta\mu_i/(\theta+\mu_i)}, working response
#'   \eqn{z_i = (\eta_i - o_i) + (y_i-\mu_i)/\mu_i},
#'   \eqn{\beta \leftarrow (X^\top WX)^{-1}X^\top Wz}; and
#' * Newton-Raphson for \eqn{\phi = \log\theta} with \eqn{\beta} fixed,
#'
#' until the log-likelihood changes by less than `control$tol`. No general
#' purpose optimizer is used. The parameterization is the one used by
#' `MASS::glm.nb()`.
#'
#' Standard errors for \eqn{\beta} are from \eqn{(X^\top\hat WX)^{-1}},
#' treating \eqn{\theta} as known (as `glm.nb()` does), and
#' \eqn{\mathrm{SE}(\hat\theta) = \{-H(\hat\theta)\}^{-1/2}}.
#'
#' @param formula A model formula, e.g. `y ~ x1 + x2`. Offsets may be given
#'   with `offset()` terms.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional offset for the linear predictor (e.g. log exposure),
#'   evaluated in `data`.
#' @param control A list from [nbreg_control()].
#' @return An object of class `"nbreg"`, a list with components
#'   `coefficients`, `vcov`, `theta`, `SE.theta`, `loglik`, `iter`,
#'   `converged`, `fitted.values` (\eqn{\hat\mu}), `linear.predictors`,
#'   `residuals` (response), `y`, `offset`, `df.residual`, `nobs`,
#'   `trace` (log-likelihood and \eqn{\theta} by iteration), `call`,
#'   `terms`, `xlevels`, `contrasts` and `control`.
#' @seealso [zinbreg()], [nbreg-methods]
#' @examples
#' # Simulated data with known beta = (1, 0.5) and theta = 2
#' set.seed(1)
#' d <- data.frame(x = rnorm(300))
#' d$y <- rnbinom(300, size = 2, mu = exp(1 + 0.5 * d$x))
#' fit <- nbreg(y ~ x, data = d)
#' summary(fit)
#'
#' # With an exposure offset
#' d$t <- runif(300, 1, 3)
#' d$y2 <- rnbinom(300, size = 2, mu = d$t * exp(1 + 0.5 * d$x))
#' coef(nbreg(y2 ~ x, data = d, offset = log(t)))
#' @export
nbreg <- function(formula, data, offset, control = nbreg_control()) {
  call <- match.call()
  if (missing(formula) || !inherits(formula, "formula"))
    stop("'formula' must be a model formula such as y ~ x", call. = FALSE)
  if (missing(data)) data <- environment(formula)
  else if (!is.data.frame(data)) stop("'data' must be a data frame", call. = FALSE)
  if (!is.list(control)) stop("'control' must be a list from nbreg_control()", call. = FALSE)
  control <- do.call(nbreg_control, control)

  # Model frame, response, design matrix and offset (as glm() builds them)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$data <- data
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y <- stats::model.response(mf)
  if (is.null(y)) stop("the formula has no response", call. = FALSE)
  check_counts(y)
  X <- stats::model.matrix(mt, mf)
  check_full_rank(X)
  n <- length(y)
  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, n)
  if (length(off) != n) stop("'offset' must have one value per observation", call. = FALSE)

  fit <- nb_fit(X, y, off, a = rep(1, n), control = control)

  structure(c(fit, list(
    y = y, offset = off, nobs = n, df.residual = n - ncol(X) - 1,
    call = call, terms = mt, xlevels = stats::.getXlevels(mt, mf),
    contrasts = attr(X, "contrasts"), control = control,
    na.action = attr(mf, "na.action"))), class = "nbreg")
}

# ---- Fitting engine ---------------------------------------------------------
# nb_fit() works on a model matrix and accepts prior weights a_i, so that
# zinbreg() can reuse it for the count-part M-step with a_i = 1 - tau_i.

nb_fit <- function(X, y, offset, a, control, beta = NULL, theta = NULL, warn = TRUE) {
  # (c) Starting values. beta: Poisson IRLS (theta = Inf gives w_i = a_i mu_i).
  #     theta: moment estimate from E[(y/mu - 1)^2] = 1/mu + 1/theta ~ 1/theta.
  if (is.null(beta)) beta <- nb_irls_beta(X, y, offset, a, theta = Inf, beta = NULL, control)$beta
  mu <- drop(exp(X %*% beta + offset))
  if (is.null(theta)) {
    theta <- sum(a) / sum(a * (y / mu - 1)^2)
    if (!is.finite(theta) || theta <= 0) theta <- 1
    theta <- min(theta, control$theta_max)
  }

  ll <- nb_loglik(y, mu, theta, a)
  trace_ll <- trace_theta <- numeric(0)
  converged <- theta_at_max <- FALSE
  for (iter in seq_len(control$maxit)) {
    # (a) IRLS for beta with theta held fixed
    beta <- nb_irls_beta(X, y, offset, a, theta, beta, control)$beta
    mu <- drop(exp(X %*% beta + offset))

    # (b) Newton-Raphson for phi = log(theta) with beta held fixed
    nr <- nb_newton_theta(y, mu, a, theta, control)
    theta <- nr$theta
    theta_at_max <- nr$at_max

    ll_new <- nb_loglik(y, mu, theta, a)
    trace_ll <- c(trace_ll, ll_new)
    trace_theta <- c(trace_theta, theta)
    if (control$trace)
      cat(sprintf("nbreg iter %3d: logLik = %.10f  theta = %.6f\n", iter, ll_new, theta))

    # (c) Stop when the log-likelihood changes by less than tol
    if (abs(ll_new - ll) < control$tol) {
      converged <- TRUE
      ll <- ll_new
      break
    }
    ll <- ll_new
  }
  if (warn && !converged)
    warning(sprintf("nbreg: no convergence after maxit = %d iterations", control$maxit),
            call. = FALSE)
  if (warn && any(a > 0 & mu < 1e-10))
    warning(paste0("nbreg: fitted means numerically 0 occurred; some coefficients may be ",
                   "infinite (e.g. a covariate pattern with only zero counts)"), call. = FALSE)
  if (warn && theta_at_max)
    warning(sprintf(paste0("nbreg: theta reached theta_max = %g (Poisson limit); ",
                           "the data show little or no overdispersion"), control$theta_max),
            call. = FALSE)

  # (d) Standard errors at convergence
  #     Var(beta) = (X'WX)^{-1} with w_i = a_i theta mu_i / (theta + mu_i)
  #     SE(theta) = {-H(theta)}^{-1/2}
  eta <- drop(X %*% beta + offset)
  w <- a * theta * mu / (theta + mu)
  vc <- xwx_inverse(X, w)
  H <- nb_hess_theta(y, mu, theta, a)
  se_theta <- if (H < 0) 1 / sqrt(-H) else NA_real_
  names(beta) <- colnames(X)

  list(coefficients = beta, vcov = vc, theta = theta, SE.theta = se_theta,
       loglik = ll, iter = iter, converged = converged, theta_at_max = theta_at_max,
       fitted.values = mu, linear.predictors = eta, residuals = y - mu,
       weights = w, trace = data.frame(iter = seq_along(trace_ll), loglik = trace_ll,
                                       theta = trace_theta))
}

# Section 3.1(a): IRLS for beta with theta fixed.
#   w_i = a_i theta mu_i / (theta + mu_i)   (= a_i mu_i when theta = Inf, i.e. Poisson)
#   z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i
#   beta <- (X'WX)^{-1} X'Wz
# If beta is NULL the first step starts from mu = y + 0.1 (as glm() does for
# the Poisson family). A step that lowers the log-likelihood is halved.
nb_irls_beta <- function(X, y, offset, a, theta, beta, control) {
  if (is.null(beta)) {
    mu <- y + 0.1
    eta <- log(mu)
    ll <- -Inf
  } else {
    eta <- drop(X %*% beta + offset)
    mu <- exp(eta)
    ll <- nb_loglik(y, mu, theta, a)
  }
  for (it in seq_len(control$inner_maxit)) {
    w <- if (is.infinite(theta)) a * mu else a * theta * mu / (theta + mu)
    z <- (eta - offset) + (y - mu) / mu
    # If mu_i underflows to 0 the weights vanish and X'WX becomes singular:
    # some coefficient is heading to -Inf (a covariate pattern with only zeros).
    beta_new <- tryCatch(wls(X, z, w), error = function(e) NULL)
    if (is.null(beta_new) || any(!is.finite(beta_new)))
      stop(paste0("IRLS failed: fitted means reached 0, so some coefficients have no finite ",
                  "estimate (e.g. a factor level or covariate pattern with only zero counts)"),
           call. = FALSE)

    eta_new <- drop(X %*% beta_new + offset)
    ll_new <- nb_loglik(y, exp(eta_new), theta, a)
    # Step-halving: never accept a step that decreases the log-likelihood
    halve <- 0
    while ((!is.finite(ll_new) || ll_new < ll) && !is.null(beta) && halve < 30) {
      beta_new <- (beta_new + beta) / 2
      eta_new <- drop(X %*% beta_new + offset)
      ll_new <- nb_loglik(y, exp(eta_new), theta, a)
      halve <- halve + 1
    }
    change <- if (is.null(beta)) Inf else max(abs(beta_new - beta))
    beta <- beta_new
    eta <- eta_new
    mu <- exp(eta)
    ll <- ll_new
    if (change < control$inner_tol) break
  }
  list(beta = beta, iter = it)
}

# Section 3.1(b): Newton-Raphson for phi = log(theta) with beta (so mu) fixed.
#   U_phi = theta U(theta),  H_phi = theta^2 H(theta) + theta U(theta)
#   phi <- phi - U_phi / H_phi
# If H_phi >= 0 (log-likelihood not concave here) a unit step in the uphill
# direction is used instead; steps that lower the log-likelihood are halved.
nb_newton_theta <- function(y, mu, a, theta, control) {
  phi <- log(theta)
  ll <- nb_loglik(y, mu, theta, a)
  at_max <- FALSE
  for (it in seq_len(control$inner_maxit)) {
    U <- nb_score_theta(y, mu, theta, a)
    H <- nb_hess_theta(y, mu, theta, a)
    U_phi <- theta * U
    H_phi <- theta^2 * H + theta * U
    step <- if (H_phi < 0) -U_phi / H_phi else sign(U_phi)

    ll_new <- nb_loglik(y, mu, exp(phi + step), a)
    halve <- 0
    while ((!is.finite(ll_new) || ll_new < ll) && halve < 30) {
      step <- step / 2
      ll_new <- nb_loglik(y, mu, exp(phi + step), a)
      halve <- halve + 1
    }
    if (halve == 30) break                  # no uphill step left: at the maximum
    phi <- phi + step
    theta <- exp(phi)
    ll <- ll_new
    if (theta >= control$theta_max) {       # Poisson limit: stop and flag it
      theta <- control$theta_max
      at_max <- TRUE
      break
    }
    if (abs(step) < control$inner_tol) break
  }
  list(theta = theta, at_max = at_max, iter = it)
}
