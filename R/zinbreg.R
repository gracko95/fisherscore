#' Control parameters for zinbreg()
#'
#' @param tol EM convergence tolerance: stop when the observed-data
#'   log-likelihood changes by less than `tol` between two EM iterations.
#' @param maxit Maximum number of EM iterations.
#' @param inner_tol Tolerance for the IRLS loops inside each M-step.
#' @param inner_maxit Maximum number of iterations of each inner loop.
#' @param theta_max Upper limit for \eqn{\theta} (see [nbreg_control()]).
#' @param hessian How to compute the observed information for the standard
#'   errors: `"analytic"` (exact second derivatives, the default) or
#'   `"numeric"` (`numDeriv::hessian()` of the package's own
#'   log-likelihood; needs the numDeriv package).
#' @param trace Logical; print the observed-data log-likelihood at each EM
#'   iteration?
#' @return A list of control parameters.
#' @examples
#' zinbreg_control(tol = 1e-10, maxit = 2000)
#' @export
zinbreg_control <- function(tol = 1e-8, maxit = 1000, inner_tol = 1e-10, inner_maxit = 100,
                            theta_max = 1e5, hessian = c("analytic", "numeric"),
                            trace = FALSE) {
  hessian <- match.arg(hessian)
  nb <- nbreg_control(tol = tol, maxit = maxit, inner_tol = inner_tol,
                      inner_maxit = inner_maxit, theta_max = theta_max, trace = trace)
  c(nb, list(hessian = hessian))
}

#' Zero-inflated negative binomial regression by EM with IRLS M-steps
#'
#' Fits the ZINB model
#' \deqn{P(Y_i = 0) = \pi_i + (1-\pi_i) f_{NB}(0;\mu_i,\theta), \qquad
#'       P(Y_i = y) = (1-\pi_i) f_{NB}(y;\mu_i,\theta),\ y \ge 1,}
#' with \eqn{\log\mu_i = x_i^\top\beta + o_i} (count model) and
#' \eqn{\mathrm{logit}\,\pi_i = z_i^\top\gamma} (zero-inflation model), by an
#' EM algorithm:
#' * **E-step**: \eqn{\tau_i = \pi_i / \{\pi_i + (1-\pi_i) f_{NB}(0;\mu_i,\theta)\}}
#'   if \eqn{y_i = 0}, else \eqn{\tau_i = 0}.
#' * **M-step, zero part**: logistic IRLS on the fractional response
#'   \eqn{\tau_i}.
#' * **M-step, count part**: the NB fit of [nbreg()] with prior weights
#'   \eqn{a_i = 1 - \tau_i}.
#'
#' EM stops when the observed-data log-likelihood changes by less than
#' `control$tol`; a warning is given if it ever decreases. Standard errors
#' come from inverting the observed information of the observed-data
#' log-likelihood in \eqn{(\beta, \gamma, \log\theta)}.
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model
#'   covariates before `|`, zero-model covariates after it. Without `|` the
#'   same covariates are used in both parts, as in `pscl::zeroinfl()`.
#'   Offsets (`offset()` terms) are allowed in the count part.
#' @param data A data frame containing the variables in `formula`.
#' @param offset Optional offset for the count model, evaluated in `data`.
#' @param control A list from [zinbreg_control()].
#' @return An object of class `"zinbreg"`, a list with components
#'   `coefficients` (a list with `count` and `zero`), `vcov` (for
#'   \eqn{(\beta,\gamma)}, named `count_*` and `zero_*` as in `zeroinfl()`),
#'   `vcov_full` (also including `log(theta)`), `theta`, `SE.theta`,
#'   `SE.logtheta`, `loglik`, `iter`, `converged`, `monotone` (did the
#'   log-likelihood never decrease?), `fitted.values` (\eqn{(1-\hat\pi)\hat\mu}),
#'   `mu`, `pi`, `tau`, `residuals` (response), `score` (gradient at the
#'   estimate), `trace`, `y`, `offset`, `nobs`, `df.residual`, `call`,
#'   `formula`, `terms`, `xlevels`, `contrasts` and `control`.
#' @seealso [nbreg()], [zinbreg-methods]
#' @examples
#' # Simulated data: count part beta = (1, 0.5), zero part gamma = (-1, 1), theta = 2
#' set.seed(1)
#' n <- 500
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' structural_zero <- rbinom(n, 1, plogis(-1 + d$z))
#' d$y <- ifelse(structural_zero == 1, 0, rnbinom(n, size = 2, mu = exp(1 + 0.5 * d$x)))
#' fit <- zinbreg(y ~ x | z, data = d)
#' summary(fit)
#' head(predict(fit, type = "zero"))
#' @export
zinbreg <- function(formula, data, offset, control = zinbreg_control()) {
  call <- match.call()
  if (missing(formula) || !inherits(formula, "formula"))
    stop("'formula' must be a model formula such as y ~ x1 | z1", call. = FALSE)
  if (missing(data)) data <- environment(formula)
  else if (!is.data.frame(data)) stop("'data' must be a data frame", call. = FALSE)
  if (!is.list(control)) stop("'control' must be a list from zinbreg_control()", call. = FALSE)
  control <- do.call(zinbreg_control, control)

  # Two-part formula. Without "|", use the same right-hand side for both parts.
  F <- Formula::as.Formula(formula)
  if (length(F)[2L] > 2L) stop("the formula can have at most two parts: y ~ x | z", call. = FALSE)
  if (length(F)[2L] == 1L)
    F <- Formula::as.Formula(stats::formula(F), stats::formula(F, lhs = 0L, rhs = 1L))

  # Model frame for both parts, then one design matrix per part
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- F
  mf$data <- data
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  mtX <- stats::terms(F, data = data, rhs = 1L)
  mtZ <- stats::delete.response(stats::terms(F, data = data, rhs = 2L))
  y <- stats::model.response(mf)
  if (is.null(y)) stop("the formula has no response", call. = FALSE)
  check_counts(y)
  if (!any(y == 0))
    stop("the response has no zeros, so the zero-inflation part cannot be estimated",
         call. = FALSE)
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)
  check_full_rank(X, "count-model matrix")
  check_full_rank(Z, "zero-model matrix")
  n <- length(y)

  # Count-model offset: offset() terms in the count formula plus 'offset ='
  off <- stats::model.offset(Formula::model.part(F, data = mf, rhs = 1L, terms = TRUE))
  if (is.null(off)) off <- rep(0, n)
  if (!is.null(mf[["(offset)"]])) {
    if (length(mf[["(offset)"]]) != n) stop("'offset' must have one value per observation", call. = FALSE)
    off <- off + mf[["(offset)"]]
  }

  fit <- zinb_fit(X, Z, y, off, control)

  structure(c(fit, list(
    y = y, offset = off, nobs = n, df.residual = n - ncol(X) - ncol(Z) - 1,
    call = call, formula = F,
    terms = list(count = mtX, zero = mtZ, full = mt),
    xlevels = list(count = stats::.getXlevels(mtX, mf), zero = stats::.getXlevels(mtZ, mf)),
    contrasts = list(count = attr(X, "contrasts"), zero = attr(Z, "contrasts")),
    control = control, na.action = attr(mf, "na.action"))), class = "zinbreg")
}

# ---- EM engine ----------------------------------------------------------------

zinb_fit <- function(X, Z, y, offset, control) {
  n <- length(y)
  inner <- control
  inner$trace <- FALSE
  inner$tol <- min(control$tol, 1e-10)      # M-steps solved more tightly than EM

  # Starting values: NB fit for (beta, theta); logistic IRLS of I(y = 0) on Z for gamma
  nb0 <- nb_fit(X, y, offset, a = rep(1, n), control = inner, warn = FALSE)
  beta <- nb0$coefficients
  theta <- nb0$theta
  gamma <- logit_irls(Z, as.numeric(y == 0), gamma = NULL, control = inner)

  mu <- drop(exp(X %*% beta + offset))
  eta_z <- drop(Z %*% gamma)
  ll <- zinb_loglik(y, mu, theta, eta_z)
  trace_ll <- numeric(0)
  monotone <- TRUE
  converged <- theta_at_max <- FALSE
  for (iter in seq_len(control$maxit)) {
    # E-step: tau_i = pi_i / (pi_i + (1 - pi_i) f_NB(0; mu_i, theta)) if y_i = 0, else 0
    tau <- zinb_estep(y, mu, theta, eta_z)

    # M-step, zero part: logistic IRLS with fractional response tau
    gamma <- logit_irls(Z, tau, gamma, inner)

    # M-step, count part: NB fit (Section 3.1) with prior weights a_i = 1 - tau_i
    nbm <- nb_fit(X, y, offset, a = 1 - tau, control = inner,
                  beta = beta, theta = theta, warn = FALSE)
    beta <- nbm$coefficients
    theta <- nbm$theta
    theta_at_max <- nbm$theta_at_max

    # Observed-data log-likelihood: EM must never decrease it
    mu <- drop(exp(X %*% beta + offset))
    eta_z <- drop(Z %*% gamma)
    ll_new <- zinb_loglik(y, mu, theta, eta_z)
    trace_ll <- c(trace_ll, ll_new)
    if (ll_new < ll - 1e-8 * (1 + abs(ll))) monotone <- FALSE
    if (control$trace)
      cat(sprintf("zinbreg EM iter %4d: logLik = %.10f  change = %+.3e  theta = %.6f\n",
                  iter, ll_new, ll_new - ll, theta))

    if (abs(ll_new - ll) < control$tol) {
      converged <- TRUE
      ll <- ll_new
      break
    }
    ll <- ll_new
  }
  if (!converged)
    warning(sprintf("zinbreg: EM did not converge in maxit = %d iterations", control$maxit),
            call. = FALSE)
  if (!monotone)
    warning("zinbreg: the observed-data log-likelihood decreased during EM", call. = FALSE)
  if (theta_at_max)
    warning(sprintf(paste0("zinbreg: theta reached theta_max = %g (Poisson limit); ",
                           "consider a zero-inflated Poisson model"), control$theta_max),
            call. = FALSE)

  # Standard errors: invert the observed information in (beta, gamma, log theta)
  names(beta) <- colnames(X)
  names(gamma) <- colnames(Z)
  par_names <- c(paste0("count_", colnames(X)), paste0("zero_", colnames(Z)), "log(theta)")
  tau <- zinb_estep(y, mu, theta, eta_z)
  if (control$hessian == "analytic") {
    H <- zinb_hessian(X, Z, y, mu, theta, eta_z, tau)
  } else {
    if (!requireNamespace("numDeriv", quietly = TRUE))
      stop("hessian = \"numeric\" needs the numDeriv package", call. = FALSE)
    p <- ncol(X)
    q <- ncol(Z)
    ll_fun <- function(par) {
      zinb_loglik(y, drop(exp(X %*% par[1:p] + offset)), exp(par[p + q + 1]),
                  drop(Z %*% par[(p + 1):(p + q)]))
    }
    H <- numDeriv::hessian(ll_fun, c(beta, gamma, log(theta)))
  }
  dimnames(H) <- list(par_names, par_names)
  vc_full <- tryCatch(solve(-H), error = function(e) {
    warning("zinbreg: the observed information is singular; standard errors are NA",
            call. = FALSE)
    matrix(NA_real_, nrow(H), ncol(H), dimnames = dimnames(H))
  })
  k <- length(par_names)
  se_logtheta <- sqrt(vc_full[k, k])

  pi <- stats::plogis(eta_z)
  list(coefficients = list(count = beta, zero = gamma),
       vcov = vc_full[-k, -k, drop = FALSE], vcov_full = vc_full,
       theta = theta, SE.logtheta = se_logtheta, SE.theta = theta * se_logtheta,
       loglik = ll, iter = iter, converged = converged, monotone = monotone,
       theta_at_max = theta_at_max,
       fitted.values = (1 - pi) * mu, mu = mu, pi = pi, tau = tau,
       residuals = y - (1 - pi) * mu,
       score = zinb_score(X, Z, y, mu, theta, eta_z, tau),
       trace = data.frame(iter = seq_along(trace_ll), loglik = trace_ll))
}

# log(exp(a) + exp(b)) without overflow/underflow
log_add_exp <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

# Observed-data log-likelihood, Section 2.2:
#   y_i = 0: log[ pi_i + (1 - pi_i) f_NB(0; mu_i, theta) ]
#   y_i > 0: log(1 - pi_i) + log f_NB(y_i; mu_i, theta)
# computed on the log scale from eta_z = logit(pi) for numerical stability.
zinb_loglik <- function(y, mu, theta, eta_z) {
  log_pi <- stats::plogis(eta_z, log.p = TRUE)       # log pi
  log_1mpi <- stats::plogis(-eta_z, log.p = TRUE)    # log(1 - pi)
  log_f <- stats::dnbinom(y, size = theta, mu = mu, log = TRUE)
  sum(ifelse(y == 0, log_add_exp(log_pi, log_1mpi + log_f), log_1mpi + log_f))
}

# E-step, Section 3.2: for y_i = 0
#   tau_i = pi_i / (pi_i + (1 - pi_i) f0_i) = 1 / (1 + exp(-eta_z_i + log f0_i))
#         = plogis(eta_z_i - log f0_i),  with log f0_i = theta log(theta / (theta + mu_i)),
# and tau_i = 0 for y_i > 0.
zinb_estep <- function(y, mu, theta, eta_z) {
  log_f0 <- stats::dnbinom(0, size = theta, mu = mu, log = TRUE)
  ifelse(y == 0, stats::plogis(eta_z - log_f0), 0)
}

# M-step for the zero part, Section 3.2: maximize
#   sum_i [ tau_i log pi_i + (1 - tau_i) log(1 - pi_i) ]   over gamma
# by IRLS: w_i = pi_i (1 - pi_i),  z_i = z_i'gamma + (tau_i - pi_i) / {pi_i (1 - pi_i)},
#          gamma <- (Z'WZ)^{-1} Z'Wz.
# With gamma = NULL the first step starts from pi = (t + 0.5) / 2 (as glm() does).
logit_irls <- function(Z, t, gamma, control) {
  q_fun <- function(eta) sum(t * stats::plogis(eta, log.p = TRUE) +
                               (1 - t) * stats::plogis(-eta, log.p = TRUE))
  if (is.null(gamma)) {
    eta <- stats::qlogis((t + 0.5) / 2)
    q <- -Inf
  } else {
    eta <- drop(Z %*% gamma)
    q <- q_fun(eta)
  }
  for (it in seq_len(control$inner_maxit)) {
    p <- stats::plogis(eta)
    w <- p * (1 - p)
    z <- eta + (t - p) / w
    # If pi_i reaches 0 or 1 the weights vanish and Z'WZ becomes singular: the
    # zero-model covariates (nearly) separate structural zeros from counts.
    gamma_new <- tryCatch(wls(Z, z, w), error = function(e) NULL)
    if (is.null(gamma_new) || any(!is.finite(gamma_new)))
      stop(paste0("the zero-model IRLS failed: the zero-model covariates (nearly) perfectly ",
                  "separate zeros from positive counts, so the zero-inflation coefficients ",
                  "have no finite estimate. Simplify the zero model (after |)."), call. = FALSE)

    eta_new <- drop(Z %*% gamma_new)
    q_new <- q_fun(eta_new)
    halve <- 0                     # step-halving: never decrease the objective
    while ((!is.finite(q_new) || q_new < q) && !is.null(gamma) && halve < 30) {
      gamma_new <- (gamma_new + gamma) / 2
      eta_new <- drop(Z %*% gamma_new)
      q_new <- q_fun(eta_new)
      halve <- halve + 1
    }
    change <- if (is.null(gamma)) Inf else max(abs(gamma_new - gamma))
    gamma <- gamma_new
    eta <- eta_new
    q <- q_new
    if (change < control$inner_tol) break
  }
  gamma
}

# ---- Analytic derivatives of the observed-data log-likelihood -----------------
# Notation for observation i: eta_c = x'beta + o (count), eta_z = z'gamma (zero),
# phi = log(theta), h = log f_NB(y; mu, theta), rho = 1 - tau.
# Derivatives of h (mu = exp(eta_c)):
#   h_c   = theta (y - mu) / (theta + mu)
#   h_cc  = -theta mu (theta + y) / (theta + mu)^2
#   h_p   = theta U_i                         (U_i: the theta score term of Section 3.1(b))
#   h_pp  = theta^2 H_i + theta U_i           (H_i: the theta Hessian term of Section 3.1(b))
#   h_cp  = theta mu (y - mu) / (theta + mu)^2
# For y = 0, l = log(pi + (1 - pi) e^h); for y > 0, l = log(1 - pi) + h. Both give
#   dl/deta_z = tau - pi,   dl/da = rho h_a        (a = eta_c or phi)
#   d2l/deta_z^2      = tau (1 - tau) - pi (1 - pi)
#   d2l/deta_z da     = -tau rho h_a
#   d2l/da db         = rho h_ab + tau rho h_a h_b
# (tau = 0, rho = 1 when y > 0, so these reduce to the NB derivatives).
zinb_derivs <- function(y, mu, theta, eta_z, tau) {
  pi <- stats::plogis(eta_z)
  rho <- 1 - tau
  U_i <- digamma(y + theta) - digamma(theta) + log(theta) + 1 -
    log(theta + mu) - (y + theta) / (theta + mu)
  H_i <- trigamma(y + theta) - trigamma(theta) + 1 / theta -
    2 / (theta + mu) + (y + theta) / (theta + mu)^2
  h_c <- theta * (y - mu) / (theta + mu)
  h_cc <- -theta * mu * (theta + y) / (theta + mu)^2
  h_p <- theta * U_i
  h_pp <- theta^2 * H_i + theta * U_i
  h_cp <- theta * mu * (y - mu) / (theta + mu)^2
  list(s_c = rho * h_c, s_z = tau - pi, s_p = rho * h_p,
       l_cc = rho * h_cc + tau * rho * h_c^2,
       l_cp = rho * h_cp + tau * rho * h_c * h_p,
       l_pp = rho * h_pp + tau * rho * h_p^2,
       l_zz = tau * (1 - tau) - pi * (1 - pi),
       l_zc = -tau * rho * h_c,
       l_zp = -tau * rho * h_p)
}

# Hessian in (beta, gamma, phi) by the chain rule: d eta_c / d beta = x_i, d eta_z / d gamma = z_i
zinb_hessian <- function(X, Z, y, mu, theta, eta_z, tau) {
  d <- zinb_derivs(y, mu, theta, eta_z, tau)
  H_bb <- crossprod(X, X * d$l_cc)
  H_bg <- crossprod(X, Z * d$l_zc)
  H_gg <- crossprod(Z, Z * d$l_zz)
  H_bp <- crossprod(X, d$l_cp)
  H_gp <- crossprod(Z, d$l_zp)
  H_pp <- sum(d$l_pp)
  rbind(cbind(H_bb, H_bg, H_bp),
        cbind(t(H_bg), H_gg, H_gp),
        cbind(t(H_bp), t(H_gp), H_pp))
}

# Score (gradient) in (beta, gamma, phi); it is ~0 at the MLE
zinb_score <- function(X, Z, y, mu, theta, eta_z, tau) {
  d <- zinb_derivs(y, mu, theta, eta_z, tau)
  c(drop(crossprod(X, d$s_c)), drop(crossprod(Z, d$s_z)), sum(d$s_p))
}
