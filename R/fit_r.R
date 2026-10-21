#' Fisher scoring in pure R
#'
#' Fits a GLM to a model matrix by Fisher scoring, also called iteratively
#' reweighted least squares (IRLS). At iteration \eqn{t}, with
#' \eqn{\eta = X\beta + \mathrm{offset}} and \eqn{\mu = g^{-1}(\eta)}:
#'
#' * working weights \eqn{w_i = a_i (d\mu_i/d\eta_i)^2 / V(\mu_i)}, where
#'   \eqn{a_i} are the prior weights
#' * working response \eqn{z_i = \eta_i - \mathrm{offset}_i + (y_i - \mu_i)\, d\eta_i/d\mu_i}
#' * update \eqn{\beta^{(t+1)} = (X^\top W X)^{-1} X^\top W z}
#'
#' This is the same as \eqn{\beta^{(t+1)} = \beta^{(t)} + \mathcal{I}(\beta^{(t)})^{-1} U(\beta^{(t)})},
#' where \eqn{U = X^\top W (z - \eta)} is the score and
#' \eqn{\mathcal{I} = X^\top W X} is the expected (Fisher) information,
#' both taken with dispersion 1.
#' The weighted least squares step is solved with a QR decomposition of
#' \eqn{W^{1/2}X}, so \eqn{X^\top W X} is never formed or inverted directly.
#'
#' Most users should call [fs_glm()]. This function is the low-level engine
#' and assumes `X` has full column rank and all `weights` are positive.
#'
#' @param X Numeric model matrix (n x p).
#' @param y Response vector. For binomial models, the proportion of successes.
#' @param weights Prior weights (for binomial models, the number of trials).
#' @param offset Offset added to the linear predictor.
#' @param family An `fs_family` object, see [fs_family()].
#' @param control A list from [fs_control()].
#' @param eta_start Starting linear predictor (including the offset).
#' @param start Optional starting coefficients. Used only as the point to
#'   fall back to if the first step needs step-halving.
#' @return A list with the coefficients, fitted values, deviance, iteration
#'   count, unscaled covariance matrix \eqn{(X^\top W X)^{-1}} and the
#'   per-iteration trace.
#' @keywords internal
#' @export
fs_fit_r <- function(X, y, weights, offset, family, control, eta_start, start = NULL) {
  nobs <- NROW(y)
  nvars <- ncol(X)
  eps <- control$epsilon

  eta <- eta_start
  mu <- family$linkinv(eta)
  if (!(family$valideta(eta) && family$validmu(mu)))
    stop("cannot find valid starting values: please specify some", call. = FALSE)
  devold <- sum(family$dev.resids(y, mu, weights))
  coefold <- start
  boundary <- conv <- FALSE
  coef <- rep(NA_real_, nvars)
  trace_dev <- numeric(0)
  trace_coef <- matrix(numeric(0), 0, nvars, dimnames = list(NULL, colnames(X)))
  qr_fit <- NULL

  for (iter in seq_len(control$maxit)) {
    # 1. Working weights and working response at the current estimate
    mu_eta <- family$mu.eta(eta)
    z <- (eta - offset) + (y - mu) / mu_eta
    sw <- sqrt(weights * mu_eta^2 / family$variance(mu))

    # 2. Weighted least squares step: solve (X'WX) beta = X'Wz by QR
    qr_fit <- qr(X * sw, tol = min(1e-7, eps / 1000))
    if (qr_fit$rank < nvars)
      stop(sprintf("X'WX is singular at iteration %d (rank %d < %d)",
                   iter, qr_fit$rank, nvars), call. = FALSE)
    start <- qr.coef(qr_fit, z * sw)

    # 3. New linear predictor, mean and deviance
    eta <- drop(X %*% start) + offset
    mu <- family$linkinv(eta)
    dev <- sum(family$dev.resids(y, mu, weights))
    if (control$trace) cat("Deviance =", format(dev, digits = 10), "Iterations -", iter, "\n")

    # 4. Step-halving if the step gave a non-finite deviance or invalid mu/eta
    if (!is.finite(dev)) {
      if (is.null(coefold))
        stop("no valid set of coefficients has been found: please supply starting values",
             call. = FALSE)
      ii <- 1
      while (!is.finite(dev)) {
        if (ii > control$maxit) stop("inner loop 1; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coefold) / 2
        eta <- drop(X %*% start) + offset
        mu <- family$linkinv(eta)
        dev <- sum(family$dev.resids(y, mu, weights))
      }
      boundary <- TRUE
    }
    if (!(family$valideta(eta) && family$validmu(mu))) {
      if (is.null(coefold))
        stop("no valid set of coefficients has been found: please supply starting values",
             call. = FALSE)
      ii <- 1
      while (!(family$valideta(eta) && family$validmu(mu))) {
        if (ii > control$maxit) stop("inner loop 2; cannot correct step size", call. = FALSE)
        ii <- ii + 1
        start <- (start + coefold) / 2
        eta <- drop(X %*% start) + offset
        mu <- family$linkinv(eta)
      }
      boundary <- TRUE
      dev <- sum(family$dev.resids(y, mu, weights))
    }

    trace_dev <- c(trace_dev, dev)
    trace_coef <- rbind(trace_coef, start)

    # 5. Convergence check on the relative change in deviance (as in glm)
    if (abs(dev - devold) / (abs(dev) + 0.1) < eps) {
      conv <- TRUE
      coef <- start
      break
    }
    devold <- dev
    coef <- coefold <- start
  }

  # Unscaled covariance (X'WX)^{-1} from the R factor of the final WLS step,
  # computed the same way as summary.glm
  cov_unscaled <- chol2inv(qr.R(qr_fit))
  cov_unscaled[qr_fit$pivot, qr_fit$pivot] <- cov_unscaled
  dimnames(cov_unscaled) <- list(colnames(X), colnames(X))
  names(coef) <- colnames(X)
  rownames(trace_coef) <- NULL

  list(coefficients = coef, linear.predictors = eta, fitted.values = mu,
       deviance = dev, iter = iter, converged = conv, boundary = boundary,
       working.weights = sw^2, cov.unscaled = cov_unscaled,
       trace = list(deviance = trace_dev, coefficients = trace_coef))
}
