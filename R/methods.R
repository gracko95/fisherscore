# S3 methods for "nbreg" and "zinbreg" objects.

# ---- nbreg -----------------------------------------------------------------------

#' Methods for nbreg fits
#'
#' Standard S3 methods for objects returned by [nbreg()].
#'
#' @param object,x An `"nbreg"` object (or its summary for `print`).
#' @param digits Number of significant digits to print.
#' @param type For `residuals()`: `"response"` (\eqn{y - \hat\mu}) or
#'   `"pearson"` (\eqn{(y - \hat\mu)/\sqrt{\hat\mu + \hat\mu^2/\hat\theta}}).
#'   For `predict()`: `"response"` (\eqn{\hat\mu}) or `"link"`
#'   (\eqn{\hat\eta = \log\hat\mu}).
#' @param newdata Optional data frame to predict for. Offsets are evaluated in
#'   `newdata`.
#' @param ... Not used.
#' @return `coef()` a named vector; `vcov()` a matrix; `logLik()` a
#'   `"logLik"` object with `df` = number of coefficients + 1 (for
#'   \eqn{\theta}), so `AIC()` and `BIC()` work; `fitted()`, `residuals()`
#'   and `predict()` numeric vectors; `nobs()` the number of observations;
#'   `summary()` an object of class `"summary.nbreg"`.
#' @examples
#' set.seed(2)
#' d <- data.frame(x = rnorm(200))
#' d$y <- rnbinom(200, size = 1.5, mu = exp(0.5 + 0.8 * d$x))
#' fit <- nbreg(y ~ x, data = d)
#' fit
#' coef(fit)
#' sqrt(diag(vcov(fit)))
#' AIC(fit); BIC(fit)
#' head(residuals(fit, type = "pearson"))
#' predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
#' @name nbreg-methods
NULL

#' @rdname nbreg-methods
#' @export
print.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Negative binomial (NB2) regression, IRLS + Newton-Raphson for theta\n\n")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L, quote = FALSE)
  cat(sprintf("\nTheta: %s (SE %s)\n", format(x$theta, digits = digits),
              format(x$SE.theta, digits = digits)))
  cat(sprintf("Log-likelihood: %s on %d Df\n", format(x$loglik, digits = digits + 2L),
              length(x$coefficients) + 1L))
  if (!x$converged) cat("Warning: the algorithm did not converge\n")
  cat("\n")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
summary.nbreg <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  structure(list(call = object$call,
                 coefficients = coef_table(object$coefficients, se),
                 theta = object$theta, SE.theta = object$SE.theta,
                 loglik = object$loglik, df = length(object$coefficients) + 1L,
                 aic = stats::AIC(object), iter = object$iter, converged = object$converged,
                 pearson = stats::residuals(object, type = "pearson")),
            class = "summary.nbreg")
}

#' @rdname nbreg-methods
#' @export
print.summary.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Pearson residuals:\n")
  q <- stats::quantile(x$pearson)
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCoefficients (negative binomial with log link):\n")
  stats::printCoefmat(x$coefficients, digits = digits)
  cat(sprintf("\nTheta = %s (SE %s)\n", format(x$theta, digits = digits + 2L),
              format(x$SE.theta, digits = digits)))
  cat(sprintf("Log-likelihood: %s on %d Df,  AIC: %s\n", format(x$loglik, digits = digits + 3L),
              x$df, format(x$aic, digits = digits + 3L)))
  cat(sprintf("Outer iterations (IRLS + Newton-Raphson): %d%s\n\n", x$iter,
              if (x$converged) "" else " (NOT converged)"))
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
coef.nbreg <- function(object, ...) object$coefficients

#' @rdname nbreg-methods
#' @export
vcov.nbreg <- function(object, ...) object$vcov

#' @rdname nbreg-methods
#' @export
logLik.nbreg <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients) + 1L,
            nobs = object$nobs, class = "logLik")
}

#' @rdname nbreg-methods
#' @export
nobs.nbreg <- function(object, ...) object$nobs

#' @rdname nbreg-methods
#' @export
fitted.nbreg <- function(object, ...) object$fitted.values

#' @rdname nbreg-methods
#' @export
residuals.nbreg <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  mu <- object$fitted.values
  r <- object$y - mu
  if (type == "pearson") r <- r / sqrt(mu + mu^2 / object$theta)   # Var(Y) = mu + mu^2/theta
  r
}

#' @rdname nbreg-methods
#' @export
predict.nbreg <- function(object, newdata = NULL, type = c("response", "link"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    eta <- object$linear.predictors
  } else {
    tt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass, xlev = object$xlevels)
    X <- stats::model.matrix(tt, mf, contrasts.arg = object$contrasts)
    eta <- drop(X %*% object$coefficients) + new_offset(mf, object$call, newdata, tt)
  }
  if (type == "link") eta else exp(eta)
}

# Offset for new data: offset() terms in the formula plus a re-evaluated 'offset =' argument
new_offset <- function(mf, call, newdata, tt) {
  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, nrow(mf))
  if (!is.null(call$offset)) off <- off + eval(call$offset, newdata, environment(tt))
  off
}

# ---- zinbreg -----------------------------------------------------------------------

#' Methods for zinbreg fits
#'
#' Standard S3 methods for objects returned by [zinbreg()].
#'
#' @param object,x A `"zinbreg"` object (or its summary for `print`).
#' @param digits Number of significant digits to print.
#' @param model For `coef()` and `vcov()`: `"full"` (both parts, named
#'   `count_*` and `zero_*` as in `pscl::zeroinfl()`), `"count"` or `"zero"`.
#' @param type For `residuals()`: `"response"` (\eqn{y - (1-\hat\pi)\hat\mu})
#'   or `"pearson"` (divided by \eqn{\sqrt{\mathrm{Var}(Y)}} with
#'   \eqn{\mathrm{Var}(Y) = (1-\pi)\mu\{1 + \mu(\pi + 1/\theta)\}}).
#'   For `predict()`: `"response"` (\eqn{(1-\hat\pi)\hat\mu}), `"count"`
#'   (\eqn{\hat\mu}), `"zero"` (\eqn{\hat\pi}) or `"prob0"`
#'   (\eqn{P(Y=0) = \hat\pi + (1-\hat\pi)f_{NB}(0;\hat\mu,\hat\theta)}).
#' @param newdata Optional data frame to predict for.
#' @param ... Not used.
#' @return As for [nbreg-methods]. `logLik()` has `df` = number of count and
#'   zero coefficients + 1.
#' @examples
#' set.seed(3)
#' n <- 400
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' d$y <- ifelse(rbinom(n, 1, plogis(-0.5 + d$z)) == 1, 0,
#'               rnbinom(n, size = 2, mu = exp(1 + 0.5 * d$x)))
#' fit <- zinbreg(y ~ x | z, data = d)
#' coef(fit)
#' coef(fit, model = "zero")
#' AIC(fit)
#' nd <- data.frame(x = 0, z = c(-1, 0, 1))
#' predict(fit, newdata = nd, type = "response")
#' predict(fit, newdata = nd, type = "zero")
#' @name zinbreg-methods
NULL

#' @rdname zinbreg-methods
#' @export
print.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Zero-inflated negative binomial regression, EM with IRLS M-steps\n\n")
  cat("Count model coefficients (negbin with log link):\n")
  print.default(format(x$coefficients$count, digits = digits), print.gap = 2L, quote = FALSE)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  print.default(format(x$coefficients$zero, digits = digits), print.gap = 2L, quote = FALSE)
  cat(sprintf("\nTheta: %s (log(theta) = %s, SE %s)\n", format(x$theta, digits = digits),
              format(log(x$theta), digits = digits), format(x$SE.logtheta, digits = digits)))
  cat(sprintf("Log-likelihood: %s on %d Df\n", format(x$loglik, digits = digits + 2L),
              attr(stats::logLik(x), "df")))
  if (!x$converged) cat("Warning: EM did not converge\n")
  cat("\n")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
summary.zinbreg <- function(object, ...) {
  se <- sqrt(diag(object$vcov_full))
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  count <- coef_table(c(object$coefficients$count, `Log(theta)` = log(object$theta)),
                      se[c(seq_len(p), p + q + 1L)])
  zero <- coef_table(object$coefficients$zero, se[p + seq_len(q)])
  structure(list(call = object$call, coefficients = list(count = count, zero = zero),
                 theta = object$theta, SE.theta = object$SE.theta,
                 loglik = object$loglik, df = p + q + 1L, aic = stats::AIC(object),
                 iter = object$iter, converged = object$converged,
                 monotone = object$monotone,
                 pearson = stats::residuals(object, type = "pearson")),
            class = "summary.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
print.summary.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Pearson residuals:\n")
  q <- stats::quantile(x$pearson)
  names(q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(q, digits = digits)
  cat("\nCount model coefficients (negbin with log link):\n")
  stats::printCoefmat(x$coefficients$count, digits = digits)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  stats::printCoefmat(x$coefficients$zero, digits = digits)
  cat(sprintf("\nTheta = %s (SE %s)\n", format(x$theta, digits = digits + 2L),
              format(x$SE.theta, digits = digits)))
  cat(sprintf("Log-likelihood: %s on %d Df,  AIC: %s\n", format(x$loglik, digits = digits + 3L),
              x$df, format(x$aic, digits = digits + 3L)))
  cat(sprintf("EM iterations: %d%s; log-likelihood never decreased: %s\n\n", x$iter,
              if (x$converged) "" else " (NOT converged)", x$monotone))
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
coef.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  cf <- object$coefficients
  switch(model,
         count = cf$count,
         zero = cf$zero,
         full = c(stats::setNames(cf$count, paste0("count_", names(cf$count))),
                  stats::setNames(cf$zero, paste0("zero_", names(cf$zero)))))
}

#' @rdname zinbreg-methods
#' @export
vcov.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  idx <- switch(model, full = seq_len(p + q), count = seq_len(p), zero = p + seq_len(q))
  v <- object$vcov[idx, idx, drop = FALSE]
  if (model != "full") {
    nm <- names(coef(object, model = model))   # drop the count_/zero_ prefix
    dimnames(v) <- list(nm, nm)
  }
  v
}

#' @rdname zinbreg-methods
#' @export
logLik.zinbreg <- function(object, ...) {
  df <- length(object$coefficients$count) + length(object$coefficients$zero) + 1L
  structure(object$loglik, df = df, nobs = object$nobs, class = "logLik")
}

#' @rdname zinbreg-methods
#' @export
nobs.zinbreg <- function(object, ...) object$nobs

#' @rdname zinbreg-methods
#' @export
fitted.zinbreg <- function(object, ...) object$fitted.values

#' @rdname zinbreg-methods
#' @export
residuals.zinbreg <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  r <- object$y - object$fitted.values
  if (type == "pearson") {
    mu <- object$mu
    pi <- object$pi
    r <- r / sqrt((1 - pi) * mu * (1 + mu * (pi + 1 / object$theta)))
  }
  r
}

#' @rdname zinbreg-methods
#' @export
predict.zinbreg <- function(object, newdata = NULL,
                            type = c("response", "count", "zero", "prob0"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    mu <- object$mu
    pi <- object$pi
  } else {
    ttX <- stats::delete.response(object$terms$count)
    mfX <- stats::model.frame(ttX, newdata, na.action = stats::na.pass,
                              xlev = object$xlevels$count)
    X <- stats::model.matrix(ttX, mfX, contrasts.arg = object$contrasts$count)
    mu <- exp(drop(X %*% object$coefficients$count) + new_offset(mfX, object$call, newdata, ttX))
    ttZ <- object$terms$zero
    mfZ <- stats::model.frame(ttZ, newdata, na.action = stats::na.pass,
                              xlev = object$xlevels$zero)
    Z <- stats::model.matrix(ttZ, mfZ, contrasts.arg = object$contrasts$zero)
    pi <- stats::plogis(drop(Z %*% object$coefficients$zero))
  }
  switch(type,
         response = (1 - pi) * mu,
         count = mu,
         zero = pi,
         prob0 = pi + (1 - pi) * stats::dnbinom(0, size = object$theta, mu = mu))
}
