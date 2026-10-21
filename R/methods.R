#' Methods for `fs_glm` objects
#'
#' Standard extractor and inference methods, which behave like their `glm`
#' counterparts.
#'
#' @param object,x An `fs_glm` object.
#' @param ... Unused, or further `fs_glm` objects for `anova()`.
#' @name fs_glm-methods
NULL

#' @rdname fs_glm-methods
#' @export
print.fs_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Fisher scoring fit (", x$engine, " engine), ", x$family$family, " family, ",
      x$family$link, " link\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", x$df.residual, "Residual\n")
  cat("Null Deviance:    ", format(signif(x$null.deviance, digits)),
      "\nResidual Deviance:", format(signif(x$deviance, digits)),
      "\tAIC:", format(signif(x$aic, digits)), "\n")
  invisible(x)
}

#' @rdname fs_glm-methods
#' @param dispersion Optional dispersion to use instead of the estimate.
#' @export
summary.fs_glm <- function(object, dispersion = NULL, ...) {
  disp <- if (is.null(dispersion)) object$dispersion else dispersion
  est <- object$coefficients[!is.na(object$coefficients)]
  se <- sqrt(diag(object$cov.unscaled) * disp)
  stat <- est / se
  # As in summary.glm(): t tests only when the dispersion is estimated
  if (object$family$fixed_dispersion || !is.null(dispersion)) {
    p <- 2 * stats::pnorm(-abs(stat))
    cn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  } else {
    p <- 2 * stats::pt(-abs(stat), object$df.residual)
    cn <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  }
  coefs <- cbind(est, se, stat, p)
  dimnames(coefs) <- list(names(est), cn)
  structure(list(call = object$call, family = object$family, engine = object$engine,
                 coefficients = coefs, aliased = object$aliased, dispersion = disp,
                 deviance = object$deviance, null.deviance = object$null.deviance,
                 df.residual = object$df.residual, df.null = object$df.null,
                 aic = object$aic, iter = object$iter, converged = object$converged,
                 deviance.resid = stats::residuals(object, "deviance"),
                 cov.scaled = object$cov.unscaled * disp),
            class = "summary.fs_glm")
}

#' @rdname fs_glm-methods
#' @param digits Number of significant digits to print.
#' @export
print.summary.fs_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Deviance Residuals:\n")
  print(stats::setNames(stats::quantile(x$deviance.resid, na.rm = TRUE),
                        c("Min", "1Q", "Median", "3Q", "Max")), digits = digits)
  if (length(x$aliased))
    cat("\nCoefficients: (", length(x$aliased),
        " not defined because of singularities)\n", sep = "")
  else cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, na.print = "NA")
  cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
      format(x$dispersion), ")\n\n", sep = "")
  cat("    Null deviance: ", format(x$null.deviance, digits = max(5L, digits + 1L)),
      " on ", x$df.null, " degrees of freedom\n", sep = "")
  cat("Residual deviance: ", format(x$deviance, digits = max(5L, digits + 1L)),
      " on ", x$df.residual, " degrees of freedom\n", sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)), "\n\n", sep = "")
  cat("Number of Fisher Scoring iterations: ", x$iter, " (", x$engine, " engine",
      if (!x$converged) ", NOT converged", ")\n\n", sep = "")
  invisible(x)
}

#' @rdname fs_glm-methods
#' @export
coef.fs_glm <- function(object, ...) object$coefficients

#' @rdname fs_glm-methods
#' @export
vcov.fs_glm <- function(object, ...) object$cov.unscaled * object$dispersion

#' @rdname fs_glm-methods
#' @export
fitted.fs_glm <- function(object, ...) stats::napredict(object$na.action, object$fitted.values)

#' @rdname fs_glm-methods
#' @export
deviance.fs_glm <- function(object, ...) object$deviance

#' @rdname fs_glm-methods
#' @export
nobs.fs_glm <- function(object, ...) sum(object$weights != 0)

#' @rdname fs_glm-methods
#' @export
formula.fs_glm <- function(x, ...) stats::formula(x$terms)

#' @rdname fs_glm-methods
#' @param type For `residuals()`: `"deviance"`, `"pearson"`, `"working"` or
#'   `"response"`. For `predict()`: `"link"` or `"response"`.
#' @export
residuals.fs_glm <- function(object, type = c("deviance", "pearson", "working", "response"),
                             ...) {
  type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  wt <- object$weights
  fam <- object$family
  res <- switch(type,
    deviance = sign(y - mu) * sqrt(pmax(fam$dev.resids(y, mu, wt), 0)),
    pearson = (y - mu) * sqrt(wt) / sqrt(fam$variance(mu)),
    working = object$residuals,
    response = y - mu)
  stats::naresid(object$na.action, res)
}

#' @rdname fs_glm-methods
#' @export
logLik.fs_glm <- function(object, ...) {
  p <- object$rank + as.integer(!object$family$fixed_dispersion)
  structure(p - object$aic / 2, df = p, nobs = stats::nobs(object), class = "logLik")
}

#' @rdname fs_glm-methods
#' @param newdata Optional data frame to predict for.
#' @param se.fit Logical; return standard errors (delta method on the
#'   response scale)?
#' @export
predict.fs_glm <- function(object, newdata = NULL, type = c("link", "response"),
                           se.fit = FALSE, ...) {
  type <- match.arg(type)
  beta <- object$coefficients
  keep <- !is.na(beta)
  if (is.null(newdata)) {
    if (is.null(object$model)) stop("refit with model = TRUE or supply 'newdata'", call. = FALSE)
    X <- stats::model.matrix(object$terms, object$model, contrasts.arg = object$contrasts)
    offset <- object$offset
  } else {
    tt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(tt, newdata, xlev = object$xlevels)
    X <- stats::model.matrix(tt, mf, contrasts.arg = object$contrasts)
    # Offsets from the formula are in the model frame; an `offset =` argument
    # is re-evaluated in newdata, as predict.lm() does.
    offset <- stats::model.offset(mf)
    if (is.null(offset)) offset <- rep.int(0, nrow(X))
    if (!is.null(object$call$offset))
      offset <- offset + eval(object$call$offset, newdata, environment(object$terms))
  }
  X <- X[, keep, drop = FALSE]
  eta <- drop(X %*% beta[keep]) + offset
  fit <- if (type == "link") eta else object$family$linkinv(eta)
  # Pad rows dropped by na.exclude back in, only when predicting for the fit data
  pad <- function(v) if (is.null(newdata)) stats::napredict(object$na.action, v) else v
  if (!se.fit) return(pad(fit))
  se <- sqrt(rowSums((X %*% stats::vcov(object)) * X))
  if (type == "response") se <- se * abs(object$family$mu.eta(eta))
  list(fit = pad(fit), se.fit = pad(se), residual.scale = sqrt(object$dispersion))
}

#' @rdname fs_glm-methods
#' @param parm Coefficients to give intervals for (default: all).
#' @param level Confidence level.
#' @export
confint.fs_glm <- function(object, parm, level = 0.95, ...) {
  s <- summary(object)$coefficients
  if (missing(parm)) parm <- rownames(s)
  q <- if (object$family$fixed_dispersion) stats::qnorm((1 + level) / 2)
       else stats::qt((1 + level) / 2, object$df.residual)
  ci <- cbind(s[parm, 1] - q * s[parm, 2], s[parm, 1] + q * s[parm, 2])
  pct <- paste(format(100 * c(1 - level, 1 + level) / 2, trim = TRUE, digits = 3), "%")
  dimnames(ci) <- list(parm, pct)
  ci
}

#' @rdname fs_glm-methods
#' @details `anova()` with two or more nested fits gives likelihood-ratio
#'   (deviance) tests, using a chi-squared reference distribution, or F tests
#'   when the dispersion is estimated (gaussian family).
#' @export
anova.fs_glm <- function(object, ...) {
  fits <- c(list(object), list(...))
  if (length(fits) < 2) stop("supply at least two nested fs_glm fits", call. = FALSE)
  if (!all(vapply(fits, inherits, logical(1), "fs_glm")))
    stop("all arguments must be fs_glm fits", call. = FALSE)
  rdf <- vapply(fits, `[[`, numeric(1), "df.residual")
  rdev <- vapply(fits, `[[`, numeric(1), "deviance")
  df <- c(NA, -diff(rdf))
  dd <- c(NA, -diff(rdev))
  tab <- data.frame(`Resid. Df` = rdf, `Resid. Dev` = rdev, Df = df, Deviance = dd,
                    check.names = FALSE)
  big <- fits[[which.min(rdf)]]
  if (big$family$fixed_dispersion) {
    tab[["Pr(>Chi)"]] <- stats::pchisq(abs(dd), abs(df), lower.tail = FALSE)
    heading <- "Analysis of Deviance Table (likelihood ratio tests)\n"
  } else {
    Fv <- (dd / df) / big$dispersion
    tab$F <- Fv
    tab[["Pr(>F)"]] <- stats::pf(abs(Fv), abs(df), big$df.residual, lower.tail = FALSE)
    heading <- "Analysis of Deviance Table (F tests)\n"
  }
  rownames(tab) <- paste("Model", seq_along(fits))
  forms <- vapply(fits, function(f) paste(deparse(formula(f)), collapse = " "), character(1))
  structure(tab, heading = c(heading, paste0("Model ", seq_along(forms), ": ", forms, "\n")),
            class = c("anova", "data.frame"))
}

#' Plot the convergence of Fisher scoring
#'
#' Plots the deviance (or the change in deviance, on a log scale) at each
#' Fisher scoring iteration.
#'
#' @param fit An `fs_glm` object.
#' @param what `"deviance"` or `"change"` (relative change in deviance,
#'   compared with the convergence tolerance).
#' @param ... Further arguments passed to [graphics::plot()].
#' @return The trace, as a data frame, invisibly.
#' @examples
#' fit <- fs_glm(cbind(ncases, ncontrols) ~ agegp + alcgp, binomial("probit"), esoph)
#' plot_convergence(fit, "change")
#' @export
plot_convergence <- function(fit, what = c("deviance", "change"), ...) {
  what <- match.arg(what)
  dev <- fit$trace$deviance
  it <- seq_along(dev)
  rel <- c(NA, abs(diff(dev)) / (abs(dev[-1]) + 0.1))
  if (what == "deviance") {
    graphics::plot(it, dev, type = "b", pch = 19, xlab = "Iteration", ylab = "Deviance",
                   main = "Fisher scoring: deviance by iteration", ...)
  } else {
    graphics::plot(it[-1], rel[-1], type = "b", pch = 19, log = "y", xlab = "Iteration",
                   ylab = "|dD| / (|D| + 0.1)", main = "Fisher scoring: convergence", ...)
    graphics::abline(h = fit$control$epsilon, lty = 2, col = "grey40")
  }
  invisible(data.frame(iteration = it, deviance = dev, rel_change = rel))
}
