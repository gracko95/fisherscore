#' Fit a generalized linear model by Fisher scoring
#'
#' `fs_glm()` is a drop-in counterpart to [stats::glm()] whose estimation
#' runs entirely through this package's own Fisher scoring code. You can
#' choose between a pure-R engine ([fs_fit_r()]) and a C++ engine
#' ([fs_fit_cpp()]). Both use the same starting values, convergence rule
#' and step-halving as `glm()`, so estimates, standard errors, deviances and
#' iteration counts can be compared exactly.
#'
#' @param formula A model formula, as for [stats::glm()]. Binomial responses
#'   may be 0/1, a factor (first level = failure), proportions with
#'   `weights` giving the number of trials, or `cbind(successes, failures)`.
#' @param family `binomial` (logit, probit, cloglog), `poisson` (log,
#'   identity, sqrt) or `gaussian` (identity). Accepts a family object,
#'   function or name.
#' @param data A data frame.
#' @param weights Optional prior weights.
#' @param subset Optional subset of observations.
#' @param na.action How to handle missing values (default `na.omit`).
#' @param offset Optional offset added to the linear predictor.
#' @param start Optional starting coefficients.
#' @param engine `"R"` (pure R, the default) or `"cpp"` (RcppArmadillo).
#' @param control A list from [fs_control()].
#' @param model Logical; keep the model frame in the result?
#' @param ... Further arguments passed to [fs_control()] when `control` is
#'   not supplied, e.g. `epsilon`, `maxit`, `trace`.
#' @return An object of class `fs_glm`. Its main components are the same
#'   as a `glm` object's: `coefficients`, `fitted.values`, `linear.predictors`,
#'   `deviance`, `null.deviance`, `aic`, `df.residual`, `df.null`, `iter`,
#'   `converged`. It also holds `cov.unscaled`, `dispersion`, `engine`,
#'   `time` (elapsed seconds in the engine) and `trace`, the deviance and
#'   coefficients at each iteration.
#' @seealso [compare_glm()] to compare a fit with [stats::glm()].
#' @examples
#' fit <- fs_glm(cbind(ncases, ncontrols) ~ agegp + alcgp, family = binomial, data = esoph)
#' summary(fit)
#'
#' # Poisson log-linear model
#' fit2 <- fs_glm(breaks ~ wool * tension, family = poisson, data = warpbreaks, engine = "cpp")
#' coef(fit2)
#' @export
fs_glm <- function(formula, family = gaussian, data, weights, subset,
                   na.action, offset, start = NULL, engine = c("R", "cpp"),
                   control = list(...), model = TRUE, ...) {
  call <- match.call()
  engine <- match.arg(engine)
  control <- do.call(fs_control, control)
  fam <- fs_family(family)

  # ---- model frame (same machinery glm uses) -----------------------------
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  Y <- stats::model.response(mf, "any")
  X <- stats::model.matrix(mt, mf)
  nobs <- NROW(Y)
  prior <- as.vector(stats::model.weights(mf))
  if (is.null(prior)) prior <- rep.int(1, nobs)
  if (any(prior < 0)) stop("negative weights not allowed", call. = FALSE)
  offset <- as.vector(stats::model.offset(mf))
  if (is.null(offset)) offset <- rep.int(0, nobs)

  # ---- response set-up (mirrors family()$initialize) ----------------------
  r <- fs_response(Y, prior, fam)
  y <- r$y
  weights <- r$weights
  n <- r$n
  if (!is.null(names(Y)) && is.null(dim(Y))) names(y) <- names(Y)
  else names(y) <- rownames(mf)

  # ---- aliased (linearly dependent) columns get NA coefficients ----------
  good <- weights > 0
  qx <- qr(X[good, , drop = FALSE])
  keep <- sort(qx$pivot[seq_len(qx$rank)])
  aliased <- setdiff(seq_len(ncol(X)), keep)
  Xk <- X[, keep, drop = FALSE]
  if (!is.null(start)) {
    if (length(start) != ncol(X))
      stop(sprintf("length of 'start' should equal %d", ncol(X)), call. = FALSE)
    start <- start[keep]
  }

  # ---- starting values ---------------------------------------------------
  eta_start <- if (!is.null(start)) drop(Xk %*% start) + offset else fam$linkfun(r$mustart)

  # ---- run the chosen engine on the rows with positive weight -------------
  timing <- system.time(
    fit <- fs_engine(engine, Xk[good, , drop = FALSE], y[good], weights[good],
                     offset[good], fam, control, eta_start[good], start)
  )[["elapsed"]]
  if (!fit$converged) warning("fs_glm: algorithm did not converge", call. = FALSE)
  if (fit$boundary) warning("fs_glm: algorithm stopped at boundary value", call. = FALSE)

  coef <- stats::setNames(rep(NA_real_, ncol(X)), colnames(X))
  coef[keep] <- fit$coefficients
  eta <- drop(Xk %*% fit$coefficients) + offset
  mu <- fam$linkinv(eta)
  names(eta) <- names(mu) <- names(y)

  eps <- 10 * .Machine$double.eps
  if (fam$family == "binomial" && (any(mu[good] > 1 - eps) || any(mu[good] < eps)))
    warning("fs_glm: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
  if (fam$family == "poisson" && any(mu[good] < eps))
    warning("fs_glm: fitted rates numerically 0 occurred", call. = FALSE)

  # ---- null deviance -----------------------------------------------------
  intercept <- attr(mt, "intercept") > 0L
  if (intercept && any(offset != 0)) {
    one <- matrix(1, sum(good), 1, dimnames = list(NULL, "(Intercept)"))
    null_fit <- fs_engine(engine, one, y[good], weights[good], offset[good], fam,
                          control, fam$linkfun(r$mustart[good]), NULL)
    null_dev <- null_fit$deviance
  } else {
    wtdmu <- if (intercept) sum(weights * y) / sum(weights) else fam$linkinv(offset)
    null_dev <- sum(fam$dev.resids(y, wtdmu, weights))
  }

  # ---- deviance, AIC, dispersion -----------------------------------------
  rank <- length(keep)
  n_ok <- sum(good)
  dev <- sum(fam$dev.resids(y, mu, weights))
  aic <- fam$aic(y, n, mu, weights, dev) + 2 * rank
  df_res <- n_ok - rank
  df_null <- n_ok - as.integer(intercept)
  pearson <- sum((weights * (y - mu)^2 / fam$variance(mu))[good])
  dispersion <- if (fam$fixed_dispersion) 1 else if (df_res > 0) pearson / df_res else NaN

  cov_unscaled <- fit$cov.unscaled
  dimnames(cov_unscaled) <- list(colnames(Xk), colnames(Xk))

  structure(list(
    coefficients = coef, fitted.values = mu, linear.predictors = eta,
    residuals = (y - mu) / fam$mu.eta(eta),
    deviance = dev, null.deviance = null_dev, aic = aic,
    rank = rank, df.residual = df_res, df.null = df_null,
    iter = fit$iter, converged = fit$converged, boundary = fit$boundary,
    prior.weights = prior, weights = weights, y = y, n = n, offset = offset,
    cov.unscaled = cov_unscaled, dispersion = dispersion, pearson.chisq = pearson,
    aliased = colnames(X)[aliased],
    family = fam, engine = engine, time = timing, control = control,
    trace = fit$trace, call = call, formula = formula, terms = mt,
    xlevels = stats::.getXlevels(mt, mf), contrasts = attr(X, "contrasts"),
    model = if (model) mf
  ), class = "fs_glm")
}

# Dispatch to the R or C++ engine with a common interface.
fs_engine <- function(engine, X, y, weights, offset, fam, control, eta_start, start) {
  if (engine == "R") {
    fs_fit_r(X, y, weights, offset, fam, control, eta_start, start)
  } else {
    fit <- fs_fit_cpp(X, y, weights, offset, eta_start,
                      if (is.null(start)) numeric(0) else start,
                      fam$family_code, fam$link_code,
                      control$epsilon, control$maxit, control$trace)
    names(fit$coefficients) <- colnames(X)
    colnames(fit$trace$coefficients) <- colnames(X)
    fit
  }
}

# Turn the raw response into (y, weights, n, mustart), as each family's
# initialize expression does in stats.
fs_response <- function(Y, prior, fam) {
  nobs <- NROW(Y)
  if (fam$family == "binomial") {
    if (NCOL(Y) == 1) {
      if (is.factor(Y)) Y <- Y != levels(Y)[1L]
      y <- as.numeric(Y)
      y[prior == 0] <- 0
      if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1", call. = FALSE)
      if (any(abs(prior * y - round(prior * y)) > 0.001))
        warning("non-integer #successes in a binomial fs_glm!", call. = FALSE)
      return(list(y = y, weights = prior, n = rep.int(1, nobs),
                  mustart = (prior * y + 0.5) / (prior + 1)))
    }
    if (NCOL(Y) == 2) {
      if (any(abs(Y - round(Y)) > 0.001))
        warning("non-integer counts in a binomial fs_glm!", call. = FALSE)
      n <- Y[, 1L] + Y[, 2L]
      y <- ifelse(n == 0, 0, Y[, 1L] / n)
      return(list(y = y, weights = prior * n, n = n, mustart = (n * y + 0.5) / (n + 1)))
    }
    stop("binomial response must be 0/1, a factor, proportions, or cbind(successes, failures)",
         call. = FALSE)
  }
  y <- as.numeric(Y)
  if (fam$family == "poisson") {
    if (any(y < 0)) stop("negative values not allowed for the 'Poisson' family", call. = FALSE)
    return(list(y = y, weights = prior, n = rep.int(1, nobs), mustart = y + 0.1))
  }
  list(y = y, weights = prior, n = rep.int(1, nobs), mustart = y)
}
