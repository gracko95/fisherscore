#' Compare Fisher scoring fits with `stats::glm()`
#'
#' Fits the same model with [stats::glm()] and with [fs_glm()] using both
#' engines, then reports the agreement in estimates, standard errors,
#' deviance and AIC, together with the iteration counts and timings.
#'
#' @param formula,family,data,... Passed to both [stats::glm()] and [fs_glm()].
#' @param engines Which `fs_glm` engines to include.
#' @param times Number of times to repeat each fit when timing it (the
#'   median time is reported).
#' @return An object of class `fs_comparison`, a list with:
#'   * `coefficients`: estimates and standard errors side by side
#'   * `summary`: one row per method with deviance, AIC, iterations, median
#'     time, and the largest absolute difference from `glm` in the
#'     coefficients and standard errors
#'   * `fits`: the fitted models
#' @examples
#' compare_glm(cbind(ncases, ncontrols) ~ agegp + alcgp, binomial, esoph)
#' @export
compare_glm <- function(formula, family = gaussian, data, ..., engines = c("R", "cpp"),
                        times = 5) {
  engines <- match.arg(engines, c("R", "cpp"), several.ok = TRUE)
  # Rebuild the call so arguments such as `weights = n` are evaluated in
  # `data`, exactly as they would be in a direct glm() call.
  cl <- match.call(expand.dots = TRUE)
  cl$engines <- cl$times <- NULL
  env <- parent.frame()
  med_time <- function(call) {
    stats::median(vapply(seq_len(times), function(i) system.time(eval(call, env))[["elapsed"]],
                         numeric(1)))
  }
  cl[[1L]] <- quote(stats::glm)
  g <- eval(cl, env)
  fits <- list(glm = g)
  t <- c(glm = med_time(cl))
  for (e in engines) {
    cl[[1L]] <- quote(nbzinbKoo::fs_glm)
    cl$engine <- e
    fits[[paste0("fs_", e)]] <- eval(cl, env)
    t[paste0("fs_", e)] <- med_time(cl)
  }

  g_coef <- stats::coef(g)
  g_se <- rep(NA_real_, length(g_coef))
  names(g_se) <- names(g_coef)
  sg <- summary(g)$coefficients
  g_se[rownames(sg)] <- sg[, 2]
  se_of <- function(f) {
    if (inherits(f, "glm")) return(g_se)
    s <- summary(f)$coefficients
    out <- stats::setNames(rep(NA_real_, length(f$coefficients)), names(f$coefficients))
    out[rownames(s)] <- s[, 2]
    out
  }

  coefs <- data.frame(term = names(g_coef), row.names = NULL)
  for (nm in names(fits)) {
    coefs[[paste0(nm, "_est")]] <- unname(stats::coef(fits[[nm]]))
    coefs[[paste0(nm, "_se")]] <- unname(se_of(fits[[nm]]))
  }

  maxdiff <- function(a, b) max(abs(a - b), na.rm = TRUE)
  summ <- data.frame(
    method = names(fits),
    deviance = vapply(fits, stats::deviance, numeric(1)),
    aic = vapply(fits, stats::AIC, numeric(1)),
    iterations = vapply(fits, function(f) as.numeric(f$iter), numeric(1)),
    converged = vapply(fits, function(f) isTRUE(f$converged), logical(1)),
    time_sec = unname(t[names(fits)]),
    max_coef_diff = vapply(fits, function(f) maxdiff(stats::coef(f), g_coef), numeric(1)),
    max_se_diff = vapply(fits, function(f) maxdiff(se_of(f), g_se), numeric(1)),
    row.names = NULL
  )
  structure(list(coefficients = coefs, summary = summ, fits = fits), class = "fs_comparison")
}

#' @export
print.fs_comparison <- function(x, digits = 6, ...) {
  cat("Coefficients and standard errors:\n")
  print(x$coefficients, digits = digits, row.names = FALSE)
  cat("\nFit summary (differences are relative to stats::glm):\n")
  s <- x$summary
  s$max_coef_diff <- format(s$max_coef_diff, digits = 3, scientific = TRUE)
  s$max_se_diff <- format(s$max_se_diff, digits = 3, scientific = TRUE)
  print(s, digits = digits, row.names = FALSE)
  invisible(x)
}
