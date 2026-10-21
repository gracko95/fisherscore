# Fit with glm() and both fs_glm() engines and check that they agree.
# Extra arguments (weights, offset, start) are passed on unevaluated, so they
# can refer to columns of `data`, as in a direct glm() call.
expect_matches_glm <- function(formula, family, data, tol = 1e-8, ...) {
  cl <- match.call(expand.dots = TRUE)
  cl$tol <- NULL
  cl[[1L]] <- quote(stats::glm)
  g <- suppressWarnings(eval(cl, parent.frame()))
  for (engine in c("R", "cpp")) {
    cl[[1L]] <- quote(fs_glm)
    cl$engine <- engine
    f <- suppressWarnings(eval(cl, parent.frame()))
    testthat::expect_equal(coef(f), coef(g), tolerance = tol, label = paste(engine, "coef"))
    testthat::expect_equal(sqrt(diag(vcov(f))), sqrt(diag(vcov(g))), tolerance = tol,
                           label = paste(engine, "SE"))
    testthat::expect_equal(deviance(f), deviance(g), tolerance = tol)
    testthat::expect_equal(f$null.deviance, g$null.deviance, tolerance = tol)
    testthat::expect_equal(AIC(f), AIC(g), tolerance = tol)
    testthat::expect_equal(f$df.residual, g$df.residual)
    testthat::expect_equal(f$iter, g$iter)
    testthat::expect_equal(unname(fitted(f)), unname(fitted(g)), tolerance = tol)
  }
  invisible(g)
}

sim_data <- function(n = 400, seed = 651) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = runif(n), g = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$yb <- rbinom(n, 1, pnorm(-0.3 + 0.8 * d$x1 - d$x2))
  d$yp <- rpois(n, (1.5 + 0.5 * d$x2 + 0.3 * (d$x1 > 0))^2)
  d$yg <- 1 + 2 * d$x1 - d$x2 + rnorm(n)
  d$exposure <- rpois(n, 5) + 1
  d
}
