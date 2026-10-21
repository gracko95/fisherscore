test_that("complete separation warns, as glm does", {
  d <- data.frame(x = 1:10, y = rep(0:1, each = 5))
  for (e in c("R", "cpp"))
    expect_warning(
      expect_warning(fs_glm(y ~ x, binomial, d, engine = e), "did not converge"),
      "numerically 0 or 1")
})

test_that("rank-deficient designs give NA coefficients like glm", {
  d <- sim_data()
  d$x3 <- 2 * d$x1
  g <- glm(yb ~ x1 + x3 + x2, binomial, d)
  f <- fs_glm(yb ~ x1 + x3 + x2, binomial, d)
  expect_true(is.na(coef(f)[["x3"]]))
  expect_equal(coef(f), coef(g), tolerance = 1e-8)
  expect_equal(f$rank, g$rank)
  expect_equal(f$aliased, "x3")
})

test_that("non-convergence is reported", {
  d <- sim_data()
  expect_warning(f <- fs_glm(yb ~ x1, binomial, d, control = fs_control(maxit = 1)),
                 "did not converge")
  expect_false(f$converged)
})

test_that("zero weights drop observations as in glm", {
  d <- sim_data()
  w <- rep(c(1, 0, 2), length.out = nrow(d))
  expect_matches_glm(yp ~ x1 + x2, poisson, d, weights = w)
})

test_that("user-supplied starting values are honoured", {
  d <- sim_data()
  g <- glm(yp ~ x1, poisson("identity"), d, start = c(3, 0))
  f <- fs_glm(yp ~ x1, poisson("identity"), d, start = c(3, 0))
  expect_equal(coef(f), coef(g), tolerance = 1e-8)
  expect_error(fs_glm(yp ~ x1, poisson, d, start = 1), "length of 'start'")
})

test_that("nearly collinear predictors are kept, as glm keeps them", {
  d <- sim_data()
  set.seed(1)
  d$z <- d$x1 + 1e-8 * rnorm(nrow(d))
  g <- glm(yg ~ x1 + z, gaussian, d)
  for (e in c("R", "cpp")) {
    f <- fs_glm(yg ~ x1 + z, gaussian, d, engine = e)
    expect_equal(f$rank, g$rank)
    expect_length(f$aliased, 0)
    expect_equal(coef(f), coef(g), tolerance = 1e-4)
  }
})

test_that("singularity checks do not depend on predictor units", {
  d <- sim_data()
  d$x_small <- d$x1 * 1e-8
  expect_matches_glm(yp ~ x_small + x2, poisson, d)
})

test_that("empty models fit with eta = offset, as in glm", {
  d <- sim_data()
  for (e in c("R", "cpp")) {
    g <- glm(yp ~ 0, poisson, d)
    f <- fs_glm(yp ~ 0, poisson, d, engine = e)
    expect_length(coef(f), 0)
    expect_equal(deviance(f), deviance(g))
    expect_equal(AIC(f), AIC(g))
    expect_equal(f$df.residual, g$df.residual)
    expect_equal(f$iter, 0L)
    g <- glm(yp ~ 0 + offset(log(exposure)), poisson, d)
    f <- fs_glm(yp ~ 0 + offset(log(exposure)), poisson, d, engine = e)
    expect_equal(unname(fitted(f)), unname(fitted(g)))
  }
})

test_that("extremely small predictor scales still fit correctly in both engines", {
  d <- data.frame(x = (1:6) * 1e-16, y = c(1, 3, 2, 5, 4, 8))
  expect_matches_glm(y ~ x, poisson, d, tol = 1e-6)
})
