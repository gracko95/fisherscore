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
