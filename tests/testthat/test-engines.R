test_that("R and C++ engines agree to machine precision", {
  d <- sim_data(1000)
  cases <- list(list(yb ~ x1 + x2 + g, binomial("probit")),
                list(yp ~ x1 + g, poisson("log")),
                list(yp ~ x1 + x2, poisson("sqrt")),
                list(yg ~ x1 + x2, gaussian()))
  for (cc in cases) {
    r <- fs_glm(cc[[1]], cc[[2]], d, engine = "R")
    cpp <- fs_glm(cc[[1]], cc[[2]], d, engine = "cpp")
    expect_equal(coef(cpp), coef(r), tolerance = 1e-10)
    expect_equal(vcov(cpp), vcov(r), tolerance = 1e-10)
    expect_equal(cpp$trace$deviance, r$trace$deviance, tolerance = 1e-10)
    expect_identical(cpp$iter, r$iter)
  }
})

test_that("the trace records one row per iteration", {
  f <- fs_glm(cbind(ncases, ncontrols) ~ agegp, binomial, esoph, engine = "cpp")
  expect_length(f$trace$deviance, f$iter)
  expect_equal(nrow(f$trace$coefficients), f$iter)
  expect_equal(unname(f$trace$coefficients[f$iter, ]), unname(coef(f)))
})
