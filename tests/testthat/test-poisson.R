test_that("poisson links match glm", {
  d <- sim_data()
  expect_matches_glm(yp ~ x1 + x2 + g, poisson("log"), d)
  expect_matches_glm(yp ~ x1 + x2, poisson("sqrt"), d)
  expect_matches_glm(yp ~ x2 + I(x1 > 0), poisson("identity"), d)
  expect_matches_glm(breaks ~ wool * tension, poisson, warpbreaks)
})

test_that("log-linear model for a 3-way contingency table matches glm", {
  ucb <- as.data.frame(UCBAdmissions)
  expect_matches_glm(Freq ~ (Admit + Gender + Dept)^2, poisson, ucb)
})

test_that("offsets (formula and argument) match glm, including the null deviance", {
  d <- sim_data()
  expect_matches_glm(yp ~ x1 + offset(log(exposure)), poisson, d)
  expect_matches_glm(yp ~ x1, poisson, d, offset = log(exposure))
})

test_that("negative counts are rejected", {
  expect_error(fs_glm(y ~ 1, poisson, data.frame(y = c(1, -1, 2))), "negative")
})
