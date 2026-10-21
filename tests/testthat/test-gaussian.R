test_that("gaussian identity matches glm and lm", {
  d <- sim_data()
  expect_matches_glm(yg ~ x1 + x2 + g, gaussian, d)
  expect_matches_glm(mpg ~ wt + hp + factor(cyl), gaussian, mtcars)
  f <- fs_glm(mpg ~ wt + hp, gaussian, mtcars)
  expect_equal(f$dispersion, summary(lm(mpg ~ wt + hp, mtcars))$sigma^2)
})

test_that("unsupported family/link combinations error clearly", {
  expect_error(fs_glm(mpg ~ wt, gaussian("log"), mtcars), "not supported")
  expect_error(fs_glm(mpg ~ wt, Gamma(), mtcars), "not supported")
})
