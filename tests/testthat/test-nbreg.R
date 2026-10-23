# nbreg() against MASS::glm.nb(), with the tolerances of Table 4:
# coefficients and SEs <= 1e-5, theta and log-likelihood <= 1e-4.

nb_reference <- function(formula, data, ...) {
  MASS::glm.nb(formula, data = data, ...,
               control = stats::glm.control(epsilon = 1e-10, maxit = 100))
}

expect_matches_glm_nb <- function(fit, ref) {
  expect_true(fit$converged)
  expect_lte(max(abs(coef(fit) - coef(ref))), 1e-5)
  expect_lte(max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(ref))))), 1e-5)
  expect_lte(abs(fit$theta - ref$theta), 1e-4)
  expect_lte(abs(fit$SE.theta - ref$SE.theta), 1e-4)
  expect_lte(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-4)
}

test_that("nbreg matches glm.nb on bioChemists", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment
  expect_matches_glm_nb(nbreg(f, bioChemists), nb_reference(f, bioChemists))
})

test_that("nbreg matches glm.nb on NMES1988", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER", envir = environment())
  f <- visits ~ hospital + health + chronic + gender + school + insurance
  expect_matches_glm_nb(nbreg(f, NMES1988), nb_reference(f, NMES1988))
})

test_that("nbreg handles an offset (argument and formula term) like glm.nb", {
  skip_if_not_installed("MASS")
  set.seed(11)
  d <- data.frame(x = rnorm(300), t = runif(300, 0.5, 4))
  d$y <- rnbinom(300, size = 1.5, mu = d$t * exp(0.3 + 0.6 * d$x))
  ref <- nb_reference(y ~ x + offset(log(t)), d)
  fit1 <- nbreg(y ~ x, data = d, offset = log(t))
  fit2 <- nbreg(y ~ x + offset(log(t)), data = d)
  expect_matches_glm_nb(fit1, ref)
  expect_matches_glm_nb(fit2, ref)
  # predict() re-evaluates the offset in newdata
  nd <- data.frame(x = c(0, 1), t = c(1, 2))
  expect_equal(predict(fit1, newdata = nd), predict(ref, newdata = nd, type = "response"),
               tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("nbreg fits a single-covariate and an intercept-only model", {
  skip_if_not_installed("MASS")
  set.seed(12)
  d <- data.frame(x = rnorm(250))
  d$y <- rnbinom(250, size = 3, mu = exp(1 - 0.4 * d$x))
  expect_matches_glm_nb(nbreg(y ~ x, d), nb_reference(y ~ x, d))
  fit0 <- nbreg(y ~ 1, d)
  expect_equal(unname(coef(fit0)), log(mean(d$y)), tolerance = 1e-8)
})

test_that("nbreg methods return the right shapes and values", {
  set.seed(13)
  d <- data.frame(x = rnorm(200), g = gl(2, 100, labels = c("a", "b")))
  d$y <- rnbinom(200, size = 2, mu = exp(0.5 + 0.5 * d$x + 0.3 * (d$g == "b")))
  fit <- nbreg(y ~ x + g, d)

  expect_s3_class(fit, "nbreg")
  expect_named(coef(fit), c("(Intercept)", "x", "gb"))
  expect_equal(dim(vcov(fit)), c(3L, 3L))
  expect_equal(nobs(fit), 200L)
  expect_equal(attr(logLik(fit), "df"), 4L)
  expect_equal(AIC(fit), -2 * fit$loglik + 2 * 4)
  expect_equal(BIC(fit), -2 * fit$loglik + log(200) * 4)
  expect_equal(fitted(fit) + residuals(fit), d$y, ignore_attr = TRUE)
  mu <- fitted(fit)
  expect_equal(residuals(fit, "pearson"), (d$y - mu) / sqrt(mu + mu^2 / fit$theta))
  expect_equal(predict(fit), mu)
  expect_equal(predict(fit, type = "link"), log(mu))
  expect_equal(predict(fit, newdata = d[1:5, ]), mu[1:5], ignore_attr = TRUE)
  expect_equal(predict(fit, newdata = data.frame(x = 0, g = "b")),
               exp(sum(coef(fit)[c(1, 3)])), ignore_attr = TRUE)

  s <- summary(fit)
  expect_s3_class(s, "summary.nbreg")
  expect_equal(colnames(s$coefficients), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_output(print(s), "Theta")
  expect_output(print(fit), "Coefficients")
})

test_that("nbreg gives informative errors for bad input", {
  d <- data.frame(x = 1:10, y = c(0, 1, 2, 3, 0, 1, 5, 2, 1, 0))
  expect_error(nbreg("y ~ x", d), "model formula")
  expect_error(nbreg(y ~ x, data = as.matrix(d)), "data frame")
  expect_error(nbreg(I(-y) ~ x, d), "non-negative integer")
  expect_error(nbreg(I(y + 0.5) ~ x, d), "non-negative integer")
  expect_error(nbreg(y ~ x + I(2 * x), d), "rank deficient")
  expect_error(nbreg(y ~ x, d, offset = rep(0, 3)), "")
  expect_error(nbreg_control(tol = -1), "finite number > 0")
  expect_error(nbreg_control(maxit = 0), "maxit")
})

test_that("nbreg warns at maxit and in the Poisson limit", {
  set.seed(14)
  d <- data.frame(x = rnorm(200))
  d$y <- rnbinom(200, size = 2, mu = exp(1 + 0.5 * d$x))
  expect_warning(fit <- nbreg(y ~ x, d, control = nbreg_control(maxit = 1)), "maxit")
  expect_false(fit$converged)

  d$y <- rpois(200, exp(1 + 0.5 * d$x))     # no overdispersion: theta -> Inf
  expect_warning(fit <- nbreg(y ~ x, d, control = nbreg_control(theta_max = 1e3)), "Poisson limit")
  expect_equal(fit$theta, 1e3)
})

test_that("trace prints the log-likelihood", {
  set.seed(15)
  d <- data.frame(x = rnorm(100))
  d$y <- rnbinom(100, size = 2, mu = exp(1 + 0.5 * d$x))
  expect_output(nbreg(y ~ x, d, control = nbreg_control(trace = TRUE)), "logLik")
})

test_that("nbreg rejects data with no finite MLE and non-finite input", {
  expect_error(nbreg(y ~ 1, data.frame(y = rep(0, 10))), "all responses are zero")
  expect_error(nbreg(y ~ 1, data.frame(y = c(1, Inf, 2))), "non-finite")
  expect_error(nbreg_control(tol = Inf), "finite")
  expect_error(nbreg_control(maxit = NA_real_), "finite")
  # a factor level with only zero counts: the coefficient goes to -Inf
  set.seed(16)
  d <- data.frame(g = gl(2, 20))
  d$y <- c(rep(0, 20), rnbinom(20, size = 2, mu = 3))
  expect_error(nbreg(y ~ g, d), "no finite")
})

test_that("predict() on training rows reproduces fitted values with poly()", {
  set.seed(17)
  d <- data.frame(x = runif(300, 0, 3))
  d$y <- rnbinom(300, size = 2, mu = exp(0.5 + 0.8 * d$x - 0.2 * d$x^2))
  fit <- nbreg(y ~ poly(x, 2), d)
  expect_equal(predict(fit, newdata = d[1:10, , drop = FALSE]), fitted(fit)[1:10],
               ignore_attr = TRUE)
})
