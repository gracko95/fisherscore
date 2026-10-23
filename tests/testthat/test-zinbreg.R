# zinbreg() against pscl::zeroinfl(dist = "negbin"), with the tolerances of
# Table 4: coefficients, SEs, log(theta) and log-likelihood all <= 1e-3.

zinb_reference <- function(formula, data, ...) {
  pscl::zeroinfl(formula, data = data, dist = "negbin", ...,
                 control = pscl::zeroinfl.control(reltol = 1e-12))
}

expect_matches_zeroinfl <- function(fit, ref) {
  expect_true(fit$converged)
  expect_true(fit$monotone)
  # zeroinfl() names its coefficients count_* and zero_*, as coef.zinbreg() does
  expect_equal(names(coef(fit)), names(coef(ref)))
  expect_lte(max(abs(coef(fit) - coef(ref))), 1e-3)
  expect_lte(max(abs(sqrt(diag(vcov(fit))) - sqrt(diag(vcov(ref))))), 1e-3)
  expect_lte(abs(log(fit$theta) - log(ref$theta)), 1e-3)
  expect_lte(abs(fit$SE.logtheta - ref$SE.logtheta), 1e-3)
  expect_lte(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), 1e-3)
}

test_that("zinbreg matches zeroinfl on bioChemists", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl", envir = environment())
  f <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  expect_matches_zeroinfl(zinbreg(f, bioChemists), zinb_reference(f, bioChemists))
})

test_that("zinbreg matches zeroinfl on NMES1988", {
  skip_on_cran()
  skip_if_not_installed("pscl")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER", envir = environment())
  f <- visits ~ hospital + health + chronic + gender + school + insurance |
    hospital + chronic + gender + school + insurance
  expect_matches_zeroinfl(zinbreg(f, NMES1988), zinb_reference(f, NMES1988))
})

sim_zinb <- function(n, seed, exposure = FALSE) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n), t = if (exposure) runif(n, 0.5, 3) else 1)
  zero <- rbinom(n, 1, plogis(-0.8 + 0.9 * d$z))
  d$y <- ifelse(zero == 1, 0, rnbinom(n, size = 1.5, mu = d$t * exp(0.8 + 0.5 * d$x)))
  d
}

test_that("a formula without | uses the same covariates in both parts", {
  skip_if_not_installed("pscl")
  d <- sim_zinb(400, 21)
  fit <- zinbreg(y ~ x + z, d)
  expect_named(coef(fit, model = "zero"), c("(Intercept)", "x", "z"))
  expect_matches_zeroinfl(fit, zinb_reference(y ~ x + z, d))
})

test_that("zinbreg handles an offset in the count model", {
  skip_if_not_installed("pscl")
  d <- sim_zinb(500, 22, exposure = TRUE)
  ref <- zinb_reference(y ~ x + offset(log(t)) | z, d)
  fit1 <- zinbreg(y ~ x | z, d, offset = log(t))
  fit2 <- zinbreg(y ~ x + offset(log(t)) | z, d)
  expect_matches_zeroinfl(fit1, ref)
  expect_matches_zeroinfl(fit2, ref)
  nd <- data.frame(x = c(0, 1), z = c(0, 0), t = c(1, 2))
  expect_equal(predict(fit1, newdata = nd, type = "count"),
               predict(fit2, newdata = nd, type = "count"), tolerance = 1e-8)
  expect_equal(predict(fit2, newdata = nd, type = "count"),
               predict(ref, newdata = nd, type = "count"), tolerance = 1e-3,
               ignore_attr = TRUE)
})

test_that("zinbreg fits single-covariate parts and an intercept-only zero part", {
  skip_if_not_installed("pscl")
  d <- sim_zinb(500, 23)
  expect_matches_zeroinfl(zinbreg(y ~ x | z, d), zinb_reference(y ~ x | z, d))
  expect_matches_zeroinfl(zinbreg(y ~ x | 1, d), zinb_reference(y ~ x | 1, d))
})

test_that("EM never decreases the observed-data log-likelihood", {
  d <- sim_zinb(500, 24)
  fit <- zinbreg(y ~ x | z, d)
  expect_true(all(diff(fit$trace$loglik) > -1e-8))
  expect_equal(fit$loglik, tail(fit$trace$loglik, 1))
})

test_that("the analytic Hessian agrees with numDeriv::hessian", {
  skip_if_not_installed("numDeriv")
  d <- sim_zinb(400, 25)
  fa <- zinbreg(y ~ x | z, d)
  fn <- zinbreg(y ~ x | z, d, control = zinbreg_control(hessian = "numeric"))
  expect_equal(fa$vcov_full, fn$vcov_full, tolerance = 1e-5)
  # and the analytic score is ~0 at the MLE
  expect_lt(max(abs(fa$score)), 1e-2)
})

test_that("zinbreg methods return the right shapes and values", {
  d <- sim_zinb(300, 26)
  fit <- zinbreg(y ~ x | z, d)
  expect_s3_class(fit, "zinbreg")
  expect_named(coef(fit), c("count_(Intercept)", "count_x", "zero_(Intercept)", "zero_z"))
  expect_named(coef(fit, model = "count"), c("(Intercept)", "x"))
  expect_equal(dim(vcov(fit)), c(4L, 4L))
  expect_equal(rownames(vcov(fit, model = "zero")), c("(Intercept)", "z"))
  expect_equal(nobs(fit), 300L)
  expect_equal(attr(logLik(fit), "df"), 5L)
  expect_equal(AIC(fit), -2 * fit$loglik + 2 * 5)
  expect_equal(BIC(fit), -2 * fit$loglik + log(300) * 5)

  mu <- predict(fit, type = "count")
  pi <- predict(fit, type = "zero")
  expect_equal(predict(fit), (1 - pi) * mu)
  expect_equal(fitted(fit), (1 - pi) * mu)
  expect_equal(fitted(fit) + residuals(fit), d$y, ignore_attr = TRUE)
  expect_equal(residuals(fit, "pearson"),
               (d$y - (1 - pi) * mu) / sqrt((1 - pi) * mu * (1 + mu * (pi + 1 / fit$theta))))
  expect_equal(predict(fit, newdata = d[1:4, ], type = "zero"), pi[1:4], ignore_attr = TRUE)
  expect_equal(predict(fit, newdata = d[1:4, ], type = "prob0"),
               pi[1:4] + (1 - pi[1:4]) * dnbinom(0, size = fit$theta, mu = mu[1:4]),
               ignore_attr = TRUE)

  s <- summary(fit)
  expect_true("Log(theta)" %in% rownames(s$coefficients$count))
  expect_output(print(s), "Zero-inflation model coefficients")
  expect_output(print(fit), "Count model coefficients")
})

test_that("zinbreg gives informative errors for bad input", {
  d <- sim_zinb(100, 27)
  expect_error(zinbreg("y ~ x | z", d), "model formula")
  expect_error(zinbreg(y ~ x | z, as.list(d)), "data frame")
  expect_error(zinbreg(I(y + 1) ~ x | z, d), "no zeros")
  expect_error(zinbreg(I(-y) ~ x | z, d), "non-negative integer")
  expect_error(zinbreg(y ~ x | z | x, d), "at most two parts")
  expect_error(zinbreg(y ~ x | z + I(2 * z), d), "zero-model matrix is rank deficient")
  expect_error(zinbreg_control(hessian = "bogus"))
})

test_that("zinbreg warns when EM hits maxit", {
  d <- sim_zinb(300, 28)
  expect_warning(fit <- zinbreg(y ~ x | z, d, control = zinbreg_control(maxit = 2)), "maxit")
  expect_false(fit$converged)
})

test_that("predict() on training subsets reproduces fitted values with poly() and factors", {
  set.seed(29)
  n <- 400
  d <- data.frame(x = runif(n, 0, 3), z = rnorm(n), f = gl(2, n / 2))
  d$y <- ifelse(rbinom(n, 1, plogis(-1 + d$z)) == 1, 0,
                rnbinom(n, size = 2, mu = exp(0.5 + 0.8 * d$x - 0.2 * d$x^2)))
  fit <- zinbreg(y ~ poly(x, 2) + f | poly(z, 2), d)
  rows <- c(1:5, 396:400)
  expect_equal(predict(fit, newdata = d[rows, ]), fitted(fit)[rows], ignore_attr = TRUE)
  expect_equal(predict(fit, newdata = d[rows, ], type = "count"), fit$mu[rows],
               ignore_attr = TRUE)
  expect_equal(predict(fit, newdata = d[rows, ], type = "zero"), fit$pi[rows],
               ignore_attr = TRUE)
})

test_that("zinbreg gives clear errors for boundary data", {
  expect_error(zinbreg(y ~ 1, data.frame(y = rep(0, 10))), "all responses are zero")
  expect_error(zinbreg(y ~ 1, data.frame(y = c(0, Inf, 2))), "non-finite")
  # zeros exactly when z < 0: complete separation in the zero model
  d <- data.frame(z = c(-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 1, 2), x = seq(-1, 1, length = 10))
  d$y <- ifelse(d$z < 0, 0, c(1, 2, 3, 1, 2, 4))
  expect_error(zinbreg(y ~ x | z, d), "separate zeros")
})
