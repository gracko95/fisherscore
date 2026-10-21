d <- sim_data()
g <- glm(yb ~ x1 + g, binomial("probit"), d)
f <- fs_glm(yb ~ x1 + g, binomial("probit"), d, engine = "cpp")

test_that("summary coefficient table matches summary.glm", {
  expect_equal(summary(f)$coefficients, summary(g)$coefficients, tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_output(print(summary(f)), "Fisher Scoring iterations")
  expect_output(print(f), "Residual Deviance")
})

test_that("residuals of every type match glm", {
  for (type in c("deviance", "pearson", "working", "response"))
    expect_equal(unname(residuals(f, type)), unname(residuals(g, type)), tolerance = 1e-8)
})

test_that("predict matches predict.glm, including newdata and se.fit", {
  nd <- data.frame(x1 = c(-1, 0, 2), g = factor(c("a", "b", "c")))
  expect_equal(unname(predict(f)), unname(predict(g)), tolerance = 1e-8)
  expect_equal(unname(predict(f, nd, type = "response")),
               unname(predict(g, nd, type = "response")), tolerance = 1e-8)
  pf <- predict(f, nd, type = "response", se.fit = TRUE)
  pg <- predict(g, nd, type = "response", se.fit = TRUE)
  expect_equal(unname(pf$se.fit), unname(pg$se.fit), tolerance = 1e-8)
})

test_that("logLik, nobs and Wald confint match glm", {
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-8)
  expect_equal(attr(logLik(f), "df"), attr(logLik(g), "df"))
  expect_equal(nobs(f), nobs(g))
  expect_equal(confint(f), suppressMessages(confint.default(g)), tolerance = 1e-8)
})

test_that("anova LRT matches anova.glm", {
  f0 <- fs_glm(yb ~ x1, binomial("probit"), d)
  g0 <- glm(yb ~ x1, binomial("probit"), d)
  a <- anova(f0, f)
  b <- anova(g0, g, test = "LRT")
  expect_equal(a$Deviance, b$Deviance, tolerance = 1e-8)
  expect_equal(a[["Pr(>Chi)"]], b[["Pr(>Chi)"]], tolerance = 1e-8)
})

test_that("compare_glm reports agreement with glm", {
  cmp <- compare_glm(yb ~ x1 + g, binomial("probit"), d, times = 1)
  expect_s3_class(cmp, "fs_comparison")
  expect_true(all(cmp$summary$max_coef_diff < 1e-8))
  expect_output(print(cmp), "Fit summary")
})

test_that("plot_convergence returns the trace", {
  pdf(NULL); on.exit(dev.off())
  tr <- plot_convergence(f, "change")
  expect_equal(nrow(tr), f$iter)
})

test_that("predict() with newdata re-evaluates an offset argument", {
  d <- sim_data()
  g <- glm(yp ~ x1, poisson, d, offset = log(exposure))
  f <- fs_glm(yp ~ x1, poisson, d, offset = log(exposure))
  expect_equal(predict(f, d), predict(f))
  nd <- d[1:5, ]
  nd$exposure <- nd$exposure * 2
  expect_equal(unname(predict(f, nd, type = "response")),
               unname(predict(g, nd, type = "response")), tolerance = 1e-8)
})

test_that("a supplied dispersion gives z tests, as in summary.glm", {
  d <- sim_data()[1:10, ]
  g <- summary(glm(yg ~ x1, gaussian, d), dispersion = 1)$coefficients
  f <- summary(fs_glm(yg ~ x1, gaussian, d), dispersion = 1)$coefficients
  expect_equal(colnames(f), colnames(g))
  expect_equal(f, g, tolerance = 1e-8)
})

test_that("na.exclude pads fitted values, residuals and predictions with NA", {
  d <- data.frame(x = c(1, 2, NA, 4, 5, 6), y = c(1, 3, 2, 5, 4, 8))
  g <- glm(y ~ x, poisson, d, na.action = na.exclude)
  f <- fs_glm(y ~ x, poisson, d, na.action = na.exclude)
  expect_equal(fitted(f), fitted(g), tolerance = 1e-8)
  expect_equal(residuals(f), residuals(g), tolerance = 1e-8)
  expect_equal(residuals(f, "pearson"), residuals(g, "pearson"), tolerance = 1e-8)
  expect_equal(predict(f), predict(g), tolerance = 1e-8)
  expect_equal(predict(f, se.fit = TRUE)$se.fit, predict(g, se.fit = TRUE)$se.fit,
               tolerance = 1e-8)
  expect_output(print(summary(f)))
})

test_that("nobs() ignores grouped binomial rows with zero trials", {
  d <- data.frame(s = c(2, 0, 3, 1), f = c(3, 0, 2, 4), x = 1:4)
  expect_equal(nobs(fs_glm(cbind(s, f) ~ x, binomial, d)),
               nobs(glm(cbind(s, f) ~ x, binomial, d)))
})

test_that("predict() with incomplete newdata returns NA rows and keeps offsets aligned", {
  d <- sim_data()
  nd <- data.frame(x1 = c(1, NA, 3), exposure = c(1, 2, 4))
  for (off in c(TRUE, FALSE)) {
    g <- if (off) glm(yp ~ x1, poisson, d, offset = log(exposure)) else glm(yp ~ x1, poisson, d)
    f <- if (off) fs_glm(yp ~ x1, poisson, d, offset = log(exposure)) else fs_glm(yp ~ x1, poisson, d)
    expect_equal(predict(f, nd), predict(g, nd), tolerance = 1e-8)
    expect_equal(predict(f, nd, se.fit = TRUE)$se.fit, predict(g, nd, se.fit = TRUE)$se.fit,
                 tolerance = 1e-8)
  }
})
