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
