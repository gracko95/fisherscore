test_that("binomial links match glm on simulated 0/1 data", {
  d <- sim_data()
  for (link in c("logit", "probit", "cloglog"))
    expect_matches_glm(yb ~ x1 + x2 + g, binomial(link), d)
})

test_that("grouped binomial data in cbind(successes, failures) form matches glm", {
  expect_matches_glm(cbind(ncases, ncontrols) ~ agegp + alcgp + tobgp, binomial, esoph)
  expect_matches_glm(cbind(ncases, ncontrols) ~ agegp + alcgp, binomial("cloglog"), esoph)
})

test_that("proportions with weights and factor responses match glm", {
  skip_if_not_installed("MASS")
  men <- MASS::menarche
  men$p <- men$Menarche / men$Total
  expect_matches_glm(p ~ Age, binomial("probit"), men, weights = Total)
  d <- sim_data()
  d$yf <- factor(ifelse(d$yb == 1, "yes", "no"))
  expect_matches_glm(yf ~ x1, binomial, d)
})

test_that("bad binomial responses are rejected", {
  d <- sim_data()
  expect_error(fs_glm(yp ~ x1, binomial, d), "0 <= y <= 1")
})
