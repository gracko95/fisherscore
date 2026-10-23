# Simulation study for zinbreg() (Section 5.1 of the project description).
#
# R = 200 data sets of size n = 500 from a ZINB model with two count-model
# covariates and one zero-model covariate, with known parameters. For each
# parameter we report the bias, the mean estimated SE versus the Monte Carlo
# SD of the estimates, and the coverage of 95% Wald intervals.
#
# Run from the package root:  Rscript inst/benchmark/simulation.R
# which also saves inst/benchmark/simulation_results.rds. The report runs this
# same script (installed as benchmark/simulation.R) in a cached chunk that is
# re-run whenever the package code changes.

library(nbzinbKoo)

set.seed(651)
R <- 200
n <- 500
beta <- c(`(Intercept)` = 1, x1 = 0.5, x2 = -0.7)   # count model: log mu = 1 + 0.5 x1 - 0.7 x2
gamma <- c(`(Intercept)` = -1, z1 = 0.8)             # zero model: logit pi = -1 + 0.8 z1
theta <- 1.5
truth <- c(setNames(beta, paste0("count_", names(beta))),
           setNames(gamma, paste0("zero_", names(gamma))),
           `log(theta)` = log(theta))

simulate_one <- function() {
  d <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.5), z1 = rnorm(n))
  mu <- exp(beta[1] + beta[2] * d$x1 + beta[3] * d$x2)
  pi <- plogis(gamma[1] + gamma[2] * d$z1)
  structural_zero <- rbinom(n, 1, pi)
  d$y <- ifelse(structural_zero == 1, 0, rnbinom(n, size = theta, mu = mu))
  d
}

fit_one <- function(d) {
  fit <- zinbreg(y ~ x1 + x2 | z1, data = d)
  est <- c(coef(fit), `log(theta)` = log(fit$theta))
  se <- c(sqrt(diag(vcov(fit))), fit$SE.logtheta)
  c(est = est, se = se, converged = fit$converged, iter = fit$iter)
}

elapsed <- system.time(
  sims <- t(replicate(R, fit_one(simulate_one())))
)[["elapsed"]]

k <- length(truth)
est <- sims[, seq_len(k)]
se <- sims[, k + seq_len(k)]
colnames(est) <- colnames(se) <- names(truth)
lower <- est - qnorm(0.975) * se
upper <- est + qnorm(0.975) * se
covered <- sweep(lower, 2, truth, "<=") & sweep(upper, 2, truth, ">=")

summary_table <- data.frame(
  parameter = names(truth),
  true = unname(truth),
  mean_estimate = colMeans(est),
  bias = colMeans(est) - truth,
  mean_se = colMeans(se),
  mc_sd = apply(est, 2, sd),
  coverage_95 = colMeans(covered),
  row.names = NULL
)
print(summary_table, digits = 3)
cat(sprintf("\n%d of %d fits converged; median EM iterations %g; total time %.1f s\n",
            sum(sims[, "converged"]), R, median(sims[, "iter"]), elapsed))

results <- list(summary = summary_table, estimates = est, se = se,
                converged = sims[, "converged"], iterations = sims[, "iter"],
                settings = list(R = R, n = n, beta = beta, gamma = gamma, theta = theta,
                                seed = 651),
                elapsed = elapsed, session = utils::sessionInfo())
# Save only when run from the package root (not when the report sources it)
if (dir.exists(file.path("inst", "benchmark")))
  saveRDS(results, file.path("inst", "benchmark", "simulation_results.rds"))
