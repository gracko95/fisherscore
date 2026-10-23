# Benchmark nbreg() / zinbreg() against MASS::glm.nb() / pscl::zeroinfl()
# on pscl::bioChemists and AER::NMES1988 (Section 5 of the project description).
#
# Run from the package root:  Rscript inst/benchmark/benchmark.R
# Saves inst/benchmark/benchmark_results.rds, which the report reads.

library(nbzinbKoo)
library(MASS)
library(pscl)
data("bioChemists", package = "pscl")
data("NMES1988", package = "AER")

# Median elapsed time of `reps` evaluations of an expression
time_it <- function(expr, reps = 5) {
  e <- substitute(expr)
  env <- parent.frame()
  stats::median(replicate(reps, system.time(eval(e, env))[["elapsed"]]))
}

# Side-by-side estimates and SEs, plus the max absolute differences of Table 4
compare_nb <- function(formula, data, label) {
  mine <- nbreg(formula, data)
  ref <- glm.nb(formula, data = data, control = glm.control(epsilon = 1e-10, maxit = 100))
  se_m <- sqrt(diag(vcov(mine)))
  se_r <- sqrt(diag(vcov(ref)))
  table <- data.frame(term = c(names(coef(mine)), "theta"),
                      nbreg = c(coef(mine), mine$theta),
                      glm.nb = c(coef(ref), ref$theta),
                      se_nbreg = c(se_m, mine$SE.theta),
                      se_glm.nb = c(se_r, ref$SE.theta), row.names = NULL)
  diffs <- data.frame(data = label, model = "NB",
                      coef = max(abs(coef(mine) - coef(ref))),
                      se = max(abs(se_m - se_r)),
                      theta = abs(mine$theta - ref$theta),
                      loglik = as.numeric(logLik(mine)) - as.numeric(logLik(ref)),
                      iterations = mine$iter,
                      time_mine = time_it(nbreg(formula, data)),
                      time_ref = time_it(glm.nb(formula, data = data,
                                                control = glm.control(epsilon = 1e-10, maxit = 100))))
  list(table = table, diffs = diffs)
}

compare_zinb <- function(formula, data, label, tol = 1e-8) {
  ctrl <- zinbreg_control(tol = tol)
  mine <- zinbreg(formula, data, control = ctrl)
  ref <- zeroinfl(formula, data = data, dist = "negbin",
                  control = zeroinfl.control(reltol = 1e-12))
  se_m <- sqrt(diag(vcov(mine)))
  se_r <- sqrt(diag(vcov(ref)))
  table <- data.frame(term = c(names(coef(mine)), "log(theta)"),
                      zinbreg = c(coef(mine), log(mine$theta)),
                      zeroinfl = c(coef(ref), log(ref$theta)),
                      se_zinbreg = c(se_m, mine$SE.logtheta),
                      se_zeroinfl = c(se_r, ref$SE.logtheta), row.names = NULL)
  diffs <- data.frame(data = label, model = sprintf("ZINB (tol %g)", tol),
                      coef = max(abs(coef(mine) - coef(ref))),
                      se = max(abs(se_m - se_r)),
                      theta = abs(log(mine$theta) - log(ref$theta)),
                      loglik = as.numeric(logLik(mine)) - as.numeric(logLik(ref)),
                      iterations = mine$iter,
                      time_mine = time_it(zinbreg(formula, data, control = ctrl), reps = 3),
                      time_ref = time_it(zeroinfl(formula, data = data, dist = "negbin",
                                                  control = zeroinfl.control(reltol = 1e-12)),
                                         reps = 3))
  list(table = table, diffs = diffs)
}

f_bio <- art ~ fem + mar + kid5 + phd + ment
f_bio_zinb <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
f_nmes <- visits ~ hospital + health + chronic + gender + school + insurance
f_nmes_zinb <- visits ~ hospital + health + chronic + gender + school + insurance |
  hospital + chronic + gender + school + insurance

results <- list(
  nb_bio = compare_nb(f_bio, bioChemists, "bioChemists"),
  nb_nmes = compare_nb(f_nmes, NMES1988, "NMES1988"),
  zinb_bio = compare_zinb(f_bio_zinb, bioChemists, "bioChemists"),
  zinb_nmes = compare_zinb(f_nmes_zinb, NMES1988, "NMES1988"),
  zinb_bio_tight = compare_zinb(f_bio_zinb, bioChemists, "bioChemists", tol = 1e-10),
  zinb_nmes_tight = compare_zinb(f_nmes_zinb, NMES1988, "NMES1988", tol = 1e-10)
)
summary_table <- do.call(rbind, lapply(results, `[[`, "diffs"))
rownames(summary_table) <- NULL
print(summary_table, digits = 3)

saveRDS(list(results = results, summary = summary_table, session = utils::sessionInfo()),
        file.path("inst", "benchmark", "benchmark_results.rds"))
