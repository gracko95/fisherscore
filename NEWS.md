# nbzinbKoo 0.1.0

* First release: `fs_glm()` fits binomial (logit, probit, cloglog), poisson
  (log, identity, sqrt) and gaussian (identity) GLMs by Fisher scoring, with
  a pure-R engine and an RcppArmadillo engine.
* `compare_glm()` compares estimates, standard errors, deviance, AIC,
  iterations and timing with `stats::glm()`.
* S3 methods: `print`, `summary`, `coef`, `vcov`, `fitted`, `residuals`,
  `predict`, `logLik`, `nobs`, `confint` (Wald), `anova` (LRT/F).
* `plot_convergence()` shows the deviance at each iteration.
