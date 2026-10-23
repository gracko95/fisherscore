# STA 651 final project: NB and ZINB regression R package (nbzinbKoo)

Individual project for STA 651 (Fall 2026, Dr. Donghyung Lee). Full spec:
project-docs/STA651_Final_Project.pdf. Setup steps: project-docs/STA651_Setup_Guide.pdf.

## Rules (must follow)
- Fitting code in R/ must NOT call optim, optimize, nlm, nlminb, uniroot,
  glm, glm.fit, lm, lm.fit, lm.wfit, MASS::glm.nb, MASS::theta.ml,
  pscl::zeroinfl, or any other optimizer package (optimx, maxLik,
  stats4::mle, bbmle, countreg, glmmTMB, gamlss, VGAM).
- Allowed in fitting code: solve, qr, chol, crossprod, %*%, lgamma,
  digamma, trigamma, dnbinom, plogis, qlogis, model.frame, model.matrix,
  model.response, model.offset, the Formula package (two-part formulas),
  and numDeriv::hessian for standard errors only.
- NB: IRLS for beta alternated with Newton-Raphson on log(theta).
- ZINB: EM algorithm; M-steps are IRLS (weighted logistic + weighted NB).
- glm, MASS::glm.nb and pscl::zeroinfl may be used ONLY in tests/,
  inst/benchmark/, and the report, for comparison.
- Follow the equations in project-docs/STA651_Final_Project.pdf, Section 3.
- Comment estimation code with the equation it implements.
- Explain your reasoning to me; I must be able to explain every line
  (graded oral walkthrough, Nov 2-6).

## Workflow
- After editing roxygen comments: Rscript -e 'devtools::document()'
- Run tests: Rscript -e 'devtools::test()'
- Full check (must be 0 errors, 0 warnings): Rscript -e 'devtools::check()'
- Report PDF: quarto render vignettes/articles/report.qmd --to pdf
- Website: Rscript -e 'pkgdown::build_site()'
- Make small git commits with clear messages after each working step.
- Before submitting, check R/ for prohibited calls:
  grep -rnE "optim|nlm|uniroot|glm|lm\\(|zeroinfl" R/
  Every match must be a comment or a use I can explain.

## Algorithms (Section 3 of the spec)
NB2: log mu = X beta + offset, Var(Y) = mu + mu^2/theta, prior weights a_i.
- IRLS for beta (theta fixed): w_i = a_i theta mu_i / (theta + mu_i),
  z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i, beta = (X'WX)^{-1} X'Wz.
- Newton-Raphson for phi = log(theta) (beta fixed): U_phi = theta U,
  H_phi = theta^2 H + theta U, with U, H from Section 3.1(b) (digamma/trigamma).
- Outer loop: start from Poisson IRLS (or log ybar) and a moment estimate of
  theta; alternate until the log-likelihood changes by < 1e-8.
- SEs: Var(beta) = (X'WX)^{-1} at convergence (theta treated as known, as
  glm.nb does); SE(theta) = {-H(theta)}^{-1/2}.

ZINB: logit pi = Z gamma (zero part), NB(mu, theta) count part.
- E-step: tau_i = pi_i / (pi_i + (1 - pi_i) f_NB(0; mu_i, theta)) if y_i = 0, else 0.
- M-step zero part: IRLS logistic regression on fractional response tau,
  w_i = pi_i (1 - pi_i), z_i = z_i'gamma + (tau_i - pi_i) / {pi_i (1 - pi_i)}.
- M-step count part: the NB fit above with prior weights a_i = 1 - tau_i
  (reuse the NB code).
- Converge on the observed-data log-likelihood (< 1e-8); it must never
  decrease. Check this every iteration.
- Starts: NB fit for (beta, theta); gamma intercept-only or logistic IRLS of
  I(y = 0) on Z.
- SEs: invert the observed information of the observed-data log-likelihood in
  (beta, gamma, log theta). Analytic Hessian earns +5 bonus; numDeriv::hessian
  of my own log-likelihood is acceptable. Report log(theta) with SE, and theta.

## Required package contents (Section 4)
- Exports: nbreg(formula, data, offset, control), zinbreg(formula, data,
  offset, control), nbreg_control(), zinbreg_control() (tolerances, maxit, trace).
- zinbreg formula: y ~ x1 + x2 | z1 + z2; without |, same covariates in both parts.
- S3 classes "nbreg" and "zinbreg", containing at least: coefficients, vcov,
  theta and its SE, log-likelihood, iterations, convergence flag, fitted
  values, matched call.
- Methods for both: print, summary (+ print.summary; estimates, SEs, z values,
  p-values), coef, vcov, logLik (with df attribute so AIC/BIC work), fitted,
  residuals (response, pearson), nobs, predict (with newdata). zinbreg predict
  types: "response" ((1 - pi) mu), "count" (mu), "zero" (pi).
- Warn if theta becomes very large (Poisson limit) or maxit is hit without converging.
- roxygen2 docs with runnable @examples for every export.
- testthat: compare with MASS::glm.nb and pscl::zeroinfl (tolerances below), plus
  an offset, a single-covariate model, and informative errors for bad input.
  Use skip_if_not_installed() for reference packages.
- DESCRIPTION: imports in Imports:; MASS, pscl, AER, testthat, etc. in Suggests:.
- Must install with devtools::install_github("gracko95/nbzinbKoo") in a fresh session.
- CI: usethis::use_github_action("check-standard"), passing badge in README.
  Explain any remaining NOTEs in the README.

## Required agreement (max absolute difference)
| Quantity                 | NB vs glm.nb() | ZINB vs zeroinfl() |
|--------------------------|----------------|--------------------|
| Coefficients             | <= 1e-5        | <= 1e-3            |
| Standard errors          | <= 1e-5        | <= 1e-3            |
| theta (NB) / log theta   | <= 1e-4        | <= 1e-3            |
| Log-likelihood           | <= 1e-4        | <= 1e-3            |
Run references with tight tolerances: glm.nb(..., control = glm.control(epsilon = 1e-10,
maxit = 100)), zeroinfl(..., control = zeroinfl.control(reltol = 1e-12)). zeroinfl
labels coefficients count_* and zero_*. A higher log-likelihood than the reference is
fine; report it. Reference: on bioChemists EM converges in about 120 iterations.

## Benchmarks, simulation, data analysis (Section 5)
- Benchmark on pscl::bioChemists and AER::NMES1988: side-by-side estimate tables,
  max abs differences, and timing (system.time or bench).
- Simulation: R >= 100 data sets, n = 500, ZINB with 1-2 covariates per part and
  known (beta, gamma, theta). Report bias, mean SE vs Monte Carlo SD, and 95% Wald
  coverage. Set a seed.
- NMES1988: response visits; count covariates hospital, health, chronic, gender,
  school, insurance; zero covariates at least hospital, chronic, gender, school,
  insurance. Fit Poisson (glm, reference), NB and ZINB (my package). Include
  overdispersion/excess-zero evidence, model comparison (logLik, AIC, BIC, observed
  vs expected zeros), IRRs and ORs with 95% CIs in plain language, and which model
  I prefer and why. Appendix A is a template only: do not copy its text.

## Report (Section 6)
- Single source vignettes/articles/report.qmd (Quarto, fully reproducible, all
  numbers from code calling the installed package), 6-10 pages.
- Sections: Introduction; Methods (equations, SEs, convergence); The R package;
  Benchmarking; Simulation study; Data analysis; Discussion; AI-use statement
  (what I asked Claude Code, what I verified, one output I corrected or
  rejected); References.
- Render the PDF to report/report.pdf (committed, and uploaded to Canvas).
- pkgdown site built into docs/, served by GitHub Pages from main /docs. Repo
  must be public. Put the site link in the README and the repo About box.
- Expected layout: R/ (nbreg.R, zinbreg.R, methods.R, utils.R), man/,
  tests/testthat/ (test-nbreg.R, test-zinbreg.R), vignettes/articles/report.qmd,
  inst/benchmark/ (benchmark.R, simulation.R + results), report/report.pdf, docs/,
  .github/workflows/R-CMD-check.yaml.

## Existing code in this repo
Project code (pure R, no compiled code): R/nbreg.R (nbreg, nbreg_control,
nb_fit engine reused by the ZINB M-step), R/zinbreg.R (zinbreg, EM engine,
analytic Hessian), R/methods.R (S3 methods), R/utils.R (log-likelihood, theta
derivatives, QR weighted least squares). Benchmark and simulation scripts and
their saved .rds results are in inst/benchmark/; the report reads
simulation_results.rds via system.file(), so reinstall after re-running it.

Left over from the package's earlier life as fisherscore (not used by
nbreg/zinbreg): the Fisher scoring GLM engine fs_glm() (R/fs_glm.R,
R/fit_r.R, R/families.R, R/control.R, R/fs_glm-methods.R, src/ C++ engine),
R/compare.R, the two vignettes and their tests. R/compare.R calls stats::glm()
inside R/ and every fs_glm line matches the prohibited-function grep: remove
these files, or move compare.R to inst/benchmark/, or be ready to justify them.

## Key dates
- Submission: Monday Nov 2, 2026, 11:59 PM Eastern (PDF report + repo link on
  Canvas, Pages URL in the comment). The repo is graded as of the last commit
  before the deadline.
- Oral walkthrough: Nov 2-6, 10-15 minutes; book the slot before Nov 2.
