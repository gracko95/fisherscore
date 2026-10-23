# nbzinbKoo

**nbzinbKoo** fits **negative binomial (NB2)** and **zero-inflated
negative binomial (ZINB)** regression models with algorithms written
from scratch, without [`optim()`](https://rdrr.io/r/stats/optim.html) or
any other general-purpose optimizer:

- [`nbreg()`](https://gracko95.github.io/nbzinbKoo/reference/nbreg.md):
  IRLS for the coefficients alternated with Newton-Raphson for
  $`\log\theta`$ (results match
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) to about
  $`10^{-8}`$);
- [`zinbreg()`](https://gracko95.github.io/nbzinbKoo/reference/zinbreg.md):
  an EM algorithm whose M-steps are a weighted logistic IRLS fit (zero
  part) and the weighted NB fit above (count part), with standard errors
  from the analytic observed information matrix (results match
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html) well
  within $`10^{-3}`$).

Website and project report: <https://gracko95.github.io/nbzinbKoo/>

STA 651 (Categorical Data Analysis) final project, Miami University,
Fall 2026.

## Installation

``` r

# install.packages("devtools")
devtools::install_github("gracko95/nbzinbKoo")
```

## Example

``` r

library(nbzinbKoo)
data("bioChemists", package = "pscl")

# Negative binomial regression
nb <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(nb)
#> 
#> Call:
#> nbreg(formula = art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.2677 -0.8754 -0.2613  0.4983  6.6573 
#> 
#> Coefficients (negative binomial with log link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.256144   0.137348   1.865 0.062191 .  
#> femWomen    -0.216418   0.072636  -2.979 0.002887 ** 
#> marMarried   0.150489   0.082097   1.833 0.066791 .  
#> kid5        -0.176415   0.052813  -3.340 0.000837 ***
#> phd          0.015271   0.035873   0.426 0.670326    
#> ment         0.029082   0.003214   9.048  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta = 2.26439 (SE 0.2712)
#> Log-likelihood: -1560.958 on 7 Df,  AIC: 3135.917
#> Outer iterations (IRLS + Newton-Raphson): 4

# Zero-inflated NB: count covariates before |, zero-inflation covariates after
zinb <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment, data = bioChemists)
summary(zinb)
#> 
#> Call:
#> zinbreg(formula = art ~ fem + mar + kid5 + phd + ment | ment, 
#>     data = bioChemists)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.3041 -0.7684 -0.2632  0.4671  6.3765 
#> 
#> Count model coefficients (negbin with log link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.404014   0.141715   2.851  0.00436 ** 
#> femWomen    -0.211906   0.071922  -2.946  0.00322 ** 
#> marMarried   0.139462   0.081192   1.718  0.08586 .  
#> kid5        -0.167624   0.052457  -3.195  0.00140 ** 
#> phd          0.001963   0.035586   0.055  0.95602    
#> ment         0.024393   0.003518   6.934 4.08e-12 ***
#> Log(theta)   1.002910   0.142830   7.022 2.19e-12 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  -0.8067     0.3532  -2.284   0.0224 *
#> ment         -0.6094     0.2458  -2.480   0.0131 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Theta = 2.7262 (SE 0.3894)
#> Log-likelihood: -1553.271 on 9 Df,  AIC: 3124.542
#> EM iterations: 112; log-likelihood never decreased: TRUE
AIC(nb, zinb)
#>      df      AIC
#> nb    7 3135.917
#> zinb  9 3124.542
head(predict(zinb, type = "zero"))   # estimated probability of a structural zero
#>            1            2            3            4            5            6 
#> 6.226167e-03 1.139290e-02 1.139290e-02 6.692070e-02 5.862284e-08 1.165473e-01
```

## Main functions

| Function | Purpose |
|----|----|
| `nbreg(formula, data, offset, control)` | NB2 regression by IRLS + Newton-Raphson for theta |
| `zinbreg(formula, data, offset, control)` | ZINB regression by EM with IRLS M-steps (`y ~ x | z`) |
| [`nbreg_control()`](https://gracko95.github.io/nbzinbKoo/reference/nbreg_control.md), [`zinbreg_control()`](https://gracko95.github.io/nbzinbKoo/reference/zinbreg_control.md) | Tolerances, maximum iterations, `trace`, Hessian type |
| `print`, `summary`, `coef`, `vcov`, `logLik`, `AIC`/`BIC`, `fitted`, `residuals`, `nobs`, `predict` | S3 methods for both classes |

## Benchmarks and simulation

- `inst/benchmark/benchmark.R`: estimates, SEs, $`\theta`$,
  log-likelihood and timing against `glm.nb()` and `zeroinfl()` on
  [`pscl::bioChemists`](https://rdrr.io/pkg/pscl/man/bioChemists.html)
  and [`AER::NMES1988`](https://rdrr.io/pkg/AER/man/NMES1988.html).
- `inst/benchmark/simulation.R`: 200 simulated ZINB data sets
  ($`n = 500`$); bias, SE accuracy and 95% Wald coverage.
- The full report is `report/report.pdf` and the [Report
  article](https://gracko95.github.io/nbzinbKoo/articles/report.html).

## R CMD check notes

[`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
gives 0 errors and 0 warnings. Any remaining NOTEs are explained here:

- TODO: list each NOTE from your final
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  run and why it is harmless (or write “None” if there are no NOTEs).
