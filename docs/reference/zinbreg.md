# Zero-inflated negative binomial regression by EM with IRLS M-steps

Fits the ZINB model \$\$P(Y_i = 0) = \pi_i + (1-\pi_i)
f\_{NB}(0;\mu_i,\theta), \qquad P(Y_i = y) = (1-\pi_i)
f\_{NB}(y;\mu_i,\theta),\\ y \ge 1,\$\$ with \\\log\mu_i =
x_i^\top\beta + o_i\\ (count model) and \\\mathrm{logit}\\\pi_i =
z_i^\top\gamma\\ (zero-inflation model), by an EM algorithm:

- **E-step**: \\\tau_i = \pi_i / \\\pi_i + (1-\pi_i)
  f\_{NB}(0;\mu_i,\theta)\\\\ if \\y_i = 0\\, else \\\tau_i = 0\\.

- **M-step, zero part**: logistic IRLS on the fractional response
  \\\tau_i\\.

- **M-step, count part**: the NB fit of
  [`nbreg()`](https://gracko95.github.io/nbzinbKoo/reference/nbreg.md)
  with prior weights \\a_i = 1 - \tau_i\\.

## Usage

``` r
zinbreg(formula, data, offset, control = zinbreg_control())
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model covariates
  before `|`, zero-model covariates after it. Without `|` the same
  covariates are used in both parts, as in
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html).
  Offsets ([`offset()`](https://rdrr.io/r/stats/offset.html) terms) are
  allowed in the count part.

- data:

  A data frame containing the variables in `formula`.

- offset:

  Optional offset for the count model, evaluated in `data`.

- control:

  A list from
  [`zinbreg_control()`](https://gracko95.github.io/nbzinbKoo/reference/zinbreg_control.md).

## Value

An object of class `"zinbreg"`, a list with components `coefficients` (a
list with `count` and `zero`), `vcov` (for \\(\beta,\gamma)\\, named
`count_*` and `zero_*` as in `zeroinfl()`), `vcov_full` (also including
`log(theta)`), `theta`, `SE.theta`, `SE.logtheta`, `loglik`, `iter`,
`converged`, `monotone` (did the log-likelihood never decrease?),
`fitted.values` (\\(1-\hat\pi)\hat\mu\\), `mu`, `pi`, `tau`, `residuals`
(response), `score` (gradient at the estimate), `trace`, `y`, `offset`,
`nobs`, `df.residual`, `call`, `formula`, `terms`, `xlevels`,
`contrasts` and `control`.

## Details

EM stops when the observed-data log-likelihood changes by less than
`control$tol`; a warning is given if it ever decreases. Standard errors
come from inverting the observed information of the observed-data
log-likelihood in \\(\beta, \gamma, \log\theta)\\.

## See also

[`nbreg()`](https://gracko95.github.io/nbzinbKoo/reference/nbreg.md),
[zinbreg-methods](https://gracko95.github.io/nbzinbKoo/reference/zinbreg-methods.md)

## Examples

``` r
# Simulated data: count part beta = (1, 0.5), zero part gamma = (-1, 1), theta = 2
set.seed(1)
n <- 500
d <- data.frame(x = rnorm(n), z = rnorm(n))
structural_zero <- rbinom(n, 1, plogis(-1 + d$z))
d$y <- ifelse(structural_zero == 1, 0, rnbinom(n, size = 2, mu = exp(1 + 0.5 * d$x)))
fit <- zinbreg(y ~ x | z, data = d)
summary(fit)
#> 
#> Call:
#> zinbreg(formula = y ~ x | z, data = d)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.0733 -0.7077 -0.3795  0.5520  5.9729 
#> 
#> Count model coefficients (negbin with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.10534    0.06122  18.056  < 2e-16 ***
#> x            0.42999    0.05244   8.200 2.41e-16 ***
#> Log(theta)   0.67132    0.19094   3.516 0.000438 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -1.0264     0.1997  -5.139 2.76e-07 ***
#> z             0.9521     0.1653   5.761 8.35e-09 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Theta = 1.95683 (SE 0.3736)
#> Log-likelihood: -949.2117 on 5 Df,  AIC: 1908.423
#> EM iterations: 34; log-likelihood never decreased: TRUE
#> 
head(predict(fit, type = "zero"))
#>         1         2         3         4         5         6 
#> 0.2783249 0.2126466 0.1040524 0.2658783 0.4794476 0.6204066 
```
