

# nbzinbKoo

<!-- badges: start -->
[![R-CMD-check](https://github.com/gracko95/nbzinbKoo/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gracko95/nbzinbKoo/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**nbzinbKoo** fits generalized linear models by **Fisher scoring**
(iteratively reweighted least squares), implemented from scratch, and lets
you check its results and speed against `stats::glm()`. It has two engines:

* `engine = "R"`: a readable pure-R implementation ([`R/fit_r.R`](R/fit_r.R))
* `engine = "cpp"`: an RcppArmadillo implementation ([`src/fisher_scoring.cpp`](src/fisher_scoring.cpp))

| Family     | Links                        |
|------------|------------------------------|
| `binomial` | `logit`, `probit`, `cloglog` |
| `poisson`  | `log`, `identity`, `sqrt`    |
| `gaussian` | `identity`                   |

Project for STA 651 (Categorical Data Analysis), Miami University.

## Installation

```r
# install.packages("remotes")
remotes::install_github("gracko95/nbzinbKoo", build_vignettes = TRUE)
```

## Example


``` r
library(nbzinbKoo)

fit <- fs_glm(cbind(ncases, ncontrols) ~ agegp + alcgp,
              family = binomial, data = esoph, engine = "cpp")
summary(fit)
#> 
#> Call:
#> fs_glm(formula = cbind(ncases, ncontrols) ~ agegp + alcgp, family = binomial, 
#>     data = esoph, engine = "cpp")
#> 
#> Deviance Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.2395 -0.7186 -0.2324  0.7930  3.3538 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -1.40326    0.20118  -6.975 3.05e-12 ***
#> agegp.L      3.68318    0.67270   5.475 4.37e-08 ***
#> agegp.Q     -1.45594    0.60521  -2.406   0.0161 *  
#> agegp.C      0.07238    0.46100   0.157   0.8752    
#> agegp^4      0.22662    0.31988   0.708   0.4787    
#> agegp^5     -0.25374    0.20994  -1.209   0.2268    
#> alcgp.L      2.59671    0.25799  10.065  < 2e-16 ***
#> alcgp.Q      0.11930    0.21915   0.544   0.5862    
#> alcgp.C      0.43863    0.17823   2.461   0.0139 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 367.95 on 87 degrees of freedom
#> Residual deviance: 105.88 on 79 degrees of freedom
#> AIC: 238.94
#> 
#> Number of Fisher Scoring iterations: 6 (cpp engine)
```

### Compare with `glm()`


``` r
compare_glm(breaks ~ wool * tension, family = poisson, data = warpbreaks)
#> Coefficients and standard errors:
#>            term   glm_est    glm_se  fs_R_est   fs_R_se fs_cpp_est fs_cpp_se
#>     (Intercept)  3.796737 0.0499375  3.796737 0.0499375   3.796737 0.0499375
#>           woolB -0.456627 0.0801920 -0.456627 0.0801920  -0.456627 0.0801920
#>        tensionM -0.618683 0.0844001 -0.618683 0.0844001  -0.618683 0.0844001
#>        tensionH -0.595799 0.0837772 -0.595799 0.0837772  -0.595799 0.0837772
#>  woolB:tensionM  0.638177 0.1221531  0.638177 0.1221531   0.638177 0.1221531
#>  woolB:tensionH  0.188363 0.1298953  0.188363 0.1298953   0.188363 0.1298953
#> 
#> Fit summary (differences are relative to stats::glm):
#>  method deviance     aic iterations converged time_sec max_coef_diff
#>     glm  182.305 468.969          4      TRUE    0.001      0.00e+00
#>    fs_R  182.305 468.969          4      TRUE    0.016      0.00e+00
#>  fs_cpp  182.305 468.969          4      TRUE    0.015      4.88e-15
#>  max_se_diff
#>     0.00e+00
#>     0.00e+00
#>     2.36e-16
```

Estimates, standard errors, deviance, AIC and iteration counts match
`glm()` to about $10^{-15}$ for every supported family and link.

## Algorithm

At each iteration, with $\eta = X\beta$, $\mu = g^{-1}(\eta)$:

$$
w_i = \frac{(\partial\mu_i/\partial\eta_i)^2}{V(\mu_i)}, \qquad
z_i = \eta_i + (y_i-\mu_i)\frac{\partial\eta_i}{\partial\mu_i}, \qquad
\beta^{(t+1)} = (X^\top W X)^{-1}X^\top W z
= \beta^{(t)} + \mathcal{I}(\beta^{(t)})^{-1}U(\beta^{(t)}).
$$

The weighted least squares step is solved by a QR decomposition, and
iteration stops when the relative change in deviance falls below
`epsilon` (default `1e-8`). See `vignette("fisher-scoring-theory")` for the
derivation and `vignette("comparison-with-glm")` for the accuracy and
speed study.

## Functions

| Function | Purpose |
|---|---|
| `fs_glm()` | Fit a GLM (formula interface, like `glm()`) |
| `fs_fit_r()`, `fs_fit_cpp()` | Low-level engines on a model matrix |
| `fs_family()` | Hand-written links, variance and deviance functions |
| `fs_control()` | Tolerance, maximum iterations, tracing |
| `compare_glm()` | Accuracy and timing compared with `glm()` |
| `plot_convergence()` | Deviance or convergence criterion by iteration |
| `summary`, `predict`, `residuals`, `anova`, `confint`, `vcov`, `logLik`, ... | S3 methods |
