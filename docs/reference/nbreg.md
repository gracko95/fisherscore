# Negative binomial regression by IRLS and Newton-Raphson

Fits the NB2 regression model \$\$\log\mu_i = x_i^\top\beta + o_i,
\qquad \mathrm{Var}(Y_i) = \mu_i + \mu_i^2/\theta\$\$ by maximum
likelihood, alternating

- IRLS for \\\beta\\ with \\\theta\\ fixed: working weights \\w_i =
  \theta\mu_i/(\theta+\mu_i)\\, working response \\z_i = (\eta_i -
  o_i) + (y_i-\mu_i)/\mu_i\\, \\\beta \leftarrow (X^\top WX)^{-1}X^\top
  Wz\\; and

- Newton-Raphson for \\\phi = \log\theta\\ with \\\beta\\ fixed,

## Usage

``` r
nbreg(formula, data, offset, control = nbreg_control())
```

## Arguments

- formula:

  A model formula, e.g. `y ~ x1 + x2`. Offsets may be given with
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms.

- data:

  A data frame containing the variables in `formula`.

- offset:

  Optional offset for the linear predictor (e.g. log exposure),
  evaluated in `data`.

- control:

  A list from
  [`nbreg_control()`](https://gracko95.github.io/nbzinbKoo/reference/nbreg_control.md).

## Value

An object of class `"nbreg"`, a list with components `coefficients`,
`vcov`, `theta`, `SE.theta`, `loglik`, `iter`, `converged`,
`fitted.values` (\\\hat\mu\\), `linear.predictors`, `residuals`
(response), `y`, `offset`, `df.residual`, `nobs`, `trace`
(log-likelihood and \\\theta\\ by iteration), `call`, `terms`,
`xlevels`, `contrasts` and `control`.

## Details

until the log-likelihood changes by less than `control$tol`. No general
purpose optimizer is used. The parameterization is the one used by
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html).

Standard errors for \\\beta\\ are from \\(X^\top\hat WX)^{-1}\\,
treating \\\theta\\ as known (as `glm.nb()` does), and
\\\mathrm{SE}(\hat\theta) = \\-H(\hat\theta)\\^{-1/2}\\.

## See also

[`zinbreg()`](https://gracko95.github.io/nbzinbKoo/reference/zinbreg.md),
[nbreg-methods](https://gracko95.github.io/nbzinbKoo/reference/nbreg-methods.md)

## Examples

``` r
# Simulated data with known beta = (1, 0.5) and theta = 2
set.seed(1)
d <- data.frame(x = rnorm(300))
d$y <- rnbinom(300, size = 2, mu = exp(1 + 0.5 * d$x))
fit <- nbreg(y ~ x, data = d)
summary(fit)
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.1750 -0.7653 -0.2194  0.3855  5.6578 
#> 
#> Coefficients (negative binomial with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.91318    0.05784  15.789  < 2e-16 ***
#> x            0.43526    0.05997   7.257 3.95e-13 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Theta = 1.78037 (SE 0.2737)
#> Log-likelihood: -626.523 on 3 Df,  AIC: 1259.046
#> Outer iterations (IRLS + Newton-Raphson): 3
#> 

# With an exposure offset
d$t <- runif(300, 1, 3)
d$y2 <- rnbinom(300, size = 2, mu = d$t * exp(1 + 0.5 * d$x))
coef(nbreg(y2 ~ x, data = d, offset = log(t)))
#> (Intercept)           x 
#>   1.0929316   0.3934869 
```
