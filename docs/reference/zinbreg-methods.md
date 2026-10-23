# Methods for zinbreg fits

Standard S3 methods for objects returned by
[`zinbreg()`](https://gracko95.github.io/nbzinbKoo/reference/zinbreg.md).

## Usage

``` r
# S3 method for class 'zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
summary(object, ...)

# S3 method for class 'summary.zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
coef(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
vcov(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
logLik(object, ...)

# S3 method for class 'zinbreg'
nobs(object, ...)

# S3 method for class 'zinbreg'
fitted(object, ...)

# S3 method for class 'zinbreg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'zinbreg'
predict(
  object,
  newdata = NULL,
  type = c("response", "count", "zero", "prob0"),
  ...
)
```

## Arguments

- digits:

  Number of significant digits to print.

- ...:

  Not used.

- object, x:

  A `"zinbreg"` object (or its summary for `print`).

- model:

  For [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html): `"full"` (both parts,
  named `count_*` and `zero_*` as in
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html)),
  `"count"` or `"zero"`.

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"response"` (\\y - (1-\hat\pi)\hat\mu\\) or `"pearson"` (divided by
  \\\sqrt{\mathrm{Var}(Y)}\\ with \\\mathrm{Var}(Y) = (1-\pi)\mu\\1 +
  \mu(\pi + 1/\theta)\\\\). For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`
  (\\(1-\hat\pi)\hat\mu\\), `"count"` (\\\hat\mu\\), `"zero"`
  (\\\hat\pi\\) or `"prob0"` (\\P(Y=0) = \hat\pi +
  (1-\hat\pi)f\_{NB}(0;\hat\mu,\hat\theta)\\).

- newdata:

  Optional data frame to predict for.

## Value

As for
[nbreg-methods](https://gracko95.github.io/nbzinbKoo/reference/nbreg-methods.md).
[`logLik()`](https://rdrr.io/r/stats/logLik.html) has `df` = number of
count and zero coefficients + 1.

## Examples

``` r
set.seed(3)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- ifelse(rbinom(n, 1, plogis(-0.5 + d$z)) == 1, 0,
              rnbinom(n, size = 2, mu = exp(1 + 0.5 * d$x)))
fit <- zinbreg(y ~ x | z, data = d)
coef(fit)
#> count_(Intercept)           count_x  zero_(Intercept)            zero_z 
#>         1.0112410         0.4266714        -0.5951456         0.8375193 
coef(fit, model = "zero")
#> (Intercept)           z 
#>  -0.5951456   0.8375193 
AIC(fit)
#> [1] 1396.044
nd <- data.frame(x = 0, z = c(-1, 0, 1))
predict(fit, newdata = nd, type = "response")
#>        1        2        3 
#> 2.219321 1.771861 1.208744 
predict(fit, newdata = nd, type = "zero")
#>         1         2         3 
#> 0.1926838 0.3554551 0.5602985 
```
