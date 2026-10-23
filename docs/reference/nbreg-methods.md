# Methods for nbreg fits

Standard S3 methods for objects returned by
[`nbreg()`](https://gracko95.github.io/nbzinbKoo/reference/nbreg.md).

## Usage

``` r
# S3 method for class 'nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
summary(object, ...)

# S3 method for class 'summary.nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
coef(object, ...)

# S3 method for class 'nbreg'
vcov(object, ...)

# S3 method for class 'nbreg'
logLik(object, ...)

# S3 method for class 'nbreg'
nobs(object, ...)

# S3 method for class 'nbreg'
fitted(object, ...)

# S3 method for class 'nbreg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'nbreg'
predict(object, newdata = NULL, type = c("response", "link"), ...)
```

## Arguments

- digits:

  Number of significant digits to print.

- ...:

  Not used.

- object, x:

  An `"nbreg"` object (or its summary for `print`).

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"response"` (\\y - \hat\mu\\) or `"pearson"` (\\(y -
  \hat\mu)/\sqrt{\hat\mu + \hat\mu^2/\hat\theta}\\). For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`
  (\\\hat\mu\\) or `"link"` (\\\hat\eta = \log\hat\mu\\).

- newdata:

  Optional data frame to predict for. Offsets are evaluated in
  `newdata`.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) a named vector;
[`vcov()`](https://rdrr.io/r/stats/vcov.html) a matrix;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) a `"logLik"` object
with `df` = number of coefficients + 1 (for \\\theta\\), so
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) work;
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) numeric vectors;
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of
observations; [`summary()`](https://rdrr.io/r/base/summary.html) an
object of class `"summary.nbreg"`.

## Examples

``` r
set.seed(2)
d <- data.frame(x = rnorm(200))
d$y <- rnbinom(200, size = 1.5, mu = exp(0.5 + 0.8 * d$x))
fit <- nbreg(y ~ x, data = d)
fit
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Negative binomial (NB2) regression, IRLS + Newton-Raphson for theta
#> 
#> Coefficients:
#> (Intercept)            x  
#>      0.4832       0.8673  
#> 
#> Theta: 1.414 (SE 0.2772)
#> Log-likelihood: -360.955 on 3 Df
#> 
coef(fit)
#> (Intercept)           x 
#>   0.4832092   0.8673129 
sqrt(diag(vcov(fit)))
#> (Intercept)           x 
#>  0.08843103  0.08225032 
AIC(fit); BIC(fit)
#> [1] 727.9105
#> [1] 737.8054
head(residuals(fit, type = "pearson"))
#>           1           2           3           4           5           6 
#> -0.69848201  0.04582069  0.43126187 -0.65218147  0.27572733 -0.89198443 
predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
#>         1         2         3 
#> 0.6810608 1.6212690 3.8594399 
```
