# Control parameters for zinbreg()

Control parameters for zinbreg()

## Usage

``` r
zinbreg_control(
  tol = 1e-08,
  maxit = 1000,
  inner_tol = 1e-10,
  inner_maxit = 100,
  theta_max = 1e+05,
  hessian = c("analytic", "numeric"),
  trace = FALSE
)
```

## Arguments

- tol:

  EM convergence tolerance: stop when the observed-data log-likelihood
  changes by less than `tol` between two EM iterations.

- maxit:

  Maximum number of EM iterations.

- inner_tol:

  Tolerance for the IRLS loops inside each M-step.

- inner_maxit:

  Maximum number of iterations of each inner loop.

- theta_max:

  Upper limit for \\\theta\\ (see
  [`nbreg_control()`](https://gracko95.github.io/nbzinbKoo/reference/nbreg_control.md)).

- hessian:

  How to compute the observed information for the standard errors:
  `"analytic"` (exact second derivatives, the default) or `"numeric"`
  ([`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
  of the package's own log-likelihood; needs the numDeriv package).

- trace:

  Logical; print the observed-data log-likelihood at each EM iteration?

## Value

A list of control parameters.

## Examples

``` r
zinbreg_control(tol = 1e-10, maxit = 2000)
#> $tol
#> [1] 1e-10
#> 
#> $maxit
#> [1] 2000
#> 
#> $inner_tol
#> [1] 1e-10
#> 
#> $inner_maxit
#> [1] 100
#> 
#> $theta_max
#> [1] 1e+05
#> 
#> $trace
#> [1] FALSE
#> 
#> $hessian
#> [1] "analytic"
#> 
```
