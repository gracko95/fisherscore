# Control parameters for nbreg()

Control parameters for nbreg()

## Usage

``` r
nbreg_control(
  tol = 1e-08,
  maxit = 100,
  inner_tol = 1e-10,
  inner_maxit = 100,
  theta_max = 1e+05,
  trace = FALSE
)
```

## Arguments

- tol:

  Outer convergence tolerance: stop when the log-likelihood changes by
  less than `tol` between two outer iterations.

- maxit:

  Maximum number of outer iterations (one IRLS fit for \\\beta\\
  followed by one Newton-Raphson fit for \\\theta\\).

- inner_tol:

  Tolerance for the inner IRLS and Newton-Raphson loops (largest
  absolute change in \\\beta\\, or in \\\log\theta\\).

- inner_maxit:

  Maximum number of iterations of each inner loop.

- theta_max:

  Upper limit for \\\theta\\. If the estimate reaches it, the data show
  no overdispersion (the Poisson limit) and a warning is given.

- trace:

  Logical; print the log-likelihood and \\\theta\\ at each outer
  iteration?

## Value

A list of control parameters.

## Examples

``` r
nbreg_control(tol = 1e-10, trace = TRUE)
#> $tol
#> [1] 1e-10
#> 
#> $maxit
#> [1] 100
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
#> [1] TRUE
#> 
```
