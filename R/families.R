#' Family definitions used by the Fisher scoring engines
#'
#' `fs_family()` turns a family specification into the pieces the Fisher
#' scoring algorithm needs: the link function \eqn{g}, its inverse, the
#' derivative \eqn{d\mu/d\eta}, the variance function \eqn{V(\mu)}, unit
#' deviances and the AIC. Every piece is written out here rather than taken
#' from [stats::family()] objects, and the clamping rules follow those in
#' `stats`, so fits agree with [stats::glm()] to numerical precision.
#'
#' Supported combinations:
#' * `binomial`: `logit`, `probit`, `cloglog`
#' * `poisson`: `log`, `identity`, `sqrt`
#' * `gaussian`: `identity`
#'
#' @param family A family object (e.g. `binomial(link = "probit")`), a family
#'   function (e.g. `poisson`), or a family name such as `"binomial"`.
#' @return An object of class `fs_family`.
#' @examples
#' fam <- fs_family(binomial("probit"))
#' fam$linkinv(0)
#' fam$variance(0.3)
#' @export
fs_family <- function(family = binomial()) {
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (inherits(family, "fs_family")) return(family)
  if (!inherits(family, "family")) stop("'family' not recognized", call. = FALSE)

  fam <- family$family
  link <- family$link
  supported <- list(binomial = c("logit", "probit", "cloglog"),
                    poisson  = c("log", "identity", "sqrt"),
                    gaussian = "identity")
  if (!fam %in% names(supported))
    stop(sprintf("family '%s' is not supported; use binomial, poisson or gaussian", fam),
         call. = FALSE)
  if (!link %in% supported[[fam]])
    stop(sprintf("link '%s' is not supported for the %s family (supported: %s)",
                 link, fam, paste(supported[[fam]], collapse = ", ")), call. = FALSE)

  lk <- fs_link(link)
  structure(
    c(list(family = fam, link = link,
           family_code = match(fam, names(supported)),
           link_code = match(link, .fs_link_names)),
      lk,
      fs_variance(fam)),
    class = "fs_family"
  )
}

# Link codes shared with the C++ engine (src/fisher_scoring.cpp).
.fs_link_names <- c("logit", "probit", "cloglog", "log", "identity", "sqrt")

fs_link <- function(link) {
  eps <- .Machine$double.eps
  switch(link,
    logit = list(
      linkfun = function(mu) log(mu / (1 - mu)),
      # same thresholds as the C code behind binomial()$linkinv
      linkinv = function(eta) {
        tmp <- ifelse(eta < -30, eps, ifelse(eta > 30, 1 / eps, exp(eta)))
        tmp / (1 + tmp)
      },
      mu.eta = function(eta) {
        ifelse(abs(eta) > 30, eps, exp(eta) / (1 + exp(eta))^2)
      },
      valideta = function(eta) TRUE
    ),
    probit = list(
      linkfun = function(mu) stats::qnorm(mu),
      linkinv = function(eta) {
        thresh <- -stats::qnorm(eps)
        stats::pnorm(pmin(pmax(eta, -thresh), thresh))
      },
      mu.eta = function(eta) pmax(stats::dnorm(eta), eps),
      valideta = function(eta) TRUE
    ),
    cloglog = list(
      linkfun = function(mu) log(-log(1 - mu)),
      linkinv = function(eta) pmax(pmin(-expm1(-exp(eta)), 1 - eps), eps),
      mu.eta = function(eta) {
        eta <- pmin(eta, 700)
        pmax(exp(eta) * exp(-exp(eta)), eps)
      },
      valideta = function(eta) TRUE
    ),
    log = list(
      linkfun = function(mu) log(mu),
      linkinv = function(eta) pmax(exp(eta), eps),
      mu.eta = function(eta) pmax(exp(eta), eps),
      valideta = function(eta) TRUE
    ),
    identity = list(
      linkfun = function(mu) mu,
      linkinv = function(eta) eta,
      mu.eta = function(eta) rep.int(1, length(eta)),
      valideta = function(eta) TRUE
    ),
    sqrt = list(
      linkfun = function(mu) sqrt(mu),
      linkinv = function(eta) eta^2,
      mu.eta = function(eta) 2 * eta,
      valideta = function(eta) all(is.finite(eta)) && all(eta > 0)
    )
  )
}

# y * log(y / mu), defined as 0 when y = 0
ylogy <- function(y, mu) ifelse(y > 0, y * log(y / mu), 0)

fs_variance <- function(family) {
  switch(family,
    binomial = list(
      variance = function(mu) mu * (1 - mu),
      validmu = function(mu) all(is.finite(mu)) && all(mu > 0 & mu < 1),
      dev.resids = function(y, mu, wt) 2 * wt * (ylogy(y, mu) + ylogy(1 - y, 1 - mu)),
      aic = function(y, n, mu, wt, dev) {
        m <- if (any(n > 1)) n else wt
        -2 * sum(ifelse(m > 0, wt / m, 0) *
                   stats::dbinom(round(m * y), round(m), mu, log = TRUE))
      },
      fixed_dispersion = TRUE
    ),
    poisson = list(
      variance = function(mu) mu,
      validmu = function(mu) all(is.finite(mu)) && all(mu > 0),
      dev.resids = function(y, mu, wt) 2 * wt * (ylogy(y, mu) - (y - mu)),
      aic = function(y, n, mu, wt, dev) -2 * sum(stats::dpois(y, mu, log = TRUE) * wt),
      fixed_dispersion = TRUE
    ),
    gaussian = list(
      variance = function(mu) rep.int(1, length(mu)),
      validmu = function(mu) TRUE,
      dev.resids = function(y, mu, wt) wt * (y - mu)^2,
      aic = function(y, n, mu, wt, dev) {
        nobs <- length(y)
        nobs * (log(2 * pi * dev / nobs) + 1) + 2 - sum(log(wt))
      },
      fixed_dispersion = FALSE
    )
  )
}

#' @export
print.fs_family <- function(x, ...) {
  cat("Family:", x$family, "\nLink:  ", x$link, "\n")
  invisible(x)
}
