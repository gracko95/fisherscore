# Internal helpers shared by nbreg() and zinbreg(). Equation numbers refer to
# Sections 2-3 of the project description (project-docs/STA651_Final_Project.pdf).

# ---- Negative binomial (NB2) building blocks ---------------------------------

# Weighted NB2 log-likelihood, Section 2.1:
#   l(beta, theta) = sum_i a_i [ log G(y_i + theta) - log G(theta) - log y_i!
#                     + theta log theta + y_i log mu_i - (theta + y_i) log(theta + mu_i) ]
# dnbinom(size = theta, mu = mu) is exactly this density. theta = Inf gives
# the Poisson log-likelihood (used for the Poisson starting fit).
nb_loglik <- function(y, mu, theta, a) {
  if (is.infinite(theta)) return(sum(a * stats::dpois(y, mu, log = TRUE)))
  sum(a * stats::dnbinom(y, size = theta, mu = mu, log = TRUE))
}

# Score for theta, Section 3.1(b):
#   U(theta) = sum_i a_i [ psi(y_i + theta) - psi(theta) + log theta + 1
#                          - log(theta + mu_i) - (y_i + theta) / (theta + mu_i) ]
nb_score_theta <- function(y, mu, theta, a) {
  sum(a * (digamma(y + theta) - digamma(theta) + log(theta) + 1 -
             log(theta + mu) - (y + theta) / (theta + mu)))
}

# Second derivative for theta, Section 3.1(b):
#   H(theta) = sum_i a_i [ psi'(y_i + theta) - psi'(theta) + 1/theta
#                          - 2 / (theta + mu_i) + (y_i + theta) / (theta + mu_i)^2 ]
nb_hess_theta <- function(y, mu, theta, a) {
  sum(a * (trigamma(y + theta) - trigamma(theta) + 1 / theta -
             2 / (theta + mu) + (y + theta) / (theta + mu)^2))
}

# ---- Weighted least squares ---------------------------------------------------

# One weighted least squares step: beta = (X'WX)^{-1} X'Wz.
# Solved through the QR decomposition of W^{1/2} X, which is numerically more
# stable than forming and inverting X'WX. Rows with w_i = 0 contribute nothing.
wls <- function(X, z, w) {
  sw <- sqrt(w)
  qx <- qr(X * sw)
  if (qx$rank < ncol(X))
    stop("X'WX is singular: the model matrix is rank deficient for the current weights",
         call. = FALSE)
  drop(qr.coef(qx, z * sw))
}

# (X'WX)^{-1}, computed from the R factor of the QR decomposition of
# W^{1/2} X: X'WX = R'R, so (X'WX)^{-1} = R^{-1} R^{-T} (chol2inv).
xwx_inverse <- function(X, w) {
  qx <- qr(X * sqrt(w))
  V <- chol2inv(qr.R(qx))
  V[qx$pivot, qx$pivot] <- V          # undo any column pivoting
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}

# ---- Input checks -------------------------------------------------------------

check_counts <- function(y) {
  if (!is.numeric(y) || !is.null(dim(y)))
    stop("the response must be a numeric vector of counts", call. = FALSE)
  if (any(y < 0) || any(abs(y - round(y)) > 1e-8))
    stop("the response must contain non-negative integer counts", call. = FALSE)
  invisible(y)
}

check_full_rank <- function(X, what = "model matrix") {
  qx <- qr(X)
  if (qx$rank < ncol(X)) {
    bad <- colnames(X)[qx$pivot[(qx$rank + 1):ncol(X)]]
    stop(sprintf("the %s is rank deficient; aliased column(s): %s", what,
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
  invisible(X)
}

# Coefficient table for summary(): estimate, SE, z = est / SE, two-sided
# p-value from the standard normal (Wald test).
coef_table <- function(est, se) {
  z <- est / se
  cbind(Estimate = est, `Std. Error` = se, `z value` = z,
        `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
}
