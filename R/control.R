#' Control parameters for Fisher scoring
#'
#' The defaults are the same as [stats::glm.control()].
#'
#' @param epsilon Positive convergence tolerance. Iterations stop when
#'   \eqn{|D^{(t)} - D^{(t-1)}| / (|D^{(t)}| + 0.1) < \epsilon}, where \eqn{D}
#'   is the deviance.
#' @param maxit Maximum number of Fisher scoring iterations.
#' @param trace Logical; print the deviance at each iteration?
#' @return A list with components `epsilon`, `maxit` and `trace`.
#' @examples
#' fs_control(epsilon = 1e-10, maxit = 50)
#' @export
fs_control <- function(epsilon = 1e-8, maxit = 25, trace = FALSE) {
  if (!is.numeric(epsilon) || epsilon <= 0) stop("'epsilon' must be > 0", call. = FALSE)
  if (!is.numeric(maxit) || maxit <= 0) stop("'maxit' must be > 0", call. = FALSE)
  list(epsilon = epsilon, maxit = as.integer(maxit), trace = isTRUE(trace))
}
