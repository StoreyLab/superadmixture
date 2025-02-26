#' Sample from Dirichlet distribution
#'
#' @description
#' `rdirichlet()` returns an \eqn{n \times m} matrix of random samples from the Dirichlet(\eqn{\boldsymbol{\alpha}}) distribution.
#'
#' @param n An integer of number of samples.
#' @param alpha A length-\eqn{m} vector of shape parameters. `alpha` should be positive.
#'
#' @returns
#' An \eqn{n \times m} matrix of random samples from the Dirichlet(\eqn{\boldsymbol{\alpha}}) distribution.
#'
#' @export
#' @keywords internal
rdirichlet <- function(n, alpha) {
  stopifnot("`alpha` must be numeric" = is.numeric(alpha))
  stopifnot("`alpha` must be a vector" = is.vector(alpha))
  stopifnot("`alpha` must be nonnegative." = all(alpha > 0))
  stopifnot("`n` must be greater than 0" = n > 0)

  # parameters in columns:
  alpha <- rbind(alpha)
  M <- ncol(alpha)

  # for a single set of parameters, extend to matrix
  if (n > nrow(alpha)) {
    alpha <- matrix(alpha, n, M, byrow = TRUE)
  }

  x <- matrix(stats::rgamma(M * n, alpha), nrow = n, ncol = M)

  # normalize to sum up to one:
  x / rowSums(x)
}
