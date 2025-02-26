D_binom <- function(X) {
  m_loci <- nrow(X)
  n_inds <- ncol(X)
  d <- apply(X, 2, function(xc) sum(2 * xc - xc^2) / m_loci)
  return(diag(d))
}


#' Estimate Individual specific Allele Frequencies
#'
#' @description
#' `est_p_indiv()` estimates the individual specific allele frequencies \eqn{\boldsymbol{\Pi}} from a genotype matrix \eqn{\boldsymbol{X}}.
#'
#' @details
#' `est_p_indiv()` adopts a three-step strategy to estimate the individual specific allele frequencies.
#' \enumerate{
#'  \item If `loci_on_cols`, \eqn{\boldsymbol{X} \leftarrow t(\boldsymbol{X})}.
#'  \item Calculates the top \eqn{K} eigen-vectors \eqn{\boldsymbol{W}} of the matrix \eqn{\boldsymbol{G} = \frac{1}{m} \boldsymbol{X}^{\prime}\boldsymbol{X} - \boldsymbol{D}}
#'  where the matrix \eqn{\boldsymbol{D}} is a diagonal matrix used to adjust for heteroskedasticity for the Binomial distribution. More specifically,
#'  \eqn{d_{jj}} is defined as \deqn{d_{jj} = \frac{1}{m}\sum_{i = 1}^m (2 x_{ij} - x_{ij}^2)} where \eqn{x_{ij}} represents the genotypic value for
#'  the individual \eqn{j} at the \eqn{i}th locus.
#'  \item Estimate the individual specific allele frequencies \eqn{\boldsymbol{\Pi}}  by projecting \eqn{\boldsymbol{X}} onto \eqn{\boldsymbol{W}}:
#'  \eqn{\boldsymbol{\Pi} = \frac{1}{2} \boldsymbol{X}\boldsymbol{W}\boldsymbol{W}^{\prime}}. If \eqn{\pi_{ij} > 1} or \eqn{\pi_{ij} < 0}, \eqn{\pi_{ij}} is set to be
#'  \eqn{\min(1, \max(0, \pi_{ij}))}.
#' }
#'
#' @param X A matrix of genotypes. `X` should be an integer matrix of 0's, 1's, 2's or `NA`.
#' @param k_antepops An integer representing the number of antecedent populations.
#' @param loci_on_cols A boolean that specifies the dimension of the input `X`. If `FALSE` (default), `X` has \eqn{m} (number of loci) \eqn{\times} \eqn{n} (number of individuals) dimension.
#' If `TRUE`, `X`'s dimension is \eqn{n \times m}.
#' @param rowspace_only,p_indiv_only Booleans for controlling the return values. If `rowspace_only = T`, an \eqn{n \times K} matrix of eigen-vectors \eqn{\boldsymbol{W}}
#' of \eqn{\boldsymbol{G}} is returned. If `p_indiv_only = T`, an \eqn{m \times n} matrix \eqn{\boldsymbol{\Pi}} of estimated individual specific allele frequencies is returned.
#' Otherwise, a list that containing \eqn{\boldsymbol{W}} and \eqn{\boldsymbol{\Pi}} is returned.
#'
#' @returns
#' + If `rowspace_only = T`, an \eqn{n \times K} matrix of eigen-vectors \eqn{\boldsymbol{W}} of \eqn{\boldsymbol{G}} is returned.
#' + If `p_indiv_only = T`,  an \eqn{m \times n} matrix \eqn{\boldsymbol{\Pi}} of estimated individual specific allele frequencies is returned.
#' + Otherwise, a named list is returned, containing:
#'    + `p_indiv`: The \eqn{m \times n} matrix of estimated IAFs,
#'    + `rowspace`: The \eqn{n \times K} matrix of eigen-vectors \eqn{\boldsymbol{W}} of \eqn{\boldsymbol{G}}.
#'
#' @examples
#' ## Estimate Individual Specific Allele Frequencies of HGDP  -----------------------------
#' ## It takes 6.451 sec to finish the following codes on a 2.3 GHz 18-Core Intel Xeon W iMac Pro.
#' data("X_hgdp", package = "superadmixture")
#' p_indiv <- est_p_indiv(X_hgdp, k_antepops = 7, loci_on_cols = TRUE, p_indiv_only = TRUE)
#' @export
est_p_indiv <- function(X, k_antepops, loci_on_cols = FALSE, rowspace_only = FALSE, p_indiv_only = FALSE) {
  # check whether X is valid
  check_geno(X)
  if (sum(is.na(X)) > 0) {
    if (loci_on_cols) {
      ind <- which(is.na(X), arr.ind=TRUE)
      X[ind] <- colMeans(X,  na.rm = TRUE)[ind[,2]]
    } else {
      ind <- which(is.na(X), arr.ind=TRUE)
      X[ind] <- rowMeans(X,  na.rm = TRUE)[ind[,1]]
    }
  }
  if (loci_on_cols) X <- t(X)
  m_loci <- nrow(X)

  # estimate the rowspace
  rowspace <- rARPACK::eigs_sym(crossprod(X) / m_loci - D_binom(X), k = k_antepops)$vectors
  if (rowspace_only) {
    return(rowspace)
  }

  # estimate IAFs
  p_indiv <- 0.5 * X %*% rowspace %*% t(rowspace)
  p_indiv <- ifelse(p_indiv > 1, 1, ifelse(p_indiv < 0, 0, p_indiv))
  if (p_indiv_only) {
    return(p_indiv)
  }

  return(list(rowspace = rowspace, p_indiv = p_indiv))
}
