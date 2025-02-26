pvals_empir <- function(stats1, stats0) {
  if (missing(stats1))
    stop("`stats1` observed statistics vector is required!")
  if (missing(stats0))
    stop("`stats0` null statistics (vector or matrix) is required!")
  stats1_order <- order(stats1)
  stats1_sorted <- stats1[stats1_order]
  stats0 <- sort(as.vector(stats0))
  m0 <- length(stats0)

  m <- length(stats1)
  pvals <- rep(NA, m)
  i0 <- 1
  for (i1 in seq_len(m)) {
    stats1_sorted_i1 <- stats1_sorted[i1]
    if (is.na(stats1_sorted_i1))
      break

    while ((i0 <= m0) && (stats1_sorted_i1 >= stats0[i0])) {
      i0 <- i0 + 1
    }
    pvals[stats1_order[i1]] <- 1 - ((i0 - 1)/m0)
  }
  return(pvals)
}


gof_stat_snp <- function(snp, p) {
  p0 <- (1 - p)^2
  p1 <- 2 * p * (1 - p)
  est <- c(sum(p0), sum(p1))
  N <- c(sum(snp == 0), sum(snp == 1))

  sigma11 <- sum(p0 * (1 - p0))
  sigma12 <- -sum(p0 * p1)
  sigma22 <- sum(p1 * (1 - p1))

  determ <- sigma11 * sigma22 - sigma12^2
  if (determ == 0)
    return(NA)

  Sigma_inv <- c(sigma22, -sigma12, -sigma12, sigma11)
  Sigma_inv <- matrix(Sigma_inv, nrow=2, ncol=2)/determ
  stat <- t(N - est) %*% Sigma_inv %*% (N - est)
  return(stat)
}


#' @title Hardy-Weinberg Equilibrium in structure populations
#'
#' @description
#' Compute structured Hardy-Weinberg Equilibrium (sHWE) p-values on a SNP-by-SNP basis.
#'
#' @inheritParams est_p_indiv
#' @param rowspace An optimal \eqn{n \times K} matrix of eigen-vectors \eqn{\boldsymbol{W}}. See details in [est_p_indiv()].
#'
#'
#' @return A vector of sHWE p-values for each SNP.
#' @examples
#' # Performing sHWE test on AMR subset of TGP for K = 3 -------------------------------------
#' data("X_amr", package = "superadmixture")
#' rowspace <- est_p_indiv(X_amr, k_antepops = 3, loci_on_cols = TRUE, rowspace_only = TRUE)
#' p_values <- sHWE(X_amr, 3, rowspace, TRUE)
#' hist(p_values)
#' @export
sHWE <- function(X, k_antepops, rowspace = NULL, loci_on_cols = FALSE) {
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

  stopifnot("nrow(rowspace) == ncol(X)" == is.null(rowspace) || (!is.null(rowspace) && nrow(rowspace) == ncol(X)))
  # calculate rowspace
  if (is.null(rowspace))
    rowspace <- est_p_indiv(X, k_antepops, loci_on_cols = FALSE, rowspace_only = TRUE)

  # calculate observed stats across matrix
  p_indiv <- 0.5 * X %*% rowspace %*% t(rowspace)
  p_indiv <- ifelse(p_indiv > 1, 1, ifelse(p_indiv < 0, 0, p_indiv))
  stats1 <- mapply(gof_stat_snp, asplit(X, 1), asplit(p_indiv, 1))

  rm(X)
  X0 <- matrix(stats::rbinom(length(p_indiv), 2, p_indiv), nrow(p_indiv), ncol(p_indiv))
  p_indiv0 <- est_p_indiv(X0, k_antepops, loci_on_cols = FALSE, p_indiv_only = TRUE)
  stats0 <- mapply(gof_stat_snp, asplit(X0, 1), asplit(p_indiv0, 1))

  # calculate empirical p-values based on these distributions
  pvals <- pvals_empir(stats1, stats0)

  return(pvals)
}
