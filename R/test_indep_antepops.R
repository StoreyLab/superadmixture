#' Compute a null test statistics for testing of no coancestry among antecedent populations
#'
#' @description
#' `compute_null()` computes a null test statistic for testing of no coancestry among antecedent populations.
#'
#' @inheritParams test_indep_antepops
#' @inheritParams dbl_admixture
#' @param inbr_antepops A length-\eqn{K} vector of inbreeding coefficient.
#'
#' @return A test statistic under the null.
#' @seealso
#' [test_indep_antepops()] for hypothesis testing of no coancestry among antecedent populations.
#' @export
#' @keywords internal
compute_null <- function(p_anc, admix_proportions, inbr_antepops, verbose = FALSE, ...) {

  # draw alelle frequencies of antecedent populations
  p_antepops <- bnpsd::draw_p_subpops(p_anc, inbr_antepops)
  p_ind     <- bnpsd::make_p_ind_admix(p_antepops, t(admix_proportions))
  X0 <- bnpsd::draw_genotypes_admix(p_ind)
  coanc_indiv0 <- popkin::inbr_diag(popkin::popkin(X0, ...))

  # estimate the coancestry among antecedent populations
  coanc_antepops_sup0 <- est_coanc(coanc_indiv0, admix_proportions)
  coanc_antepops_std0 <- est_coanc(coanc_indiv0, admix_proportions, model = "standard")

  # estimate the test statistics
  test_stat0 <- norm(coanc_antepops_sup0 - coanc_antepops_std0, "F")
  if (verbose) print(paste0("T0: ", test_stat0))
  return(test_stat0)
}


#' Test of no coancestry among antecedent populations
#'
#' @description
#' `test_indep_antepops()` performs a hypothesis testing on whether the assumption of no coancestry
#' among antecedent populations holds.
#'
#' @inheritParams est_p_indiv
#' @inheritParams est_coanc
#' @param B An optional integer (defaults to 100) of number of bootstrap replications.
#' @param ... Other options used to control the estimation of coancestry in [popkin::popkin()]. Passed on to [popkin::popkin()].
#'
#' @return A named list is returned, containing:
#' - `test_stat1`: A scalar corresponding to the observed test statistics
#' - `test_stats0`: A length-\eqn{B} vector corresponding to the test statistics
#'   under the null
#' - `p_value`: A scalar that may be interpreted as the bootstrap p-value
#'
#' @export
test_indep_antepops <- function(X, admix_proportions, loci_on_cols = TRUE, B = 100, verbose = FALSE, ...) {

  # check whether X is valid
  check_geno(X)
  if (loci_on_cols) X <- t(X)

  # check if `admix_proportions` is valid
  check_admix_proportions(admix_proportions)
  if (ncol(admix_proportions) != ncol(X)) {
      admix_proportions <- t(admix_proportions)
  }

  # calculate the ancestral allele frequencies
  p_anc <- rowMeans(X, na.rm=TRUE) / 2

  # estimate the matrix of the coancestry
  coanc_indiv <- popkin::inbr_diag(popkin::popkin(X, ...))

  # estimate the coancestry among antecedent populations
  coanc_antepops_sup <- est_coanc(coanc_indiv, admix_proportions)
  coanc_antepops_std <- est_coanc(coanc_indiv, admix_proportions, model = "standard")

  # estimate the test statistics
  test_stat1 <- norm(coanc_antepops_sup - coanc_antepops_std, "F")
  if (verbose) print(paste0("T: ", test_stat1))

  # generate the null distributions
  test_stats0 <- c()
  for (b in 1:B) {
    test_stats0 <- c(test_stats0, compute_null(p_anc, admix_proportions, diag(coanc_antepops_std), verbose, ...))
  }

  p_value <- sum(test_stat1 < test_stats0) / B
  return(list(test_stat1=test_stat1, test_stats0=test_stats0, p_value=p_value))
}
