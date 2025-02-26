#' Calculate the negative entropy of sHWE p-values
#'
#' @description
#' `calc_neg_entropy()` groups the sHWE p-values into equal-sized bins to quantify how Uniform distributed they are over a given range. The number of bins is denoted by C in the algorithm; note that while we found \eqn{C = 150} to be sufficient for analyses with \eqn{10^5 < m < 10^6}, it may be helpful to choose a higher value if there are many more SNPs. This algorithm first divide the range of P-values into C equal-sized bins, and drop the bin with the lower bound of zero. Then for each of the remaining bins, compute the proportion of p-values in each bin. Finally, compute the entropy using the formula \deqn{\sum_{c=2}^C f_c\log f_c}
#'
#' @param pvalues A length-\eqn{m} (number of loci) vector of sHWE p-values
#' @param C An interger of number of bins
#' @return Negative entropy
#'
#' @examples
#' ## Calculating the entropy of sHWE p-values on AMR subset of TGP for K = 3 ---------------------
#' data("X_amr", package = "superadmixture")
#' rowspace <- est_p_indiv(X_amr, k_antepops = 3, loci_on_cols = TRUE, rowspace_only = TRUE)
#' p_values <- sHWE(X_amr, 3, rowspace, TRUE)
#' calc_neg_entropy(p_values)
#' @export
calc_neg_entropy <- function(pvalues, C=150) {

  # Divide the range of P-values into C equal-sized bins.
  breaks <- seq(0, 1, by = 1/C)
  p <- table(cut(pvalues, breaks = breaks))

  # Drop the bin with the lower bound of zero,
  # since this bin should contain the most significant P-values.
  p <- p[-1]

  # For each of the remaining  bins, compute the proportion of
  # P-values in each bin. These proportions should sum to one.
  p <- p / sum(p)

  # sing these proportions, compute the entropy using the formula
  return(sum(p * log10(p)))
}
