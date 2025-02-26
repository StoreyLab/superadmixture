test_that("sHWE() works", {
  m_loci <- 10000
  n_inds <- 200
  k_antepops <- 3
  X <- matrix(rbinom(m_loci * n_inds, size = 2, prob = runif(m_loci * n_inds, 0.01, 0.5)), m_loci, n_inds)
  rowspace <- est_p_indiv(X, k_antepops, rowspace_only = TRUE)
  pval1 <- sHWE(X, k_antepops, rowspace = rowspace)
  pval2 <- lfa::sHWE(X, rowspace, B = 1)
  expect_gt(cor(pval1, pval2), 0.95)


  k_antepops <- 2
  inbr_subpops <- c(0.1, 0.3)
  admix_proportions <- bnpsd::admix_prop_1d_linear(n_inds, k_antepops, sigma = 1)
  coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_subpops)
  out <- bnpsd::draw_all_admix(admix_proportions, inbr_subpops, m_loci)
  X <- out$X
  rowspace <- est_p_indiv(X, k_antepops, rowspace_only = TRUE)
  pval1 <- sHWE(X, k_antepops, rowspace = rowspace)
  pval2 <- lfa::sHWE(X, rowspace, B = 1)
  expect_gt(cor(pval1, pval2), 0.95)
})
