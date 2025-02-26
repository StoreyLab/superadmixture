test_that("est_p_indiv works", {
  m_loci <- 10000
  n_inds <- 200
  X <- matrix(rbinom(m_loci * n_inds, size = 2, prob = runif(m_loci * n_inds, 0.01, 0.5)), m_loci, n_inds)
  expect_no_condition(check_geno(X))

  # check if D_binom() works
  expect_equal(D_binom(X), D_binomial(X))

  # check if est_p_indiv() works
  expect_equal(est_p_indiv(X, k_antepops = 3, p_indiv_only = TRUE), estimate_F(X, d = 3)$F_hat)


  k_subpops <- 2
  inbr_subpops <- c(0.1, 0.3)
  admix_proportions <- bnpsd::admix_prop_1d_linear(n_inds, k_subpops, sigma = 1)
  coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_subpops)
  out <- bnpsd::draw_all_admix(admix_proportions, inbr_subpops, m_loci)
  X <- out$X
  expect_no_condition(check_geno(X))

  # check if D_binom() works
  expect_equal(D_binom(X), D_binomial(X))

  # check if est_p_indiv() works
  expect_equal(est_p_indiv(X, k_antepops = 2, p_indiv_only = TRUE), estimate_F(X, d = 2)$F_hat)
})
