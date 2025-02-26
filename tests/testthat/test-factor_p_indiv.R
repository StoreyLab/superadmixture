test_that("factor_p_indiv works", {

  EPSILON <- 0.1

  m_loci <- 10000
  n_inds <- 200
  X <- matrix(rbinom(m_loci * n_inds, size = 2, prob = runif(m_loci * n_inds, 0.01, 0.5)), m_loci, n_inds)
  expect_no_condition(check_geno(X))

  obj1 <- est_p_indiv(X, k_antepops = 3)
  obj2 <- factor_p_indiv(obj1$p_indiv, 3, init_method = "rowspace", rowspace = obj1$rowspace, verbose = FALSE)
  p_indiv <- obj1$p_indiv
  P_hat <- obj2$P_hat
  Q_hat <- obj2$Q_hat
  expect_lt(norm(p_indiv - P_hat %*% Q_hat, "F") / norm(p_indiv, "F"), EPSILON)

  k_antepops <- 2
  inbr_antepops <- c(0.1, 0.3)
  admix_proportions <- bnpsd::admix_prop_1d_linear(n_inds, k_antepops, sigma = 1)
  coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_antepops)
  out <- bnpsd::draw_all_admix(admix_proportions, inbr_antepops, m_loci)
  X <- out$X
  expect_no_condition(check_geno(X))

  obj1 <- est_p_indiv(X, k_antepops = 2)
  obj2 <- factor_p_indiv(obj1$p_indiv, 2, init_method = "rowspace", rowspace = obj1$rowspace, verbose = FALSE)
  p_indiv <- obj1$p_indiv
  P_hat <- obj2$P_hat
  Q_hat <- obj2$Q_hat
  expect_lt(norm(p_indiv - P_hat %*% Q_hat, "F") / norm(p_indiv, "F"), EPSILON)

})
