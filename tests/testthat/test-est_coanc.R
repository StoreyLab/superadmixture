test_that("coanc_est() can find the optimal solution when `model` = 'super'", {

  EPSILON <- 0.1

  for (k_antepops in c(3, 6, 9)) {
    for (alpha_type in c("1", "2", "spatial")) {
      coanc_antepops <- construct_coanc(k_antepops)
      admix_props <- construct_admix_props(k_antepops, alpha_type = alpha_type)
      coanc_indiv <- t(admix_props) %*% coanc_antepops %*% admix_props
      coanc_antepops_est <- est_coanc(coanc_indiv, admix_props, init_method = "random", tol = 1e-6, verbose = FALSE)
      expect_lt(norm(coanc_antepops_est - coanc_antepops, "F") / norm(coanc_antepops, "F"), EPSILON)
    }
  }
})


test_that("coanc_est() can find the optimal solution when `model` = 'standard'", {

  EPSILON <- 0.1

  for (k_antepops in c(3, 6, 9)) {
    for (alpha_type in c("1", "2", "spatial")) {
      coanc_antepops <- diag(runif(k_antepops))
      admix_props <- construct_admix_props(k_antepops, alpha_type = alpha_type)
      coanc_indiv <- t(admix_props) %*% coanc_antepops %*% admix_props
      coanc_antepops_est <- est_coanc(coanc_indiv, admix_props, init_method = "random", model = "standard", tol = 1e-6, verbose = FALSE)
      expect_lt(norm(coanc_antepops_est - coanc_antepops, "F") / norm(coanc_antepops, "F"), EPSILON)
    }
  }

})
