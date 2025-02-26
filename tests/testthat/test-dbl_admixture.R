test_that("dbl_admixture() works", {

  skip_on_ci()

  EPSILON <- 5e-2

  draw_p_antepops <- function(S, k_antepops, B = 100000) {
    p_anc <- stats::runif(1)

    Gamma <- diag(stats::runif(S))
    W <- t(superadmixture:::rdirichlet(n = k_antepops, alpha = rep(1, S)))
    coanc_antepops <- t(W) %*% Gamma %*% W

    admix_props <- construct_admix_props(k_antepops, alpha_type = "1")

    # draw allele frequencies of antecedent populations
    p_antepops <- dbl_admixture(rep(p_anc, B), admix_proportions = admix_props, coanc_antepops = coanc_antepops, p_antepops_only = TRUE, S = S)
    p_tru <- rep(p_anc, k_antepops)
    p_est <- colMeans(p_antepops)
    covar_tru <- p_anc * (1 - p_anc) * coanc_antepops
    covar_est <- cov(p_antepops)

    return(list(p_tru=p_tru, p_est=p_est, covar_tru=covar_tru, covar_est=covar_est))
  }

  obj <- draw_p_antepops(S = 6,  k_antepops = 3)
  expect_lt(norm(obj$p_est - obj$p_tru, "2") / norm(obj$p_tru, "2"), EPSILON)
  expect_lt(norm(obj$covar_est - obj$covar_tru, "F") / norm( obj$covar_tru, "F"), EPSILON)

  obj <- draw_p_antepops(S = 10, k_antepops = 6)
  expect_lt(norm(obj$p_est - obj$p_tru, "2") / norm(obj$p_tru, "2"), EPSILON)
  expect_lt(norm(obj$covar_est - obj$covar_tru, "F") / norm( obj$covar_tru, "F"), EPSILON)

})


