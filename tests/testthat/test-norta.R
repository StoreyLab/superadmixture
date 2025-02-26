test_that("norta() works", {

  skip_on_ci()

  EPSILON <- 5e-1

  draw_p_antepops <- function(S, k_antepops) {
    Gamma <- diag(stats::runif(S))
    p_anc <- stats::runif(1)
    p_anc <- rep(p_anc, 100)
    W <- t(superadmixture:::rdirichlet(n = k_antepops, alpha = rep(1, S)))
    coanc_antepops <- t(W) %*% Gamma %*% W
    p_antepops <- norta(p_anc, coanc_antepops, parallel = TRUE, p_antepops_only = TRUE)

    p_tru <- rep(p_anc[1], k_antepops)
    p_est <- colMeans(p_antepops)
    covar_tru <-  p_anc[1] * (1 - p_anc[1]) * coanc_antepops
    covar_est <- cov(p_antepops)
    return(list(p_antepops=p_antepops, p_tru=p_tru, p_est=p_est, covar_tru=covar_tru, covar_est=covar_est))
  }

  obj <- draw_p_antepops(S=6, k=3)
  expect_lt(norm(obj$p_tru - obj$p_est, "2") / norm(obj$p_tru, "2"), EPSILON)
  expect_lt(norm(obj$covar_tru - obj$covar_est, "F") / norm(obj$covar_tru, "F"), EPSILON)
  expect_true(all(obj$samples <= 1) & all(obj$samples >= 0))

})
