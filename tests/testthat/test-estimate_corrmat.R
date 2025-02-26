test_that("estimate_corrmat() works", {

  skip_on_ci()

  EPSILON <- 5e-2

  draw_p_antepops <- function(S, k_antepops, B = 100000) {
    Gamma <- diag(stats::runif(S))
    p_anc <- stats::runif(1)
    W <- t(superadmixture:::rdirichlet(n = k_antepops, alpha = rep(1, S)))
    coanc_antepops <- t(W) %*% Gamma %*% W
    corrmat <- estimate_corrmat(p_anc, coanc_antepops)
    samples <- matrix(nrow=B, ncol=k_antepops)
    for (iter in 1:B) {
      z <- MASS::mvrnorm(n = 1, mu = rep(0, k_antepops), Sigma = corrmat)
      a <- (1 - diag(coanc_antepops)) / diag(coanc_antepops) * p_anc
      b <- (1 - diag(coanc_antepops)) / diag(coanc_antepops) * (1 - p_anc)
      samples[iter, ] <- qbeta(pnorm(z), a, b)
    }
    p_tru <- rep(p_anc, k_antepops)
    p_est <- colMeans(samples)
    covar_tru <-  p_anc * (1 - p_anc) * coanc_antepops
    covar_est <- cov(samples)
    return(list(samples=samples, p_tru=p_tru, p_est=p_est, covar_tru=covar_tru, covar_est=covar_est))
  }

  obj <- draw_p_antepops(S=6, k=3)
  expect_lt(norm(obj$p_tru - obj$p_est, "2") / norm(obj$p_tru, "2"), EPSILON)
  expect_lt(norm(obj$covar_tru - obj$covar_est, "F") / norm(obj$covar_tru, "F"), EPSILON)
  expect_true(all(obj$samples <= 1) & all(obj$samples >= 0))

})

