test_that("est_paras_dbladmix() errors if the input `coanc_antepops` is a scalar", {
  expect_error(est_paras_dbladmix(3))
})

test_that("est_paras_dbladmix() works", {

  EPSILON <- 1e-1

  # a wrapper function
  est_paras <- function(S, k_antepops) {
    Gamma <- diag(stats::runif(S))
    W <- t(superadmixture:::rdirichlet(n = k_antepops, alpha = rep(1, S)))
    coanc_antepops <- t(W) %*% Gamma %*% W
    paras <- est_paras_dbladmix(coanc_antepops, S, verbose = FALSE)
    return(list(W=paras$W, Gamma=paras$Gamma, coanc_antepops=coanc_antepops))
  }

  obj <- est_paras(S = 6, k_antepops = 3)
  expect_lt(norm(t(obj$W) %*% obj$Gamma %*% obj$W - obj$coanc_antepops, "F") / norm(obj$coanc_antepops, "F"), EPSILON)
  expect_no_condition(check_admix_proportions(obj$W))
  expect_true(all(diag(obj$Gamma) < 1) && all(diag(obj$Gamma) > 0))

  obj <- est_paras(S = 10, k_antepops = 6)
  expect_lt(norm(t(obj$W) %*% obj$Gamma %*% obj$W - obj$coanc_antepops, "F") / norm(obj$coanc_antepops, "F"), EPSILON)
  expect_no_condition(check_admix_proportions(obj$W))
  expect_true(all(diag(obj$Gamma) < 1) && all(diag(obj$Gamma) > 0))
})


