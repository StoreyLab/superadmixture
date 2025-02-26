init_coanc <- function(coanc_indiv, admix_proportions, init_method) {
  stopifnot('`init_method` %in% c("random", "LS")' = init_method %in% c("random", "LS"))
  if (init_method == "random") {
    return(diag(stats::runif(nrow(admix_proportions))))
  }
  if (init_method == "LS") {
    qqt <- tcrossprod(admix_proportions)
    q_theta_qt <- admix_proportions %*% coanc_indiv %*% t(admix_proportions)
    return(MASS::ginv(qqt) %*% q_theta_qt %*% MASS::ginv(qqt))
  }
}


calc_err <- function(coanc_indiv, coanc_antepops, admix_proportions) {
  return(norm(coanc_indiv - t(admix_proportions) %*% coanc_antepops %*% admix_proportions, "F"))
}


#' Estimate the coancestry among antecedent populations
#'
#' @description
#' `est_coanc()` infers the coancestry among antecedent populations \eqn{\boldsymbol{\Lambda}}.
#'
#' @details
#' When the super-admixture model is assumed, `est_coanc()` seeks a \eqn{K \times K} matrix
#' \eqn{\boldsymbol{\Lambda}} that solves this optimization problem:
#' \deqn{\arg\min\|\boldsymbol{\Theta} - \boldsymbol{Q}' \boldsymbol{\Lambda} \boldsymbol{Q}\|_F}
#' subject to the constraints:
#' \enumerate{
#'     \item \eqn{\boldsymbol{\Lambda}} is symmetric
#'     \item \eqn{0 \leq \lambda_{uv} \leq 1} for all \eqn{u,v}
#' }
#' When the standard admixture model is assumed, `est_coanc()` finds a diagonal \eqn{\boldsymbol{\Lambda}}
#' that minimizes:
#' \deqn{\arg\min\|\boldsymbol{\Theta} - \boldsymbol{Q}' \boldsymbol{\Lambda} \boldsymbol{Q}\|_F}
#' subject to the constraints:
#' \enumerate{
#'     \item \eqn{0 \leq \lambda_{uu} \leq 1} for all \eqn{u}
#'     \item \eqn{\lambda_{uv} = 0} for all \eqn{u \neq v}
#' }
#'
#' @param coanc_indiv An \eqn{n \times n} matrix representing coancestry among individuals.
#' @param admix_proportions A \eqn{K} \eqn{\times} \eqn{n} matrix representing admixture proportions.
#' @param tol A scalar for convergence criterion. If the relative change of
#' \eqn{\|\boldsymbol{\Theta} - \boldsymbol{Q}' \boldsymbol{\Lambda} \boldsymbol{Q}\|_F} is less than `tol`,
#' the algorithm halts.
#' @param max_iters Maximum number of iterations.
#' @param init_method One of "random" or "LS" (default). If `init_method` is "random",
#' the initial \eqn{\boldsymbol{\Lambda}} is set to be a diagonal matrix whose diagonal elements
#' are sampled from the Uniform(0,1) distribution. Otherwise, the initial \eqn{\boldsymbol{\Lambda}} is set to be
#' \eqn{(\boldsymbol{Q}\boldsymbol{Q}^{\prime})^{-1} \boldsymbol{\Theta} (\boldsymbol{Q}\boldsymbol{Q}^{\prime})^{-1}}.
#' @param model One of "standard" or "super" (default). If `model` is "standard", the
#' standard admixture model is fitted. If `model` is "super", the the super-admixture model is fitted.
#' @param verbose A Boolean that controls the returned values. If `verbose` is FALSE (default), only the estimated \eqn{\boldsymbol{\Lambda}}
#' is returned. If TRUE, `est_coanc` returns a list consisting of the estimated \eqn{\boldsymbol{\Lambda}},
#' a list of \eqn{\boldsymbol{\Lambda}}s over iterations, a list of errors over iterations.
#'
#' @returns If `verbose` is FALSE, the estimated \eqn{\boldsymbol{\Lambda}} is returned.
#'
#' If `verbose` is TRUE, a list with the following elements is returned:
#' + `coanc_antepops`: The estimated \eqn{\boldsymbol{\Lambda}},
#' + `coanc_list`: A list of \eqn{\boldsymbol{\Lambda}}'s over iterations,
#' + `err_list`: A list of \eqn{\|\boldsymbol{\Theta} - \boldsymbol{Q}' \boldsymbol{\Lambda} \boldsymbol{Q}\|_F}'s over iterations.
#'
#' @examples
#' ## Estimate Coancestry among antecedent populations of HGDP ----------------------------------------
#' data("X_hgdp", package = "superadmixture")
#' data("admix_props_hgdp", package = "superadmixture")
#' ## Estimate kinship coefficients
#' kinship <- popkin::popkin(t(X_hgdp))
#' coanc_indiv <- popkin::inbr_diag(kinship)
#' coanc_indiv <- ifelse(coanc_indiv < 0, 0, coanc_indiv)
#'
#' ## Estimate coancestry among populations
#' coanc_pops_hgdp <- est_coanc(coanc_indiv, admix_props_hgdp, model = "super")
#
#' @export
est_coanc <- function(coanc_indiv, admix_proportions, tol = 1e-8, max_iters = 10000, init_method = "LS", model = "super", verbose = FALSE) {
  # check whether the parameter `init_method` is valid
  stopifnot('`init_method` %in% c("random", "LS")' = init_method %in% c("random", "LS"))
  # check whether the parameter `model` is valid
  stopifnot('`model` %in% c("super", "standard")' = model %in% c("super", "standard"))
  # check if `coanc_indiv` is valid
  check_coancestry(coanc_indiv)
  # check if `admix_proportions` is valid
  check_admix_proportions(admix_proportions)

  # define `k_antepops` and `n_inds`
  if (nrow(admix_proportions) == nrow(coanc_indiv)) {
    admix_proportions <- t(admix_proportions)
  }
  k_antepops <- nrow(admix_proportions)
  n_inds <- ncol(admix_proportions)

  # pre-calculate the variables that keep fixed through the optimization
  q_theta_qt <- admix_proportions %*% coanc_indiv %*% t(admix_proportions)
  qqt <- tcrossprod(admix_proportions)
  if (nrow(admix_proportions) == 2) {
    lps <- 2 * max(svd(admix_proportions)$d)^ 4
  } else {
    lps <- 2 * rARPACK::svds(admix_proportions, k = 1, nu = 0, nv = 0)$d^4
  }
  coanc_antepops <- init_coanc(coanc_indiv, admix_proportions, init_method)

  # run optimization to estimate coancestry
  converged <- FALSE
  previous_err <- calc_err(coanc_indiv, coanc_antepops, admix_proportions)
  coanc_list <- list(length = 0)
  err_list <- list(length = 0)
  for (i in 1:max_iters) {
    grd <- -2 * (q_theta_qt - qqt %*% coanc_antepops %*% qqt)
    coanc_antepops <- coanc_antepops - 1 / lps * grd
    idx_offdiag <- col(coanc_antepops) != row(coanc_antepops)
    if (model == "standard") {
      coanc_antepops[idx_offdiag] <- 0
    }
    diag(coanc_antepops) <- ifelse(diag(coanc_antepops) < 0.01, 0.01, diag(coanc_antepops))
    coanc_antepops[idx_offdiag] <- ifelse(coanc_antepops[idx_offdiag] < 0, 0, coanc_antepops[idx_offdiag])
    coanc_antepops <- ifelse(coanc_antepops > 0.99, 0.99, coanc_antepops)
    current_err <- calc_err(coanc_indiv, coanc_antepops, admix_proportions)
    coanc_list[[i]] <- coanc_antepops
    err_list[[i]] <- current_err

    converged <- (previous_err - current_err) / previous_err < tol
    previous_err <- current_err
    if (converged) break
  }

  # symmetrize the result
  coanc_antepops <- (coanc_antepops + t(coanc_antepops)) / 2
  if (verbose) {
    return(list(coanc_antepops = coanc_antepops, err_list = err_list, coanc_list = coanc_list))
  } else {
    return(coanc_antepops)
  }
}
