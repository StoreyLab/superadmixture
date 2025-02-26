#' Estimating parameters used in the double-admixture simulation
#'
#' @description
#' `est_paras_dbladmix()` identifies a \eqn{K \times S} matrix \eqn{\boldsymbol{W}} and a diagonal \eqn{S \times S}
#' matrix \eqn{\boldsymbol{\Gamma}} such that \eqn{\boldsymbol{\Lambda} \approx \boldsymbol{W}^{\prime}\boldsymbol{\Gamma}\boldsymbol{W}}.
#'
#' @details
#' This function solves the following optimization problem based on Proximal Alternating Linearized Minimization (PALM),
#' and it is guaranteed to converge to a stationary point:
#' \deqn{\arg\min \|\boldsymbol{\Lambda} - \boldsymbol{W}^{\prime}\boldsymbol{\Gamma}\boldsymbol{W}\|_F}
#' subjects to those constraints:
#' \enumerate{
#'     \item All entries of \eqn{\boldsymbol{W}} are greater or equal to 0
#'     \item \eqn{\sum_{s = 1}^S w_{su} = 1} for all \eqn{u}
#'     \item \eqn{0.01 \leq \gamma_{ss} \leq 0.99} for all \eqn{s}
#'     \item \eqn{\gamma_{st} = 0} for all \eqn{s \neq t}
#' }
#'
#' @param coanc_antepops A \eqn{K \times K} matrix representing coancestry among antecedent populations.
#' @param S An integer representing the number of independent Balding-Nichols distributions.
#' @param max_iters Maximum number of iterations. `max_iters` should be set to a large number (100,000 by default) to ensure convergence.
#' @param tau1,tau2 Scalars controlling the step size of the optimization iteration.
#' They should be greater than 1.0 to ensure the convergence (defaults to 1.1).
#' @param tol A Scalar for convergence criterion. If the relative change of
#' \eqn{\|\boldsymbol{\Lambda} - \boldsymbol{W}^{\prime}\boldsymbol{\Gamma}\boldsymbol{W}\|_F} is less than `tol`,
#' the algorithm halts.
#' @param verbose A Boolean that controls the returned values. If `verbose = F` (default),
#' only the estimated \eqn{\boldsymbol{W}} and \eqn{\boldsymbol{\Gamma}}  are returned. If `T`, this function returns
#' a list consisting of the estimated \eqn{\boldsymbol{W}}, the estimated \eqn{\boldsymbol{\Gamma}} and a list of errors over iterations.
#'
#' @returns
#' + If `verbose = T`, a list with the following elements is returned:
#'    + `W`: The estimated \eqn{\boldsymbol{W}}
#'    + `Gamma`: The estimated \eqn{\boldsymbol{\Gamma}}
#'    + `err_list`: A list of \eqn{\|\boldsymbol{\Lambda} - \boldsymbol{W}^{\prime}\boldsymbol{\Gamma}\boldsymbol{W}\|_F}'s over iterations.
#' + If `verbose = F`, a list containing the following elements is returned:
#'    + `W`: The estimated \eqn{\boldsymbol{W}}
#'    + `Gamma`: The estimated \eqn{\boldsymbol{\Gamma}}
#'
#' @examples
#' # estimate parameters  --------------------------------------------------------------------
#' data("coanc_pops_amr",  package = "superadmixture")
#' paras <- est_paras_dbladmix(coanc_pops_amr, verbose = FALSE)
#'
#' @export
#' @keywords internal
est_paras_dbladmix <- function(coanc_antepops, S = 2 * nrow(coanc_antepops), tau1 = 1.1, tau2 = 1.1, max_iters = 100000, tol = 1e-6, verbose = FALSE) {
  # check if `coanc_antepops` is valid
  check_coancestry(coanc_antepops)
  k_antepops <- nrow(coanc_antepops)
  stopifnot("`k_antepops` > 1" = k_antepops > 1)

  # calculate the L2-norm of `coanc_antepops`
  l2_coanc_antepops <- norm(coanc_antepops, "2")

  # initialize W and gamma
  W <- t(rdirichlet(n = k_antepops, alpha = rep(1, S)))
  Gamma <- diag(stats::runif(S))

  # run optimization to estimate parameters for dbl-admixture
  converged <- FALSE
  previous_err <- calc_err(coanc_antepops, Gamma, W)
  err_list <- list(length = 0)

  for (i in 1:max_iters) {
    # update Gamma
    lipschitz <- 2 * norm(W, "2")^4
    grd <- -2 * (W %*% (coanc_antepops - t(W) %*% Gamma %*% W) %*% t(W))
    Gamma <- Gamma - 1 / (tau1 * lipschitz) * grd
    Gamma[col(Gamma) != row(Gamma)] <- 0
    diag(Gamma) <- ifelse(diag(Gamma) > 0.99, 0.99, diag(Gamma))
    diag(Gamma) <- ifelse(diag(Gamma) < 0.01, 0.01, diag(Gamma))

    # update W
    l2_Gamma <- norm(Gamma, "2")
    lipschitz <- 4 * (l2_coanc_antepops * l2_Gamma + 3 * k_antepops * l2_Gamma^2)
    grd <- -4 * Gamma %*% W %*% (coanc_antepops - t(W) %*% Gamma %*% W)
    W <- W - 1 / (tau2 * lipschitz) * grd
    for (j in 1:ncol(W)) {
      W[, j] <- projsplx(W[, j])
    }

    # check convergence
    current_err <- calc_err(coanc_antepops, Gamma, W)
    if (verbose && i %% 100 == 1) print(paste0("Iter = ", i, ", current_error = ", round(current_err, digits = 3)))
    err_list[[i]] <- current_err
    converged <- (previous_err - current_err) / previous_err < tol
    previous_err <- current_err
    if (converged) break
  }

  if (verbose) {
    return(list(W = W, Gamma = Gamma, err_list = err_list))
  } else {
    return(list(W = W, Gamma = Gamma))
  }
}


#' Simulating allele frequencies and genotypes by the double-admixture method
#'
#' @description
#' `dbl_admixture()` generates random samples of antecedent allele frequencies \eqn{\boldsymbol{P}} and random samples of genotypes \eqn{\boldsymbol{X}}
#' according to the double-admixture method.
#'
#' @param p_anc A length-\eqn{m} vector representing allele frequencies at the ancestral population.
#' @inheritParams est_coanc
#' @inheritParams est_paras_dbladmix
#' @param paras An optional named list containing `Gamma` and `W`. If provided, `Gamma` and `W` won't be estimated. See details in [est_paras_dbladmix].
#' @param paras_only,p_antepops_only,geno_only Booleans that control the returned values. If `paras_only = T`, the estimated \eqn{\boldsymbol{W}}
#' and \eqn{\boldsymbol{\Gamma}} are returned. If `p_antepops_only = T`, only the \eqn{m \times K} matrix of antecedent allele frequencies
#' \eqn{\boldsymbol{P}} is returned. If `geno_only = T`, only the \eqn{m \times n} matrix of genotypes \eqn{\boldsymbol{X}} is returned. Otherwise,
#' a list containing \eqn{\boldsymbol{W}}, \eqn{\boldsymbol{\Gamma}}, \eqn{\boldsymbol{P}} and \eqn{\boldsymbol{X}} is returned.
#' @param ... Other options used to control parameter estimation. Passed on to [est_paras_dbladmix].
#'
#' @returns
#' + If `paras_only = T`, returns parameters estimated from [est_paras_dbladmix].
#' + If `p_antepops = T`, an \eqn{m \times K} matrix of allele frequencies of antecedent populations \eqn{\boldsymbol{P}} is returned.
#' + If `geno_only = T`, an \eqn{m \times n} matrix of genotypes is returned.
#' + Otherwise, a named list is returned, containing:
#'    + `W` and `Gamma`: the parameters estimated from [est_paras_dbladmix],
#'    + `p_antepops`: the simulated \eqn{m \times K} matrix of antecedent allele frequencies,
#'    + `X`: the simulated \eqn{m \times n} matrix of genotypes.
#'
#' @examples
#' ## Generate a genotype from HGDP -------------------------------------------------------------------
#' ## It takes 87.842 sec to finish the following codes on a 2.3 GHz 18-Core Intel Xeon W iMac Pro.
#' \dontrun{
#' data("fam_hgdp",         package = "superadmixture")
#' data("p_anc_hgdp",       package = "superadmixture")
#' data("admix_props_hgdp", package = "superadmixture")
#' data("coanc_pops_hgdp",  package = "superadmixture")
#' X_hgdp <- dbl_admixture(
#'       p_anc = p_anc_hgdp,
#'       coanc_antepops = coanc_pops_hgdp,
#'       admix_proportions = admix_props_hgdp,
#'       geno_only = TRUE)
#'
#' ## Estimate kinship from the simulated genotypes
#' kinship <- popkin::popkin(X_hgdp)
#' coanc_indiv <- popkin::inbr_diag(kinship)
#' coanc_indiv <- ifelse(coanc_indiv < 0, 0, coanc_indiv)
#'
#' ## Visualize kinship
#' popkin::plot_popkin(coanc_indiv, labs = fam_hgdp$subpop, labs_las = 2)}
#' @seealso
#' [est_paras_dbladmix] for details about estimating parameters of the Double-Admixture simulation.
#'
#'
#' @export
dbl_admixture <- function(p_anc, coanc_antepops, admix_proportions, paras = NULL, paras_only = FALSE, p_antepops_only = FALSE, geno_only = FALSE, ...) {
  # check if `p_anc` is valid
  check_p_anc(p_anc)
  # check if `coanc_antepops` is valid
  check_coancestry(coanc_antepops)

  # check if `paras` is valid
  stopifnot("paras should contain `Gamma` and `W`" = is.null(paras) || (!is.null(paras) && all(names(paras) %in% c("W", "Gamma"))))

  # estimate W and Gamma
  if (is.null(paras)) {
    paras <- est_paras_dbladmix(coanc_antepops = coanc_antepops, ...)
    if (paras_only) {
      return(paras)
    }
  }

  # sample allele frequencies among antecedent populations
  p_antepops <- bnpsd::draw_p_subpops(p_anc, diag(paras$Gamma)) %*% paras$W
  if (p_antepops_only) {
    return(p_antepops)
  }

  # check if `admix_proportions` is valid
  check_admix_proportions(admix_proportions)
  # check if the dimension of `coanc_antepops` and `admix_proportions` match
  stopifnot("`nrow(coanc_antepops)` and `nrow(admix_proportions)` should match" = nrow(admix_proportions) == nrow(coanc_antepops))

  # draw genotypes
  p_indiv <- p_antepops %*% admix_proportions
  p_indiv <- ifelse(p_indiv > 1, 1, ifelse(p_indiv < 0, 0, p_indiv))
  X <- bnpsd::draw_genotypes_admix(p_indiv)
  if (geno_only) {
    return(X)
  }

  # returns
  return(list(W = paras$W, Gamma = paras$Gamma, p_antepops = p_antepops, X = X))
}
