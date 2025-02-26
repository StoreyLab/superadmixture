g_cubature <- function(rho, au, bu, av, bv) {
  bivnorm_corrmat <- diag(2)
  bivnorm_corrmat[1,2] <- bivnorm_corrmat[2, 1] <- rho

  # do numeric integral
  integrand <- function(z, au, bu, av, bv, bivnorm_corrmat) {
    zu <- z[1]
    zv <- z[2]
    pu <- stats::qbeta(stats::pnorm(zu), au, bu)
    pv <- stats::qbeta(stats::pnorm(zv), av, bv)
    pu * pv * mvtnorm::dmvnorm(z, sigma = bivnorm_corrmat)
  }

  return(cubature::cuhre(
    f = integrand, lowerLimit = rep(-5, 2), upperLimit = rep(5, 2),
    relTol = 1e-3, absTol = 1e-6,
    au = au, bu = bu, av = av, bv = bv, bivnorm_corrmat = bivnorm_corrmat)$integral)
}


#g_monteCarlo <- function(rho, au, bu, av, bv, n_mc = 10000) {
#  bivnorm_corrmat <- diag(2)
#  bivnorm_corrmat[1,2] <- bivnorm_corrmat[2, 1] <- rho
#  z <- mvtnorm::rmvnorm(n_mc, sigma = bivnorm_corrmat)
#  zu <- z[, 1]
#  zv <- z[, 2]
#  pu <- stats::qbeta(stats::pnorm(zu), au, bu)
#  pv <- stats::qbeta(stats::pnorm(zv), av, bv)
#  g <- mean(zu * zv)
#  gprime <- mean(pu * pv * (rho / (1 - rho^2) + (zu * zv - rho * (zu^2 - rho * zu * zv + zv^2)) / (1 - rho^2)^2))
#  return(list(g = g, gprime = gprime))
#}


gprime_cubature <- function(rho, au, bu, av, bv) {
  bivnorm_corrmat <- diag(2)
  bivnorm_corrmat[1,2] <- bivnorm_corrmat[2, 1] <- rho

  # do numeric integral
  integrand <- function(z, au, bu, av, bv, bivnorm_corrmat) {
    zu <- z[1]
    zv <- z[2]
    pu <- stats::qbeta(stats::pnorm(zu), au, bu)
    pv <- stats::qbeta(stats::pnorm(zv), av, bv)
    rho <- bivnorm_corrmat[1, 2]
    pu * pv * mvtnorm::dmvnorm(z, sigma = bivnorm_corrmat) *
      (rho / (1 - rho^2) + (zu * zv - rho * (zu^2 - rho * zu * zv + zv^2)) / (1 - rho^2)^2)
  }

  return(cubature::cuhre(
    f = integrand, lowerLimit = rep(-5, 2), upperLimit = rep(5, 2),
    relTol = 1e-3, absTol = 1e-6,
    au = au, bu = bu, av = av, bv = bv, bivnorm_corrmat = bivnorm_corrmat)$integral)
}


#gprime_monteCarlo <- function(rho, au, bu, av, bv, n_mc = 10000) {
#  bivnorm_corrmat <- diag(2)
#  bivnorm_corrmat[1,2] <- bivnorm_corrmat[2, 1] <- rho
#  z <- mvtnorm::rmvnorm(n_mc, sigma = bivnorm_corrmat)
#  zu <- z[, 1]
#  zv <- z[, 2]
#  pu <- stats::qbeta(stats::pnorm(zu), au, bu)
#  pv <- stats::qbeta(stats::pnorm(zv), av, bv)
#  return(mean(pu * pv * (rho / (1 - rho^2) + (zu * zv - rho * (zu^2 - rho * zu * zv + zv^2)) / (1 - rho^2)^2)))
#}


estimate_rho <- function(p_anc, coanc_u, coanc_v, coanc_uv, max_iters = 100, tol = 1e-4, verbose = FALSE) {
  au <- (1 - coanc_u) / coanc_u * p_anc
  bu <- (1 - coanc_u) / coanc_u * (1 - p_anc)
  av <- (1 - coanc_v) / coanc_v * p_anc
  bv <- (1 - coanc_v) / coanc_v * (1 - p_anc)

  # initialize rho
  rho <- coanc_uv / sqrt(coanc_u * coanc_v)
  if (rho > 0.99) rho <- 0.99
  if (rho < 0.01) rho <- 0.01
  # initialize optimization
  converged <- FALSE
  if (verbose) print("Estimating rho....")

  # run Newton-Raphson
  for (i in 1:max_iters) {
    #  if (method == "numeric") {
    g <- g_cubature(rho, au, bu, av, bv) - p_anc^2 - p_anc * (1 - p_anc) * coanc_uv
    gprime <- gprime_cubature(rho, au, bu, av, bv)
   #  } else {
    #   obj <- g_monteCarlo(rho, au, bu, av, bv, n_mc = n_mc)
    #   g <- obj$g -  - p_anc^2 - p_anc * (1 - p_anc) * coanc_uv
    #  gprime <- obj$gprime
    #gprime <- gprime_monteCarlo(rho, au, bu, av, bv, n_mc = n_mc)
    #}

    # update
    rho_next <- rho - g / gprime
    if (rho_next < 0.01) rho_next <- 0.01
    if (rho_next > 0.99) rho_next <- 0.99
    if (verbose) print(paste0("Current g = ", round(g, digits = 6), "; current rho = ", round(rho, digits = 3), "; current rho_next = ", round(rho_next, digits = 3)))
    converged <- abs(g) < tol || rho == rho_next
    if (converged) break
    rho <- rho_next
  }
  return(rho)
}

#' Estimate the correlation matrix \eqn{\boldsymbol{\Sigma}_{\boldsymbol{z}}} used in the NORTA method
#'
#' @description
#' `estimate_corrmat()` adopts the Newton-Raphson method to estimate the the correlation matrix \eqn{\boldsymbol{\Sigma}_{\boldsymbol{z}}} used in [norta].
#'
#' @details
#' This function finds the \eqn{\boldsymbol{\Sigma}_{\boldsymbol{z}}} used in [norta]. This problem can be reduced to \eqn{K \times (K - 1) / 2} independent problems:
#' for each \eqn{u < v}, we need to find such \eqn{\rho} that
#' \deqn{\mathbb{E}[X_u X_v] = a^2 + a(1 - a)\lambda_{uv}}
#' We note that there is no close form soluation for \eqn{\rho}. We adopt the numeric integration to perform a numerical search for \eqn{\rho}.
#' Define
#' \deqn{g(\rho) = \mathbb{E}[F_u^{-1}(\Phi(z_u))F_v^{-1}(\Phi(z_v))] - a^2 - a(1 - a)\lambda_{uv}}
#' \deqn{g^{\prime}(\rho) = \mathbb{E}[F_u^{-1}(\Phi(z_u))F_v^{-1}(\Phi(z_v))\left(\frac{\rho}{1 - \rho^2} + \frac{z_u z_v - \rho(z_u^2 - \rho z_u z_v + z_v^2)}{(1 - \rho^2)^2}\right)]}
#' where \eqn{\Phi} denotes the univariate standard normal cumulative density function (cdf); \eqn{F_u} denotes the Balding-Nichols(\eqn{a_i}, \eqn{\lambda_{uu}});
#' \eqn{F_u^{-1}} denotes the inverse cdf of \eqn{F_u}.The iteration of Newton's method is
#' \deqn{\rho \leftarrow \rho \leftarrow \frac{g(\rho)}{g^{\prime}(\rho)}}
#'
#' @param p_anc A scalar between 0 and 1.
#' @inheritParams dbl_admixture
#' @param max_iters Maximum number of iterations.
#' @param tol A scalar for stopping criteria. If \eqn{|g(\rho)|} is less than `tol`, the Newton's method halts.
#'
#' @returns
#' A \eqn{K \times K} matrix \eqn{\boldsymbol{\Sigma}_{\boldsymbol{z}}} of the correlation matrix standard multivariate Normal distribution that will be used by [norta()].
#'
#' @seealso
#' [norta] for details about simulating allele frequencies and genotypes.
#'
#'
#' @export
#' @keywords internal
estimate_corrmat <- function(p_anc, coanc_antepops, max_iters = 100, tol = 1e-4, verbose = FALSE) {
  k_antepops <- nrow(coanc_antepops)
  corrmat <- diag(1, k_antepops, k_antepops)
  for (u in 1:(k_antepops - 1)) {
    for (v in (u + 1):k_antepops) {
      coanc_u  <- coanc_antepops[u, u]
      coanc_v  <- coanc_antepops[v, v]
      coanc_uv <- coanc_antepops[u, v]
      sigma_uv <- estimate_rho(p_anc, coanc_u, coanc_v, coanc_uv, max_iters, tol, verbose)
      corrmat[u, v] <- corrmat[v, u] <- sigma_uv
    }
  }
  return(Matrix::nearPD(corrmat)$mat)
}


#' Simulate allele frequencies and genotypes using the NORTA method
#'
#' @description
#' `norta()` generates random samples of antecedent allele frequencies \eqn{\boldsymbol{P}} and random samples of genotypes \eqn{\boldsymbol{X}}
#' according to NORmal-To-Anything method.
#'
#' @details
#' Here we consider the problem of simulating a \eqn{K\times 1} dimensional random vector \eqn{\boldsymbol{p}} given a scalar \eqn{a} and a \eqn{K \times K} matrix
#' \eqn{\boldsymbol{\Lambda}} that satisfies the following moment constraints:
#' \deqn{\mathbb{E}[p_{u}] = a \quad \forall u \in 1, \ldots, K}
#' \deqn{\text{Cov}(p_u, p_v) = a(1 - a)\lambda_{uv} \quad \forall (u,v) \in 1, \ldots, K}
#' Let \eqn{\Phi} denote the univariate standard normal cumulative density function (cdf). Let \eqn{F_u} denote the Balding-Nichols(\eqn{a}, \eqn{\lambda_{uu}}).
#' Let \eqn{F_u^{-1}} denotes the inverse cdf of \eqn{F_u}. Based on the NORmal-To-Anything method, \eqn{\boldsymbol{p}} can be generated in two steps.
#' \enumerate{
#'     \item Simulate a length-\eqn{K} random vector \eqn{\boldsymbol{z}} from the standard multivariate Normal distribution \eqn{\Sigma_{\boldsymbol{z}}}.
#'     \item Let \eqn{p_u} be \eqn{p_u = F_u^{-1}(\Phi(z_u))} for \eqn{u = 1, \ldots, K}.
#' }
#' The method for identifying  \eqn{\Sigma_{\boldsymbol{z}}} is described in [estimate_corrmat]. To simulate a matrix of antecedent populations allele frequencies,
#' we adopt the strategy described above to each locus and simulate a random vector \eqn{\boldsymbol{p}_i} at each locus given \eqn{a_i} and  \eqn{\boldsymbol{\Lambda}}.
#
#' @inheritParams dbl_admixture
#' @param p_antepops_only,geno_only Booleans used to control the returned values. If `p_antepops_only = T`, only the \eqn{m \times K} matrix of
#' antecedent allele frequencies \eqn{\boldsymbol{P}} is returned. If `geno_only = T`, only the \eqn{m \times n} matrix of genotypes \eqn{\boldsymbol{X}} is returned.
#' Otherwise, a list containing \eqn{\boldsymbol{P}} and \eqn{\boldsymbol{X}} is returned.
#' @param parallel A boolean indicating whether `mclapply` is used.
#' @param max_iters Maximum number of iterations for Newton-Raphson method in [estimate_corrmat].
#' @param tol A scalar for stopping criteria of Newton-Raphson method.
#' @param mc_cores An optional integer of number of cores for `mclapply`. If not available, `mc_cores` will be determined by `parallel::detectCores()`.
#'
#' @returns
#' + If `p_antepops = T`, an \eqn{m \times K} matrix of allele frequencies of antecedent populations \eqn{\boldsymbol{P}} is returned.
#' + If `geno_only = T`, an \eqn{m \times n} matrix of genotypes is returned.
#' + Otherwise, a named list is returned, containing:
#'    + `p_antepops`: the simulated \eqn{m \times K} matrix of antecedent allele frequencies,
#'    + `X`: the simulated \eqn{m \times n} matrix of genotypes.
#'
#' @examples
#' ## Generate a genotype from HGDP -------------------------------------------------------------------
#' ## It takes several hours to run the following codes.
#' ## For a quick implementation of NORTA method, see `norta_approx` function.
#' \dontrun{
#' data("fam_hgdp",         package = "superadmixture")
#' data("p_anc_hgdp",       package = "superadmixture")
#' data("admix_props_hgdp", package = "superadmixture")
#' data("coanc_pops_hgdp",  package = "superadmixture")
#' X_hgdp <- norta(
#'     p_anc = p_anc_hgdp,
#'     coanc_antepops = coanc_pops_hgdp,
#'     admix_proportions = admix_props_hgdp,
#'     parallel = TRUE,
#'     geno_only = TRUE)
#'
#' ## Estimate kinship from the simulated genotypes
#' kinship <- popkin::popkin(X_hgdp)
#' coanc_indiv <- popkin::inbr_diag(kinship)
#' coanc_indiv <- ifelse(coanc_indiv < 0, 0, coanc_indiv)
#'
#' ## Visualize kinship
#' popkin::plot_popkin(coanc_indiv, labs = fam_hgdp$subpop, labs_las = 2)}
#'
#' @seealso
#' [estimate_corrmat] for details about finding the correlation matrix \eqn{\Sigma_{\boldsymbol{z}}}.
#'
#'
#' @export
norta <- function(p_anc, coanc_antepops, admix_proportions, p_antepops_only = FALSE, geno_only = FALSE, parallel = FALSE,
                  max_iters = 100, tol = 1e-4, mc_cores = parallel::detectCores()) {
  check_coancestry(coanc_antepops)
  coanc_antepops <- as.matrix(Matrix::nearPD(coanc_antepops)$mat)

  # a helper function for drawing a length-K random vector of allele frequencies
  draw_p_antepops <- function(p_anc, coanc_antepops, max_iters = 100, tol = 1e-4) {
      if (p_anc == 0) return(rep(0, k_antepops))
      if (p_anc == 1) return(rep(1, k_antepops))
      k_antepops <- nrow(coanc_antepops)
      corrmat <- estimate_corrmat(p_anc, coanc_antepops, max_iters, tol)
      z <- MASS::mvrnorm(n = 1, mu = rep(0, k_antepops), Sigma = corrmat)
      a <- (1 - diag(coanc_antepops)) / diag(coanc_antepops) * p_anc
      b <- (1 - diag(coanc_antepops)) / diag(coanc_antepops) * (1 - p_anc)
      return(stats::qbeta(stats::pnorm(z), a, b))
  }

  # draw allele frequencies
  if (parallel) {
    p_antepops <- parallel::mclapply(p_anc, draw_p_antepops,
                         coanc_antepops=coanc_antepops,
                         max_iters=max_iters,
                         tol=tol,
                         mc.cores = mc_cores)
    p_antepops <- purrr::reduce(p_antepops, rbind)
    } else {
    p_antepops <- sapply(p_anc, draw_p_antepops,
                         coanc_antepops=coanc_antepops,
                         max_iters=max_iters,
                         tol=tol)
    p_antepops <- t(p_antepops)
  }

  if (p_antepops_only) {
    return(p_antepops)
  }

  check_admix_proportions(admix_proportions)
  # define `k_antepops` and `n_inds`
  if (ncol(admix_proportions) == nrow(coanc_antepops)) {
    admix_proportions <- t(admix_proportions)
  }
  k_antepops <- nrow(admix_proportions)
  n_inds <- ncol(admix_proportions)

  # draw genotypes
  p_indiv <- p_antepops %*% admix_proportions
  p_indiv <- ifelse(p_indiv > 1, 1, ifelse(p_indiv < 0, 0, p_indiv))
  X <- bnpsd::draw_genotypes_admix(p_indiv)
  if (geno_only) {
    return(X)
  }

  # returns
  return(list(p_antepops = p_antepops, X = X))

}


#' An accelerated NORTA implementation for simulating allele frequencies and genotypes
#'
#' @description
#' `norta()` generates random samples of antecedent allele frequencies \eqn{\boldsymbol{P}} and random samples of genotypes \eqn{\boldsymbol{X}}
#' according to NORmal-To-Anything method. Instead of estimating \eqn{\boldsymbol{\Sigma}_{\boldsymbol{z}}} for each locus, this method only estimates
#' \eqn{\boldsymbol{\Sigma}_{\boldsymbol{z}}} for \eqn{a_i = 0.01, 0.02, \ldots, 0.99}. See [norta] for details about NORTA simulation.
#'
#' @inheritParams norta
#'
#' @returns
#' + If `p_antepops = T`, an \eqn{m \times K} matrix of allele frequencies of antecedent populations \eqn{\boldsymbol{P}} is returned.
#' + If `geno_only = T`, an \eqn{m \times n} matrix of genotypes is returned.
#' + Otherwise, a named list is returned, containing:
#'    + `p_antepops`: the simulated \eqn{m \times K} matrix of antecedent allele frequencies,
#'    + `X`: the simulated \eqn{m \times n} matrix of genotypes.
#'
#' @examples
#' ## Generate a genotype from HGDP -------------------------------------------------------------------
#' ## It takes 230.851 sec to finish the following codes on a 2.3 GHz 18-Core Intel Xeon W iMac Pro.
#' \dontrun{
#' data("fam_hgdp",         package = "superadmixture")
#' data("p_anc_hgdp",       package = "superadmixture")
#' data("admix_props_hgdp", package = "superadmixture")
#' data("coanc_pops_hgdp",  package = "superadmixture")
#' X_hgdp <- norta_approx(
#'       p_anc = p_anc_hgdp,
#'       coanc_antepops = coanc_pops_hgdp,
#'       admix_proportions = admix_props_hgdp,
#'       parallel = TRUE,
#'       geno_only = TRUE)
#'
#' ## Estimate kinship from the simulated genotypes
#' kinship <- popkin::popkin(X_hgdp)
#' coanc_indiv <- popkin::inbr_diag(kinship)
#' coanc_indiv <- ifelse(coanc_indiv < 0, 0, coanc_indiv)
#'
#' ## Visualize kinship
#' popkin::plot_popkin(coanc_indiv, labs = fam_hgdp$subpop, labs_las = 2)}
#'
#' @seealso
#' [norta] for details about finding the correlation matrix \eqn{\Sigma_{\boldsymbol{z}}}.
#'
#'
#' @export
norta_approx <- function(p_anc, coanc_antepops, admix_proportions, p_antepops_only = FALSE, geno_only = FALSE, parallel = FALSE,
                         max_iters = 100, tol = 1e-4, mc_cores = parallel::detectCores()) {
  check_coancestry(coanc_antepops)
  k_antepops <- nrow(coanc_antepops)
  coanc_antepops <- as.matrix(Matrix::nearPD(coanc_antepops)$mat)

  # draw allele frequencies
  if (parallel) {
    corrmat_list <- parallel::mclapply(as.list(seq(0.01, 0.99, length.out = 99)), function(a) estimate_corrmat(a,
                                       coanc_antepops=coanc_antepops,
                                       max_iters=max_iters,
                                       tol=tol),
                                       mc.cores = mc_cores)
  } else {
    corrmat_list <- lapply(as.list(seq(0.01, 0.99, length.out = 99)), function(a) estimate_corrmat(a,
                           coanc_antepops=coanc_antepops,
                           max_iters=max_iters,
                           tol=tol))
  }

  k_antepops <- nrow(coanc_antepops)
  p_antepops <- sapply(p_anc, function(p_anc_one) {
    p_anc_one <- round(p_anc_one, digits = 2)
    if (p_anc_one == 0) return(rep(0, k_antepops))
    if (p_anc_one == 1) return(rep(1, k_antepops))
    corrmat <- corrmat_list[[p_anc_one * 100]]
    z <- MASS::mvrnorm(n = 1, mu = rep(0, k_antepops), Sigma = corrmat)
    a <- (1 - diag(coanc_antepops)) / diag(coanc_antepops) * p_anc_one
    b <- (1 - diag(coanc_antepops)) / diag(coanc_antepops) * (1 - p_anc_one)
    return(stats::qbeta(stats::pnorm(z), a, b))
  })

  p_antepops <- t(p_antepops)

  if (p_antepops_only) {
    return(p_antepops)
  }

  check_admix_proportions(admix_proportions)
  # define `k_antepops` and `n_inds`
  if (ncol(admix_proportions) == nrow(coanc_antepops)) {
    admix_proportions <- t(admix_proportions)
  }
  n_inds <- ncol(admix_proportions)

  # draw genotypes
  p_indiv <- p_antepops %*% admix_proportions
  p_indiv <- ifelse(p_indiv > 1, 1, ifelse(p_indiv < 0, 0, p_indiv))
  X <- bnpsd::draw_genotypes_admix(p_indiv)
  if (geno_only) {
    return(X)
  }

  # returns
  return(list(p_antepops = p_antepops, X = X))
}

