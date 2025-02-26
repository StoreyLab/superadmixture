#' Factor Individual specific Allele Frequencies
#'
#' @description
#' `factor_p_indiv()` factors the individual specific allele frequencies \eqn{\boldsymbol{\Pi}} into two matrices \eqn{\boldsymbol{P}} and \eqn{\boldsymbol{Q}} such that
#' \eqn{\boldsymbol{\Pi} = \boldsymbol{P}\boldsymbol{Q}}.
#'
#' @details
#' `factor_p_indiv()` finds matrices \eqn{\boldsymbol{P}} and \eqn{\boldsymbol{Q}} such that \eqn{\boldsymbol{\Pi} = \boldsymbol{P}\boldsymbol{Q}}, by solving the following
#' optimization problem
#' \deqn{\arg\min_{\boldsymbol{P}, \boldsymbol{Q}} \|\boldsymbol{\Pi} - \boldsymbol{P}\boldsymbol{Q}\|_F}
#' where  \eqn{\boldsymbol{P}} and \eqn{\boldsymbol{Q}} obeys these conditions:
#' \enumerate{
#'  \item \eqn{p_{iu} \in [0,1]  ~ \forall (i,u)}
#'  \item \eqn{q_{uj} \ge 0  ~ \forall (u,j)} and \eqn{\sum_u q_{uj} = 1 ~ \forall j}
#' }
#'
#' @param p_indiv  A matrix of individual specific allele frequencies \eqn{\boldsymbol{\Pi}}. Each \eqn{\pi_{ij}} should be an numeric matrix each element lies between 0 and 1.
#' @inheritParams est_p_indiv
#' @param rowspace An optional \eqn{n \times K} dimensional rowspace of \eqn{\boldsymbol{Q}}.
#' @param init_method One of "random" or "rowspace" (default). If `init_method` is "random",
#' the columns of \eqn{\boldsymbol{Q}} are drawn from Dirichlet distribution. If `init_method` is "rowspace",
#' \eqn{\boldsymbol{Q}} are initialized by projecting `rowspace` to the simplex.
#' @param max_iters Maximum number of iterations.
#' @param tol A scalar for convergence criterion. If the relative change of \eqn{\|\boldsymbol{\Pi} - \boldsymbol{P}\boldsymbol{Q}\|_F} is less than `tol`,
#' the algorithm halts.
#' @param verbose  A Boolean that controls the returned values. If `verbose = F` (default), only the estimated \eqn{\boldsymbol{P}} and \eqn{\boldsymbol{Q}}
#' are returned. If `T`, this function returns a list consisting of the estimated \eqn{\boldsymbol{P}},
#' the estimated \eqn{\boldsymbol{Q}} and a list of errors over iterations.
#'
#' @returns
#' + If `verbose = T`, a named list is returned, containing:
#'    + `P_hat`: The \eqn{m \times K} matrix of estimated \eqn{\boldsymbol{P}},
#'    + `Q_hat`: The \eqn{K \times n} matrix of estimated \eqn{\boldsymbol{Q}},
#'    + `err_list`: A list of \eqn{\|\boldsymbol{\Pi} - \boldsymbol{P}\boldsymbol{Q}\|_F}'s over iterations.
#' + Otherwise, a named list containing `P_hat` and `Q_hat` is returned.
#'
#' @examples
#' ## Estimate Admixture Proportions of HGDP  ----------------------------------------
#' ## It takes 14.469 sec to finish the following codes on a 2.3 GHz 18-Core Intel Xeon W iMac Pro.
#' data("X_hgdp", package = "superadmixture")
#' obj <- est_p_indiv(X_hgdp, k_antepops = 7, loci_on_cols = TRUE)
#' rowspace <- obj$rowspace
#' p_indiv <- obj$p_indiv
#' admix_props <- factor_p_indiv(p_indiv, k_antepops = 7,
#'           rowspace = rowspace, tol = 1e-2, max_iters = 200)$Q_hat
#' @export
factor_p_indiv <- function(p_indiv, k_antepops, init_method = "rowspace", rowspace = NULL, loci_on_cols = FALSE, max_iters = 1000000, tol = 1e-4, verbose = FALSE) {
  # check whether the parameter `init_method` is valid
  stopifnot('`init_method` %in% c("random", "rowspace")' = init_method %in% c("random", "rowspace"))
  stopifnot('`init_method` = "rowspace" && `!is.null(rowspace)`' = (init_method == "rowspace" && !is.null(rowspace)) || init_method == "random")

  if (loci_on_cols) p_indiv <- t(p_indiv)
  m_loci <- nrow(p_indiv)
  n_inds <- ncol(p_indiv)

  stopifnot('`init_method` = "rowspace" && `nrow(rowspace) == n_inds`' = (init_method == "rowspace" && nrow(rowspace) == n_inds) || init_method == "random")
  stopifnot('`init_method` = "rowspace" && `ncol(rowspace) == k_antepops`' = (init_method == "rowspace" && ncol(rowspace) == k_antepops) || init_method == "random")

  # initialization
  if (init_method == "random") Q <- t(rdirichlet(n = n_inds, alpha = rep(1, k_antepops)))
  if (init_method == "rowspace") {
    Q <- t(rowspace)
    for (j in 1:n_inds) {
      Q[, j] <- projsplx(Q[, j])
    }
  }
  P <- p_indiv %*% t(Q) %*% MASS::ginv(tcrossprod(Q))
  P <- ifelse(P < 0, 0, ifelse(P > 1, 1, P))

  # run optimization to estimate parameters for dbl-admixture
  converged <- FALSE
  first_order_optim <- FALSE
  previous_err <- norm(p_indiv - P %*% Q, "F")
  err_list <- list()

  for (i in 1:max_iters) {
    if (!first_order_optim) {
      # use newton method for updating
      Q <- MASS::ginv(crossprod(P)) %*% t(P) %*% p_indiv
      for (j in 1:n_inds) {
        Q[, j] <- projsplx(Q[, j])
      }
      P <- p_indiv %*% t(Q) %*% MASS::ginv(tcrossprod(Q))
      P <- ifelse(P < 0, 0, ifelse(P > 1, 1, P))
    } else {
      # use first-order gradient descent for updating
      Q <- Q - 1 / norm(crossprod(P), "F") * t(P) %*% (P %*% Q - p_indiv)
      for (j in 1:n_inds) {
        Q[, j] <- projsplx(Q[, j])
      }
      P <- P - 1 / norm(tcrossprod(Q), "F") * (P %*% Q - p_indiv) %*% t(Q)
      P <- ifelse(P < 0, 0, ifelse(P > 1, 1, P))
    }

    # check convergence
    current_err <- norm(p_indiv - P %*% Q, "F")
    if (verbose) {
      if (first_order_optim) {
        print(paste0("Iter = ", i, ", current_error = ", round(current_err, digits = 3), " (first order iteration)"))
      } else {
        print(paste0("Iter = ", i, ", current_error = ", round(current_err, digits = 3), " (Newton iteration)"))
      }
    }

    err_list[[i]] <- current_err
    converged <- first_order_optim && (previous_err - current_err) / previous_err < tol
    first_order_optim <- first_order_optim || ((previous_err - current_err) / previous_err < 1e-2)
    previous_err <- current_err
    if (converged) break
  }

  if (verbose) {
    return(list(P_hat = P, Q_hat = Q, err_list = err_list))
  } else {
    return(list(P_hat = P, Q_hat = Q))
  }
}
