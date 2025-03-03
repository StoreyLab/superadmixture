% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbl_admixture.R
\name{est_paras_dbladmix}
\alias{est_paras_dbladmix}
\title{Estimating parameters used in the double-admixture simulation}
\usage{
est_paras_dbladmix(
  coanc_antepops,
  S = 2 * nrow(coanc_antepops),
  tau1 = 1.1,
  tau2 = 1.1,
  max_iters = 1e+05,
  tol = 1e-06,
  verbose = FALSE
)
}
\arguments{
\item{coanc_antepops}{A \eqn{K \times K} matrix representing coancestry among antecedent populations.}

\item{S}{An integer representing the number of independent Balding-Nichols distributions.}

\item{tau1, tau2}{Scalars controlling the step size of the optimization iteration.
They should be greater than 1.0 to ensure the convergence (defaults to 1.1).}

\item{max_iters}{Maximum number of iterations. \code{max_iters} should be set to a large number (100,000 by default) to ensure convergence.}

\item{tol}{A Scalar for convergence criterion. If the relative change of
\eqn{\|\boldsymbol{\Lambda} - \boldsymbol{W}^{\prime}\boldsymbol{\Gamma}\boldsymbol{W}\|_F} is less than \code{tol},
the algorithm halts.}

\item{verbose}{A Boolean that controls the returned values. If \code{verbose = F} (default),
only the estimated \eqn{\boldsymbol{W}} and \eqn{\boldsymbol{\Gamma}}  are returned. If \code{T}, this function returns
a list consisting of the estimated \eqn{\boldsymbol{W}}, the estimated \eqn{\boldsymbol{\Gamma}} and a list of errors over iterations.}
}
\value{
\itemize{
\item If \code{verbose = T}, a list with the following elements is returned:
\itemize{
\item \code{W}: The estimated \eqn{\boldsymbol{W}}
\item \code{Gamma}: The estimated \eqn{\boldsymbol{\Gamma}}
\item \code{err_list}: A list of \eqn{\|\boldsymbol{\Lambda} - \boldsymbol{W}^{\prime}\boldsymbol{\Gamma}\boldsymbol{W}\|_F}'s over iterations.
}
\item If \code{verbose = F}, a list containing the following elements is returned:
\itemize{
\item \code{W}: The estimated \eqn{\boldsymbol{W}}
\item \code{Gamma}: The estimated \eqn{\boldsymbol{\Gamma}}
}
}
}
\description{
\code{est_paras_dbladmix()} identifies a \eqn{K \times S} matrix \eqn{\boldsymbol{W}} and a diagonal \eqn{S \times S}
matrix \eqn{\boldsymbol{\Gamma}} such that \eqn{\boldsymbol{\Lambda} \approx \boldsymbol{W}^{\prime}\boldsymbol{\Gamma}\boldsymbol{W}}.
}
\details{
This function solves the following optimization problem based on Proximal Alternating Linearized Minimization (PALM),
and it is guaranteed to converge to a stationary point:
\deqn{\arg\min \|\boldsymbol{\Lambda} - \boldsymbol{W}^{\prime}\boldsymbol{\Gamma}\boldsymbol{W}\|_F}
subjects to those constraints:
\enumerate{
\item All entries of \eqn{\boldsymbol{W}} are greater or equal to 0
\item \eqn{\sum_{s = 1}^S w_{su} = 1} for all \eqn{u}
\item \eqn{0.01 \leq \gamma_{ss} \leq 0.99} for all \eqn{s}
\item \eqn{\gamma_{st} = 0} for all \eqn{s \neq t}
}
}
\examples{
# estimate parameters  --------------------------------------------------------------------
data("coanc_pops_amr",  package = "superadmixture")
paras <- est_paras_dbladmix(coanc_pops_amr, verbose = FALSE)

}
\keyword{internal}
