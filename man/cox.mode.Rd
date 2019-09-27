\name{cox_mode}
\alias{cox_mode}

\title{
  Cox Model Posterior Mode Estimation
}

\description{
  This function computes posterior mode estimates of the parameters of a flexible Cox model
  with structured additive predictors using a Newton-Raphson algorithm. Integrals are solved
  numerically. Moreover, optimum smoothing variances are computed using a stepwise optimization,
  see also the details section of function \code{\link{bfit}}.
}

\usage{
cox_mode(x, y, start, weights, offset,
  criterion = c("AICc", "BIC", "AIC"),
  nu = 0.1, update.nu = TRUE,
  eps = .Machine$double.eps^0.25, maxit = 400,
  verbose = TRUE, digits = 4, ...)
}

\arguments{
 \item{x}{The \code{x} list, as returned from function
    \code{\link{bamlss.frame}} and transformed by function \code{\link{surv_transform}},
    holding all model matrices and other information that is used for
    fitting the model.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{weights}{Prior weights on the data, as returned from function \code{\link{bamlss.frame}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{criterion}{Set the information criterion that should be used, e.g., for smoothing
    variance selection. Options are the corrected AIC \code{"AICc"}, the \code{"BIC"} and
    \code{"AIC"}.}
  \item{nu}{Calibrates the step length of parameter updates of one Newton-Raphson update.}
  \item{update.nu}{Should the updating step length be optimized in each iteration
    of the backfitting algorithm.}
  \item{eps}{The relative convergence tolerance of the backfitting algorithm.}
  \item{maxit}{The maximum number of iterations for the backfitting algorithm}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{\dots}{Currently not used.}
}

\value{
  A list containing the following objects:
  \item{fitted.values}{A named list of the fitted values of the modeled parameters
    of the selected distribution.}
  \item{parameters}{The estimated set regression coefficients and smoothing variances.}
  \item{edf}{The equivalent degrees of freedom used to fit the model.}
  \item{logLik}{The value of the log-likelihood.}
  \item{logPost}{The value of the log-posterior.}
  \item{hessian}{The Hessian matrix evaluated at the posterior mode.}
  \item{converged}{Logical, indicating convergence of the backfitting algorithm.}
  \item{time}{The runtime of the algorithm.}
}

\references{
  Umlauf N, Klein N, Zeileis A (2016). Bayesian Additive Models for Location
    Scale and Shape (and Beyond). \emph{(to appear)}
}

\seealso{
  \code{\link{cox_mcmc}}, \code{\link{cox_bamlss}}, \code{\link{surv_transform}},
  \code{\link{simSurv}}, \code{\link{bamlss}}
}

\examples{
## Please see the examples of function cox_mcmc()!
}

\keyword{regression}
\keyword{survival}

