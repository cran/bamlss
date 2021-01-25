\name{simSurv}
\alias{simSurv}
\alias{rSurvTime2}

\title{
  Simulate Survival Times
}

\description{
  Function \code{simSurv()} and \code{rSurvtime2()} simulate arbitrary (Cox-type) survival times
  using the inversion method. Function \code{simSurv()} is a simple wrapper that calls
  \code{rSurvtime2()}. The functions are based on the \R package \pkg{CoxFlexBoost} implementation
  \code{rSurvtime()} and only slightly modify the code.
}

\usage{
## Simulate a pre-specified survival times data set.
simSurv(n = 300)

## Simulate arbitrary survival times.
rSurvTime2(lambda, x, cens_fct, upper = 1000, ...,
  file = NULL, subdivisions = 1000)
}

\arguments{
  \item{n}{The number of individuals for which survival times should be simulated.}
  \item{lambda}{function. Baseline hazard \eqn{\lambda(t, x)}{lambda(t,x)}
    where time must be first argument. }
  \item{x}{matrix. (Sampled) values for covariates (without time). }
  \item{cens_fct}{ function. Function to compute (random) censoring. }
  \item{upper}{ upper boundary of the interval the random survival times fall into. }
  \item{\dots}{further arguments to be passed to \code{lambda} or \code{cens_fct}. }
  \item{file}{character. name of the data file the generated data set should be stored into
    (e.g., "survtimes.RData") or NULL if the dataset should directly be returned in \R. }
  \item{subdivisions}{The maximum number of subintervals for the integration.}
}

\details{
  This is basically a slight modification according the computation of the integral,
  see the manual page of function \code{rSurvtime()} of package \code{CoxFlexBoost} for
  details.
}

\value{
  A \code{data.frame} consisting of the observed survival time
  (\code{time}), the non-censoring indicator (\code{event}) and further
  covariates \code{x} is returned. If \code{file} is specified, the
  data.frame is additionally stored on the disc.
}

\references{
  Benjamin Hofner (2016). \pkg{CoxFlexBoost}: Boosting Flexible Cox Models
    (with Time-Varying Effects). R package version 0.7-0.

  Ralph Bender and Thomas Augustin and Maria Blettner (2005),
    Generating Survival Times to Simulate Cox Proportional Hazards Models.
    \emph{Statistics in Medicine}, \bold{24}, 1713-1723.
}

\seealso{
  \code{\link{cox_bamlss}}, \code{\link{opt_Cox}}, \code{\link{sam_Cox}}, \code{\link{bamlss}}
}

\examples{
## The following shows the code of the
## wrapper function simSurv().
set.seed(111)
n <- 100
X <- matrix(NA, nrow = n, ncol = 3)
X[, 1] <- runif(n, -1, 1)
X[, 2] <- runif(n, -3, 3)
X[, 3] <- runif(n, -1, 1)

## Specify censoring function.
cens_fct <- function(time, mean_cens) {
  ## Censoring times are independent exponentially distributed.
  censor_time <- rexp(n = length(time), rate = 1 / mean_cens)
  event <- (time <= censor_time)
  t_obs <- apply(cbind(time, censor_time), 1, min)
  ## Return matrix of observed survival times and event indicator.
  return(cbind(t_obs, event))
}

## log(time) is the baseline hazard.
lambda <-  function(time, x) {
  exp(log(time) + 0.7 * x[1] + sin(x[2]) + sin(time * 2) * x[3])
}

## Simulate data with lambda() and cens_fct().
d <- rSurvTime2(lambda, X, cens_fct, mean_cens = 5)
}

\keyword{survival}
\keyword{distribution}

