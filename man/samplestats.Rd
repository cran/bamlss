\name{samplestats}
\alias{samplestats}

\title{
  Sampling Statistics
}

\description{
  The function computes the average the log-likelihood, log-posterior, the deviance
  information criterion and estimated degrees of freedom from samples of, e.g., a
  \code{\link{bamlss}} object.
}

\usage{
samplestats(samples, x = NULL, y = NULL,
  family = NULL, logLik = FALSE, ...)
}

\arguments{
  \item{samples}{An object of class \code{"mcmc.list"} or \code{"bamlss"} which contains MCMC
    samples.}
  \item{x}{The \code{x} list as returned by function \code{\link{bamlss.frame}}.}
  \item{y}{The model response, as returned by function \code{\link{bamlss.frame}}}
  \item{family}{A \code{\link{bamlss.family}} object.}
  \item{logLik}{Logical, should the log-likelihood be computed, may take some time!}
  \item{\dots}{Currently not used.}
}

\details{
  If the log-likelihood is not available in the samples, the function tries to compute
  the information. Depending on the complexity of the model, this may take some time.
  Computations are based on the \code{$d()} or \code{$loglik()} function of the
  \code{\link{bamlss.family}} object.

  If a \code{\link{bamlss.family}} object contains a function \code{$p2d()} or \code{$p2loglik()},
  which computes the log-likelihood from parameters, these functions are used for computation.
}

\value{
  A list with the following entries (if available):
    \item{logLik}{The average log-likelihood.}
    \item{logPost}{The average log-posterior.}
    \item{DIC}{The deviance information criterion.}
    \item{pd}{The estimated degrees of freedom.}
}

\seealso{
  \code{\link{bamlss}}
}

\examples{
\dontrun{## Generate some data.
d <- GAMart()

## Estimate model without sampling statistics
b <- bamlss(num ~ s(x1) + s(x2) + s(x3) + te(lon,lat),
  data = d, samplestats = FALSE)

## Note: needs the $d() or $loglik() function in the family!
names(family(b))

## Compute sampling statistics.
samplestats(b)
}
}

\keyword{regression}

