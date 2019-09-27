\name{smooth_check}
\alias{smooth_check}

\title{
  MCMC Based Simple Significance Check for Smooth Terms
}

\description{
  For each smooth term estimated with MCMC, the function computes 95% credible
  intervals and simply computes the fraction of the cases where the interval does not contain zero.
}
 
\usage{
smooth_check(object, newdata = NULL, model = NULL, term = NULL, ...)
}

\arguments{
  \item{object}{A fitted model object which contains MCMC samples.}
  \item{newdata}{Optionally, use new data for computing the check.}
  \item{model}{Character, for which model should the check be computed?}
  \item{term}{Character, for which term should the check be computed?}
  \item{\dots}{Arguments passed to \code{\link{predict.bamlss}}.}
}

\examples{
\dontrun{## Simulate some data.
d <- GAMart()

## Model formula.
f <- list(
  num ~ s(x1) + s(x2) + s(x3),
  sigma ~ s(x1) + s(x2) + s(x3)
)

## Estimate model with MCMC.
b <- bamlss(f, data = d)

## Run the check, note that all variables
## for sigma should have no effect.
smooth_check(b)
}
}

\keyword{regression}
