\name{rb}
\alias{rb}
\alias{smooth.construct.randombits.smooth.spec}

\title{Random Bits for BAMLSS}

\description{
  This smooth constructor implements random bits model terms. Note that this is experimental.
}

\usage{
## Linear smooth constructor.
rb(..., k = 50)

## For mgcv.
\method{smooth.construct}{randombits.smooth.spec}(object, data, knots, ...)
}

\arguments{
  \item{\dots}{For function \code{rb()} a formula of the type \code{~x1+x2+x3} that specifies
    the covariates that should be modeled.}
  \item{k}{Integer, number of radnom bit columns in the design matrix.}
  \item{object, data, knots}{See \code{\link[mgcv]{smooth.construct}}.}
}

\value{
  Function \code{rb()}, similar to function \code{\link[mgcv]{s}} a simple smooth specification
  object.
}


\seealso{
\code{\link{bamlss}}, \code{\link{predict.bamlss}}, \code{\link{bfit}}, \code{\link{boost}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart()

## Estimate model.
f <- num ~ rb(x1) + rb(x2) + rb(x3) + rb(~lon+lat)

b <- bamlss(f, data = d)

plot(b)
}
}

\keyword{regression}

