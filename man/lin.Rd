\name{lin}
\alias{lin}
\alias{smooth.construct.linear.smooth.spec}

\title{Linear Effects for BAMLSS}

\description{
  This smooth constructor implements simple linear effects. The columns of the design matrix
  are automatically scaled. The main advantage of this constructor is speed when used in
  the BAMLSS boosting algorithm \code{\link{boost}}. Optionally, a ridge penalty can be
  added, please see the example.
}

\usage{
## Linear smooth constructor.
lin(...)

## For mgcv.
\method{smooth.construct}{linear.smooth.spec}(object, data, knots, ...)
}

\arguments{
  \item{\dots}{For function \code{lin()} a formula of the type \code{~x1+x2+x3} that specifies
    the covariates that should be modeled.}
  \item{object, data, knots}{See \code{\link[mgcv]{smooth.construct}}.}
}

\value{
  Function \code{lin()}, similar to function \code{\link[mgcv]{s}} a simple smooth specification
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
f <- num ~ lin(~x1+x2+x3+fac,ridge=TRUE)

b <- bamlss(f, data = d)

summary(b)
}
}

\keyword{regression}

