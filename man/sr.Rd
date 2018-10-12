\name{smooth.construct.sr.smooth.spec}
\alias{smooth.construct.sr.smooth.spec}

\title{Random Effects P-Spline}

\description{
  This smooth constructor implements the random effects representation of a P-spline.
}

\usage{
\method{smooth.construct}{sr.smooth.spec}(object, data, knots, ...)
}

\arguments{
  \item{object, data, knots}{See \code{\link[mgcv]{smooth.construct}}.}
  \item{\dots}{Currently not used.}
}

\value{
  See \code{\link[mgcv]{smooth.construct}}
}


\seealso{
\code{\link{bamlss}}, \code{\link{predict.bamlss}}, \code{\link{bfit}}, \code{\link{boost}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart()

## Estimate model.
f <- num ~ x1 + x2 + x3 + s2(x1,bs="sr") + s2(x2,bs="sr") + s2(x3,bs="sr")

b <- bamlss(f, data = d, optimizer = boost, sampler = FALSE)

plot(b)
}
}

\keyword{regression}

