\name{smooth.construct.kr.smooth.spec}
\alias{smooth.construct.kr.smooth.spec}
\alias{Predict.matrix.kriging.smooth}

\title{Kriging Smooth Constructor}

\description{
  This smooth constructor implements a kriging based model term.
}

\usage{
\method{smooth.construct}{kr.smooth.spec}(object, data, knots, ...)
\method{Predict.matrix}{kriging.smooth}(object, data)
}

\arguments{
  \item{object, data, knots}{See \code{\link[mgcv]{smooth.construct}}.}
  \item{\dots}{Currently not  used.}
}

\details{
  This smooth constructor implements univariate and bivariate Kriging terms. The
  basis functions are based on the \code{\link[geoR]{matern}} function. For finding
  knots, a space filling algorithm is used, see \code{\link[fields]{cover.design}}.
}

\value{
  A smooth specification object, see also \code{\link[mgcv]{smooth.construct}.}
}

\references{
  Fahrmeir, L., Kneib, T. , Lang, S., Marx, B. (2013): Regression. Models, Methods and Applications,
  Springer Verlag. \url{http://www.regressionbook.org/}
}

\seealso{
\code{\link{bamlss}}, \code{\link[mgcv]{smooth.construct}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart()

## Estimate model.
f <- num ~ s(x1,bs="kr") + s(x2,bs="kr") + s(x3,bs="kr") + s(lon,lat,bs="kr",k=30)

## Set the seed, estimate model.
set.seed(111)
b <- bamlss(f, data = d)

## Plot estimated effects.
plot(b)
}
}

\keyword{regression}

