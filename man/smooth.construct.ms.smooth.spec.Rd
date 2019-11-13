\name{smooth.construct.ms.smooth.spec}
\alias{smooth.construct.ms.smooth.spec}

\title{Smooth constructor for monotonic P-splines}

\description{The function sets up a smooth term for shape constraint estimation of
  P-spline model terms. Note that this currently only works using boosting.}

\usage{
\method{smooth.construct}{ms.smooth.spec}(object, data, knots, ...)
}

\arguments{
  \item{object}{Either a smooth specification object, or object of class \code{"bamlss"},
    \code{"bamlss.frame"}, \code{"bamlss.formula"} or \code{"bamlss.terms"}. For smooth
    specification objects, see function \code{\link[mgcv]{smooth.construct}}.}
  \item{data}{A data frame or list, see also see function \code{\link[mgcv]{smooth.construct}}.}
  \item{knots}{See function \code{\link[mgcv]{smooth.construct}}.}
  \item{\dots}{Arguments passed to the smooth term constructor functions.}
}

\value{
See function see \code{\link[mgcv]{smooth.construct}}.
}

\seealso{
  \code{\link{bamlss.frame}}, \code{\link{bamlss.formula}},
  \code{\link{bamlss}}, \code{\link[mgcv]{smooth.construct}}.
}

\examples{
\dontrun{## Generate some data.
d <- GAMart()

## Increasing: constr = 1.
## Decreasing: constr = 2.
b <- bamlss(num ~ s2(x1,bs="ms",xt=list(constr=2)),
  data = d, sampler = FALSE, optimizer = boost)

plot(b)
}}

\keyword{models}
\keyword{smooth}
\keyword{regression}

