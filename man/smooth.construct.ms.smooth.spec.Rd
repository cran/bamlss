\name{smooth.construct.ms.smooth.spec}
\alias{smooth.construct.ms.smooth.spec}

\title{Smooth constructor for monotonic P-splines}

\description{The function sets up a smooth term for shape constraint estimation of
  P-spline model terms. Note that this currently only works using boosting and backfitting.}

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
set.seed(123)

n <- 300
x <- runif(n, -2, 3)
y <- sin(x) + rnorm(n, sd = 0.1)

d <- data.frame("y" = y, "x" = x)

## Increasing: constr = 1.
## Decreasing: constr = 2.
b <- bamlss(y ~ s2(x,bs="ms",xt=list(constr=1)),
  data = d, optimizer = opt_bfit, sampler = sam_MVNORM)

## Predict and plot.
p <- predict(b, model = "mu", FUN = c95)
plot(y ~ x)
plot2d(p ~ x, add = TRUE, col.lines = 4, lwd = 2)
}}

\keyword{models}
\keyword{smooth}
\keyword{regression}

