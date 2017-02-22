\name{randomize}
\alias{randomize}

\title{
  Transform Smooth Constructs to Random Effects
}

\description{
  The transformer function takes a \code{\link{bamlss.frame}} object and transforms all
  \code{\link[bamlss]{smooth.construct}}s into a random effects representation. Note that this
  is only possible for smooth terms with a single smoothing variance. The function
  is based on function \code{\link[mgcv]{smooth2random}}.
}

\usage{
randomize(x)
}

\arguments{
  \item{x}{Object returned from function \code{\link{bamlss.frame}}.}
}

\details{
  The decomposition is achieved by a spectral decomposition of the penalty and design matrix by
  finding a basis of the null space of the penalty matrix. This feature is used, e.g., for the
  \code{\link{JAGS}} sampler function. For more details see also \code{\link[mgcv]{jagam}}.
}

\value{
  A transformed \code{\link{bamlss.frame}}. To each \code{\link[bamlss]{smooth.construct}} model term
  an element named \code{"Xf"}, the fixed effects design matrix, and an element \code{"rand$Xr"},
  the random effects design matrix, is added. In addition, for re-transforming parameters elements
  \code{"trans.U"} and \code{"trans.D"} are supplied. See also function
  \code{\link[mgcv]{smooth2random}}.
}

\references{
  Fahrmeir L, Kneib T, Lang S, Marx B (2013). Regression - Models, Methods and Applications.
    Springer-Verlag, Berlin. ISBN 978-3-642-34332-2.

  Wood S.N. (2006). Generalized Additive Models: An Introduction with
    \R. Chapman and Hall/CRC.
}

\seealso{
  \code{\link{bamlss.frame}}, \code{\link{bamlss}}, \code{\link[mgcv]{smooth2random}}.
}

\examples{
## Simulate data.
d <- GAMart()

## Create a "bamlss.frame".
bf <- bamlss.frame(num ~ s(x1) + s(x2) + s(x3) + s(lon,lat), data = d)

## Structure of the "s(x1)" smooth.construct.
str(bf$x$mu$smooth.construct[["s(x1)"]])

## Transform.
bf <- randomize(bf)

## New structure adding fixed
## and random effect matrices.
str(bf$x$mu$smooth.construct[["s(x1)"]])
}

\keyword{regression}
\keyword{smooth}

