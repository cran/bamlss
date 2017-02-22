\name{smooth.construct}
\alias{smooth.construct}
\alias{smooth.construct.bamlss.frame}
\alias{smooth.construct.bamlss.formula}
\alias{smooth.construct.bamlss.terms}

\title{Constructor Functions for Smooth Terms in BAMLSS}

\description{The generic function is only a copy of \code{\link[mgcv]{smooth.construct}} adding
  a \code{\dots} argument. For objects of class \code{"bamlss.frame"} and \code{"bamlss"} the
  method extracts all smooth model terms, see function \code{\link{bamlss.frame}} for details
  on the setup of BAMLSS.}

\usage{
## Function as in package mgcv
## but with additional dots argument.
smooth.construct(object, data, knots, ...)

## For 'bamlss.frame's.
\method{smooth.construct}{bamlss.frame}(object, data = NULL, knots = NULL,
  model = NULL, drop = TRUE, ...)

\method{smooth.construct}{bamlss.formula}(object, data = NULL, knots = NULL,
  model = NULL, drop = TRUE, ...)

\method{smooth.construct}{bamlss.terms}(object, data = NULL, knots = NULL,
  model = NULL, drop = TRUE, ...)
}

\arguments{
  \item{object}{Either a smooth specification object, or object of class \code{"bamlss"},
    \code{"bamlss.frame"}, \code{"bamlss.formula"} or \code{"bamlss.terms"}. For smooth
    specification objects, see function \code{\link[mgcv]{smooth.construct}}.}
  \item{data}{A data frame or list, see also see function \code{\link[mgcv]{smooth.construct}}.}
  \item{knots}{See function \code{\link[mgcv]{smooth.construct}}.}
  \item{model}{Character, specifies for which model parameter the smooth constructs should be
    created.}
  \item{drop}{If there is only one model parameter the returned named list is simplified.}
  \item{\dots}{Arguments passed to the smooth term constructor functions.}
}

\value{
For smooth specification objects see function see \code{\link[mgcv]{smooth.construct}}. For
objects of class \code{"bamlss.frame"} or \code{"bamlss"} the list of smooth constructs,
see function \code{\link{bamlss.frame}} for more details.
}

\seealso{
  \code{\link{bamlss.frame}}, \code{\link{bamlss.formula}},
  \code{\link{bamlss}}, \code{\link[mgcv]{smooth.construct}}.
}

\examples{
## Generate some data.
d <- GAMart()

## Create a "bamlss.frame".
bf <- bamlss.frame(num ~ s(x1) + s(x2), data = d)

## Extract the smooth construct.
sc <- smooth.construct(bf)
str(sc)

## Also possible with formulas.
f <- bamlss.formula(list(
  num ~ s(x1) + te(lon,lat),
  sigma ~ s(x2)
), family = "gaussian")

sc <- smooth.construct(f, data = d)
str(sc)
}

\keyword{models}
\keyword{smooth}
\keyword{regression}

