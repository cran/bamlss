\name{response_name}
\alias{response_name}

\title{
  Extract the reponse name of a \code{\link{bamlss.frame}} object.
}

\description{
  This is a small helper function to quickly extract the response name(s) of an object of
  class \code{"bamlss.frame"} or \code{"bamlss"}.
}

\usage{
response_name(object, ...)
}

\arguments{
  \item{object}{An object of class \code{"bamlss.frame"} or \code{"bamlss"}.}
  \item{\dots}{Not used.}
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}}
}

\examples{
## Simulate some data.
d <- GAMart()

## Create a bamlss.frame.
bf <- bamlss.frame(num ~ s(x1) + s(x2) + s(x3), data = d)

## Extract the response name.
response_name(bf)
}

\keyword{regression}

