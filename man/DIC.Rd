\name{DIC}
\alias{DIC}

\title{
  Deviance Information Criterion
}

\description{
  Generic function returning the deviance information criterion (DIC) of a fitted model object.
}
 
\usage{
DIC(object, ...)
}

\arguments{
  \item{object}{A fitted model object for which there exists a \code{DIC} method.}
  \item{\dots}{Optionally more fitted model objects.}
}

\examples{
\dontrun{d <- GAMart()
b1 <- bamlss(num ~ s(x1), data = d)
b2 <- bamlss(num ~ s(x1) + s(x2), data = d)
DIC(b1, b2)
}
}

\keyword{regression}
