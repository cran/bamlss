\name{Volcano}
\alias{Volcano}

\title{
  Artificial Data Set based on Auckland's Maunga Whau Volcano
}

\description{
  This function creates a data set based on the \code{\link[datasets]{volcano}}
  data by adding normal errors to the topographic information.
}

\usage{
Volcano(sd = 0.3)
}

\arguments{
  \item{sd}{The standard deviation of the normal errors.}
}

\value{
  A data frame with coordinates and noisy elevation.
}

\seealso{
\code{\link[datasets]{volcano}}
}

\examples{
d <- Volcano()
head(d)

\dontrun{b <- bamlss(y ~ te(lon,lat,k=10), data = d)
plot(b, theta = -130)
}
}

\keyword{datasets}

