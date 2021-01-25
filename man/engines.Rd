\name{engines}
\alias{engines}

\title{
  Show Available Engines for a Family Object
}

\description{
  The function shows available optimizer and sampling engines for a given family object.
}
 
\usage{
engines(family, ...)
}

\arguments{
  \item{family}{A family object or the name of the family.}
  \item{\dots}{Further family objects or names.}
}

\examples{
engines(gaussian_bamlss, "gamma", cox_bamlss)
}

\keyword{regression}
\keyword{models}
\keyword{distribution}
