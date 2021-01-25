\name{gamlss_distributions}
\alias{gamlss_distributions}

\title{
  Extract Distribution families of the \pkg{gamlss.dist} Package
}

\description{
  The functions searches in the \pkg{gamlss.dist} namespace for available
  distributions. It returns a named list of faily functions which can be used
  with \code{\link{bamlss}}.
}

\usage{
gamlss_distributions(type = c("continuous", "discrete"))
}

\arguments{
  \item{type}{Character specifying the type of distribution to be extracted.}
}

\seealso{
  \code{\link{bamlss}}
}

\examples{\dontrun{
dists <- gamlss_distributions(type = "continuous")
print(dists)
}}


\keyword{regression}

