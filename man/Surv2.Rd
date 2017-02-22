\name{Surv2}
\alias{Surv2}

\title{
  Create a Survival Object for Joint Models
}
\description{
  This function is only a slight extension of \code{\link[survival]{Surv}} for
  joint models.
}
\usage{
Surv2(..., obs = NULL)
}

\arguments{
  \item{\dots}{Arguments passed to function \code{\link[survival]{Surv}}.}
  \item{obs}{The observed longitudinal response.}
}

\value{
  An object of class \code{"Surv2"} and \code{"matrix"}.
}

\seealso{
  \code{\link{jm.mode}}, \code{\link{jm.mcmc}}, \code{\link{bamlss}}
}

\examples{
## Surv2(time, event, obs = y)
## See the examples of jm.mode() and jm.mcmc()!
}

\keyword{regression}

