\name{mvn_modchol}
\alias{mvn_modchol}
\title{Modified Cholesky MVN}
\usage{
mvn_modchol(k = 2L, ...)
}
\arguments{
\item{k}{integer. The dimension of the multivariate distribution.}

\item{...}{not used.}
}
\value{
a bamlss family.
}
\description{
BAMLSS Family for MVN with Modified Cholesky Parameterization
}
\details{
This is a prototype implementation of a BAMLSS family that models
a multivariate Normal (Gaussian) distribution by a modified Cholesky
decomposition of the covariance matrix.
}

\seealso{
  \code{\link{mvnchol_bamlss}}
}
