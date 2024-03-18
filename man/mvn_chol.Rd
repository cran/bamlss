\name{mvn_chol}
\alias{mvn_chol}
\title{Cholesky MVN}
\usage{
mvn_chol(k = 2L, ...)
}
\arguments{
\item{k}{integer. The dimension of the multivariate distribution.}

\item{...}{not used.}
}
\value{
a bamlss family.
}
\description{
BAMLSS Family for MVN with Cholesky Parameterization
}
\details{
This is a prototype implementation of a BAMLSS family that models
a multivariate Normal (Gaussian) distribution by a Cholesky
decomposition of the covariance matrix.
}

\seealso{
  \code{\link{mvnchol_bamlss}}
}
