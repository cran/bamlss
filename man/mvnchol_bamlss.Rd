\name{mvnchol_bamlss}
\alias{mvnchol_bamlss}
\title{Cholesky MVN}
\usage{
mvnchol_bamlss(k, type = c("basic", "modified", "chol"), ...)
}
\arguments{
\item{k}{integer. The dimension of the multivariate distribution.}

\item{type}{character. Choose \code{"basic"} Cholesky decomposition or \code{"modified"}
Cholesky decomposition. (For back compatibility \code{"chol"} is identical to \code{"basic"}.)}

\item{...}{not used.}
}
\value{
a bamlss family.
}
\description{
BAMLSS Families for MVN with Cholesky Parameterization
}
\details{
BAMLSS families that models a multivariate Normal (Gaussian)
distribution by (modified) Cholesky decomposition of the covariance
matrix.

For examples see \code{\link{TempIbk}}.
}

\seealso{
  \code{\link{simdata}}, \code{\link{TempIbk}}
}
