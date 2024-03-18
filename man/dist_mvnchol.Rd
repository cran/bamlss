\name{dist_mvnchol}
\alias{dist_mvnchol}
\title{Cholesky MVN (disttree)}
\usage{
dist_mvnchol(k, r = k - 1L, type = c("basic", "modified", "chol"), ...)
}
\arguments{
\item{k}{integer. The dimension of the multivariate distribution.}

\item{r}{Integer, the number of off-diagonals to model (AD-r covariance).}

\item{type}{character. Choose \code{"basic"} Cholesky decomposition or \code{"modified"}
Cholesky decomposition. (For back compatibility \code{"chol"} is identical to \code{"basic"}.)}

\item{...}{not used.}
}
\value{
a bamlss family.
}
\description{
disttree Families for MVN with Cholesky Parameterization
}
\details{
NOTE: These functions are under development!! 
disttree families that models a multivariate Normal (Gaussian)
distribution by (modified) Cholesky decomposition of the covariance
matrix.
}
