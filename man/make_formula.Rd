\name{make_formula}
\alias{make_formula}
\title{Formula Generator}
\usage{
make_formula(formula, type = "basic")
}
\arguments{
\item{formula}{formula.}

\item{type}{character. Type of Cholesky decomposition.}
}
\description{
Generate a formula for a MVN Cholesky model.
}
\details{
This is a helper to generate a formula for a bamlss model with \code{k}-dimensional
multi-variate normal distribution and Cholesky decomposed variance-covariance matrix.
It is helpful if one formula should be used for means, another for all diagonal
entries of the Cholesky factor, and a third one for all lower triangular entries
of the Cholesky factor.
The left hand side has \code{k} elements separated by \code{|}.
The right hand side has one to three elements separated by \code{|} specifying
the formulas used for all means, diagonal entries of the Cholesky factor and
lower triangular entries of the Cholesky factor, respectively.
}
\examples{
f <- O | C | E | A | N ~ s(x1) + s(x2) | s(y) | z 
f2 <- make_formula(f)
f2
}

\seealso{
  \code{\link{mvnchol_bamlss}}
}
