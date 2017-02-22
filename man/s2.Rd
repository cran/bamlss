\name{s2}
\alias{s2}

\title{
  Special Smooths in BAMLSS Formulae
}

\description{
  This is a simple wrapper function to define special smooth terms
  in BAMLSS formulae. The function calls the smooth term constructor function
  \code{\link[mgcv]{s}}. The return value of \code{\link[mgcv]{s}} is only slightly
  modified, such that function \code{\link{bamlss.frame}} identifies this term as a special
  term and uses the appropriate (internal) infrastructures.

  This structure is useful when the model term structure is, e.g., not a
  linear combination of a design matrix and coefficients. See the example section of
  function \code{\link{bfit}} on how to use this setup.
}

\usage{
s2(...)
}

\arguments{
  \item{\dots}{Arguments passed to function \code{\link[mgcv]{s}}.}
}

\value{
  Slightly modified return value of function \code{\link[mgcv]{s}}.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}}, \code{\link{bamlss.formula}}, \code{\link{bfit}}
}

\examples{
print(names(s(x)))
print(names(s2(x)))
}

\keyword{regression}

