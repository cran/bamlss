\name{WAIC}
\alias{WAIC}

\title{
  Watanabe-Akaike Information Criterion (WAIC)
}

\description{
  Function returning the Watanabe-Akaike Information Criterion (WAIC) of a fitted model object.
}
 
\usage{
WAIC(object, ..., newdata = NULL)
}

\arguments{
  \item{object}{A fitted model object which contains MCMC samples.}
  \item{\dots}{Optionally more fitted model objects.}
  \item{newdata}{Optionally, use new data for computing the WAIC.}
}

\value{
  A data frame containing the WAIC and estimated number of parameters.
}

\references{
  Watanabe S. (2010). Asymptotic Equivalence of {B}ayes Cross Validation and Widely
  Applicable Information Criterion in Singular Learning Theory. \emph{The Journal of Machine
  Learning Research}, \bold{11}, 3571--3594.
  \url{https://jmlr.org/papers/v11/watanabe10a.html}
}

\examples{
\dontrun{d <- GAMart()
b1 <- bamlss(num ~ s(x1), data = d)
b2 <- bamlss(num ~ s(x1) + s(x2), data = d)
WAIC(b1, b2)
}
}

\keyword{regression}
