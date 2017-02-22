\name{scale2}
\alias{scale2}

\title{
  Scaling Vectors and Matrices
}

\description{
  The function scales numeric objects to specific ranges.
}

\usage{
scale2(x, lower = -1.5, upper = 1.5)
}

\arguments{
  \item{x}{Numeric, vector or matrix.}
  \item{lower}{The upper range.}
  \item{upper}{The lower range.}
}

\value{
  A scaled numeric vector or matrix, scaled to the range provided in \code{lower}
  and \code{upper}.
}

\examples{
set.seed(123)
x <- runif(5)
scale2(x, -1, 1)
scale2(x, 0, 10)
}

\keyword{manip}

