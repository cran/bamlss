\name{gF}
\alias{gF}

\title{
  Get a BAMLSS Family
}

\description{
  Function to get a \code{\link{family.bamlss}} object to be used for fitting. The main
  purpose of this function is to ease the handling of extra arguments to the family object.
}

\usage{
gF(x, ...)
}

\arguments{
  \item{x}{The name of the \code{\link{family.bamlss}} without the \code{".bamlss"} extension.}
  \item{\dots}{Arguments passed to the family object.}
}

\value{
  A \code{\link{family.bamlss}} object.
}

\seealso{
  \code{\link{family.bamlss}}.
}

\examples{
f <- gF(gaussian)
print(f)
}

\keyword{misc}

