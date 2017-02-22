\name{GAMart}
\alias{GAMart}

\title{GAM Artificial Data Set}

\description{
  This function creates artificial GAM-type \code{\link{data.frame}}s.
  The function is mainly used for testing purposes.
}

\usage{
GAMart(n = 500, sd = 0.1, seed = FALSE,
  ti = c("none", "vcm", "main", "both"))
}

\arguments{
  \item{n}{The number of observations.}
  \item{sd}{Standard deviation of the normal errors.}
  \item{seed}{Sets the \code{seed} to \code{111}.}
  \item{ti}{For tensor product interaction term, the type of interaction.}
}

\examples{
d <- GAMart()
head(d)

\dontrun{b <- bamlss(num ~ s(x1) + s(x2) + s(x3) + te(lon,lat), data = d)
plot(b)
}
}

\keyword{datasets}
