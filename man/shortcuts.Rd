\name{boost2}
\alias{boost2}
\alias{lasso2}
\alias{bayesx2}

\title{Some Shortcuts}

\description{
  Some simple shortcuts to model fitting engines.
}

\usage{
## BayesX.
bayesx2(...)

## Gradient boosting.
boost2(...)

## Lasso.
lasso2(...)
}

\arguments{
  \item{\dots}{Arguments passed to \code{\link{bamlss}} and \code{\link{predict.bamlss}}.}
}

\value{
  A \code{bamlss} object.
}


\seealso{
\code{\link{bamlss}}, \code{\link{boost}}, \code{\link{lasso}}, \code{\link{BayesX}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart()

## Estimate model.
f <- num ~ s(x1) + s(x2) + s(x3) + s(lon,lat)

## Boosted model.
b <- boost2(f, data = d)

## Plot estimated effects.
plot(b)
}
}

\keyword{regression}

