\name{CRPS}
\alias{CRPS}

\title{
  Continuous Rank Probability Score
}

\description{
  The function computes the continuous rank probability score (CRPS). Note that the function
  uses numerical integration, for highly efficient computation please see
  the \pkg{scoringRules} package.
}
 
\usage{
CRPS(object, newdata = NULL,
  interval = c(-Inf, Inf), FUN = mean,
  term = NULL, ...)
}

\arguments{
  \item{object}{An object returned from \code{\link{bamlss}}.}
  \item{newdata}{Optional new data that should be used for calculation.}
  \item{interval}{The interval that should be used for numerical integration}
  \item{FUN}{Function to be applied on the CRPS scores.}
  \item{term}{If required, specify the model terms that should be used within the
    \code{\link{predict.bamlss}} function.}
  \item{\dots}{Arguments passed to function \code{FUN}.}
}

\references{
  Gneiting T, Raftery AE (2007). Strictly Proper Scoring Rules, Prediction, and
  Estimation." Journal of the American Statistical Association, 102(477), 359--378.
  \doi{10.1198/016214506000001437}cd ..

  Gneiting T, Balabdaoui F, Raftery AE (2007). Probabilistic Forecasts, Calibration and
  Sharpness. Journal of the Royal Statistical Society B, 69(2), 243--268.
  \doi{10.1111/j.1467-9868.2007.00587.x}
}

\examples{\dontrun{## Simulate data.
d <- GAMart()

## Model only including covariate x1.
b1 <- bamlss(num ~ s(x1), data = d)

## Now, also including x2 and x2.
b2 <- bamlss(num ~ s(x1) + s(x2) + s(x3), data = d)

## Compare using the CRPS score.
CRPS(b1)
CRPS(b2)
}}

\keyword{regression}
\keyword{models}
\keyword{distribution}
