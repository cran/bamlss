\name{summary.bamlss}
\alias{summary.bamlss}
\alias{print.summary.bamlss}

\title{
Summary for BAMLSS
}

\description{
  The function takes an object of class \code{"bamlss"} and produces summaries
  of optimizer and sampler function outputs.
}

\usage{
\method{summary}{bamlss}(object, model = NULL,
  FUN = NULL, parameters = TRUE, ...)

\method{print}{summary.bamlss}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object}{An object of class "bamlss".}
  \item{x}{An oject of class \code{"summary.bamlss"}.}
  \item{model}{Character or integer, specifies the model for which a summary should be computed.}
  \item{FUN}{Function that should be applied on samples, see also function \code{\link{coef.bamlss}}.}
  \item{parameters}{If an optimizer function is applied within the \code{\link{bamlss}} call, should
    the values of the estimated parameters be part of the summary?}
  \item{digits}{Controls number of digits printed in output.}
  \item{\dots}{Other arguments.}
}

\details{
  If the fitted model contains samples, summaries according to the supplied function can be
  computed, e.g., different quantiles of samples. See also function \code{\link{coef.bamlss}}
  that extracts the coefficient summaries.

  If an optimizer function was used within the \code{\link{bamlss}} call, estimated parameters
  will be included per default into the summary.

  Note that summaries not based on samples can be user defined, e.g., as returned from function
  \code{\link{samplestats}} or the return values of optimizer function, e.g., see function
  \code{\link{bfit}}.
}

\value{
\code{summary.bamlss} produces the following summary:
  \item{call}{The initial \code{\link{bamlss}} call.}
  \item{family}{The family that is used for modeling.}
  \item{formula}{The model formula.}
  \item{model.matrix}{Summary of parameteric terms.}
  \item{smooth.construct}{Summary of smooth terms.}
  \item{model.stats}{Other model statistics, e.g., as returned from
    optimizer functions and/or produces by function \code{\link{samplestats}}.}
}

\seealso{
  \code{\link{bamlss}}
}

\examples{
\dontrun{## Generate some data.
d <- GAMart()

## Model formula.
f <- list(
  num ~ s(x1) + s(x2),
  sigma ~ s(x3) + te(lon,lat)
)

## Estimate model.
b <- bamlss(f, data = d)

## Print the summary.
print(summary(b))
}
}

\keyword{models}
\keyword{regression}

