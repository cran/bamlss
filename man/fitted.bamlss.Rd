\name{fitted.bamlss}
\alias{fitted.bamlss}

\title{
  BAMLSS Fitted Values
}

\description{
  Function to compute fitted values for \code{\link{bamlss}} models. The function calls
  \code{\link{predict.bamlss}} to compute fitted values from samples.
}

\usage{
\method{fitted}{bamlss}(object, model = NULL, term = NULL,
  type = c("link", "parameter"), samples = TRUE,
  FUN = c95, nsamps = NULL, ...)
}

\arguments{
  \item{object}{An object of class \code{"bamlss"}}
  \item{model}{Character or integer, specifies the model for which fitted values should be computed.}
  \item{term}{Character or integer, specifies the model terms for which fitted values are required.
    Note that if \code{samples = TRUE}, e.g., \code{term = c("s(x1)", "x2")} will compute the
    combined fitted values \code{s(x1) + x2}.}
  \item{type}{If \code{type = "link"} the predictor of the corresponding \code{model}
    is returned. If \code{type = "parameter"} fitted values on the distributional parameter scale
    are returned.}
  \item{samples}{Should fitted values be computed using samples of parameters or estimated parameters
    as returned from optimizer functions (e.g., function \code{\link{bfit}} returns
    \code{"fitted.values"}). The former results in a call to \code{\link{predict.bamlss}}, the
    latter simply extracts the \code{"fitted.values"} of the \code{\link{bamlss}} object and
    is not model term specific.}
  \item{FUN}{A function that should be applied on the samples of predictors or
    parameters, depending on argument \code{type}.}
  \item{nsamps}{If the fitted \code{\link{bamlss}} object contains samples of parameters,
    computing fitted values may take quite some time. Therefore, to get a first feeling it can
    be useful to compute fitted values only based on \code{nsamps} samples, i.e., \code{nsamps}
    specifies the number of samples which are extracted on equidistant intervals.}
  \item{\dots}{Arguments passed to function \code{\link{predict.bamlss}}.}
}

\value{
  Depending on arguments \code{model}, \code{FUN} and the structure of the \code{\link{bamlss}}
  model, a list of fitted values or simple vectors or matrices of fitted values.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{predict.bamlss}}.
}

\examples{
\dontrun{## Generate some data.
d <- GAMart()

## Model formula.
f <- list(
  num ~ s(x1) + s(x2) + s(x3) + te(lon,lat),
  sigma ~ s(x1) + s(x2) + s(x3) + te(lon,lat)
)

## Estimate model.
b <- bamlss(f, data = d)

## Fitted values returned from optimizer.
f1 <- fitted(b, model = "mu", samples = FALSE)

## Fitted values returned from sampler.
f2 <- fitted(b, model = "mu", samples = TRUE, FUN = mean)

plot(f1, f2)
}
}

\keyword{regression}
\keyword{models}

