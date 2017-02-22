\name{predict.bamlss}
\alias{predict.bamlss}

\title{
  BAMLSS Prediction 
}

\description{
  Takes a fitted \code{\link{bamlss}} object and computes predictions. Predictions can be based on
  estimated parameters of optimizer functions or on samples returned from sampler functions.
  It is possible to compute predictions on multiple cores using the
  \code{\link[parallel]{parallel}} and to chunk predictions to save computation time and memory
  storage. Predictions can be computed for full distributional parameters or specific model
  terms. If a \code{link{bamlss}} model was fitted on multiple cores, i.e., the samples
  are provided as \code{link{mcmc.list}} where each list entry represents samples from one
  core, function \code{predict.bamlss()} computes combined predictions based on samples of all
  cores.
}

\usage{
\method{predict}{bamlss}(object, newdata, model = NULL, term = NULL,
  match.names = TRUE, intercept = TRUE, type = c("link", "parameter"),
  FUN = function(x) { mean(x, na.rm = TRUE) }, trans = NULL,
  what = c("samples", "parameters"), nsamps = NULL,
  verbose = FALSE, drop = TRUE,
  cores = NULL, chunks = 1, ...)
}

\arguments{
  \item{object}{An object of class \code{"bamlss"}}
  \item{newdata}{A data frame or list containing the values of the model
    covariates at which predictions are required. Note that depending
    on argument \code{term}, only covariates that are needed by the
    corresponding model terms need to be supplied.}
  \item{model}{Character or integer, specifies the model for which predictions should be computed.}
  \item{term}{Character or integer, specifies the model terms for which predictions are required.
    Note that, e.g., \code{term = c("s(x1)", "x2")} will compute the combined prediction
    \code{s(x1) + x2}.}
  \item{match.names}{Should partial string matching be used to select the \code{term}s for
    prediction. Note that, e.g., \code{term = "x1"} will select all terms including \code{"x1"} if
    \code{match.names = TRUE}.}
  \item{intercept}{Should the intercept be included?}
  \item{type}{If \code{type = "link"} the predictor of the corresponding \code{model}
    is returned. If \code{type = "parameter"} predictions on the distributional parameter scale
    are returned.}
  \item{FUN}{A function that should be applied on the samples of predictors or
    parameters, depending on argument \code{type}.}
  \item{trans}{A transformer function or named list of transformer functions that computes
    transformed predictions. If \code{trans} is a list, the list names must match the names
    of the parameters of the \code{\link{bamlss.family}}.}
  \item{what}{Predictions can be computed from samples or estimated parameters of optimizer
    functions. If no samples are available the default is to use estimated parameters.}
  \item{nsamps}{If the fitted \code{\link{bamlss}} object contains samples of parameters,
    computing predictions may take quite some time. Therefore, to get a first feeling it can
    be useful to compute predictions only based on \code{nsamps} samples, i.e., \code{nsamps}
    specifies the number of samples which are extracted on equidistant intervals.}
  \item{verbose}{If predictions are chunked, information on the prediction process can be printed.}
  \item{drop}{If predictions for only one \code{model} are returned, the list structure is dropped.}
  \item{cores}{Specifies the number of cores that should be used for prediction. Note that
    this functionality is based on the \code{\link[parallel]{parallel}} package.}
  \item{chunks}{Should computations be split into \code{chunks}? Prediction is then processed
    sequentially.}
  \item{\dots}{Arguments passed to prediction functions that are part of a \code{\link{bamlss.family}}
    object, i.e., the objects has a \code{$predict()} function that should be used instead.}
}

\value{
  Depending on arguments \code{model}, \code{FUN} and the structure of the \code{\link{bamlss}}
  model, a list of predictions or simple vectors or matrices of predictions.
}

\seealso{
  \code{link{bamlss}}, \code{\link{fitted.bamlss}}.
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

## Predictions.
p <- predict(b)
str(b)

## Prediction for "mu" model and term "s(x2)".
p <- predict(b, model = "mu", term = "s(x2)")

## Plot effect
plot2d(p ~ x2, data = d)

## Same for "sigma" model.
p <- predict(b, model = "sigma", term = "s(x2)")
plot2d(p ~ x2, data = d)

## Prediction for "mu" model and term "s(x1)" + "s(x2)"
## without intercept.
p <- predict(b, model = "mu", term = c("s(x1)", "s(x2)"),
  intercept = FALSE)

## Prediction based on quantiles.
p <- predict(b, model = "mu", term = "s(x2)", FUN = c95)
plot2d(p ~ x2, data = d)

## Extract samples of predictor for "s(x2)".
p <- predict(b, model = "mu", term = "s(x2)",
  intercept = FALSE, FUN = function(x) { x })
print(dim(p))
plot2d(p ~ x2, data = d, col.lines = rgb(0.1, 0.1, 0.1, alpha = 0.1))

## Or using specific combinations of terms.
p <- predict(b, model = "mu", term = c("s(x2)", "te(lon,lat)"),
  intercept = FALSE, FUN = function(x) { x })
head(p)

## Prediction using new data.
## Only need x3 data when predicting
## for s(x3).
nd <- data.frame("x3" = seq(0, 1, length = 100))
nd <- cbind(nd, predict(b, newdata = nd, term = "s(x3)"))
print(head(nd))
plot2d(mu ~ x3, data = nd)
plot2d(sigma ~ x3, data = nd)
}
}

\keyword{regression}
\keyword{models}

