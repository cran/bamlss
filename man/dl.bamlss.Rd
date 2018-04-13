\name{dl.bamlss}
\alias{dl.bamlss}
\alias{predict.dl.bamlss}

\title{Deep Learning BAMLSS}

\description{
  This function interfaces \pkg{keras} infrastructures for high-level neural networks. The function
  can be used as a standalone model fitting engine such as \code{\link{bamlss}} or as an on top
  model engine to capture special features in the data that could not be captures by other
  model fitting engines.
}

\usage{
## Deep learning bamlss.
dl.bamlss(object, offset = NULL, weights = NULL,
  eps = .Machine$double.eps^0.25, maxit = 100,
  force.stop = TRUE, epochs = 30, optimizer = NULL,
  batch_size = NULL, keras.model = NULL,
  verbose = TRUE, digits = 4, ...)

## Predict method.
\method{predict}{dl.bamlss}(object, newdata, model = NULL, drop = TRUE, ...)
}

\arguments{
  \item{object}{An object of class \code{"bamlss"} or a \code{\link{bamlss.formula}}.}
  \item{offset}{A \code{\link[base]{list}} or \code{\link[base]{data.frame}}. Can be used to supply
    model offsets for use in fitting, e.g., fitted values from an initial model fit.}
  \item{weights}{Prior weights on the data.}
  \item{eps}{The relative convergence tolerance of the algorithm.}
  \item{maxit}{Integer, maximum number of iterations of the algorithm.}
  \item{force.stop}{Logical. should the algorithm stop if relative change is smaller than
    \code{eps}.}
  \item{epochs}{Number of times to iterate over the training data arrays, see
    \code{\link[keras]{fit}}.}
  \item{optimizer}{Character or call to optimizer functions to be used within \code{\link[keras]{fit}}.
    For character, options are: \code{"adam"} \code{"sgd"}, \code{"rmsprop"}, \code{"adagrad"},
    \code{"adadelta"}, \code{"adamax"}, \code{"nadam"}. The default is
    \code{\link[keras]{optimizer_rmsprop}} with learning rate set to \code{1e-04}.}
  \item{batch_size}{Number of samples per gradient update, see \code{\link[keras]{fit}}.}
  \item{keras.model}{A compiled model using \pkg{keras}, e.g., using
    \code{\link[keras]{keras_model_sequential}} and \code{\link[keras]{compile}}. Note
    that the last layer only has one unit with linear activation function.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{newdata}{A \code{\link[base]{list}} or \code{\link[base]{data.frame}} that should
    be used for prediction.}
  \item{model}{Character or integer specifying for which distributional parameter predictions should
    be computed.}
  \item{drop}{If predictions for only one \code{model} are returned, the list structure is dropped.}
  \item{\dots}{For function \code{dl.boost()}, arguments passed to \code{\link{bamlss.frame}}.}
}

\details{
  The default \pkg{keras} model is a sequential model with two hidden layers with \code{"relu"}
  activation function and 100 units in each layer. Between each layer is a dropout layer with
  0.1 dropout rate.
}

\value{
  For function \code{dl.bamlss()} an object of class \code{"dl.bamlss"}. Note that extractor
  functions \code{\link[stats]{fitted}} and \code{\link{residuals.bamlss}} can be applied.
  For function \code{predict.dl.bamlss()} a list or vector of predicted values.
}

\section{WARNINGS}{
  The BAMLSS deep learning infrastructure is still experimental!
}

\seealso{
\code{\link{bamlss.frame}}, \code{\link{bamlss}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
n <- 300
x <- runif(n, -3, 3)
fsigma <- -2 + cos(x)
y <- sin(x) + rnorm(n, sd = exp(fsigma))

## Setup model formula.
f <- list(
  y ~ x,
  sigma ~ x
)

## Fit neural network.
library("keras")
b <- dl.bamlss(f)

## Plot estimated functions.
par(mfrow = c(1, 2))
plot(x, y)
plot2d(fitted(b)$mu ~ x, add = TRUE)
plot2d(fitted(b)$sigma ~ x,
  ylim = range(c(fitted(b)$sigma, fsigma)))
plot2d(fsigma ~ x, add = TRUE, col.lines = "red")

## Another example identifying structures that are
## not captured by the initial model.
set.seed(123)
d <- GAMart()
b1 <- bamlss(num ~ s(x1) + s(x2) + s(x3), data = d, sampler = FALSE)
b2 <- dl.bamlss(num ~ lon + lat, data = d, offset = fitted(b1))
p <- predict(b2, model = "mu")
par(mfrow = c(1, 1))
plot3d(p ~ lon + lat, data = d, symmetric = FALSE)
}
}

\keyword{regression}

