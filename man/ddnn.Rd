\name{ddnn}
\alias{ddnn}
\alias{predict.ddnn}
\alias{cv_ddnn}

\title{Deep Distributional Neural Network}

\description{
  This function interfaces \pkg{keras} infrastructures for high-level neural networks. The function
  can be used as a standalone model fitting engine such as \code{\link{bamlss}} or as an on top
  model engine to capture special features in the data that could not be captures by other
  model fitting engines.
}

\usage{
## Deep distributional neural net.
ddnn(object, optimizer = "adam",
  learning_rate = 0.01,
  epochs = 100, batch_size = NULL,
  nlayers = 2, units = 100, activation = "relu",
  l1 = NULL, l2 = NULL,
  validation_split = 0.2, early_stopping = TRUE, patience = 50,
  verbose = TRUE, ...)

## Predict method.
\method{predict}{ddnn}(object, newdata,
  model = NULL, type = c("link", "parameter"),
  drop = TRUE, ...)

## CV method for optimizing
## the number of epochs using
## the CRPS.
cv_ddnn(formula, data, folds = 10,
  min_epochs = 300, max_epochs = 400,
  interval = c(-Inf, Inf), ...)
}

\arguments{
  \item{object}{An object of class \code{"bamlss"} or a \code{\link{bamlss.formula}}.}
  \item{optimizer}{Character or call to optimizer functions to be used within \code{\link[keras]{fit}}.
    For character, options are: \code{"adam"} \code{"sgd"}, \code{"rmsprop"}, \code{"adagrad"},
    \code{"adadelta"}, \code{"adamax"}, \code{"adam"}. The default is
    \code{\link[keras]{optimizer_rmsprop}} with learning rate set to \code{1e-04}.}
  \item{learning_rate}{The learning rate of the optimizer.}
  \item{epochs}{Number of times to iterate over the training data arrays, see
    \code{\link[keras]{fit}}.}
  \item{batch_size}{Number of samples per gradient update, see \code{\link[keras]{fit}}.}
  \item{nlayers}{Number of hidden layers.}
  \item{units}{Number of nodes per hidden layer, can be a vector.}
  \item{activation}{Activation functions used for the hidden layers, can be a vector.}
  \item{l1}{Shrinkage parameter for L1 penalty.}
  \item{l2}{Shrinkage parameter for L2 penalty.}
  \item{validation_split}{Proportion of data that should be used for validation.}
  \item{early_stopping}{Logical, should early stopping of the optimizer be applied?}
  \item{patience}{Integer, number of iterations the optimizer waits until early stopping is applied
    after changes get small in validation data set.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{newdata}{A \code{\link[base]{list}} or \code{\link[base]{data.frame}} that should
    be used for prediction.}
  \item{model}{Character or integer specifying for which distributional parameter predictions should
    be computed.}
  \item{type}{If \code{type = "link"} the predictor of the corresponding \code{model}
    is returned. If \code{type = "parameter"} predictions on the distributional parameter scale
    are returned.}
  \item{drop}{If predictions for only one \code{model} are returned, the list structure is dropped.}
  \item{formula}{The model formula.}
  \item{data}{The data used for estimation.}
  \item{folds}{The number of folds that should be generated.}
  \item{min_epochs, max_epochs}{Defines the minimum and maximum epochs thet should be used.}
  \item{interval}{Response interval, see function \code{\link{CRPS}}.}
  \item{\dots}{Arguments passed to \code{\link{bamlss.frame}}.}
}

\details{
  The default \pkg{keras} model is a sequential model with two hidden layers with \code{"relu"}
  activation function and 100 units in each layer. Between each layer is a dropout layer with
  0.1 dropout rate.
}

\value{
  For function \code{ddnn()} an object of class \code{"ddnn"}. Note that extractor
  functions \code{\link[stats]{fitted}} and \code{\link{residuals.bamlss}} can be applied.
  For function \code{predict.ddnn()} a list or vector of predicted values.
}

\section{WARNINGS}{
  The deep learning infrastructure is experimental!
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
b <- ddnn(f, epochs = 2000)

## Plot estimated functions.
par(mfrow = c(1, 2))
plot(x, y)
plot2d(fitted(b)$mu ~ x, add = TRUE)
plot2d(fitted(b)$sigma ~ x,
  ylim = range(c(fitted(b)$sigma, fsigma)))
plot2d(fsigma ~ x, add = TRUE, col.lines = "red")

## Predict with newdata.
nd <- data.frame(x = seq(-6, 6, length = 100))
nd$p <- predict(b, newdata = nd, type = "link")

par(mfrow = c(1, 2))
plot(x, y, xlim = c(-6, 6), ylim = range(c(nd$p$mu, y)))
plot2d(p$mu ~ x, data = nd, add = TRUE)
plot2d(p$sigma ~ x, data = nd,
  ylim = range(c(nd$p$sigma, fsigma)))
plot2d(fsigma ~ x, add = TRUE, col.lines = "red")

## Plot quantile residuals.
e <- residuals(b)
plot(e)
}
}

\keyword{regression}

