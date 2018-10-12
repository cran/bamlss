\name{n}
\alias{n}
\alias{n.weights}
\alias{make_weights}
\alias{predictn}

\title{Neural Networks for BAMLSS}

\description{
  This smooth constructor implements single hidden layer neural networks.
}

\usage{
## The neural network smooth constructor.
n(..., k = 10, type = 2)

## Initialize weights.
n.weights(nodes, k, r = NULL, s = NULL,
  type = c("sigmoid", "gauss", "softplus", "cos", "sin"),
  x = NULL, ...)

## Second weights initializer, internally calls n.weights.
make_weights(object, data, dropout = 0.2)

## Boosted neural net predictions.
predictn(object, newdata, model = NULL,
  mstop = NULL, type = c("link", "parameter"))
}

\arguments{
  \item{\dots}{For function \code{n()} a formula of the type \code{~x1+x2+x3} that specifies
    the covariates that should be modeled by the neural network. For function \code{predictn()},
    arguments to be passed to \code{\link{predict.bamlss}}.}
  \item{k}{For function \code{n()}, the number of hidden nodes of the network. Note that one can set
    an argument \code{split = TRUE} to split up the neural network into, e.g., \code{nsplit = 5}
    parts with \code{k} nodes each. For function \code{n.weights()}, argument \code{k}
    is the number of input variables of the network (number of covariates).}
  \item{type}{Integer. Type \code{1} fits a complete network in each boosting iteration, \code{type = 2} selects
    the best fitting node in each boosting iteration. for function \code{n.weights()}, the type of
    activation function that should be used. For function \code{predictn()}, the type of prediction
    that should be computed.}
  \item{nodes}{Number of nodes for each layer, i.e., can also be a vector.}
  \item{r, s}{Parameters controlling the shape of the activation functions.}
  \item{x}{A scaled covariate matrix, the data will be used to identify the range of the weights.}
  \item{object, data}{See \code{\link[mgcv]{smooth.construct}}. For function \code{predictn()},
    a boosted \code{"bamlss"} object.}
  \item{dropout}{The fraction of inner weights that should be set to zero.}
  \item{newdata}{The data frame that should be used for prediction.}
  \item{model}{For which parameter of the distribution predictions should be computed.}
  \item{mstop}{The stopping iteration for which predictions should be computed. The default
    is to return a matrix of predictions, each column represents the prediction of one
    boosting iteration.}
}

\value{
  Function \code{n()}, similar to function \code{\link[mgcv]{s}} a simple smooth specification
  object.
}


\seealso{
\code{\link{bamlss}}, \code{\link{predict.bamlss}}, \code{\link{bfit}}, \code{\link{boost}}
}

\examples{
## ... coming soon ...!
}

\keyword{regression}

