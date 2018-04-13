\name{n}
\alias{n}
\alias{smooth.construct.nn.smooth.spec}
\alias{n.weights}

\title{Neural Networks for BAMLSS}

\description{
  This smooth constructor implements single hidden layer neural networks with sigmoidal
  activation functions.
}

\usage{
## The neural network smooth constructor.
n(..., k = 10)

## Initialize weights.
n.weights(nodes, k, r = NULL, s = NULL, type = c("sigmoid", "gauss"), ...)

## For mgcv.
\method{smooth.construct}{nn.smooth.spec}(object, data, knots, ...)
}

\arguments{
  \item{\dots}{For function \code{n()} a formula of the type \code{~x1+x2+x3} that specifies
    the covariates that should be modeled by the neural network. For function
    \code{boost.nnet.predict()} arguments passed to function \code{\link{predict.bamlss}}.}
  \item{k}{For function \code{n()}, the number of hidden nodes of the network. Note that one can set
    an argument \code{split = TRUE} to split up the neural network into, e.g., \code{nsplit = 5}
    parts with \code{k} nodes each. For function \code{nnet.weights()}, argument \code{k}
    is the number of input variables of the network (number of covariates).}
  \item{nodes}{Number of nodes for each layer, i.e., can also be a vector.}
  \item{r, s}{Parameters controlling the shape of the activation functions.}
  \item{type}{The type of activation function that should be used.}
  \item{object, data, knots}{See \code{\link[mgcv]{smooth.construct}}.}
}

\value{
  Function \code{n()}, similar to function \code{\link[mgcv]{s}} a simple smooth specification
  object.
}


\seealso{
\code{\link{bamlss}}, \code{\link{predict.bamlss}}, \code{\link{bfit}}, \code{\link{boost}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart()

## Estimate model.
f <- num ~ n(x1,k=50) + n(x2,k=50) + n(x3,k=50) + n(~lon+lat,k=50,split=TRUE)

## Set the seed, weights are initialized randomly.
set.seed(111)
b <- bamlss(f, data = d, sampler = FALSE)

## Plot estimated effects.
plot(b)

## Boosted version.
set.seed(111)
b <- bamlss(f, data = d, sampler = FALSE, optimizer = boost,
  nu = 0.05, maxit = 4000)

plot(b)

## Second example with MCMC.
set.seed(222)
n <- 2000
d <- data.frame("x1" = runif(n, 0, pi), "x2" = runif(n, 0, pi))
d$y <- sin(d$x1) * sin(d$x2) + rnorm(n, sd = exp(-2 + sin(2 * d$x2)))

f <- list(
  y ~ n(~x1+x2,k=50),
  sigma ~ n(x2)
)

b <- bamlss(f, data = d)
plot(b, ask = FALSE)

## Initialize some weights.
set.seed(111)
n.weights(10, 3)
}
}

\keyword{regression}

