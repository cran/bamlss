\name{bboost}
\alias{bboost}
\alias{bboost_plot}
\alias{predict.bboost}

\title{Bootstrap Boosting}

\description{
  Wrapper function for applying bootstrap estimation using gradient boosting.
}

\usage{
## Bootstrap boosting.
bboost(..., data, type = 1, cores = 1,
  n = 2, prob = 0.623, fmstop = NULL,
  trace = TRUE, drop = FALSE, replace = FALSE)

## Plotting function.
bboost_plot(object, col = NULL)

## Predict method.
\method{predict}{bboost}(object, newdata, ..., cores = 1, pfun = NULL)
}

\arguments{
  \item{\dots}{Arguments passed to \code{\link{bamlss}} and \code{\link{predict.bamlss}}.}
  \item{data}{The data frame to be used for modeling.}
  \item{type}{Type of algorithm, \code{type = 1} uses all observations and samples with
    replacement, \code{type = 2} uses only a fraction specified in \code{prob} and samples
    with replacement.}
  \item{cores}{The number of cores to be used.}
  \item{n}{The number of bootstrap iterations.}
  \item{prob}{The fraction that should be used to fit the model in each bootstrap iteration.}
  \item{fmstop}{The function that should return the optimum stopping iteration. The function must
    have two arguments: (1) the \code{model} end (2) the \code{data}. The function must return
    a list with two named arguments: (1) \code{"mstop"} the optimum stopping iteration and (2)
    a vector of the objective criterion that should be evaluated by the hold out sample data during
    each bootstrap iteration. See the examples.}
  \item{trace}{Prints out the current state of the bootstrap algorithm.}
  \item{drop}{Should only the best set of parameters be saved?}
  \item{replace}{Sampling with replacement, or sampling \code{ceiling(nobs * prob)}
    rows of the data for fitting the \code{n} models.}
  \item{object}{The \code{"bboost"} object used for prediction and plotting.}
  \item{col}{The color that should be used for plotting.}
  \item{newdata}{The data frame predictions should be made for.}
  \item{pfun}{The prediction function that should be used, for example \code{\link{predictn}} could
    be used, too. Note that this is experimental.}
}

\value{
  A list of \code{bamlss} objects.
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

## Function for evaluation of hold out sample
## criterion to find the optimum mstop.
fmstop <- function(model, data) {
  p <- predict(model, newdata = data, model = "mu")
  mse <- NULL
  for(i in 1:nrow(model$parameters))
    mse <- c(mse, mean((data$num - p[, i])^2))
  list("MSE" = mse, "mstop" = which.min(mse))
}

## Bootstrap boosted models.
b <- bboost(f, data = d, n = 50, cores = 3, fmstop = fmstop)

## Plot hold out sample MSE.
bboost_plot(b)

## Predict for each bootstrap sample.
nd <- data.frame("x2" = seq(0, 1, length = 100))
p <- predict(b, newdata = nd, model = "mu", term = "x2")
plot2d(p ~ x2, data = nd)
}
}

\keyword{regression}

