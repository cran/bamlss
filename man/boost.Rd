\name{boost}
\alias{boost}
\alias{boostm}
\alias{boost_summary}
\alias{boost_plot}
\alias{print.boost_summary}
\alias{plot.boost_summary}
\alias{boost_frame}

\title{Boosting BAMLSS}

\description{
  Optimizer functions for gradient and likelihood boosting with \code{\link{bamlss}}. In each
  boosting iteration the function selects the model term with the largest contribution to the
  log-likelihood, AIC or BIC.
}

\usage{
## Gradient boosting optimizer.
boost(x, y, family, weights = NULL,
  offset = NULL, nu = 0.1, nu.adapt = TRUE, df = 4, maxit = 400,
  mstop = NULL, maxq = NULL, qsel.splitfactor = FALSE,
  verbose = TRUE, digits = 4, flush = TRUE,
  eps = .Machine$double.eps^0.25,
  nback = NULL, plot = TRUE, initialize = TRUE,
  stop.criterion = NULL, select.type = 1,
  force.stop = TRUE, hatmatrix = !is.null(stop.criterion),
  reverse.edf = FALSE, approx.edf = FALSE,
  always = FALSE, ...)

## Modified likelihood based boosting.
boostm(x, y, family, offset = NULL,
  nu = 0.1, df = 3, maxit = 400, mstop = NULL,
  verbose = TRUE, digits = 4, flush = TRUE,
  eps = .Machine$double.eps^0.25, plot = TRUE,
  initialize = TRUE, stop.criterion = "BIC",
  force.stop = !is.null(stop.criterion),
  do.optim = TRUE, always = FALSE, ...)

## Boosting summary extractor.
boost_summary(object, ...)

## Plot all boosting paths.
boost_plot(x, which = c("loglik", "loglik.contrib", "parameters",
  "aic", "bic", "user"), intercept = TRUE, spar = TRUE, mstop = NULL,
  name = NULL, drop = NULL, labels = NULL, color = NULL, ...)

## Boosting summary printing and plotting.
\method{print}{boost_summary}(x, summary = TRUE, plot = TRUE,
  which = c("loglik", "loglik.contrib"), intercept = TRUE,
  spar = TRUE, ...)
\method{plot}{boost_summary}(x, ...)

## Model frame for out-of-sample selection.
boost_frame(formula, train, test, family = "gaussian", ...)
}

\arguments{
  \item{x}{For function \code{boost()} the \code{x} list, as returned from function
    \code{\link{bamlss.frame}}, holding all model matrices and other information that is used for
    fitting the model. For the plotting function the corresponding \code{\link{bamlss}} object
    fitted with the \code{boost()} optimizer.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.}
  \item{weights}{Prior weights on the data, as returned from function \code{\link{bamlss.frame}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{nu}{Numeric, between [0, 1], controls the step size, i.e., the amount
    that should be added to model term parameters.}
  \item{nu.adapt}{Logical. If set to TRUE (default) step size \code{nu} is divided by 2,
    if current boosting iteration did not improve the loglikelihood.}
  \item{df}{Integer, defines the initial degrees of freedom that should be assigned
    to each smooth model term. May also be a named vector, the names must match the model term
    labels, e.g., as provided in \code{\link{summary.bamlss}}.}
  \item{maxit}{Integer, the maximum number of boosting iterations.}
  \item{mstop}{For convenience, overwrites \code{maxit}.}
  \item{maxq}{Integer, defines the maximum number of selected base-learners. The algorithm stops
    if this numer is exceeded.}
  \item{qsel.splitfactor}{Logical, if set to \code{TRUE} dummy variables of categorical predictors are counted individually.}
  \item{name}{Character, the name of the coefficient (group) that should be plotted. Note that
    the string provided in \code{name} will be removed from the labels on the 4th axis.}
  \item{drop}{Character, the name of the coefficient (group) that should not be plotted.}
  \item{labels}{A character string of labels that should be used on the 4 axis.}
  \item{color}{Colors or color function that creates colors for the (group) paths.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{flush}{use \code{\link{flush.console}} for displaying the current output in the console.}
  \item{eps}{The tolerance used as stopping mechanism, see argument \code{nback}.}
  \item{nback}{Integer. If \code{nback} is not \code{NULL}, then the algorithm stops if the
    the change in the log-likelihood of the last \code{nback} iterations is smaller or
    equal to \code{eps}. If \code{maxit = NULL} the maximum number of iterations is set to 10000.}
  \item{plot}{Should the boosting summary be printed and plotted?}
  \item{initialize}{Logical, should intercepts be initialized?}
  \item{stop.criterion}{Character, selects the information criterion that should be used
    to determine the optimum number of boosting iterations. Either \code{"AIC"} or \code{"BIC"}
    is possible. Note that this feature requires to compute hat-matrices for each distributional
    parameter, therefore, the routine may be slow and computer storage intensive.}
  \item{select.type}{Should model terms be selected by the log-likelihood contribution,
    \code{select.type = 1}, or by the corresponding \code{stop.criterion}, \code{select.type = 2}.}
  \item{force.stop}{Logical, should the algorithm stop if the information criterion increases?}
  \item{do.optim}{Logical. Should smoothing parameters be optimized in each boosting iteration?}
  \item{hatmatrix}{Logical, if set to \code{TRUE} the hat-matrices for each distributional parameter
    will be computed. The hat-matrices are used to determine the effective (equivalent) degrees of
    freedom in each boosting iteration, i.e., it is possible to compute information criteria
    like the AIC or BIC for selecting the optimum number of boosting iterations.}
  \item{reverse.edf}{Logical. Instead of computing degrees of freedom with hat-matrices, the actual
    smoothing parameters are reverse engineered to compute the corresponding actual smoother matrix.
    Note that this option is still experimental.}
  \item{approx.edf}{Logical. Another experimental and fast approximation of the degrees of freedom.}
  \item{always}{Logical or character. Should the intercepts forced to be updated in each boosting iteration?
    If \code{always = TRUE} each intercept of each distributional parameter is updated,
    if \code{always = "best"} only the intercept corresponding to the distributional of the best
    fitting model term is updated.}
  \item{object}{A \code{\link{bamlss}} object that was fitted using \code{boost()}.}
  \item{summary}{Should the summary be printed?}
  \item{which}{Which of the three provided plots should be created?}
  \item{intercept}{Should the coefficient paths of intercepts be dropped in the plot?}
  \item{spar}{Should graphical parmeters be set with \code{\link[graphics]{par}}?}
  \item{formula}{See \code{\link{bamlss.frame}}.}
  \item{train, test}{Data frames used for training and testing the model..}
  \item{\dots}{For function \code{boost()}, arguments passed to \code{\link{bamlss.engine.setup}}.
    for function \code{boost_summary()} arguments passed to function \code{print.boost_summary()}.}
}

\value{
  For function \code{boost_summary()} a list containing information on selection frequencies etc.
  For function \code{boost()} and \code{boostm()} a list containing the following objects:
  \item{fitted.values}{A named list of the fitted values based on the last boosting iteration
    of the modeled parameters of the selected distribution.}
  \item{parameters}{A matrix, each row corresponds to the parameter values of one boosting iteration.}
  \item{boost_summary}{The boosting summary which can be printed and plotted.}
}

\section{WARNINGS}{
  The function does not take care of variable scaling for the linear parts! This must be done by the
  user, e.g., one option is to use argument \code{scale.d} in function \code{\link{bamlss.frame}},
  which uses \code{\link[base]{scale}}.

  Function \code{boost()} does not select the optimum stopping iteration! The modified likelihood
  based algorithm implemented in function \code{boostm()} is still experimental!
}

\seealso{
\code{\link{bamlss.frame}}, \code{\link{bamlss}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart()

## Estimate model.
f <- num ~ x1 + x2 + x3 + lon + lat +
  s(x1) + s(x2) + s(x3) + s(lon) + s(lat) + te(lon,lat)

b <- bamlss(f, data = d, optimizer = boost,
  sampler = FALSE, scale.d = TRUE, nu = 0.01,
  maxit = 1000, plot = FALSE)

## Plot estimated effects.
plot(b)

## Print and plot the boosting summary.
boost_summary(b, plot = FALSE)
boost_plot(b, which = 1)
boost_plot(b, which = 2)
boost_plot(b, which = 3, name = "mu.s.te(lon,lat).")

## Extract estimated parameters for certain
## boosting iterations.
parameters(b, mstop = 1)
parameters(b, mstop = 100)

## Also works with predict().
head(do.call("cbind", predict(b, mstop = 1)))
head(do.call("cbind", predict(b, mstop = 100)))

## Another example using the modified liklihood
## bootsing algorithm.
f <- list(
  num ~ x1 + x2 + x3 + lon + lat +
    s(x1) + s(x2) + s(x3) + s(lon) + s(lat) + te(lon,lat),
  sigma ~ x1 + x2 + x3 + lon + lat +
    s(x1) + s(x2) + s(x3) + s(lon) + s(lat) + te(lon,lat)
)

b <- bamlss(f, data = d, optimizer = boostm,
  sampler = FALSE, scale.d = TRUE, nu = 0.05,
  maxit = 400, stop.criterion = "AIC", force.stop = FALSE)

## Plot estimated effects.
plot(b)

## Plot AIC and log-lik contributions.
boost_plot(b, "AIC")
boost_plot(b, "loglik.contrib")

## Out-of-sample selection of model terms.
set.seed(123)
d <- GAMart(n = 5000)

## Split data into training and testing
i <- sample(1:2, size = nrow(d), replace = TRUE)
d.test <- subset(d, i == 1)
d.train <- subset(d, i == 2)

## Model formula
f <- list(
  num ~ s(x1) + s(x2) + s(x3),
  sigma ~ s(x1) + s(x2) + s(x3)
)

## Create model frame for out-of-sample selection.
sm <- boost_frame(f, train = d.train, test = d.test, family = "gaussian")

## Out-of-sample selection function.
sfun <- function(parameters) {
  sm$parameters <- parameters
  p <- predict(sm, type = "parameter")
  -1 * sum(sm$family$d(d.test$num, p, log = TRUE))
}

## Start boosting with out-of-sample negative
## log-likelihood selection of model terms.
b <- bamlss(f, data = d.train, sampler = FALSE, optimizer = boost,
  selectfun = sfun, always = "best")

## Plot curve of negative out-of-sample log-likelihood.
boost_plot(b, which = "user")
}
}

\keyword{regression}

