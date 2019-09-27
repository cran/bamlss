\name{bfit}
\alias{bfit}
\alias{bfit_iwls}
\alias{bfit_iwls_Matrix}
\alias{bfit_iwls_lm}
\alias{bfit_iwls_optim}
\alias{bfit_lm}
\alias{bfit_optim}
\alias{bfit_glmnet}

\title{
  Fit BAMLSS with Backfitting
}

\description{
  This optimizer function is a generic tool for fitting BAMLSS using a backfitting algorithm.
  The backfitting procedure is based on iteratively weighted least squares (IWLS) for finding
  posterior mode estimates, however, the updating methods for model terms can be more general, see
  the details section. In addition, the default IWLS updating scheme implements optimum smoothing
  variance selection based on information criteria using a stepwise approach.
}

\usage{
## Optimizer function:
bfit(x, y, family, start = NULL, weights = NULL, offset = NULL,
  update = "iwls", criterion = c("AICc", "BIC", "AIC"),
  eps = .Machine$double.eps^0.25, maxit = 400,
  outer = NULL, inner = FALSE, mgcv = FALSE,
  verbose = TRUE, digits = 4, flush = TRUE,
  nu = NULL, stop.nu = NULL, ...)

## Model term updating functions:
bfit_iwls(x, family, y, eta, id, weights, criterion, ...)
bfit_iwls_Matrix(x, family, y, eta, id, weights, criterion, ...)
bfit_lm(x, family, y, eta, id, weights, criterion, ...)
bfit_optim(x, family, y, eta, id, weights, criterion, ...)
bfit_glmnet(x, family, y, eta, id, weights, criterion, ...)
}

\arguments{
  \item{x}{For function \code{bfit()} the \code{x} list, as returned from function
    \code{\link{bamlss.frame}}, holding all model matrices and other information that is used for
    fitting the model. For the updating functions an object as returned from function
    \code{\link[mgcv]{smooth.construct}} or \code{\link[mgcv]{smoothCon}}.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{weights}{Prior weights on the data, as returned from function \code{\link{bamlss.frame}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{update}{Sets the updating function for model terms, e.g. for a term \code{s(x)} in the model
    formula. Per default this is set to \code{"iwls"}, a character pointing to the set of updating
    functions, see above. Other options are \code{"optim"} and \code{"lm"} etc., however, this is
    more experimental and should not be set by the user. Another option is to pass a full updating
    function which should be used for each model term, the structure of updating functions is
    described in the details below. Model terms may also have different updating functions, see
    the example section implementing a new model term constructor for Gompertz growth curves
    using this feature.}
  \item{criterion}{Set the information criterion that should be used, e.g., for smoothing
    variance selection. Options are the corrected AIC \code{"AICc"}, the \code{"BIC"} and
    \code{"AIC"}.}
  \item{eps}{The relative convergence tolerance of the backfitting algorithm.}
  \item{maxit}{The maximum number of iterations for the backfitting algorithm}
  \item{outer}{Should the current working observations and weights be computed in one outer
    iteration, otherwise the working observations are computed anew for each model term
    updating step. The default will run one outer iteration first, afterwards model weights
    are computed for each term anew.}
  \item{inner}{Should the model terms for one parameter of the modeled distribution be
    fully updated until convergence in an inner iteration, i.e., the algorithm waits until
    coefficients for the current distribution parameter do not change anymore before updating
    the next parameter.}
  \item{mgcv}{Should the \pkg{mgcv} \code{\link[mgcv]{gam}} function be used for computing updates
    in an \code{inner} iteration with working observations provided in an \code{outer} iteration.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{flush}{use \code{\link{flush.console}} for displaying the current output in the console.}
  \item{nu}{Numeric or logical. Function \code{bfit} uses step length optimization of parameters when 
    updating a model term, useful when the algorithm encounters convergence problems. If
    \code{nu = TRUE} the step length parameter is optimized for each model term in each iteration of
    the backfitting algorithm. If \code{nu} is numeric, e.g. \code{nu = 1}, then \code{nu} is
    halfed until an improvement in the log-posterior is obtained or nu is smaller than
    \code{.Machine$double.eps}.}
  \item{stop.nu}{Integer. Should step length reduction be stopped after \code{stop.nu} iterations
    of the backfitting algorithm?}
  \item{eta}{The current value of the predictors, provided as a named list, one list entry
    for each parameter. The names correspond to the parameter names in the family object,
    see \code{\link{family.bamlss}}. E.g., when using the \code{\link{gaussian_bamlss}} family
    object, the current values for the mean can be extracted by \code{eta\$mu} and for the
    standard deviation by \code{eta\$sigma}.}
  \item{id}{Character, the name of the current parameter for which the model term should
    be updated.}
  \item{\dots}{For function \code{bfit()}, arguments passed to function
    \code{\link{bamlss.engine.setup}}. For updating functions, within the dots argument the
    actual \code{iteration} number of the backfitting algorithm, the actual total number of
    equivalent degrees of freedom \code{edf} and vectors \code{z} and \code{hess}
    only if argument \code{outer = TRUE} are provided.}
}

\details{
  This algorithm is based on iteratively weighted least squares (IWLS) for BAMLSS, i.e.,
  a Newton-Raphson or Fisher scoring algorithm is applied, similar to
  Rigby~and~Stasinopoulos~(2005). The algorithm utilizes the chain rule for computing derivatives of
  the log-posterior w.r.t. regression coefficients, therefore, to compute the working observations
  and weights only the derivatives of the log-likelihood w.r.t. the predictors are required.

  It is assumed that the provided \code{family} object holds functions for computing the first
  and second order derivatives of the log-likelihood w.r.t. the predictors. These Functions
  are provided within the named lists \code{"score"} and \code{"hess"} within the \code{family}
  object. See the documentation of \code{\link{family.bamlss}} and the code of the provided
  families, e.g. \code{\link{gaussian_bamlss}}, for examples of the required structure.

  The algorithm either updates each model term over all distributional parameters sequentially,
  or does a full update until convergence for model terms for one distributional parameter before
  updating the next parameter, see argument \code{inner}. Additionally, working observations and
  weights can be computed only once in an \code{outer} iteration.

  Starting values of regression coefficients and smoothing variances can be supplied, moreover,
  if a family object holds functions for initializing the distributional parameters, see also
  \code{\link{family.bamlss}}, starting values are based on the initialize functions.

  The default updating function for model terms is based on IWLS, which is assigned by function
  \code{\link{bamlss.engine.setup}}, however, special updating functions can be used.
  This is achieved by providing an updating function to argument
  \code{update}, which should be used for all model terms. Another option is to set the updating
  function within the \code{xt} argument of the \pkg{mgcv} smooth term constructor functions, see
  e.g. function \code{\link[mgcv]{s}}. If the \code{xt} list then holds an element named \code{"update"},
  which is a valid updating function, this updating function is used for the corresponding model
  term. This way it is possible to call different (special) updating functions for specific terms,
  e.g., that do not fit in the IWLS scheme. See the examples below. Note that this does not work if
  \code{mgcv = TRUE}, since the \code{\link[mgcv]{gam}} function assumes a strict linear
  representation of smooth terms.
  
  A model term updating function has the following arguments:

  \code{update(x, family, y, eta, id, weights, criterion, ...)}

  Here \code{x} is an object as returned from function \code{\link[mgcv]{smooth.construct}}
  or \code{\link[mgcv]{smoothCon}}. The \code{x} object is preprocessed by function
  \code{\link{bamlss.engine.setup}}, i.e., an element called \code{"state"} is assigned. The state
  element represents the current state of the model term holding the current values of the
  parameters with corresponding fitted values, as well as equivalent degrees of freedom, see
  also the values that are returned by such functions below. The backfitting algorithm uses the
  state of a model term for generating updates of the parameters. Note that for special model
  terms the state list should already be provided within the call to the corresponding
  smooth constructor function, see the growth curve example below.

  In addition, for special model terms the fitted values may not be computed by a linear combination
  of the design matrix and the coefficients. Therefore, the \code{x} object should hold an element
  named \code{"fit.fun"} which is a function for computing the fitted values.
  See also \code{\link{smooth.construct.bamlss.frame}} and \code{\link{predict.bamlss}} that use
  this setup. The arguments of fitting functions are

  \code{fit.fun(X, b, ...)}

  where \code{X} is the design matrix and \code{b} is the vector of coefficients. Hence, for
  usual IWLS updating the fitted values are computed by \code{X \%*\% b}. For special terms like
  nonlinear growth curves this may not be the case, see the example below. The fitting functions
  are assigned by \code{\link{bamlss.engine.setup}}, unless the function is already provided
  after calling the constructor function \code{\link[mgcv]{smooth.construct}} or
  \code{\link[mgcv]{smoothCon}}. Note that the dots argument is usually not needed by the user.

  The default updating function is \code{bfit_iwls()}. Function \code{bfit_iwls_Matrix()} uses the
  sparse matrix infrastructures of package \pkg{Matrix}. The \pkg{Matrix} package and
  \code{bfit_iwls_Matrix()} is used for model terms where
  the maximum number of non-zero entries in the design matrix is less than half of the total number
  of columns, if an additional argument \code{force.Matrix} is set to \code{TRUE} in the \code{bfit()}
  call.

  The IWLS updating functions find optimum smoothing variances according to an information criterion
  using a stepwise approach, i.e., in each iteration and for each model term update the updating
  functions try to find a better smoothing variance to control the trade-off between over-smoothing
  and nonlinear functional estimation. The search interval is centered around the current state
  of the smoothing variances, hence, in each iteration only a slight improvement is achieved.
  This algorithm is based on Belitz~and~Lang~(2008) and can also be viewed as a boosting approach
  for optimization.
}

\value{
  For function \code{bfit()} a list containing the following objects:
  \item{fitted.values}{A named list of the fitted values of the modeled parameters
    of the selected distribution.}
  \item{parameters}{The estimated set regression coefficients and smoothing variances.}
  \item{edf}{The equivalent degrees of freedom used to fit the model.}
  \item{logLik}{The value of the log-likelihood.}
  \item{logPost}{The value of the log-posterior.}
  \item{IC}{The value of the information criterion.}
  \item{converged}{Logical, indicating convergence of the backfitting algorithm.}
  \item{For updating functions a list providing the current state}{}
  \item{fitted.values}{The resulting fitted values after updating.}
  \item{parameters}{The resulting named numeric vector of updated model term parameters.
    Coefficients should be named with \code{"b1"}, \dots, \code{"bk"}, where \code{k} is
    the total number of coefficients. Smoothing variances should be named with
    \code{"tau21"}, \dots, \code{"tau2m"}, where \code{m} is the total number of smoothing
    variances assigned to the model term.}
  \item{edf}{The equivalent degrees of freedom used to produce the fitted values.}
  \item{hessian}{Optional, the coefficient Hessian information}
  \item{log.prior}{Optional, the value of the log-prior of the model term.}
}

\references{
  Belitz C, Lang S (2008). Simultaneous Selection of Variables and Smoothing Parameters in
  Structured Additive Regression Models. \emph{Computational Statistics \& Data Analysis},
  \bold{53}, pp 61-81.

  Umlauf N, Klein N, Zeileis A (2016). Bayesian Additive Models for Location
    Scale and Shape (and Beyond). \emph{(to appear)}

  Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location,
  scale and shape, (with discussion), \emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}},
  \code{\link{bamlss.engine.setup}}, \code{\link{set.starting.values}}, \code{\link{s2}}
}

\examples{
\dontrun{## Simulated data example illustrating
## how to call the optimizer function.
## This is done internally within
## the setup of function bamlss().
d <- GAMart(n = 200)
f <- num ~ s(x1) + s(x2) + s(x3)
bf <- bamlss.frame(f, data = d, family = "gaussian")
opt <- with(bf, bfit(x, y, family))
print(str(opt))

## Same with bamlss().
b <- bamlss(f, data = d, family = "gaussian", sampler = FALSE)
plot(b)
summary(b)

## Use of different updating function.
b <- bamlss(f, data = d, family = "gaussian",
  sampler = FALSE, update = bfit_lm)
plot(b)

## Use mgcv gam() function for updating.
b <- bamlss(f, data = d, family = "gaussian",
  sampler = FALSE, mgcv = TRUE)
plot(b)

## Special smooth constructor including updating/sampler
## function for nonlinear Gompertz curves.
## Note: element special.npar is needed here since this
##       function has 3 parameters but the design matrix only
##       one column!
smooth.construct.gc.smooth.spec <- function(object, data, knots) 
{
  object$X <- matrix(as.numeric(data[[object$term]]), ncol = 1)
  center <- if(!is.null(object$xt$center)) {
    object$xt$center
  } else TRUE
  object$by.done <- TRUE
  if(object$by != "NA")
    stop("by variables not supported!")
  object$fit.fun <- function(X, b, ...) {
    f <- b[1] * exp(-b[2] * exp(-b[3] * drop(X)))
    if(center)
      f <- f - mean(f)
    f
  }
  object$update <- bfit_optim
  object$propose <- GMCMC_slice
  object$prior <- function(b) { sum(dnorm(b, sd = 1000, log = TRUE)) }
  object$fixed <- TRUE
  object$state$parameters <- c("b1" = 0, "b2" = 0.5, "b3" = 0.1)
  object$state$fitted.values <- rep(0, length(object$X))
  object$state$edf <- 3
  object$special.npar <- 3 ## Important!
  class(object) <- c("gc.smooth", "no.mgcv", "special")
  object
}

## Work around for the "prediction matrix" of a growth curve.
Predict.matrix.gc.smooth <- function(object, data, knots) 
{
  X <- matrix(as.numeric(data[[object$term]]), ncol = 1)
  X
}

## Heteroscedastic growth curve data example.
set.seed(111)

d <- data.frame("time" = 1:30)
d$y <- 2 + 1 / (1 + exp(0.5 * (15 - d$time))) +
  rnorm(30, sd = exp(-3 + 2 * cos(d$time/30 * 6 - 3)))

## Special model terms must be called with s2()!
f <- list(
  y ~ s2(time, bs = "gc"),
  sigma ~ s(time)
)

## Fit model with special model term.
b <- bamlss(f, data = d,
  optimizer = bfit, sampler = GMCMC)

## Plot the fitted curves.
plot(b)

## Predict with special model term.
nd <- data.frame("time" = seq(1, 30, length = 100))
p <- predict(b, newdata = nd, model = "mu", FUN = c95)
plot(d, ylim = range(c(d$y, p)))
matplot(nd$time, p, type = "l",
  lty = c(2, 1, 2), col = "black", add = TRUE)
}
}

\keyword{regression}

