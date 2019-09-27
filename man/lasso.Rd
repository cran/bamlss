\name{la}
\alias{la}
\alias{lasso}
\alias{lasso_plot}
\alias{lasso_stop}
\alias{lasso_coef}
\alias{lasso_transform}

\title{Lasso Smooth Constructor}

\description{
  Smooth constructors and optimizer for Lasso penalization with \code{\link{bamlss}}. The
  penalization is based on a Taylor series approximation of the Lasso penalty.
}

\usage{
## Smooth constructor function.
la(formula, type = c("single", "multiple"), ...)

## Single Lasso smoothing parameter optimizer.
lasso(x, y, start = NULL, adaptive = TRUE, lower = 0.001, upper = 1000,
  nlambda = 100, lambda = NULL,  multiple = FALSE, verbose = TRUE,
  digits = 4, flush = TRUE, nu = NULL, stop.nu = NULL,
  ridge = .Machine$double.eps^0.5, zeromodel = NULL, ...)

## Lasso transformation function to set
## adaptive weights from an unpenalized model.
lasso_transform(x, zeromodel, nobs = NULL, ...)

## Plotting function for lasso() optimizer.
lasso_plot(x, which = c("criterion", "parameters"),
  spar = TRUE, model = NULL, name = NULL, mstop = NULL,
  retrans = FALSE, color = NULL, show.lambda = TRUE,
  labels = NULL, digits = 2, ...)

## Extract optimum stopping iteration for lasso() optimizer.
## Based on the minimum of the information criterion.
lasso_stop(x)

## Extract retransformed Lasso coefficients.
lasso_coef(x, ...)
}

\arguments{
  \item{formula}{A formula like \code{~ x1 + x2 + ... + xk} of variables which should be
    penalized with Lasso.}
  \item{type}{Should one single penalty parameter be used or multiple parameters, one for each
    covariate in \code{formula}.}
  \item{x}{For function \code{lasso()} and \code{lasso_transform()} the \code{x} list, as returned
    from function \code{\link{bamlss.frame}}, holding all model matrices and other information that
    is used for fitting the model. For the plotting function and
    \code{lasso_stop()}/\code{lasso_coef()} the
    corresponding \code{\link{bamlss}} object fitted with the \code{lasso()} optimizer.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{start}{A vector of starting values. Note, Lasso smoothing parameters will be dropped.}
  \item{adaptive}{Should adaptive weights be used for fused Lasso terms?}
  \item{lower}{Numeric. The minimum lambda value.}
  \item{upper}{Numeric. The maximum lambda value.}
  \item{nlambda}{Integer. The number of smoothing parameters for which coefficients should be
    estimated, i.e., the vector of smoothing parameters is build up as a sequence from
    \code{lower} to \code{upper} with length \code{nlambda}.}
  \item{lambda}{Numeric. A sequence/vector of lambda parameters that should be used.}
  \item{multiple}{Logical. Should the lambda grid be exapnded to search for multiple lambdas, one for
    each distributional parameter.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}. If the optimum lambda value
    is plotted, the number of decimal decimal places to be used within \code{lasso_plot()}.}
  \item{flush}{use \code{\link{flush.console}} for displaying the current output in the console.}
  \item{nu}{Numeric or logical. Defines the step length for parameter updating of a model term,
    useful when the algorithm encounters convergence problems. If \code{nu = TRUE} the step length
    parameter is optimized for each model term in each iteration of the backfitting algorithm.}
  \item{stop.nu}{Integer. Should step length reduction be stopped after \code{stop.nu} iterations
    of the Lasso algorithm?}
  \item{ridge}{A ridge penalty parameter that should be used when finding adaptive weights, i.e.,
    parameters from an unpenalized model. The ridge penalty is used to stabilize the estimation
    in complex models.}
  \item{zeromodel}{A model containing the unpenalized parameters, e.g., for each \code{la()}
    terms one can place a simple ridge penalty with \code{la(x, ridge = TRUE, sp = 0.1)}. This
    way it is possible to find the unpenalized parameters that can be used as adaptive
    weights for fusion penalties.}
  \item{nobs}{Integer, number of observations of the data used for modeling. If not supplied
    \code{nobs} is taken from the number of rows from the model term design matrices.}
  \item{which}{Which of the two provided plots should be created, character or integer \code{1} and \code{2}.}
  \item{spar}{Should graphical parameters be set by the plotting function?}
  \item{model}{Character selecting for which model the plot shpuld be created.}
  \item{name}{Character, the name of the coefficient group that should be plotted. Note that
    the string provided in \code{name} will be removed from the labels on the 4th axis.}
  \item{mstop}{Integer vector, defines the path length to be plotted.}
  \item{retrans}{Logical, should coefficients be re-transformed before plotting?}
  \item{color}{Colors or color function that creates colors for the group paths.}
  \item{show.lambda}{Logical. Should the optimum value of the penalty parameter lambda be shown?}
  \item{labels}{A character string of labels that should be used on the 4 axis.}
  \item{\dots}{Arguments passed to the subsequent smooth constructor function.
    \code{lambda} controls the starting value of the penalty parameter, \code{const} the constant
    that is added within the penalty approximation. Moreover, \code{fuse = 1} enforces nominal
    fusion of categorical variables and \code{fuse = 2} ordered fusion within \code{la()} Note that
    \code{la()} terms with and without fusion should not be mixed when using the \code{lasso()}
    optimizer function.
    For the optimizer \code{lasso()} arguments passed to function \code{\link{bfit}}.}
}

\value{
  For function \code{la()}, similar to function \code{\link[mgcv]{s}} a simple smooth
  specification object.

  For function \code{lasso()} a list containing the following objects:
  \item{fitted.values}{A named list of the fitted values based on the last lasso iteration
    of the modeled parameters of the selected distribution.}
  \item{parameters}{A matrix, each row corresponds to the parameter values of one boosting iteration.}
  \item{lasso.stats}{A matrix containing information about the log-likelihood, log-posterior
    and the information criterion for each lambda.}
}

\references{
  Andreas Groll, Julien Hambuckers, Thomas Kneib, and Nikolaus Umlauf (2019). Lasso-type penalization in
  the framework of generalized additive models for location, scale and shape.
  \emph{Computational Statistics \& Data Analysis}.
  \url{https://doi.org/10.1016/j.csda.2019.06.005}

  Oelker Margreth-Ruth and Tutz Gerhard (2015). A uniform framework for combination of
  penalties in generalized structured models. \emph{Adv Data Anal Classif}.
  \url{http://dx.doi.org/10.1007/s11634-015-0205-y}
}

\seealso{
\code{\link[mgcv]{s}}, \code{\link[mgcv]{smooth.construct}}
}

\examples{
\dontrun{## Simulated fusion Lasso example.
bmu <- c(0,0,0,2,2,2,4,4,4)
bsigma <- c(0,0,0,-2,-2,-2,-1,-1,-1)
id <- factor(sort(rep(1:length(bmu), length.out = 300)))

## Response.
set.seed(123)
y <- bmu[id] + rnorm(length(id), sd = exp(bsigma[id]))

## Estimate model:
## fuse=1 -> nominal fusion,
## fuse=2 -> ordinal fusion,
## first, unpenalized model to be used for adaptive fusion weights.
f <- list(y ~ la(id,fuse=2,fx=TRUE), sigma ~ la(id,fuse=1,fx=TRUE))
b0 <- bamlss(f, sampler = FALSE)

## Model with single lambda parameter.
f <- list(y ~ la(id,fuse=2), sigma ~ la(id,fuse=1))
b1 <- bamlss(f, sampler = FALSE, optimizer = lasso,
  criterion = "BIC", zeromodel = b0)

## Plot information criterion and coefficient paths.
lasso_plot(b1, which = 1)
lasso_plot(b1, which = 2)
lasso_plot(b1, which = 2, model = "mu", name = "mu.s.la(id).id")
lasso_plot(b1, which = 2, model = "sigma", name = "sigma.s.la(id).id")

## Extract coefficients for optimum Lasso parameter.
coef(b1, mstop = lasso_stop(b1))

## Predict with optimum Lasso parameter.
p1 <- predict(b1, mstop = lasso_stop(b1))

## Full MCMC, needs lasso_transform() to assign the
## adaptive weights from unpenalized model b0.
b2 <- bamlss(f, optimizer = FALSE, transform = lasso_transform,
  zeromodel = b0, nobs = length(y), start = coef(b1, mstop = lasso_stop(b1)),
  n.iter = 4000, burnin = 1000)
summary(b2)
plot(b2)

ci <- confint(b2, model = "mu", pterms = FALSE, sterms = TRUE)
lasso_plot(b1, which = 2, model = "mu", name = "mu.s.la(id).id", spar = FALSE)
for(i in 1:8) {
  abline(h = ci[i, 1], lty = 2, col = "red")
  abline(h = ci[i, 2], lty = 2, col = "red")
}
}
}

\keyword{regression}

