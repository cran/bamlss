\name{sam_GMCMC}
\alias{sam_GMCMC}
\alias{GMCMC}
\alias{GMCMC_iwls}
\alias{GMCMC_iwlsC}
\alias{GMCMC_iwlsC_gp}
\alias{GMCMC_slice}

\title{
  General Markov Chain Monte Carlo for BAMLSS
}

\description{
  These functions provide a quite general infrastructure for sampling BAMLSS. The default
  proposal function is based on iteratively weighted least squares (IWLS), however,
  each model term may have a different updating function, see the details.
}

\usage{
## Sampler functions:
sam_GMCMC(x, y, family, start = NULL, weights = NULL, offset = NULL,
  n.iter = 1200, burnin = 200, thin = 1, verbose = TRUE,
  step = 20, propose = "iwlsC_gp", chains = NULL, ...)

GMCMC(x, y, family, start = NULL, weights = NULL, offset = NULL,
  n.iter = 1200, burnin = 200, thin = 1, verbose = TRUE,
  step = 20, propose = "iwlsC_gp", chains = NULL, ...)

## Propose functions:
GMCMC_iwls(family, theta, id, eta, y, data,
  weights = NULL, offset = NULL, ...)
GMCMC_iwlsC(family, theta, id, eta, y, data,
  weights = NULL, offset = NULL, zworking, resids, rho, ...)
GMCMC_iwlsC_gp(family, theta, id, eta, y, data,
  weights = NULL, offset = NULL, zworking, resids, rho, ...)
GMCMC_slice(family, theta, id, eta, y, data, ...)
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
  \item{n.iter}{Sets the number of MCMC iterations.}
  \item{burnin}{Sets the burn-in phase of the sampler, i.e., the number of starting samples that
    should be removed.}
  \item{thin}{Defines the thinning parameter for MCMC simulation. E.g., \code{thin = 10} means,
    that only every 10th sampled parameter will be stored.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{step}{How many times should algorithm runtime information be printed, divides \code{n.iter}.}
  \item{propose}{Sets the propose function for model terms, e.g. for a term \code{s(x)} in the model
    formula. Per default this is set to \code{"iwlsC"}, a character pointing to the set of propose
    functions, see above. Other options are \code{"iwls"} and \code{"slice"}, however, this is
    more experimental and should not be set by the user. Another option is to pass a full propose
    function which should be used for each model term, the structure of propose functions is
    described in the details below. Model terms may also have different propose functions, see
    the example section.}
  \item{chains}{How many chains should be started? Chains a sampled sequentially!}
  \item{theta}{The current state of parameters, provided as a named list. The first level
    represents the parameters of the distribution, the second level the parameters of
    the model terms. E.g., using the \code{\link{gaussian_bamlss}} family object
    \code{theta[["mu"]][["s(x)"]]} extracts the current state of a model term \code{"s(x)"}
    of the \code{"mu"} parameter. Extraction is done with the \code{id} argument.}
  \item{id}{The parameter identifier, a character vector of length 2. The first character
    specifies the current distributional parameter, the second the current model term.}
  \item{eta}{The current value of the predictors, provided as a named list, one list entry
    for each parameter. The names correspond to the parameter names in the family object,
    see \code{\link{family.bamlss}}. E.g., when using the \code{\link{gaussian_bamlss}} family
    object, the current values for the mean can be extracted by \code{eta\$mu} and for the
    standard deviation by \code{eta\$sigma}.}
  \item{data}{An object as returned from function
    \code{\link[mgcv]{smooth.construct}} or \code{\link[mgcv]{smoothCon}}. The object is
    preprocessed by function \code{\link{bamlss.engine.setup}}.}
  \item{zworking}{Preinitialized numeric vector of \code{length(y)}, only for internal usage.}
  \item{resids}{Preinitialized numeric vector of \code{length(y)}, only for internal usage.}
  \item{rho}{An environment, only for internal usage.}
  \item{\dots}{Arguments passed to function \code{\link{bamlss.engine.setup}} and to the
    propose functions.}
}

\details{
  The sampler function \code{sam_GMCMC()} cycles through all distributional parameters and corresponding
  model terms in each iteration of the MCMC chain. Samples of the parameters of a model term (e.g.,
  \code{s(x)}) are generated by proposal functions, e.g. \code{GMCMC_iwls()}.

  The default proposal function that should be used for all model terms is set with argument
  \code{propose}. For smooth terms, e.g. terms created with function \code{\link[mgcv]{s}}, if
  a valid propose function is supplied within the extra \code{xt} list, this propose function will
  be used. This way each model term may have its own propose function for creating samples
  of the parameters. See the example section.

  The default proposal function \code{GMCMC_iwlsC_gp} allows for general priors for the
  smoothing variances and general penalty functions. Samples of smoothing variances are computed
  using slice sampling. Function \code{GMCMC_iwlsC} samples smoothing variances of univariate terms
  assuming an inverse gamma prior. Terms of higher dimensions use again slice sampling for creating
  samples of smoothing variances.

  Function \code{GMCMC_iwls} is similar to function \code{GMCMC_iwlsC} but uses plain \R code.

  Function \code{GMCMC_slice} applies slice sampling also for the regression coefficients and
  is therefore relatively slow.
}

\value{
  The function returns samples of parameters, depending on the return value of the propose functions
  other quantities can be returned. The samples are provided as a \code{\link[coda]{mcmc}} matrix.
  If \code{chains > 1}, the samples are provided as a \code{\link[coda]{mcmc.list}}.
}


\references{
  Umlauf N, Klein N, Zeileis A (2016). Bayesian Additive Models for Location
    Scale and Shape (and Beyond). \emph{(to appear)}
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}},
  \code{\link{bamlss.engine.setup}}, \code{\link{set.starting.values}}, \code{\link{s2}}
}

\examples{
\dontrun{## Simulated data example illustrating
## how to call the sampler function.
## This is done internally within
## the setup of function bamlss().
d <- GAMart()
f <- num ~ s(x1, bs = "ps")
bf <- bamlss.frame(f, data = d, family = "gaussian")

## First, find starting values with optimizer.
opt <- with(bf, bfit(x, y, family))

## Sample.
samps <- with(bf, sam_GMCMC(x, y, family, start = opt$parameters))
plot(samps)
}
}

\keyword{regression}

