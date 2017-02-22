\name{JAGS}
\alias{JAGS}
\alias{BUGSeta}
\alias{BUGSmodel}

\title{
  Markov Chain Monte Carlo for BAMLSS using JAGS
}

\description{
  This sampler function for BAMLSS is an interface to the JAGS library
  using package \code{\link[rjags]{rjags}}. The function basically interprets the
  \code{\link{bamlss.frame}} into BUGS code, similar to the \code{\link[mgcv]{jagam}} function of
  package \code{\link[mgcv]{mgcv}}. I.e., the function uses the random effects representation of
  smooth terms, see the transformer function \code{\link{randomize}} to generate the BUGS code.

  Note that estimating BAMLSS with JAGS is not very efficient.
  Also note that this function is more experimental and support is only provided for a small
  number of \code{\link{bamlss.family}} objects.

  Function \code{BUGSeta()} therefore computes the code and data for one parameter of
  the modeled distribution. Function \code{BUGSmodel()} then collects all parameter model code
  and data, which can be send to JAGS.
}

\usage{
## Sampler function:
JAGS(x, y, family, start = NULL,
  tdir = NULL, n.chains = 1, n.adapt = 100,
  n.iter = 4000, thin = 2, burnin = 1000,
  seed = NULL, verbose = TRUE, set.inits = TRUE,
  save.all = FALSE, modules = NULL, ...)

## Function to interpret an additive predictor into BUGS code:
BUGSeta(x, id = NULL, ...)

## Function to interpret the full BAMLSS:
BUGSmodel(x, family, is.stan = FALSE, reference = NULL, ...)
}

\arguments{
  \item{x}{For function \code{JAGS()} and \code{BUGSmodel()} the \code{x} list, as returned from
    function \code{\link{bamlss.frame}}, holding all model matrices and other information that is
    used for fitting the model. For function \code{BUGSeta()} argument \code{x} is one element of
    the \code{x} object, i.e., one parameter.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{tdir}{The path to the temporary directory that should be used.}
  \item{n.chains}{Specifies the number of sequential MCMC chains that should be run with JAGS.}
  \item{n.adapt}{Specifies the number of iterations that should be used as an initial adaptive phase.}
  \item{n.iter}{Sets the number of MCMC iterations.}
  \item{thin}{Defines the thinning parameter for MCMC simulation. E.g., \code{thin = 10} means,
    that only every 10th sampled parameter will be stored.}
  \item{burnin}{Sets the burn-in phase of the sampler, i.e., the number of starting samples that
    should be removed.}
  \item{seed}{Sets the seed.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{set.inits}{Should initial values of BAMLSS \code{\link{parameters}} be provided to JAGS,
    if available. Set in argument \code{start}.}
  \item{save.all}{Should all JAGS files be saved in \code{tdir}.}
  \item{modules}{Specify additional modules that should be loaded, see function \code{\link[rjags]{load.module}}.}
  \item{id}{Character, the current parameter name for which the BUGS code should be produced.}
  \item{is.stan}{Should the BUGS code be translated to STAN code. Note that this is only experimental.}
  \item{reference}{A \code{character} specifying a reference category, e.g., when
    fitting a multinomial model.}
  \item{\dots}{Currently not used.}
}

\value{
  Function \code{JAGS()} returns samples of parameters. The samples are provided as a \code{\link[coda]{mcmc}}
  matrix. If \code{n.chains > 1}, the samples are provided as a \code{\link[coda]{mcmc.list}}.

  Function \code{BUGSeta()} returns the BUGS model code and preprocessed data for one
  additive predictor. Function \code{BUGSmodel()} then combines all single BUGS code chunks and
  the data and creates the final BUGS model code that can be send to JAGS.
}

\note{
  Note that for setting up a new family object to be used with \code{JAGS()} additional
  information needs to be supplied. The extra information must be placed within the
  family object in an element named \code{"bugs"}. The following entries should be supplied
  within the \code{..$bugs} list:
  \itemize{
    \item \code{"dist"}. The name of the distribution in BUGS/JAGS model language.
    \item \code{"eta"}. The function that computes the BUGS code for one structured additive predictor.
      Function \code{BUGSeta()} is used per default.
    \item \code{"model"}. The function that merges all single predictor BUGS model code and data.
      The default function is \code{BUGSmodel()}.
    \item \code{"reparam"}. A named vector of character strings that specify a re-parametrization.
  }
  See also the example code of \code{\link{family.bamlss}}.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}},
  \code{\link{bamlss.engine.setup}}, \code{\link{set.starting.values}}, \code{\link{bfit}},
  \code{\link{GMCMC}}
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

## Sample with JAGS.
if(require("rjags")) {
  samps <- with(bf, JAGS(x, y, family, start = opt$parameters))
  plot(samps)
}
}
}

\keyword{regression}

