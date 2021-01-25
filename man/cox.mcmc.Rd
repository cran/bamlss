\name{sam_Cox}
\alias{sam_Cox}
\alias{cox_mcmc}

\title{
  Cox Model Markov Chain Monte Carlo
}

\description{
  This sampler function implements a derivative based MCMC algorithm for
  flexible Cox models with structured additive predictors.
}

\usage{
sam_Cox(x, y, family, start, weights, offset,
  n.iter = 1200, burnin = 200, thin = 1,
  verbose = TRUE, digits = 4, step = 20, ...)

cox_mcmc(x, y, family, start, weights, offset,
  n.iter = 1200, burnin = 200, thin = 1,
  verbose = TRUE, digits = 4, step = 20, ...)
}

\arguments{
 \item{x}{The \code{x} list, as returned from function
    \code{\link{bamlss.frame}} and transformed by function \code{\link{surv_transform}},
    holding all model matrices and other information that is used for
    fitting the model.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.
    In this case this is the \code{\link{cox_bamlss}} family object.}
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
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{step}{How many times should algorithm runtime information be printed, divides \code{n.iter}.}
  \item{\dots}{Currently not used.}
}

\details{
  The sampler uses derivative based proposal functions to create samples of parameters.
  For time-dependent functions the proposals are based on one Newton-Raphson iteration centered
  at the last state, while for the time-constant functions proposals can be based
  on iteratively reweighted least squares (IWLS), see also function \code{\link{GMCMC}}.
  The integrals that are part of the time-dependent function updates are solved numerically.
  In addition, smoothing variances are sampled using slice sampling.
}

\value{
  The function returns samples of parameters. The samples are provided as a
  \code{\link[coda]{mcmc}} matrix.
}

\references{
  Umlauf N, Klein N, Zeileis A (2016). Bayesian Additive Models for Location
    Scale and Shape (and Beyond). \emph{(to appear)}
}

\seealso{
  \code{\link{opt_Cox}}, \code{\link{cox_bamlss}}, \code{\link{surv_transform}},
  \code{\link{simSurv}}, \code{\link{bamlss}}
}

\examples{
\dontrun{library("survival")
set.seed(123)

## Simulate survival data.
d <- simSurv(n = 500)

## Formula of the survival model, note
## that the baseline is given in the first formula by s(time).
f <- list(
  Surv(time, event) ~ s(time) + s(time, by = x3),
  gamma ~ s(x1) + s(x2)
)

## Cox model with continuous time.
## Note the the family object cox_bamlss() sets
## the default optimizer and sampler function!
## First, posterior mode estimates are computed
## using function opt_Cox(), afterwards the
## sampler sam_Cox() is started.
b <- bamlss(f, family = "cox", data = d)

## Plot estimated effects.
plot(b)
}
}

\keyword{regression}
\keyword{survival}

