\name{MVNORM}
\alias{MVNORM}

\title{
  Create Samples for BAMLSS by Multivariate Normal Approximation
}

\description{
  This sampler function for BAMLSS uses estimated \code{\link{parameters}} and the Hessian
  information to create samples from a multivariate normal distribution. Note that smoothing
  variance uncertainty is not accounted for, therefore, the resulting credible intervals
  are most likely too narrow. 
}

\usage{
MVNORM(x, y = NULL, family = NULL, start = NULL,
  n.samples = 500, hessian = NULL, ...)
}

\arguments{
  \item{x}{The \code{x} list, as returned from
    function \code{\link{bamlss.frame}}, holding all model matrices and other information that is
    used for fitting the model. Or an object returned from function \code{\link{bamlss}}.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{n.samples}{Sets the number of samples that should be generated.}
  \item{hessian}{The Hessian matrix that should be used. Note that the row and column names
    must be the same as the names of the \code{\link{parameters}}. If \code{hessian = NULL}
    the function uses \code{\link{optim}} to compute the Hessian if it is not provided
    within \code{x}.}
  \item{\dots}{Arguments passed to function \code{\link{optim}}.}
}

\value{
  Function \code{MVNORM()} returns samples of parameters. The samples are provided as a
  \code{\link[coda]{mcmc}} matrix.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}},
  \code{\link{bamlss.engine.setup}}, \code{\link{set.starting.values}}, \code{\link{bfit}},
  \code{\link{GMCMC}}
}

\examples{
## Simulated data example illustrating
## how to call the sampler function.
## This is done internally within
## the setup of function bamlss().
d <- GAMart()
f <- num ~ s(x1, bs = "ps")
bf <- bamlss.frame(f, data = d, family = "gaussian")

## First, find starting values with optimizer.
o <- with(bf, bfit(x, y, family))

## Sample.
samps <- with(bf, MVNORM(x, y, family, start = o$parameters))
plot(samps)
}

\keyword{regression}

