\name{continue}
\alias{continue}

\title{
  Continue Sampling
}

\description{
  This function takes a \code{\link{bamlss}} object which was created using a sampler function
  and continues sampling from the last state of the MCMC chain.
}

\usage{
continue(object, cores = NULL, combine = TRUE,
  sleep = NULL, results = TRUE, ...)
}

\arguments{
  \item{object}{A \code{\link{bamlss}} object which contains samples.}
  \item{cores}{An integer specifying the number of cores that should be used for the sampler
    function. This is based on function \code{\link[parallel]{mclapply}} of the \pkg{parallel}
    package.}
  \item{combine}{Should the new samples be combined with the old samples into one
    \code{\link[coda]{mcmc}} matrix? Note that if \code{combine = FALSE} and the number of
    iterations differ from one sampling step to the other there will be an error since
    the start and end points in the samples \code{\link[coda]{mcmc}} objects are different!}
  \item{sleep}{Time the system should sleep before the next core is started.}
  \item{results}{If a results function was used to create the \code{\link{bamlss}} object,
    should the results function be applied using the new samples?}
  \item{\dots}{Arguments passed to the sampler function.}
}

\value{
A \code{\link{bamlss}} object.
}

\seealso{
\code{\link{bamlss}}
}

\examples{
\dontrun{## Simulate and run model with MCMC.
set.seed(123)
d <- GAMart()
b <- bamlss(num ~ s(x1) + s(x2) + s(x3) + te(lon,lat), data = d)

## Continue sampling.
a <- continue(b)

## Plot all samples.
plot(samples(a, combine = FALSE))
}
}

\keyword{regression}

