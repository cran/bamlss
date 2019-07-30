\name{bbfit}
\alias{bbfit}

\title{Batch-Wise Backfitting}

\description{
  Batch-wise backfitting estimation engine for GAMLSS using very large data sets.
}

\usage{
bbfit(x, y, family, shuffle = TRUE, start = NULL, offset = NULL,
  epochs = 1, nbatch = 10, maxit = Inf, verbose = TRUE, ...)
}

\arguments{
   \item{x}{For function \code{bfit()} the \code{x} list, as returned from function
    \code{\link{bamlss.frame}}, holding all model matrices and other information that is used for
    fitting the model. For the updating functions an object as returned from function
    \code{\link[mgcv]{smooth.construct}} or \code{\link[mgcv]{smoothCon}}.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.}
  \item{shuffle}{Should observations be shuffled?}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{epochs}{For how many epochs should the algorithm run?}
  \item{nbatch}{Number of batches.}
  \item{maxit}{Maximum number of iterations of the inner loop.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{\dots}{Not used.}
}

\details{
  The algorithm uses batch-wise estimation of smoothing variances, which are estimated on an
  hold-out batch. This way, models for very large data sets can be estimated. Note, the algorithm
  only works in combination withe the \pkg{ff} and \pkg{ffbase} package. The data needs to be stored
  as comma separated file on disc, see the example.
}

\value{
  For function \code{bbfit()} a list containing the following objects:
  \item{fitted.values}{A named list of the fitted values of the modeled parameters
    of the selected distribution.}
  \item{parameters}{The estimated set regression coefficients and smoothing variances.}
  \item{shuffle}{Logical}
  \item{runtime}{The runtime of the algorithm.}
}

\seealso{
\code{\link{bamlss}}, \code{\link{bfit}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart(n = 27000, sd = -1)

## Write data to disc.
tf <- tempdir()
write.table(d, file.path(tf, "d.raw"), quote = FALSE, row.names = FALSE, sep = ",")

## Estimation using batch-wise backfitting.
f <- list(
  num ~ s(x1,k=40) + s(x2,k=40) + s(x3,k=40) + te(lon,lat,k=10),
  sigma ~ s(x1,k=40) + s(x2,k=40) + s(x3,k=40) + te(lon,lat,k=10)
)

b <- bamlss(f, data = file.path(tf, "d.raw"), optimizer = bbfit,
  sampler = FALSE, nbatch = 10, epochs = 2, loglik = TRUE)

## Show estimated effects.
plot(b)
}
}

\keyword{regression}

