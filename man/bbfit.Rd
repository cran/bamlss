\name{bbfit}
\alias{bbfit}
\alias{bbfitp}
\alias{contribplot}

\title{Batchwise Backfitting}

\description{
  Batchwise backfitting estimation engine for GAMLSS using very large data sets.
}

\usage{
## Batchwise backfitting engine.
bbfit(x, y, family, shuffle = TRUE, start = NULL, offset = NULL,
  epochs = 1, nbatch = 10, verbose = TRUE, ...)

## Parallel version.
bbfitp(x, y, family, mc.cores = 1, ...)

## Loglik contribution plot.
contribplot(x, ...)
}

\arguments{
   \item{x}{For function \code{bfit()} the \code{x} list, as returned from function
    \code{\link{bamlss.frame}}, holding all model matrices and other information that is used for
    fitting the model. For the updating functions an object as returned from function
    \code{\link[mgcv]{smooth.construct}} or \code{\link[mgcv]{smoothCon}}. For function
    \code{contribplot()}, a \code{"bamlss"} object using \code{bbfit()} with argument
    \code{select = TRUE}.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.}
  \item{shuffle}{Should observations be shuffled?}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{epochs}{For how many epochs should the algorithm run?}
  \item{nbatch}{Number of batches. Can also be a number between 0 and 1, i.e., determining
    the fraction of observations that should be used for fitting.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{mc.cores}{On how many cores should estimation be started?}
  \item{\dots}{For \code{bbfitp()} all arguments to be passed to \code{bbfit()}.}
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

