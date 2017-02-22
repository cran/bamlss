\name{results.bamlss.default}
\alias{results.bamlss.default}

\title{
  Compute BAMLSS Results for Plotting and Summaries
}

\description{
  The results function combines estimated parameters and/or samples with the
  \code{\link{bamlss.frame}} and computes the data that can be used, e.g., for
  creating effect plots or summary statistics. The function is usually used
  internally within \code{\link{bamlss}}. The object returned is of class "bamlss.results",
  which has a plotting method, see \code{\link{plot.bamlss.results}}.
}

\usage{
results.bamlss.default(x, what = c("samples", "parameters"),
  grid = -1, nsamps = NULL, burnin = NULL, thin = NULL, ...)
}

\arguments{
  \item{x}{A \code{\link{bamlss.frame}} which has estimated parameters
    or samples. See also \code{\link{bfit}} and \code{\link{GMCMC}}.}
  \item{what}{Should the results data be prepared using estimated parameters or samples?}
  \item{grid}{Integer, sets the number of grid points for univariate functions to be used for
    creating results data, e.g., for plotting. This is more efficient when using data sets with
    a large number of unique covariate values. If negative suitable defaults are chosen.}
  \item{nsamps}{Integer, if results are computed using parameter samples, this argument controls
    the number of samples that should be used, e.g., if \code{nsamps = 100} only 100 samples
    with equidistant intervals are selected. Basically similar to argument \code{thin}.}
  \item{burnin}{Integer, sets the number of samples that should be dropped from the beginning
    of the MCMC chain when creating results.}
  \item{thin}{Integer, should the MCMC chain be thinned additionally?}
  \item{\dots}{Currently not used.}
}

\value{
  An object of class \code{"bamlss.results"}.
}

\seealso{
  \code{\link{plot.bamlss.results}}, \code{\link{bamlss}}.
}

\examples{
\dontrun{## Simulate data.
d <- GAMart()

## Estimate model with no results.
b <- bamlss(num ~ s(x1) + s(x2) + s(x3),
  data = d, results = FALSE)

## Compute model results
a <- results.bamlss.default(b)

## Plot results for smooth terms.
plot(a)
}
}

\keyword{dplot}

