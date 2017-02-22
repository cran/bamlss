\name{plot.bamlss}
\alias{plot.bamlss}
\alias{plot.bamlss.results}

\title{
  Plotting BAMLSS
}

\description{
  Plotting methods for objects of class \code{"bamlss"} and \code{"bamlss.results"}, which can be
  used for producing effect plots of model terms, trace plots of samples or residual plots.
  Note that effect plots of model terms with more than two covariates are not supported, for this
  purpose use function \code{\link{predict.bamlss}}.
}

\usage{
\method{plot}{bamlss}(x, model = NULL, term = NULL,
  which = "effects", parameters = FALSE,
  ask = dev.interactive(), spar = TRUE, ...)

\method{plot}{bamlss.results}(x, model = NULL, term = NULL,
  ask = dev.interactive(), scale = 1, spar = TRUE, ...)
}

\arguments{
  \item{x}{An object of class \code{"bamlss"} or \code{"bamlss.results"}.}
  \item{model}{Character or integer. For which model should the plots be created?}
  \item{term}{Character or integer. For which model term should a plot be created?}
  \item{which}{Character or integer, selects the type of plot: \code{"effects"} produces
    effect plots of smooth model terms, \code{"samples"} shows trace plots of samples,
    \code{"hist-resid"} shows a histogram of residuals (see also
    \code{\link{residuals.bamlss}} for the different available types), \code{"qq-resid"} shows
    a quantile-quantile plot of residuals, \code{"scatter-resid"} shows a scatter plot of
    residuals with fitted values for the distribution mean (if available in the family object),
    \code{"max-acf"} shows an \code{\link[stats]{acf}} plot of the maximum autocorrelation of
    all parameter samples.}
  \item{parameters}{For trace plots of parameters, should corresponding parameter values as returned
    from an optimizer function (e.g., \code{\link{bfit}}) be added as horizontal lines?}
  \item{ask}{For multiple plots, the user is asked to show the next plot.}
  \item{scale}{If set to \code{1}, effect plots all have the same scale on the y-axis. If set to
    \code{0} each effect plot has its own scale for the y-axis.}
  \item{spar}{Should graphical parameters be set?}
  \item{\dots}{Arguments to be passed to \code{\link{plot2d}}, \code{\link{plot3d}},
    \code{\link{sliceplot}}, \code{\link{plotblock}}, \code{\link{plotmap}} and
    \code{\link{residuals.bamlss}}.}
}

\seealso{
  \code{\link{bamlss}}, \code{\link{results.bamlss.default}}, \code{\link{residuals.bamlss}}.
}

\examples{
\dontrun{## Generate some data.
d <- GAMart()

## Model formula.
f <- list(
  num ~ s(x1) + s(x2) + s(x3) + te(lon,lat),
  sigma ~ s(x2) + te(lon,lat)
)

## Estimate model.
b <- bamlss(f, data = d)

## Effect plots
plot(b, ask = FALSE)
plot(b, model = "mu")
plot(b, model = "sigma", term = "te(lon,lat)")

## Trace plots.
plot(b, which = "samples")

## Residual plots.
plot(b, which = 3:4)
}
}

\keyword{aplot}
