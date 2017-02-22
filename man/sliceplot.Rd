\name{sliceplot}
\alias{sliceplot}

\title{
  Plot Slices of Bivariate Functions 
}

\description{
  This function plots slices from user defined values of bivariate surfaces.
}
 
\usage{
sliceplot(x, y = NULL, z = NULL, view = 1, c.select = NULL,
  values = NULL, probs = c(0.1, 0.5, 0.9), grid = 100,
  legend = TRUE, pos = "topright", digits = 2, data = NULL,
  rawdata = FALSE, type = "mba", linear = FALSE,
  extrap = FALSE, k = 40, rug = TRUE, rug.col = NULL,
  jitter = TRUE, ...)
}

\arguments{
  \item{x}{A matrix or data frame, containing the covariates for which the effect should be plotted 
    in the first and second column and at least a third column containing the effect. Another
    possibility is to specify the plot via a \code{formula}, 
    e.g., for simple plotting of bivariate surfaces \code{z ~ x + y}, see the examples.}
  \item{y}{If \code{x} is a vector the argument \code{y} and \code{z} must also be supplied as
    vectors.}
  \item{z}{If \code{x} is a vector the argument \code{y} and \code{z} must also be supplied as
    vectors, \code{z} defines the surface given by \eqn{z = f(x, y)}.}
  \item{view}{Which variable should be used for the x-axis of the plot, the other variable will be
    used to compute the slices. May also be a \code{character} with the name of the corresponding
    variable.}
  \item{c.select}{Integer, selects the column that is used in the resulting matrix to be
    used as the \code{z} argument.}
  \item{values}{The values of the \code{x} or \code{y} variable that should be used for computing
    the slices, if set to \code{NULL}, slices will be constructed according to the quantiles, see
    also argument \code{probs}.}
  \item{probs}{Numeric vector of probabilities with values in [0,1] to be used within function
    \code{\link[stats]{quantile}} to compute the \code{values} for plotting the slices.}
  \item{grid}{The grid size of the surface where the slices are generated from.}
  \item{legend}{If set to \code{TRUE}, a legend with the \code{values} that where used for slicing
    will be added.}
  \item{pos}{The position of the legend, see also function \code{\link[graphics]{legend}}.}
  \item{digits}{The decimal place the legend values should be rounded.}
  \item{data}{If \code{x} is a \code{formula}, a \code{data.frame} or \code{list}. By default the 
    variables are taken from \code{environment(x)}: typically the environment from which 
    \code{plot3d} is called.}
  \item{rawdata}{If set to \code{TRUE}, the data will not be interpolated, only raw data will be
    used. This is useful when displaying data on a regular grid.}
  \item{type}{Character, which type of interpolation method should be used. The default is
    \code{type = "akima"}, see function \code{\link[akima]{interp}}. The two other options are
    \code{type = "mba"}, which calls function \code{\link[MBA]{mba.surf}} of package \pkg{MBA}, or
    \code{type = "mgcv"}, which uses a spatial smoother withing package \pkg{mgcv} for interpolation.
    The last option is definitely the slowest, since a full regression model needs to be estimated.}
  \item{linear}{Logical, should linear interpolation be used withing function
     \code{\link[akima]{interp}}?}
  \item{extrap}{Logical, should interpolations be computed outside the observation area
    (i.e., extrapolated)?}
  \item{k}{Integer, the number of basis functions to be used to compute the interpolated surface
    when \code{type = "mgcv"}.}
  \item{rug}{Add a \code{\link[graphics]{rug}} to the plot.}
  \item{jitter}{If set to \code{TRUE} a \code{\link[base]{jitter}}ed 
      \code{\link[graphics]{rug}} plot is added.}
  \item{rug.col}{Specify the color of the rug representation.}
  \item{\dots}{Parameters passed to \code{\link{matplot}} and \code{\link[graphics]{legend}}.}
}

\details{
  Similar to function \code{\link{plot3d}}, this function first applies bivariate interpolation
  on a regular \code{grid}, afterwards the slices are computed from the resulting surface. 
}

\note{
  Function \code{sliceplot} can use the \pkg{akima} package to construct smooth interpolated
  surfaces, therefore, package \pkg{akima} needs to be installed. The \pkg{akima} package has an ACM
  license that restricts applications to non-commercial usage, see

  \url{http://www.acm.org/publications/policies/softwarecrnotice}
  
  Function \code{sliceplot} prints a note referring to the ACM license. This note can be suppressed by
  setting

  \code{options("use.akima" = TRUE)} 
}

\seealso{
  \code{\link{plot2d}}, \code{\link{plot3d}},
  \code{\link{plotmap}}, \code{\link{plotblock}}.
}

\examples{
## Generate some data.
set.seed(111)
n <- 500

## Regressors.
d <- data.frame(z = runif(n, -3, 3), w = runif(n, 0, 6))

## Response.
d$y <- with(d, 1.5 + cos(z) * sin(w) + rnorm(n, sd = 0.6))

\dontrun{## Estimate model.
b <- bamlss(y ~ te(z, w), data = d)
summary(b)

## Plot estimated effect.
plot(b, term = "te(z,w)", sliceplot = TRUE)
plot(b, term = "te(z,w)", sliceplot = TRUE, view = 2)
plot(b, term = "te(z,w)", sliceplot = TRUE, view = "w")
plot(b, term = "te(z,w)", sliceplot = TRUE, probs = seq(0, 1, length = 10))
}

## Variations.
d$f1 <- with(d, sin(z) * cos(w))
sliceplot(cbind(z = d$z, w = d$w, f1 = d$f1))

## Same with formula. 
sliceplot(sin(z) * cos(w) ~ z + w, ylab = "f(z)", data = d)

## Compare with plot3d().
plot3d(sin(z) * 1.5 * w ~ z + w, zlab = "f(z,w)", data = d)
sliceplot(sin(z) * 1.5 * w ~ z + w, ylab = "f(z)", data = d)
sliceplot(sin(z) * 1.5 * w ~ z + w, view = 2, ylab = "f(z)", data = d)
}

\keyword{hplot}

