\name{plotblock}
\alias{plotblock}

\title{
  Factor Variable and Random Effects Plots
}

\description{
  Function to plot effects for model terms including factor, or group variables for random effects.
}
 
\usage{
plotblock(x, residuals = FALSE, range = c(0.3, 0.3), 
  col.residuals = "black", col.lines = "black", c.select = NULL, 
  fill.select = NULL , col.polygons = NULL, data = NULL,
  shift = NULL, trans = NULL, labels = NULL, ...)
}

\arguments{
  \item{x}{Either a \code{list} of length of the unique factors, where each \code{list} element 
    contains the estimated effects for one factor as a matrix, or 
    one data matrix with first column as the group or factor variable. Also formulas are accepted, 
    e.g it is possible to specify the plot with \code{f ~ x} or \code{f1 + f2 ~ x}. By convention, 
    the covariate for which effects should be plotted, is always in the first column in the 
    resulting data matrix, that is used for plotting, i.e. in the second formula example, the data 
    matrix is \code{cbind(x, f1, f2)}, also see argument \code{c.select} and \code{fill.select}.}
  \item{residuals}{If set to \code{TRUE}, residuals will be plotted if available. Residuals may be
    set as an \code{\link[base]{attr}}ibute of \code{x} named 
    \code{"residuals"}, where the residuals must be a matrix with first column 
    specifying the covariate, and second column the residuals that should be plotted.}
  \item{range}{Numeric vector, specifying the left and right bound of the block.}
  \item{col.residuals}{The color of the partial residuals.}
  \item{col.lines}{Vector of maximum length of columns of \code{x} minus 1, specifying the color of 
    the lines.}
  \item{c.select}{Integer vector of maximum length of columns of \code{x}, selects the 
    columns of the resulting data matrix that should be used for plotting. E.g. if \code{x} has 5 
    columns, then \code{c.select = c(1, 2, 5)} will select column 1, 2 and 5 for plotting. Note that 
    first element of \code{c.select} should always be 1, since this is the column of the covariate 
    the effect is plotted for.}
  \item{fill.select}{Integer vector, select pairwise the columns of the resulting data matrix 
    that should form one polygon with a certain background color specified in argument \code{col}. 
    E.g. \code{x} has three columns, or is specified with formula \code{f1 + f2 ~ x}, then setting 
    \code{fill.select = c(0, 1, 1)} will draw a polygon with \code{f1} and \code{f2} as boundaries. 
    If \code{x} has five columns or the formula is e.g. \code{f1 + f2 + f3 + f4 ~ x}, then setting
    \code{fill.select = c(0, 1, 1, 2, 2)}, the pairs \code{f1}, \code{f2} and \code{f3}, \code{f4}
    are selected to form two polygons.}
  \item{col.polygons}{Specify the background color for the upper and lower confidence bands, e.g. 
      \code{col = c("green", "red")}.}
  \item{data}{If \code{x} is a formula, a \code{data.frame} or \code{list}. By default the variables 
    are taken from \code{environment(x)}: typically the environment from which \code{plotblock} is 
    called.}
  \item{shift}{Numeric constant to be added to the smooth before plotting.}
  \item{trans}{Function to be applied to the smooth before plotting, e.g., to transform the
    plot to the response scale.}
  \item{labels}{Character, labels for the factor levels.}
  \item{\dots}{Graphical parameters, please see the details.}
}

\details{
  Function \code{plotblock} draws for every factor or group the effect as a "block" in one graphic, 
  i.e., similar to boxplots, estimated fitted effects, e.g., containing quantiles of MCMC samples,
  are drawn as one block, where the upper lines represent upper quantiles, the 
  middle line the mean or median, and lower lines lower quantiles, also see the examples. The 
  following graphical parameters may be supplied additionally:
  \itemize{
    \item \code{cex}: Specify the size of partial residuals,
    \item \code{lty}: The line type for each column that is plotted, e.g. \code{lty = c(1, 2)},
    \item \code{lwd}: The line width for each column that is plotted, e.g. \code{lwd = c(1, 2)},
    \item \code{poly.lty}: The line type to be used for the polygons,
    \item \code{poly.lwd}: The line width to be used for the polygons,
    \item \code{density} \code{angle}, \code{border}: See \code{\link[graphics]{polygon}},
    \item \code{\dots}: Other graphical parameters, see function \code{\link[graphics]{plot}}.
    }
}

\seealso{
  \code{\link{plot2d}}, \code{\link{plot3d}},
  \code{\link{plotmap}}, \code{\link{sliceplot}}.
}

\examples{
## Generate some data.
set.seed(111)
n <- 500

## Regressors.
d <- data.frame(fac = factor(rep(1:10, n/10)))

## Response.
d$y <- with(d, 1.5 + rnorm(10, sd = 0.6)[fac] + 
  rnorm(n, sd = 0.6))

\dontrun{## Estimate model.
b <- bamlss(y ~ s(fac,bs="re"), data = d)
summary(b)

## Plot random effects.
plot(b)

## Extract fitted values.
f <- fitted(b, model = "mu", term = "fac")
f <- cbind(d["fac"], f)

## Now use plotblock.
plotblock(f)

## Variations.
plotblock(f, fill.select = c(0, 1, 0, 1), col.poly = "red")
plotblock(f, fill.select = c(0, 1, 0, 1), col.poly = "lightgray",
  lty = c(2, 1, 2), lwd = c(2, 1, 2))
}

## More examples.
plotblock(y ~ fac, data = d, range = c(0.45, 0.45))

d <- data.frame(fac = factor(rep(1:10, n/10)))
d$y <- with(d, c(2.67, 5, 6, 3, 4, 2, 6, 7, 9, 7.5)[fac])
plotblock(y ~ fac, data = d)
plotblock(cbind(y - 0.1, y + 0.1) ~ fac, data = d)
}

\keyword{hplot}

