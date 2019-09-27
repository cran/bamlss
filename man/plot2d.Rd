\name{plot2d}
\alias{plot2d}
\alias{plotnonp}
  
\title{
  Plot 2D Effects
}
  
\description{
  Function to plot simple 2D graphics for univariate effects/functions.
}
 
\usage{
plot2d(x, residuals = FALSE, rug = FALSE, jitter = TRUE, 
  col.residuals = NULL, col.lines = NULL, col.polygons = NULL, 
  col.rug = NULL, c.select = NULL, fill.select = NULL, 
  data = NULL, sep = "", month = NULL, year = NULL,
  step = 12, shift = NULL, trans = NULL,
  scheme = 2, s2.col = NULL, grid = 50, ...)
}

\arguments{
  \item{x}{A matrix or data frame, containing the covariate for which the effect should be plotted 
    in the first column and at least a second column containing the effect. Another possibility is
    to specify the plot via a \code{formula}, e.g. \code{y ~ x}, see the examples.
    \code{x} may also be a character file path to the data to be used for plotting.}
  \item{residuals}{If set to \code{TRUE}, residuals may also be plotted if available.
    Residuals must be supplied as an \code{\link{attr}}ibute named \code{"residuals"}, which
    is a matrix or data frame where the first column is the covariate and the second column
    the residuals.}
  \item{rug}{Add a \code{\link[graphics]{rug}} to the plot.}
  \item{jitter}{If set to \code{TRUE} a \code{\link[base]{jitter}}ed 
      \code{\link[graphics]{rug}} plot is added.}
  \item{col.residuals}{The color of the partial residuals.}
  \item{col.lines}{The color of the lines.}
  \item{col.polygons}{Specify the background color of polygons, if \code{x} has at least 3 columns, 
    i.e. column 2 and 3 can form one polygon.}
  \item{col.rug}{Specify the color of the rug representation.}
  \item{c.select}{Integer vector of maximum length of columns of \code{x}, selects the 
    columns of the resulting data matrix that should be used for plotting. E.g. if \code{x} has 5
    columns, then \code{c.select = c(1, 2, 5)} will select column 1, 2 and 5 for plotting. Note that
    first element of \code{c.select} should always be the column that holds the variable for the
    x-axis.}
  \item{fill.select}{Integer vector, select pairwise the columns of the resulting data matrix 
    that should form one polygon with a certain background color specified in argument \code{col}. 
    E.g. \code{x} has three columns, or is specified with formula \code{f1 + f2 ~ x}, then setting 
    \code{fill.select = c(0, 1, 1)} will draw a polygon with \code{f1} and \code{f2} as boundaries. 
    If \code{x} has five columns or the formula is e.g. \code{f1 + f2 + f3 + f4 ~ x}, then setting
    \code{fill.select = c(0, 1, 1, 2, 2)}, the pairs \code{f1}, \code{f2} and \code{f3}, \code{f4}
    are selected to form two polygons.}
  \item{data}{If \code{x} is a formula, a \code{data.frame} or \code{list}. By default the variables 
    are taken from \code{environment(x)}: typically the environment from which \code{plot2d} is 
    called. Note that \code{data} may also be a character file path to the data.}
  \item{sep}{The field separator character when \code{x} or \code{data} is a character, see function
    \code{\link{read.table}}.}
  \item{month, year, step}{Provide specific annotation for plotting estimation results for temporal
    variables. \code{month} and \code{year} define the minimum time point whereas step specifies the
    type of temporal data with \code{step = 4}, \code{step = 2} and \code{step = 1} corresponding to
    quarterly, half yearly and yearly data.}
  \item{shift}{Numeric constant to be added to the smooth before plotting.}
  \item{trans}{Function to be applied to the smooth before plotting, e.g., to transform the
    plot to the response scale.}
  \item{scheme}{Sets the plotting scheme for polygons, possible values are \code{1} and \code{2}.}
  \item{s2.col}{The color for the second plotting scheme.}
  \item{grid}{Integer, specifies the number of polygons for the second plotting scheme.}
  \item{\dots}{Other graphical parameters, please see the details.}
}

\details{
  For 2D plots the following graphical parameters may be specified additionally: 
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
  \code{\link{plot3d}}, \code{\link{plotmap}},
  \code{\link{plotblock}}, \code{\link{sliceplot}}.
}

\examples{
## Generate some data.
set.seed(111)
n <- 500
## Regressor.
d <- data.frame(x = runif(n, -3, 3))

## Response. 
d$y <- with(d, 10 + sin(x) + rnorm(n, sd = 0.6))

\dontrun{## Estimate model.
b <- bamlss(y ~ s(x), data = d)
summary(b)

## Plot estimated effect.
plot(b)
plot(b, rug = FALSE)

## Extract fitted values.
f <- fitted(b, model = "mu", term = "s(x)")
f <- cbind(d["x"], f)

## Now use plot2d.
plot2d(f)
plot2d(f, fill.select = c(0, 1, 0, 1))
plot2d(f, fill.select = c(0, 1, 0, 1), lty = c(2, 1, 2))
plot2d(f, fill.select = c(0, 1, 0, 1), lty = c(2, 1, 2),
  scheme = 2)

## Variations.
plot2d(sin(x) ~ x, data = d)
d$f <- with(d, sin(d$x))
plot2d(f ~ x, data = d)
d$f1 <- with(d, f + 0.1)
d$f2 <- with(d, f - 0.1)
plot2d(f1 + f2 ~ x, data = d)
plot2d(f1 + f2 ~ x, data = d, fill.select = c(0, 1, 1), lty = 0)
plot2d(f1 + f2 ~ x, data = d, fill.select = c(0, 1, 1), lty = 0,
  density = 20, poly.lty = 1, poly.lwd = 2)
plot2d(f1 + f + f2 ~ x, data = d, fill.select = c(0, 1, 0, 1), 
  lty = c(0, 1, 0), density = 20, poly.lty = 1, poly.lwd = 2)
}}

\keyword{hplot}

