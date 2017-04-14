\name{plot3d}
\alias{plot3d}

\title{
  Plot 3D Effects
}

\description{
  Function to plot 3D graphics or image and/or contour plots for bivariate effects/functions.
}
 
\usage{
plot3d(x, residuals = FALSE, col.surface = NULL, 
  ncol = 99L, swap = FALSE, col.residuals = NULL, col.contour = NULL, 
  c.select = NULL, grid = 30L, image = FALSE, contour = FALSE, 
  legend = TRUE, cex.legend = 1, breaks = NULL, range = NULL, 
  digits = 2L, d.persp = 1L, r.persp = sqrt(3), outscale = 0,
  data = NULL, sep = "", shift = NULL, trans = NULL,
  type = "mba", linear = FALSE, extrap = FALSE,
  k = 40, ...)
}

\arguments{
  \item{x}{A matrix or data frame, containing the covariates for which the effect should be plotted 
    in the first and second column and at least a third column containing the effect. Another
    possibility is to specify the plot via a \code{formula}, 
    e.g. for simple plotting of bivariate surfaces \code{z ~ x + y}, see the examples.
    \code{x} may also be a character file path to the data to be used for plotting.}
  \item{residuals}{If set to \code{TRUE}, residuals may also be plotted if available. Residuals must
    be supplied as an \code{\link{attr}}ibute named \code{"residuals"}, which
    is a matrix or data frame where the first two columns are covariates and the third column
    the residuals.}
  \item{col.surface}{The color of the surface, may also be a function, e.g. 
    \code{col.surface = heat.colors}.}
  \item{ncol}{the number of different colors that should be generated, if \code{col.surface} is a 
    function.}
  \item{swap}{If set to \code{TRUE} colors will be represented in reverse order.}
  \item{col.residuals}{The color of the partial residuals, or if \code{contour = TRUE} the color of the
    contour lines.}
  \item{col.contour}{The color of the contour lines.}
  \item{c.select}{Integer vector of maximum length of columns of \code{x}, selects the 
    columns of the resulting data matrix that should be used for plotting. E.g. if \code{x} has 5
    columns, then \code{c.select = c(1, 2, 5)} will select column 1, 2 and 5 for plotting. 
    If \code{c.select = 95} or \code{c.select = 80}, function \code{plot3d} will search for the
    corresponding columns to plot a 95\eqn{\%} or 80\eqn{\%} confidence surfaces respectively. Note
    that if e.g. \code{c.select = c(1, 2)}, \code{plot3d} will use columns 1 + 2 and 2 + 2 for 
    plotting.}
  \item{grid}{The grid size of the surface(s).}
  \item{image}{If set to \code{TRUE}, an \code{\link[fields]{image.plot}} is drawn.}
  \item{contour}{If set to \code{TRUE}, a \code{\link[graphics]{contour}} plot is drawn.}
  \item{legend}{If \code{image = TRUE} an additional legend may be added to the plot.}
  \item{cex.legend}{The expansion factor for the legend text, see \code{\link[graphics]{text}}.}
  \item{breaks}{A set of breakpoints for the colors: must give one more breakpoint than 
    \code{ncol}.}
  \item{range}{Specifies a certain range values should be plotted for.}
  \item{digits}{Specifies the legend decimal places.}
  \item{d.persp}{See argument \code{d} in function \code{\link[graphics]{persp}}.}
  \item{r.persp}{See argument \code{r} in function \code{\link[graphics]{persp}}.}
  \item{outscale}{Scales the outer ranges of \code{x} and \code{z} limits used for interpolation.}
  \item{data}{If \code{x} is a \code{formula}, a \code{data.frame} or \code{list}. By default the 
    variables are taken from \code{environment(x)}: typically the environment from which 
    \code{plot3d} is called. Note that \code{data} may also be a character file path to the data.}
  \item{sep}{The field separator character when \code{x} or \code{data} is a character, see function
    \code{\link{read.table}}.}
  \item{shift}{Numeric constant to be added to the smooth before plotting.}
  \item{trans}{Function to be applied to the smooth before plotting, e.g., to transform the
    plot to the response scale.}
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
  \item{\dots}{Parameters passed to \code{\link{colorlegend}} if an image plot with legend is drawn,
    also other graphical parameters, please see the details.}
}

\details{
  For 3D plots the following graphical parameters may be specified additionally: 
  \itemize{
    \item \code{cex}: Specify the size of partial residuals,
    \item \code{col}: It is possible to specify the color for the surfaces if \code{se > 0}, then 
      e.g. \code{col = c("green", "black", "red")},
    \item \code{pch}: The plotting character of the partial residuals,
    \item \code{\dots}: Other graphical parameters passed functions \code{\link[graphics]{persp}},
      \code{\link[fields]{image.plot}} and \code{\link[graphics]{contour}}.
    }
}

\note{
  Function \code{plot3d} can use the \pkg{akima} package to construct smooth interpolated
  surfaces, therefore, package \pkg{akima} needs to be installed. The \pkg{akima} package has an ACM
  license that restricts applications to non-commercial usage, see

  \url{http://www.acm.org/publications/policies/softwarecrnotice}
  
  Function \code{plot3d} prints a note referring to the ACM license. This note can be suppressed by
  setting

  \code{options("use.akima" = TRUE)} 
}

\seealso{
  \code{\link{colorlegend}}, \code{\link{plot2d}}, \code{\link{plotmap}},
  \code{\link{plotblock}}, \code{\link{sliceplot}}.
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
b <- bamlss(y ~ s(z,w), data = d)
summary(b)

## Plot estimated effect.
plot(b, model = "mu", term = "s(z,w)")

## Extract fitted values.
f <- fitted(b, model = "mu", term = "s(z,w)", intercept = FALSE)
f <- cbind(d[, c("z", "w")], f)

## Now use plot3d().
plot3d(f)
plot3d(f, swap = TRUE)
plot3d(f, grid = 100, border = NA)

## Only works if columns are named with
## '2.5%' and '97.5%'.
plot3d(f, c.select = 95, border = c("red", NA, "green"),
  col.surface = c(1, NA, 1), resid = TRUE, cex.resid = 0.2)

## Now some image and contour.
# plot3d(f, image = TRUE, legend = FALSE)
# plot3d(f, image = TRUE, legend = TRUE)
# plot3d(f, image = TRUE, contour = TRUE)
# plot3d(f, image = TRUE, contour = TRUE, swap = TRUE)
# plot3d(f, image = TRUE, contour = TRUE, col.contour = "white")
# plot3d(f, contour = TRUE)
# plot3d(f, image = TRUE, contour = TRUE, c.select = 3)
# plot3d(f, image = TRUE, contour = TRUE, c.select = "Mean")
# plot3d(f, image = TRUE, contour = TRUE, c.select = "97.5%")
}

## Variations.
d$f1 <- with(d, sin(z) * cos(w))
with(d, plot3d(cbind(z, w, f1)))

## Same with formula. 
plot3d(sin(z) * cos(w) ~ z + w, zlab = "f(z,w)", data = d)
plot3d(sin(z) * cos(w) ~ z + w, zlab = "f(z,w)", data = d, 
  ticktype = "detailed")

## Play with palettes.
plot3d(sin(z) * cos(w) ~ z + w, col.surface = heat.colors, data = d)
plot3d(sin(z) * cos(w) ~ z + w, col.surface = topo.colors, data = d)
plot3d(sin(z) * cos(w) ~ z + w, col.surface = cm.colors, data = d)
plot3d(sin(z) * cos(w) ~ z + w, col.surface = rainbow, data = d)
plot3d(sin(z) * cos(w) ~ z + w, col.surface = terrain.colors, data = d)

plot3d(sin(z) * cos(w) ~ z + w, col.surface = rainbow_hcl, data = d)
plot3d(sin(z) * cos(w) ~ z + w, col.surface = diverge_hcl, data = d)
plot3d(sin(z) * cos(w) ~ z + w, col.surface = sequential_hcl, data = d)

plot3d(sin(z) * cos(w) ~ z + w, 
  col.surface = rainbow_hcl(n = 99, c = 300, l = 80, start = 0, end = 100), 
  data = d)
# plot3d(sin(z) * cos(w) ~ z + w, 
#  col.surface = rainbow_hcl(n = 99, c = 300, l = 80, start = 0, end = 100), 
#  image = TRUE, grid = 200, data = d)
}

\keyword{hplot}

