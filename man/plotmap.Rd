\name{plotmap}
\alias{plotmap}

\title{
  Plot Maps
}

\description{
  The function takes a \code{\link[base]{list}} polygons and draws the corresponding map. Different 
  colors for each polygon can be used.
}
 
\usage{
plotmap(map, x = NA, id = NULL, select = NULL,
  legend = TRUE, names = FALSE, values = FALSE, ...)
}

\arguments{
  \item{map}{The map to be plotted, usually an object that inherits from
    \code{\link[sp]{SpatialPolygons}}, but may also be a list of polygons, i.e., each list
    entry is a matrix with x- and y-coordinates.}
  \item{x}{A vector, data.frame or matrix. In the latter case \code{x} should have two columns,
    one column that is the region identifier and one that contains the values to be plotted.
    In case \code{x} is a data.frame, the function searches for \code{"character"} or
    \code{"factor"} columns for the region identifier. If \code{x} is a matrix, the second
    column of \code{x} is supposed to be the region identifier. If \code{x = NA} and
    \code{map} is an object of class \code{"SpatialPolygonsDataFrame"} only the polygons
    will be plotted without using the data.}
  \item{id}{If argument \code{x} is a vector, argument \code{id} should contain the region identifier
    vector.}
  \item{select}{Select the column of the data in \code{x} which should be used for plotting, may
    be an integer or character with the corresponding column name.}
  \item{legend}{Should a color legend be added to the plot, see also function \code{\link{colorlegend}}.}
  \item{names}{If set to \code{TRUE} the name for each polygon will also be plotted at the centroids
    of the corresponding polygons.}
  \item{values}{If set to \code{TRUE} the corresponding values for each polygon will also be plotted 
    at the centroids of the polygons.}
  \item{\dots}{Arguments to be passed to \code{\link{colorlegend}} and others, e.g. change the 
    \code{border} of the polygons and \code{density}, see \code{\link[graphics]{polygon}}.}
}

\seealso{
  \code{\link{plot2d}}, \code{\link{plot3d}},
  \code{\link{sliceplot}}, \code{\link{plotblock}}.
}

\examples{
## Example from mgcv ?mrf.
## Load Columbus Ohio crime data (see ?columbus for details and credits).
data("columb", package = "mgcv")
data("columb.polys", package = "mgcv")

## Plot the map.
plotmap(columb.polys)

## Plot aggregated data.
a <- with(columb, aggregate(crime,
  by = list("district" = district), FUN = mean))

plotmap(columb.polys, x = a$x, id = a$district)
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft")
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft", side.legend = 2)
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft", side.legend = 2, side.ticks = 2)
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft", side.legend = 2, side.ticks = 2,
  col = heat_hcl, swap = TRUE)
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft", side.legend = 2, side.ticks = 2,
  col = heat_hcl, swap = TRUE, range = c(10, 50))
plotmap(columb.polys, x = a$x, id = a$district,
  pos = "topleft", side.legend = 2, side.ticks = 2,
  col = heat_hcl(5), swap = TRUE, range = c(10, 50),
  lrange = c(0, 60))

\dontrun{## Estimate model.
xt <- list("polys" = columb.polys)
b <- bamlss(crime ~ s(district,bs="mrf",xt=xt) + s(income), data = columb)

## Plot effects.
plot(b, map = columb.polys, pos = "top")
}}

\keyword{hplot}

