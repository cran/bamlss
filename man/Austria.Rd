\name{Austria}
\alias{Austria}
\alias{AustriaTopo}

\title{Austria States and Topography}

\description{
  These two data sets provide the states of Austria as a
  \code{\link[sp]{SpatialPolygonsDataFrame}} and the topography of Austria as
  a \code{\link[raster]{raster}} object.
}

\usage{data("Austria")}

\references{
  Austria (2016). Nuts Einheiten, \emph{STATISTIK AUSTRIA}.
    \url{http://www.statistik.at/web_de/klassifikationen/regionale_gliederungen/nuts_einheiten/index.html}

  Austria Topography (2016). \emph{National Centers for Environmental Information},
    \emph{National Oceanic and Atmospheric Administration (NOAA)}.
    \url{https://www.ngdc.noaa.gov/mgg/global/global.html}
}

\examples{
require("raster")
data("Austria")
plot(AustriaTopo, col = terrain.colors(10))
plot(Austria, add = TRUE)
}

\keyword{datasets}

