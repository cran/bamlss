\name{Austria}
\alias{Austria}
\alias{AustriaTopo}

\title{Austria States and Topography}

\description{
  These two data sets \code{Austria} and \code{AustriaTopo} provide the states of Austria as a
  \code{\link[sp]{SpatialPolygonsDataFrame}} and the topography of Austria as
  data frame.
}

\usage{data("Austria")}

\format{
The \code{AustriaTopo} data contains the following variables:
\describe{
  \item{lon:}{Numeric, longitude coordinates.}
  \item{lat:}{Numeric, latitude coordinates.}
  \item{alt:}{Numeric, altitude.}
}
}

\references{
  Austria (2016). Nuts Einheiten, \emph{STATISTIK AUSTRIA}.
    \code{http://www.statistik.at/web_de/klassifikationen/regionale_gliederungen/nuts_einheiten/index.html}

  Austria Topography (2016). \emph{National Centers for Environmental Information},
    \emph{National Oceanic and Atmospheric Administration (NOAA)}.
    \code{https://www.ngdc.noaa.gov/mgg/global/global.html}
}

\examples{
data("Austria")
plot(Austria)
print(head(AustriaTopo))
}

\keyword{datasets}

