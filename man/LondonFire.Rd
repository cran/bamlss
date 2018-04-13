\name{LondonFire}
\alias{LondonFire}
\alias{LondonBoroughs}
\alias{LondonBoundaries}
\alias{LondonFStations}

\title{London Fire Data}

\description{
  Provides the compiled dwelling fire data of London in 2015. The data is provided
  as a \code{\link[sp]{SpatialPointsDataFrame}} in object \code{LondonFire}. In addition
  the boundary and borough information is provided in objects \code{LondonBoundaries} and
  \code{LondonBoroughs}. Locations of all fire stations in London of 2015 are provided in
  object \code{LondonFStations}.
}

\usage{data("LondonFire")}

\format{
The \code{LondonFire} data contains the following variables:
\describe{
  \item{arrivaltime:}{Numeric, the time after the emergency call until the first fire engine arrived.}
  \item{daytime:}{Numeric, The time of day at which the emergency call was received.}
  \item{fsintens:}{Numeric, the fire station intensity at the location of the fire scene. The
    intensity is measured using a kernel density estimate of the \code{LondonFStations} using the
    \pkg{spatstat} package.}
}
}

\references{
  London Fire (2015). London Fire Brigade Incident Records.
    \emph{London Data Store}, UK Open Government Licence (OGL v2),
    \url{https://data.london.gov.uk/dataset/london-fire-brigade-incident-records}.

  London Boroughs/Boundaries (2015). Statistical GIS Boundary Files for London.
    \emph{London Data Store}, UK Open Government Licence (OGL v2),
    \url{https://data.london.gov.uk/dataset/statistical-gis-boundary-files-london}.

  London Fire Stations (2015). \emph{London Fire Brigade},
    \url{https://www.london-fire.gov.uk/},
    (old url \code{http://www.london-fire.gov.uk/A-ZFireStations.asp}).

  Taylor BM (2016). Spatial Modelling of Emergency Service Response Times.
    \emph{Journal of the Royal Statistical Society}: Series A (Statistics in Society),
    pp. n/a-n/a. URL \url{http://dx.doi.org/10.1111/rssa.12192}.
}

\examples{
data("LondonFire")
plot(LondonFire, col = "red")
plot(LondonFStations, col = "blue", cex = 1.5, pch = 16, add = TRUE)
plot(LondonBoroughs, lwd = 0.5, add = TRUE)
plot(LondonBoundaries, add = TRUE, lwd = 1.5)
head(LondonFire@data)
}

\keyword{datasets}

