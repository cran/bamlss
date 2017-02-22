\name{homstart_data}
\alias{homstart_data}

\title{
HOMSTART Precipitation Data
}

\description{
This function downloads and compiles the HOMSTART-project data set. The data is downloaded from
the Zentralanstalt fuer Meteorologie und Geodynamik (ZAMG, \url{http://www.zamg.ac.at}) and
funded by the Austrian Climate Research Programme (ACRP) and is free for research purposes.
}

\usage{
homstart_data(dir = NULL, load = TRUE, tdir = NULL)
}

\arguments{
  \item{dir}{The directory where the \code{homstart.rda} file should be stored.}
  \item{load}{Should the \code{homstart} data be loaded?}
  \item{tdir}{An optional temporary directory where all downloaded files are processed.}
}

\value{
A data frame containing the following variables:
  \item{raw}{The daily precipitation observations.}
  \item{cens}{Precipitation observations censored at 0.}
  \item{bin}{Factor with levels \code{"yes"} or \code{"no"} indicating precipitation.}
  \item{cat}{Factor with levels \code{"none"}, \code{"low"}, \code{"medium"} and \code{"high"}
    indicating the amount of precipitation.}
  \item{trend}{A numeric time trend}
  \item{month}{Month of of the observation.}
  \item{year}{Year of the observation.}
  \item{day}{Day of the year.}
  \item{lon}{The longitude coordinate of the corresponding meteorological station.}
  \item{lat}{The latitude coordinate of the corresponding meteorological station.}
  \item{id}{Factor, meteorological station identifier.}
  \item{cos1, cos2, sin1, sin2}{Transformed time trend for harmonic regression.}
  \item{weekend}{Factor, with levels \code{"yes"} and \code{"no"} indication if the observation
    was measured on a weekend.}
  \item{elevation}{Numeric, the elevation of the meteorological station.}
}

\references{
Nemec J, Gruber C, Chimani B, Auer I (2012).
  Trends in extreme temperature indices in {A}ustria based on a new homogenised dataset.
  \emph{International Journal of Climatology}. DOI 10.1002/joc.3532.

Nemec J, Chimani B, Gruber C, Auer I (2011).
  Ein neuer {D}atensatz homogenisierter {T}agesdaten.
  \emph{OEGM Bulletin}, \bold{1}, 19--20.
  \url{http://www.meteorologie.at/docs/OEGM_bulletin_2011_1.pdf}

Umlauf N, Mayr G, Messner J, Zeileis A (2012).
  Why does it always rain on me? A spatio-temporal analysis of precipitation in Austria.
  \emph{Austrian Journal of Statistics}, \bold{41}(1), 81--92.
  \url{http://www.stat.tugraz.at/AJS/ausg121/121Zeileis.pdf}
}

\examples{
\dontrun{homstart_data(load = TRUE)
head(homstart)
}
}

\keyword{datasets}

