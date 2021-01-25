\name{fatalities}
\alias{fatalities}

\title{Weekly Number of Fatalities in Austria}

\description{
  This data set includes weekly fatalities in Austria from 2000 to 46 weeks in 2020. The data
  is taken from the Eurostat data base.
}

\usage{data("fatalities")}

\format{
The \code{fatalities} data contains the following variables:
\describe{
  \item{num:}{Integer, the number of fatalities.}
  \item{year:}{Integer, the corresponding year fatalities are recorded.}
  \item{week:}{Integer, the corresponding week fatalities are recorded..}
}
}

\references{
  Eurostat Database (2020). \emph{Population and social conditions, demography and migration,
    mortality, weekly deaths, deaths by week and NUTS 3 region, Austria}
    \url{https://ec.europa.eu/eurostat/}
}

\examples{
data("fatalities")
plot(num ~ week, data = fatalities)
}

\keyword{datasets}

