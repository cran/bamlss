\name{GEFCom2014}
\alias{GEFCom2014}

\title{GEFCom2014 Electricity Load Data}

\description{
  This data set is a subset of the electricity load data of the GEFCom2014 challenge. The full
  data is published in the Appendix of Hong T. et al. (2016).
}

\usage{
  data("GEFCom2014")
}

\format{
  The \code{GEFCom2014} contains the following variables:
  \describe{
    \item{date}{The date of recording.}
    \item{load}{The corresponding electricty load.}
  }
}

\references{
  Hong T., Pinson P, Fan S., Zareipour H., Troccoli A., and Hyndman R.J. (2016).
  Electric load forecasting, Electricity price forecasting, Wind power forecasting, Solar power forecasting,
  Probabilistic forecasting, Forecasting competition. \emph{International Journal of Forecasting},
  \bold{32}(3), 896--913. \doi{10.1016/j.ijforecast.2016.02.001}
}

\examples{
  data("GEFCom2014", package = "bamlss")
  GEFCom2014$yday <- as.POSIXlt(GEFCom2014$date)$yday
  plot(load ~ yday, data = GEFCom2014)
}

\keyword{datasets}

