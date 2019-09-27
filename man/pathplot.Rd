\name{pathplot}
\alias{pathplot}

\title{Plot Coefficients Paths}

\description{
  This is a simple wrapper function to plot coefficients paths obtained from
  the boosting optimizer function \code{\link{boost}} and the LASSO optimizer
  \code{\link{lasso}}.
}

\usage{
pathplot(object, ...)
}

\arguments{
  \item{object}{An object of class \code{"bamlss"}.}
  \item{\dots}{Arguments passed to \code{\link{boost_plot}} or \code{\link{lasso_plot}}.}
}

\seealso{
  \code{\link{boost_plot}}, \code{\link{lasso_plot}}
}

\keyword{hplot}

