\name{residuals.bamlss}
\alias{residuals.bamlss}
\alias{plot.bamlss.residuals}

\title{
  Compute BAMLSS Residuals
}

\description{
  Function to compute quantile and response residuals.
}

\usage{
\method{residuals}{bamlss}(object, type = c("quantile", "response"),
  nsamps = NULL, ...)

\method{plot}{bamlss.residuals}(x, which = c("hist-resid", "qq-resid"),
  spar = TRUE, ...)
}

\arguments{
  \item{object}{An object of class \code{"bamlss"}.}
  \item{type}{The type of residuals wanted, possible types are
    \code{"quantile"} residuals and \code{"response"} residuals.}
  \item{nsamps}{If the fitted \code{\link{bamlss}} object contains samples of parameters,
    computing residuals may take quite some time. Therefore, to get a first feeling it can
    be useful to compute residuals only based on \code{nsamps} samples, i.e., \code{nsamps}
    specifies the number of samples which are extracted on equidistant intervals.}
  \item{x}{Object returned from function \code{residuals.bamlss()}.}
  \item{which}{Should a histogram with kernel density estimates be plotted or a qq-plot?}
  \item{spar}{Should graphical parameters be set by the plotting function?}
  \item{\dots}{For function \code{residuals.bamlss()} arguments passed to possible
    \code{$residuals()} functions that may be part of a \code{\link{bamlss.family}}. For function
    \code{plot.bamlss.residuals()} arguments passed to function
    \code{\link[graphics]{hist.default}} and \code{\link[stats]{qqnorm.default}}.}
}

\details{
  Response residuals are the raw residuals, i.e., the response data minus the fitted distributional
  mean. If the \code{\link{bamlss.family}} object contains a function \code{$mu(par, \dots)}, then
  raw residuals are computed with \code{y - mu(par)} where \code{par} is the named list of fitted
  values of distributional parameters. If \code{$mu(par, ...)} is missing, then the fitted values
  of the first distributional parameter are used.

  Randomized quantile residuals are based on the cumulative distribution function of the
  \code{\link{bamlss.family}} object, i.e., the \code{$p(y, par, ...)} function.
}

\value{
  A vector of residuals.
}

\references{
  Dunn P. K., and Smyth G. K. (1996). Randomized Quantile Residuals.
    \emph{Journal of Computational and Graphical Statistics} \bold{5}, 236-244. 
}

\seealso{
  \code{\link{bamlss}}, \code{\link{predict.bamlss}}, \code{\link{fitted.bamlss}}.
}

\examples{
\dontrun{## Generate data.
d <- GAMart()

## Estimate models.
b1 <- bamlss(num ~ s(x1), data = d)
b2 <- bamlss(num ~ s(x1) + s(x2) + s(x3), data = d)

## Extract quantile residuals.
e1 <- residuals(b1, type = "quantile")
e2 <- residuals(b2, type = "quantile")

## Plots.
plot(e1)
plot(e2)
}
}

\keyword{regression}
\keyword{models}

