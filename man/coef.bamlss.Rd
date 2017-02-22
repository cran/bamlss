\name{coef.bamlss}
\alias{coef.bamlss}
\alias{confint.bamlss}

\title{
  Extract BAMLSS Coefficients
}

\description{
  Methods to extract coefficients of fitted \code{\link{bamlss}} objects, either coefficients
  returned from optimizer functions, or samples from a sampler functions.

  Method \code{confint.bamlss()} produces credible intervals or parameter samples
  using quantiles.
}

\usage{
\method{coef}{bamlss}(object, model = NULL, term = NULL,
  FUN = NULL, parameters = NULL,
  pterms = TRUE, sterms = TRUE,
  hyper.parameters = TRUE, list = FALSE,
  full.names = TRUE, rescale = FALSE, ...)

\method{confint}{bamlss}(object, parm, level = 0.95,
  model = NULL, pterms = TRUE, sterms = FALSE,
  full.names = FALSE, hyper.parameters = FALSE, ...)
}

\arguments{
  \item{object}{An object of class \code{"bamlss"}}
  \item{model}{Character or integer. For which model should coefficients be extracted?}
  \item{term}{Character or integer. For which term should coefficients be extracted?}
  \item{FUN}{A function that is applied on the parameter samples.}
  \item{parameters}{If is set to \code{TRUE}, additionally adds estimated parameters
    returned from an optimizer function (if available).}
  \item{pterms}{Should coefficients of parametric terms be included?}
  \item{sterms}{Should coefficients of smooths terms be included?}
  \item{hyper.parameters}{For smooth terms, should hyper parameters such as smoothing variances
    be included?}
  \item{list}{Should the returned object have a list structure for each distribution parameter?}
  \item{full.names}{Should full names be assigned, indicating whether a term is parametric "p" or
    smooth "s".}
  \item{rescale}{Should parameters of the linear parts be rescaled if the \code{scale.d} argument
    in \code{\link{bamlss.frame}} is set to \code{TRUE}.}
  \item{parm}{Character or integer. For which term should coefficients intervals be extracted?}
  \item{level}{The credible level which defines the lower and upper quantiles that should
    be computed from the samples.}
  \item{\dots}{Arguments to be passed to \code{FUN} and function \code{\link{samples}}.}
}

\value{
  Depending on argument \code{list} and the number of distributional parameters, either a
  \code{list} or vector/matrix of model coefficients.
}

\seealso{
  \code{\link{bamlss}}.
}

\examples{
\dontrun{## Simulate data.
d <- GAMart()

## Model formula.
f <- list(
  num ~ s(x1) + s(x2) + s(x3),
  sigma ~ s(x1) + s(x2) + s(x3)
)

## Estimate model.
b <- bamlss(f, data = d)

## Extract coefficients based on MCMC samples.
coef(b)

## Now only the mean.
coef(b, FUN = mean)

## As list without the full names.
coef(b, FUN = mean, list = TRUE, full.names = FALSE)

## Coefficients only for "mu".
coef(b, model = "mu")

## And "s(x2)".
coef(b, model = "mu", term = "s(x2)")

## With optimizer parameters.
coef(b, model = "mu", term = "s(x2)", parameters = TRUE)

## Only parameteric part.
coef(b, sterms = FALSE, hyper.parameters = FALSE)

## For sigma.
coef(b, model = "sigma", sterms = FALSE,
  hyper.parameters = FALSE)

## 95 perc. credible interval based on samples.
confint(b)
}
}

\keyword{regression}
\keyword{models}

