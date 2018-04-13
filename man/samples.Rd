\name{samples}
\alias{samples}
\alias{samples.bamlss}
\alias{samples.bamlss.frame}

\title{
  Extract Samples
}

\description{
  Generic function to extract samples from objects.
}

\usage{
## Generic.
samples(object, ...)

## Method for "bamlss" objects.
\method{samples}{bamlss}(object, model = NULL, term = NULL,
  combine = TRUE, drop = TRUE, burnin = NULL,
  thin = NULL, coef.only = FALSE, ...)
}

\arguments{
  \item{object}{An object for which samples should be extracted.}
  \item{model}{Character or integer, specifies the model for which samples should be extracted.}
  \item{term}{Character or integer, specifies the term for which samples should be extracted.}
  \item{combine}{Samples stored as a \code{\link[coda]{mcmc.list}}, e.g., when a model is estimated
    on multiple cores, can be combined into one large sample matrix.}
  \item{drop}{If there is only one \code{model} for which samples should be extracted, should
    the list structure be dropped?}
  \item{burnin}{Integer, specifies the number of samples that should be withdrawn as a burn-in phase.}
  \item{thin}{Integer, specifies the step length of samples that should be extracted, e.g.,
    \code{thin = 10} mean that only every 10th sample is returned.}
  \item{coef.only}{Logical, should only samples of model coefficients be returned?}
  \item{\dots}{Other arguments.}
}

\seealso{
  \code{\link{bamlss}}.
}

\examples{
\dontrun{## Generate data.
d <- GAMart()

## Estimate model.
b <- bamlss(num ~ s(x1) + s(x2) + s(x3), data = d)

## Extract samples for "s(x2)".
sa <- samples(b, term = "s(x2)")
head(sa)

## Trace plot.
plot(sa)
}
}

\keyword{regression}
\keyword{models}

