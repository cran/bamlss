\name{stabsel}
\alias{stabsel}
\alias{plot.stabsel}
\title{Stability selection.}
\description{
Performs stability selection based on gradient boosting.
}
\usage{
stabsel(formula, data, family = "gaussian",
  q, maxit, B = 100, thr = .9, fraction = 0.5, seed = NULL, ...)

## Plot selection frequencies.
\method{plot}{stabsel}(x, show = NULL,
  pal = function(n) gray.colors(n, start = 0.9, end = 0.3), ...)
}
\arguments{
\item{formula}{ A formula or extended formula. }
\item{data}{ A \code{\link{data.frame}}. }
\item{family}{ A \code{\link{bamlss.family}} object. }
\item{q}{ An integer specifying how many terms to select in each boosting run. }
\item{maxit}{ An integer specifying the maximum number of boosting iterations.
  See \code{\link{boost}}. Either choose \code{q} or \code{maxit} as hyper-parameter
  for regularization.}
\item{B}{ An integer. The boosting is run B times. }
\item{thr}{ Cut-off threshold of relative frequencies (between 0 and 1) for selection. }
\item{fraction}{ Numeric between 0 and 1. The fraction of data to be used in each
  boosting run. }
\item{seed}{A seed to be set before the stability selection.}
\item{x}{ A object of class stabsel. }
\item{show}{ Number of terms to be shown. }
\item{pal}{ Color palette for different model terms. }
\item{\dots}{ Not used yet in \code{stabsel}. }
}
\value{
A object of class stabsel.
}
\details{
\code{stabsel} performs stability selection based on gradient
boosting (\code{\link{boost}}): The boosting algorithm is run
\code{B} times on a randomly drawn \code{fraction} of the \code{data}.
Each boosting run is stopped either when \code{q} terms have been selected,
or when \code{maxit} iterations have been performed, i.e. either \code{q}
or \code{maxit} can be used to tune the regularization of the boosting.
After the boosting the relative selection frequencies are evaluated.
Terms with a relative selection frequency larger then \code{thr}
are suggested for a final regression model.

If neither \code{q} nor \code{maxit} has been specified, \code{q}
will be set to the square root of the number of columns in \code{data}.

Gradient boosting does not depend on random numbers. Thus, the
individual boosting runs differ only in the subset of data which
is used.
}
\author{Thorsten Simon}
\examples{
\dontrun{## Simulate some data.
set.seed(111)
d <- GAMart()
n <- nrow(d)

## Add some noise variables.
for(i in 4:9)
  d[[paste0("x",i)]] <- rnorm(n)

f <- paste0("~ ", paste("s(x", 1:9, ")", collapse = "+", sep = ""))
f <- paste(f, "+ te(lon,lat)")
f <- as.formula(f)
f <- list(update(f, num ~ .), f)

## Run stability selection.
sel <- stabsel(f, data = d, q = 6, B = 10)
plot(sel)

## Estimate selected model.
nf <- formula(sel)
b <- bamlss(nf, data = d)
plot(b)
}
}
\keyword{model selection}

