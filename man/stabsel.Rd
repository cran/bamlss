\name{stabsel}
\alias{stabsel}
\alias{plot.stabsel}
\title{Stability selection.}
\description{
Performs stability selection based on gradient boosting.
}
\usage{
stabsel(formula, data, family = "gaussian",
  q, B = 100, thr = .9, ...)

## Plot selection frequencies.
\method{plot}{stabsel}(x, show = NULL,
  pal = function(n) gray.colors(n, start = 0.9, end = 0.3), ...)
}
\arguments{
\item{formula}{ A formula or extended formula. }
\item{data}{ A \code{\link{data.frame}}. }
\item{family}{ A \code{\link{bamlss.family}} object. }
\item{q}{ An integer specifying how many terms to select in each boosting run. }
\item{B}{ An integer. The boosting is run B times. }
\item{thr}{ Cut-off threshold for selection. }
\item{x}{ A object of class stabsel. }
\item{show}{ Number of terms to be shown. }
\item{pal}{ Color palette for different model terms. }
\item{\dots}{ Not used yet in \code{stabsel}. }
}
\value{
A object of class stabsel.
}
\author{Thorsten Simon}
\examples{
\dontrun{## Simulate some data.
d <- GAMart()
n <- nrow(d)

## Add some noise variables.
for(i in 4:9)
  d[[paste0("x",i)]] <- rnorm(n)

f <- paste0("~ ", paste("s(x", 1:9, ")", collapse = "+", sep = ""))
f <- paste(f, "+ te(lon,lat)")
f <- eval(parse(text = f))
f <- list(update(f, num ~ .), f)

## Run stability selection.
sel <- stabsel(f, data = d, q = 6, B = 10)
plot(sel)

## Estimate selected model.
b <- bamlss(sel$formula.new, data = d)
plot(b)
}
}
\keyword{model selection}

