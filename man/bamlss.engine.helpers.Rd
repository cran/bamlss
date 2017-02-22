\name{bamlss.engine.helpers}
\alias{bamlss.engine.helpers}
\alias{get.par}
\alias{set.par}
\alias{get.state}
\alias{set.starting.values}

\title{
  BAMLSS Engine Helper Functions
}

\description{
  These functions can be useful when setting up new model fitting engines that are based
  on the setup function \code{\link{bamlss.engine.setup}}. See the examples.
}

\usage{
## Functions to extract parameter states.
get.par(x, what = NULL)
get.state(x, what = NULL)
set.par(x, replacement, what)

## Function for setting starting values.
set.starting.values(x, start)
}

\arguments{
  \item{x}{For function \code{get.par()} and \code{set.par()} argument \code{x} is a
    named numeric vector. For function \code{get.state()} argument \code{x} is an object
    of the \code{\link{smooth.construct}} list that is processed by function
    \code{\link{bamlss.engine.setup}}, i.e., which has a \code{"state"} object.
    For function \code{set.starting.values()} argument \code{x} is the \code{x} list, as
    returned from function \code{\link{bamlss.frame}}.}
  \item{what}{The name of the parameter(s) that should be extracted or replaced.}
  \item{replacement}{The value(s) that should be used for replacement.}
  \item{start}{The named numeric vector of starting values. The name convention is based
    on function \code{\link{parameters}}.}
}

\seealso{
  \code{\link{bamlss.engine.setup}}
}

\examples{
## Create a bamlss.frame.
d <- GAMart()
bf <- bamlss.frame(num ~ s(x1) + s(x2) + te(lon,lat), data = d, family = "gaussian")
names(bf$x$mu$smooth.construct)

## Use the setup function for
## adding state elements.
bf$x <- bamlss.engine.setup(bf$x, df = c("s(x1)" = 1, "s(x2)" = 3))
names(bf$x$mu$smooth.construct)

## Extract regression coefficients.
get.state(bf$x$mu$smooth.construct[["te(lon,lat)"]], "b")

## Extract smoothing variances.
get.state(bf$x$mu$smooth.construct[["te(lon,lat)"]], "tau2")

## More examples.
state <- bf$x$mu$smooth.construct[["te(lon,lat)"]]$state
get.par(state$parameters, "b")
get.par(state$parameters, "tau2")

state$parameters <- set.par(state$parameters, c(0.1, 0.5), "tau2")
get.par(state$parameters, "tau2")

## Setting starting values.
start <- c("mu.s.s(x1).b" = 1:9, "mu.s.s(x1).tau2" = 0.1)
bf$x <- set.starting.values(bf$x, start = start)
get.state(bf$x$mu$smooth.construct[["s(x1)"]], "b")
get.state(bf$x$mu$smooth.construct[["s(x1)"]], "tau2")
}

\keyword{regression}

