\name{parameters}
\alias{parameters}

\title{Extract or Initialize Parameters for BAMLSS}

\description{
  The function either sets up a list of all parameters of a \code{\link{bamlss.frame}}, which
  can be used for setting up models, or extracts the estimated parameters of a \code{\link{bamlss}}
  object.
}

\usage{
parameters(x, model = NULL, start = NULL,
  fill = c(0, 1e-04), list = FALSE,
  simple.list = FALSE, extract = FALSE,
  ...)
}

\arguments{
  \item{x}{A \code{\link{bamlss.frame}} or \code{\link{bamlss}} object.}
  \item{model}{The model name for which parameters should be initialized or extracted.}
  \item{start}{A named numeric vector which should be used when creating the parameter list.
    See also function \code{link{set.starting.values}}}
  \item{fill}{Numeric, when setting up a parameter list, the values the should be used for regression
    coefficients (first element of \code{fill}) and for smoothing variances
    (second element of \code{fill}).}
  \item{list}{Should the function return a list of all parameters?}
  \item{simple.list}{Should the names of parameter vectors be dropped?}
  \item{extract}{Should parameters of a \code{\link{bamlss.frame}} be extracted or initialized?}
  \item{\dots}{Currently not used.}
}

\details{
  Parameters for BAMLSS are used for optimizer functions in function \code{\link{bamlss}}.
  The function is useful for initializing all parameters given a \code{\link{bamlss.frame}}
  (which is done internally in function \code{\link{bamlss}}), but also for extracting all
  estimated parameters of some optimizer.

  The naming convention of the parameter list is used by a couple of functions in this package.
  For each parameter of the modeled distribution, e.g., \code{\link{gaussian_bamlss}} has
  parameters \code{"mu"} and \code{"sigma"}, a list element is created. These elements the contain
  the list of all model term parameters. Parametric model terms are indicated with \code{"p"} and
  smooth model terms with \code{"s"}. If the design matrix of a model term in the \code{x} list
  of a \code{\link{bamlss.frame}} does not contain any columns names, then the parameters are named
  with a leading \code{"b"}, otherwise the column names of the design matrix are used. Smoothing
  variances parameter vectors are named with a leading \code{"tau2"}.

  The naming convention is useful when setting up new model fitting engines for \code{\link{bamlss}}
  and is used, e.g., by \code{\link{bfit}} and \code{\link{GMCMC}}, which are based on parameter
  state list objects as provided by function \code{\link{bamlss.engine.setup}}.
}

\value{
A named list of all parameters of a \code{\link{bamlss.frame}} or \code{\link{bamlss}}
object.
}

\seealso{
\code{\link{bamlss.frame}}, \code{\link{bamlss}}, \code{\link{bfit}}, \code{\link{GMCMC}},
\code{\link{get.par}}, \code{\link{set.par}}
}

\examples{
## Create a "bamlss.frame"
set.seed(123)
d <- GAMart()
bf <- bamlss.frame(num ~ s(x1) + te(lon,lat), data = d)

## Create list of all parameters from "bamlss.frame".
p <- parameters(bf, list = TRUE)
str(p)

\dontrun{## Estimate model.
f <- list(num ~ s(x1) + te(lon,lat), sigma ~ s(x1))
b <- bamlss(f, data = d, sampler = FALSE)

## Extract estimated parameters.
parameters(b)
parameters(b, list = TRUE)
}
}

\keyword{regression}

