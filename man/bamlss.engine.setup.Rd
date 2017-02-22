\name{bamlss.engine.setup}
\alias{bamlss.engine.setup}

\title{
  BAMLSS Engine Setup Function
}

\description{
  This function takes the \code{x} object of a \code{\link{bamlss.frame}} and adds
  additional objects that are useful for model fitting engines. This is applied only
  for 'regular' terms, e.g., as created by \code{\link[mgcv]{s}} and \code{\link[mgcv]{te}}.
  For special model terms the corresponding \code{\link[mgcv]{smooth.construct}} method is in charge
  of this (see also the examples for function \code{\link{bfit}}).
}

\usage{
bamlss.engine.setup(x, update = "iwls", propose = "iwlsC_gp",
  do.optim = NULL, df = NULL, parametric2smooth = TRUE, ...)
}

\arguments{
  \item{x}{The \code{x} list, as returned from function \code{\link{bamlss.frame}},
    holding all model matrices and other information that is used for
    fitting the model.}
  \item{update}{Sets the updating function for model terms, see function \code{\link{bfit}}.}
  \item{propose}{Sets the propose function for model terms, see function \code{\link{GMCMC}}.}
  \item{do.optim}{Adds list element \code{"do.optim"} in the \code{"state"} element, see the
    details below.}
  \item{df}{The initial degrees of freedom that should be assigned to a smooth model term,
    based on the trace of the smoother matrix. Note that \code{df} can be a named numeric vector.
    If the names match the labels of the model terms, the corresponding \code{df} are used, e.g.,
    \code{df = c("s(x1)" = 1, "s(x2)" = 2)} sets different \code{df}s for each term.}
  \item{parametric2smooth}{Should parametric model terms be transformed into an artificial
    smooth model term and be added to the \code{"smooth.construct"} object within the \code{x} list?
    This feature is handy, since algorithms can then cycle over the \code{"smooth.construct"}
    object, only.}
  \item{\dots}{Currently not used.}
}

\details{
  For each model term in the \code{"smooth.construct"} object of the \code{x} list (as returned
  from \code{\link{bamlss.frame}}), this function adds a named list called \code{"state"} with the
  following entries:
  \itemize{
    \item \code{"parameters"}: A numeric vector. Regression coefficients are named with \code{"b"},
      smooth variances are named with \code{"tau2"}.
    \item \code{"fitted.values"}: Given the \code{"parameters"}, the actual fitted values of
      the model term.
    \item \code{"edf"}: Given the smoothing variances, the actual equivalent degrees of freedom
      (edf) of the model term.
    \item \code{"do.optim"}: Should an optimizer function try to find optimum smoothing variances?
  }
  The state will be changed in each iteration and can be passed outside an updating function.

  Additionally, if missing in the \code{xt} argument of a model term (see, e.g., function
  \code{\link[mgcv]{s}} for \code{xt}) the function adds the corresponding log-prior and its first
  and second order derivatives w.r.t. regression coefficients in functions \code{grad()} and
  \code{hess()}.

  Also, objects named \code{"lower"} and \code{"upper"} are added to each model term. These
  indicate the lower and upper boundaries of the parameter space.
}

\value{
  A transformed \code{x} list, as returned from function \code{\link{bamlss.frame}}.
}

\seealso{
  \code{\link{bamlss.frame}}, \code{\link{bfit}}, \code{\link{GMCMC}},
  \code{\link{get.par}}, \code{\link{set.par}}, \code{\link{get.state}}
}

\examples{
d <- GAMart()
bf <- bamlss.frame(num ~ s(x1) + s(x2), data = d, family = "gaussian")
names(bf$x$mu$smooth.construct)
bf$x <- bamlss.engine.setup(bf$x, df = c("s(x1)" = 1, "s(x2)" = 3))
names(bf$x$mu$smooth.construct)
names(bf$x$mu$smooth.construct[["s(x1)"]])
names(bf$x$mu$smooth.construct[["s(x1)"]]$state)
sapply(bf$x$mu$smooth.construct, function(x) {
  c(x$state$edf, get.state(x, "tau2"))
})
}

\keyword{regression}
\keyword{models}

