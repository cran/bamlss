\name{model.frame.bamlss}
\alias{model.frame.bamlss}
\alias{model.frame.bamlss.frame}
\alias{bamlss.model.frame}

\title{
  BAMLSS Model Frame
}

\description{
  Extracts the model frame of a \code{\link{bamlss}} or \code{\link{bamlss.frame}}
  object.
}

\usage{
\method{model.frame}{bamlss}(formula, ...)

\method{model.frame}{bamlss.frame}(formula, ...)

## Internal helper function for creating
## the model frame.
bamlss.model.frame(formula, data, family = gaussian_bamlss(),
  weights = NULL, subset = NULL, offset = NULL,
  na.action = na.omit, specials = NULL, contrasts.arg = NULL,
  drop.unused.levels = TRUE, ...)
}

\arguments{
  \item{formula}{An object of class \code{"bamlss"} or \code{"bamlss.frame"}.}
  \item{data}{A \code{\link{data.frame}} or \code{\link{list}} containing the model
    response variable(s) and covariates specified in the \code{formula}.
    By default the variables are taken from \code{environment(formula)}.}
  \item{family}{A \code{\link{bamlss.family}} object, see function \code{\link{bamlss}}.}
  \item{weights}{Prior weights on the data.}
  \item{subset}{An optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{offset}{Can be used to supply model offsets for use in fitting.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}'s. The default is set by the \code{na.action}
    setting of \code{\link{options}}, and is \code{\link{na.omit}} if set to \code{NULL}.}
  \item{specials}{Special function in formulas, see \code{\link[stats]{terms.object}}.}
  \item{contrasts.arg}{An optional list. See the \code{contrasts.arg} of 
    \code{\link[stats]{model.matrix.default}}.}
  \item{drop.unused.levels}{Should factors have unused levels dropped?}
  \item{\dots}{Arguments to be passed to \code{bamlss.model.frame()} and others.}
}

\value{
  The \code{\link[base]{data.frame}} containing the variables used for modeling
  with function \code{\link{bamlss}}.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}},
  \code{\link{model.matrix.bamlss.frame}}.
}

\examples{
\dontrun{## Generate some data.
d <- GAMart()

## Model formula.
f <- list(
  log(pnum) ~ s(x1) + log(x2),
  sigma ~ s(x3)
)

## Estimate model.
b <- bamlss(f, data = d)

## Extract the model frame.
head(model.frame(b))
}
}

\keyword{regression}
\keyword{models
}
