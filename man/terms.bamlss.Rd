\name{terms.bamlss}
\alias{terms.bamlss}
\alias{terms.bamlss.frame}
\alias{terms.bamlss.formula}

\title{
BAMLSS Model Terms
}

\description{
  Extract \code{\link[stats]{terms.object}}s for BAMLSS.
}

\usage{
\method{terms}{bamlss}(x, specials = NULL, data = NULL,
  model = NULL, pterms = TRUE, sterms = TRUE,
  drop = TRUE, ...)

\method{terms}{bamlss.frame}(x, specials = NULL, data = NULL,
  model = NULL, pterms = TRUE, sterms = TRUE,
  drop = TRUE, ...)

\method{terms}{bamlss.formula}(x, specials = NULL, data = NULL,
  model = NULL, pterms = TRUE, sterms = TRUE,
  drop = TRUE, ...)
}

\arguments{
  \item{x}{An \code{link{bamlss}}, \code{\link{bamlss.frame}} or
    \code{\link{bamlss.formula}} object.}
  \item{specials}{See \code{\link[stats]{terms.object}}.}
  \item{data}{Data passed to \code{\link[stats]{terms.formula}}.}
  \item{model}{Character or integer, specifies the model for which terms should be returned.}
  \item{pterms}{Should parametric terms be part of the object?}
  \item{sterms}{Should smooth terms be part of the object?}
  \item{drop}{If terms for only one \code{model} are returned, the list structure is dropped.}
  \item{\dots}{Arguments passed to \code{\link{bamlss.formula}}.}
}

\value{
  Object of class \code{"bamlss.terms"}, a list of \code{\link[stats]{terms.object}}s, depending on
  the structure of the \code{\link{bamlss.formula}} object.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}}, \code{\link{bamlss.formula}}.
}

\examples{
## Model formula.
f <- list(
  num ~ x1 + x2 + id,
  sigma ~ x3 + fac + lon + lat
)

## Create the terms object.
terms(bamlss.formula(f))
}

\keyword{regression}
\keyword{models}
