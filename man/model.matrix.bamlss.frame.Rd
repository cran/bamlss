\name{model.matrix.bamlss.frame}
\alias{model.matrix.bamlss.frame}
\alias{model.matrix.bamlss.formula}
\alias{model.matrix.bamlss.terms}

\title{
Construct/Extract BAMLSS Design Matrices
}

\description{
  The function creates design (or model) matrices for BAMLSS, i.e., for
  each parameter of a \code{\link{bamlss.family}} object.
}

\usage{
\method{model.matrix}{bamlss.frame}(object, data = NULL, model = NULL,
  drop = TRUE, scale.x = FALSE, ...)

\method{model.matrix}{bamlss.formula}(object, data = NULL, model = NULL,
  drop = TRUE, scale.x = FALSE, ...)

\method{model.matrix}{bamlss.terms}(object, data = NULL, model = NULL,
  drop = TRUE, scale.x = FALSE, ...)
}

\arguments{
  \item{object}{A \code{\link{bamlss.frame}}, \code{\link{bamlss.formula}} or
    \code{\link{terms.bamlss}} object.}
  \item{data}{A data frame or list.}
  \item{model}{Character or integer, specifies the model for which design matrices should be returned.}
  \item{drop}{If model matrices for only one \code{model} are returned, the list structure is dropped.}
  \item{scale.x}{Logical, should the model matrices of the linear parts be scaled?}
  \item{\dots}{Not used.}
}

\value{
  Depending on the type of model a named list of model matrices or a single
  model matrix.
}

\seealso{
  \code{\link[stats]{model.matrix}}, \code{\link{bamlss.frame}},
  \code{\link{bamlss.formula}}, \code{\link{terms.bamlss}}.
}

\examples{
## Generate some data.
d <- GAMart()

## Model formula.
f <- list(
  num ~ x1 + x2 + id,
  sigma ~ x3 + fac + lon + lat
)

## Create a "bamlss.frame".
bf <- bamlss.frame(f, data = d)

## Get the model matrices.
X <- model.matrix(bf)
head(X$sigma)

## Same with "bamlss.formula".
X <- model.matrix(bamlss.formula(f), data = d)
head(X$sigma)
}

\keyword{manip}
\keyword{regression}
\keyword{models}

