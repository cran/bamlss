\name{bamlss.formula}
\alias{bamlss.formula}

\title{
  Formulae for BAMLSS
}

\description{
  This function creates an extended BAMLSS \code{\link{formula}}. In combination with a
  \code{\link{bamlss.family}} object, each parameter of the response distribution is linked
  to a single formula. If no formula is supplied for a parameter, a simple intercept only
  model is created. Moreover, the function identifies hierarchical structures, see the examples.
  This function is useful for creating complex \code{\link{model.frame}}s for (hierarchical)
  multi-parameter models and is used by function \code{\link{bamlss.frame}}. 
}

\usage{
bamlss.formula(formula, family = NULL, specials = NULL, ...)
}

\arguments{
  \item{formula}{A simple \code{\link{formula}}, or a \code{\link{list}} of simple formulae,
    or an extended \code{\link[Formula]{Formula}}. For formula lists or extended
    \code{\link[Formula]{Formula}}e, each single formula represents one model for the
    respective parameter as specified in the \code{family} object.}
  \item{family}{A \code{\link{bamlss.family}} object.}
  \item{specials}{A character vector specifying special functions to be used within
    formulae, see \code{\link{terms.object}}.}
  \item{\dots}{Arguments passed to the \code{family}.}
}

\value{
  A named list of class \code{"bamlss.formula"}. Each list entry specifies a model, e.g., for
  one parameter of a provided \code{\link{bamlss.family}} object. Each entry (parameter model) then
  holds:
  \item{formula}{A simple \code{\link{formula}}.}
  \item{fake.formula}{A formula with all function calls being dropped, e.g., the formula
    \code{y ~ s(x1) + s(x2)} is represented in the \code{fake.formula} entry as
    \code{y ~ x1 + x2}. The \code{fake.formula} is useful for creating model frames.}
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}}, \code{\link{bamlss.family}}
}

\examples{
## Simple formula without family object.
f <- bamlss.formula(y ~ x1 + s(x2))
print(f)
print(str(f))

## Complex formula with list of formulae.
f <- list(
  y1 ~ x1 + s(x2),
  y2 ~ x3 + te(lon, lat),
  u ~ x4 + x1
)

f <- bamlss.formula(f)
print(f)
print(names(f))

## Same formula but using extended formulae
## of package Formula.
f <- y1|y2|u ~ x1 + s(x2)|x3 + te(lon,lat)|x4 + x1
f <- bamlss.formula(f)
print(f)
print(names(f))

## Using a bamlss family object, e.g., gaussian_bamlss().
## The family has two parameters, mu and sigma, for
## each parameter one formula is returned. If no
## formula is specified an intercept only model is
## generated for the respective parameter.
f <- bamlss.formula(y ~ x1 + s(x2), family = gaussian_bamlss)

## Note, same as:
f <- bamlss.formula(y ~ x1 + s(x2), family = "gaussian")
print(f)

## Specify model for parameter sigma, too
f <- list(
  y ~ x1 + s(x2),
  sigma ~ x2 + te(lon, lat)
)
f <- bamlss.formula(f, family = "gaussian")
print(f)

## With complex hierarchical structures,
## each parameter is another list of formulae,
## indicated by the h1,...,hk, names.
f <- list(
  y ~ x1 + s(x2) + id1,
  sigma ~ x2 + te(lon, lat) + id2,
  id1 ~ s(x3) + x4 + s(id3),
  id3 ~ x5 + s(x5, x6),
  id2 ~ x7
)
f <- bamlss.formula(f, family = "gaussian")
print(f)
}

\keyword{regression}

