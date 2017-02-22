\name{bamlss.frame}
\alias{bamlss.frame}

\title{
  Create a Model Frame for BAMLSS
}

\description{
  This function parses the \code{data} and the model \code{formula}, or
  extended \code{\link{bamlss.formula}}, as well as the \code{\link{bamlss.family}}
  into a \code{\link{bamlss.frame}} object. The \code{\link{bamlss.frame}} then holds all model
  matrices and information that is needed for setting up estimation engines.
}

\usage{
bamlss.frame(formula, data = NULL, family = "gaussian",
  weights = NULL, subset = NULL, offset = NULL,
  na.action = na.omit, contrasts = NULL,
  knots = NULL, specials = NULL, reference = NULL,
  model.matrix = TRUE, smooth.construct = TRUE,
  ytype = c("matrix", "vector", "integer"),
  scale.x = FALSE, scale.d = FALSE, ...)
}

\arguments{
  \item{formula}{A formula or extended formula, i.e., the \code{formula} can be a
    \code{\link{list}} of formulas where each list entry specifies the details of one parameter
    of the modeled response distribution, see \code{\link{bamlss.formula}}. For incorporating
    smooth terms, all model term constructors implemented in \code{\link{mgcv}} such as
    \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}} and \code{\link[mgcv]{ti}} can be used, amongst others.}
  \item{data}{A \code{\link{data.frame}} or \code{\link{list}} containing the model
    response variable(s) and covariates specified in the \code{formula}.
    By default the variables are taken from \code{environment(formula)}:
    typically the environment from which \code{bamlss} is called.}
  \item{family}{A \code{\link{bamlss.family}} object, specifying the details of the modeled
    distribution such as the parameter names, the density function, link functions, etc.}
  \item{weights}{Prior weights on the data.}
  \item{subset}{An optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{offset}{Can be used to supply model offsets for use in fitting.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}'s. The default is set by the \code{na.action}
    setting of \code{\link{options}}, and is \code{\link{na.omit}} if set to \code{NULL}.}
  \item{contrasts}{An optional list. See the \code{contrasts.arg} of 
    \code{\link[stats]{model.matrix.default}}.}
  \item{knots}{An optional list containing user specified knots, see the documentation of
    function \code{\link[mgcv]{gam}}.}
  \item{specials}{Specify new special terms here to be used with the \code{\link{bamlss.formula}},
    see also \code{\link{terms.object}}.}
  \item{reference}{A \code{character} specifying a reference category, e.g., when
    fitting a multinomial model.}
  \item{model.matrix}{Logical, should model matrices for linear parts be returned?}
  \item{smooth.construct}{Logical, should model matrices, e.g., as returned from
    \code{\link[mgcv]{smooth.construct}} and \code{\link[mgcv]{smoothCon}} be part
    of returned \code{bamlss.frame}?.}
  \item{ytype}{For categorical responses, should the response be a vector or matrix. If
    \code{ytype == "matrix"} code{bamlss.frame()} uses function \code{\link{model.matrix}}
    to construct the response matrix from levels. If the response is a factor
    \code{ytype == "integer"} will create an integer response.}
  \item{scale.x}{Logical, should the model matrices of the linear parts be scaled?}
  \item{scale.d}{Logical, should the numeric variables in the model frame be scaled?}
  \item{\dots}{Arguments passed to function \code{\link{smooth.construct.bamlss.frame}}.}
}

\details{
  The function parses the \code{data}, the \code{formula} or the extended
  \code{\link{bamlss.formula}} as well as the \code{\link{bamlss.family}} into a model frame
  like object, the \code{\link{bamlss.frame}}. This object holds all necessary model matrices
  and information that is needed for model fitting engines. Per default,
  all package \code{\link{mgcv}} smooth term constructor functions like
  \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}}, \code{\link[mgcv]{t2}} and
  \code{\link[mgcv]{ti}} can be used (see also function \code{\link[mgcv]{smooth.construct}}),
  however, even special user defined constructors can be included, see the examples below.

  Function \code{bamlss.frame()} uses function \code{\link{model.matrix.bamlss.frame}} to
  compute all design matrices for simple linear parts, all smooth terms are parsed with
  function \code{\link{smooth.construct.bamlss.frame}}.

  It is also possible to create a \code{"bamlss.frame"} using hierarchical formulae, see the
  example below.
}

\value{
  An list of class \code{"bamlss.frame"} with the following elements:
  \item{call}{The initial call.}
  \item{model.frame}{The \code{\link{model.frame}} used to compute all design matrices.}
  \item{formula}{The \code{\link{bamlss.formula}}.}
  \item{family}{The \code{\link{bamlss.family}} object.}
  \item{terms}{The \code{\link{terms.bamlss}} object.}
  \item{x}{A named list, the elements correspond to the parameters that are specified
    within the \code{\link{bamlss.family}} object. For each parameter the corresponding
    \code{formula}, a \code{fake.formula} only holding the covariate names, a \code{\link{terms}}
    object, a \code{\link{model.matrix}} for the linear part and a list \code{smooth.construct}
    holding all information for smooth terms as returned from function
    \code{link{smooth.construct.bamlss.frame}} is created.}
  \item{y}{The response data.}
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.formula}}, \code{\link{bamlss.family}},
  \code{\link{smooth.construct.bamlss.frame}}, \code{\link{model.matrix.bamlss.frame}}
}

\examples{
## Create a 'bamlss.frame'.
d <- GAMart()
f <- list(
  num ~ fac + s(x1) + s(x2) + te(lon, lat),
  sigma ~ id + s(x2) + s(x3)
)
bf <- bamlss.frame(f, data = d, family = "gaussian")

## Show parts of the 'bamlss.frame'.
print(bf)

## Categorical responses.
f <- list(
  cat ~ fac + s(x1) + s(x2)
)

bf <- bamlss.frame(f, data = d, family = "multinomial", reference = "low")
print(bf)

## The response is a matrix per default.
head(bf$y)

## 0/1 responses.
d <- cbind(d, model.matrix(~ -1 + cat, data = d))

f <- list(
  catnone ~ fac + s(x1),
  catlow ~ s(x2),
  catmedium ~ s(x3)
)

bf <- bamlss.frame(f, data = d, family = "multinomial")
print(bf)

## Hierarchical structures.
f <- list(
  num ~ s(x1) + s(x2) + id,
  id ~ te(lon, lat),
  sigma ~ s(x1) + fac
)

bf <- bamlss.frame(f, data = d, family = "gaussian")
print(bf)

## Special model term constructors,
## set up "new" constructor function and eval
## with bamlss.frame().
s77 <- function(...) {
  sm <- s(...)
  sm$label <- paste("s77(", paste(sm$term, collapse = ","), ")", sep = "")
  sm
}

f <- list(
  num ~ s77(x1) + s(x2) + id,
  sigma ~ s77(x1)
)

bf <- bamlss.frame(f, data = d, family = "gaussian", specials = "s77")
print(bf)
names(bf$x$mu$smooth.construct)
}

\keyword{regression}
\keyword{MCMC}

