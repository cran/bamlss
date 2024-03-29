\name{family.bamlss}
\alias{family.bamlss}
\alias{family.bamlss.frame}
\alias{bamlss.family}
\alias{beta_bamlss}
\alias{binomial_bamlss}
\alias{dw_bamlss}
\alias{cnorm_bamlss}
\alias{cox_bamlss}
\alias{gaussian_bamlss}
\alias{gaussian2_bamlss}
\alias{Gaussian_bamlss}
\alias{gamma_bamlss}
\alias{multinomial_bamlss}
\alias{mvnorm_bamlss}
\alias{mvnormAR1_bamlss}
\alias{poisson_bamlss}
\alias{gpareto_bamlss}
\alias{glogis_bamlss}
\alias{AR1_bamlss}
\alias{beta1_bamlss}
\alias{ztnbinom_bamlss}
\alias{lognormal_bamlss}
\alias{ALD_bamlss}
\alias{nbinom_bamlss}
\alias{weibull_bamlss}
\alias{dirichlet_bamlss}
\alias{logNN_bamlss}
\alias{Sichel_bamlss}
\alias{GEV_bamlss}
\alias{gumbel_bamlss}
\alias{mix_bamlss}
\alias{ELF_bamlss}
\alias{DGP_bamlss}
\alias{ZANBI_bamlss}

\title{
  Distribution Families in \pkg{bamlss}
}

\description{
  Family objects in \pkg{bamlss} specify the information that is needed for using (different)
  model fitting engines, e.g., the parameter names and corresponding link functions, the
  density function, derivatives of the log-likelihood w.r.t. the predictors, and so
  forth. The optimizer or sampler functions that are called by \code{\link{bamlss}} must know how
  much information is needed to interpret the model since the family objects are simply
  passed through. Family objects are also used for computing post-modeling statistics, e.g., for
  residual diagnostics or random number generation. See the details and examples.
}

\usage{
## Family objects in bamlss:
ALD_bamlss(..., tau = 0.5, eps = 0.01)
beta_bamlss(...)
binomial_bamlss(link = "logit", ...)
cnorm_bamlss(...)
cox_bamlss(...)
dw_bamlss(...)
DGP_bamlss(...)
dirichlet_bamlss(...)
ELF_bamlss(..., tau = 0.5)
gaussian_bamlss(...)
gaussian2_bamlss(...)
Gaussian_bamlss(...)
gamma_bamlss(...)
logNN_bamlss(...)
multinomial_bamlss(...)
mvnorm_bamlss(k = 2, ...)
mvnormAR1_bamlss(k = 2, ...)
poisson_bamlss(...)
gpareto_bamlss(...)
glogis_bamlss(...)
AR1_bamlss(...)
beta1_bamlss(ar.start, ...)
nbinom_bamlss(...)
ztnbinom_bamlss(...)
lognormal_bamlss(...)
weibull_bamlss(...)
Sichel_bamlss(...)
GEV_bamlss(...)
gumbel_bamlss(...)
mix_bamlss(f1, f2, ...)
ZANBI_bamlss(...)

## Extractor functions:
\method{family}{bamlss}(object, ...)
\method{family}{bamlss.frame}(object, ...)
}

\arguments{
  \item{object}{An object of class \code{"bamlss"} or \code{"bamlss.frame"}, see function
    \code{\link{bamlss}} and \code{\link{bamlss.frame}}.}
  \item{k}{The dimension of the multivariate normal. Note, if \code{k = 1} function
    \code{gaussian_bamlss()} is called.}
  \item{ar.start}{Logical vector of length equal to the number of rows of the full data set
    used for modeling. Must hold entries \code{TRUE} indicating the start of a time series of
    a section. If \code{ar.start = NULL} lagged residuals are computed by simple shifting.
    See also \code{\link[mgcv]{bam}}.}
  \item{link}{Possible link functions.}
  \item{tau}{The quantile the should be fitted.}
  \item{eps}{Constant to be used for the approximation of the absolute function.}
  \item{f1, f2}{A family of class \code{"gamlss.family"}, see package \pkg{gamlss.dist}.}
  \item{\dots}{Arguments passed to functions that are called within the family object.}
}

\details{
  The following lists the minimum requirements on a \pkg{bamlss} family object to be used with
  \code{\link{bamlss}} and \code{\link{bamlss.frame}}:
  \itemize{
    \item The family object must return a \code{\link{list}} of class \code{"family.bamlss"}.
    \item The object must contain the family name as a character string. 
    \item The object must contain the names of the parameters as a character string, as well as
      the corresponding link functions as character string.
  }
  For most optimizer and sampling functions at least the density function, including a \code{log}
  argument, should be provided.
  When using generic model fitting engines like \code{\link{opt_bfit}} or \code{\link{sam_GMCMC}}, as
  well as for computing post-modeling statistics with function \code{\link{samplestats}}, and others,
  it is assumed that the density function in a family object has the following arguments:

  \code{d(y, par, log = FALSE, ...)}

  where argument \code{y} is the response (possibly a matrix) and \code{par} is a named list
  holding the evaluated parameters of the distribution, e.g., using a normal distribution \code{par}
  has two elements, one for the mean \code{par$mu} and one for the standard deviation
  \code{par$sigma}. The dots argument is for passing special internally used objects, depending
  on the type of model this feature is usually not needed.

  Similarly, for derivative based algorithms, e.g. using iteratively weighted least squares
  (IWLS, see function \code{\link{opt_bfit}}, the family object holds derivative functions evaluating
  derivatives of the log-likelihood w.r.t. the predictors (or expectations of derivatives).
  For each parameter, these functions also hold the following arguments:
 
  \code{score(y, par, ...)}

  for computing the first derivative of the log-likelihood w.r.t. a predictor and

  \code{hess(y, par, ...)}

  for computing the negative second derivatives. Within the family object these functions
  are organized in a named list, see the examples below.

  In addition, for the cumulative distribution function (\code{p(y, par, ...)}), for the quantile
  function (\code{q(y, par, ...)}) or for creating random numbers (\code{r(n, par, ...)}) the same
  structure is assumed. See, e.g., the code of function \code{gaussian.bamlss()}.

  Some model fitting engines can initialize the distributional parameters which oftentimes
  leads to much faster convergence. The initialize functions are again organized within a named list,
  one entry for each parameter, similar to the \code{score} and \code{hess} functions, e.g.,
  see the code of family object \code{gaussian.bamlss}.

  Using function \code{\link{bamlss}}, \code{\link{residuals.bamlss}} and
  \code{\link{predict.bamlss}} the family objects may
  also specify the \code{transform()}er, \code{optimizer()}, \code{sampler()},
  \code{samplestats()}, \code{results()}, \code{residuals()} and \code{predict()} function that
  should be used with this family. See for example the setup of \code{\link{cox_bamlss}}.

  For using specialized estimation engines like \code{\link{sam_JAGS}} it is recommended to supply
  any extra arguments needed by those engines with an additional list entry within the family
  object, e.g., using \code{gaussian_bamlss()} with \code{\link{sam_JAGS}} the family objects holds
  special details in an element named \code{"bugs"}.

  The examples below illustrate this setup. See also the code of the \pkg{bamlss} family
  functions.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}}
}

\examples{
## New family object for the normal distribution,
## can be used by function opt_bfit() and sam_GMCMC().
normal_bamlss <- function(...) {
  f <- list(
    "family" = "normal",
    "names" = c("mu", "sigma"),
    "links" = c("identity", "log"),
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        drop((y - par$mu) / (par$sigma^2))
      },
      "sigma" = function(y, par, ...) {
        drop(-1 + (y - par$mu)^2 / (par$sigma^2))
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        drop(1 / (par$sigma^2))
      },
      "sigma" = function(y, par, ...) { 
        rep(2, length(y))
      }
    )
  )
  class(f) <- "family.bamlss"
  return(f)
}

\dontrun{## Test on simulated data.
d <- GAMart()
b <- bamlss(num ~ s(x1) + s(x2) + s(x3),
  data = d, family = "normal")
plot(b)

## Compute the log-likelihood using the family object.
f <- family(b)
sum(f$d(y = d$num, par = f$map2par(fitted(b)), log = TRUE))

## For using JAGS() more details are needed.
norm4JAGS_bamlss <- function(...) {
  f <- normal_bamlss()
  f$bugs <- list(
    "dist" = "dnorm",
    "eta" = BUGSeta,
    "model" = BUGSmodel,
    "reparam" = c(sigma = "1 / sqrt(sigma)")
  )
  return(f)
}

## Now with opt_bfit() and sam_JAGS().
b <- bamlss(num ~ s(x1) + s(x2) + s(x3), data = d,
  optimizer = opt_bfit, sampler = sam_JAGS, family = "norm4JAGS")
plot(b)
}
}

\keyword{regression}
\keyword{models}
\keyword{distribution}

