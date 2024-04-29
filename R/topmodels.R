#' Create a BAMLSS Distribution
#'
#' A single class and corresponding methods encompassing all \code{bamlss.family}
#' distributions (from the \pkg{bamlss} package) using the workflow from the
#' \pkg{distributions3} package.
#' 
#' The constructor function \code{BAMLSS} sets up a distribution
#' object, representing a distribution from the BAMLSS (Bayesian additive
#' model of location, scale, and shape) framework by the corresponding parameters
#' plus a \code{family} attribute, e.g., \code{\link[bamlss]{gaussian_bamlss}} for the
#' normal distribution or \code{\link[bamlss]{binomial_bamlss}} for the binomial
#' distribution. The parameters employed by the family vary across the families
#' but typically capture different distributional properties (like location, scale,
#' shape, etc.).
#'
#' All parameters can also be vectors, so that it is possible to define a vector
#' of BAMLSS distributions from the same family with potentially different parameters.
#' All parameters need to have the same length or must be scalars (i.e.,
#' of length 1) which are then recycled to the length of the other parameters.
#' 
#' For the \code{BAMLSS} distribution objects there is a wide range
#' of standard methods available to the generics provided in the \pkg{distributions3}
#' package: \code{\link[distributions3]{pdf}} and \code{\link[distributions3]{log_pdf}}
#' for the (log-)density (PDF), \code{\link[distributions3]{cdf}} for the probability
#' from the cumulative distribution function (CDF), \code{quantile} for quantiles,
#' \code{\link[distributions3]{random}} for simulating random variables,
#' and \code{\link[distributions3]{support}} for the support interval
#' (minimum and maximum). Internally, these methods rely on the usual d/p/q/r
#' functions provided in \pkg{bamlss}, see the manual pages of the individual
#' families. The methods \code{\link[distributions3]{is_discrete}} and
#' \code{\link[distributions3]{is_continuous}} can be used to query whether the
#' distributions are discrete on the entire support or continuous on the entire
#' support, respectively.
#'
#' See the examples below for an illustration of the workflow for the class and methods.
#' 
#' @param family object. BAMLSS family specifications recognized by
#'   \code{\link[bamlss]{bamlss.family}}, including \code{family.bamlss} objects,
#'   family-generating functions (e.g., \code{\link[bamlss]{gaussian_bamlss}}),
#'   or characters with family names (e.g., \code{"gaussian"} or \code{"binomial"}).
#' @param \dots further arguments passed as parameters to the BAMLSS family.
#'   Can be scalars or vectors.
#' 
#' @return A \code{BAMLSS} distribution object.
#' 
#' @seealso \code{\link[bamlss]{bamlss.family}}
#' 
#' @examples
#' \dontshow{ if(!requireNamespace("bamlss")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("not all packages required for the example are installed")
#'   } else q() }
#' }
#' ## package and random seed
#' library("distributions3")
#' set.seed(6020)
#' 
#' ## three Weibull distributions
#' X <- BAMLSS("weibull", lambda = c(1, 1, 2), alpha = c(1, 2, 2))
#' X
#' 
#' ## moments (FIXME: mean and variance not provided by weibull_bamlss)
#' ## mean(X)
#' ## variance(X)
#' 
#' ## support interval (minimum and maximum)
#' support(X)
#' is_discrete(X)
#' is_continuous(X)
#' 
#' ## simulate random variables
#' random(X, 5)
#' 
#' ## histograms of 1,000 simulated observations
#' x <- random(X, 1000)
#' hist(x[1, ], main = "Weibull(1,1)")
#' hist(x[2, ], main = "Weibull(1,2)")
#' hist(x[3, ], main = "Weibull(2,2)")
#' 
#' ## probability density function (PDF) and log-density (or log-likelihood)
#' x <- c(2, 2, 1)
#' pdf(X, x)
#' pdf(X, x, log = TRUE)
#' log_pdf(X, x)
#' 
#' ## cumulative distribution function (CDF)
#' cdf(X, x)
#' 
#' ## quantiles
#' quantile(X, 0.5)
#' 
#' ## cdf() and quantile() are inverses
#' cdf(X, quantile(X, 0.5))
#' quantile(X, cdf(X, 1))
#' 
#' ## all methods above can either be applied elementwise or for
#' ## all combinations of X and x, if length(X) = length(x),
#' ## also the result can be assured to be a matrix via drop = FALSE
#' p <- c(0.05, 0.5, 0.95)
#' quantile(X, p, elementwise = FALSE)
#' quantile(X, p, elementwise = TRUE)
#' quantile(X, p, elementwise = TRUE, drop = FALSE)
#' 
#' ## compare theoretical and empirical mean from 1,000 simulated observations
#' ## (FIXME: mean not provided by weibull_bamlss)
#' ## cbind(
#' ##   "theoretical" = mean(X),
#' ##   "empirical" = rowMeans(random(X, 1000))
#' ## )
#' @export
BAMLSS <- function(family, ...) {
  stopifnot(requireNamespace("bamlss"))
  ## get family object
  f <- bamlss::bamlss.family(family)

  ## collect parameters in data.frame
  d <- data.frame(...)

  ## check whether all necessary parameters are specified
  n <- setdiff(f$names, names(d))
  if(length(n) > 0L) stop(sprintf("the parameter%s %s of the '%s' family %s not specified",
    if(length(n) > 1L) "s" else "", paste("'", n, "'", collapse = ", "), f$family, if(length(n) > 1L) "are" else "is"))
  n <- setdiff(names(d), f$names)
  if(length(n) > 0L) warning(sprintf("%s %s of the %s family",
    paste0("'", n, "'", collapse = ", "), if(length(n) > 1L) "are not parameters" else "is not a parameter", f$family))
  
  ## set up distribution
  d <- d[, f$names, drop = FALSE]
  class(d) <- c("BAMLSS", "distribution")
  attr(d, "family") <- f
  return(d)
}

#' @rdname BAMLSS
#' @method family BAMLSS
#' @export
#' @usage NULL
#' @importFrom stats family
family.BAMLSS <- function(object, ...) {
  attr(object, "family")
}

#' @rdname BAMLSS
#' @method mean BAMLSS
#' @export
#' @usage NULL
#' @importFrom stats setNames
mean.BAMLSS <- function(x, ...) {
  f <- family(x)
  if(!("mean" %in% names(f))) stop(sprintf("no mean() function provided by '%s' family", f$family))
  m <- f$mean(x)
  setNames(m, names(x))
}

#' @rdname BAMLSS
#' @method variance BAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 variance
#' @importFrom stats setNames
variance.BAMLSS <- function(x, ...) {
  f <- family(x)
  if(!("variance" %in% names(f))) stop(sprintf("no variance() function provided by '%s' family", f$family))
  m <- f$variance(x)
  setNames(m, names(x))
}

#' @rdname BAMLSS
#' @method skewness BAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 skewness
skewness.BAMLSS <- function(x, ...) {
  f <- family(x)
  if(!("skewness" %in% names(f))) stop(sprintf("no skewness() function provided by '%s' family", f$family))
  m <- f$skewness(x)
  setNames(m, names(x))
}

#' @rdname BAMLSS
#' @method kurtosis BAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 kurtosis
kurtosis.BAMLSS <- function(x, ...) {
  f <- family(x)
  if(!("kurtosis" %in% names(f))) stop(sprintf("no kurtosis() function provided by '%s' family", f$family))
  m <- f$kurtosis(x)
  setNames(m, names(x))
}

#' @rdname BAMLSS
#' @method random BAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 random apply_dpqr make_positive_integer
random.BAMLSS <- function(x, n = 1L, drop = TRUE, ...) {
  f <- family(x)
  if(!("r" %in% names(f))) stop(sprintf("no r() function provided by '%s' family", f$family))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) f$r(n = at, par = d, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

#' @rdname BAMLSS
#' @importFrom distributions3 pdf apply_dpqr
#' @method pdf BAMLSS
#' @export
#' @usage NULL
pdf.BAMLSS <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  f <- family(d)
  if(!("d" %in% names(f))) stop(sprintf("no d() function provided by '%s' family", f$family))
  FUN <- function(at, d) f$d(y = at, par = d, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

#' @rdname BAMLSS
#' @method log_pdf BAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 log_pdf apply_dpqr
log_pdf.BAMLSS <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  f <- family(d)
  if(!("d" %in% names(f))) stop(sprintf("no d() function provided by '%s' family", f$family))
  FUN <- function(at, d) f$d(y = at, par = d, log = TRUE, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

#' @rdname BAMLSS
#' @method cdf BAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 cdf apply_dpqr
cdf.BAMLSS <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  f <- family(d)
  if(!("p" %in% names(f))) stop(sprintf("no p() function provided by '%s' family", f$family))
  FUN <- function(at, d) f$p(y = at, par = d, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

#' @rdname BAMLSS
#' @method quantile BAMLSS
#' @export
#' @usage NULL
#' @importFrom stats quantile
#' @importFrom distributions3 apply_dpqr
quantile.BAMLSS <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  f <- family(x)
  if(!("q" %in% names(f))) stop(sprintf("no q() function provided by '%s' family", f$family))
  FUN <- function(at, d) f$q(p = at, par = d, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

#' @rdname BAMLSS
#' @method support BAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 support make_support
support.BAMLSS <- function(d, drop = TRUE, ...) {
  s <- quantile(d, probs = c(0, 1), elementwise = FALSE)
  distributions3::make_support(s[, 1L], s[, 2L], d, drop = drop)
}

#' @rdname BAMLSS
#' @method is_discrete BAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 is_discrete
#' @importFrom stats setNames
is_discrete.BAMLSS <- function(d, ...) {
  f <- family(d)
  if(!("type" %in% names(f))) {
    warning(sprintf("no 'type' information provided by '%s' family", f$family))
    f$type <- "unknown"
  }
  setNames(rep.int(f$type == "discrete", length(d)), names(d))
}

#' @rdname BAMLSS
#' @method is_continuous BAMLSS
#' @export
#' @usage NULL
#' @importFrom distributions3 is_continuous
#' @importFrom stats setNames
is_continuous.BAMLSS <- function(d, ...) {
  f <- family(d)
  if(!("type" %in% names(f))) {
    warning(sprintf("no 'type' information provided by '%s' family", f$family))
    f$type <- "unknown"
  }
  setNames(rep.int(f$type == "continuous", length(d)), names(d))
}

#' @rdname BAMLSS
#' @method format BAMLSS
#' @export
#' @usage NULL
format.BAMLSS <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
  class(x) <- c(paste("BAMLSS", family(x)$family), "distribution")
  NextMethod()
}

#' @rdname BAMLSS
#' @method print BAMLSS
#' @export
#' @usage NULL
print.BAMLSS <- function(x, digits = pmax(3L, getOption("digits") - 3L), ...) {
  class(x) <- c(paste("BAMLSS", family(x)$family), "distribution")
  NextMethod()
}



#' Extracting Fitted or Predicted Probability Distributions from bamlss Models
#' 
#' Methods for \pkg{bamlss} model objects for extracting
#' fitted (in-sample) or predicted (out-of-sample) probability distribution
#' objects.
#' 
#' To facilitate making probabilistic forecasts based on \code{\link[bamlss]{bamlss}}
#' model objects, the \code{\link[distributions3]{prodist}} method extracts fitted or
#' predicted probability \code{distribution} objects. Internally, the
#' \code{\link[bamlss]{predict.bamlss}} method from the \pkg{bamlss} package is
#' used first to obtain the distribution parameters. Subsequently, the corresponding \code{distribution}
#' object is set up using \code{\link{BAMLSS}}, enabling the workflow provided by
#' the \pkg{distributions3} package.
#' 
#' Note that these probability distributions only reflect the random variation in
#' the dependent variable based on the model employed (and its associated
#' distributional assumption for the dependent variable). This does not capture
#' the uncertainty in the parameter estimates.
#' 
#' @param object A model object of class \code{\link[bamlss]{bamlss}}.
#' @param ... Arguments passed on to \code{\link[bamlss]{predict.bamlss}}, 
#' e.g., \code{newdata}.
#' @param distributions3 logical. If a dedicated \pkg{distributions3} object
#' is available (e.g., such as \code{\link[distributions3]{Normal}}) and uses
#' the same parameterization, should this be used instead of the general
#' \code{BAMLSS} distribution?
#' 
#' @return An object inheriting from \code{distribution}.
#' 
#' @seealso \code{\link{BAMLSS}} \code{\link[bamlss]{predict.bamlss}}
#' 
#' @keywords distribution
#' 
#' @examples
#' \dontshow{ if(!requireNamespace("bamlss")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("not all packages required for the example are installed")
#'   } else q() }
#' }
#' ## packages, code, and data
#' library("bamlss")
#' library("distributions3")
#' data("cars", package = "datasets")
#' 
#' ## fit heteroscedastic normal BAMLSS model
#' f <- list(dist ~ s(speed), sigma ~ s(speed))
#' m <- bamlss(f, data = cars, family = "gaussian")
#' 
#' ## obtain predicted distributions for three levels of speed
#' d <- prodist(m, newdata = data.frame(speed = c(10, 20, 30)))
#' print(d)
#' 
#' ## obtain quantiles (works the same for any distribution object 'd' !!)
#' quantile(d, 0.5)
#' quantile(d, c(0.05, 0.5, 0.95), elementwise = FALSE)
#' quantile(d, c(0.05, 0.5, 0.95), elementwise = TRUE)
#' 
#' ## visualization
#' plot(dist ~ speed, data = cars)
#' nd <- data.frame(speed = 0:240/4)
#' nd$dist <- prodist(m, newdata = nd)
#' nd$fit <- quantile(nd$dist, c(0.05, 0.5, 0.95))
#' matplot(nd$speed, nd$fit, type = "l", lty = 1, col = "slategray", add = TRUE)
#' 
#' ## moments
#' mean(d)
#' variance(d)
#' 
#' ## simulate random numbers
#' random(d, 5)
#' 
#' ## density and distribution
#' pdf(d, 50 * -2:2)
#' cdf(d, 50 * -2:2)
#' 
#' ## further diagnostics: graphical and scores
#' pithist(m)
#' qqrplot(m)
#' proscore(m, type = c("LogLik", "CRPS", "MAE", "MSE"), aggregate = TRUE)
#' 
#' ## note that proscore can replicate logLik() value
#' proscore(m, aggregate = sum)
#' logLik(m)
#' 
#' ## Poisson example
#' data("FIFA2018", package = "distributions3")
#' m2 <- bamlss(goals ~ s(difference), data = FIFA2018, family = "poisson")
#' d2 <- prodist(m2, newdata = data.frame(difference = 0))
#' print(d2)
#' quantile(d2, c(0.05, 0.5, 0.95))
#' @export
prodist.bamlss <- function(object, ..., distributions3 = FALSE) {
  d <- predict(object, type = "parameter", drop = FALSE, ...)
  ## check whether distributions3 object is already available
  if(distributions3) {
    rval <- switch(object$family$family,
      "binomial"  = distributions3::Binomial(p = d$pi, size = 1),
      "gaussian"  = distributions3::Normal(mu = d$mu, sigma = d$sigma),
      "lognormal" = distributions3::LogNormal(log_mu = d$mu, log_sigma = d$sigma),
      "gpareto"   = distributions3::GP(xi = d$xi, sigma = d$sigma),
      "GEV"       = distributions3::GEV(mu = d$mu, sigma = d$sigma, xi = d$xi),
      "weibull"   = distributions3::Weibull(shape = d$alpha, scale = d$lambda),
      "poisson"   = distributions3::Poisson(lambda = d$lambda),
      "nbinom"    = distributions3::NegativeBinomial(mu = d$mu, size = d$theta),
      "ztnbinom"  = distributions3::ZTNegativeBinomial(mu = d$mu, theta = d$theta),
      NULL
    )
    if(is.null(rval)) {
      warning(sprintf("no dedicated distributions3 object available for '%s' family, using general BAMLSS distribution object instead", object$family$family))
    } else {
      return(rval)
    }
  }
  ## otherwise use general BAMLSS distributions3 object
  do.call("BAMLSS", c(list(family = object$family), d))
}

#' @rdname crps.distribution
#' @exportS3Method scoringRules::crps BAMLSS
crps.BAMLSS <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  if(requireNamespace("scoringRules")) {
    ## manually match bamlss family distribution names with scoringRules, if possible
    f <- family(y)$family
    FUN <- switch(EXPR = f,
      "binomial"  = function(at, d) scoringRules::crps_binom(y = at, prob = d$pi, size = 1L),
      "gaussian"  = function(at, d) scoringRules::crps_norm(y = at, mean = d$mu, sd = d$sigma),
      "lognormal" = function(at, d) scoringRules::crps_lnorm(y = at, meanlog = d$mu, sdlog = d$sigma),
      "GEV"       = function(at, d) scoringRules::crps_gev(y = at, location = d$mu, scale = d$sigma, shape = d$xi),
      "gpareto"   = function(at, d) scoringRules::crps_gpd(y = at, scale = d$sigma, shape = d$xi),
      "poisson"   = function(at, d) scoringRules::crps_pois(y = at, lambda = d$lambda),
      "nbinom"    = function(at, d) scoringRules::crps_nbinom(y = at, size = d$theta, mu = d$mu),
      NULL ## FIXME: bamlss and scoringRules can also be matched with different parameterizations for beta, gamma, Gumbel, etc.
    )
  } else {
    FUN <- NULL
  }
  if(is.null(FUN)) {
    ## use crps.distribution() if no closed-form solution from scoringRules is available
    NextMethod()
  } else {
    ## use apply_dpqr() with scoringRules::crps_*() function
    distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
  }
}
