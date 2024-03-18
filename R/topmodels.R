BAMLSS <- function(family, ...) {
  stopifnot(requireNamespace("bamlss"))
  f <- bamlss::bamlss.family(family)
  d <- data.frame(...)
  n <- setdiff(f$names, names(d))
  if (length(n) > 0) 
 stop(sprintf("the parameter%s %s of the '%s' family %s not specified", if (length(n) > 1) 
   "s"
 else "", paste("'", n, "'", collapse = ", "), f$family, if (length(n) > 1) 
   "are"
 else "is"))
  n <- setdiff(names(d), f$names)
  if (length(n) > 0) 
 warning(sprintf("%s %s of the %s family", paste0("'", n, "'", collapse = ", "), if (length(n) > 1) 
   "are not parameters"
 else "is not a parameter", f$family))
  d <- d[, f$names, drop = FALSE]
  class(d) <- c("BAMLSS", "distribution")
  attr(d, "family") <- f
  return(d)
}
family.BAMLSS <- function(object, ...) {
  attr(object, "family")
}
mean.BAMLSS <- function(x, ...) {
  f <- family(x)
  if (!("mean" %in% names(f))) 
 stop(sprintf("no mean() function provided by '%s' family", f$family))
  m <- f$mean(x)
  setNames(m, names(x))
}
variance.BAMLSS <- function(x, ...) {
  f <- family(x)
  if (!("variance" %in% names(f))) 
 stop(sprintf("no variance() function provided by '%s' family", f$family))
  m <- f$variance(x)
  setNames(m, names(x))
}
skewness.BAMLSS <- function(x, ...) {
  f <- family(x)
  if (!("skewness" %in% names(f))) 
 stop(sprintf("no skewness() function provided by '%s' family", f$family))
  m <- f$skewness(x)
  setNames(m, names(x))
}
kurtosis.BAMLSS <- function(x, ...) {
  f <- family(x)
  if (!("kurtosis" %in% names(f))) 
 stop(sprintf("no kurtosis() function provided by '%s' family", f$family))
  m <- f$kurtosis(x)
  setNames(m, names(x))
}
random.BAMLSS <- function(x, n = 1, drop = TRUE, ...) {
  f <- family(x)
  if (!("r" %in% names(f))) 
 stop(sprintf("no r() function provided by '%s' family", f$family))
  n <- distributions3::make_positive_integer(n)
  if (n == 0) 
 return(numeric(0))
  FUN <- function(at, d) f$r(n = at, par = d, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}
pdf.BAMLSS <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  f <- family(d)
  if (!("d" %in% names(f))) 
 stop(sprintf("no d() function provided by '%s' family", f$family))
  FUN <- function(at, d) f$d(y = at, par = d, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}
log_pdf.BAMLSS <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  f <- family(d)
  if (!("d" %in% names(f))) 
 stop(sprintf("no d() function provided by '%s' family", f$family))
  FUN <- function(at, d) f$d(y = at, par = d, log = TRUE, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}
cdf.BAMLSS <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  f <- family(d)
  if (!("p" %in% names(f))) 
 stop(sprintf("no p() function provided by '%s' family", f$family))
  FUN <- function(at, d) f$p(y = at, par = d, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}
quantile.BAMLSS <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  f <- family(x)
  if (!("q" %in% names(f))) 
 stop(sprintf("no q() function provided by '%s' family", f$family))
  FUN <- function(at, d) f$q(p = at, par = d, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}
support.BAMLSS <- function(d, drop = TRUE, ...) {
  s <- quantile(d, probs = c(0, 1), elementwise = FALSE)
  distributions3::make_support(s[, 1], s[, 2], d, drop = drop)
}
is_discrete.BAMLSS <- function(d, ...) {
  f <- family(d)
  if (!("type" %in% names(f))) {
 warning(sprintf("no 'type' information provided by '%s' family", f$family))
 f$type <- "unknown"
  }
  setNames(rep.int(f$type == "discrete", length(d)), names(d))
}
is_continuous.BAMLSS <- function(d, ...) {
  f <- family(d)
  if (!("type" %in% names(f))) {
 warning(sprintf("no 'type' information provided by '%s' family", f$family))
 f$type <- "unknown"
  }
  setNames(rep.int(f$type == "continuous", length(d)), names(d))
}
format.BAMLSS <- function(x, digits = pmax(3, getOption("digits") - 3), ...) {
  class(x) <- c(paste("BAMLSS", family(x)$family), "distribution")
  NextMethod()
}
print.BAMLSS <- function(x, digits = pmax(3, getOption("digits") - 3), ...) {
  class(x) <- c(paste("BAMLSS", family(x)$family), "distribution")
  NextMethod()
}
prodist.bamlss <- function(object, ..., distributions3 = FALSE) {
  d <- predict(object, type = "parameter", drop = FALSE, ...)
  if (distributions3) {
 rval <- switch(object$family$family, binomial = distributions3::Binomial(p = d$pi, size = 1), gaussian = distributions3::Normal(mu = d$mu, sigma = d$sigma), lognormal = distributions3::LogNormal(log_mu = d$mu, log_sigma = d$sigma), gpareto = distributions3::GP(xi = d$xi, sigma = d$sigma), GEV = distributions3::GEV(mu = d$mu, sigma = d$sigma, xi = d$xi), weibull = distributions3::Weibull(shape = d$alpha, scale = d$lambda), poisson = distributions3::Poisson(lambda = d$lambda), nbinom = distributions3::NegativeBinomial(mu = d$mu, 
   size = d$theta), ztnbinom = distributions3::ZTNegativeBinomial(mu = d$mu, theta = d$theta), NULL)
 if (is.null(rval)) {
   warning(sprintf("no dedicated distributions3 object available for '%s' family, using general BAMLSS distribution object instead", object$family$family))
 }
 else {
   return(rval)
 }
  }
  do.call("BAMLSS", c(list(family = object$family), d))
}
crps.BAMLSS <- function(y, x, drop = TRUE, elementwise = NULL, ...) {
  if (requireNamespace("scoringRules")) {
 f <- family(y)$family
 FUN <- switch(EXPR = f, binomial = function(at, d) scoringRules::crps_binom(y = at, prob = d$pi, size = 1), gaussian = function(at, d) scoringRules::crps_norm(y = at, mean = d$mu, sd = d$sigma), lognormal = function(at, d) scoringRules::crps_lnorm(y = at, meanlog = d$mu, sdlog = d$sigma), GEV = function(at, d) scoringRules::crps_gev(y = at, location = d$mu, scale = d$sigma, shape = d$xi), gpareto = function(at, d) scoringRules::crps_gpd(y = at, scale = d$sigma, shape = d$xi), poisson = function(at, 
   d) scoringRules::crps_pois(y = at, lambda = d$lambda), nbinom = function(at, d) scoringRules::crps_nbinom(y = at, size = d$theta, mu = d$mu), NULL)
  }
  else {
 FUN <- NULL
  }
  if (is.null(FUN)) {
 NextMethod()
  }
  else {
 distributions3::apply_dpqr(d = y, FUN = FUN, at = x, type = "crps", drop = drop, elementwise = elementwise)
  }
}
