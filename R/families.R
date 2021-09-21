######################
## BAMLSS Families. ##
######################
## Print method.
print.family.bamlss <- function(x, full = TRUE, ...)
{
  cat("Family:", x$family, if(!is.null(x$full.name)) paste0("(", x$full.name, ")") else NULL,  "\n")
  links <- paste(names(x$links), x$links, sep = " = ")
  links <- paste(links, collapse = ", ")
  cat(if(length(links) > 1) "Link functions:" else "Link function:", links, sep = " ")
  cat("\n")
  if(full) {
    nfun <- names(x[c("transform", "optimizer", "sampler", "results", "predict")])
    if(!all(is.na(nfun))) {
      nfun <- nfun[!is.na(nfun)]
      cat("---\nFamily specific functions:\n")
      for(j in nfun)
        cat(" ..$ ", j, "\n", sep = "")
    }
    nfun <- names(x[c("score", "hess")])
    if(!all(is.na(nfun))) {
      nfun <- nfun[!is.na(nfun)]
      cat("---\nDerivative functions:\n")
      for(j in nfun) {
        cat(" ..$ ", j, "\n", sep = "")
        for(i in names(x[[j]]))
          cat(" .. ..$ ", i, "\n", sep = "")
      }
    }
  }
}


## Second make.link function.
make.link2 <- function(link)
{
  if(is.null(link))
    link <- "identity"
  link0 <- link
  if(link0 == "tanhalf"){
    rval <- list(
      "linkfun" = function (mu) {
        tan(mu/2)},
      "linkinv" = function(eta) {
        2 * atan(eta)},
      "mu.eta" = function(eta) {
        2 / (eta^2 + 1)},
      "mu.eta2" = function(eta) {
        (-4 * eta ) / (eta^2 + 1)^2},
      "valideta" = function(eta) TRUE,
      "name" = "tanhalf"
      )
  } else {
    mu.eta2 <- function(x) {
      if(link0 == "identity") {
        x$mu.eta2 <- function(eta) rep.int(0, length(eta))
        return(x)
      }
      if(link0 == "log") {
        x$mu.eta2 <- function(eta) exp(eta)
        return(x)
      }
      if(link0 == "logit") {
        x$mu.eta2 <- function(eta) {
          eta <- exp(eta)
          return(-eta * (eta - 1) / (eta + 1)^3)
        }
        return(x)
      }
      if(link0 == "probit") {
        x$mu.eta2 <- function(eta) {
          -eta * dnorm(eta, mean = 0, sd = 1)
        }
        return(x)
      }
      if(link0 == "inverse") {
        x$mu.eta2 <- function(eta) {
          2 / (eta^3)
        }
        return(x)
      }
      if(link0 == "1/mu^2") {
        x$mu.eta2 <- function(eta) {
          0.75 / eta^(2.5)
        }
        return(x)
      }
      if(link0 == "sqrt") {
        x$mu.eta2 <- function(eta) { rep(2, length = length(eta)) }
        return(x)
      }
      x$mu.eta2 <- function(eta) rep.int(0, length(eta))
      ## warning(paste('higher derivatives of link "', link, '" not available!', sep = ''))
      return(x)
    }

    if(link %in% c("logit", "probit", "cauchit", "cloglog", "identity",
                   "log", "sqrt", "1/mu^2", "inverse")) {
      rval <- make.link(link)
    } else {
      rval <- switch(link,
        "rhogit" = list(
          "linkfun" = function(mu) { mu / sqrt(1 - mu^2) },
          "linkinv" = function(eta) {
              rval <- eta / sqrt(1 + eta^2)
              rval <- (abs(rval) - .Machine$double.eps) * sign(rval)
              rval
          },
          "mu.eta" = function(eta) { 1 / (1 + eta^2)^1.5 }
        ),
        "cloglog2" = list(
          "linkfun" = function(mu) { log(-log(mu)) },
          "linkinv" = function(eta) {
            pmax(pmin(1 - expm1(-exp(eta)), .Machine$double.eps), .Machine$double.eps)
          },
          "mu.eta" = function(eta) {
            eta <- pmin(eta, 700)
            pmax(-exp(eta) * exp(-exp(eta)), .Machine$double.eps)
          }
        ),
        "sigmoid" = list(
          "linkfun" = function(mu) {
            i <- mu <= -1
            if(any(i))
              mu[i] <- mu[i] <- -0.9999
            i <- mu >= 1
            if(any(i))
              mu[i] <- mu[i] <- 0.9999 
            -log(2/(mu + 1) - 1)
          },
          "linkinv" = function(eta) {
            tanh(eta/2)
          },
          "mu.eta" = function(eta) {
            0.5 / cosh(eta * 0.5)^2
          },
          "mu.eta2" = function(eta) {
            eta2 <- eta * 0.5
            -(0.5 * (2 * (sinh(eta2) * 0.5 * cosh(eta2)))/(cosh(eta2)^2)^2)
          }
        )
      )
    }

    rval <- mu.eta2(rval)
  }
  rval$name <- link
  rval
}

parse.links <- function(links, default.links, ...)
{
  dots <- list(...)
  nl <- names(default.links)
  if(length(dots))
    links <- as.character(dots)
  if(is.null(names(links)))
    names(links) <- rep(nl, length.out = length(links))
  links <- as.list(links)
  for(j in nl) {
    if(is.null(links[[j]]))
      links[[j]] <- default.links[j]
  }
  links <- links[nl]
  links <- as.character(links)
  names(links) <- nl
  links
}


## http://stats.stackexchange.com/questions/41536/how-can-i-model-a-proportion-with-bugs-jags-stan
beta_bamlss <- function(...)
{
  links <- c(mu = "logit", sigma2 = "logit")

  rval <- list(
    "family" = "beta",
    "names" = c("mu", "sigma2"),
    "links" =  parse.links(links, c(mu = "logit", sigma2 = "logit"), ...),
    "valid.response" = function(x) {
      if(ok <- !all(x > 0 & x < 1)) stop("response values not in (0, 1)!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("beta", "mu"),
      "sigma2" = c("beta", "sigma2")
    ),
    "bugs" = list(
      "dist" = "dbeta",
      "eta" = BUGSeta,
      "model" = BUGSmodel,
      "reparam" = c(
        mu = "mu * (1 / sigma2)",
        sigma2 = "(1 - mu) * (1 / sigma2)"
      )
    ),
    "score" = list(
      "mu" = function(y, par, ...) {
        a <- par$mu
        b <- par$sigma2
        h1 <- a * (1 - b) / b
        h2 <- (1 - a) * (1 - b) / b
        drop(a * h2 * log(y) - a * h2 * log(1 - y) + ((1 - b) / b) * a * (1 - a) * (-digamma(h1) + digamma(h2)))
      },
      "sigma2" = function(y, par, ...) {
        a <- par$mu
        b <- par$sigma2
        h1 <- a*(1-b)/b
        h2 <- (1-a)*(1-b)/b
        drop(-(1 - b) / (b) * ( -a * digamma(h1) - (1 - a) * digamma(h2) + digamma((1 - b) / (b)) + a * log(y) + (1 - a) * log(1 - y)))
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        a <- par$mu
        b <- par$sigma2
        h1 <- a * (1 - b) / b
        h2 <- (1 - a) * (1 - b) / b
        drop(((1 - b) / b)^2 * a^2 * (1 - a)^2 * (trigamma(h1) + trigamma(h2)))
      },
      "sigma2" = function(y, par, ...) {
        a <- par$mu
        b <- par$sigma2
        h1 <- a * (1 - b) / b
        h2 <- (1 - a) * (1 - b) / b
        drop(((1 - b) / b)^2 * (a^2 * trigamma(h1) + (1 - a)^2 * trigamma(h2) - trigamma((1 - b) / (b))))
      }
    ),
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
       mu <- par$mu
       sigma2 <- par$sigma2
       a <- mu * (1 - sigma2) / (sigma2)
       b <- a * (1 - mu) / mu
       dbeta(y, shape1 = a, shape2 = b, log = log)
    },
    "p" = function(y, par, ...) {
       mu <- par$mu
       sigma2 <- par$sigma2
       a <- mu * (1 - sigma2) / (sigma2)
       b <- a * (1 - mu) / mu
       pbeta(y, shape1 = a, shape2 = b, ...)
    },
    "r" = function(n, par) {
       mu <- par$mu
       sigma2 <- par$sigma2
       a <- mu * (1 - sigma2) / (sigma2)
       b <- a * (1 - mu) / mu
       rbeta(n, shape1 = a, shape2 = b)
    },
    "q" = function(p, par) {
      mu <- par$mu
      sigma2 <- par$sigma2
      a <- mu * (1 - sigma2) / (sigma2)
      b <- a * (1 - mu) / mu
      qbeta(p, shape1 = a, shape2 = b)
    },
    "mean" = function(par) {
      a <- par$mu * (1 - par$sigma2) / (par$sigma2)
      b <- a * (1 - par$mu) / par$mu
      ex <- (1) / (1 + (b / a))
      return(ex)
    },
    "variance" = function(par) {
      a <- par$mu * (1 - par$sigma2) / (par$sigma2)
      b <- a * (1 - par$mu) / par$mu
      vx <- (a * b) / (((a + b)^2) * (a + b + 1))
      return(vx)
    }
  )
  class(rval) <- "family.bamlss"
  rval
}


betazoi_bamlss <- function(...)
{
  links <- c(mu = "logit", sigma2 = "logit", nu = "log", tau = "log")

  rval <- list(
    "family" = "betazoi",
    "names" = c("mu", "sigma2", "nu", "tau"),
    "links" =  parse.links(links, c(mu = "logit", sigma2 = "logit", nu = "log", tau = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x >= 0 & x <= 1)) stop("response values not in [0, 1]!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("betainf", "mu"),
      "sigma2" = c("betainf", "sigma2"),
      "nu" = c("betainf", "nu"),
      "tau" = c("betainf", "tau"),
      "order" = 1:4,
      "hess" = list(
        "mu" = function(x) { 1 * ((x != 1) & (x != 0)) },
        "sigma2" = function(x) { 1 * ((x != 1) & (x != 0)) }
      )
    ),
    "mu" = function(par, ...) {
       par$mu * (1 - (par$nu + par$tau) / (1 + par$nu + par$tau)) + par$tau / (1 + par$nu + par$tau)
    },
    "d" = function(y, par, log = FALSE) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      d <- ifelse(y == 0, par$nu / (1 + par$nu + par$tau), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + par$nu + par$tau))
      ifelse (y==1, par$tau / (1 + par$nu + par$tau), d)
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      h1 <- par$nu / (1 + par$nu + par$tau)
      h2 <- par$tau / (1 + par$nu + par$tau)
      cdf <- ifelse(y == 0, h1, h1 + (1 - (h1 + h2)) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
      ifelse(y == 1, 1, cdf)
    }
  )
  class(rval) <- "family.bamlss"
  rval
}


betazi_bamlss <- function(...)
{
  links <- c(mu = "logit", sigma2 = "logit", nu = "log")

  rval <- list(
    "family" = "betazi",
    "names" = c("mu", "sigma2", "nu"),
    "links" =  parse.links(links, c(mu = "logit", sigma2 = "logit", nu = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x >= 0 & x < 1)) stop("response values not in [0, 1)!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("betainf0", "mu"),
      "sigma2" = c("betainf0", "sigma2"),
      "nu" = c("betainf0", "nu"),
      "order" = 1:3,
      "hess" = list(
        "mu" = function(x) { 1 * ((x != 0)) },
        "sigma2" = function(x) { 1 * ((x != 0)) }
      )
    ),
    "mu" = function(par, ...) {
      par$mu * (1 - (par$nu) / (1 + par$nu))
    },
    "d" = function(y, par, log = FALSE) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      d <- ifelse(y == 0, par$nu / (1 + par$nu), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + par$nu))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      h1 <- par$nu / (1 + par$nu)
      ifelse(y == 0, h1, h1 + (1 - h1) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
    }
  )
  class(rval) <- "family.bamlss"
  rval
}


betaoi_bamlss <- function(...)
{
  links <- c(mu = "logit", sigma2 = "logit", tau = "log")

  rval <- list(
    "family" = "betazi",
    "names" = c("mu", "sigma2", "tau"),
    "links" =  parse.links(links, c(mu = "logit", sigma2 = "logit", tau = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0 & x <= 1)) stop("response values not in (0, 1]!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("betainf1", "mu"),
      "sigma2" = c("betainf1", "sigma2"),
      "tau" = c("betainf1", "tau"),
      "order" = 1:3,
      "hess" = list(
        "mu" = function(x) { 1 * ((x != 1)) },
        "sigma2" = function(x) { 1 * ((x != 1)) }
      )
    ),
    "mu" = function(par, ...) {
      par$mu * (1 - par$tau / (1 + par$tau)) +
        par$tau / (1 + par$tau)
    },
    "d" = function(y, par, log = FALSE) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      d <- ifelse(y == 1, par$tau / (1 + par$tau), dbeta(y, shape1 = a, shape2 = b, ncp = 0) / (1 + par$tau))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <- par$mu
      sigma <- par$sigma
      a <- mu * (1 - sigma) / (sigma)
      b <- a * (1 - mu) / mu
      h1 <- par$tau / (1 + par$tau)
      ifelse(y == 1, h1, h1 + (1 - h1) * pbeta(y, shape1 = a, shape2 = b, ncp = 0))
    }
  )
  class(rval) <- "family.bamlss"
  rval
}


process_factor_response <- function(x)
{
   if(!is.null(attr(x, "contrasts"))) {
     if(!is.null(dim(x)))
       x <- x[, ncol(x)]
   }
   if(is.factor(x)) {
     x <- as.integer(x) - 1L
   } else {
     if(!is.null(dim(x)))
         x <- x[, ncol(x)]
   }
   as.integer(x)
}


binomial_bamlss <- function(link = "logit", ...)
{
  if(link != "logit")
    return(binomial2_bamlss(...))

  rval <- list(
    "family" = "binomial",
    "names" = "pi",
    "links" = c(pi = "logit"),
    "valid.response" = function(x) {
      if(!is.factor(x)) {
        if(length(unique(x)) > 2)
          stop("response has more than 2 levels!", call. = FALSE)
      } else {
        if(nlevels(x) > 2)
          stop("more than 2 levels in factor response!", call. = FALSE)
      }
      TRUE
    },
    "bayesx" = list(
      "pi" = c("binomial_logit", "mu")
    ),
    "bugs" = list(
      "dist" = "dbern",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "mu" = function(par, ...) {
      par$pi
    },
    "d" = function(y, par, log = FALSE) {
      y <- process_factor_response(y)
      dbinom(y, size = 1, prob = par$pi, log = log)
    },
    "p" = function(y, par, ...) {
      y <- process_factor_response(y)
      pbinom(y, size = 1, prob = par$pi, ...)
    },
    "r" = function(n, par) {
      rbinom(n, size = 1, prob = par$pi)
    },
    "q" = function(p, par) {
      qbinom(p, size = 1, prob = par$pi)
    },
    "score" = list(
      "pi" = function(y, par, ...) {
        y <- process_factor_response(y)
        y - par$pi
      }
    ),
    "hess" = list(
      "pi" = function(y, par, ...) {
        par$pi * (1 - par$pi)
      }
    ),
    "initialize" = list("pi" = function(y, ...) {
      y <- process_factor_response(y)
      (y + 0.5) / 2
    }),
    "mean" = function(par) par$pi,
    "variance" = function(par) par$pi * (1 - par$pi)
  )

  class(rval) <- "family.bamlss"
  rval
}


binomial2_bamlss <- function(...)
{
  rval <- list(
    "family" = "binomial",
    "names" = "pi",
    "links" = c(pi = "probit"),
    "valid.response" = function(x) {
      if(!is.factor(x)) {
        if(length(unique(x)) > 2)
          stop("response has more than 2 levels!", call. = FALSE)
      } else {
        if(nlevels(x) > 2)
          stop("more than 2 levels in factor response!", call. = FALSE)
      }
      TRUE
    },
    "bayesx" = list(
      "pi" = c("binomial_probit", "mu")
    ),
    "mu" = function(par, ...) {
      par$pi
    },
    "d" = function(y, par, log = FALSE) {
      y <- process_factor_response(y)
      i <- y < 1
      par$pi[i] <- 1 - par$pi[i]
      if(log)
        par$pi <- log(par$pi)
      par$pi
    },
    "initialize" = list("pi" = function(y, ...) {
      y <- process_factor_response(y)
      (y + 0.5) / 2
    })
  )

  class(rval) <- "family.bamlss"
  rval
}

cloglog_bamlss <- function(...)
{
  link <- "cloglog"

  rval <- list(
    "family" = "cloglog",
    "names" = "pi",
    "links" = parse.links(link, c(pi = "cloglog"), ...),
    "valid.response" = function(x) {
      if(!is.factor(x)) stop("response must be a factor!", call. = FALSE)
      if(nlevels(x) > 2) stop("more than 2 levels in factor response!", call. = FALSE)
      TRUE
    },
    "bayesx" = list(
      "pi" = c(paste("cloglog", link, sep = "_"), "mean")
    ),
    "mu" = function(par, ...) {
      par$pi
    },
    "d" = function(y, par, log = FALSE) {
      dbinom(y, size = 1, prob = par$pi, log = log)
    },
    "p" = function(y, par, ...) {
      pbinom(y, size = 1, prob = par$pi, ...)
    },
    "r" = function(n, par) {
      pbinom(n, size = 1, prob = par$pi)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


coxph_bamlss <- function(...)
{
  rval <- list(
    "family" = "coxph",
    "names" = "gamma",
    "links" = c(gamma = "identity"),
    "transform" = function(x, ...) {
      get_ties <- function(x) {
        n <- length(x)
        ties <- vector(mode = "list", length = n)
        for(i in 1:n) {
          ties[[i]] <- which(x == x[i])
        }
        ties
      }
      y <- x$y
      attr(y[[1]], "ties") <- get_ties(y[[1]][, "time"])

      return(list("y" = y))
    },
    "d" = function(y, par, log = FALSE) {
      status <- y[, "status"]
      time <- y[, "time"]
      ties <- attr(y, "ties")
      d <- sapply(ties, length)
      n <- length(time)
      d <- rep(0, n)
      for(i in 1:n) {
        if(status[i] > 0) {
          d[i] <- sum(par$gamma[ties[[i]]]) - log(sum(exp(par$gamma[time >= time[i]]))^d[i])
        }
      }
      if(!log)
        d <- exp(d)
      d
    },
    "score" = list(
      "gamma" = function(y, par, ...) {
        status <- y[, "status"]
        time <- y[, "time"]
        ties <- attr(y, "ties")
        d <- sapply(ties, length)
        n <- length(time)
        risk <- rep(0, n)
        for(i in 1:n) {
          C <- which(time >= time[i])
          risk[i] <- sum(d[i] / sum(exp(par$gamma[C])))
        }
        score <- status - exp(par$gamma) * risk
        score
      }
    ),
    "hess" = list(
      "gamma" = function(y, par, ...) {
        status <- y[, "status"]
        time <- y[, "time"]
        ties <- attr(y, "ties")
        d <- sapply(ties, length)
        n <- length(time)
        risk <- risk2 <- rep(0, n)
        for(i in 1:n) {
          C <- which(time >= time[i])
          risk[i] <- sum(d[i] / sum(exp(par$gamma[C])))
          risk2[i] <- sum(d[i] / (sum(exp(par$gamma[C]))^2))
        }
        hess <- - exp(par$gamma) * risk + exp(2 * par$gamma) * risk2
        -hess
      }
    ),
    "initialize" = list(
      "gamma" = function(y, ...) { rep(0, nrow(y)) }
    )
  )

  class(rval) <- "family.bamlss"
  rval
}


if(FALSE) {
  library(survival)
  col1 <- colon[colon$etype==1, ]
  col1$differ <- as.factor(col1$differ)
  col1$sex <- as.factor(col1$sex)

  b1 <- bamlss(Surv(time, status) ~ s(nodes), family = "coxph", data = col1, sampler = FALSE)

  b2 <- gam(time ~ perfor + rx + obstruct + adhere + sex + s(age,by=sex) + s(nodes), family = cox.ph, data = col1)
}


gaussian_bamlss <- function(...)
{
  links <- c(mu = "identity", sigma = "log")

  rval <- list(
    "family" = "gaussian",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "bayesx" = list(
      "mu" = switch(links["mu"],
        "identity" = c("normal2", "mu"),
        "inverse" = c("normal_mu_inv", "mean")
      ),
      "sigma" = switch(links["sigma"],
        "log" = c("normal2", "sigma"),
        "logit" = c("normal_sigma_logit", "scale")
      )
    ),
    "bugs" = list(
      "dist" = "dnorm",
      "eta" = BUGSeta,
      "model" = BUGSmodel,
      "reparam" = c(sigma = "1 / sqrt(sigma)")
    ),
    "keras" = list(
      "nloglik" = function(y_true, y_pred) {
        K = keras::backend()

        mu = y_pred[, 1]
        sigma = K$exp(y_pred[,2])
        sigma2 = K$pow(sigma, 2)

        ll = -0.5 * K$log(6.28318530717959 * sigma2) - 0.5 * K$pow((y_true[,1] - mu), 2) / sigma2
        ll = K$sum(ll)

        return(-1 * ll)
      }
    ),
    "score" = list(
      "mu" = function(y, par, ...) { drop((y - par$mu) / (par$sigma^2)) },
      "sigma" = function(y, par, ...) { drop(-1 + (y - par$mu)^2 / (par$sigma^2)) }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { drop(1 / (par$sigma^2)) },
      "sigma" = function(y, par, ...) { rep(2, length(y)) }
    ),
    "loglik" = function(y, par, ...) {
      sum(dnorm(y, par$mu, par$sigma, log = TRUE))
    },
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "p" = function(y, par, ...) {
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    },
    "r" = function(n, par) {
      rnorm(n, mean = par$mu, sd = par$sigma)
    },
    "q" = function(p, par) {
      qnorm(p, mean = par$mu, sd = par$sigma)
    },
    "crps" = function(y, par, ...) {
      scoringRules::crps_norm(y, mean = par$mu, sd = par$sigma)
    },
    "initialize" = list(
      "mu"    = function(y, ...) { (y + mean(y)) / 2 },
      "sigma" = function(y, ...) { rep(sd(y), length(y)) }
    ),
    "mean"      = function(par) par$mu,
    "variance"  = function(par) par$sigma^2,
    "valid.response" = function(x) {
      if(is.factor(x) | is.character(x))
        stop("the response should be numeric!")
      return(TRUE)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


Gaussian_bamlss <- function(...)
{
  rval <- list(
    "family" = "gaussian",
    "names" = "mu",
    "links" = c(mu = "identity"),
    "bayesx" = list(
      "mu" = c("gaussian", "mu")
    ),
    "score" = list(
      "mu" = function(y, par, ...) { y - par$mu }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { rep(1, length(par$mu)) }
    ),
    "loglik" = function(y, par, ...) {
      sum(dnorm(y, par$mu, 1, log = TRUE))
    },
    "initialize" = list(
      "mu" = function(y, ...) { (y + mean(y)) / 2 }
    ),
    "valid.response" = function(x) {
      if(is.factor(x) | is.character(x))
        stop("the response should be numeric!")
      return(TRUE)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


AR0_bamlss <- function(ar.start = NULL, ...)
{
  rval <- list(
    "family" = "AR1",
    "names" = c("mu", "alpha", "sigma"),
    "links" = c(mu = "identity", alpha = "rhogit", sigma = "log"),
    "initialize" = list(
      "mu" = function(y, ...) { (y + mean(y)) / 2 },
      "alpha" = function(y, ...) {
        rep(acf(y - mean(y), lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2], length(y))
      },
      "sigma" = function(y, ...) { rep(sd(y), length(y)) }
    )
  )
  if(is.null(ar.start)) {
    rval$d <- function(y, par, log = FALSE, ...) {
      e <- y - par$mu
      e1 <- c(0, e[-length(e)])
      mu <- par$mu + par$alpha * e1
      dnorm(y, mu, par$sigma, log = log)
    }
    rval$p <- function(y, par, ...) {
      e <- y - par$mu
      e1 <- c(0, e[-length(e)])
      mu <- par$mu + par$alpha * e1
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    }
#    rval$score <- list(
#      "mu" = function(y, par, ...) {
#        e <- y - par$mu
#        e1 <- c(0, e[-length(e)])
#        mu <- par$mu + par$alpha * e1
#        drop((y - mu) / (par$sigma^2))
#      },
#      "sigma" = function(y, par, ...) {
#        e <- y - par$mu
#        e1 <- c(0, e[-length(e)])
#        mu <- par$mu + par$alpha * e1
#        drop(-1 + (y - mu)^2 / (par$sigma^2))
#      }
#    )
#    rval$hess <- list(
#      "mu" = function(y, par, ...) { drop(1 / (par$sigma^2)) },
#      "sigma" = function(y, par, ...) { rep(2, length(y)) }
#    )
  } else {
    nobs <- length(ar.start)
    i <- which(ar.start)
    rval$d <- function(y, par, log = FALSE, ...) {
      e <- y - par$mu
      e <- c(0, e[-nobs])
      e[i] <- 0
      mu <- par$mu + par$alpha * e
      dnorm(y, mu, par$sigma, log = log)
    }
    rval$p <- function(y, par, ...) {
      e <- y - par$mu
      e <- c(0, e[-nobs])
      e[i] <- 0
      mu <- par$mu + par$alpha * e
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    }
#    rval$score <- list(
#      "mu" = function(y, par, ...) {
#        e <- e0 <- y - par$mu
#        e <- c(0, e[-nobs])
#        e[i] <- 0
#        mu <- par$mu + par$alpha * e
#        drop((y - mu) / (par$sigma^2))
#      },
#      "sigma" = function(y, par, ...) {
#        e <- y - par$mu
#        e <- c(0, e[-nobs])
#        e[i] <- 0
#        mu <- par$mu + par$alpha * e
#        drop(-1 + (y - mu)^2 / (par$sigma^2))
#      }
#    )
#    rval$hess <- list(
#      "mu" = function(y, par, ...) { drop(1 / (par$sigma^2)) },
#      "sigma" = function(y, par, ...) { rep(2, length(y)) }
#    )
  }

  rval$r = function(n, par) {
    rnorm(n, mean = par$mu, sd = par$sigma)
  }
  rval$q = function(p, par) {
    qnorm(p, mean = par$mu, sd = par$sigma)
  }

  class(rval) <- "family.bamlss"
  rval
}


AR1_bamlss <- function(...)
{
  ## https://arxiv.org/pdf/1711.05204.pdf
  rval <- list(
    "family" = "AR1",
    "names" = c("mu", "sigma", "rho"),
    "links" = c(mu = "identity", sigma = "log", rho = "sigmoid"),
    "d" = function(y, par, log = FALSE, ...) {
      if(!is.null(dim(y))) {
        i <- which(y[, 2L] > 0)
        y <- y[, 1L]
      } else {
        i <- 1L
      }
      e <- y - par$mu
      e <- c(0, e[-length(e)])
      e[i] <- 0
      mu <- par$mu + par$rho * e
      d <- dnorm(y, mu, par$sigma, log = log)
      if(!log)
        d <- exp(d)
      return(d)
    },
    "p" = function(y, par, ...) {
      if(!is.null(dim(y))) {
        i <- which(y[, 2L] > 0)
        y <- y[, 1L]
      } else {
        i <- 1L
      }
      e <- y - par$mu
      e <- c(0, e[-length(e)])
      e[i] <- 0
      mu <- par$mu + par$rho * e
      d <- pnorm(y, mu, par$sigma)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        if(!is.null(dim(y))) {
          i <- which(y[, 2L] > 0)
          y <- y[, 1L]
        } else {
          i <- 1L
        }
        e <- y - par$mu
        e <- c(0, e[-length(e)])
        e[i] <- 0
        mu <- par$mu + par$rho * e
        (y - mu) / par$sigma^2
      },
      "sigma" = function(y, par, ...) {
        if(!is.null(dim(y))) {
          i <- which(y[, 2L] > 0)
          y <- y[, 1L]
        } else {
          i <- 1L
        }
        e <- y - par$mu
        e <- c(0, e[-length(e)])
        e[i] <- 0
        mu <- par$mu + par$rho * e
        (y - mu)^2 / par$sigma^2 - 1
      },
      "rho" = function(y, par, ...) {
        if(!is.null(dim(y))) {
          i <- which(y[, 2L] > 0)
          y <- y[, 1L]
        } else {
          i <- 1L
        }
        e <- y - par$mu
        e <- c(0, e[-length(e)])
        e[i] <- 0
        mu <- par$mu + par$rho * e
        eta_rho <- tanh(par$rho / 2)
        (y - mu) / par$sigma^2 * 2 * exp(-eta_rho) * e / (exp(-eta_rho) + 1)^2
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        1 / par$sigma^2
      },
      "sigma" = function(y, par, ...) {
        rep(2, length(par[[1L]]))
      },
      "rho" = function(y, par, ...) {
        if(!is.null(dim(y))) {
          i <- which(y[, 2L] > 0)
          y <- y[, 1L]
        } else {
          i <- 1L
        }
        e <- y - par$mu
        e <- c(0, e[-length(e)])
        e[i] <- 0
        eta_rho <- tanh(par$rho / 2)
        4 * e^2 * exp(2 * (eta_rho - log(par$sigma))) / (exp(eta_rho) + 1)^4
      }
    ),
    "crps" = function(y, par, ...) {
      if(!is.null(dim(y))) {
        i <- which(y[, 2L] > 0)
        y <- y[, 1L]
      } else {
        i <- 1L
      }
      e <- y - par$mu
      e <- c(0, e[-length(e)])
      e[i] <- 0
      mu <- par$mu + par$rho * e
      scoringRules::crps_norm(y, mean = mu, sd = par$sigma)
    },
    "initialize" = list(
      "mu"    = function(y, ...) {
        start <- if(is.null(dim(y))) {
          (y + mean(y, na.rm = TRUE)) / 2
        } else {
          (y[, 1L] + mean(y[, 1L], na.rm = TRUE)) / 2
        }
        start
      },
      "sigma" = function(y, ...) {
        start <- if(is.null(dim(y))) {
          rep(sd(y), length(y))
        } else {
          rep(sd(y[, 1L], na.rm = TRUE), nrow(y))
        }
        start
      },
      "rho" = function(y, ...) {
        rep(1e-05, if(is.null(dim(y))) length(y) else nrow(y))
      }
    )
  )

  class(rval) <- "family.bamlss"
  rval
}


if(FALSE) {
  time <- seq(-3, 3, length = 5000)
  y <- rep(0, length(time))
  rho <- bamlss:::make.link2("sigmoid")$linkinv

  for(i in 2:length(time)) {
    y[i] <- sin(time[i]) + rho(0.6 * cos(time[i] * 4)) * (y[i - 1] - sin(time[i - 1])) + rnorm(1, sd = 0.2)
  }
  time <- time[-1]
  y <- y[-1]

  plot(time, y)

  f <- list(
    Y ~ s(time),
    sigma ~ s(time),
    rho ~ s(lag1time)
  )

  d <- data.frame("time" = time[-1], "lag1time" = time[-length(time)])
  d$Y <- cbind("y" = y[-1], "ar.start" = c(1, rep(0, nrow(d) - 1)))

  b <- bamlss(f, family = AR1_bamlss, data = d)

  p <- predict(b, model = "rho")
  plot(d$lag1time, p, type = "l")
  lines(0.6 * cos(d$lag1time * 4) ~ d$lag1time, col = 2, lwd = 2)
}


AR00_bamlss <- function(ar.start = NULL, ...)
{
  if(is.null(list(...)$mu)) {
    rval <- list(
      "family" = "AR1",
      "names" = c("alpha", "sigma"),
      "links" = c(alpha = "rhogit", sigma = "log"),
      "initialize" = list(
        "alpha" = function(y, ...) {
          rep(acf(y, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2], length(y))
        },
        "sigma" = function(y, ...) { rep(sd(y), length(y)) }
      )
    )
    if(is.null(ar.start)) {
      rval$d <- function(y, par, log = FALSE, ...) {
        y1 <- c(0, y[-length(y)])
        mu <- par$alpha * y1
        dnorm(y, mu, par$sigma, log = log)
      }
      rval$p <- function(y, par, ...) {
        y1 <- c(0, y[-length(y)])
        mu <- par$alpha * y1
        pnorm(y, mean = par$mu, sd = par$sigma, ...)
      }
    } else {
      nobs <- length(ar.start)
      i <- which(ar.start)
      rval$d <- function(y, par, log = FALSE, ...) {
        y1 <- c(0, y[-nobs])
        y1[i] <- 0
        mu <- par$alpha * y1
        dnorm(y, mu, par$sigma, log = log)
      }
      rval$p <- function(y, par, ...) {
        y1 <- c(0, y[-nobs])
        y1[i] <- 0
        mu <- par$alpha * y1
        pnorm(y, mean = par$mu, sd = par$sigma, ...)
      }
    }
  } else {
    rval <- list(
      "family" = "AR1",
      "names" = c("mu", "sigma", "alpha"),
      "links" = c(mu = "identity", sigma = "log", alpha = "rhogit"),
      "initialize" = list(
        "mu" = function(y, ...) { (y + mean(y)) / 2 },
        "sigma" = function(y, ...) { rep(sd(y), length(y)) },
        "alpha" = function(y, ...) {
          rep(0, length(y))
        }
      )
    )
    if(is.null(ar.start)) {
      rval$d <- function(y, par, log = FALSE, ...) {
        e <- y - par$mu
        e1 <- c(0, e[-length(e)])
        mu <- par$mu + par$alpha * e1
        dnorm(y, mu, par$sigma, log = log)
      }
      rval$p <- function(y, par, ...) {
        e <- y - par$mu
        e1 <- c(0, e[-length(e)])
        mu <- par$mu + par$alpha * e1
        pnorm(y, mean = mu, sd = par$sigma, ...)
      }
    } else {
      nobs <- length(ar.start)
      i <- which(ar.start)
      rval$d <- function(y, par, log = FALSE, ...) {
        e <- y - par$mu
        e <- c(0, e[-nobs])
        e[i] <- 0
        mu <- par$mu + par$alpha * e
        dnorm(y, mu, par$sigma, log = log)
      }
      rval$p <- function(y, par, ...) {
        e <- y - par$mu
        e <- c(0, e[-nobs])
        e[i] <- 0
        mu <- par$mu + par$alpha * e
        pnorm(y, mean = par$mu, sd = par$sigma, ...)
      }
    }
  }

  class(rval) <- "family.bamlss"
  rval
}

beta1_bamlss <- function(ar.start, ...)
{
  rval <- list(
    "family" = "beta1",
    "names" = c("mu", "phi", "rho"),
    "links" = c(mu = "logit", phi = "log", rho = "rhogit"),
    "valid.response" = function(x) {
      if(ok <- !all(x > 0 & x < 1)) stop("response values not in (0, 1)!", call. = FALSE)
      ok
    }
  )

  lmu <- make.link2(rval$links["mu"])
  lphi <- make.link2(rval$links["phi"])
  lrho <- make.link2(rval$links["rho"])

  nobs <- length(ar.start)
  i <- which(ar.start)

  rval$d <- function(y, par, log = FALSE, ...) {
    mu <- lmu$linkfun(par$mu)
    e <- lmu$linkfun(y) - mu
    e <- c(0, e[-nobs])
    e[i] <- 0
    mu <- lmu$linkinv(mu + par$rho * e)
    a <- mu * par$phi
    b <- (1 - mu) * par$phi
    return(dbeta(y, a, b, log = log))
  }

  rval$score <- list(
    "mu" = function(y, par, ...) {
      mu <- lmu$linkfun(par$mu)
      e <- lmu$linkfun(y) - mu
      e <- c(0, e[-nobs])
      e[i] <- 0
      eta <- mu + par$rho * e
      mu <- lmu$linkinv(eta)
      ystar <- qlogis(y)
      mustar <- digamma(mu * par$phi) - digamma((1 - mu) * par$phi)
      return(par$phi * (ystar - mustar) * lmu$mu.eta(eta))
    },
    "phi" = function(y, par, ...) {
      mu <- lmu$linkfun(par$mu)
      e <- lmu$linkfun(y) - mu
      e <- c(0, e[-nobs])
      e[i] <- 0
      eta <- mu + par$rho * e
      mu <- lmu$linkinv(eta)
      ystar <- qlogis(y)
      mustar <- digamma(mu * par$phi) - digamma((1 - mu) * par$phi)
      (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*par$phi) + digamma(par$phi)) * lphi$mu.eta(lphi$linkfun(par$phi))
    },
    "rho" = function(y, par, ...) {
      mu <- lmu$linkfun(par$mu)
      e <- lmu$linkfun(y) - mu
      e <- c(0, e[-nobs])
      e[i] <- 0
      eta <- mu + par$rho * e
      mu <- lmu$linkinv(eta)
      ystar <- qlogis(y)
      mustar <- digamma(mu * par$phi) - digamma((1 - mu) * par$phi)
      sa <- par$phi * (ystar - mustar) * lmu$mu.eta(eta)
      return(sa * e * lrho$mu.eta(lrho$linkfun(par$rho)))
    }
  )

  rval$p <- function(y, par, ...) {
    mu <- lmu$linkfun(par$mu)
    e <- lmu$linkfun(y) - mu
    e <- c(0, e[-nobs])
    e[i] <- 0
    mu <- lmu$linkinv(mu + par$rho * e)
    a <- mu * par$phi
    b <- (1 - mu) * par$phi
    pbeta(y, shape1 = a, shape2 = b, ...)
  }

  class(rval) <- "family.bamlss"
  rval
}


gpareto_bamlss <- function(...)
{
  rval <- list(
    "family" = "gpareto",
    "names" = c("xi", "sigma"),
    "links" = c(xi = "log", sigma = "log"),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x >= 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "d" = function(y, par, log = FALSE) {
      d <- -log(par$sigma) - (1 / par$xi + 1) * log(1 + par$xi * y / par$sigma)
      if(!log)
        d <- exp(d)
      d
    },
    "mu" = function(par, ...) {
      par$sigma / (1 - par$xi)
    },
    "p" = function(y, par, ...) {
      ##p <- 1 - (1 + par$xi * y / par$sigma)^(-1 / par$xi)
      y <- pmax(y, 0) / par$sigma
      p <- pmax(1 + par$xi * y, 0)
      p <- 1 - p^(-1 / par$xi)
      p
    },
    "q" = function(p, par, ...) {
      ##(par$sigma / par$xi) * ((1 - p)^(-par$xi) - 1)
      par$sigma * ((1 - p)^(-par$xi) - 1) / par$xi
    },
    "r" = function(n, par) {
      nobs <- length(par[[1]])
      samps <- matrix(NA, nrow = nobs, ncol = n)
      for(i in 1:n)
        samps[, i] <- par$sigma * ((1 - runif(nobs))^(-par$xi) - 1) / par$xi
      if(n < 2)
        samps <- as.vector(samps)
      samps
    },
    "mean" = function(par) {
#      if (par$xi >= 1) {
#        return(NA)
#      } else {
#        return(par$sigma / (1 - par$xi))
#      }
      ## TODO: NA vs Inf ???
      ifelse(par$xi >= 1, NA, par$sigma / (1 - par$xi))
    },
    "variance" = function(par) {
#      if (par$xi >= 0.5) {
#        return(NA)
#      } else {
#        return((par$sigma^2) / (1 - par$xi)^2 * (1 - 2 * par$xi))
#      }
      ifelse(par$xi >= 0.5, NA, (par$sigma^2) / (1 - par$xi)^2 * (1 - 2 * par$xi))
    }
  )

  rval$score <- list(
    "xi" = function(y, par, ...) {
      .Call("gpareto_score_xi", as.numeric(y), as.numeric(par$xi),
        as.numeric(par$sigma), PACKAGE = "bamlss")
    },
    "sigma" = function(y, par, ...) {
      .Call("gpareto_score_sigma", as.numeric(y), as.numeric(par$xi),
        as.numeric(par$sigma), PACKAGE = "bamlss")
    }
  )

  rval$hess <- list(
    "xi" = function(y, par, ...) {
      .Call("gpareto_hess_xi", as.numeric(y), as.numeric(par$xi),
        as.numeric(par$sigma), PACKAGE = "bamlss")
    },
    "sigma" = function(y, par, ...) {
      .Call("gpareto_hess_sigma", as.numeric(y), as.numeric(par$xi),
        as.numeric(par$sigma), PACKAGE = "bamlss")
    }
  )

#  rval$score <- list(
#    "xi" = function(y, par, ...) {
#      ys <- y / par$sigma
#      xi1 <- 1 / par$xi
#      xi1ys <- 1 + par$xi * ys
#      -((xi1 + 1) * (par$xi * ys/xi1ys) - xi1 * log(xi1ys))
#    },
#    "sigma" = function(y, par, ...) {
#      ys <- y / par$sigma
#      -(1 - (1/par$xi + 1) * (par$xi * ys /(1 + par$xi * ys)))
#    }
#  )

#  rval$hess <- list(
#    "xi" = function(y, par, ...) {
#      ys <- y / par$sigma
#      xi1 <- 1 / par$xi
#      xi2 <- par$xi^2
#      xiys <- par$xi * ys
#      xi1ys <- 1 + xiys
#      xiysxi1ys <- xiys / xi1ys
#      xi1xiysxi1ys <- xi1 * xiysxi1ys
#      (xi1 + 1) * (xiysxi1ys - xiys^2/xi1ys^2) - xi1xiysxi1ys - ((xi1 - par$xi * (2 * xi2)/xi2^2) * log(xi1ys) + xi1xiysxi1ys)
#    },
#    "sigma" = function(y, par, ...) {
#      s1 <- 1 / par$sigma
#      ys <- y / par$sigma
#      xiys1 <- par$xi * y * s1
#      xiys <- par$xi * ys
#      xi1ys <- 1 + xiys
#      -1 * (1/par$xi + 1) * ((xiys1 - par$xi * y * par$sigma * 2 * s1^2)/xi1ys + xiys1 * xiys1/xi1ys^2)
#    }
#  )

  rval$initialize <- list(
    "xi" = function(y, ...) { rep(mean(y) + 0.5, length(y)) },
    "sigma" = function(y, ...) { rep(sd(y), length(y)) }
  )

  class(rval) <- "family.bamlss"
  rval
}


dpareto_bamlss <- function(...)
{
  rval <- list(
    "family" = "dpareto",
    "names" = c("alpha", "sigma"),
    "links" = c(alpha = "log", sigma = "log"),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x >= 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "d" = function(y, par, log = FALSE) {
      d <- log((1 / (1 + (y - 1) / par$sigma))^par$alpha - (1 / (1 + y/par$sigma))^par$alpha)
      if(!log)
        d <- exp(d)
      d
    }
  )

  rval$initialize <- list(
    "alpha" = function(y, ...) { rep(mean(y) + 0.5, length(y)) },
    "sigma" = function(y, ...) { rep(sd(y), length(y)) }
  )

  class(rval) <- "family.bamlss"
  rval
}


gaussian1_bamlss <- function(...)
{
  links <- c(mu = "identity", sigma = "log")

  rval <- list(
    "family" = "gaussian",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "bayesx" = list(
      "mu" = switch(links["mu"],
        "identity" = c("normal2", "mu"),
        "inverse" = c("normal_mu_inv", "mean")
      ),
      "sigma" = switch(links["sigma"],
        "log" = c("normal2", "sigma"),
        "logit" = c("normal_sigma_logit", "scale")
      )
    ),
    "bugs" = list(
      "dist" = "dnorm",
      "eta" = BUGSeta,
      "model" = BUGSmodel,
      "reparam" = c(sigma = "1 / sqrt(sigma)")
    ),
    "score" = list(
      "mu" = function(y, par, ...) { drop((y - par$mu) / (par$sigma^2)) },
      "sigma" = function(y, par, ...) { drop(-1 + (y - par$mu)^2 / (par$sigma^2)) }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { drop(1 / (par$sigma^2)) },
      "sigma" = function(y, par, ...) { rep(2, length(y)) }
    ),
    "gradient" = list(
      "mu" = function(g, y, par, x, ...) {
        par$mu <- par$mu + x$get.mu(x$X, g)
        crossprod(if(!is.null(x$xbin.ind)) x$X[x$xbin.ind, , drop = FALSE] else x$X, drop((y - par$mu) / (exp(par$sigma)^2)))
      },
      "sigma" = function(g, y, par, x, ...) {
        par$sigma <- par$sigma + x$get.mu(x$X, g)
        crossprod(if(!is.null(x$xbin.ind)) x$X[x$xbin.ind, , drop = FALSE] else x$X, drop(-1 + (y - par$mu)^2 / (exp(par$sigma)^2)))
      }
    ),
    "hessian" = list(
      "mu" = function(g, y, par, x, ...) {
        if(!is.null(x$xbin.ind)) {
          return(crossprod(x$X[x$xbin.ind, , drop = FALSE], x$X[x$xbin.ind, , drop = FALSE] * drop(1 / exp(par$sigma)^2)))
        } else {
          return(crossprod(x$X, x$X * drop(1 / exp(par$sigma)^2)))
        }
      },
      "sigma" = function(g, y, par, x, ...) {
        if(is.null(x$XX)) return(crossprod(x$X)) else return(x$XX)
      }
    ),
    "loglik" = function(y, par, ...) {
      sum(dnorm(y, par$mu, par$sigma, log = TRUE))
    },
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "p" = function(y, par, ...) {
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    },
    "initialize" = list(
      "mu" = function(y, ...) { (y + mean(y)) / 2 },
      "sigma" = function(y, ...) { rep(sd(y), length(y)) }
    )
  )

  class(rval) <- "family.bamlss"
  rval
}


gaussian2_bamlss <- function(...)
{
  links <- c(mu = "identity", sigma2 = "log")

  rval <- list(
    "family" = "gaussian2",
    "names" = c("mu", "sigma2"),
    "links" = parse.links(links, c(mu = "identity", sigma2 = "log"), ...),
    "bayesx" = list(
      "mu" = switch(links["mu"],
        "identity" = c("normal", "mu"),
        "inverse" = c("normal_mu_inv", "mean")
      ),
      "sigma2" = switch(links["sigma2"],
        "log" = c("normal", "sigma2"),
        "logit" = c("normal_sigma2_logit", "scale")
      )
    ),
    "bugs" = list(
      "dist" = "dnorm",
      "eta" = BUGSeta,
      "model" = BUGSmodel,
      "reparam" = c(sigma2 = "1 / sigma2")
    ),
    "score" = list(
      "mu" = function(y, par, ...) { drop((y - par$mu) / par$sigma2) },
      "sigma2" = function(y, par, ...) { drop(-0.5 + (y - par$mu)^2 / (2 * par$sigma2)) }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { drop(1 / par$sigma2) },
      "sigma2" = function(y, par, ...) { rep(0.5, length(y)) }
    ),
    "loglik" = function(y, par, ...) {
      sum(dnorm(y, par$mu, sqrt(par$sigma2), log = TRUE))
    },
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = sqrt(par$sigma2), log = log)
    },
    "p" = function(y, par, ...) {
      pnorm(y, mean = par$mu, sd = sqrt(par$sigma2), ...)
    },
    "q" = function(p, par) {
      qnorm(p, mean = par$mu, sd = sqrt(par$sigma2))
    },
    "mean"     = function(par) par$mu,
    "variance" = function(par) par$sigma2
  )

  class(rval) <- "family.bamlss"
  rval
}


truncgaussian2_bamlss <- function(...)
{
  links <- c(mu = "identity", sigma2 = "log")

  rval <- list(
    "family" = "truncgaussian2",
    "names" = c("mu", "sigma2"),
    "links" = parse.links(links, c(mu = "identity", sigma2 = "log"), ...),
    "bayesx" = list(
      "mu" =  c("truncnormal", "mu"),
      "sigma2" = c("truncnormal", "sigma2")
    ),
    "mu" = function(par, ...) {
      mu <-  par$mu
      sigma <- sqrt(par$sigma2)
      arg <- - mu / sigma
      mu + sigma * dnorm(arg) / (1 - pnorm(arg))
    },
    "d" = function(y, par, log = FALSE) {
      sigma <- sqrt(par$sigma2)
      arg <- - par$mu / sigma
      d <- dnorm(y / sigma + arg) / (1 - pnorm(arg))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      sigma <- sqrt(par$sigma2)
      arg <- - par$mu / sigma
      2 * (pnorm(y / sigma + arg) - pnorm(arg))
    }
  )

  class(rval) <- "family.bamlss"
  rval
}

truncgaussian_bamlss <- function(...)
{
  links <- c(mu = "identity", sigma = "log")

  rval <- list(
    "family" = "truncgaussian",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "bayesx" = list(
      "mu" =  c("truncnormal2", "mu"),
      "sigma" = c("truncnormal2", "sigma")
    ),
    "mu" = function(par, ...) {
      mu <-  par$mu
      sigma <- par$sigma
      arg <- - mu / sigma
      mu + sigma * dnorm(arg) / (1 - pnorm(arg))
    },
    "d" = function(y, par, log = FALSE) {
      mu <-  par$mu
      sigma <- par$sigma
      arg <- - mu / sigma
      d <- dnorm(y / sigma + arg) / (1 - pnorm(arg))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <-  par$mu
      sigma <- par$sigma
      arg <- - mu / sigma
      2 * (pnorm(y / sigma + arg) - pnorm(arg))
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        rval <- with(par, (y - mu) / (sigma^2) - (1 / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)))
        return(drop(rval))
      },
      "sigma" = function(y, par, ...) {
        rval <- with(par, -1 + (y - mu)^2 / (sigma^2) + (mu / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)))
        return(drop(rval))
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        rval <- with(par, 1 / (sigma^2) - (mu / sigma^2) * (1 / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma))
          - ((1 / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)))^2)
        return(drop(rval))
      },
      "sigma" = function(y, par, ...) {
        rval <- with(par, 2 - (mu / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma)) *
          (1 + (mu / sigma)^2 + (mu / sigma)*(dnorm(mu / sigma) / pnorm(mu / sigma))))
        return(drop(rval))
      }
    ),
    "loglik" = function(y, par, ...) {
      rval <- with(par, sum(-0.5 * log(2 * pi) - log(sigma) - (y - mu)^2 / (2*sigma^2) - log(pnorm(mu / sigma))))
      return(rval)
    },
    "valid.response" = function(x) {
      if(is.factor(x) | is.character(x))
        stop("the response should be numeric!")
      return(TRUE)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


#trunc_bamlss <- function(direction = "left", point = 0, ...)
#{
#  links <- c(mu = "identity", sigma = "log")

#  tgrad <- function(y, par, what = "mu") {
#    par$sigma[par$sigma < 1] <- 1
#    par$mu[par$mu < -10] <- -10
#    par$mu[par$mu > 10] <- 10
#    resid <- y - par$mu
#    if(direction == "left") {
#      trunc <- par$mu - point
#      sgn <- 1
#    } else {
#      trunc <- point - par$mu
#      sgn <- -1
#    }
#    mills <- dnorm(trunc / par$sigma) / pnorm(trunc / par$sigma)
#    g <- if(what == "mu") {
#      resid / par$sigma^2 - sgn / par$sigma * mills
#    } else {
#      resid^2 / par$sigma^3 - 1 / par$sigma + trunc / par$sigma^2 * mills
#    }
#    return(g)
#  }

#  rval <- list(
#    "family" = "truncreg",
#    "names" = c("mu", "sigma"),
#    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
#    "d" = function(y, par, log = FALSE) {
#      par$sigma[par$sigma < 1] <- 1
#      par$mu[par$mu < -10] <- -10
#      par$mu[par$mu > 10] <- 10
#      resid <- y - par$mu
#      if(direction == "left") {
#        trunc <- par$mu - point
#      } else {
#        trunc <- point - par$mu
#      }
#      ll <- log(dnorm(resid / par$sigma)) - log(par$sigma) - log(pnorm(trunc / par$sigma))
#      if(!log) ll <- exp(ll)
#      ll
#    },
#    "p" = function(y, par, ...) {
#      ## Needs msm package!
#      ptnorm(y, mean = par$mu, sd = par$sigma,
#        lower = if(direction == "left") point else -Inf,
#        upper = if(direction == "right") point else Inf)
#    },
#    "score" = list(
#      "mu" = function(y, par, ...) {
#        as.numeric(tgrad(y, par, what = "mu"))
#      },
#      "sigma" = function(y, par, ...) {
#        as.numeric(tgrad(y, par, what = "sigma"))
#      }
#    )
#  )
#
#  class(rval) <- "family.bamlss"
#  rval
#}


cnorm_bamlss <- function(...)
{
  f <- list(
    "family" = "cnorm",
    "names" = c("mu", "sigma"),
    "links" = c("mu" = "identity", "sigma" = "log")
  )
  f$transform = function(x, ...) {
    attr(x$y[[1]], "check") <- as.integer(x$y[[1]] <= 0)
    list("y" = x$y)
  }
  f$bayesx = list(
    "mu" =  c("cnormal", "mu"),
    "sigma" = c("cnormal", "sigma")
  )
  f$score <- list(
    "mu" = function(y, par, ...) {
      .Call("cnorm_score_mu",
        as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")), PACKAGE = "bamlss")
    },
    "sigma" = function(y, par, ...) {
      .Call("cnorm_score_sigma",
        as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")), PACKAGE = "bamlss")
    }
  )
  f$hess <- list(
    "mu" = function(y, par, ...) {
      .Call("cnorm_hess_mu",
        as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")), PACKAGE = "bamlss")
    },
    "sigma" = function(y, par, ...) {
      .Call("cnorm_hess_sigma",
        as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")), PACKAGE = "bamlss")
    }
  )
  f$loglik <- function(y, par, ...) {
    .Call("cnorm_loglik",
      as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma),
      as.integer(attr(y, "check")), PACKAGE = "bamlss")
  }
  f$d <- function(y, par, log = FALSE) {
    ifelse(y <= 0, pnorm(-par$mu / par$sigma, log.p = log),
      dnorm((y - par$mu) / par$sigma, log = log) / par$sigma^(1 - log) - log(par$sigma) * log)
  }
  f$p <- function(y, par, log = FALSE) {
    return(ifelse(y == 0, runif(length(y), 0, pnorm(0, par$mu, par$sigma)), pnorm(y - par$mu, 0, par$sigma)))
  }
  f$q <- function(p, par, ...) {
    rval <- qnorm(p) * par$sigma + par$mu
    pmax(pmin(rval, Inf), 0)
  }
  f$r <- function(n, par, ...) {
    rval <- rnorm(n) * par$sigma + par$mu
    pmax(pmin(rval, Inf), 0)
  }
  f$initialize = list(
    "mu" = function(y, ...) { (y + mean(y)) / 2 },
    "sigma" = function(y, ...) { rep(sd(y), length(y)) }
  )
  f$mean <- function(par) {
    mu <- par$mu
    sigma <- par$sigma
    ex <- pnorm(mu / sigma) * (mu + sigma * ((dnorm(mu / sigma) / pnorm(mu / sigma))))
    return(ex)
  }
  f$variance <- function(par) {
    mu <- par$mu
    sigma <- par$sigma
    X <- pnorm(-mu / sigma)
    Y <- dnorm(mu / sigma) / (1 - X)
    Z <- Y^2 - Y * (-mu / sigma)
    vx <- (sigma^2) * (1 - X) * ((1 - Z) + (((-mu/sigma) - Y)^2)* X)
    return(vx)
  }
  class(f) <- "family.bamlss"
  f
}


pcnorm_bamlss <- function(start = 2, update = FALSE, ...)
{
  ## ll_yp <- "-1/2*(log(2*pi) + 2*log(sigma) + (y^(1/exp(lambda)) - mu)^2/sigma^2"
  ## derivs_yp <- compute_derivatives(ll, "lambda")
  f <- list(
    "family" = "pcnorm",
    "names" = c("mu", "sigma", "lambda"),
    "links" = c(mu = "identity", sigma = "log", lambda = "log")
  )
  f$score <- list(
    "mu" = function(y, par, ...) {
      .Call("cnorm_score_mu",
        as.numeric(y^(1 / par$lambda)), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")), PACKAGE = "bamlss")
    },
    "sigma" = function(y, par, ...) {
      .Call("cnorm_score_sigma",
        as.numeric(y^(1 / par$lambda)), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")), PACKAGE = "bamlss")
    },
    "lambda" = function(y, par, ...) {
      lambda <- 1/par$lambda
      yl <- y^lambda
      score <- 1/2 * (2 * (yl * (log(y) * lambda) * (yl - par$mu)))/par$sigma^2 - 1 - lambda * log(y)
      ifelse(y <= 0, 0, score)
    }
  )
  f$hess <- list(
    "mu" = function(y, par, ...) {
      .Call("cnorm_hess_mu",
        as.numeric(y^(1 / par$lambda)), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")), PACKAGE = "bamlss")
    },
    "sigma" = function(y, par, ...) {
      .Call("cnorm_hess_sigma",
        as.numeric(y^(1 / par$lambda)), as.numeric(par$mu), as.numeric(par$sigma),
        as.integer(attr(y, "check")), PACKAGE = "bamlss")
    },
    "lambda" = function(y, par, ...) {
      hess <- 1/2 * (2 * ((y^(1/par$lambda) * (log(y) * (par$lambda/par$lambda^2 -
        par$lambda * (2 * (par$lambda * par$lambda))/(par$lambda^2)^2)) -
        y^(1/par$lambda) * (log(y) * (par$lambda/par$lambda^2)) *
          (log(y) * (par$lambda/par$lambda^2))) * (y^(1/par$lambda) -
        par$mu) - y^(1/par$lambda) * (log(y) * (par$lambda/par$lambda^2)) *
        (y^(1/par$lambda) * (log(y) * (par$lambda/par$lambda^2)))))/par$sigma^2 -
        (par$lambda/par$lambda^2 - par$lambda * (2 * (par$lambda *
          par$lambda))/(par$lambda^2)^2) * log(y)
      ifelse(y <= 0, 0, -hess)
    }
  )
#  f$loglik <- function(y, par, ...) {
#    .Call("cnorm_power_loglik",
#      as.numeric(y), as.numeric(par$mu), as.numeric(par$sigma), as.numeric(par$lambda),
#      as.integer(attr(y, "check")), PACKAGE = "bamlss")
#  }
  f$d <- function(y, par, log = FALSE) {
    dy <- ifelse(y <= 0, pnorm(-par$mu / par$sigma, log.p = log),
      -0.918938533204673 - log(par$sigma) - 0.5 * (y^(1/par$lambda) - par$mu)^2/par$sigma^2)
    if(!log)
      dy <- exp(dy)
    dy
  }
  f$p <- function(y, par, log = FALSE) {
    y <- y^(1 / par$lambda)
    return(ifelse(y == 0, runif(length(y), 0, pnorm(0, par$mu, par$sigma)), pnorm(y - par$mu, 0, par$sigma)))
  }
  f$q <- function(y, par, ...) {
    rval <- qnorm(y^(1 / par$lambda)) * par$sigma + par$mu
    pmax(pmin(rval, Inf), 0)
  }
  f$initialize = list(
    "mu" = function(y, ...) {
      (y^(1 / start) + mean(y^(1 / start))) / 2
    },
    "sigma" = function(y, ...) {
      rep(sd(y^(1 / start)), length(y))
    },
    "lambda" = function(y, ...) { rep(start, length(y)) }
  )
  class(f) <- "family.bamlss"
  f
}


cens_bamlss <- function(links = c(mu = "identity", sigma = "log", df = "log"),
  left = 0, right = Inf, dist = "gaussian", ...)
{
  dist <- match.arg(dist, c("student", "gaussian", "logistic"))

  ddist <- switch(dist,
    "student"  = function(x, location, scale, df, log = TRUE)
      dt((x - location)/scale, df = df, log = log)/(scale^(1-log)) -
      log*log(scale),
    "gaussian" = function(x, location, scale, df, log = TRUE)
      dnorm((x - location)/scale, log = log)/(scale^(1-log)) -
      log*log(scale),
    "logistic" = function(x, location, scale, df, log = TRUE)
      dlogis((x - location)/scale, log = log)/(scale^(1-log)) -
      log*log(scale)
  )
  pdist <- switch(dist,
    "student"  = function(x, location, scale, df, lower.tail = TRUE,
      log.p = TRUE) pt((x - location)/scale, df = df, lower.tail = lower.tail,
      log.p = log.p),
    "gaussian" = function(x, location, scale, df, lower.tail = TRUE,
      log.p = TRUE) pnorm((x - location)/scale, lower.tail = lower.tail,
      log.p = log.p),
    "logistic" = function(x, location, scale, df, lower.tail = TRUE,
      log.p = TRUE) plogis((x - location)/scale, lower.tail = lower.tail,
      log.p = log.p)
  )

  gradfun <- function(y, par, type = "gradient", name = "mu") {
    ## functions used to evaluate gradient and hessian
    mills <- function(y, lower.tail = TRUE) {
      with(par, sigma * ddist(y, mu, sigma, df, log = FALSE)/
      pdist(y, mu, sigma, df, log.p = FALSE, lower.tail = lower.tail))
    }

    ## ddensity/dmu
    d1 <- with(par, switch(dist,
      "student"  = function(x)
        (x - mu)/sigma^2 * (df + 1) / (df + (x - mu)^2/sigma^2),
      "gaussian" = function(x)
        (x - mu)/sigma^2,
      "logistic" = function(x)
        (1 - 2 * pdist(-x, - mu, sigma, log.p = FALSE))/sigma)
    )

    ## ddensity/dsigma
    d2 <- function(x) with(par, d1(x) * (x-mu))

    ## d^2density/dmu^2
    d3 <- with(par, switch(dist,
      "student"  = function(x)
        (df + 1)*((x - mu)^2 - df*sigma^2) / (df*sigma^2 + (x - mu)^2)^2,
      "gaussian" = function(x)
        - 1/sigma^2,
      "logistic" = function(x)
        - 2/sigma * ddist(x, mu, sigma, log = FALSE))
    )

    ## d^2density/dsigma^2
    d5 <- with(par, switch(dist,
      "student"  = function(x)
        - (x - mu)^2 * (df + 1) / (df*sigma^2 + (x - mu)^2)^2*2*df*sigma^2,
      "gaussian" = function(x)
        2 * d3(x) * (x-mu)^2,
      "logistic" = function(x)
        - d2(x) - 2*(x-mu)^2/sigma*ddist(x,mu,sigma, log = FALSE)
    ))

    ## d^2density/dmudsigma
    d4 <- with(par, switch(dist,
      "student"  = function(x)
          d5(x) / (x - mu),
      "gaussian" = function(x)
        2 * d3(x) * (x-mu),
      "logistic" = function(x)
        - d1(x) + (x-mu)*d3(x)
    ))

    ## compute gradient
    if(type == "gradient") {
      if(name == "mu") {
        rval <- with(par, ifelse(y <= left,
          - mills(left)/sigma,
          ifelse(y >= right,
            mills(right, lower.tail = FALSE)/sigma,
            d1(y)
          )))
      } else {
        rval <- with(par, ifelse(y <= left,
          - mills(left) * (left - mu)/sigma,
          ifelse(y >= right,
            mills(right, lower.tail = FALSE) * (right - mu)/sigma,
            d2(y) - 1
          )))
      }

    ## compute hessian
    } else {
      if(name == "mu") {
        rval <- with(par, ifelse(y <= left,
          -d1(left)/sigma * mills(left) - mills(left)^2/sigma^2,
          ifelse(y >= right,
            d1(right)/sigma * mills(right, lower.tail = FALSE) -
              mills(right, lower.tail = FALSE)^2/sigma^2,
            d3(y)
          )))
      } else {
        rval <- with(par, ifelse(y <= left,
          ((left-mu)/sigma - (left-mu)*d2(left))*mills(left) -
            (left - mu)^2/sigma^2 * mills(left)^2,
          ifelse(y >= right,
            (-(right-mu)/sigma + (right-mu)*d2(right))*
              mills(right, lower.tail = FALSE)
              - (right - mu)^2/sigma^2 * mills(right, lower.tail = FALSE)^2,
            d5(y)
          )))
      }
    }
    return(rval)
  }

  score <- list(
    "mu" =  function(y, par, ...) {
      gradmu <- gradfun(y, par, type = "gradient", name = "mu")
      return(drop(gradmu))
    },
    "sigma" =  function(y, par, ...) {
      gradsigma <- gradfun(y, par, type = "gradient", name = "sigma")
      return(drop(gradsigma))
    }
  )

  hess <- list(
    "mu" =  function(y, par, ...) {
      wmu <- -1 * gradfun(y, par, type = "hess", name = "mu")
      return(drop(wmu))
    },
    "sigma" =  function(y, par, ...) {
      wsigma <- -1 * gradfun(y, par, type = "hess", name = "sigma")
      return(drop(wsigma))
    }
  )

  names <- switch(dist,
    "student" = c("mu", "sigma", "df"),
    "gaussian" = c("mu", "sigma"),
    "logistic" = c("mu", "sigma")
  )

  i <- 1:length(names)

  rval <- list(
    "family" = "cens",
    "names" = names,
    "links" = parse.links(links[i], c(mu = "identity", sigma = "log", df = "log")[i], ...),
    "d" = function(y, par, log = FALSE, ...) {
      ll <- with(par, ifelse(y <= left,
        pdist(left, mu, sigma, df, lower.tail = TRUE, log = TRUE),
        ifelse(y >= right,
          pdist(right, mu, sigma, df, lower.tail = FALSE, log = TRUE),
          ddist(y, mu, sigma, df, log = TRUE))))
      if(!log) ll <- exp(ll)
      return(ll)
    },
    "p" = function(y, par, log = FALSE, ...) {
      with(par, pdist(y, mu, sigma, df, lower.tail = TRUE, log = log))
    },
    "score" = score,
    "hess" = hess
  )

  class(rval) <- "family.bamlss"
  rval
}


## Truncated distributions.
dtrunc <- function(x, spec, a = 1, b = Inf, ...) {
  tt <- rep(0, length(x))
  g <- get(paste("d", spec, sep = ""), mode = "function")
  G <- get(paste("p", spec, sep = ""), mode = "function")
  tt[x >= a & x <= b] <- g(x[x >= a & x <= b], ...) / (G(b, ...) - G(a, ...))
  return(tt)
}

ptrunc <- function(x, spec, a = -Inf, b = Inf, ...)
{
  tt <- x
  aa <- rep(a, length(x))
  bb <- rep(b, length(x))
  G <- get(paste("p", spec, sep = ""), mode = "function")
  tt <- G(apply(cbind(apply(cbind(x, bb), 1, min), aa), 1, max), ...)
  tt <- tt - G(aa, ...)
  tt <- tt / (G(bb, ...) - G(aa, ...))
  return(tt)
}

qtrunc <- function(p, spec, a = -Inf, b = Inf, ...)
{
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(a, ...) + p * (G(b, ...) - G(a, ...)), ...)
  return(tt)
}

trunc2_bamlss <- function(links = c(mu = "identity", sigma = "log"),
  name = "norm", a = -Inf, b = Inf, ...)
{
  rval <- list(
    "family" = "trunc2",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "d" = function(y, par, log = FALSE) {
      if(name == "gamma") {
		    a2 <- par$sigma
		    s2 <- par$mu / par$sigma
        d <- dtrunc(y, name, a = a, b = b, shape = s2, scale = s2, log = log)
      } else {
        d <- dtrunc(y, name, a = a, b = b, mean = par$mu, sd = par$sigma, log = log)
      }
      return(d)
    },
    "p" = function(y, par, ...) {
      if(name == "gamma") {
		    a2 <- par$sigma
		    s2 <- par$mu / par$sigma
        p <- ptrunc(y, name, a = a, b = b, shape = a2, scale = s2, ...)
      } else {
        p <- ptrunc(y, name, a = a, b = b, mean = par$mu, sd = par$sigma, ...)
      }
      return(p)
    },
    "q" = function(y, par, ...) {
      q <- qtrunc(y, name, a = a, b = b, mean = par$mu, sd = par$sigma, ...)
      return(q)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


t_bamlss <- function(...)
{
  links <- c(mu = "identity", sigma2 = "log", df = "log")

  rval <- list(
    "family" = "t",
    "names" = c("mu", "sigma2", "df"),
    "links" = parse.links(links, c(mu = "identity", sigma2 = "log", df = "log"), ...),
    "bayesx" = list(
      "mu" = c("t", "mu"),
      "sigma2" =  c("t", "sigma2"),
	    "df" = c("t", "df")
    ),
    "bugs" = list(
      "dist" = "dt",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "mu" = function(par, ...) {
      rval <- par$mu
      rval[par$df <= 1] <- 0
      rval
    },
    "d" = function(y, par, log = FALSE) {
      arg <- (y - par$mu) / sqrt(par$sigma2)
      dt(arg, df = par$df, log = log)
    },
    "p" = function(y, par, ...) {
      arg <- (y - par$mu) / sqrt(par$sigma2)
      pt(arg, df = par$df, ...)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        ((par$df + 1) * (y - par$mu)) / (par$sigma2 * par$df + (y - par$mu)^2)
      },
      "sigma2" = function(y, par, ...) {
        -0.5 + ((par$df + 1) * (y - par$mu)^2) / (2 * (par$sigma2 * par$df + (y - par$mu)^2))
      },
      "df" = function(y, par, ...) {
        a <- par$df/2 * (digamma((par$df + 1) / 2) - digamma(par$df / 2) -
          log(1 + (y - par$mu)^2/(par$sigma2 * par$df)))
        b <- -0.5 + ((par$df + 1) * (y - par$mu)^2) / (2 * (par$sigma2 * par$df + (y - par$mu)^2))
        a + b
      }
    )
  )

  class(rval) <- "family.bamlss"
  rval
}


invgaussian_bamlss <- function(...)
{
  links <- c(mu = "log", sigma2 = "log")

  rval <- list(
    "family" = "invgaussian",
    "names" = c("mu", "sigma2"),
    "links" = parse.links(links, c(mu = "log", sigma2 = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu"  = c("invgaussian", "mu"),
      "sigma2" = c("invgaussian", "sigma2")
    ),
    "score" = list(
      "mu" = function(y, par, ...) {
        mu <- par$mu
        (y - mu) / (mu^2 * par$sigma2)
      },
      "sigma2" = function(y, par, ...) {
        mu <- par$mu
        -0.5 + (y - mu)^2 / (2 * y * mu^2 * par$sigma2)
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { 1 / (par$mu * par$sigma2) },
      "sigma2" = function(y, par, ...) { rep(0.5, length(y)) }
    ),
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      mu <- par$mu
      sigma <- sqrt(par$sigma2)
      d <- exp( -0.5 * log(2 * pi) - log(sigma) - (3 / 2) * log(y) - ((y - mu)^2) / (2 * sigma^2 * (mu^2) * y))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <- par$mu
      lambda <- 1 / sqrt(par$sigma2)
      lq <- sqrt(lambda / y)
      qm <- y / mu
      pnorm(lq * (qm - 1)) + exp(2 * lambda / mu) * pnorm(-lq * (qm + 1), ...)
    },
    "valid.response" = function(x) {
      if(is.factor(x) | is.character(x))
        stop("the response should be numeric!")
      return(TRUE)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}

weibull_bamlss <- function(...)
{
  links <- c(lambda = "log", alpha = "log")

  rval <- list(
    "family" = "weibull",
    "names" = c("lambda", "alpha"),
    "links" = parse.links(links, c(lambda = "log", alpha = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "lambda"  = c("weibull", "lambda"),
      "alpha" = c("weibull", "alpha")
    ),
    "bugs" = list(
      "dist" = "dweib",
      "eta" = BUGSeta,
      "model" = BUGSmodel,  ## FIXME: order of parameters?
      "order" = 2:1
    ),
    "mu" = function(par, ...) {
      alpha <-  par$alpha
      lambda <- par$lambda
      alpha * gamma(1 + 1 / lambda)
    },
    "d" = function(y, par, log = FALSE) {
      alpha <-  par$alpha
      lambda <- par$lambda
      dweibull(y, scale = lambda, shape = alpha, log = log)
    },
    "p" = function(y, par, ...) {
      alpha <- par$alpha
      lambda <- par$lambda
      pweibull(y, scale = lambda, shape = alpha, ...)
    },
    "q" = function(p, par, ...) {
      alpha <- par$alpha
      lambda <- par$lambda
      qweibull(p, scale = lambda, shape = alpha, ...)
    },
    "r" = function(n, par) {
      alpha <-  par$alpha
      lambda <- par$lambda
      rweibull(n, scale = lambda, shape = alpha)
    },
    "score" = list(
      "lambda" = function(y, par, ...) {
        sc <- par$lambda
        sh <- par$alpha
        -sh * (1 - (y/sc)^sh)
      },
      "alpha" = function(y, par, ...) {
        sc <- par$lambda
        sh <- par$alpha
        1 + sh*log(y/sc) - sh*((y/sc)^sh)*log(y/sc)
      }
    ),
    "hess" = list(
      "lambda" = function(y, par, ...) {
        par$alpha^2
      },
      "alpha" = function(y, par, ...) {
        ## page 16 in
        ## https://projecteuclid.org/download/suppdf_2/euclid.aoas/1437397122
        ## 1 + trigamma(2) + digamma(2)^2
        rep(1.823681, length(y))
      }
    ),
    "loglik" = function(y, par, ...) {
       sum(dweibull(y, shape = par$alpha, scale = par$lambda,
                    log = TRUE), na.rm = TRUE)
    },
    "initialize" = list(
      "lambda" = function(y, ...) {
        k <- 1.283 / var(log(y))
        exp(log(y) + 0.5772 / k)
      },
      "alpha" = function(y, ...) {
        k <- 1.283 / var(log(y))
        rep(k, length(y))
      }
    )
  )

  class(rval) <- "family.bamlss"
  rval
}


dw_bamlss <- function(...)
{
  links <- c(lambda = "logit", alpha = "log")
  
  trunc <- list(...)$trunc
  if(is.null(trunc))
    trunc <- FALSE

  link_lambda <- make.link2(links["lambda"])
  link_alpha <- make.link2(links["alpha"])

  rval <- list(
    "family" = "dw",
    "full.name" = "Discrete Weibull",
    "names" = c("lambda", "alpha"),
    "links" = links,
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x >= 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "d" = function(y, par, log = FALSE) {
      d <- par$lambda^(y^par$alpha) - par$lambda^((y + 1)^par$alpha)
      if(trunc) {
        d / par$lambda
      }
      if(log)
        d <- log(d)
      d[is.na(d) | !is.finite(d)] <- 1.490116e-08
      return(d)
    },
    "p" = function(y, par, log = FALSE) {
      d <- 1 - par$lambda^((y + 1)^par$alpha)
      if(log)
        d <- log(d)
      return(d)
    },
    "score" = list(
      "lambda" = function(y, par, ...) {        
        ya <- y^par$alpha
        y1a <- ((y + 1)^par$alpha)

        sw <- (par$lambda^(ya - 1) * ya - par$lambda^(y1a - 1) * y1a)/(par$lambda^ya - par$lambda^y1a)
        sw <- sw * link_lambda$mu.eta(link_lambda$linkfun(par$lambda))

        if(trunc)
          sw <- sw - 1/par$lambda * link_lambda$mu.eta(link_lambda$linkfun(par$lambda))

        return(sw)
      },
      "alpha" = function(y, par, ...) {
        ya <- y^par$alpha
        y1 <- y + 1
        ll <- log(par$lambda)
        lya <- par$lambda^(y1^par$alpha)
        lya2 <- par$lambda^ya

        sw <- (lya2 * (ll * (ya * log(y))) - lya * (ll * (y1^par$alpha * log(y1))))/(lya2 - lya)
        sw <- sw * link_alpha$mu.eta(link_alpha$linkfun(par$alpha))

        return(sw)
      }
    ),
    "hess" = list(
      "lambda" = function(y, par, ...) {
        d1 <- link_lambda$mu.eta(link_lambda$linkfun(par$lambda))
        d2 <- link_lambda$mu.eta2(link_lambda$linkfun(par$lambda))

        y1 <- y + 1
        ya <- y^par$alpha
        y1a <- y1^par$alpha
        ya1 <- (ya - 1)
        ly1 <- par$lambda^ya1
        ly1ya <- ly1 * ya
        lya <- par$lambda^ya
        lya1 <- par$lambda^y1a
        y1a1 <- (y1a - 1)
        lyalya1 <- (lya - lya1)
        ly1a1 <- par$lambda^y1a1
        a <- (ly1ya - ly1a1 * y1a)

        dll1 <- a/lyalya1
        dll2 <- (par$lambda^(ya1 - 1) * ya1 * ya - par$lambda^(y1a1 - 1) * y1a1 * y1a)/lyalya1 - a^2/lyalya1^2

        hw <- d1^2 * dll2 + d2 * dll1

        if(trunc) {
          hw <- hw + d1^2 * 1 / par$lambda^2 - d2 * 1 / par$lambda
        }

        return(-hw)
      },
      "alpha" = function(y, par, ...) {
        lambda <- par$lambda
        alpha <- par$alpha

        d1 <- link_alpha$mu.eta(link_alpha$linkfun(par$alpha))
        d2 <- link_alpha$mu.eta2(link_alpha$linkfun(par$alpha))

        ya <- y^par$alpha
        lla <- log(par$lambda)
        y1 <- y + 1
        ly <- log(y)
        ly1 <- log(y1)
        y1a <- y1^par$alpha
        ly1a <- par$lambda^y1a
        y1aly1 <- y1a * ly1
        yaly <- ya * ly
        llayaly <- lla * yaly
        lya <- lambda^ya
        lyaly1a <- lya - ly1a
        llay1aly1 <- lla * y1aly1
        lyallayaly <- lya * llayaly
        a <- ly1a * llay1aly1
        b <- lyallayaly - a

        dll1 <- b/lyaly1a

        dll2 <- (lyallayaly * llayaly + lya * (lla * (yaly * ly)) -
          (a * llay1aly1 + ly1a * (lla * (y1aly1 * ly1))))/lyaly1a - b^2/lyaly1a^2

        hw <- d1^2 * dll2 + d2 * dll1

        return(-hw)
      }
    ),
    "initialize" = list(
      "lambda" = function(y, ...) {
        m1 <- mean(y)
        m2 <- mean(y^2)
        rep(ifelse(any(y < 1), m1 / (m1 + 1), (m1 - 1) / m1), length(y))
      },
      "alpha" = function(y, ...) {
        rep(1, length(y))
      }
    )
  )

  class(rval) <- "family.bamlss"
  rval
}


gamma_bamlss <- function(...)
{
  rval <- list(
    "family" = "gamma",
    "names" = c("mu", "sigma"),
    "links" = c(mu = "log", sigma = "log"),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("gamma", "mu"),
      "sigma" = c("gamma", "sigma")
    ),
    "bugs" = list(
      "dist" = "dgamma",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "keras" = list(
      "nloglik" = function(y_true, y_pred) {
        K = keras::backend()

        mu = K$exp(y_pred[, 1])
        sigma = K$exp(y_pred[,2])

        stop("not working yet!")
      }
    ),
    "score" = list(
      "mu" = function(y, par, ...) {
        sigma <- par$sigma
        sigma * (-1 + y / par$mu)
      },
      "sigma" = function(y, par, ...) {
        mu <- par$mu
        sigma <- par$sigma
        sigma * (log(sigma) + 1 - log(mu) - digamma(sigma) + log(y) - y / mu)
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { par$sigma },
      "sigma" = function(y, par, ...) {
        sigma <- par$sigma
        sigma^2 * trigamma(sigma) - sigma
      }
    ),
    "loglik" = function(y, par, ...) {
       a <- par$sigma
       s <- par$mu / par$sigma
       sum(dgamma(y, shape = a, scale = s, log = TRUE), na.rm = TRUE)
    },
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      a <- par$sigma
      s <- par$mu / par$sigma
      dgamma(y, shape = a, scale = s, log = log)
    },
    "q" = function(p, par, lower.tail = TRUE, log.p = FALSE) {
      a <- par$sigma
      s <- par$mu / par$sigma
      qgamma(p, shape = a, scale = s, lower.tail = lower.tail, log.p = log.p)
    },
    "p" = function(y, par, lower.tail = TRUE, log.p = FALSE) {
      a <- par$sigma
      s <- par$mu / par$sigma
      pgamma(y, shape = a, scale = s, lower.tail = lower.tail, log.p = log.p)
    },
    "r" = function(n, par) {
      a <- par$sigma
      s <- par$mu / par$sigma
      rgamma(n, shape = a, scale = s)
    },
    "initialize" = list(
      "mu" = function(y, ...) { (y + mean(y)) / 2 },
      "sigma" = function(y, ...) { rep(1, length(y)) }
    ),
    "mean" = function(par) {
      a <- par$sigma
      s <- par$mu / par$sigma
      ex <- a * s
      return(ex)
    },
    "variance" = function(par) {
      a <- par$sigma
      s <- par$mu / par$sigma
      vx <- a * s^2
      return(vx)
    },
    "valid.response" = function(x) {
      if(is.factor(x) | is.character(x))
        stop("the response should be numeric!")
      return(TRUE)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


gengamma_bamlss <- function(...)
{
  links <- c(mu = "log", sigma = "log", tau = "log")

  rval <- list(
    "family" = "gengamma",
    "names" = c("mu", "sigma", "tau"),
    "links" = parse.links(links, c(mu = "log", sigma = "log", tau = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("gengamma", "mu"),
      "sigma" = c("gengamma", "sigma"),
      "tau" = c("gengamma", "tau")
    )
  )

  class(rval) <- "family.bamlss"
  rval
}


lognormal_bamlss <- function(...)
{
  links <- c(mu = "identity", sigma = "log")

  rval <- list(
    "family" = "lognormal",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("lognormal2", "mu"),
      "sigma" = c("lognormal2", "sigma")
    ),
    "score" = list(
      "mu" = function(y, par, ...) { (log(y) - par$mu) / (par$sigma^2) },
      "sigma" = function(y, par, ...) { -1 + (log(y) - par$mu)^2 / (par$sigma^2) }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { 1 / (par$sigma^2) },
      "sigma" = function(y, par, ...) { rep(2, length(y)) }
    ),
    "mu" = function(par, ...) {
      exp(par$mu + 0.5 * (par$sigma)^2)
    },
    "d" = function(y, par, log = FALSE) {
      dlnorm(y, meanlog = par$mu, sdlog = par$sigma, log = log)
    },
    "p" = function(y, par, ...) {
      plnorm(y, meanlog = par$mu, sdlog = par$sigma, ...)
    },
    "q" = function(p, par, ...) {
      qlnorm(p, meanlog = par$mu, sdlog = par$sigma, ...)
    },
    "r" = function(n, par) {
      rlnorm(n, meanlog = par$mu, sdlog = par$sigma)
    },
    "initialize" = list(
      "mu"    = function(y, ...) { (log(y) + mean(log(y))) / 2 },
      "sigma" = function(y, ...) { rep(sd(log(y)), length(y)) }
    ),
    "mean" = function(par) exp(par$mu + 0.5 * par$sigma^2),
    "variance" = function(par) exp(2 * par$mu + par$sigma^2) * (exp(par$sigma^2) - 1)
  )

  class(rval) <- "family.bamlss"
  rval
}


lognormal2_bamlss <- function(...)
{
  links <- c(mu = "identity", sigma2 = "log")

  rval <- list(
    "family" = "lognormal2",
    "names" = c("mu", "sigma2"),
    "links" = parse.links(links, c(mu = "identity", sigma2 = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("lognormal", "mu"),
      "sigma2" = c("lognormal", "sigma2")
    ),
	  "score" = list(
      "mu" = function(y, par, ...) { (log(y) - par$mu) / (par$sigma2) },
      "sigma2" = function(y, par, ...) { -0.5 + (log(y) - par$mu)^2 / (2 * par$sigma2) }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { 1 / (par$sigma2) },
      "sigma2" = function(y, par, ...) { rep(0.5, length(y)) }
    ),
    "mu" = function(par, ...) {
      exp(par$mu + 0.5 * (par$sigma2))
    },
    "d" = function(y, par, log = FALSE) {
      dlnorm(y, meanlog = par$mu, sdlog = sqrt(par$sigma2), log = log)
    },
    "p" = function(y, par, ...) {
      plnorm(y, meanlog = par$mu, sdlog = sqrt(par$sigma2), ...)
    },
    "r" = function(n, par) {
      rlnorm(n, meanlog = par$mu, sdlog = sqrt(par$sigma2))
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


gumbel_bamlss <- function(...)
{
 links <- c(mu = "identity", sigma = "log")

  rval <- list(
    "family" = "gumbel",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "keras" = list(
      "nloglik" = function(y_true, y_pred) {
        K = keras::backend()

        mu = y_pred[, 1]
        sigma = K$exp(y_pred[,2])

        ll = -1 * y_pred[, 2] + ((y_true[,1] - mu)/sigma) - K$exp((y_true[,1] - mu)/sigma)
        ll = K$sum(ll)

        return(-1 * ll)
      }
    ),
    "score" = list(
      "mu" = function(y, par, ...) {
        s1 <- 1/par$sigma
        return(-1 * (s1 - exp((y - par$mu)/par$sigma) * s1))
      },
      "sigma" = function(y, par, ...) {
        ymus <- (y - par$mu)/par$sigma
        return(-1 * (ymus + 1 - exp(ymus) * ymus))
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { 
        return(exp((y - par$mu)/par$sigma) * 1/par$sigma^2)
      },
      "sigma" = function(y, par, ...) {
        ymus <- (y - par$mu)/par$sigma
        h <- ymus - (exp(ymus) * ymus + exp(ymus) * ymus^2)
        return(-1 * h)
      }
    ),
    "d" = function(y, par, log = FALSE) {
      d <- -log(par$sigma) + ((y - par$mu)/par$sigma) - exp((y - par$mu)/par$sigma)
      if(!log)
        d <- exp(d)
      return(d)
    },
    "p" = function(y, par, ...) {
      1 - exp(-exp((y - par$mu)/par$sigma))
    },
    "r" = function(n, par) {
      rnorm(n, mean = par$mu, sd = par$sigma)
    },
    "q" = function(p, par) {
      q <- par$mu + par$sigma * log(-log(1 - p))
    },
    "initialize" = list(
      "mu"    = function(y, ...) { (y + mean(y)) / 2 },
      "sigma" = function(y, ...) { rep((sqrt(6) * sd(y))/pi, length(y)) }
    ),
    "valid.response" = function(x) {
      if(is.factor(x) | is.character(x))
        stop("the response should be numeric!")
      return(TRUE)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


dagum_bamlss <- function(...)
{
  links <- c(a = "log", b = "log", p = "log")

  rval <- list(
    "family" = "dagum",
    "names" = c("a", "b", "p"),
    "links" = parse.links(links, c(a = "log", b = "log", p = "log"), ...),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "a" = c("dagum", "a"),
      "b" = c("dagum", "b"),
      "p" = c("dagum", "p")
    ),
    "mu" = function(par, ...) {
      a <- par$a
      b <- par$b
      p <- par$p
      -(b/a) * (gamma(- 1 / a) * gamma(p + 1 / a)) / (gamma(p))
    },
    "d" = function(y, par, log = FALSE) {
      a <- par$a
      b <- par$b
      p <- par$p
      ap <- a * p
      d <- ap * y^(ap - 1) / (b^ap * (1 + (y / b)^a)^(p + 1))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      a <- par$a
      b <- par$b
      p <- par$p
      (1 + (y / b)^(-a))^(-p)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


BCCG2_bamlss <- function(...)
{
  links <- c(mu = "log", sigma = "log", nu = "identity")

  rval <- list(
    "family" = "BCCG",
    "names" = c("mu", "sigma", "nu"),
    "links" = parse.links(links, c(mu = "log", sigma = "log", nu = "identity"), ...),
	  "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x > 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "bayesx" = list(
      "mu" = c("BCCG", "mu"),
      "sigma" =  c("BCCG", "sigma"),
      "nu" = c("BCCG", "nu"),
      "order" = c("nu", "sigma", "mu")
    ),
    "d" = function(y, par, log = FALSE) {
      mu <- par$mu
      sigma <- par$sigma
      nu <- par$nu
      z <- ifelse(nu == 0, log(y/mu)/sigma, (((y/mu)^nu - 1)/(nu * sigma)))
      d <- (1 / (sqrt(2 * pi) * sigma)) * (y^(nu - 1) / mu^nu) * exp(-z^2 / 2)
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      mu <- par$mu
      sigma <- par$sigma
      nu <- par$nu
      z <- ifelse(nu == 0, log(y/mu)/sigma, (((y/mu)^nu - 1)/(nu * sigma)))
      FYy1 <- pnorm(z, ...)
      FYy2 <- ifelse(nu > 0, pnorm(-1/(sigma * abs(nu))), 0)
      FYy3 <- pnorm(1/(sigma * abs(nu)), ...)
      (FYy1 - FYy2)/FYy3
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


bivnorm_bamlss <- function(...)
{
  rval <- list(
    "family" = "mvnorm",
    "names" = c("mu1", "mu2", "sigma1", "sigma2", "rho"),
    "links" = c("identity", "identity", "log", "log", "rhogit"),
    "bayesx" = list(
      "mu1" = c("bivnormal", "mu"),
      "mu2" = c("bivnormal", "mu"),
      "sigma1" = c("bivnormal", "sigma"),
      "sigma2" = c("bivnormal", "sigma"),
      "rho" = c("bivnormal", "rho"),
      "order" = 5:1,
      "rm.number" = TRUE
    ),
    "d" = function(y, par, log = FALSE) {
      d <- .Call("bivnorm_loglik", as.numeric(y[, 1]), as.numeric(y[, 2]),
        as.numeric(par$mu1), as.numeric(par$mu2), as.numeric(par$sigma1),
        as.numeric(par$sigma2), as.numeric(par$rho), PACKAGE = "bamlss")
      if(!log)
        d <- exp(d)
      return(d)
    },
    "loglik" = function(y, par, log = FALSE) {
      ll <- .Call("bivnorm_loglik", as.numeric(y[, 1]), as.numeric(y[, 2]),
        as.numeric(par$mu1), as.numeric(par$mu2), as.numeric(par$sigma1),
        as.numeric(par$sigma2), as.numeric(par$rho), PACKAGE = "bamlss")
      return(ll)
    },
    "score" = list(
      "mu1" = function(y, par, ...) {
        1 / ((1 - par$rho^2) * par$sigma1) * ((y[, 1] - par$mu1) / par$sigma1 - par$rho * ((y[, 2] - par$mu2) / par$sigma2))
      },
      "mu2" = function(y, par, ...) {
        1 / ((1 - par$rho^2) * par$sigma2) * ((y[, 2] - par$mu2) / par$sigma2 - par$rho * ((y[, 1] - par$mu1) / par$sigma1))
      },
      "sigma1" = function(y, par, ...) {
        -1 + 1 / (1 - par$rho^2) * (y[, 1] - par$mu1) / par$sigma1 * ((y[, 1] - par$mu1) / par$sigma1 - par$rho * (y[, 2] - par$mu2) / par$sigma2)
      },
      "sigma2" = function(y, par, ...) {
        -1 + 1 / (1 - par$rho^2) * (y[, 2] - par$mu2) / par$sigma2 * ((y[, 2] - par$mu2) / par$sigma2 - par$rho * (y[, 1] - par$mu1) / par$sigma1)
      },
      "rho" = function(y, par, ...) {
        etarho <- par$rho / sqrt(1 - par$rho^2)
        sval <- 1 / (1 - par$rho^2) * par$rho / ((1 + etarho^2)^(1.5)) - etarho * (((y[, 1] - par$mu1) / par$sigma1)^2 +
          ((y[, 2] - par$mu2) / par$sigma2)^2) + (1 + 2 * etarho^2) / sqrt(1 + etarho^2) *
          (y[, 1] - par$mu1) / par$sigma1 * (y[, 2] - par$mu2) / par$sigma2
        sval
      }
    ),
    "hess" = list(
      "mu1" = function(y, par, ...) {
        1 / ((1 - par$rho^2) * par$sigma1^2)
      },
      "mu2" = function(y, par, ...) {
        1 / ((1 - par$rho^2) * par$sigma2^2)
      },
      "sigma1" = function(y, par, ...) {
        1 + 1 / (1 - par$rho^2)
      },
      "sigma2" = function(y, par, ...) {
        1 + 1 / (1 - par$rho^2)
      },
      "rho" = function(y, par, ...) {
        1 - par$rho^4
      }
    ),
    "p" = function(y, par, ...) {
      p <- cbind(
        pnorm(y[, 1], mean = par$mu1, sd = par$sigma1, ...),
        pnorm(y[, 2], mean = par$mu2, sd = par$sigma2, ...)
      )
      colnames(p) <- colnames(y)
      p
    },
    "mu" = function(y, par, ...) {
      cbind(
        "mu1" = par$mu1,
        "mu2" = par$mu2
      )
    },
    "initialize" = list(
      "mu1" = function(y, ...) {
        (y[, 1] + mean(y[, 1])) / 2
      },
      "mu2" = function(y, ...) {
        (y[, 2] + mean(y[, 2])) / 2
      },
      "sigma1" = function(y, ...) {
        rep(sd(y[, 1]), length(y[, 1]))
      },
      "sigma2" = function(y, ...) {
        rep(sd(y[, 2]), length(y[, 2]))
      },
      "rho" = function(y, ...) {
        rep(0, length(y[, 1]))
      }
    ),
    "mean" = function(par) {
      means_in_par <- grepl("mu", names(par))
      return(unlist(par[means_in_par]))
    },
    "variance" = function(par) {
      sigmas_in_par <- grepl("sigma", names(par))
      sigmas <- par[sigmas_in_par]
      vars <- sapply(sigmas, FUN = "^", e2 = 2)
      return(vars)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}

################################################################################
##  MULTIVARIATE NORMAL
################################################################################

dev_mvnorm_bamlss <- function(k = 2, cov_type = "cor", geo_fun = NULL,
                              dist_mat = NULL, ...)
{
  ## --- jump to univariate gaussian for k = 1 ---
  if(k == 1)
    return(gaussian_bamlss())

  ## --- jump to bivariate gaussion for k = 2 ---
  if(k == 2)
    return(bivnorm_bamlss())

  ## --- Set up family w/o correlations ---  

  ## --- set names of distributional  parameters ---
  mu <- paste0("mu", seq_len(k))
  sigma <- paste0("sigma", seq_len(k))

  ## --- set names of link functions ---
  links <- c(
    rep("identity", k),
    rep("log", k)
  )
  names(links) <- c(mu, sigma)

  ## --- family list ---
  rval <- list(
    "family" = "mvnorm",
    "names" = c(mu, sigma),
    "links" = links,
    "d" = function(y, par, log = FALSE) {     ## This density needs all correlation
      d <- log_dmvnorm(y, par)                ##   parameters
      if(!log)
        d <- exp(d)
      return(d)
    },
    "p" = function(y, par, ...) {
      p <- NULL
      for(j in 1:k) {
        p <- cbind(p, pnorm(y[, j],
          mean = par[[paste("mu", j, sep = "")]],
          sd = par[[paste("sigma", j , sep = "")]], ...))
      }
      colnames(p) <- colnames(y)
      p
    },
    "mu" = function(y, par, ...) {
      do.call("cbind", par[grep("mu", names(par))])
    }
  )

  ## --- score functions ---
  mu_score_calls <- paste0(
    "function(y, par, ...) {mu_score_mvnorm(y, par, j = ", seq_len(k),")}"
  )
  sigma_score_calls <- paste0(
    "function(y, par, ...) {sigma_score_mvnorm(y, par, j = ", seq_len(k),")}"
  )

  scores <- list()
  for(j in seq_along(mu)) {
    scores[[mu[j]]] <- eval(parse(text = mu_score_calls[j]))
  }
  for(j in seq_along(sigma)) {
    scores[[sigma[j]]] <- eval(parse(text = sigma_score_calls[j]))
  }
  rval$score <- scores

  ## --- initialize function ---
  mu_calls <- paste0(
    "function(y, ...) {(y[,", seq_len(k), "] + mean(y[,", seq_len(k), "])) / 2}"
  )
  sigma_calls <- paste0(
    "function(y, ...) {rep(sd(y[,", seq_len(k), "]), length(y[,", seq_len(k), "]))}"
  )

  init <- list()
  for(j in seq_along(mu))
    init[[mu[j]]] <- eval(parse(text = mu_calls[j]))
  for(j in seq_along(sigma))
    init[[sigma[j]]] <- eval(parse(text = sigma_calls[j]))
  rval$initialize <- init

  class(rval) <- "family.bamlss"

  ## --- add correlation model ---
  if(cov_type == "cor") {
    ## --- names ---
    rho <- combn(seq_len(k), 2, function(x) paste0("rho", x[1], x[2]))
    k_rho <- k * (k-1) / 2    ## number of rho parameters
    
    ## --- links ---
    nl <- names(rval$links)
    rval$links <- c(rval$links, rep("rhogit", k_rho))
    names(rval$links) <- c(nl, rho)
  
    ## --- add score functions ---
    rho_score_calls <- combn(seq_len(k), 2, function(x) {paste0(
      "function(y, par, ...) {rho_score_mvnorm(y, par, i = ", x[1],", j = ", x[2],")}"
    )})
    for(j in seq_along(rho)) {
      rval$score[[rho[j]]] <- eval(parse(text = rho_score_calls[j]))
    }

    ## --- init functions ---
    rho_calls <- rep("function(y, ...) { rep(0, length(y[,1])) }", k_rho)
    for(j in seq_along(rho)) {
      rval$initialize[[rho[j]]] <- eval(parse(text = rho_calls[j]))
    }

  }

  rval
}

mvnorm_bamlss <- function(k = 2, ...)
{
  ## --- jump to univariate gaussian for k = 1 ---
  if(k == 1)
    return(gaussian_bamlss())

  ## --- jump to bivariate gaussion for k = 2 ---
  if(k == 2)
    return(bivnorm_bamlss())

  ## --- set names of distributional  parameters ---
  mu <- paste0("mu", seq_len(k))
  sigma <- paste0("sigma", seq_len(k))
  rho <- combn(seq_len(k), 2, function(x) paste0("rho", x[1], x[2]))
  k_rho <- k * (k-1) / 2    ## number of rho parameters

  ## --- set names of link functions ---
  links <- c(
    rep("identity", k),
    rep("log", k),
    rep("rhogit", k_rho)
  )
  names(links) <- c(mu, sigma, rho)

  ## --- family list ---
  rval <- list(
    "family" = "mvnorm",
    "names" = c(mu, sigma, rho),
    "links" = links,
    "d" = function(y, par, log = FALSE) {
      d <- log_dmvnorm(y, par)
      if(!log)
        d <- exp(d)
      return(d)
    },
    "p" = function(y, par, ...) {
      p <- NULL
      for(j in 1:k) {
        p <- cbind(p, pnorm(y[, j],
          mean = par[[paste("mu", j, sep = "")]],
          sd = par[[paste("sigma", j , sep = "")]], ...))
      }
      colnames(p) <- colnames(y)
      p
    },
    "mu" = function(y, par, ...) {
      do.call("cbind", par[grep("mu", names(par))])
    }
  )

  ## --- score functions ---
  mu_score_calls <- paste0(
    "function(y, par, ...) {mu_score_mvnorm(y, par, j = ", seq_len(k),")}"
  )
  sigma_score_calls <- paste0(
    "function(y, par, ...) {sigma_score_mvnorm(y, par, j = ", seq_len(k),")}"
  )
  rho_score_calls <- combn(seq_len(k), 2, function(x) {paste0(
    "function(y, par, ...) {rho_score_mvnorm(y, par, i = ", x[1],", j = ", x[2],")}"
  )})

  scores <- list()
  for(j in seq_along(mu)) {
    scores[[mu[j]]] <- eval(parse(text = mu_score_calls[j]))
  }
  for(j in seq_along(sigma)) {
    scores[[sigma[j]]] <- eval(parse(text = sigma_score_calls[j]))
  }
  for(j in seq_along(rho)) {
    scores[[rho[j]]] <- eval(parse(text = rho_score_calls[j]))
  }
  rval$score <- scores

  ## --- initialize function ---
  mu_calls <- paste0(
    "function(y, ...) {(y[,", seq_len(k), "] + mean(y[,", seq_len(k), "])) / 2}"
  )
  sigma_calls <- paste0(
    "function(y, ...) {rep(sd(y[,", seq_len(k), "]), length(y[,", seq_len(k), "]))}"
  )
  rho_calls <- rep("function(y, ...) { rep(0, length(y[,1])) }", k_rho)

  init <- list()
  for(j in seq_along(mu)) {
    init[[mu[j]]] <- eval(parse(text = mu_calls[j]))
  }
  for(j in seq_along(sigma)) {
    init[[sigma[j]]] <- eval(parse(text = sigma_calls[j]))
  }
  for(j in seq_along(rho)) {
    init[[rho[j]]] <- eval(parse(text = rho_calls[j]))
  }
  rval$initialize <- init

  class(rval) <- "family.bamlss"
  rval
}


log_dmvnorm <- function(y, par) {
  npar <- names(par)
  par <- do.call("cbind", par)
  y <- as.matrix(y)
  sj <- grep("sigma", npar)
  mj <- grep("mu", npar)
  rj <- as.integer(min(grep("rho", npar)))
  return(.Call("log_dmvnorm", y, par, nrow(y), ncol(y), mj, sj, rj, PACKAGE = "bamlss"))
}

mu_score_mvnorm <- function(y, par, j) {
  par <- do.call("cbind", par)
  y <- as.matrix(y)
  cn <- colnames(par)
  mj <- grep("mu", cn)
  sj <- grep("sigma", cn)
  rj <- as.integer(min(grep("rho", cn)))
  kj <- as.integer(j-1)
  return(.Call("mu_score_mvnorm", y, par, nrow(y), ncol(y), mj, sj, rj, kj, PACKAGE = "bamlss"))
}

mu_score_mvnormR <- function(y, par, j) {
  n <- nrow(y)
  k <- ncol(y)
  par <- do.call("cbind", par)
  y <- as.matrix(y)
  cn <- colnames(par)
  sj <- grep("sigma", cn)
  mj <- grep("mu", cn)
  rj <- as.integer(min(grep("rho", cn)))

  rval <- NULL
  for (kk in seq(n)) {
    ## Fill Sigma
    Sigma <- matrix(NA, ncol=k, nrow=k)
    l <- 0
    for ( ii in seq(k) ) {
      Sigma[ii,ii] <- par[kk,k+ii]^2
      for ( jj in seq(k) ) {
        if ( ii<jj ) {
          Sigma[ii,jj] <- par[kk,k+ii] * par[kk,k+jj] * par[kk,rj+l]
          Sigma[jj,ii] <- Sigma[ii,jj]
          l <- l + 1
        }
      }
    }
    ## invert Sigma
    InvSig <- chol2inv(chol(Sigma))

    m <- drop(y[kk,] - par[kk,mj])
    rval[kk] <- sum( InvSig[j,] * m )
  }

  return(rval)
}

sigma_score_mvnorm <- function(y, par, j) {
  par <- do.call("cbind", par)
  y <- as.matrix(y)
  cn <- colnames(par)
  sj <- grep("sigma", cn)
  mj <- grep("mu", cn)
  rj <- as.integer(min(grep("rho", cn)))
  kj <- as.integer(j-1)
  return(.Call("sigma_score_mvnorm", y, par, nrow(y), ncol(y), mj, sj, rj, kj, PACKAGE = "bamlss"))
}

sigma_score_mvnormR <- function(y, par, j) {
 n <- nrow(y)
 k <- ncol(y)
 par <- do.call("cbind", par)
 y <- as.matrix(y)
 cn <- colnames(par)
 sj <- grep("sigma", cn)
 mj <- grep("mu", cn)
 rj <- as.integer(min(grep("rho", cn)))

 rval <- NULL
 for ( kk in seq(n) ) {
 ## Fill Rho
   Rho <- matrix(0, ncol=k, nrow=k)
   l <- 0
   for ( ii in seq(k) ) {
     Rho[ii,ii] <- 1
     for ( jj in seq(k) ) {
       if ( ii<jj ) {
         Rho[ii,jj] <- par[kk,rj+l]
         Rho[jj,ii] <- Rho[ii,jj]
         l <- l + 1
       }
     }
   }
   ## invert Rho
   InvRho <- chol2inv(chol(Rho))

   m <- drop((y[kk,] - par[kk,mj])/par[kk,sj])
   rval[kk] <- -1 + m[j]*sum(m*InvRho[j,])
 }
 return(rval)
}

rho_score_mvnorm <- function(y, par, i, j)
{
  par <- do.call("cbind", par)
  y <- as.matrix(y)
  cn <- colnames(par)
  sj <- grep("sigma", cn)
  mj <- grep("mu", cn)
  rj <- as.integer(min(grep("rho", cn)))
  kj <- as.integer(j-1)
  lj <- as.integer(i-1)
  rval <- .Call("rho_score_mvnorm", y, par, nrow(y), ncol(y), mj, sj, rj, kj, lj, PACKAGE = "bamlss")
  return(rval)
}

rho_score_mvnormR <- function(y, par, i, j)
{
 n <- nrow(y)
 k <- ncol(y)
 par <- do.call("cbind", par)
 y <- as.matrix(y)
 cn <- colnames(par)
 sj <- grep("sigma", cn)
 mj <- grep("mu", cn)
 rj <- as.integer(min(grep("rho", cn)))

 rval <- NULL
 for ( kk in seq(n) ) {
   ## Fill Rho
   Rho <- matrix(0, ncol=k, nrow=k)
   l <- 0
   for ( ii in seq(k) ) {
     Rho[ii,ii] <- 1
     for ( jj in seq(k) ) {
       if ( ii<jj ) {
         Rho[ii,jj] <- par[kk,rj+l]
         Rho[jj,ii] <- Rho[ii,jj]
         l <- l + 1
       }
     }
   }
   ## compute deriv
   mu <- Rho[i,j]
   eta <- mu / sqrt(1 - mu^2)
   deriv <- 1 / (1 + eta^2)^1.5

   ## invert Rho
   InvRho <- chol2inv(chol(Rho))

   m <- drop((y[kk,] - par[kk,mj])/par[kk,sj])
   rval[kk] <- drop(-.5*InvRho[j,i] + .5*sum(InvRho[j,]*m)*sum(InvRho[,i]*m))*deriv
 }
 return(rval)
}


bivprobit_bamlss <- function(...)
{
  links <- c(mu1 = "identity", mu2 = "identity", rho = "rhogit")

  rval <- list(
    "family" = "bivprobit",
    "names" = c("mu1", "mu2", "rho"),
    "links" = parse.links(links, c(mu1 = "identity", mu2 = "identity", rho = "rhogit"), ...),
    "bayesx" = list(
      "mu1" = c("bivprobit", "mu"),
      "mu2" = c("bivprobit", "mu"),
      "rho" = c("bivprobit", "rho"),
      "order" = 3:1,
      "rm.number" = TRUE
    ),
    "mu" = function(par, ...) {
      c(par$mu1, par$mu2)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


bivlogit_bamlss <- function(...)
{
  links <- c(p1 = "logit", p2 = "logit", psi = "log")

  rval <- list(
    "family" = "bivlogit",
    "names" = c("p1", "p2", "psi"),
    "links" = parse.links(links, c(p1 = "logit", p2 = "logit", psi = "log"), ...),
    "bayesx" = list(
      "mu1" = c("bivlogit", "mu"),
      "mu2" = c("bivlogit", "mu"),
      "psi" = c("bivlogit", "oddsratio"),
      "order" = 3:1,
      "rm.number" = TRUE
    ),
    "mu" = function(par, ...) {
      c(par$mu1, par$mu2)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


mvt_bamlss <- function(...) 
{
  links <- c(mu1 = "identity", mu2 = "identity",
    sigma1 = "log", sigma2 = "log", rho = "rhogit", df = "log")
  rval <- list(
    "family" = "mvt",
    "names" = c("mu1", "mu2", "sigma1", "sigma2", "rho", "df"),
    "links" = parse.links(links, c(mu1 = "identity", mu2 = "identity",
      sigma1 = "log", sigma2 = "log", rho = "rhogit", df = "log"), ...),
    "bayesx" = list(
      "mu1" = c("bivt", "mu"),
      "mu2" = c("bivt", "mu"),
      "sigma1" = c("bivt", "sigma"),
      "sigma2" = c("bivt", "sigma"),
      "rho" = c("bivt", "rho"),
      "df" = c("bivt", "df"),
      "order" = 6:1,
      "rm.number" = TRUE
    ),
    "mu" = function(par, ...) {
      c(par$mu1, par$mu2)
    },
    "d" = function(y, par, log = FALSE, ...) {
      c1 <- 1 + 1/(par$df*(1 - par$rho)) * ( ((y[,1] - par$mu1)/par$sigma1)^2 -
        2*par$rho*((y[,1] - par$mu1)/par$sigma1) * ((y[,2] - par$mu2)/par$sigma2) +
        ((y[,2] - par$mu2)/par$sigma2)^2 )
      d <- lgamma((par$df + 2)*0.5) - lgamma(par$df*0.5) - 1.1447298858494 - log(par$df) -
        log(par$sigma1) - log(par$sigma2) - 0.5*log(1 - par$rho^2) -
        (par$df + 2)*0.5 * log(c1)
      if(!log)
        d <- exp(d)
      return(d)
    },
   "hess" = list(
      "mu1" = function(y, par, ...) {
        1 / ((1 - par$rho^2) * par$sigma1^2)
      },
      "mu2" = function(y, par, ...) {
        1 / ((1 - par$rho^2) * par$sigma2^2)
      },
      "sigma1" = function(y, par, ...) {
        1 + 1 / (1 - par$rho^2)
      },
      "sigma2" = function(y, par, ...) {
        1 + 1 / (1 - par$rho^2)
      },
      "rho" = function(y, par, ...) {
        1 - par$rho^4
      }
    )
  )
  class(rval) <- "family.bamlss"
  rval
}


dirichlet_bamlss <- function(...)
{
  k <- list(...)$k
  if(is.null(k))
    stop("argument k is missing, number of parameters need to be specified!")
  rval <- list(
    "family" = "dirichlet",
    "names" = paste0("alpha", 1:k)
  )
  rval$links <- rep("log", k)
  names(rval$links) <- rval$names

  rval$d <- function(y, par, log = FALSE, ...) {
    par <- do.call("cbind", par)
    s1 <- -1 * rowSums(lgamma(par))
    s2 <- lgamma(rowSums(par))
    s3 <- rowSums((par - 1) * log(y))
    ll <- s1 + s2 + s3
    if(!log)
      ll <- exp(ll)
    ll
  }
  score <- function(y, par, id, ...) {
    rval <- -1 * par[[id]] * digamma(par[[id]]) +
      par[[id]] * digamma(rowSums(do.call("cbind", par))) +
      par[[id]] * log(y[, grep(id, names(par))])
    rval
  }
  rval$score <- rep(list(score), k)
  names(rval$score) <- rval$names
  hess <- function(y, par, id, ...) {
    rval <- par[[id]]^2 * (-1 * trigamma(rowSums(do.call("cbind", par))) +
      trigamma(par[[id]]))
    rval
  }
  rval$hess <- rep(list(hess), k)
  names(rval$hess) <- rval$names

  rval$initialize <- list()
  for(j in 1:k) {
    ft <- c(
      "function(y, ...) {",
         paste0("mj <- mean(y[, ", j, "]);"),
         "return(rep(mj, nrow(y)))",
      "}"
    )
    rval$initialize[[paste0("alpha", j)]] <- eval(parse(text = paste(ft, collapse = "")))
  }

  rval$expectation <- function(par, log = FALSE) {
    if (log) {
      par <- do.call("cbind", par)
      alpha0 <- rowSums(par)
      par <- digamma(par) - digamma(alpha0)
      as.list(as.data.frame(par))
    } else {
      par <- do.call("cbind", par)
      alpha0 <- rowSums(par)
      par <- par / alpha0
      as.list(as.data.frame(par))
    }
  }

  class(rval) <- "family.bamlss"
  rval
}

#if(FALSE) {
#  data("ArcticLake", package = "DirichletReg")
# 
#  AL <- as.matrix(ArcticLake[, 1:3])
# 
#  d <- data.frame("depth" = ArcticLake$depth)
#  d$y <- AL
# 
#  f <- list(
#    y ~ s(depth),
#      ~ s(depth),
#      ~ s(depth)
#  )
# 
#  b <- bamlss(f, data = d, family = dirichlet_bamlss(k = 3))
# 
#  p <- predict(b, type = "parameter")
#  p <- as.data.frame(p)
#  p <- p / rowSums(p)
#  par(mfrow = c(1, 3))
#  for(j in 1:3) {
#    plot(d$y[, j] ~ d$depth, ylab = colnames(d$y)[j], xlab = "depth")
#    plot2d(p[, j] ~ d$depth, lwd = 2, add = TRUE)
#  }

#  ###
#  p <- predict(b, type = "parameter", FUN = function(x) x)
#  p <- lapply(p, as.matrix)
#  alpha_0 <- purrr::reduce(p, `+`)
#  props <- lapply(p, function(x) x / alpha_0)
#  props_95 <- lapply(props, function(x) t(apply(x, 1L, c95)))

#  cols <- hcl.colors(3, "Dark 2")
#  for(j in 1:3) {
#    plot2d(props_95[[j]] ~ d$depth, lwd = 2,
#	   ylab = colnames(d$y)[j], xlab = "depth")
#    points(d$y[, j] ~ d$depth, pch = 19, col = cols[j])
#    points(d$y[, j] ~ d$depth)
#  }

#  
#}


multinomial_bamlss <- multinom_bamlss <- function(...)
{
  link <- "log"

  rval <- list(
    "family" = "multinomial",
    "names" = "pi",
    "links" = link,
    "bugs" = list(
      "dist" = "dcat",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "bayesx" = list(
      "pi" = c(paste("multinom", "probit", sep = "_"), "main", "servant")
    ),
    "score" = function(y, par, id, ...) {
      pi <- par[[id]] / (1 + rowSums(do.call("cbind", par)))
      if(is.factor(y))
        1 * (y == id)
      else
        return(y[, id] - pi)
    },
    "hess" = function(y, par, id, ...) {
      pi <- par[[id]] / (1 + rowSums(do.call("cbind", par)))
      return(pi * (1 - pi))
    },
    "d" = function(y, par, log = FALSE) {
      if(is.factor(y))
        y <- model.matrix(~ y - 1)
      par <- cbind(do.call("cbind", par), 1)
      d1 <- rowSums(y * log(par))
      d2 <- log(rowSums(par))
      d <- d1 - d2
      if(!log)
        d <- exp(d)
      return(d)
    },
    "loglik" = function(y, par, ...) {
      if(is.factor(y))
        y <- model.matrix(~ y - 1)
      par <- cbind(do.call("cbind", par), 1)
      d1 <- rowSums(y * log(par))
      d2 <- log(rowSums(par))
      d <- d1 - d2
      return(sum(d, na.rm = TRUE))
    },
    "cat" = TRUE
  )

  rval$probabilities <- function(par, numeric = TRUE, ...) {
    pi <- exp(do.call("cbind", par))
    pi <- pi / (rowSums(pi) + 1)
    pi <- cbind(1 - rowSums(pi), pi)

    if(!numeric) {
      pi <- t(apply(pi, 1, function(x) {
        y <- rep(0,length(x))
        y[which.max(x)] <- 1
        y
      }))
    }

    colnames(pi) <- c("reference", names(par))

    return(as.data.frame(pi))
  }

  class(rval) <- "family.bamlss"
  rval
}


## Count Data distributions.
poisson_bamlss <- function(...)
{
  rval <- list(
    "family" = "poisson",
    "names" = "lambda",
    "links" = c(lambda = "log"),
    "bayesx" = list(
      "lambda" = c("poisson", "lambda")
    ),
    "bugs" = list(
      "dist" = "dpois",
      "eta" = BUGSeta,
      "model" = BUGSmodel
    ),
    "mu" = function(par, ...) {
       par$lambda
    },
    "d" = function(y, par, log = FALSE) {
      dpois(y, lambda = par$lambda, log = log)
    },
    "p" = function(y, par, ...) {
      ppois(y, lambda = par$lambda, ...)
    },
    "r" = function(n, par) {
      rpois(n, lambda = par$lambda)
    },
    "q" = function(p, par) {
      qpois(p, lambda = par$lambda)
    },
    "score" = list(
      "lambda" = function(y, par, ...) {
        y - par$lambda
      }
    ),
    "hess" = list(
      "lambda" = function(y, par, ...) {
        par$lambda
      }
    ),
    initialize = list(
      "lambda" = function(y, ...) {
        (y + mean(y)) / 2
      }
    ),
    "mean" = function(par) par$lambda,
    "variance" = function(par) par$lambda,
    "valid.response" = function(x) {
      if(is.factor(x) | is.character(x))
        stop("the response should be numeric!")
      if(any(x < 0))
        stop("response values < 0 are not allowed!")
      return(TRUE)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


zip_bamlss <- function(...)
{
  links <- c(lambda = "log", pi = "logit")

  rval <- list(
    "family" = "zip",
    "names" = c("lambda", "pi"),
    "links" = parse.links(links, c(lambda = "log", pi = "logit"), ...),
    "bayesx" = list(
      "lambda" = c("zip", "lambda"),
      "pi" = switch(links["pi"],
        "logit" = c("zip", "pi"),
        "cloglog2" = c("zip", "pi")
      )
    ),
	  "mu" = function(par, ...) {
      par$lambda * (1 - par$pi)
    },
    "d" = function(y, par, log = FALSE) {
      d <- ifelse(y == 0, par$pi + (1 - par$pi) * dpois(y, lambda = par$lambda),
				(1 - par$pi) * dpois(y, lambda = par$lambda))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      ifelse(y < 0, 0, par$pi + (1 - par$pi) * ppois(y, lambda = par$lambda))
    }
  )
  if(rval$bayesx[[2]][[1]] == "zip_pi_cloglog")
    rval$bayesx[[1]][[1]] <- "zip_lambda_cloglog"

  class(rval) <- "family.bamlss"
  rval
}


hurdleP_bamlss <- function(...)
{
  links <- c(lambda = "log", pi = "logit")

  rval <- list(
    "family" = "hurdle",
    "names" = c("lambda", "pi"),
    "links" = parse.links(links, c(lambda = "log", pi = "logit"), ...),
    "bayesx" = list(
      "lambda" = c("hurdle", "lambda"),
      "pi" = c("hurdle", "pi"),
      "hess" = list(
        "lambda" = function(x) { 1 * (x != 0)}
      )
    ),
    "mu" = function(par, ...) {
      (1 - par$pi) * par$lambda / (1 - exp(-par$lambda))
    },
    "d" = function(y, par, log = FALSE) {
      d <- ifelse(y == 0, par$pi,
        (1 - par$pi) * dpois(y, lambda = par$lambda) / (1 - exp(-par$lambda)))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      cdf1 <- ppois(y, lambda = par$lambda)
      cdf2 <- ppois(0, lambda = par$lambda)
      cdf3 <- par$pi + ((1 - par$pi) * (cdf1 - cdf2)/(1 - cdf2))
      cdf <- ifelse((y == 0), par$pi, cdf3)
      cdf
    }
  )

  class(rval) <- "family.bamlss"
  rval
}

negbin_bamlss <- function(...)
{
  links <- c(mu = "log", delta = "log")

  rval <- list(
    "family" = "negbin",
    "names" = c("mu", "delta"),
    "links" = parse.links(links, c(mu = "log", delta = "log"), ...),
    "bayesx" = list(
     "mu" = c("negbin", "mu"),
      "delta" = c("negbin", "delta")
    ),
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      dnbinom(y, mu = par$mu, size = par$delta, log = log)
    },
    "p" = function(y, par, ...) {
      pnbinom(y, mu = par$mu, size = par$delta)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


zinb_bamlss <- function(...)
{
  links <- c(mu = "log", pi = "logit", delta = "log")

  rval <- list(
    "family" = "zinb",
    "names" = c("mu", "pi", "delta"),
    "links" = parse.links(links, c(mu = "log", pi = "logit", delta = "log"), ...),
    "bayesx" = list(
      "mu" = c("zinb", "mu"),
      "pi" = c("zinb", "pi"),
      "delta" = c("zinb", "delta")
    ),
    "mu" = function(par, ...) {
      par$mu * (1 - par$pi)
    },
    "d" = function(y, par, log = FALSE) {
      d <- ifelse(y == 0, par$pi + (1 - par$pi) * dnbinom(y, mu = par$mu, size = par$delta),
				(1 - par$pi) * dnbinom(y, mu = par$mu, size = par$delta))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      ifelse(y<0, 0, par$pi + (1 - par$pi) * pnbinom(y, size = par$delta, mu = par$mu))
    }
  )

  class(rval) <- "family.bamlss"
  rval
}

hurdleNB_bamlss <- function(...)
{
  links <- c(mu = "log", pi = "logit", delta = "log")

  rval <- list(
    "family" = "hurdleNB",
    "names" = c("mu", "delta", "pi"),
    "links" = parse.links(links, c(mu = "log", delta = "log", pi = "logit"), ...),
    "bayesx" = list(
      "mu" = c("hurdle", "mu"),
      "delta" = c("hurdle", "delta"),
      "pi" = c("hurdle", "pi"),
      "hess" = list(
         "mu" = function(x) { 1 * (x != 0)},
         "delta" = function(x) { 1 * (x != 0)}
      )
    ),
    "mu" = function(par, ...) {
      (1 - par$pi) * par$mu / (1 - (par$delta) / (par$delta + par$mu)^par$delta)
    },
    "d" = function(y, par, log = FALSE) {
      d <- ifelse(y == 0, par$pi + (1 - par$pi) * dnbinom(y, mu = par$mu, size = par$delta),
        (1 - par$pi) * dnbinom(y, mu = par$mu, size = par$delta))
      if(log) d <- log(d)
      d
    },
    "p" = function(y, par, ...) {
      cdf1 <- pnbinom(y, size = par$delta, mu = par$mu)
      cdf2 <- pnbinom(0, size = par$delta, mu = par$mu)
      cdf3 <- par$pi + ((1 - par$pi) * (cdf1 - cdf2)/(1 - cdf2))
      cdf <- ifelse((y == 0), par$pi, cdf3)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}

ztnbinom_bamlss <- function(...) {
### Zero-truncated neg bin
### Author: Thorsten Simon
### Date:   2018 Nov
  rval <- list(
    "family" = "ztnbinom",
    "names"  = c("mu", "theta"),
    "links"  = c(mu = "log", theta = "log"),
    "mu" = function(par, ...) {
      par$mu / stats::pnbinom(0, mu = par$mu, size = par$theta, lower.tail = FALSE)
    },
    "loglik" = function(y, par, ...) {
        rval <- stats::dnbinom(y, mu = par$mu, size = par$theta, log = TRUE) -
                stats::pnbinom(0, mu = par$mu, size = par$theta, lower.tail = FALSE, log.p = TRUE)
        sum(rval) 
    },
    "d" = function(y, par, log = FALSE) {
        rval <- stats::dnbinom(y, mu = par$mu, size = par$theta, log = TRUE) -
                stats::pnbinom(0, mu = par$mu, size = par$theta, lower.tail = FALSE, log.p = TRUE)
        if(log) rval else exp(rval)
    },
    "p" = function(y, par, ...) {
        rval <- log(stats::pnbinom(y, mu = par$mu, size = par$theta,
                                   lower.tail = TRUE, log.p = FALSE) -
                    stats::dnbinom(0, mu = par$mu, size = par$theta)) -
                    stats::pnbinom(0, mu = par$mu, size = par$theta,
                                   lower.tail = FALSE, log.p = TRUE)
        exp(rval)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        .Call("ztnbinom_score_mu", as.numeric(y), as.numeric(par$mu),
              as.numeric(par$theta), PACKAGE = "bamlss")
      },
      "theta" = function(y, par, ...) {
        .Call("ztnbinom_score_theta", as.numeric(y), as.numeric(par$mu),
              as.numeric(par$theta), PACKAGE = "bamlss")
      }
    ),
    "initialize" = list(
      "mu" = function(y, ...) {
##        rep(mean(y) - 1, length(y))
        (y + mean(y) - 1) / 2
      },
      "theta" = function(y, ...) {
##        rep( mean(y)^2 / (var(y) - mean(y)), length(y))
        rep( 1, length(y))
      }
    )
  )

  rval$rps <- function(y, par, ymin = 1L, ymax = max(max(y), 100L)) {
    K <- seq(ymin, ymax, by = 1L)
    rps <- rep(0, length(y))
    for (k in K) {
      P <- rval$p(k, par)
      O <- y <= k
      rps <- rps + (P - O)^2
    }
    return(rps)
  }
  rval$type <- "discrete"

  class(rval) <- "family.bamlss"
  rval
}

ztnbinomVAR_bamlss <- function(...) {
### Zero-truncated neg bin with alternative parameterization
###   gamma models the VAR of the un-truncated neg-bin, and replaces theta
### Author: Thorsten Simon
### Date:   2018 Nov
  theta_from_gamma <- function(par) {
    par$mu / sqrt(par$gamma - par$mu)  ## Thus gamma has to be greater than mu! Do we have to enforce this?
  }
  rval <- list(
    "family" = "ztnbinom",
    "names"  = c("mu", "gamma"),
    "links"  = c(mu = "log", gamma = "log"),
    "mu" = function(par, ...) {
        theta <- theta_from_gamma(par)
        par$mu / stats::pnbinom(0, mu = par$mu, size = theta, lower.tail = FALSE)
    },
    "loglik" = function(y, par, ...) {
        theta <- theta_from_gamma(par)
        rval <- stats::dnbinom(y, mu = par$mu, size = theta, log = TRUE) -
                stats::pnbinom(0, mu = par$mu, size = theta, lower.tail = FALSE, log.p = TRUE)
        sum(rval) 
    },
    "d" = function(y, par, log = FALSE) {
        theta <- theta_from_gamma(par)
        rval <- stats::dnbinom(y, mu = par$mu, size = theta, log = TRUE) -
                stats::pnbinom(0, mu = par$mu, size = theta, lower.tail = FALSE, log.p = TRUE)
        if(log) rval else exp(rval)
    },
    "p" = function(y, par, ...) {
        theta <- theta_from_gamma(par)
        rval <- log(stats::pnbinom(y, mu = par$mu, size = theta,
                                   lower.tail = TRUE, log.p = FALSE) -
                    stats::dnbinom(0, mu = par$mu, size = theta)) -
                    stats::pnbinom(0, mu = par$mu, size = theta,
                                   lower.tail = FALSE, log.p = TRUE)
        exp(rval)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        theta <- theta_from_gamma(par)
        .Call("ztnbinom_score_mu", as.numeric(y), as.numeric(par$mu),
              as.numeric(theta), PACKAGE = "bamlss")
      },
      "gamma" = function(y, par, ...) {
        theta <- theta_from_gamma(par)
        .Call("ztnbinom_score_theta", as.numeric(y), as.numeric(par$mu),
              as.numeric(theta), PACKAGE = "bamlss")
      }
    ),
    "initialize" = list(
      "mu" = function(y, ...) {
        (y + mean(y) - 1) / 2
      },
      "gamma" = function(y, ...) {
        (y + mean(y) - 1) / 2
      }
    )
  )

  class(rval) <- "family.bamlss"
  rval
}

## http://stats.stackexchange.com/questions/17672/quantile-regression-in-jags
quant_bamlss <- function(prob = 0.5)
{
  rval <- list(
    "family" = "quant",
    "names" = "mean",
    "links" = c("mean" = "identity"),
    "bayesx" = list(
      "mean" = c("quantreg", "mean"),
      "quantile" = prob
    )
  )
  class(rval) <- "family.bamlss"
  rval
}


quant2_bamlss <- function(prob = 0.5, ...)
{
  links <- c(mu = "identity", sigma = "log")
  rval <- list(
    "family" = "quant2",
    "names" = c("mu", "sigma"),
    "links" = parse.links(links, c(mu = "identity", sigma = "log"), ...),
    "bugs" = list(
      "dist" = "dnorm",
      "eta" = BUGSeta,
      "model" = BUGSmodel,
      "reparam" = c(
        "mu" = "(1 - 2 * prop) / (prop * (1 - prop)) * w[i] + mu",
        "sigma" = "(prop * (1 - prop) * (1 / sigma)) / (2 * w[i])"
      ),
      "addparam" = list("w[i] ~ dexp(1 / sigma[i])"),
      "addvalues" = list("prop" = prob)
    )
  )
  class(rval) <- "family.bamlss"
  rval
}

kde_bamlss <- function(..., err)
{
  rval <- list(
    "family" = "kde",
    "names" = "mu",
    "links" = c("mu" = "identity"),
    "d" = function(y, par, log = FALSE) {
      n <- length(y)
      K <- function(x) { (1/sqrt(2 * pi) * exp(-0.5 * (x)^2)) }
      r <- function(x) { exp(x) / (1 + exp(x)) }
      d <- rep(0, n)
      e <- y - par$mu
      h <- 1.06 * sd(y) * n^(-1/5)
      for(i in 1:n) {
        z <- e[i] - err[-i]
        d[i] <- log(1/n * sum(1/h * K(z / h)))
      }
      if(!log)
        d <- exp(d)
      d
    }
  )
  class(rval) <- "family.bamlss"
  rval
}


## Zero adjusted families.
#zero_bamlss <- function(g = invgaussian)
#{
#  pi <- "logit"
#  gg <- try(inherits(g, "family.bamlss"), silent = TRUE)
#  if(inherits(gg, "try-error")) {
#    g <- deparse(substitute(g), backtick = TRUE, width.cutoff = 500)
#  } else {
#    if(is.function(g)) {
#      if(inherits(try(g(), silent = TRUE), "try-error"))
#        g <- deparse(substitute(g), backtick = TRUE, width.cutoff = 500)
#    }
#  }
#  g <- bamlss.family(g)
#  g[c("mu", "map2par", "loglik")] <- NULL
#  g0 <- g
#  np <- g$names
#  g$links <- c(g$links, "pi" = pi)
#  g$family <- paste("zero-adjusted | ", g$family, ", ", "binomial", sep = "")
#  g$names <- c(g$names, "pi")
#  g$valid.response <- function(x) {
#    if(any(x < 0)) stop("response includes values smaller than zero!")
#    TRUE
#  }
#  dfun <- g0$d
#  pfun <- g0$p
#  if(is.function(dfun)) {
#    g$d <- function(y, par, log = TRUE) {
#      d <- dfun(y, par, log = FALSE) * par$pi * (y > 0) + (1 - par$pi) * (y == 0)
#      if(log)
#        d <- log(d)
#      d
#    }
#  }
#  if(is.function(pfun)) {
#    g$p <- function(y, par, log = FALSE) {
#      pfun(y, par, log = log) * par$pi + (1 - par$pi)
#    }
#  }
#  g$bayesx <- c(g$bayesx, list("pi" = c(paste("binomial", pi, sep = "_"), "meanservant")))
#  g$bayesx$hess <- list()
#  for(j in np) {
#    g$bayesx$hess[[j]] <- function(x) { 1 * (x > 0) }
#    if(grepl("mean", g$bayesx[[j]][2]))
#      g$bayesx[[j]][2] <- "meanservant"
#  }
#  g$bayesx$zero <- TRUE
#  class(g) <- "family.bamlss"
#  g
#}


## General bamlss family creator.
gF <- function(x, ...) {
  if(is.function(x)) {
    if(inherits(x(), "gamlss.family"))
      return(tF(x, ...))
  }

  if(!is.character(x))
    x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
  fam <- get(paste(x, "bamlss", sep = "_"), mode = "function")

  return(fam(...))
}

gF2 <- function(x, ...) {
  if(!is.character(x))
    x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
  F <- get(paste(x, "bamlss", sep = "_"), mode = "function")
  x <- F(...)
  linkinv <- vector(mode = "list", length = length(x$names))
  for(j in x$names)
    linkinv[[j]] <- make.link2(x$links[j])$linkinv
  x$map2par <- function(eta) {
    for(j in names(eta)) {
      eta[[j]] <- linkinv[[j]](eta[[j]])
      eta[[j]][is.na(eta[[j]])] <- 0
      if(any(jj <- eta[[j]] == Inf))
        eta[[j]][jj] <- 10
      if(any(jj <- eta[[j]] == -Inf))
        eta[[j]][jj] <- -10
    }
    return(eta)
  }
  if(is.null(x$loglik))
    x$loglik <- function(y, par, ...) { sum(x$d(y, par, log = TRUE), na.rm = TRUE) }
  x
}


## Function to transform gamlss.family objects.
tF <- function(x, ...)
{
  if(is.function(x)) x <- x()
  if(!inherits(x, "gamlss.family")) stop('only "gamlss.family" objects can be transformed!')

  args <- list(...)
  bd <- if(is.null(args$bd)) 1 else args$bd
  args$bd <- NULL
  pr <- args$range
  check_range <- function(par) {
    for(j in names(par)) {
      if(!is.null(pr[[j]])) {
        par[[j]][par[[j]] < min(pr[[j]])] <- min(pr[[j]])
        par[[j]][par[[j]] > max(pr[[j]])] <- max(pr[[j]])
      }
    }
    par
  }
  nx <- names(x$parameters)
  score <- hess <- initialize <- list()

  make_call <- function(fun) {
    fn <- deparse(substitute(fun), backtick = TRUE, width.cutoff = 500)
    nf <- names(formals(fun))
    if(length(nf) < 1) {
      call <- paste(fn, "()", sep = "")
    } else {
      call <- paste(fn, "(", if("y" %in% nf) "y," else "", sep = "")
      np <- nx[nx %in% nf]
      call <- paste(call, paste(np, '=', 'par$', np, sep = '', collapse = ','), sep = "")
      if("bd" %in% nf) {
        call <- paste(call, ",bd=", bd, sep = "")
      }
    }
    call <- parse(text = paste(call, ")", sep = ""))
    return(call)
  }

  if("mu" %in% nx) {
    mu.link <- make.link2(x$mu.link)
    mu.cs <- make_call(x$dldm)
    mu.hs <- make_call(x$d2ldm2)
    score$mu  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(mu.cs) * mu.link$mu.eta(mu.link$linkfun(par$mu))
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    hess$mu <- function(y, par, ...) {
      par <- check_range(par)
      score <- eval(mu.cs)
      hess <- -1 * eval(mu.hs)
      eta <- mu.link$linkfun(par$mu)
      res <- drop(score * mu.link$mu.eta2(eta) + hess * mu.link$mu.eta(eta)^2)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    if(!is.null(x$mu.initial)) {
      initialize$mu <- function(y, ...) {
        if(!is.null(attr(y, "contrasts"))) {
          if(!is.null(dim(y)))
            y <- y[, ncol(y)]
        }
        if(!is.null(bd))
          bd <- rep(bd, length.out = if(!is.null(dim(y))) nrow(y) else length(y))
        res <- eval(x$mu.initial)
        if(!is.null(dim(res)))
          res <- res[, 1]
        res
      }
    }
  }

  if("sigma" %in% nx) {
    sigma.link <- make.link2(x$sigma.link)
    sigma.cs <- make_call(x$dldd)
    sigma.hs <- make_call(x$d2ldd2)
    score$sigma  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(sigma.cs) * sigma.link$mu.eta(sigma.link$linkfun(par$sigma))
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    hess$sigma <- function(y, par, ...) {
      par <- check_range(par)
      score <- eval(sigma.cs)
      hess <- -1 * eval(sigma.hs)
      eta <- sigma.link$linkfun(par$sigma)
      res <- drop(score * sigma.link$mu.eta2(eta) + hess * sigma.link$mu.eta(eta)^2)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    if(!is.null(x$sigma.initial)) {
      initialize$sigma <- function(y, ...) {
        if(!is.null(bd))
          bd <- rep(bd, length.out = if(!is.null(dim(y))) nrow(y) else length(y))
        res <- eval(x$sigma.initial)
        if(!is.null(dim(res)))
          res <- res[, 1]
        res
      }
    }
  }

  if("nu" %in% nx) {
    nu.link <- make.link2(x$nu.link)
    nu.cs <- make_call(x$dldv)
    nu.hs <- make_call(x$d2ldv2)
    score$nu  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(nu.cs) * nu.link$mu.eta(nu.link$linkfun(par$nu))
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    hess$nu <- function(y, par, ...) {
      par <- check_range(par)
      score <- eval(nu.cs)
      hess <- -1 * eval(nu.hs)
      eta <- nu.link$linkfun(par$nu)
      res <- drop(score * nu.link$mu.eta2(eta) + hess * nu.link$mu.eta(eta)^2)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    if(!is.null(x$nu.initial)) {
      initialize$nu <- function(y, ...) {
        if(!is.null(bd))
          bd <- rep(bd, length.out = if(!is.null(dim(y))) nrow(y) else length(y))
        res <- eval(x$nu.initial)
        if(!is.null(dim(res)))
          res <- res[, 1]
        res
      }
    }
  }

  if("tau" %in% nx) {
    tau.link <- make.link2(x$tau.link)
    tau.cs <- make_call(x$dldt)
    tau.hs <- make_call(x$d2ldt2)
    score$tau  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(tau.cs) * tau.link$mu.eta(tau.link$linkfun(par$tau))
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    hess$tau <- function(y, par, ...) {
      par <- check_range(par)
      score <- eval(tau.cs)
      hess <- -1 * eval(tau.hs)
      eta <- tau.link$linkfun(par$tau)
      res <- drop(score * tau.link$mu.eta2(eta) + hess * tau.link$mu.eta(eta)^2)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
    if(!is.null(x$tau.initial)) {
      initialize$tau <- function(y, ...) {
        if(!is.null(bd))
          bd <- rep(bd, length.out = if(!is.null(dim(y))) nrow(y) else length(y))
        res <- eval(x$tau.initial)
        if(!is.null(dim(res)))
          res <- res[, 1]
        res
      }
    }
  }

  dfun <- get(paste("d", x$family[1], sep = ""))
  pfun <- try(get(paste("p", x$family[1], sep = "")), silent = TRUE)
  qfun <- try(get(paste("q", x$family[1], sep = "")), silent = TRUE)
  rfun <- try(get(paste("r", x$family[1], sep = "")), silent = TRUE)

  nf <- names(formals(dfun))
  bdc <- "bd" %in% nf

  dc <- parse(text = paste('dfun(y,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',log=log,...',
    if(bdc) paste0(",bd=", bd) else NULL, ")", sep = ""))
  pc <- parse(text = paste('pfun(q,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',log=log,...',
    if(bdc) paste0(",bd=", bd) else NULL, ")", sep = ""))
  qc <- parse(text = paste('qfun(p,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',log=log,...',
    if(bdc) paste0(",bd=", bd) else NULL, ")", sep = ""))
  rc <- parse(text = paste('rfun(n,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',...',
    if(bdc) paste0(",bd=", bd) else NULL, ")", sep = ""))

  rval <- list(
    "family" = x$family[1],
    "names" = nx,
    "links" = unlist(x[paste(nx, "link", sep = ".")]),
    "score" = score,
    "hess" = hess,
    "d" = function(y, par, log = FALSE, ...) {
       par <- check_range(par)
       d <- try(eval(dc), silent = TRUE)
       if(inherits(d, "try-error"))
         d <- rep(NA, length(par[[1L]]))
       return(d)
    },
    "p" = if(!inherits(pfun, "try-error")) function(q, par, log = FALSE, ...) {
      par <- check_range(par)
      eval(pc)
    } else NULL,
    "q" = if(!inherits(qfun, "try-error")) function(p, par, log = FALSE, ...) {
      par <- check_range(par)
      eval(qc)
    } else NULL,
    "r" = if(!inherits(rfun, "try-error")) function(n, par, ...) {
      par <- check_range(par)
      eval(rc)
    } else NULL
  )
  names(rval$links) <- nx
  rval$valid.response <- x$y.valid
  rval$initialize <- initialize
  rval$type <- tolower(x$type)

  if(!is.null(x$mean)) {
    meanc <- make_call(x$mean)
    rval$mean  <- function(par, ...) {
      par <- check_range(par)
      res <- eval(meanc)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
  } else {
    rval$mean <- function(par, ...) { par[[1L]] }
  }

  if(!is.null(x$variance)) {
    varc <- make_call(x$variance)
    rval$variance  <- function(par, ...) {
      par <- check_range(par)
      res <- eval(varc)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
  } else {
    rval$variance <- function(par, ...) { par[[2L]] }
  }

  class(rval) <- "family.bamlss"
  rval
}


### Ordered logit.
### http://staff.washington.edu/lorenc2/bayesian/ologit.R

## Categorical distribution.
dcat <- function(x, p, log=FALSE)
{
  if(is.vector(x) & !is.matrix(p))
    p <- matrix(p, length(x), length(p), byrow = TRUE)
  if(is.matrix(x) & !is.matrix(p))
    p <- matrix(p, nrow(x), length(p), byrow = TRUE)
  if(is.vector(x) & {length(x) == 1}) {
    temp <- rep(0, ncol(p))
    temp[x] <- 1
    x <- t(temp)
  } else if(is.vector(x) & (length(x) > 1)) x <- as.indicator.matrix(x)
  if(!identical(nrow(x), nrow(p))) stop("The number of rows of x and p differ.")
  if(!identical(ncol(x), ncol(p))) {
    x.temp <- matrix(0, nrow(p), ncol(p))
    x.temp[,as.numeric(colnames(x))] <- x
    x <- x.temp
  }
  dens <- x * p
  if(log) dens <- x * log(p)
  dens <- as.vector(rowSums(dens))
  return(dens)
}

qcat <- function(pr, p, lower.tail = TRUE, log.pr = FALSE)
{
  if(!is.vector(pr)) pr <- as.vector(pr)
  if(!is.vector(p)) p <- as.vector(p)
  if(log.pr == FALSE) {
    if(any(pr < 0) | any(pr > 1))
      stop("pr must be in [0,1].")
  } else if(any(!is.finite(pr)) | any(pr > 0)) stop("pr, as a log, must be in (-Inf,0].")
  if(sum(p) != 1) stop("sum(p) must be 1.")
  if(lower.tail == FALSE) pr <- 1 - pr
  breaks <- c(0, cumsum(p))
  if(log.pr == TRUE) breaks <- log(breaks)
  breaks <- matrix(breaks, length(pr), length(breaks), byrow = TRUE)
  x <- rowSums(pr > breaks)
  return(x)
}

rcat <- function(n, p)
{
  if(is.vector(p)) {
    x <- as.vector(which(rmultinom(n, size = 1, prob = p) == 1, arr.ind = TRUE)[, "row"])
  } else {
    n <- nrow(p)
    x <- apply(p, 1, function(x) {
      as.vector(which(rmultinom(1, size = 1, prob = x) == 1, arr.ind = TRUE)[, "row"])
    })
  }
  return(x)
}

as.indicator.matrix <- function(x)
{
  n <- length(x)
  x <- as.factor(x)
  X <- matrix(0, n, length(levels(x)))
  X[(1:n) + n*(unclass(x) - 1)] <- 1
  dimnames(X) <- list(names(x), levels(x))
  X
}


#####################################
## New family specification setup. ##
#####################################
##############################################
## (1) Score, hessian and fisher functions. ##
##############################################
##########################
## Normal distribution. ##
##########################
snorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"))
{
  if(!is.character(which))
    which <- c("mu", "sigma")[as.integer(which)]
  which <- tolower(which)
  score <- NULL
  dxm <- x - mean
  sd2 <- sd^2
  for(w in which) {
    if(w == "mu")
      score <- cbind(score, dxm / sd2)
    if(w == "sigma")
      score <- cbind(score, (dxm^2 - sd2) / sd^3)
  }
  if(is.null(dim(score)))
    score <- matrix(score, ncol = 1)
  colnames(score) <- paste("d", which, sep = "")
  score
}

hnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"))
{
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  n <- length(x)
  hess <- list()
  sd2 <- sd^2
  for(w in which) {
    if(w == "mu")
      hess[[w]] <- rep(-1 / sd2, length.out = n)
    if(w == "sigma")
      hess[[w]] <- (sd2 - 3 * (x - mean)^2) / sd2^2
    if(w == "mu.sigma")
      hess[[w]] <- -2 * (x - mean) / sd^3
    if(w == "sigma.mu")
      hess[[w]] <- 2 * (mean - x) / sd^3
  }
  hess <- do.call("cbind", hess)
  colnames(hess) <- gsub("mu", "dmu", colnames(hess))
  colnames(hess) <- gsub("sigma", "dsigma", colnames(hess))
  colnames(hess)[colnames(hess) == "dmu"] <- "d2mu"
  colnames(hess)[colnames(hess) == "dsigma"] <- "d2sigma"
  hess
}

fnorm <- function(x, mean = 0, sd = 1, which = c("mu", "sigma"))
{
  if(!is.character(which))
    which <- c("mu", "sigma", "mu.sigma", "sigma.mu")[as.integer(which)]
  which <- tolower(which)
  n <- length(x)
  fish <- list()
  sd2 <- sd^2
  for(w in which) {
    if(w == "mu")
      fish[[w]] <- rep(1 / sd2, length.out = n)
    if(w == "sigma")
      fish[[w]] <- rep(2 / sd2, length.out = n)
    if(w %in% c("mu.sigma", "sigma.mu"))
      fish[[w]] <- 0
  }
  fish <- do.call("cbind", fish)
  colnames(fish) <- gsub("mu", "dmu", colnames(fish))
  colnames(fish) <- gsub("sigma", "dsigma", colnames(fish))
  colnames(fish)[colnames(fish) == "dmu"] <- "d2mu"
  colnames(fish)[colnames(fish) == "dsigma"] <- "d2sigma"
  fish
}


###################################
## (2) Family creator functions. ##
###################################
get.dist <- function(distribution = "norm")
{
  ddist <- get(paste("d", distribution, sep = ""))
  pdist <- try(get(paste("p", distribution, sep = "")), silent = TRUE)
  if(inherits(pdist, "try-error"))
    pdist <- NULL
  qdist <- try(get(paste("q", distribution, sep = "")), silent = TRUE)
  if(inherits(qdist, "try-error"))
    qdist <- NULL
  rdist <- try(get(paste("r", distribution, sep = "")), silent = TRUE)
  if(inherits(pdist, "try-error"))
    rdist <- NULL
  sdist <- try(get(paste("s", distribution, sep = "")), silent = TRUE)
  if(inherits(sdist, "try-error"))
    sdist <- NULL
  hdist <- try(get(paste("h", distribution, sep = "")), silent = TRUE)
  if(inherits(hdist, "try-error"))
    hdist <- NULL
  fdist <- try(get(paste("f", distribution, sep = "")), silent = TRUE)
  if(inherits(fdist, "try-error"))
    fdist <- NULL
  dist <- list("d" = ddist, "p" = pdist, "q" = qdist, "r" = rdist,
    "s" = sdist, "h" = hdist, "f" = fdist)
  return(dist)
}


gaussian5_bamlss <- function(links = c(mu = "identity", sigma = "log"), ...)
{
  links <- parse.links(links, c(mu = "identity", sigma = "log"), ...)
  lfun <- list()
  for(j in names(links))
    lfun[[j]] <- make.link2(links[j])

  rval <- list(
    "family" = "gaussian",
    "names" = c("mu", "sigma"),
    "links" = links,
    "score" = list(
      "mu" = function(y, par, ...) {
        mu <- lfun$mu$linkfun(par$mu)
        drop(snorm(y, mean = par$mu, sd = par$sigma, which = 1) * lfun$mu$mu.eta(mu))
      },
      "sigma" = function(y, par, ...) {
        sigma <- lfun$sigma$linkfun(par$sigma)
        drop(snorm(y, mean = par$mu, sd = par$sigma, which = 2) * lfun$sigma$mu.eta(sigma))
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        mu <- lfun$mu$linkfun(par$mu)
        w <- -1 * drop(snorm(y, mean = par$mu, sd = par$sigma, which = 1) * lfun$mu$mu.eta2(mu) +
          hnorm(y, mean = par$mu, sd = par$sigma, which = 1) * lfun$mu$mu.eta(mu)^2)
        w
      },
      "sigma" = function(y, par, ...) {
        sigma <- lfun$sigma$linkfun(par$sigma)
        w <- -1 * drop(snorm(y, mean = par$mu, sd = par$sigma, which = 2) * lfun$sigma$mu.eta2(sigma) +
          hnorm(y, mean = par$mu, sd = par$sigma, which = 2) * lfun$sigma$mu.eta(sigma)^2)
        w
      }
    ),
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "p" = function(y, par, ...) {
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


mvnormAR1_bamlss <- function(k = 2, ...)
{
  if(k == 1)
    return(gaussian_bamlss())

  if(k == 2)
    return(bivnorm_bamlss())

  mu <- paste("mu", 1:k, sep = "")
  sigma <- paste("sigma", 1:k, sep = "")
  rho <- "rho"

  links <- c(rep("identity", length(mu)), rep("log", length(sigma)), rep("rhogit", length(rho)))
  names(links) <- c(mu, sigma, rho)

  rval <- list(
    "family" = "mvnormAR1",
    "names" = c(mu, sigma, rho),
    "links" = links,
    "d" = function(y, par, log = FALSE) {
      d <- log_dmvnormAR1(y, par)
      if(!log)
        d <- exp(d)
      return(d)
    },
    "p" = function(y, par, ...) {
      p <- NULL
      for(j in 1:k) {
        p <- cbind(p, pnorm(y[, j],
          mean = par[[paste("mu", j, sep = "")]],
          sd = par[[paste("sigma", j , sep = "")]], ...))
      }
      colnames(p) <- colnames(y)
      p
    },
    "mu" = function(y, par, ...) {
      do.call("cbind", par[grep("mu", names(par))])
    }
  )

  mu_score_calls <- sigma_score_calls <- NULL
  for(j in seq_along(mu))
    mu_score_calls <- c(mu_score_calls, paste("function(y, par, ...) {mu_score_mvnormAR1(y, par, j=",j,")}", sep=""))
  for(j in seq_along(sigma))
    sigma_score_calls <- c(sigma_score_calls, paste("function(y, par, ...) {sigma_score_mvnormAR1(y, par, j=",j,")}", sep=""))
  scores <- list()
  for(j in seq_along(mu))
    scores[[mu[j]]] <- eval(parse(text = mu_score_calls[j]))
  for(j in seq_along(sigma))
    scores[[sigma[j]]] <- eval(parse(text = sigma_score_calls[j]))
  scores[["rho"]] <- function(y, par, ...) { rho_score_mvnormAR1(y, par) }
  rval$score <- scores

  mu_calls <- sigma_calls <- NULL
  for(j in seq_along(mu))
    mu_calls <- c(mu_calls, paste("function(y, ...) { (y[,", j, "] + mean(y[,", j , "])) / 2 }"))
  for(j in seq_along(sigma))
    sigma_calls <- c(sigma_calls, paste("function(y, ...) { rep(sd(y[,", j , "]), length(y[,", j, "])) }"))
  init <- list()
  for(j in seq_along(mu))
    init[[mu[j]]] <- eval(parse(text = mu_calls[j]))
  for(j in seq_along(sigma))
    init[[sigma[j]]] <- eval(parse(text = sigma_calls[j]))
  init[["rho"]] <- function(y, ...) { rep(0, length(y[, 1])) }
  rval$initialize <- init

  class(rval) <- "family.bamlss"
  rval
}

log_dmvnormAR1 <- function(y, par)
{
  par <- do.call("cbind", par)
  y <- as.matrix(y)
  cn <- colnames(par)
  mj <- grep("mu", cn)
  sj <- grep("sigma", cn)
  rj <- as.integer(grep("rho", cn))
  return(.Call("log_dmvnormAR1", y, par, nrow(y), ncol(y), mj, sj, rj, PACKAGE = "bamlss"))
}

log_dmvnormAR1_R <- function(y, par)
{
  par <- do.call("cbind", par)
  y <- as.matrix(y)
  cn <- colnames(par)
  mj <- grep("mu", cn)
  sj <- grep("sigma", cn)
  rj <- as.integer(grep("rho", cn))
  k <- ncol(y)
  n <- nrow(y)

  rval <- NULL
  for (kk in seq(n)) {
    term1 <- - k/2 * log(2*pi)
    term2 <- - sum( log( par[kk,sj] ))
    term3 <- - (k-1)/2 * log( 1 - par[kk,rj]^2 )
    ytilde <- ( y[kk,] - par[kk,mj] ) / par[kk,sj]
    term4 <- - 1/( 1 - par[kk,rj]^2 ) * 1/2 *
               ( sum(ytilde^2) - 2 * par[kk,rj] * sum( ytilde[-1]*ytilde[-k] ) +
                 par[kk,rj]^2 * sum( ytilde[-c(1,k)]^2 ) )

    rval[kk] <- term1 + term2 + term3 + term4
  }

  return(rval)
}

mu_score_mvnormAR1 <- function(y, par, j)
{
  par <- do.call("cbind", par)
  y <- as.matrix(y)
  cn <- colnames(par)
  mj <- grep("mu", cn)
  sj <- grep("sigma", cn)
  rj <- as.integer(grep("rho", cn))
  kj <- as.integer(j-1)
  return(.Call("mu_score_mvnormAR1", y, par, nrow(y), ncol(y), mj, sj, rj, kj, PACKAGE = "bamlss"))
}

sigma_score_mvnormAR1 <- function(y, par, j)
{
  par <- do.call("cbind", par)
  y <- as.matrix(y)
  cn <- colnames(par)
  sj <- grep("sigma", cn)
  mj <- grep("mu", cn)
  rj <- as.integer(grep("rho", cn))
  kj <- as.integer(j-1)
  return(.Call("sigma_score_mvnormAR1", y, par, nrow(y), ncol(y), mj, sj, rj, kj, PACKAGE = "bamlss"))
}

rho_score_mvnormAR1 <- function(y, par)
{
  par <- do.call("cbind", par)
  y <- as.matrix(y)
  cn <- colnames(par)
  sj <- grep("sigma", cn)
  mj <- grep("mu", cn)
  rj <- as.integer(grep("rho", cn))
  return(.Call("rho_score_mvnormAR1", y, par, nrow(y), ncol(y), mj, sj, rj, PACKAGE = "bamlss"))
}


# -------------------------------------------------------------------
# Generalized logistic type I distribution with gradients.
# -------------------------------------------------------------------
glogis_bamlss <- function(...) {

   requireNamespace('glogis')

   links <- c(mu="identity",sigma="log",alpha="log")

   rval <- list(
      "family" = "glogis", # Generalized Logistic Distribution Type I (a.k.a. skewed logistic distribution)
      "names"  = c("mu","sigma","alpha"),
      "score" = list(
         "mu" = function(y,par,...) {
            .Call("bamlss_glogis_score",as.integer(1),as.numeric(y),
                  as.numeric(par$mu),as.numeric(par$sigma),as.numeric(par$alpha))
         },
         "sigma" = function(y,par,...) {
            .Call("bamlss_glogis_score",as.integer(2),as.numeric(y),
                  as.numeric(par$mu),as.numeric(par$sigma),as.numeric(par$alpha))
         },
         "alpha" = function(y,par,...) {
            .Call("bamlss_glogis_score",as.integer(3),as.numeric(y),
                  as.numeric(par$mu),as.numeric(par$sigma),as.numeric(par$alpha))
         }
      ),
      "hess" = list(
         "mu" = function(y,par,...) {
            .Call("bamlss_glogis_hesse",as.integer(1),as.numeric(y),
                  as.numeric(par$mu),as.numeric(par$sigma),as.numeric(par$alpha))
         },
         "sigma" = function(y,par,...) {
            .Call("bamlss_glogis_hesse",as.integer(2),as.numeric(y),
                  as.numeric(par$mu),as.numeric(par$sigma),as.numeric(par$alpha))
         },
         "alpha" = function(y,par,...) {
            .Call("bamlss_glogis_hesse",as.integer(3),as.numeric(y),
                  as.numeric(par$mu),as.numeric(par$sigma),as.numeric(par$alpha))
         }
      ),
      "links"  = parse.links(links,c(mu="identity",sigma="log",alpha="log"),...),
      "loglik" = function(y, par, ... ) {
         .Call("bamlss_glogis_loglik",as.numeric(y),
                  as.numeric(par$mu),as.numeric(par$sigma),as.numeric(par$alpha))
      },
      "d" = function(y,par,log=FALSE) {
        sapply(y, FUN = function(x) {
          .Call("bamlss_glogis_density",as.numeric(x),
                as.numeric(par$mu),as.numeric(par$sigma),as.numeric(par$alpha),
                as.integer(log), PACKAGE = "bamlss")
          })
      },
      "p" = function(y,par,...) {
        sapply(y, FUN = function(x) {
          .Call("bamlss_glogis_distr",as.numeric(x),
                as.numeric(par$mu),as.numeric(par$sigma),as.numeric(par$alpha))
        })
      },
      "q" = function(p,par,...) {
        sapply(p, FUN = function(x) {
          .Call("bamlss_glogis_quantile",as.numeric(p),
                as.numeric(par$mu),as.numeric(par$sigma),as.numeric(par$alpha))
        })
      },
      "r" = function(y,par,...) {
         glogis::rglogis(y,par$mu,par$sigma,par$alpha)
      },
      "initialize" = list(
         "mu"    = function(y, ...) { (y + mean(y)) / 2 },
         "sigma" = function(y, ...) { rep(sd(y), length(y)) },
         "alpha" = function(y, ...) { rep(1,length(y)) }
      ),
      "mean"      = function(par) as.vector(par$mu + (digamma(par$alpha) - digamma(1)) * par$sigma),
      "variance"  = function(par) as.vector((psigamma(par$alpha, deriv = 1) + psigamma(1, deriv = 1)) * par$sigma),
      "skewness"  = function(par) as.vector((psigamma(par$alpha, deriv = 2) - psigamma(1, deriv = 2)) /
                                              (psigamma(par$alpha, deriv = 1) + psigamma(1, deriv = 1))^(3/2))
   )

   # Return family object
   class(rval) <- "family.bamlss"
   return(rval)
}


## Most likely transformations.
mlt_bamlss <- function(todistr = "Normal")
{
  rval <- list(
    "family" = "mlt",
    "names" = "mu",
    "links" = c("mu" = "identity"),
    "optimizer" = opt_mlt,
    "sampler" = FALSE
  )
  rval$distr <- mlt_distr(todistr)
  class(rval) <- "family.bamlss"
  rval
}

mlt_Normal <- function() {
    list(p = pnorm, d = dnorm, q = qnorm,
         ### see also MiscTools::ddnorm
         dd = function(x) -dnorm(x = x) * x,
         ddd = function(x) dnorm(x = x) * (x^2 - 1),
         name = "normal")
}

mlt_Logistic <- function() {
    list(p = plogis, d = dlogis, q = qlogis,
         dd = function(x) {
             ex <- exp(x)
             (ex - exp(2 * x)) / (1 + ex)^3
         },
         ddd = function(x) {
             ex <- exp(x)
             (ex - 4*(exp(2 * x)) + exp(3 * x)) / (1 + ex)^4
         },
         name = "logistic")
}

mlt_MinExtrVal <- function() {
    list(p = function(x) 1 - exp(-exp(x)),
         q = function(p) log(-log(1 - p)),
         d = function(x, log = FALSE) {
             ret <- x - exp(x)
             if (!log) return(exp(ret))
             ret
         },
         dd = function(x) {
             ex <- exp(x)
             (ex - exp(2 * x)) / exp(ex)
         },
         ddd = function(x) {
             ex <- exp(x)
             (ex - 3*exp(2 * x) + exp(3 * x)) / exp(ex)
         },
         name = "minimum extreme value")
}

mlt_distr <- function(which = c("Normal", "Logistic", "MinExtrVal")) {
    which <- match.arg(which)
    do.call(paste("mlt_", which, sep = ""), list())
}


## Family object for the nested multinomial model.
pgev <- function(x, xi = 1)
{
  if(any(i <- xi < 0)) {
    tx <- x[i]
    txi <- xi[i]
    if(any(j <- tx > -1/xi[i])) {
      tx[j] <- -1/xi[i][j] - 0.01
      x[i] <- tx
    }
  }
  if(any(i <- xi > 0)) {
    tx <- x[i]
    txi <- xi[i]
    if(any(j <- tx < -1/xi[i])) {
      tx[j] <- -1/xi[i][j] + 0.01
      x[i] <- tx
    }
  }
  i <- xi == 0
  p <- x
  p[i] <- exp(-exp(-x[i]))
  p[!i] <- exp(-pmax(1 + xi[!i] * x[!i], 0)^(-1/xi[!i]))
  if(any(i <- p == 0))
    p[i] <- 0.0001
  p
}

pgev_dw <- function(x, xi = 1)
{
  if(any(i <- xi < 0)) {
    tx <- x[i]
    txi <- xi[i]
    if(any(j <- tx > -1/xi[i])) {
      tx[j] <- -1/xi[i][j] - 0.01
      x[i] <- tx
    }
  }
  if(any(i <- xi > 0)) {
    tx <- x[i]
    txi <- xi[i]
    if(any(j <- tx < -1/xi[i])) {
      tx[j] <- -1/xi[i][j] + 0.01
      x[i] <- tx
    }
  }
  i <- xi == 0
  p <- x
  p[i] <- exp(-x[i])
  p[!i] <- pmax(1 + xi[!i] * x[!i], 0)^(-1/xi[!i] - 1)
  p
}

pgev_dxi <- function(x, xi = 1)
{
  if(any(i <- xi < 0)) {
    tx <- x[i]
    txi <- xi[i]
    if(any(j <- tx > -1/xi[i])) {
      tx[j] <- -1/xi[i][j] - 0.01
      x[i] <- tx
    }
  }
  if(any(i <- xi > 0)) {
    tx <- x[i]
    txi <- xi[i]
    if(any(j <- tx < -1/xi[i])) {
      tx[j] <- -1/xi[i][j] + 0.01
      x[i] <- tx
    }
  }
  i <- xi != 0
  p <- x
  p[!i] <- 0
  pm <- pmax(1 + xi[i] * x[i], 0)
  p[i] <- pm^(-1/xi[i]) * log(pm) * 1/(xi[i]^2) + pm^(-1/xi[i] - 1) * x[i]/xi[i]
  p
}

nmult_bamlss <- function(K)
{
  links <- c(rep("identity", K - 1L), rep("identity", K), rep("identity", K))
  names(links) <- c(paste0("alpha", 1L:(K - 1L)), paste0("w", 1L:K), paste0("xi", 1L:K))

  h <- make.link("logit")$linkinv

  rval <- list(
    "family" = "Nested Multinomial",
    "names"  = names(links),
    "links" = links,
    "d" = function(y, par, log = FALSE) {
      if(is.null(dim(y))) {
        y <- model.matrix(~ y - 1)
      }

      alpha <- exp(do.call("cbind", par[grep("alpha", names(par))]))
      w <- do.call("cbind", par[grep("w", names(par))])
      xi <- do.call("cbind", par[grep("xi", names(par))])

      probs_alpha <- alpha / (rowSums(alpha) + 1)
      probs_alpha <- cbind(probs_alpha, 1 - rowSums(probs_alpha))

      probs_w <- 1 - pgev(w, xi)

      d <- NULL
      for(j in 1:K) {
        d <- cbind(d, (1 - probs_w[, j]) * probs_alpha[, j])
        d <- cbind(d, probs_w[, j] * probs_alpha[, j])
      }

      d[d < 1e-10] <- 1e-10

      d <- rowSums(y * log(d))
      if(!log)
        d <- exp(d)

      return(d)
    }
  )

  score_alpha <- function(y, par, id, ...) {
    if(is.null(dim(y))) {
      y <- model.matrix(~ y - 1)
    }
    i <- as.integer(gsub("alpha", "", id))
    ind <- sort(rep(1:K, 2))
    j <- which(ind == i)
    alpha <- exp(do.call("cbind", par[grep("alpha", names(par))]))
    pi <- alpha[, id] / (rowSums(alpha) + 1)
    rowSums(y[, j]) - pi
  }

  salpha <- rep(list(score_alpha), K - 1L)
  names(salpha) <- paste0("alpha", K - 1L)

  rval$score <- salpha

#  score_w <- function(y, par, id, ...) {
#    if(is.null(dim(y))) {
#      y <- model.matrix(~ y - 1)
#    }
#    i <- as.integer(gsub("w", "", id))
#    ind <- sort(rep(1:K, 2))
#    j <- which(ind == i)
#    j <- j[2]
#    y[, j] * pgev_dw(par[[id]], par[[paste0("xi", i)]]) 
#  }

#  sw <- rep(list(score_w), K - 1L)
#  names(sw) <- paste0("w", K - 1L)

#  rval$score <- c(rval$score, sw)

#  score_xi <- function(y, par, id, ...) {
#    if(is.null(dim(y))) {
#      y <- model.matrix(~ y - 1)
#    }
#    i <- as.integer(gsub("xi", "", id))
#    ind <- sort(rep(1:K, 2))
#    j <- which(ind == i)
#    j <- j[2]
#    y[, j] * pgev_dxi(par[[paste0("w", i)]], par[[id]]) 
#  }

#  sxi <- rep(list(score_xi), K - 1L)
#  names(sxi) <- paste0("xi", K - 1L)

#  rval$score <- c(rval$score, sxi)

  hess_alpha <- function(y, par, id, ...) {
    alpha <- exp(do.call("cbind", par[grep("alpha", names(par))]))
    pi <- alpha[, id] / (rowSums(alpha) + 1)
    pi * (1 - pi)
  }

  halpha <- rep(list(hess_alpha), K - 1L)
  names(halpha) <- paste0("alpha", K - 1L)

  rval$hess <- halpha

  rval$initialize <- list()
  init_xi <- rep(list(function(y, ...) {
    rep(1, if(is.null(dim(y))) length(y) else nrow(y))
  }), K)
  names(init_xi) <- paste0("xi", 1:K)
  rval$initialize <- c(rval$initialize, init_xi)

  rval$nocat <- TRUE

  rval$probabilities <- function(par, y = NULL, numeric = TRUE, ...) {
    if(!is.null(y)) {
      if(!is.character(y)) {
        if(is.null(dim(y))) {
          y <- model.matrix(~ y - 1)
        }
        y <- colnames(y)
      }
    } else {
      y <- paste0("p", 1:(2 * K))
    }

    alpha <- exp(do.call("cbind", par[grep("alpha", names(par))]))
    w <- do.call("cbind", par[grep("w", names(par))])
    xi <- do.call("cbind", par[grep("xi", names(par))])

    probs_alpha <- alpha / (rowSums(alpha) + 1)
    probs_alpha <- cbind(probs_alpha, 1 - rowSums(probs_alpha))

    probs_w <- 1 - pgev(w, xi)

    d <- NULL
    for(j in 1:K) {
      d <- cbind(d, (1 - probs_w[, j]) * probs_alpha[, j])
      d <- cbind(d, probs_w[, j] * probs_alpha[, j])
    }

    if(!numeric) {
      d <- t(apply(d, 1, function(x) {
        y <- rep(0,length(x))
        y[which.max(x)] <- 1
        y
      }))
    }

    colnames(d) <- y

    return(as.data.frame(d))
  }

  ## Return family object
  class(rval) <- "family.bamlss"
  return(rval)
}

#library(numDeriv)

#foo <- function(xi, w = 0.5) {
#  xi[abs(xi) < 0.00001] <- 0.00001
#  xi[xi <= -0.5] <- -0.4999
#  a <- (1 - xi * w)
#  a[a <= 0] <- 1e-05
#  1 - exp(-a^(-1/xi))
#}

#bar <- function(xi, w = 0.5) {
#  xi[abs(xi) < 0.00001] <- 0.00001
#  xi[xi <= -0.5] <- -0.4999
#  a <- (1 - xi * w)
#  a[a <= 0] <- 1e-05
#  exp(-a^(-1/xi)) * (a^(-1/xi) * (log(a) * (1/xi^2)) - a^((-1/xi) - 1) * ((-1/xi) * w))
#}

#curve(foo, -0.5, 10)
#numDeriv::grad(foo, 0.2)
#bar(0.2)

#if(FALSE) {
#foo1 <- function(x, p = 0.7) {
#  log(p * (1 - p)) - x * (p - 1 * (x < 0))
#}

#absx <- function(x) { x * tanh(x/0.1) }

#foo2 <- function(x, p = 0.7) {
#  i <- x > 0
#  rval <- rep(0, length(x))
#  rval[i] <- p * absx(x[i])
#  rval[!i] <- (1 - p) * absx(x[!i])
#  log(p * (1 - p)) - rval
#}

#curve(foo1, -1, 1)
#curve(foo2, -1, 1, col = 2, add = TRUE)

#library(numDeriv)

#numDeriv::grad(foo3, -0.01)
#foo3prime(-0.01)

#n <- 300
#x <- runif(n, -3, 3)
#y <- sin(x) + rnorm(n, sd = exp(-1 + 0.5 * x))

#data("mcycle", package = "MASS")
#x <- mcycle$times
#y <- mcycle$accel / max(abs(mcycle$accel))

#f <- list(
#  y ~ s(x),
#  sigma ~ s(x),
#  lambda ~ s(x)
#)

#maxit = 500

#b1 <- bamlss(f, family = gF("ELF", tau = 0.1), optimizer = boost, initialize = FALSE)
#b2 <- bamlss(f, family = gF("ELF", tau = 0.5), optimizer = boost, initialize = FALSE)
#b3 <- bamlss(f, family = gF("ELF", tau = 0.9), optimizer = boost, initialize = FALSE)

#p1 <- predict(b1, model = "mu")
#p2 <- predict(b2, model = "mu")
#p3 <- predict(b3, model = "mu")

#plot(x, y)
#plot2d(cbind(p1, p2, p3) ~ x, col.lines = 4, lty = c(2, 1, 2), add = TRUE)
#}

ALD_bamlss <- function(..., tau = 0.5, eps = 0.01)
{
  if(tau < 0.001 | tau > (1 - 0.001))
    stop("tau must be between 0 and 1!")

  absx <- function(x) { sqrt(x^2 + eps) }

  rval <- list(
    "family" = "Asymmetric Laplace",
    "names" = c("mu", "sigma"),
    "links" = c(mu = "identity", sigma = "log"),
    "d" = function(y, par, log = FALSE) {
      ## d <- logpc - r * (p - 1 * (r < 0))
      ## e1 <- expression(  -2 * log(s) - tau * sqrt((y-mu)^2)/s^2  )
      ## e2 <- expression(  -2 * log(s) - (1-tau) * sqrt((y-mu)^2)/s^2  )
      ## Deriv(e1, cache.exp = FALSE, "mu")
      r <- y - par$mu
      i <- r > 0
      d <- rep(0, length(r))
      d[i] <- tau * absx(r[i]) / par$sigma[i]^2
      d[!i] <- (1 - tau) * absx(r[!i]) / par$sigma[!i]^2
      d <- log(tau * (1 - tau)) - 2 * log(par$sigma) - d
      if(!log)
        d <- exp(d)
      return(d)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        r <- (y - par$mu)
        i <- r > 0
        score <- r/(par$sigma^2 * absx(r))
        score[i] <- tau * score[i]
        score[!i] <- (1 - tau) * score[!i]
        return(score)
      },
      "sigma" = function(y, par, ...) {
        r <- (y - par$mu)
        i <- r > 0
        score <- rep(0, length(r))
        score[i] <- (2 * (tau * absx(r[i])/par$sigma[i]^2) - 2)/par$sigma[i]
        score[!i] <- (2 * ((1 - tau) * absx(r[!i])/par$sigma[!i]^2) - 2)/par$sigma[!i]
        return(score)
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        r <- (y - par$mu)
        i <- r > 0
        hess <- rep(0, length(r))
        hess[i] <- tau * (par$sigma[i]^2 * sign(r[i]) * r[i]/(par$sigma[i]^2 * absx(r[i]))^2 -
          1/(par$sigma[i]^2 * absx(r[i])))
        hess[!i] <- (1 - tau) * (par$sigma[!i]^2 * sign(r[!i]) * (r[!i])/(par$sigma[!i]^2 * absx(r[!i]))^2 -
          1/(par$sigma[!i]^2 * absx(r[!i])))
        return(-hess)
      },
      "sigma" = function(y, par, ...) {
        r <- (y - par$mu)
        i <- r > 0
        hess <- rep(0, length(r))
        hess[i] <- -((6 * (tau * absx(r[i])/par$sigma[i]^2) - 2)/par$sigma[i]^2)
        hess[!i] <- -((6 * ((1 - tau) * absx(r[!i])/par$sigma[!i]^2) - 2)/par$sigma[!i]^2)
        return(-hess)
      }
    ),
    "initialize" = list(
      "mu"    = function(y, ...) { (y + quantile(y, prob = tau)) / 2 },
      "sigma" = function(y, ...) { rep(sd(y), length(y)) }
    )
  )

  class(rval) <- "family.bamlss"
  rval
}


log1pexp <- function(x)
{
  indx <- .bincode(x, c(-Inf, -37, 18, 33.3, Inf), T)
  
  kk <- which(indx==1)
  if( length(kk) ){  x[kk] <- exp(x[kk])  }
  
  kk <- which(indx==2)
  if( length(kk) ){  x[kk] <- log1p( exp(x[kk]) ) }
  
  kk <- which(indx==3)
  if( length(kk) ){  x[kk] <- x[kk] + exp(-x[kk]) }
  
  return(x)
}

ELF_bamlss <- function(..., tau = 0.5)
{
  if(tau < 0.001 | tau > (1 - 0.001))
    stop("tau must be between 0 and 1!")

  rval <- list(
    "family" = "Elf density",
    "names" = c("mu", "sigma", "lambda"),
    "links" = c(mu = "identity", sigma = "log", lambda = "log"),
    "d" = function(y, par, log = FALSE) {
      r <- y - par$mu
#      d <- (1 - tau) * r/par$sigma - par$lambda * log(1 + exp(r/(par$lambda * par$sigma))) -
#        log(par$lambda * par$sigma * beta(par$lambda *(1 - tau), par$lambda*tau))

## library("Deriv")
## e <- expression((1 - tau) * (y-mu)/s - lam * log(1 + exp((y-mu)/(lam*s))) - log(lam * s * beta(lam*(1 - tau), lam*tau)))
## Deriv(e, cache.exp = FALSE, "mu")

      d <- (1 - tau) * r/par$sigma - par$lambda * log1pexp(r/(par$lambda * par$sigma)) -
        log(par$lambda * par$sigma * beta(par$lambda *(1 - tau), par$lambda*tau))

      if(!log)
        d <- exp(d)

      return(d)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        .e2 <- exp((y - par$mu)/(par$lambda * par$sigma))
        score <- (.e2/(1 + .e2) + tau - 1)/par$sigma
        return(score)
      },
      "sigma" = function(y, par, ...) {
        .e1 <- par$lambda * par$sigma
        .e2 <- y - par$mu
        .e4 <- exp(.e2/.e1)
        score <- (par$lambda^2 * .e4/((1 + .e4) * .e1^2) - (1 - tau)/par$sigma^2) * .e2 - 1/par$sigma
        return(score)
      },
      "lambda" = function(y, par, ...) {
        .e1 <- par$lambda * par$sigma
        .e2 <- y - par$mu
        .e4 <- exp(.e2/.e1)
        .e5 <- 1 - tau
        .e6 <- digamma(par$lambda)
        score <- -((1 + par$lambda * (.e5 * (digamma(par$lambda * .e5) - .e6) + tau * (digamma(par$lambda * tau) - .e6)))/par$lambda + log1p(.e4) - .e1 * .e4 * .e2/((1 + .e4) * .e1^2))
        return(score)
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        .e2 <- exp((y - par$mu)/(par$lambda * par$sigma))
        .e3 <- 1 + .e2
        hess <- .e2 * (.e2/.e3 - 1)/(par$lambda * par$sigma^2 * .e3)
        return(-hess)
      },
      "sigma" = function(y, par, ...) {
        .e1 <- par$lambda * par$sigma
        .e2 <- y - par$mu
        .e4 <- exp(.e2/.e1)
        .e5 <- 1 + .e4
        hess <- (2 * ((1 - tau)/par$sigma^3) - par$lambda^3 * ((2 * (.e1 * .e5) - .e4 * .e2)/(.e5 * .e1^2)^2 + .e2/(.e5 * .e1^4)) * .e4) * .e2 + 1/par$sigma^2
        return(-hess)
      },
      "lambda" = function(y, par, ...) {
        .e1 <- par$lambda * par$sigma
        .e2 <- y - par$mu
        .e3 <- 1 - tau
        .e5 <- exp(.e2/.e1)
        .e6 <- digamma(par$lambda)
        .e7 <- .e1^2
        .e8 <- 1 + .e5
        .e9 <- par$lambda * .e3
        .e10 <- par$lambda * tau
        .e11 <- .e3 * (digamma(.e9) - .e6)
        .e12 <- .e8 * .e7
        .e15 <- tau * (digamma(.e10) - .e6)
        .e16 <- trigamma(par$lambda)
        hess <- -((.e11 + par$lambda * ((.e3 * trigamma(.e9) - .e16) * .e3 + tau * 
          (tau * trigamma(.e10) - .e16)) + .e15 - (1 + par$lambda * (.e11 + .e15))/par$lambda)/par$lambda -
          par$sigma * ((2 - .e1 * .e2/.e7)/.e12 - .e1 * (2 * (.e1 * .e8) - .e5 * .e2)/.e12^2) * .e5 * .e2)
        return(-hess)
      }
    ),
    "initialize" = list(
      "mu"    = function(y, ...) { (y + quantile(y, prob = tau)) / 2 },
      "sigma"    = function(y, ...) { sd((y - quantile(y, prob = tau)) / 2) },
      "lambda"    = function(y, ...) { rep(4, length(y)) }
    )
  )

  class(rval) <- "family.bamlss"
  rval
}

nbinom_bamlss <- function(...) {
### neg bin
### Author: Thorsten Simon
### Date:   2019 Feb
  rval <- list(
    "family" = "nbinom",
    "names" = c("mu", "theta"),
    "links" = c(mu = "log", theta = "log"),
    "mu" = function(par, ...) { par$mu },
    "loglik" = function(y, par, ...) {
      sum(stats::dnbinom(y, mu = par$mu, size = par$theta, log = TRUE))
    },
    "d" = function(y, par, log = FALSE) {
      stats::dnbinom(y, mu = par$mu, size = par$theta, log = log)
    },
    "p" = function(y, par, ...) {
      stats::pnbinom(y, mu = par$mu, size = par$theta)
    },
    "q" = function(p, par, ...) {
      stats::qnbinom(p, mu = par$mu, size = par$theta)
    },
    "r" = function(n, par, ...) {
      stats::rnbinom(n, mu = par$mu, size = par$theta)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        y - (y + par$theta)*par$mu/(par$mu + par$theta)
      },
      "theta" = function(y, par, ...) {
        (digamma(y + par$theta) - digamma(par$theta) + log(par$theta) + 1 -
          log(par$mu + par$theta) - (y + par$theta)/(par$mu + par$theta)) *
          par$theta
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        par$mu * par$theta / (par$mu + par$theta)
      },
      "theta" = function(y, par, ...) {
        (digamma(y + par$theta) - digamma(par$theta) + log(par$theta) + 1 -
          log(par$mu + par$theta) - (y + par$theta)/(par$mu + par$theta)) *
          -par$theta - (trigamma(y + par$theta) - trigamma(par$theta) +
          1/par$theta - 2/(par$mu + par$theta) + (y + par$theta) /
          (par$mu + par$theta)^2) * par$theta^2
      }
    ),
    "initialize" = list(
      "mu" = function(y, ...) { (y + mean(y)) / 2 },
      "theta" = function(y, ...) {
        rep( min( mean(y)^2 / (var(y) - mean(y)), 10), length(y))
      }
    )
  )

#  dnbinom2 <- function(x, mu, size, log = FALSE) {
#    prob <- size / (size + mu)

#    d <- lgamma(x + size) - lgamma(size) - lfactorial(x) +
#      size * log(prob) + x * log(1-prob)
#    if(!log)
#      d <- exp(d)

#    d
#  }

#  dnbinom(1:10, mu = 3, size=2)
#  dnbinom2(1:10, mu = 3, size=2)


  rval$keras <- list(
    "nloglik" = function(y_true, y_pred) {
      K = keras::backend()
      tf = tensorflow::tf

      mu = K$exp(y_pred[, 1])
      size = K$exp(y_pred[, 2])
      x = y_true[, 1]

      prob = size / (size + mu)

      ll = tf$math$lgamma(x + size) - tf$math$lgamma(size) -
        tf$math$lgamma(x + 1) + size * K$log(prob) + x * K$log(1 - prob)
      ll = K$sum(ll)

      return(-1 * ll)
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


gpareto2_bamlss <- function(...)
{
  rval <- list(
    "family" = "gpareto",
    "names" = c("xi", "sigma", "lambda", "gamma"),
    "links" = c(xi = "log", sigma = "log", lambda = "log", gamma = "logit"),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x >= 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "d" = function(y, par, log = FALSE) {
      yt <- (y + par$lambda - 1)^(par$gamma)
      d <- -log(par$sigma) - (1 / par$xi + 1) * log(1 + par$xi * yt / par$sigma)
      if(!log)
        d <- exp(d)
      d
    }
  )

  class(rval) <- "family.bamlss"
  rval
}


#ipp_bamlss <- function(...)
#{
#  rval <- list(
#    "family" = "IPP",
#    "names" = "lambda",
#    "links" = c(lambda = "identity"),
#    "d" = function(y, par, ...) {
#      f <- function(x,y) {
#        exp(par$lambda)
#      }
#      intL <- rmutil::int2(f, a=c(0,0), b=c(1,1))
#      return(par$lambda - intL)
#    }
#  )

#  class(rval) <- "family.bamlss"
#  rval
#}


## logNN model.
logNN_bamlss <- function(...)
{
  stopifnot(requireNamespace("statmod"))

  N <- list(...)$N
  if(is.null(N))
    N <- 100
  check <- list(...)$check
  if(is.null(check))
    check <- FALSE

  gq <- statmod::gauss.quad(N, kind = "hermite")

  gq$weights <- gq$weights * exp(gq$nodes^2)

  A <- function(t, y, mu, sigma, lambda) {
    exp(-(((t - mu) / sigma)^2 + ((y - exp(t))/lambda)^2 )/2 ) / (2*pi*sigma*lambda)
  }

  rval <- list(
    "family" = "logNormal-Normal Convolution",
    "names" = c("mu", "sigma", "lambda"),
    "links" = c("identity", "log", "log"),
    "d" = function(y, par, log = FALSE, ...) {
      if(check) {
        d <- .Call("logNN_dens", gq$nodes, gq$weights, y, par$mu,
          par$sigma, par$lambda, package = "bamlss")
        foo <- function(t, y, mu, sigma, lambda) {
          A(t, y, mu, sigma, lambda)
        }
        n <- length(y)
        di <- rep(0, n)
        for(i in 1:n) {
          di[i] <- integrate(foo, -Inf, Inf, y = y[i], mu = par$mu[i],
            sigma = par$sigma[i], lambda = par$lambda[i],
            rel.tol = .Machine$double.eps^.75, subdivisions = 500L)$value
        }
        plot(d, di,
          main = paste("Mean absolute difference:", round(mean(abs(d - di)), 4)),
          xlab = "Gauss quadrature", ylab = "integrate()", ...)
        abline(0, 1)
      } else {
        d <- .Call("logNN_dens", gq$nodes, gq$weights, y, par$mu,
          par$sigma, par$lambda, package = "bamlss")
      }
      if(log)
        d <- log(d)
      return(d)
    },
    "p" = function(y, par, lower.tail = TRUE, log.p = FALSE, ...) {
      foo <- function(t, y, mu, sigma, lambda) {
        dnorm(t, mean = mu, sd = sigma) * pnorm(y - exp(t), sd = lambda)
      }
      n <- length(y)
      p <- rep(0, n)
      for(i in 1:n) {
        p[i] <- integrate(foo, -Inf, Inf,
          y = y[i], mu = par$mu[i], sigma = par$sigma[i], lambda = par$lambda[i],
          rel.tol = .Machine$double.eps^.75)$value
      }
      if(!lower.tail)
        p <- 1 - p
     if(log.p)
       return(log(p))
     else
      return(p)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
#        d <- .Call("logNN_dens", ba2, nodes, gq$weights, y, par$mu,
#          par$sigma, par$lambda, package = "bamlss")
#        foo <- function(t, y, mu, sigma, lambda) {
#          A(t, y, mu, sigma, lambda) * (t - mu)
#        }
#        n <- length(y)
#        rval <- rep(0, n)
#        for(i in 1:n) {
#          fx <- foo(nodes, y = y[i], mu = par$mu[i],
#            sigma = par$sigma[i], lambda = par$lambda[i])
#          rval[i] <- sum(gq$weights * fx) * ba2
#        }
#        rval <- 1/d * rval/par$sigma^2
        rval <- .Call("logNN_score_mu", gq$nodes, gq$weights, y, par$mu,
          par$sigma, par$lambda, package = "bamlss")
        return(rval)
      },
      "sigma" = function(y, par, ...) {
#        d <- .Call("logNN_dens", ba2, nodes, gq$weights, y, par$mu,
#          par$sigma, par$lambda, package = "bamlss")
#        foo <- function(t, y, mu, sigma, lambda) {
#          A(t, y, mu, sigma, lambda) * ((t - mu)^2 - sigma^2)
#        }
#        n <- length(y)
#        rval <- rep(0, n)
#        for(i in 1:n) {
#          fx <- foo(nodes, y = y[i], mu = par$mu[i],
#            sigma = par$sigma[i], lambda = par$lambda[i])
#          rval[i] <- sum(gq$weights * fx) * ba2
#        }
#        rval <- 1/d * rval/par$sigma^2
        rval <- .Call("logNN_score_sigma", gq$nodes, gq$weights, y, par$mu,
          par$sigma, par$lambda, package = "bamlss")
        rval
      },
      "lambda" = function(y, par, ...) {
#        d <- .Call("logNN_dens", ba2, nodes, gq$weights, y, par$mu,
#          par$sigma, par$lambda, package = "bamlss")
#        foo <- function(t, y, mu, sigma, lambda) {
#          A(t, y, mu, sigma, lambda) * ((y - exp(t))^2 - lambda^2)
#        }
#        n <- length(y)
#        rval <- rep(0, n)
#        for(i in 1:n) {
#          fx <- foo(nodes, y = y[i], mu = par$mu[i],
#            sigma = par$sigma[i], lambda = par$lambda[i])
#          rval[i] <- sum(gq$weights * fx) * ba2
#        }
#        rval <- 1/d * rval/par$lambda^2
        rval <- .Call("logNN_score_lambda", gq$nodes, gq$weights, y, par$mu,
          par$sigma, par$lambda, package = "bamlss")
        rval
      }
    ),
    "initialize" = list(
      "mu"    = function(y, ...) {
        rep(mean(log(y[y > 0])), length(y))
      },
      "sigma" = function(y, ...) {
        rep(sd(log(y[y > 0])), length(y))
      },
      "lambda" = function(y, ...) {
        e <- y - exp(mean(log(y[y > 0])))
        rep(sd(e), length(y))
      }
    )
  )

  class(rval) <- "family.bamlss"
  rval
}

if(FALSE) {
  ## Simulate logNN data.
  set.seed(123)

  n <- 3000

  ## Round observations for using the binning option.
  x <- round(runif(n, -3, 3), 2)
  mu <- 0.1 + sin(x)
  Z <- rlnorm(n, mu, sd = exp(-2 + cos(x)))
  y <- Z + rnorm(n, 0, 0.3)

  f <- list(
    y ~ s(x),
    sigma ~ s(x),
    lambda ~ 1
  )

  b <- bamlss(f, family = logNN_bamlss(N=200,check=FALSE), binning = TRUE, eps = 0.01)

  mu <- as.matrix(predict(b, model = "mu", FUN = identity))
  sigma <- as.matrix(predict(b, model = "sigma", type = "parameter", FUN = identity))
  fit <- t(apply(exp(mu + sigma^2/2), 1, FUN = c95))
  fsigma <- predict(b, model = "sigma", FUN = c95)

  par(mfrow = c(1, 3))
  plot(y ~ x, ylim = range(c(fit, y)), col = rgb(0.1, 0.1, 0.1, alpha = 0.3),
    main = "Data and fit")
  plot2d(fit ~ x, add = TRUE, col.lines = 4, lwd = 3)
  fmu <- t(apply(mu, 1, c95))
  plot2d(fmu ~ x, main = "True mu and fit")
  plot2d(I(0.1 + sin(x)) ~ x, col.lines = 2, add = TRUE)
  plot2d(fsigma ~ x, main = "True log(sigma) and fit")
  plot2d(I(-2 + cos(x)) ~ x, col.lines = 2, add = TRUE)
}

pgev <- function(q, loc = 0, scale = 1, shape = 0, log = FALSE) 
{
  q <- (q - loc)/scale
  rval <- rep(0, length(q))
  i <- shape == 0
  if(log) {
    rval[i] <- -exp(-q[i])
    rval[!i] <- -pmax(1 + shape[!i] * q[!i], 0)^(-1/shape[!i])
  } else {
    rval[i] <- exp(-exp(-q[i]))
    rval[!i] <- exp(-pmax(1 + shape[!i] * q[!i], 0)^(-1/shape[!i]))
  }
  return(rval)
}

pgev_bamlss <- function(...)
{
  rval <- list(
    "family" = "pgev",
    "names" = c("mu", "lambda"),
    "links" = c(mu = "identity", lambda = "identity"),
    "valid.response" = function(x) {
      if(length(unique(x)) > 2)
        stop("only binary responses!") 
    },
    "d" = function(y, par, log = FALSE) {
      p <- pgev(par$mu, loc = 0, scale = 1, shape = par$lambda)
      d <- y * log(p) + (1 - y) * log(1 - p)
      if(!log)
        d <- exp(d)
      d
    },
    "initialize" = list(
      "mu" = function(y, ...) {
        y <- process_factor_response(y)
        rep(log(-log(mean(y))), length = length(y))
      },
      "lambda" = function(y) { rep(0, length(y)) }
    )
  )

  rval$probabilities <- function(par, ...) {
    pgev(par$mu, loc = 0, scale = 1, shape = par$lambda)
  }

  class(rval) <- "family.bamlss"
  rval
}

eel <- function(x, mu = 0, sigma = 1, lambda = 0.5, alpha = 0.5) {
  (1 - (1 + exp((x - mu)/sigma))^(-lambda))^alpha
}

eel_bamlss <- function(...)
{
  rval <- list(
    "family" = "eel",
    "names" = c("mu", "lambda", "alpha"),
    "links" = c(mu = "identity", lambda = "identity", alpha = "identity"),
    "valid.response" = function(x) {
      if(length(unique(x)) > 2)
        stop("only binary responses!") 
    },
    "d" = function(y, par, log = FALSE) {
      p <- eel(par$mu, mu = 0, sigma = 1, lambda = par$lambda, alpha = par$alpha)
      p[p < 1e-10] <- 1e-10
      p[p > 0.9999999] <- 0.9999999
      d <- y * log(p) + (1 - y) * log(1 - p)
      if(!log)
        d <- exp(d)
      d
    },
    "initialize" = list(
      "mu" = function(y, ...) {
        y <- process_factor_response(y)
        rep(mean(y), length = length(y))
      },
      "lambda" = function(y) { rep(1, length(y)) },
      "alpha" = function(y) { rep(1, length(y)) }
    )
  )

  rval$probabilities <- function(par, ...) {
    eel(par$mu, mu = 0, sigma = 1, lambda = par$lambda, alpha = par$alpha)
  }

  class(rval) <- "family.bamlss"
  rval
}


Sbqr <- function(u, tau = 0.5, alpha = 0.05) {
  tau * u + alpha * log(1 + exp(-u / alpha))
}

bqr_bamlss <- function(...)
{
  tau <- list(...)$tau
  if(is.null(tau))
    tau <- 0.5

  rval <- list(
    "family" = "bqr",
    "names" = "mu",
    "links" = c(mu = "identity"),
    "valid.response" = function(x) {
      if(length(unique(x)) > 2)
        stop("only binary responses!") 
    },
    "d" = function(y, par, ...) {
      -1 * Sbqr(y - pnorm(par$mu), tau = tau)^2
    },
    "loglik" = function(y, par, ...) {
      -1 * sum(Sbqr(y - pnorm(par$mu), tau = tau)^2, na.rm = TRUE)
    },
    "initialize" = list(
      "mu" = function(y, ...) {
        y <- process_factor_response(y)
        rep(qnorm(tau), length = length(y))
      }
    )
  )

  rval$probabilities <- function(par, ...) {
    pnorm(par$mu)
  }

  class(rval) <- "family.bamlss"
  rval
}


dgp_bamlss <- function(...)
{
  fam <- gpareto_bamlss()

  rval <- list(
    "family" = "discrete generalized pareto",
    "names" = c("xi", "sigma"),
    "links" = c(xi = "log", sigma = "log"),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x >= 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "d" = function(y, par, log = FALSE, ...) {
      d <- fam$p(y + 1, par) - fam$p(y, par)
      if(log)
        d <- log(d)
      return(d)
    },
    "p" = function(y, par, log = FALSE, ...) {
      par <- as.data.frame(par)
      n <- length(y)
      p <- rep(0, n)
      for(i in 1:n) {
        dy <- fam$p((y[i] + 1):1, par[i, ]) - fam$p((y[i]):0, par[i, ])
        p[i] <- sum(dy)
      }
      return(p)
    }
  )
  rval$type <- "discrete"

  class(rval) <- "family.bamlss"
  rval
}


discretize <- function(family)
{
  family <- bamlss.family(family)
  rval <- list()
  rval$family <- paste("discrete", family$family)
  rval$names <- family$names
  rval$links <- family$links
  rval$valid.response <- function(x) {
    if(is.factor(x)) return(FALSE)
    if(ok <- !all(x >= 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
    ok
  }
  rval$d <- function(y, par, log = FALSE, ...) {
    d <- family$p(y + 1, par) - family$p(y, par)
    if(log)
      d <- log(d)
    return(d)
  }
  rval$p <- function(y, par, log = FALSE, ...) {
    par <- as.data.frame(par)
    n <- length(y)
    p <- rep(0, n)
    for(i in 1:n) {
      dy <- family$p((y[i] + 1):1, par[i, ]) - family$p((y[i]):0, par[i, ])
      p[i] <- sum(dy)
    }
    return(p)
  }
  rval$q <- function(p, par, ...) {
    par <- as.data.frame(par)
    n <- nrow(par)
    x <- rep(NA, n)
    p <- rep(p, length.out = n)
    for(i in 1:n) {
      alpha <- 0
      y <- 0
      while(alpha < p[i]) {
        alpha <- rval$p(y, par[i, ])
        y <- y + 1
      }
      x[i] <- y
    }
    return(x)
  }
  rval$initialize <- family$initialize
  rval$type <- "discrete"
  class(rval) <- "family.bamlss"
  rval
}


## https://www.jstor.org/stable/2288808?seq=1#metadata_info_tab_contents
## From gamlss.dist.
dSichel <- function(x, mu = 1, sigma = 1, nu = 1, log = FALSE, ...)
{
  alpha <- sqrt(1/sigma^2 + 2 * mu/sigma)

  a <- x * log(mu) + log(besselK(alpha, x + nu))
  b <- lfactorial(x) + (x + nu) * (log(alpha) + log(sigma)) + log(besselK(1/sigma, nu))
  d <- a - b

  if(!log)
    d <- exp(d)

  return(d)
}

if(FALSE) {
  x <- 0:100
  d <- dSichel(x)
  d2 <- dSICHEL(x, mu = 1, sigma = 1, nu = 1)
  plot(d ~ x, type = "h")
  print(sum(d))
}

Sichel_bamlss <- function(...)
{
  rval <- list(
    "family" = "Sichel",
    "names" = c("mu", "sigma", "nu"),
    "links" = c(mu = "log", sigma = "log", nu = "identity"),
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x >= 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "d" = function(y, par, log = FALSE) {
      dSichel(y, mu = par$mu, sigma = par$sigma, nu = par$nu, log = log)
    }
  )

  rval$p <- function(y, par, log = FALSE, ...) {
    par <- as.data.frame(par)
    n <- length(y)
    p <- rep(0, n)
    for(i in 1:n) {
      dy <- dSichel(0:y[i], mu = par$mu[i], sigma = par$sigma[i], nu = par$nu[i])
      p[i] <- sum(dy)
    }
    return(p)
  }

  rval$rps <- function(y, par, ymin = 1L, ymax = max(max(y), 100L)) {
    K <- seq(ymin, ymax, by = 1L)
    rps <- rep(0, length(y))
    for (k in K) {
      P <- rval$p(k, par)
      O <- y <= k
      rps <- rps + (P - O)^2
    }
    return(rps)
  }

  rval$mean <- function(par, ...) {
    s1 <- 1/par$sigma
    k1 <- besselK(s1, par$nu + 1)
    k2 <- besselK(s1, par$nu)
    a <- k1/k2
    return(par$mu * a)
  }

  rval$type <- "discrete"

  class(rval) <- "family.bamlss"
  rval
}


exp2 <- function(x) {
  x[x > 20] <- 20
  x <- exp(x)
  x
}


GEV_bamlss <- function(...)
{
   ## log(1/s) - (1+xi*((y-m)/s))^(-1/xi) - (1/xi+1) * log(1 + xi*((y-m)/s))
   ## log(1/exp(s)) - (1+xi*((y-m)/exp(s)))^(-1/xi) - (1/xi+1) * log(1 + xi*((y-m)/exp(s)))

   eps <- 1e-05

   rval <- list(
    "family" = "GEV",
    "names" = c("mu", "sigma", "xi"),
    "links" = c(mu = "identity", sigma = "log", xi = "identity"),
    "d" = function(y, par, log = FALSE, ...) {
      par$xi[par$xi >= 0 & par$xi < eps] <- eps
      par$xi[par$xi < 0 & par$xi > -eps] <- -eps
      par$sigma[par$sigma < 0.001] <- 0.001
      x <- (y - par$mu) / par$sigma
      xx <- 1 + par$xi * x
      xx[xx < eps] <- eps
      d <- -log(par$sigma) - xx^(-1/par$xi) - (1 / par$xi + 1) * log(xx)
      if(!log)
        d <- exp2(d)
      return(d)
    },
    "p" = function(y, par, ...) {
      par$xi[par$xi >= 0 & par$xi < eps] <- eps
      par$xi[par$xi < 0 & par$xi > -eps] <- -eps
      par$sigma[par$sigma < 0.001] <- 0.001
      x <- (y - par$mu) / par$sigma
      xx <- 1 + par$xi * x
      xx[xx < eps] <- eps
      p <- exp2(-xx^(-1/par$xi))
      return(p)
    },
    "q" = function(p, par, ...) {
      par$xi[par$xi >= 0 & par$xi < eps] <- eps
      par$xi[par$xi < 0 & par$xi > -eps] <- -eps
      par$sigma[par$sigma < 0.001] <- 0.001
      par$mu + par$sigma / par$xi * ((-log(p))^(-par$xi) - 1)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        m <- par$mu
        s <- par$sigma
        xi <- par$xi

        xi[xi >= 0 & xi < eps] <- eps
        xi[xi < 0 & xi > -eps] <- -eps
        s[s < 0.001] <- 0.001

        yms <- (y - m)/s
        xi1 <- 1/xi
        xiyms <- 1 + xi * yms
        xiyms[xiyms < eps] <- eps
        s1 <- 1/s
        smu <- xiyms^(-xi1 - 1) * -xi1 * xi * s1 + (xi1 + 1) * xi * s1/xiyms

        return(smu)
      },
      "sigma" = function(y, par, ...) {
        m <- par$mu
        s <- par$sigma
        xi <- par$xi

        xi[xi >= 0 & xi < eps] <- eps
        xi[xi < 0 & xi > -eps] <- -eps
        s[s < 0.001] <- 0.001

        yms <- (y - m)/s
        xiyms <- 1 + xi * yms
        xiyms[xiyms < eps] <- eps
        xi1 <- 1/xi

        ssigma <- -(s/s^2/(1/s) - xiyms^(-xi1 - 1) * -xi1 * xi * (y - m) * s/s^2 -
          (xi1 + 1) * (xi * ((y - m) * s/s^2)/xiyms))

        return(ssigma)
      },
      "xi" = function(y, par, ...) {
        m <- par$mu
        s <- par$sigma
        xi <- par$xi

        xi[xi >= 0 & xi < eps] <- eps
        xi[xi < 0 & xi > -eps] <- -eps
        s[s < 0.001] <- 0.001

        xi1 <- 1/xi
        yms <- (y - m)/s
        xiyms <- 1 + xi * yms
        xiyms[xiyms < eps] <- eps
        lxiyms <- log(xiyms)

        sxi <- -(xiyms^((-xi1) - 1) * ((-xi1) * yms) +
          xiyms^(-xi1) * (lxiyms *  (xi1^2)) +
          ((xi1 + 1) * (yms/xiyms) -
          xi1^2 * lxiyms))

        return(sxi)
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        m <- par$mu
        s <- par$sigma
        xi <- par$xi

        xi[xi >= 0 & xi < eps] <- eps
        xi[xi < 0 & xi > -eps] <- -eps
        s[s < 0.001] <- 0.001

        yms <- (y - m)/s
        xiyms <- 1 + xi * yms
        xiyms[xiyms < eps] <- eps
        xi1 <- 1/xi
        s1 <- 1/s
        xis1 <- xi * s1

        hmu <- (xi1 + 1) * xis1^2/xiyms^2 - xiyms^(-xi1 - 1 - 1) * ((-xi1 - 1) * xis1) * -xi1 * xis1

        return(-hmu)
      },
      "sigma" = function(y, par, ...) {
        m <- par$mu
        s <- par$sigma
        xi <- par$xi

        xi[xi >= 0 & xi < eps] <- eps
        xi[xi < 0 & xi > -eps] <- -eps
        s[s < 0.001] <- 0.001

        ym <- y - m
        yms <- ym/s
        s2 <- s^2
        xiyms <- 1 + xi * yms
        xiyms[xiyms < eps] <- eps
        xi1 <- 1/xi
        yms2 <- ym * s
        yms2s2 <- yms2/s2

        hsigma <- -((s/s2 - s * (2 * s2)/s2^2)/(1/s) + 
          s/s2 * (s/s2)/(1/s)^2 - (xiyms^(-xi1 - 1) *
          (-xi1 * (xi * (yms2s2 - yms2 * (2 * (s2))/s2^2))) -
          xiyms^(-xi1 - 1 - 1) * ((-xi1 - 1) *
          (xi * yms2s2)) * (-xi1 * xi * yms2s2)) -
          (xi1 + 1) * (xi * (yms2s2 - yms2 * 
          (2 * s2)/s2^2)/xiyms + 
          xi * yms2s2 * (xi * yms2s2)/xiyms^2))

        return(-hsigma)
      },
      "xi" = function(y, par, ...) {
        m <- par$mu
        s <- par$sigma
        xi <- par$xi

        xi[xi >= 0 & xi < eps] <- eps
        xi[xi < 0 & xi > -eps] <- -eps
        s[s < 0.001] <- 0.001

        yms <- (y - m)/s
        xiyms <- 1 + xi * yms
        xiyms[xiyms < eps] <- eps
        xi1 <- 1 / xi
        a <- -xi1 - 1
        lxiyms <- log(xiyms)

        hxi <- -((xiyms^(a - 1) * (a * yms) +
          xiyms^a * (lxiyms * xi1^2)) * (-xi1 * yms) + 
          xiyms^a * (xi1^2 * yms) + 
          ((xiyms^a * (-xi1 * yms) +
          xiyms^(-xi1) * (lxiyms * xi1^2)) * (lxiyms * (xi1^2)) +
          xiyms^(-xi1) * (yms/xiyms * xi1^2 - lxiyms * (2 * xi/(xi^2)^2))) -
          ((xi1 + 1) * yms^2/xiyms^2 + xi1^2 * yms/xiyms +
          (xi1^2 * (yms/xiyms) - 2 * xi/(xi^2)^2 * lxiyms)))

        return(-hxi)
      }
    )
  )

  rval$initialize <- list(
    "mu" = function(y, ...) { (y + mean(y)) / 2 },
    "sigma" = function(y, ...) { rep(sd(y), length(y)) },
    "xi" = function(y, ...) { rep(1e-02, length(y)) }
  )

  rval$mean <- function(par, ...) {
    par$mu + par$sigma * (gamma(1 - par$xi) - 1) / par$xi
  }

  class(rval) <- "family.bamlss"
  rval
}

gamlss_distributions <- function(type = c("continuous", "discrete"))
{
  stopifnot(requireNamespace("gamlss.dist"))
  if(!("package:gamlss.dist" %in% search()))
    attachNamespace("gamlss.dist")
  funs <- ls("package:gamlss.dist")
  type <- match.arg(type)
  d <- list()
  tf <- tempfile()
  tf2 <- tempfile()
  warn <- getOption("warn")
  options("warn" = -1)
  for(j in seq_along(funs)) {
    fj0 <- get(funs[j])
    png(tf)
    capture.output(fj <- try(fj0(), silent = TRUE), file = tf2)
    dev.off()
    if(!inherits(fj, "try-error")) {
      if(inherits(fj, "gamlss.family")) {
        if(tolower(fj$type) == tolower(type)) {
          d[[funs[j]]] <- fj0
        }
      }
    }
  }
  unlink(tf)
  unlink(tf2)
  options("warn" = warn)
  return(d)
}

mix_bamlss <- function(f1, f2, ...)
{
  if(is.function(f1))
    f1 <- f1()
  if(is.function(f2))
    f2 <- f2()
  if(!inherits(f1, "gamlss.family"))
    stop("all families need to be of class 'gamlss.family'!")
  if(!inherits(f2, "gamlss.family"))
    stop("all families need to be of class 'gamlss.family'!")

  npar <- f1$nopar + f2$nopar + 1
  pn0 <- c("mu", "sigma", "nu", "tau")
  links <- c(f1[paste0(pn0[1:f1$nopar], ".link")],
    f2[paste0(pn0[1:f2$nopar], ".link")])
  links <- unlist(links)
  links <- c("logit", links)
  names(links) <- c("alpha", paste0("beta", 1:(npar - 1)))

  fn1 <- f1$family[1]
  fn2 <- f2$family[1]

  fn1 <- paste0("d", fn1, "(y")
  fn2 <- paste0("d", fn2, "(y")

  k <- 1
  for(j in 1:f1$nopar) {
    fn1 <- paste0(fn1, ",", pn0[j], "=par$beta", k)
    k <- k + 1
  }
  for(j in 1:f2$nopar) {
    fn2 <- paste0(fn2, ",", pn0[j], "=par$beta", k)
    k <- k + 1
  }
  fn1 <- paste0(fn1, ")")
  fn2 <- paste0(fn2, ")")

  dfun <- paste0("dy <- par$alpha*", fn1, "+(1-par$alpha)*", fn2)
  dfun <- c("function(y,par,log=FALSE,...) {", dfun, "if(log)\n  dy <- log(dy)", "return(dy)", "}")
  dfun <- paste0(dfun, collapse ="\n")
  dfun <- eval(parse(text = dfun)) 

  rval <- list(
    "family" = paste("mix", f1$family[1], f2$family[1]),
    "names" = names(links),
    "links" = links,
    "valid.response" = function(x) {
      if(is.factor(x)) return(FALSE)
      if(ok <- !all(x >= 0)) stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    },
    "d" = dfun,
    "p" = function(y, par, log = FALSE, ...) {
      par <- as.data.frame(par)
      n <- length(y)
      p <- rep(0, n)
      for(i in 1:n) {
        p[i] <- sum(dfun(0:y[i], par[i, ]), na.rm = TRUE)
      }
      p[p > 1] <- 1
      return(p)
    }
  )
  rval$type <- unique(tolower(c(f1$type, f2$type)))
  if(length(rval$type) > 1)
    stop("types of mixing distributions must be identical!")

  class(rval) <- "family.bamlss"
  rval
}

