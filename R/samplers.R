#MCMCpack <- function(x, y, family, start = NULL,
#  n.iter = 1200, burnin = 200, thin = 1, verbose = 100, ...)
#{
#  par <- make_par(x, type = 2)

#  post.samp <- MCMCmetrop1R(log_posterior, theta.init = par$par,
#    x = x, logfun = TRUE, V = attr(x, "hessian"),
#    mcmc = n.iter, burnin = burnin, thin = thin, verbose = verbose)

#  colnames(post.samp) <- names(par$par)

#  post.samp
#}


sam_GMCMC <- GMCMC <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
  n.iter = 1200, burnin = 200, thin = 1, verbose = TRUE, step = 20,
  propose = "iwlsC_gp", chains = NULL, ...)
{
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  np <- length(nx)

  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, propose = propose, ...)

  nobs <- nrow(y)
  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }

  if(!is.null(start))
    x <- set.starting.values(x, start)

  if(is.character(propose)) {
    propose <- if(grepl("GMCMC", propose, fixed = TRUE)) {
      get(propose)
    } else {
      get(paste("GMCMC", propose, sep = "_"))
    }
  }

  theta <- propose2 <- fitfun <- list()
  for(i in nx) {
    theta[[i]] <- propose2[[i]] <- fitfun[[i]] <- list()
    nt <- NULL
    x[[i]] <- x[[i]]$smooth.construct
    for(j in names(x[[i]])) {
      theta[[i]][[j]] <- x[[i]][[j]]$state$parameters
      attr(theta[[i]][[j]], "fitted.values") <- x[[i]][[j]]$fit.fun(x[[i]][[j]]$X, x[[i]][[j]]$state$parameters)
      attr(theta[[i]][[j]], "hess") <- x[[i]][[j]]$state$hessian
      x[[i]][[j]]$state$fitted.values <- NULL
      if(!inherits(x[[i]][[j]], "special")) {
        x[[i]][[j]]$XW <- t(x[[i]][[j]]$X)
        x[[i]][[j]]$XWX <- crossprod(x[[i]][[j]]$X)
      }
      x[[i]][[j]]$dmvnorm_log <- dmvnorm_log
      if(is.null(x[[i]][[j]]$fxsp))
        x[[i]][[j]]$fxsp <- FALSE
      nt <- c(nt, if(j == "model.matrix") "p" else paste("s", j, sep = "."))
      if(is.null(x[[i]][[j]]$propose)) {
        if(!is.null(x[[i]][[j]]$xt$propose))
          x[[i]][[j]]$propose <- x[[i]][[j]]$xt$propose
      }
      if(is.null(family$propose)) {
        propose2[[i]][[j]] <- if(is.null(x[[i]][[j]]$propose)) propose else x[[i]][[j]]$propose
      } else {
        propose2[[i]][[j]] <- family$propose[[i]]
      }
      fitfun[[i]][[j]] <- function(x, p) {
        attr(p, "fitted.values")
      }
      x[[i]][[j]]$penaltyFunction <- as.integer(sapply(x[[i]][[j]]$S, is.function))
#      if(!is.null(x[[i]][[j]]$sparse.setup$block.index)) {
#        propose2[[i]][[j]] <- GMCMC_iwls
#      }
    }
    names(theta[[i]]) <- names(propose2[[i]]) <- names(fitfun[[i]]) <- names(x[[i]]) <- nt
  }

  logLik <- function(eta) {
    family$loglik(y, family$map2par(eta))
  }

  zworking <- as.numeric(rep(0, length = nobs))
  resids <- as.numeric(rep(0, length = nobs))
  if(!is.null(weights))
    weights <- as.data.frame(weights)
  if(!is.null(offset))
    offset <- as.data.frame(offset)

  samps <- gmcmc(fun = family, theta = theta, fitfun = fitfun, data = x,
    propose = propose2, logLik = logLik, n.iter = n.iter, burnin = burnin, thin = thin,
    y = y, simplify = FALSE, zworking = zworking, resids = resids,
    cores = 1, chains = chains, combine = FALSE, weights = weights, offset = offset,
    verbose = verbose, step = step, ...)

  samps
}


gmcmc <- function(fun, theta, priors = NULL, propose = NULL,
  fitfun = NULL, logLik = NULL, data = NULL, attr.copy = "hess",
  n.iter = 12000, burnin = 2000, thin = 10, verbose = TRUE, step = 20,
  simplify = TRUE, chains = NULL, cores = NULL,
  combine = TRUE, sleep = 1, ...)
{
  if(!is.list(theta)) {
    theta <- list(theta)
    names(theta) <- names(formals(fun))[1]
  }
  if(is.null(names(theta)))
    names(theta) <- paste("theta[", 1:length(theta), "]", sep = "")
  ntheta <- names(theta)
  k <- length(theta)
  for(i in 1:k) {
    if(is.list(theta[[i]])) {
      if(is.null(names(theta[[i]])) & length(theta[[i]]))
        names(theta[[i]]) <- paste("term[", 1:length(theta[[i]]), "]", sep = "")
      if(length(theta[[i]])) {
        for(j in seq_along(theta[[i]])) {
          if(is.null(names(theta[[i]][[j]])) & length(theta[[i]][[j]]))
            names(theta[[i]][[j]]) <- paste("[", 1:length(theta[[i]][[j]]), "]", sep = "")
        }
      }
    }
  }
  class(theta) <- c("gmcmc.theta", "list")

  plot <- list(...)$plot
  if(is.null(plot))
    plot <- FALSE

  parse_input <- function(input = NULL, default, type) {
    ninput <- deparse(substitute(input), backtick = TRUE, width.cutoff = 500)
    if(is.null(input))
      input <- default
    if(!is.list(input))
      input <- rep(list(input), length.out = k)
    if(is.null(names(input)))
      names(input) <- names(theta)
    if(!all(names(input) %in% names(theta)))
      stop(paste("the '", ninput,"' list() names are different from theta!", sep = ""))
    for(i in names(theta)) {
      if(!is.list(input[[i]]) & is.list(theta[[i]]))
        input[[i]] <- rep(list(input[[i]]), length.out = length(theta[[i]]))
      if(is.list(theta[[i]])) {
        if(is.null(names(input[[i]])))
          names(input[[i]]) <- names(theta[[i]])
        if(!all(names(input[[i]]) %in% names(theta[[i]])))
          stop(paste(ninput, "list() names for parameter", i, "different from the theta list() names!"))
        if(type != "NA") {
          for(j in names(input[[i]])) {
            if(!inherits(input[[i]][[j]], type)) {
              stop(paste("the '", ninput, "' object for block [", i, "][", j,
                "] is not a ", type,"!", sep = ""), call. = FALSE)
            }
          }
        }
      } else {
        if(type != "NA") {
          if(!inherits(input[[i]], type)) {
            stop(paste("the '", ninput, "' object for block [", i,
              "] is not a ", type,"!", sep = ""), call. = FALSE)
          }
        }
      }
    }
    return(input)
  }

  if(is.null(fitfun)) {
    fitfun <- function(x, p) {
      if(is.null(x)) {
        return(sum(p))
      } else return(drop(x %*% p))
    }
  }

  propose <- parse_input(propose, gmcmc_mvnorm, "function")
  fitfun <- parse_input(fitfun, NULL, "function")
  if(!is.null(data))
    data <- parse_input(data, NULL, "NA")
  if(!is.null(priors))
    priors <- parse_input(priors, NULL, "function")

  if(burnin < 1) burnin <- 1
  if(burnin > n.iter) burnin <- floor(n.iter * 0.1)
  if(thin < 1)
    thin <- 1
  thin <- as.integer(thin)
  iterthin <- as.integer(seq(burnin, n.iter, by = thin))

  rho <- new.env()
  theta.save <- vector(mode = "list", length = length(theta))
  names(theta.save) <- names(theta)
  eta <- theta.save
  for(i in ntheta) {
    eta[[i]] <- 0
    if(is.list(theta[[i]])) {
      for(j in names(theta[[i]]))
        eta[[i]] <- eta[[i]] + fitfun[[i]][[j]](data[[i]][[j]], theta[[i]][[j]])
    } else {
      eta[[i]] <- eta[[i]] + fitfun[[i]](data[[i]], theta[[i]])
    }
  }

  for(i in ntheta) {
    if(is.list(theta[[i]])) {
      theta.save[[i]] <- vector(mode = "list", length = length(theta[[i]]))
      names(theta.save[[i]]) <- names(theta[[i]])
      for(j in names(theta[[i]])) {
        p0 <- propose[[i]][[j]](fun, theta, id = c(i, j),
          prior = priors[[i]][[j]], data = data[[i]][[j]], eta = eta,
          iteration = 1, n.iter = n.iter, burnin = burnin, rho = rho, ...)
        if(!is.list(p0)) {
          stop(paste("the propose() function for block [", i, "][", j,
            "] must return a named list()!", sep = ""))
        }
        if(is.null(p0$alpha)) {
          stop(paste("the propose() function for block [", i, "][", j,
            "] must return the acceptance probability 'alpha'!", sep = ""))
        }
        if(is.null(p0$parameters)) {
          stop(paste("the propose() function for block [", i, "][", j,
            "] must return a vector 'parameters'!", sep = ""))
        }
        if(length(p0$parameters) != length(theta[[i]][[j]])) {
          stop(paste("the propose() function for block [", i, "][", j,
            "] must return a vector 'parameters' with the same length of initial parameters in theta!",
            sep = ""))
        }
        if(is.null(names(p0$parameters)))
          names(p0$parameters) <- paste("[", 1:length(p0$parameters), "]", sep = "")
        if(length(p0$parameters) < 2 & FALSE) ## FIXME!!!
          names(p0$parameters) <- NULL
        theta.save[[i]][[j]] <- list(
          "samples" = matrix(NA, nrow = length(iterthin), ncol = length(p0$parameters)),
          "alpha" = rep(NA, length = length(iterthin)),
          "accepted" = rep(NA, length = length(iterthin))
        )
        colnames(theta.save[[i]][[j]]$samples) <- names(p0$parameters)
        if(!is.null(p0$extra)) {
          theta.save[[i]][[j]]$extra <- matrix(NA, nrow = length(iterthin), ncol = length(p0$extra))
          colnames(theta.save[[i]][[j]]$extra) <- names(p0$extra)
        }
        names(theta[[i]][[j]]) <- names(p0$parameters)
        if(!is.null(attr.copy)) {
          for(a in attr.copy)
            attr(theta[[i]][[j]], a) <- attr(p0$parameters, a)
        }
      }
    } else {
      p0 <- propose[[i]](fun, theta, id = i,
        prior = priors[[i]], data = data[[i]], eta = eta,
        iteration = 1, n.iter = n.iter, burnin = burnin, rho = rho, ...)
      if(!is.list(p0)) {
        stop(paste("the propose() function for block [", i,
          "] must return a named list()!", sep = ""))
      }
      if(is.null(p0$alpha)) {
        stop(paste("the propose() function for block [", i,
          "] must return the acceptance probability 'alpha'!", sep = ""))
      }
      if(is.null(p0$parameters)) {
        stop(paste("the propose() function for block [", i,
          "] must return a vector 'parameters'!", sep = ""))
      }
      if(length(p0$parameters) != length(theta[[i]])) {
        stop(paste("the propose() function for block [", i,
          "] must return a vector 'parameters' with the same length of initial parameters in theta!",
          sep = ""))
      }
      if(is.null(names(p0$parameters)))
        names(p0$parameters) <- paste("[", 1:length(p0$parameters), "]", sep = "")
      if(length(p0$parameters) < 2 & FALSE) ## FIXME!!!
        names(p0$parameters) <- NULL
      theta.save[[i]] <- list(
        "samples" = matrix(NA, nrow = length(iterthin), ncol = length(p0$parameters)),
        "alpha" = rep(NA, length = length(iterthin)),
        "accepted" = rep(NA, length = length(iterthin))
      )
      colnames(theta.save[[i]]$samples) <- names(p0$parameters)
      if(!is.null(p0$extra)) {
        theta.save[[i]]$extra <- matrix(NA, nrow = length(iterthin), ncol = length(p0$extra))
        colnames(theta.save[[i]]$extra) <- names(p0$extra)
      }
      names(theta[[i]]) <- names(p0$parameters)
      if(!is.null(attr.copy)) {
        for(a in attr.copy)
          attr(theta[[i]], a) <- attr(p0$parameters, a)
      }
    }
  }

  ll <- rep(NA, length = length(iterthin))

  nstep <- step
  step <- floor(n.iter / step)

  sampler <- function(...) {
    if(verbose) {
      cat2("Starting the sampler...")
      if(!interactive())
        cat("\n")
    }
    ptm <- proc.time()
    for(iter in 1:n.iter) {
      if(save <- iter %in% iterthin)
        js <- which(iterthin == iter)
      for(i in ntheta) {
        if(is.list(theta[[i]])) {
          for(j in names(theta[[i]])) {
            ## Get proposed states.
            state <- try(propose[[i]][[j]](fun, theta, id = c(i, j),
              prior = priors[[i]][[j]], data = data[[i]][[j]], eta = eta,
              iteration = iter, n.iter = n.iter, burnin = burnin, rho = rho, ...), silent = TRUE)

            if(inherits(state, "try-error")) {
              state <- list("alpha" = NA)
            }

            ## If accepted, set current state to proposed state.
            accepted <- if(is.na(state$alpha)) FALSE else log(runif(1)) <= state$alpha

            if(accepted) {
              if(is.null(attr(theta[[i]][[j]], "fitted.values"))) {
                eta[[i]] <- eta[[i]] - fitfun[[i]][[j]](data[[i]][[j]], theta[[i]][[j]])
              } else {
                eta[[i]] <- eta[[i]] - attr(theta[[i]][[j]], "fitted.values")
              }

              theta[[i]][[j]] <- state$parameters

              if(is.null(attr(theta[[i]][[j]], "fitted.values"))) {
                eta[[i]] <- eta[[i]] + fitfun[[i]][[j]](data[[i]][[j]], theta[[i]][[j]])
              } else {
                eta[[i]] <- eta[[i]] + attr(theta[[i]][[j]], "fitted.values")
              }
            }

            ## Save the samples.
            if(save) {
              theta.save[[i]][[j]]$samples[js, ] <- theta[[i]][[j]]
              theta.save[[i]][[j]]$alpha[js] <- min(c(exp(state$alpha), 1), na.rm = TRUE)
              theta.save[[i]][[j]]$accepted[js] <- accepted
              if(!is.null(state$extra)) {
                theta.save[[i]][[j]]$extra[js, ] <- state$extra
              }
              ll[js] <- if(!is.null(logLik)) {
                logLik(eta)
              } else sum(do.call(fun, c(theta, list(...))[names(formals(fun))]), na.rm = TRUE)

              if(plot) {
                if(js > 1) {
                  plot(ll[1:js], type = "l",
                    xlab = "Iterations", ylab = "logLik")
                }
              }
            }
          }
        } else {
          ## Get proposed states.
          state <- try(propose[[i]](fun, theta, id = i,
            prior = priors[[i]], data = data[[i]], eta = eta,
            iteration = iter, n.iter = n.iter, burnin = burnin, rho = rho, ...), silent = TRUE)

          if(inherits(state, "try-error")) {
            state <- list("alpha" = NA)
          }

          ## If accepted, set current state to proposed state.
          accepted <- if(is.na(state$alpha)) FALSE else log(runif(1)) <= state$alpha

          if(accepted) {
            if (is.null(attr(theta[[i]], "fitted.values"))) {
              eta[[i]] <- eta[[i]] - fitfun[[i]](data[[i]], theta[[i]])
            } else {
              eta[[i]] <- eta[[i]] - attr(theta[[i]], "fitted.values")
            }

            theta[[i]] <- state$parameters

            if (is.null(attr(theta[[i]], "fitted.values"))) {
              eta[[i]] <- eta[[i]] + fitfun[[i]](data[[i]], theta[[i]])
            } else {
              eta[[i]] <- eta[[i]] + attr(theta[[i]], "fitted.values")
            }
          }

          ## Save the samples.
          if(save) {
            theta.save[[i]]$samples[js, ] <- theta[[i]]
            theta.save[[i]]$alpha[js] <- min(c(exp(state$alpha), 1), na.rm = TRUE)
            theta.save[[i]]$accepted[js] <- accepted
            if(!is.null(state$extra)) {
              theta.save[[i]]$extra[js, ] <- state$extra
            }
            ll[js] <- if(!is.null(logLik)) {
              logLik(eta)
            } else sum(do.call(fun, c(theta, list(...))[names(formals(fun))]), na.rm = TRUE)

            if(plot) {
              if(js > 1) {
                plot(ll[1:js], type = "l",
                  xlab = "Iterations", ylab = "logLik")
              }
            }
          }
        }
      }

      if(verbose) barfun(ptm, n.iter, iter, step, nstep)
    }

    if(verbose) cat("\n")

    for(i in ntheta) {
      theta.save[[i]] <- if(is.list(theta[[i]])) {
        lapply(theta.save[[i]], function(x) { do.call("cbind", x) })
      } else do.call("cbind", theta.save[[i]])
      if(is.list(theta[[i]])) {
        for(j in names(theta.save[[i]])) {
          cn <- i
          if(length(theta.save[[i]]) > 1 | !simplify)
            cn <- paste(cn, j, sep = ".")
          cn2 <- colnames(theta.save[[i]][[j]])
          for(k in seq_along(cn2)) {
            cn2[k] <- if(grepl("[", cn2[k], fixed = TRUE)) {
              paste(if(strsplit(cn2[k], "[", fixed = TRUE)[[1]][1] != "") "." else NULL, cn2[k], sep = "")
            } else paste(if(cn2[k] != "") "." else NULL, cn2[k], sep = "")
          }
          cn <- paste(cn, cn2, sep = "")
          colnames(theta.save[[i]][[j]]) <- cn
        }
      } else {
        cn <- colnames(theta.save[[i]])
        sep <- c(".", "")[as.integer(grepl("[", cn, fixed = TRUE)) + 1L]
        for(j in seq_along(cn)) {
          cn[j] <- paste(if(length(theta) > 1) i else NULL,
            if(length(theta[[i]]) > 1) cn[j] else NULL,
            sep = if(length(theta) > 1) sep[j] else "")
        }
        colnames(theta.save[[i]]) <- cn
      }
      theta.save[[i]] <- if(length(theta.save[[i]]) < 2 & simplify) {
        if(is.list(theta[[i]])) theta.save[[i]][[1]] else theta.save[[i]]
      } else {
        if(is.list(theta[[i]])) do.call("cbind", theta.save[[i]]) else theta.save[[i]]
      }
    }
    theta.save <- if(length(theta.save) < 2) theta.save[[1]] else do.call("cbind", theta.save)
    theta.save <- cbind(theta.save, "logLik" = ll)

    theta.save <- mcmc(theta.save, start = burnin, end = n.iter, thin = thin)
    class(theta.save) <- c("gmcmc", class(theta.save))

    return(theta.save)
  }

  sampling <- function(chains, ...) {
    if(is.null(chains))  chains <- 1
    if(chains < 2) {
      return(sampler(...))
    } else {
      samps <- list()
      for(j in 1:chains)
        samps[[j]] <- sampler(...)
      samps <- as.mcmc.list(samps)
      return(samps)
    }
  }

  if(is.null(cores)) cores <- 1
  if(cores < 2) {
    rval <- sampling(chains, ...)
  } else {
    parallel_fun <- function(j) {
      if(j > 1 & !is.null(sleep)) Sys.sleep(sleep)
      sampling(chains, ...)
    }
    rval <- parallel::mclapply(1:cores, parallel_fun, mc.cores = cores)
    if(!inherits(rval, "mcmc.list") & inherits(rval[[1L]], "mcmc.list"))
      rval <- as.mcmc.list(do.call("c", rval))
    if(!inherits(rval, "mcmc.list"))
      rval <- as.mcmc.list(rval)
    if(combine)
      rval <- process.chains(rval)
  }
  
  return(rval)
}


## Print info.
barfun <- function(ptm, n.iter, i, step, nstep, start = TRUE)
{
  ia <- interactive()
  if(i == 10 & start) {
    cat(if(ia) "\r" else "\n")
    elapsed <- c(proc.time() - ptm)[3]
    rt <- elapsed / i * (n.iter - i)
    rt <- if(rt > 60) {
      paste(formatC(format(round(rt / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(rt, 2), nsmall = 2), width = 5), "sec", sep = "")
    cat("|", rep(" ", nstep), "|   0% ", rt, sep = "")
    if(.Platform$OS.type != "unix" & ia) flush.console()
  }
  istep <- i %% step
  if(is.na(istep))
    istep <- 0
  if(istep == 0) {
    cat(if(ia) "\r" else "\n")
    p <- i / n.iter
    p <- paste("|", paste(rep("*", round(nstep * p)), collapse = ""),
      paste(rep(" ", round(nstep * (1 - p))), collapse = ""), "| ",
      formatC(round(p, 2) * 100, width = 3), "%", sep = "")
    elapsed <- c(proc.time() - ptm)[3]
    rt <- elapsed / i * (n.iter - i)
    rt <- if(rt > 60) {
      paste(formatC(format(round(rt / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(rt, 2), nsmall = 2), width = 5), "sec", sep = "")
    elapsed <- if(elapsed > 60) {
      paste(formatC(format(round(elapsed / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(elapsed, 2), nsmall = 2), width = 5), "sec", sep = "")
    cat(p, rt, elapsed, sep = " ")
    if(.Platform$OS.type != "unix" & ia) flush.console()
  }
}


DIC <- function(object, ...)
{
  UseMethod("DIC")
}

DIC.gmcmc <- function(object, ...)
{
  object <- list(object, ...)
  rval <- NULL
  for(j in seq_along(object)) {
    dev <- -2 * object[[j]][, "logLik"]
    pd <- var(dev, na.rm = TRUE) / 2
    rval <- rbind(rval, c("DIC" = mean(dev, na.rm = TRUE), "pd" = pd))
  }
  Call <- match.call()
  row.names(rval) <- if(nrow(rval) > 1) as.character(Call[-1L]) else ""
  rval
}

plot.gmcmc <- function(x, ...)
{
  x <- na.omit(x)
  class(x) <- "mcmc"
  plot(x, ...)
}

"[.gmcmc.theta" <- function(x, i) {
  if(length(i) < 2)
    return(x[[i]])
  else
    return(x[[i[1]]][[i[2]]])
}

"[<-.gmcmc.theta" <- function(x, i, value) {
  if(length(i) < 2)
    x[[i]] <- value
  else
    x[[i[1]]][[i[2]]] <- value
  return(x)
}

gmcmc_mvnorm <- function(fun, theta, id, prior, ...)
{
  args <- list(...)
  iteration <- args$iteration

  if(is.null(attr(theta[id], "scale")))
    attr(theta[id], "scale") <- 1
  if(is.null(attr(theta[id], "P")))
    attr(theta[id], "sigma") <- diag(length(theta[id]))

  if(iteration < floor(0.15 * args$n.iter)) {
    scale <- attr(theta[id], "scale")
    sigma <- attr(theta[id], "sigma")

    k <- 1
    do <- TRUE
    while(do & k < 100) {
      ll0 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
      p0 <- if(is.null(prior)) {
        sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
      } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

      theta[id] <- drop(mvtnorm::rmvnorm(n = 1, mean = theta[id], sigma = scale * sigma))

      ll1 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
      p1 <- if(is.null(prior)) {
        sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
      } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

      alpha <- drop((ll1 + p1) - (ll0 + p0))

      if(!is.finite(alpha)) next

      accepted <- if(is.na(alpha)) FALSE else log(runif(1)) <= alpha

      scale <- if(alpha < log(0.23)) {
        scale - 0.1 * scale
      } else {
        scale + 0.1 * scale
      }

      if(accepted) do <- FALSE

      k <- k + 1
    }

    attr(theta[id], "scale") <- scale
    attr(theta[id], "sigma") <- sigma
  }

  scale <- attr(theta[id], "scale")
  sigma <- attr(theta[id], "sigma")

  ll0 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p0 <- if(is.null(prior)) {
    sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  theta[id] <- drop(mvtnorm::rmvnorm(n = 1, mean = theta[id], sigma = scale * sigma))

  ll1 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p1 <- if(is.null(prior)) {
    sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  alpha <- drop((ll1 + p1) - (ll0 + p0))

  rval <- list("parameters" = theta[id], "alpha" = alpha)

  attr(rval$parameters, "scale") <- scale
  attr(rval$parameters, "sigma") <- sigma

  rval
}


gmcmc_unislice <- function(fun, theta, id, prior, j, ...,
  w = 1, m = 30, lower = -Inf, upper = +Inf)
{
  args <- list(...)

  x0 <- theta[id][j]
  gL <- gR <- theta

  fun2 <- function(theta) {
    ll <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
    lp <- if(is.null(prior)) {
      sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
    } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)
    return(ll + lp)
  }

  gx0 <- fun2(theta)

  ## Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)

  ## Find the initial interval to sample from.
  ## w <- w * abs(x0) ## FIXME???
  u <- runif(1, 0, w)
  gL[id][j] <- theta[id][j] - u
  gR[id][j] <- theta[id][j] + (w - u)  ## should guarantee that g[j] is in [L, R], even with roundoff

  ## Expand the interval until its ends are outside the slice, or until
  ## the limit on steps is reached.
  if(is.infinite(m)) {
    repeat {
      if(gL[id][j] <= lower) break
      if(fun2(gL) <= logy) break
      gL[id][j] <- gL[id][j] - w
    }
    repeat {
      if(gR[id][j] >= upper) break
      if(fun2(gR) <= logy) break
      gR[id][j] <- gR[id][j] + w
    }
  } else {
    if(m > 1) {
      J <- floor(runif(1, 0, m))
      K <- (m - 1) - J
      while(J > 0) {
        if(gL[id][j] <= lower) break
        if(fun2(gL) <= logy) break
        gL[id][j] <- gL[id][j] - w
        J <- J - 1
      }
      while(K > 0) {
        if(gR[id][j] >= upper) break
        if(fun2(gR) <= logy) break
        gR[id][j] <- gR[id][j] + w
        K <- K - 1
      }
    }
  }

  ## Shrink interval to lower and upper bounds.
  if(gL[id][j] < lower) {
    gL[id][j] <- lower
  }
  if(gR[id][j] > upper) {
    gR[id][j] <- upper
  }

  ## Sample from the interval, shrinking it on each rejection.
  repeat {
    theta[id][j] <- runif(1, gL[id][j], gR[id][j])

    gx1 <- fun2(theta)

    if(gx1 >= logy) break

    if(theta[id][j] > x0) {
      gR[id][j] <- theta[id][j]
    } else {
      gL[id][j] <- theta[id][j]
    }
  }

  ## Return the point sampled
  return(theta[id][j])
}


gmcmc_slice <- function(fun, theta, id, prior, ...)
{
  for(j in seq_along(theta[id]))
    theta[id][j] <- gmcmc_unislice(fun, theta, id, prior, j, ...)
  return(list("parameters" = theta[id], "alpha" = log(1)))
}


dmvnorm_log <- function(x, mean, sigma)
{
  d <- try(dmvnorm(x, mean = mean, sigma = sigma, log = TRUE), silent = TRUE)
  if(inherits(d, "try-error"))
    d <- NA
  return(d)
}


GMCMC_iwlsC <- function(family, theta, id, eta, y, data,
  weights = NULL, offset = NULL, zworking, resids, rho, ...)
{
  if(!is.null(offset)) {
    for(j in names(offset))
      eta[[j]] <- eta[[j]] + offset[[j]]
  }

  W <- if(is.null(weights[[id[1]]])) 1.0 else weights[[id[1]]]
  rval <- .Call("gmcmc_iwls", family, theta, id, eta, y, data,
    zworking, resids, id[1], W, rho, PACKAGE = "bamlss")

  ## Sample variance parameter.
  if(!data$fixed & !data$fxsp & length(data$S)) {
    if(length(data$S) > 1) {
      i <- grep("tau2", names(rval$parameters))
      for(j in i) {
        rval$parameters <- uni.slice(rval$parameters, data, family, NULL,
          NULL, id[1], j, logPost = gmcmc_logPost, lower = 0, ll = rval$loglik)
      }
    }
  }

  return(list("parameters" = rval$parameters, "alpha" = rval$alpha, "extra" = c("edf" = rval$edf)))
}


GMCMC_iwlsC_gp <- function(family, theta, id, eta, y, data,
  weights = NULL, offset = NULL, zworking, resids, rho, ...)
{
  if(!is.null(offset)) {
    for(j in names(offset))
      eta[[j]] <- eta[[j]] + offset[[j]]
  }

  W <- if(is.null(weights[[id[1]]])) 1.0 else weights[[id[1]]]
  rval <- .Call("gmcmc_iwls_gp", family, theta, id, eta, y, data,
    zworking, resids, id[1], W, data$sparse.setup$block.index, data$sparse.setup$is.diagonal, rho)

  ## Sample variance parameter.
  if(!data$fixed & !data$fxsp & length(data$S)) {
    if((length(data$S) < 2) & (attr(data$prior, "var_prior") == "ig")) {
      g <- get.par(rval$parameters, "b")
      if(is.function(data$S[[1]])) {
        K <- data$S[[1]](g)
        data$rank <- ncol(K)
      } else K <- data$S[[1]]
      a <- data$rank / 2 + data$a
      b <- 0.5 * crossprod(g, K) %*% g + data$b
      tau2 <- 1 / rgamma(1, a, b)
      rval$parameters <- set.par(rval$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(rval$parameters))
      for(j in i) {
        rval$parameters <- uni.slice(rval$parameters, data, family, NULL,
          NULL, id[1], j, logPost = gmcmc_logPost, lower = 0, ll = rval$loglik)
      }
    }
  }

  return(list("parameters" = rval$parameters, "alpha" = rval$alpha, "extra" = c("edf" = rval$edf)))
}


GMCMC_iwlsC_gp_diag_lasso <- function(family, theta, id, eta, y, data,
  weights = NULL, offset = NULL, zworking, resids, rho, ...)
{
  if(!is.null(offset)) {
    for(j in names(offset))
      eta[[j]] <- eta[[j]] + offset[[j]]
  }

  W <- if(is.null(weights[[id[1]]])) 1.0 else weights[[id[1]]]
  rval <- .Call("gmcmc_iwls_gp_diag_lasso", family, theta, id, eta, y, data,
    zworking, resids, id[1], W, rho, PACKAGE = "bamlss")

  ## Sample variance parameter.
  if(!data$fixed & !data$fxsp & length(data$S)) {
    if((length(data$S) < 2) & (attr(data$prior, "var_prior") == "ig")) {
      g <- get.par(rval$parameters, "b")
      if(is.function(data$S[[1]])) {
        K <- data$S[[1]](g)
        data$rank <- ncol(K)
      } else K <- data$S[[1]]
      a <- data$rank / 2 + data$a
      b <- 0.5 * crossprod(g, K) %*% g + data$b
      tau2 <- 1 / rgamma(1, a, b)
      rval$parameters <- set.par(rval$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(rval$parameters))
      for(j in i) {
        rval$parameters <- uni.slice(rval$parameters, data, family, NULL,
          NULL, id[1], j, logPost = gmcmc_logPost, lower = 0, ll = rval$loglik)
      }
    }
  }

  return(list("parameters" = rval$parameters, "alpha" = rval$alpha, "extra" = c("edf" = rval$edf)))
}


process.derivs <- function(x, is.weight = FALSE)
{
  .Call("process_derivs", as.numeric(x), as.logical(is.weight), PACKAGE = "bamlss")
}


GMCMC_iwls <- function(family, theta, id, eta, y, data, weights = NULL, offset = NULL, ...)
{
  theta <- theta[[id[1]]][[id[2]]]

  if(is.null(attr(theta, "fitted.values")))
    attr(theta, "fitted.values") <- data$fit.fun(data$X, theta)

  if(!is.null(offset)) {
    for(j in names(offset))
      eta[[j]] <- eta[[j]] + offset[[j]]
  }

  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute weights.
  hess <- family$hess[[id[1]]](y, peta, id = id[1])

  if(!is.null(weights[[id[1]]]))
    hess <- hess * weights[[id[1]]]

  hess <- process.derivs(hess, is.weight = TRUE)

  ## Score.
  score <- process.derivs(family$score[[id[1]]](y, peta, id = id[1]), is.weight = FALSE)

  ## Compute working observations.
  z <- eta[[id[1]]] + 1 / hess * score

  ## Compute old log likelihood and old log coefficients prior.
  pibeta <- family$loglik(y, peta)

  if(inherits(data, "nnet0.smooth")) {
    data$X <- data$getZ(data$X, theta)
  }

  ## Sample variance parameter.
  if(!data$fixed & !data$fxsp & length(data$S)) {
    if((length(data$S) < 2) & (attr(data$prior, "var_prior") == "ig")) {
      g <- get.par(theta, "b")
      if(is.function(data$S[[1]])) {
        K <- data$S[[1]](c(g, data$fixed.hyper))
        data$rank <- ncol(K)
      } else K <- data$S[[1]]
      a <- data$rank / 2 + data$a
      b <- 0.5 * crossprod(g, K) %*% g + data$b
      tau2 <- 1 / rgamma(1, a, b)
      theta <- set.par(theta, tau2, "tau2")
    } else {
      i <- grep("tau2", names(theta))
      for(j in i) {
        theta <- uni.slice(theta, data, family, NULL,
          NULL, id[1], j, logPost = gmcmc_logPost, lower = 0, ll = pibeta)
      }
    }
  }

  p1 <- data$prior(theta)

  ## Compute partial predictor.
  eta2 <- eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values")

  ## Compute reduced residuals.
  e <- z - eta2
  xbin.fun(data$binning$sorted.index, hess, e, data$weights, data$rres, data$binning$order)

  ## Compute mean and precision.
  XWX <- do.XWX(data$X, 1 / data$weights, data$sparse.setup$matrix)

  S <- 0
  P <- if(data$fixed) {
    if((k <- ncol(data$X)) < 2) {
      1 / XWX
    } else matrix_inv(XWX + if(!is.null(data$xt[["pS"]])) data$xt[["pS"]] else 0, data$sparse.setup)
  } else {
    tau2 <- get.par(theta, "tau2")
    for(j in seq_along(data$S))
      S <- S + 1 / tau2[j] * if(is.function(data$S[[j]])) data$S[[j]](c(theta, data$fixed.hyper)) else data$S[[j]]
    matrix_inv(XWX + S + if(!is.null(data$xt[["pS"]])) data$xt[["pS"]] else 0, data$sparse.setup)
  }
  P[P == Inf] <- 0
  if(is.null(data$xt[["pm"]])) {
    M <- P %*% crossprod(data$X, data$rres)
  } else {
    pS <- if(!is.null(data$xt[["pS"]])) {
      data$xt[["pS"]]
    } else {
      if(!is.null(data$xt[["pSa"]])) {
        1 / tau2[length(tau2)] * data$xt[["pSa"]]
      } else 0
    }
    M <- P %*% (crossprod(data$X, data$rres) + pS %*% data$xt[["pm"]])
  }

  ## Degrees of freedom.
  edf <- sum_diag(XWX %*% P)

  ## Save old coefficients
  g0 <- drop(get.par(theta, "b"))

  ## Sample new parameters.
  g <- try(drop(rmvnorm(n = 1, mean = M, sigma = P, method = "chol")), silent = TRUE)
  if(inherits(g, "try-error")) {
    return(list("parameters" = theta, "alpha" = -Inf, "extra" = c("edf" = NA)))
  }
  if(!is.null(data$doCmat)) {
    V <- P %*% t(data$C)
    W <- data$C %*% V
    U <- chol2inv(chol(W)) %*% t(V)
    g <- drop(g - t(U) %*% data$C %*% g)
  }

  ## Compute log priors.
  p2 <- data$prior(c("b" = g, get.par(theta, "tau2")))
  qbetaprop <- try(dmvnorm(g, mean = M, sigma = P, log = TRUE), silent = TRUE)
  if(inherits(qbetaprop, "try-error")) {
    return(list("parameters" = theta, "alpha" = -Inf, "extra" = c("edf" = NA)))
  }

  ## Compute fitted values.        
  attr(theta, "fitted.values") <- data$fit.fun(data$X, g)

  ## Set up new predictor.
  eta[[id[1]]] <- eta[[id[1]]] + attr(theta, "fitted.values")

  ## Map predictor to parameter scale.
  peta <- family$map2par(eta)

  ## Compute new log likelihood.
  pibetaprop <- family$loglik(y, peta)

  ## Compute new weights
  hess <- family$hess[[id[1]]](y, peta, id = id[1])

  if(!is.null(weights[[id[1]]]))
    hess <- hess * weights[[id[1]]]

  hess <- process.derivs(hess, is.weight = TRUE)

  ## New score.
  score <- process.derivs(family$score[[id[1]]](y, peta, id = id[1]), is.weight = FALSE)

  ## New working observations.
  z <- eta[[id[1]]] + 1 / hess * score

  ## Compute reduced residuals.
  e <- z - eta2
  xbin.fun(data$binning$sorted.index, hess, e, data$weights, data$rres, data$binning$order)

  ## New penalty.
  S <- 0
  if(!data$fixed) {
    tau2 <- get.par(theta, "tau2")
    for(j in seq_along(data$S))
      S <- S + 1 / tau2[j] * if(is.function(data$S[[j]])) data$S[[j]](c(g, data$fixed.hyper)) else data$S[[j]]
  }

  ## Compute mean and precision.
  XWX <- do.XWX(data$X, 1 / data$weights, data$sparse.setup$matrix)
  P2 <- if(data$fixed) {
    if(k < 2) {
      1 / (XWX)
    } else matrix_inv(XWX + if(!is.null(data$xt[["pS"]])) data$xt[["pS"]] else 0, data$sparse.setup)
  } else {
    matrix_inv(XWX + S + if(!is.null(data$xt[["pS"]])) data$xt[["pS"]] else 0, data$sparse.setup)
  }
  P2[P2 == Inf] <- 0
  if(is.null(data$xt[["pm"]])) {
    M2 <- P2 %*% crossprod(data$X, data$rres)
  } else {
    pS <- if(!is.null(data$xt[["pS"]])) {
      data$xt[["pS"]]
    } else {
      if(!is.null(data$xt[["pSa"]])) {
        1 / tau2[length(tau2)] * data$xt[["pSa"]]
      } else 0
    }
    M2 <- P2 %*% (crossprod(data$X, data$rres) + pS %*% data$xt[["pm"]])
  }

  ## Get the log prior.
  qbeta <- try(dmvnorm(g0, mean = M2, sigma = P2, log = TRUE), silent = TRUE)
  if(inherits(qbeta, "try-error")) {
    return(list("parameters" = theta, "alpha" = -Inf, "extra" = c("edf" = NA)))
  }

  theta <- set.par(theta, g, "b")
  data$state$parameters <- as.numeric(theta)
  names(data$state$parameters) <- names(theta)

  ## Compute acceptance probablity.
  alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

#cat("\n-----------\n")
#cat("pibetaprop", pibetaprop, "\n")
#cat("qbeta", qbeta, "\n")
#cat("p2", p2, "\n")
#cat("pibeta", pibeta, "\n")
#cat("qbetaprop", qbetaprop, "\n")
#cat("p1", p1, "\n")
#cat("alpha", exp(alpha), "\n")
#print(data$label)

  return(list("parameters" = theta, "alpha" = alpha, "extra" = c("edf" = edf)))
}


GMCMC_mvn <- function(family, theta, id, prior, eta, y, data, ...)
{
  theta <- theta[[id[1]]][[id[2]]]

  if(is.null(attr(theta, "fitted.values")))
    attr(theta, "fitted.values") <- data$fit.fun(data$X, theta)

  ll1 <- family$loglik(y, family$map2par(eta))
  p1 <- data$prior(theta)

  eta2 <- eta
  eta2[[id[1]]] <- eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values")
  if(!is.null(data$state)) {
    if(!is.null(data$state$hessian))
      attr(theta, "hess") <- data$state$hessian
    if(is.null(attr(theta, "hess"))) {
      g <- get.par(data$state$parameters, "b")
      tau2 <- if(!data$fixed) get.par(data$state$parameters, "tau2") else NULL

      lp <- function(g) {
        eta[[id[1]]] <- eta[[id[1]]] + data$fit.fun(data$X, g)
        family$loglik(y, family$map2par(eta)) + data$prior(c(g, tau2))
      }

      gfun <- hfun <- NULL

      g.grad <- grad(fun = lp, theta = g, id = id[1], prior = NULL,
        args = list("gradient" = gfun, "x" = data, "y" = y, "eta" = eta))

      g.hess <- hess(fun = lp, theta = g, id = id[1], prior = NULL,
        args = list("gradient" = gfun, "hessian" = hfun, "x" = data, "y" = y, "eta" = eta))

      attr(theta, "hess") <- matrix_inv(g.hess)
    }
  }

  g <- drop(rmvnorm(n = 1, mean = drop(get.par(theta, "b")), sigma = attr(theta, "hess")))
  theta <- set.par(theta, g, "b")

  attr(theta, "fitted.values") <- data$fit.fun(data$X, theta)

  eta[[id[1]]] <- eta[[id[1]]] + attr(theta, "fitted.values")
  ll2 <- family$loglik(y, family$map2par(eta))
  p2 <- data$prior(theta)

  ## Sample variance parameter.
  if(!data$fixed & !data$fxsp & length(data$S)) {
    if(length(data$S) < 2) {
      g <- get.par(theta, "b")
      a <- data$rank / 2 + data$a
      b <- 0.5 * crossprod(g, data$S[[1]]) %*% g + data$b
      tau2 <- 1 / rgamma(1, a, b)
      theta <- set.par(theta, tau2, "tau2")
    } else {
      i <- grep("tau2", names(theta))
      for(j in i) {
        theta <- uni.slice(theta, data, family, NULL,
          NULL, id[1], j, logPost = gmcmc_logPost, lower = 0, ll = ll1)
      }
    }
  }

  ## Compute acceptance probablity.
  alpha <- drop((ll2 + p2) - (ll1 + p1))

  data$state$parameters <- as.numeric(theta)
  names(data$state$parameters) <- names(theta)

  rval <- list("parameters" = theta, "alpha" = alpha)
  rval
}


GMCMC_newton <- function(family, theta, id, prior, eta, y, data, ...)
{
  theta <- theta[id]
  g <- get.par(theta, "b")
  tau2 <- if(!data$fixed) get.par(theta, "tau2") else NULL
  nu <- if(is.null(data$nu)) 0.1 else data$nu

  pibeta <- family$loglik(y, family$map2par(eta))
  p1 <- data$prior(theta)

  if(is.null(attr(theta, "fitted.values")))
    attr(theta, "fitted.values") <- data$fit.fun(data$X, theta)

  eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values")

  lp <- function(g) {
    eta[[id[1]]] <- eta[[id[1]]] + data$fit.fun(data$X, g)
    family$loglik(y, family$map2par(eta)) + data$prior(c(g, tau2))
  }

  if(is.null(family$gradient[[id[1]]])) {
    gfun <- NULL
  } else {
    gfun <- list()
    gfun[[id[1]]] <- function(g, y, eta, x, ...) {
      gg <- family$gradient[[id[1]]](g, y, eta, x)
      if(!is.null(data$grad)) {
        gg <- gg + data$grad(score = NULL, c(g, tau2), full = FALSE)
      }
      drop(gg)
    }
  }

  if(is.null(family$hessian[[id[1]]])) {
    hfun <- NULL
  } else {
    hfun <- list()
    hfun[[id[1]]] <- function(g, y, eta, x, ...) {
      hg <- family$hessian[[id[1]]](g, y, eta, x)
      if(!is.null(data$hess)) {
        hg <- hg + data$hess(score = NULL, c(g, tau2), full = FALSE)
      }
      hg
    }
  }

  g.grad <- grad(fun = lp, theta = g, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "x" = data, "y" = y, "eta" = eta))

  g.hess <- hess(fun = lp, theta = g, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "hessian" = hfun, "x" = data, "y" = y, "eta" = eta))

  Sigma <- matrix_inv(g.hess)
  mu <- drop(g + nu * Sigma %*% g.grad)

  g2 <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma))
  names(g2) <- names(g)

  p2 <- data$prior(c("b" = g2, tau2))
  qbetaprop <- dmvnorm(g2, mean = mu, sigma = Sigma, log = TRUE)

  g.grad2 <- grad(fun = lp, theta = g2, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "x" = data, "y" = y, "eta" = eta))

  g.hess2 <- hess(fun = lp, theta = g2, id = id[1], prior = NULL,
    args = list("gradient" = gfun, "hessian" = hfun, "x" = data, "y" = y, "eta" = eta))

  Sigma2 <- matrix_inv(g.hess2)
  mu2 <- drop(g2 + nu * Sigma2 %*% g.grad2)

  theta <- set.par(theta, g2, "b")

  attr(theta, "fitted.values") <- data$fit.fun(data$X, g2)
  eta[[id[1]]] <- eta[[id[1]]] + attr(theta, "fitted.values")

  pibetaprop <- family$loglik(y, family$map2par(eta))

  qbeta <- dmvnorm(matrix(g, nrow = 1), mean = mu2, sigma = Sigma2, log = TRUE)

  alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

  ## Sample variance parameter.
  if(!data$fixed & !data$fxsp) {
    if(!data$fixed & !data$fxsp) {
      tau2 <- NULL
      for(j in seq_along(data$S)) {
        a <- data$rank[j] / 2 + data$a
        b <- 0.5 * crossprod(g2, data$S[[j]]) %*% g2 + data$b
        tau2 <- c(tau2, 1 / rgamma(1, a, b))
      }
      theta <- set.par(theta, tau2, "tau2")
    }
  }

  data$state$parameters <- as.numeric(theta)
  names(data$state$parameters) <- names(theta)

  rval <- list("parameters" = theta, "alpha" = alpha, "extra" = c("edf" = data$edf(data)))
}


gmcmc_logPost <- function(g, x, family, y = NULL, eta = NULL, id, ll = NULL)
{
  if(is.null(ll)) {
    eta[[id[1]]] <- eta[[id[1]]] + x$fit.fun(x$X, g)
    ll <- family$loglik(y, family$map2par(eta))
  }
  lp <- x$prior(g)

  return(ll + lp)
}


GMCMC_slice <- function(family, theta, id, eta, y, data, ...)
{
  theta <- theta[[id[1]]][[id[2]]]

  if(is.null(attr(theta, "fitted.values")))
    attr(theta, "fitted.values") <- data$fit.fun(data$X, theta)

  ## Remove fitted values.
  eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values")

  attr(theta, "fitted.values") <- NULL

  ## Sample coefficients.
  i <- 1:length(theta)
  if(!data$fixed)
    i <- i[-grep("tau2", names(theta))]
  for(j in i) {
    theta <- uni.slice(theta, data, family, y,
      eta, id[1], j, logPost = gmcmc_logPost)
  }

  ## New fitted values.
  fit <- data$fit.fun(data$X, theta)

  ## Sample variance parameter.
  if(!data$fixed & !data$fxsp & length(data$S)) {
    ## New logLik.
    eta[[id[1]]] <- eta[[id[1]]] + fit
    ll <- family$loglik(y, family$map2par(eta))
    i <- grep("tau2", names(theta))
    for(j in i) {
      theta <- uni.slice(theta, data, family, NULL,
        NULL, id[1], j, logPost = gmcmc_logPost, lower = 0, ll = ll)
    }
  }

  ## New theta.
  attr(theta, "fitted.values") <- fit

  return(list("parameters" = theta, "alpha" = log(1)))
}


gmcmc_mvnorm2 <- function(fun, theta, id, prior, ...)
{
  args <- list(...)
  iteration <- args$iteration

  ll0 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p0 <- if(is.null(prior)) {
    sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  adapt <- if(is.null(args$adapt)) args$burnin else args$adapt 
  if(iteration <= adapt) {
    eps <- attr(theta[id], "eps")
    if(is.null(eps))
      eps <- 1
    if(eps > 0.001) {
      objfun <- function(par) {
        theta[id] <- par
        ll <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
        lp <- if(is.null(prior)) {
          sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
        } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)
        return(-1 * (ll + lp))
      }
      start <- attr(theta[id], "mode")
      if(is.null(start))
        start <- theta[id]
      opt <- optim(start, fn = objfun, method = "L-BFGS-B", hessian = TRUE,
        control = list(maxit = 100))
      theta[id] <- opt$par
      attr(theta[id], "sigma") <- matrix_inv(diag(diag(opt$hessian)))
      attr(theta[id], "scale") <- 1
      attr(theta[id], "eps") <- mean(abs((start - opt$par) / start), na.rm = TRUE)
      attr(theta[id], "mode") <- opt$par
    }
  }

  theta.attr <- attributes(theta[id])

  scale <- theta.attr$scale
  sigma <- theta.attr$sigma
  sigma <- scale * sigma

  theta[id] <- drop(rmvnorm(n = 1, mean = theta[id], sigma = sigma))

  ll1 <- sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  p1 <- if(is.null(prior)) {
    sum(dnorm(theta[id], sd = 1000, log = TRUE), na.rm = TRUE)
  } else sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)

  alpha <- drop((ll1 + p1) - (ll0 + p0))

  rval <- list("parameters" = theta[id], "alpha" = alpha)
  attributes(rval$parameters) <- theta.attr

  rval
}


eval_fun <- function(fun, theta, prior, id, args)
{
  ll <- if(is.list(theta)) {
    sum(do.call(fun, c(theta, args)[names(formals(fun))]), na.rm = TRUE)
  } else sum(fun(theta), na.rm = TRUE)
  lp <- if(is.null(prior)) {
    if(is.list(theta)) {
      if(is.list(theta[[id[1]]])) {
        sum(dnorm(theta[[id[1]]][[id[2]]], sd = 1000, log = TRUE), na.rm = TRUE)
      } else {
        sum(dnorm(theta[[id[1]]], sd = 1000, log = TRUE), na.rm = TRUE)
      }
    } else {
      sum(dnorm(theta, sd = 1000, log = TRUE), na.rm = TRUE)
    }
  } else {
    if(is.list(theta)) {
      sum(do.call(prior, c(theta, args)[names(formals(prior))]), na.rm = TRUE)
    } else {
      sum(prior(theta), na.rm = TRUE)
    }
  }
  ll + lp
}


grad <- function(fun, theta, prior, id, args, eps = .Machine$double.eps^0.25)
{
  if(!is.null(args$gradient[[id[1]]])) {
    gfun <- args$gradient[[id[1]]]
    args[[names(formals(gfun))[1]]] <- theta
    grad.theta <- do.call(gfun, args[names(formals(gfun))])
  } else {
    gfun <- function(theta, j) {
      theta2 <- theta
      if(is.list(theta)) {
        if(is.list(theta[[id[1]]])) {
          theta2[[id[1]]][[id[2]]][j] <- theta2[[id[1]]][[id[2]]][j] + eps
          theta[[id[1]]][[id[2]]][j] <- theta[[id[1]]][[id[2]]][j] - eps
        } else {
          theta2[[id[1]]][j] <- theta2[[id[1]]][j] + eps
          theta[[id[1]]][j] <- theta[[id[1]]][j] - eps
        }
      } else {
        theta2[j] <- theta2[j] + eps
        theta[j] <- theta[j] - eps
      }
      lp1 <- eval_fun(fun, theta2, prior, id, args)
      lp2 <- eval_fun(fun, theta, prior, id, args)
      drop((lp1 - lp2) / (2 * eps))
    }
    k <- if(is.list(theta)) {
      if(is.list(theta[[id[1]]])) {
        length(theta[[id[1]]][[id[2]]])
      } else length(theta[[id[1]]])
    } else length(theta)
    grad.theta <- rep(NA, k)
    for(j in 1:k) {
      grad.theta[j] <- gfun(theta, j)
    }
  }
  names(grad.theta) <- if(is.list(theta)) {
    if(is.list(theta[[id[1]]])) {
      names(theta[[id[1]]][[id[2]]])
    } else names(theta[[id[1]]])
  } else names(theta)
  grad.theta
}

hess <- function(fun, theta, prior, id, args, diag = FALSE)
{
  if(!is.null(args$hessian)) {
    hfun <- args$hessian[[id[1]]]
    args[[names(formals(hfun))[1]]] <- theta
    H <- do.call(hfun, args[names(formals(hfun))])
  } else {
    theta2 <- theta
    objfun <- function(par) {
      if(is.list(theta)) {
        if(is.list(theta[[id[1]]])) {
          theta2[[id[1]]][[id[2]]] <- par
        } else theta2[[id[1]]] <- par
        return(eval_fun(fun, theta2, prior, id, args))
      } else {
        return(eval_fun(fun, par, prior, id, args))
      }
    }
    gfun <- NULL
    if(!is.null(args$gradient)) {
      gfun2 <- args$gradient[[id[1]]]
      gfun <- function(par) {
        if(is.list(theta)) {
          if(is.list(theta[[id[1]]])) {
            theta2[[id[1]]][[id[2]]] <- par
          } else theta2[[id[1]]] <- par
          return(do.call(gfun2, c(theta2, args)[names(formals(gfun2))]))
        } else {
          args[[names(formals(gfun2))[1]]] <- par
          return(do.call(gfun2, args[names(formals(gfun2))]))
        }
      }
    }
    start <- if(is.list(theta)) {
      if(is.list(theta[[id[1]]])) {
        theta[[id[1]]][[id[2]]]
      } else theta[[id[1]]]
    } else theta
    H <- -1 * optimHess(start, fn = objfun, gr = gfun, control = list(fnscale = -1))
  }
  if(diag)
    H <- diag(diag(H))
  nh <- if(is.list(theta)) {
    if(is.list(theta[[id[1]]])) {
      names(theta[[id[1]]][[id[2]]])
    } else names(theta[[id[1]]])
  } else names(theta)
  rownames(H) <- colnames(H) <- nh
  H
}


gmcmc_newton <- function(fun, theta, id, prior, ...)
{
  args <- list(...)

  adapt <- if(is.null(args$adapt)) {
    if(args$burnin < 1) 1 else args$burnin
  } else args$adapt
  eps0 <- if(is.null(args$eps0)) 1e-04 else args$eps0
  if(args$iteration <= adapt) {
    eps <- attr(theta[id], "eps")
    if(is.null(eps))
      eps <- 1
    if(is.na(eps))
      eps <- 1
    if(eps > eps0) {
      theta2 <- theta
      objfun <- function(par, ...) {
        theta2[id] <- par
        -1 * eval_fun(fun, theta2, prior, id, args)
      }
      gfun <- NULL
      if(!is.null(args$gradient)) {
        gfun2 <- args$gradient[[id[1]]]
        gfun <- function(par, ...) {
          theta2[id] <- par
          -1 * do.call(gfun2, c(theta2, args)[names(formals(gfun2))])
        }
      }

      start <- attr(theta[id], "mode")
      if(is.null(start))
        start <- theta[id]

      opt <- optim(start, fn = objfun, gr = gfun, method = "L-BFGS-B")
      
      rval <- list("parameters" = opt$par, "alpha" = log(10))
      attr(rval$parameters, "eps") <- mean(abs((start - opt$par) / start), na.rm = TRUE)
      attr(rval$parameters, "mode") <- opt$par
      mode(rval$parameters) <- "numeric"
    } else {
      rval <- list("parameters" = attr(theta[id], "mode"), "alpha" = log(10))
    }
  } else {
    p.old <- eval_fun(fun, theta, prior, id, args)
    theta2 <- theta

    grad.theta <- grad(fun, theta, prior, id, args)
    hess.theta <- hess(fun, theta, prior, id, args, diag = FALSE)

    Sigma <- matrix_inv(hess.theta)
    mu <- drop(theta[id] + Sigma %*% grad.theta)

    q.prop <- dmvnorm(matrix(theta[id], nrow = 1), mean = mu, sigma = Sigma, log = TRUE)

    theta2[id] <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma))

    grad.theta2 <- grad(fun, theta2, prior, id, args)
    hess.theta2 <- hess(fun, theta2, prior, id, args, diag = FALSE)

    Sigma2 <- matrix_inv(hess.theta2)
    mu2 <- drop(theta2[id] + Sigma2 %*% grad.theta2)

    p.prop <- eval_fun(fun, theta2, prior, id, args)

    q.old <- dmvnorm(matrix(theta[id], nrow = 1), mean = mu2, sigma = Sigma2, log = TRUE)

    alpha <- (p.prop - p.old) + (q.old - q.prop)

    rval <- list("parameters" = theta2[id], "alpha" = alpha)
  }

  rval
}


## Naive sampler.
sam_MVNORM <- MVNORM <- function(x, y = NULL, family = NULL, start = NULL, n.samples = 500, hessian = NULL, ...)
{
  if(!is.null(start))
    start <- unlist(start)
  if(is.null(hessian)) {
    if(inherits(x, "bamlss")) {
      if(is.null(x$hessian)) {
        hessian <- opt_optim(x = x$x, y = x$y, family = x$family,
          start = start, hessian = TRUE, ...)
      } else {
        hessian <- x$hessian
      }
    } else {
      hessian <- opt_optim(x = x, y = y, family = family, start = start, hessian = TRUE, ...)
    }
  }

  if(is.null(hessian))
    stop("need the hessian for sampling!")

  par <- if(inherits(x, "bamlss")) {
    unlist(x$parameters)
  } else {
    if(is.null(start)) get.all.par(x, list = FALSE, drop = TRUE) else start
  }
  if(length(par) < 1)
    par <- start
  if(is.null(par))
    stop("need parameters for sampling!")

  if(length(i <- grep2(c(".alpha", ".edf", ".tau2", ".accepted"), names(par), fixed = TRUE)))
    par <- par[-i]

  npar <- names(par)
  hessian <- hessian[npar, npar]
  Sigma <- matrix_inv(-1 * hessian)
  if(ncol(Sigma) != length(par)) {
    for(i in c(1e-05, 1e-04, 1e-03, 1e-02)) {
      if(ncol(Sigma) != length(par)) {
        hessian2 <- hessian + diag(i, ncol(hessian))
        Sigma <- matrix_inv(-1 * hessian2)
      }
    }
  }

  samps <- rmvnorm(n.samples, mean = par, sigma = Sigma)
  colnames(samps) <- npar

  as.mcmc(samps)
}


cat2 <- function(x, sleep = 0.03) {
  if(interactive()) {
    x <- strsplit(x, "")[[1]]
    add <- ""
    for(i in seq_along(x)) {
      if(i > 1) {
        Sys.sleep(sleep)
        cat("\r")
      }
      cat(paste(x[1:i], collapse = ""))
    }
  } else cat(x)
  invisible(NULL)
}


uni.slice <- function(g, x, family, response, eta, id, j, ...,
  w = 1, m = 30, lower = -Inf, upper = +Inf, logPost)
{
  x0 <- g[j]
  gL <- gR <- g

  gx0 <- logPost(g, x, family, response, eta, id, ...)

  ## Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)

  ## Find the initial interval to sample from.
  ## w <- w * abs(x0) ## FIXME???
  u <- runif(1, 0, w)
  gL[j] <- g[j] - u
  gR[j] <- g[j] + (w - u)  ## should guarantee that g[j] is in [L, R], even with roundoff

  ## Expand the interval until its ends are outside the slice, or until
  ## the limit on steps is reached.
  if(is.infinite(m)) {
    repeat {
      if(gL[j] <= lower) break
      if(logPost(gL, x, family, response, eta, id, ...) <= logy) break
      gL[j] <- gL[j] - w
    }
    repeat {
      if(gR[j] >= upper) break
      if(logPost(gR, x, family, response, eta, id, ...) <= logy) break
      gR[j] <- gR[j] + w
    }
  } else {
    if(m > 1) {
      J <- floor(runif(1, 0, m))
      K <- (m - 1) - J
      while(J > 0) {
        if(gL[j] <= lower) break
        if(logPost(gL, x, family, response, eta, id, ...) <= logy) break
        gL[j] <- gL[j] - w
        J <- J - 1
      }
      while(K > 0) {
        if(gR[j] >= upper) break
        if(logPost(gR, x, family, response, eta, id, ...) <= logy) break
        gR[j] <- gR[j] + w
        K <- K - 1
      }
    }
  }

  ## Shrink interval to lower and upper bounds.
  if(gL[j] < lower) {
    gL[j] <- lower
  }
  if(gR[j] > upper) {
    gR[j] <- upper
  }

  ## Sample from the interval, shrinking it on each rejection.
  repeat {
    g[j] <- runif(1, gL[j], gR[j])

    gx1 <- logPost(g, x, family, response, eta, id, ...)

    if(gx1 >= logy) break

    if(g[j] > x0) {
      gR[j] <- g[j]
    } else {
      gL[j] <- g[j]
    }
  }

  ## Return the point sampled
  return(g)
}


uni.slice2 <- function(g, x, y, j, ...,
  w = 1, m = 30, lower = -Inf, upper = +Inf, logPost)
{
  x0 <- g[j]
  gL <- gR <- g

  gx0 <- logPost(g, x, y, ...)

  ## Determine the slice level, in log terms.
  logy <- gx0 - rexp(1)

  ## Find the initial interval to sample from.
  ## w <- w * abs(x0) ## FIXME???
  u <- runif(1, 0, w)
  gL[j] <- g[j] - u
  gR[j] <- g[j] + (w - u)  ## should guarantee that g[j] is in [L, R], even with roundoff

  ## Expand the interval until its ends are outside the slice, or until
  ## the limit on steps is reached.
  if(is.infinite(m)) {
    repeat {
      if(gL[j] <= lower) break
      if(logPost(gL, x, y, ...) <= logy) break
      gL[j] <- gL[j] - w
    }
    repeat {
      if(gR[j] >= upper) break
      if(logPost(gR, x, y, ...) <= logy) break
      gR[j] <- gR[j] + w
    }
  } else {
    if(m > 1) {
      J <- floor(runif(1, 0, m))
      K <- (m - 1) - J
      while(J > 0) {
        if(gL[j] <= lower) break
        if(logPost(gL, x, y, ...) <= logy) break
        gL[j] <- gL[j] - w
        J <- J - 1
      }
      while(K > 0) {
        if(gR[j] >= upper) break
        if(logPost(gR, x, y, ...) <= logy) break
        gR[j] <- gR[j] + w
        K <- K - 1
      }
    }
  }

  ## Shrink interval to lower and upper bounds.
  if(gL[j] < lower) {
    gL[j] <- lower
  }
  if(gR[j] > upper) {
    gR[j] <- upper
  }

  ## Sample from the interval, shrinking it on each rejection.
  repeat {
    g[j] <- runif(1, gL[j], gR[j])

    gx1 <- logPost(g, x, y, ...)

    if(gx1 >= logy) break

    if(g[j] > x0) {
      gR[j] <- g[j]
    } else {
      gL[j] <- g[j]
    }
  }

  ## Return the point sampled
  return(g)
}
