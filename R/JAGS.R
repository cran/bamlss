######################################################################
## (1) General JAGS setup functions, model code, initials and data. ##
######################################################################
## Setup the model structure needed for the sampler.
## These functions create the model code, initials and data
## to run a JAGS sampler.
## Examples: http://sourceforge.net/projects/mcmc-jags/files/
##           http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/Programs/
## Families: http://www.jrnold.me/categories/jags.html
## Default linear predictor and model setup functions.
BUGSeta <- function(x, id = NULL, ...)
{
  setup <- list()
  setup$inits <- setup$data <- list()

  ## Parametric part.
  if(FALSE) {
    id2 <- if(is.null(id)) 1 else id
    setup$data[[paste("X", id2, sep = "")]] <- x$X
    for(j in 1:k) {
      setup$param <- c(setup$param, paste(if(k < 2) paste("beta", id2, sep = "") else paste("beta", id2, "[", j, "]",
        sep = ""), "*X", id2, "[i, ", j, "]", sep = ""))
    }
    setup$param <- paste(setup$param, collapse = " + ")
    setup$param <- paste("    param", id2, "[i] <- ", setup$param, sep = "")
    setup$loops <- k
    setup$priors.coef <- if(k > 1) {
      paste("    beta", id2, "[j] ~ dnorm(0, 1.0E-6)", sep = "")
    } else paste("  beta", id2, " ~ dnorm(0, 1.0E-6)", sep = "")
    setup$inits[[paste("beta", id2, sep = "")]] <- runif(k)
    setup$psave <- c(setup$psave, paste("beta", id2, sep = ""))
    setup$eta <- paste("param", id2, "[i]", sep = "")
  }

  ## Cycle through all terms.
  if(m <- length(x$smooth.construct)) {
    for(i in 1:m) {
      setup <- if(!is.null(x$smooth.construct[[i]]$special)) {
        buildBUGS.smooth.special(x$smooth.construct[[i]], setup, paste(i, id, sep = ""))
      } else {
        buildBUGS.smooth(x$smooth.construct[[i]], setup, paste(i, id, sep = ""))
      }
    }
  }

  setup
}


## Get link functions.
BUGSlinks <- function(x)
{
  switch(x,
    "identity" = "eta",
    "log" = "exp(eta)",
    "exp" = "log(eta)",
    "inverse" = "1 / (eta^2)",
    "logit" = "1 / (1 + exp(-(eta)))",
    "probit" = "phi(eta)",
    "cloglog" = "log(-log(1 - eta))",
    "pow" = "pow(eta, -2)"
  )
}

## Construct the final model code.
BUGSmodel <- function(x, family, is.stan = FALSE, reference = NULL, ...) {
  if(is.function(family))
    family <- family()
  k <- if(all(c("inits", "data", "psave") %in% names(x))) {
    x <- list(x)
    1
  } else length(x)
  if(k > length(family$names)) {
    stop(paste("more parameters specified than existing in family ",
      family$family, ".bamlss()!", sep = ""), call. = FALSE)
  }
  model <- "model {"
  for(j in 1:k) {
    model <- c(model, x[[j]]$start)
  }
  pn <- family$names
  if(!is.null(family$bugs$reparam))
    pn[repi <- match(names(family$bugs$reparam), pn)] <- paste("rp", 1:length(family$bugs$reparam), sep = "")
  if(is.null(pn)) pn <- paste("theta", 1:k, sep = "")
  if(length(pn) < 2 & length(pn) != k)
    pn <- paste(pn, 1:k, sep = "")
  if(!is.null(family$bugs$order))
    pn <- pn[family$bugs$order]

  pn[1:k] <- paste(pn[1:k], "[i]", sep = "")
  on <- NULL
  links <- family[[grep("links", names(family), fixed = TRUE, value = TRUE)]]
  links <- rep(sapply(links, BUGSlinks), length.out = k)
  model2 <- c(
    paste("    y[i] ~ ", family$bugs$dist, "(",
      paste(if(!is.null(reference)) paste("probs[i, 1:", k + 1, "]", sep = "") else pn, collapse = ", "), ")", sep = ""))
  if(!is.null(reference)) {
    for(j in 1:k) {
      model2 <- c(model2,
        paste("    probs[i, ", j, "] <- ", pn[j], sep = "")
      )
    }
    model2 <- c(model2, paste("    probs[i, ", k + 1, "] <- ", 1, sep = ""))
  }
  if(!is.null(family$bugs$reparam)) {
    reparam <- NULL
    for(j in seq_along(family$bugs$reparam))
      reparam <- c(reparam, paste("    rp", j, "[i] <- ", family$bugs$reparam[j], sep = ""))
    for(j in family$names)
      reparam <- gsub(j, paste(j, "[i]", sep = ""), reparam)
    model2 <- c(model2, reparam)
    pn[repi] <- paste(family$names[repi], "[i]", sep = "")
  }
  if(!is.null(family$bugs$addparam)) {
    for(j in family$bugs$addparam)
      model2 <- c(model2, paste("   ", j))
  }
  if(!is.null(family$bugs$addvalues)) {
    for(j in names(family$bugs$addvalues))
      model2 <- gsub(j, family$bugs$addvalues[[j]], model2)
  }

  for(j in 1:k) {
    model2 <- c(model2, paste("    ", if(is.null(on)) pn[j] else paste(on, "[i, ", j, "]", sep = ""),
      " <- ", gsub("eta", x[[j]]$eta, links[[j]]), sep = ""))
  }
  for(j in 1:k)
    model2 <- c(model2, x[[j]]$adds)
  for(j in 1:k)
    model2 <- c(model2, x[[j]]$param, x[[j]]$smooth.construct)
  model2 <- rev(model2)
  model <- c(model, "  for(i in 1:n) {", model2)
  model <- c(model, "  }")

  for(i in 1:k)
    model <- c(model, x[[i]]$close1)

  for(i in 1:k) {
    lp <- list()
    if(!is.null(x[[i]]$loops)) {
      for(j in 1:length(x[[i]]$loops))
        lp[[paste(x[[i]]$loops[j])]] <- c(lp[[paste(x[[i]]$loops[j])]], x[[i]]$priors.coef[j])
    }
    if(length(lp)) {
      for(j in names(lp)) {
        if(j != 1)
          tmp <- c(paste("  for(j in 1:", j, ") {", sep = ""), lp[[j]], "  }")
        else
          tmp <- lp[[j]]
        model <- c(model, tmp)
      }
    }
    model <- c(model, x[[i]]$priors.scale, x[[i]]$close2, x[[i]]$close3)
  }
  model <- c(model, "}")

  model
}

## Create final model setup.
setupJAGS <- function(x, y, family, is.stan = FALSE)
{
  nx <- names(x)
  if(is.null(nx))
    nx <- 1:length(x)
  rval <- list()
  fn <- family$names
  if(length(fn) < length(x))
    fn <- paste(fn, 1:length(nx), sep = "")
  for(i in seq_along(nx)) {
    rval[[nx[i]]] <- family$bugs$eta(x[[i]], i)
  }

  y <- y[[1]]

  reference <- NULL
  if(is.factor(y)) {
    nl <- nlevels(y)
    y <- as.integer(y)
    if(nl < 3)
      y <- y - 1
  }
  if(is.matrix(y)) {
    if(length(nx) < ncol(y)) {
      cn <- colnames(y)
      reference <- cn[!(cn %in% nx)]
      y <- drop(apply(y, 1, function(x) {
        which(x != 0)
      }))
      y <- as.integer(y)
    }
  }

  ## Create model code.
  model <- family$bugs$model(rval, family, is.stan, reference)

  ## Collect data.
  if(all(c("inits", "data", "psave") %in% names(rval)))
    rval <- list(rval)
  data <- inits <- psave <- NULL
  for(j in seq_along(rval)) {
    data <- c(data, rval[[j]]$data)
    inits <- c(inits, rval[[j]]$inits)
    psave <- c(psave, rval[[j]]$psave)
  }
  data <- data[unique(names(data))]
  inits <- inits[unique(names(inits))]
  psave <- unique(psave)
  data$n <- length(y)
  psave <- c(psave, attr(model, "psave"))
  data$y <- y

  rval <- list("model" = model, "data" = data,
    "inits" = inits, "psave" = psave)
  
  return(rval)
}


## Build the JAGS model code for a smooth term. 
buildBUGS.smooth <- function(smooth, setup, i)
{
  fall <- NULL
  kr <- if(is.null(smooth$rand$Xr)) 0 else ncol(smooth$rand$Xr)
  kx <- if(is.null(smooth$Xf)) 0 else ncol(smooth$Xf)
  if(kr < 1 & kx < 1) {
    if(!is.null(smooth$X)) {
      smooth$Xf <- smooth$X
      kx <- ncol(smooth$Xf)
    }
  }
  hcauchy <- if(is.null(smooth$xt$hcauchy)) FALSE else smooth$xt$hcauchy
  ig <- itau2 <- NULL
  if(!is.null(smooth$state)) {
    ig <- get.par(smooth$state$parameters, "b")
    itau2 <- get.par(smooth$state$parameters, "tau2")
  }
  if(kx > 0) {
    fall <- c(fall, paste("b", i, if(kx > 1) paste("[", 1:kx, "]", sep = ""),
      "*Xf", i, "[i, ", 1:kx, "]", sep = ""))
    setup$data[[paste("Xf", i, sep = "")]] <- smooth$Xf
    tmp <- if(kx > 1) {
        paste("    b", i, "[j] ~ dnorm(0, 1.0E-6)", sep = "")
    } else paste("  b", i, " ~ dnorm(0, 1.0E-6)", sep = "")
    setup$priors.coef <- c(setup$priors.coef, tmp)
    setup$loops <- c(setup$loops, kx)
    setup$inits[[paste("b", i, sep = "")]] <- if(!is.null(ig)) {
      if(kr > 0)
        ig[-1 * (1:kr)]
      else
        ig
    } else runif(kx, 0.1, 0.2)
    setup$psave <- c(setup$psave, paste("b", i, sep = ""))
  }
  if(kr > 0) {
    fall <- c(fall, paste("g", i, if(kr > 1) paste("[", 1:kr, "]", sep = ""), "*Xr",
      i, "[i, ", 1:kr, "]", sep = ""))
    setup$data[[paste("Xr", i, sep = "")]] <- smooth$rand$Xr
    taug <- paste("taug", if(is.null(smooth$id)) i else smooth$id, sep = "")
    tmp <- if(kr > 1) {
      paste("    g", i, paste("[j] ~ dnorm(0, ", taug, ")", sep = ""), sep = "")
    } else paste("g", i, paste(" ~ dnorm(0, ", taug, ")", sep = ""), sep = "")
    setup$priors.coef <- c(setup$priors.coef, tmp)
    setup$loops <- c(setup$loops, kr)
    if(is.null(setup$priors.scale) || !any(grepl(taug, setup$priors.scale))) {
      if(!hcauchy) {
        setup$priors.scale <- c(setup$priors.scale, paste("  ", taug,
          " ~ dgamma(1.0E-4, 1.0E-4)", sep = ""))
      } else {
        setup$priors.scale <- c(setup$priors.scale, paste("  ", taug,
          " <- abs(", taug, 0, ")", sep = ""),
          paste("  ", taug, 0, "~ dt(0, 10, 1)", sep = ""))
      }
      setup$inits[[taug]] <- if(!is.null(itau2)) {
        itau2[1]
      } else runif(1, 0.1, 0.2)
      setup$psave <- c(setup$psave, taug)
    }
    setup$inits[[paste("g", i, sep = "")]] <- if(!is.null(ig)) {
      ig[1:kr]
    } else runif(kr, 0.1, 0.2)
    setup$psave <- c(setup$psave, paste("g", i, sep = ""))
  }

  setup$smooth.construct <- c(setup$smooth.construct, paste("    sm", i, "[i] <- ",
    paste(fall, collapse = " + ", sep = ""), sep = ""))
  setup$eta <- paste(setup$eta, paste("sm", i, "[i]", sep = ""),
    sep = if(length(setup$eta)) " + " else "")

  setup
}


## For special terms, e.g. growth curves, this function
## builds the model code.
buildBUGS.smooth.special <- function(smooth, setup, i)
{
  UseMethod("buildBUGS.smooth.special")
}


## Default special model term builder.
buildBUGS.smooth.special.default <- function(smooth, setup, i)
{
  buildBUGS.smooth(smooth, setup, i)
}


## Special code builder for growth curve terms.
buildBUGS.smooth.special.gc.smooth <- function(smooth, setup, i, zero)
{
  center <- if(is.null(smooth$xt$center)) TRUE else smooth$xt$center
  pn <- paste("g", i, sep = "")

  setup$data[[paste("X", pn, sep = "")]] <- as.numeric(smooth$X)
  setup$inits[[paste(pn, sep = "")]] <- runif(3, 0.1, 0.2)
  setup$psave <- c(setup$psave, pn)

  setup$close2 <- c(setup$close2,
    "  for(j in 1:3) {",
    paste("    ", pn, "[j] ~ dnorm(0, 1.0E-6)", sep = ""),
    "  }"
  )

  fall <- paste(pn, "[1] / (1 + exp(", pn, "[2]) * (exp(",
    pn, "[3]) / (1 + exp(", pn, "[3])))^(X", pn, "[i]))", sep = "")
  ##fall <- paste(pn, "[1] * X[i]", sep = "")

  if(!center) {
    setup$smooth.construct <- c(setup$smooth.construct, paste("    sm", i, "[i] <- ",
      paste(fall, collapse = " + ", sep = ""), sep = ""))
  } else {
    setup$close1 <- c(setup$close1, paste("  sm", i,  1, " <- sm", i, 0, " - mean(sm", i, 0, ")", sep = ""))
    setup$close1 <- c(setup$close1,
      paste("  for(i in 1:n) {", sep = ""),
      paste("    sm", i, 0, "[i] <- ",
        paste(fall, collapse = " + ", sep = ""), sep = ""), "  }")
  }
  setup$eta <- paste(setup$eta, paste("sm", i, if(center) 1 else NULL, "[i]", sep = ""),
    sep = if(length(setup$eta)) " + " else "")

  setup
}


########################################
## (3) Interface to the JAGS sampler. ##
########################################
samplerJAGS <- function(x, tdir = NULL,
  n.chains = 1, n.adapt = 100,
  n.iter = 4000, thin = 2, burnin = 1000,
  seed = NULL, verbose = TRUE, set.inits = TRUE,
  save.all = FALSE, modules = NULL)
{
  ## Temporary directory handling.
  if(is.null(tdir)) {
    dir.create(tdir <- tempfile())
    on.exit(unlink(tdir))
  } else tdir <- path.expand(tdir)
  if(!file.exists(tdir))
    dir.create(tdir)

  ## Write the model code.
  writeLines(paste(x$model, collapse = "\n"), mfile <- file.path(tdir, "model.txt"))

  ## Set the seed of the random number generator.
  if(is.null(seed))
    seed <- floor(runif(n.chains) * .Machine$integer.max)
  inits <- rep(list(x$inits), n.chains)
  for(j in seq_along(inits)) {
    inits[[j]][[".RNG.name"]] <- "base::Super-Duper"
    inits[[j]][[".RNG.seed"]] <- seed[j]
  }

  ## Sampling.
  rjags::load.module("dic")
  if(!is.null(modules)) {
    for(m in modules)
      rjags::load.module(m)
  }
  
  if(verbose) writeLines(x$model)
  
  if(save.all) {
    mdata <- x$data
    vnames <- x$psave
    save(mdata, vnames, file = file.path(tdir, "msetup.rda"))
    writeLines(c(
        'library("rjags")',
        'load("msetup.rda")',
        'm <- jags.model("model.txt", data = mdata, inits = inits)',
        paste('samples <- coda.samples(m, variable.names = vnames, n.iter = ',
          n.iter, ', thin = ', thin, ')', sep = '')
      ), con = file.path(tdir, "model.R")
    )
  }
  
  if(set.inits) {
    jmodel <- rjags::jags.model(mfile, data = x$data, inits = inits,
      n.chains = n.chains, n.adapt = n.adapt)
  } else {
    jmodel <- rjags::jags.model(mfile, data = x$data,
      n.chains = n.chains, n.adapt = n.adapt)
  }
  jsamples <- rjags::coda.samples(jmodel, variable.names = c(x$psave, "deviance"),
    n.iter = n.iter, thin = thin)

  ## Remove burnin.
  if(is.null(burnin))
    burnin <- floor(n.iter * 0.2)
  jsamples <- window(jsamples, start = burnin)

  jsamples
}


## Main JAGS function.
sam_JAGS <- JAGS <- function(x, y, family, start = NULL,
  tdir = NULL, n.chains = 1, n.adapt = 100,
  n.iter = 4000, thin = 2, burnin = 1000,
  seed = NULL, verbose = TRUE, set.inits = TRUE,
  save.all = FALSE, modules = NULL, ...)
{
  stopifnot(requireNamespace("rjags"))

  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, ...)
  if(!is.null(start))
    x <- set.starting.values(x, start)
  x <- randomize(x)
  ms <- setupJAGS(x, y, family)
  samps <- samplerJAGS(ms, tdir, n.chains, n.adapt,
    n.iter, thin, burnin, seed, verbose, set.inits,
    save.all, modules)
  return(transform.bugs.samples(x, samps))
}


transform.bugs.samples <- function(x, samples)
{
  nx <- names(x)
  snames <- colnames(samples[[1]])
  samps.out <- list()
  for(j in seq_along(samples)) {
    rval <- NULL
    for(i in seq_along(x)) {
      if(has_pterms(x[[i]]$terms)) {
        id <- which(names(x[[i]]$smooth.construct) == "model.matrix")
        sn <- paste("b", id, i, sep = "")
        sn <- grep(sn, snames, value = TRUE)
        psamples <- as.matrix(samples[, grep(sn, snames), drop = FALSE])
        colnames(psamples) <- paste(nx[i], "p", colnames(x[[i]]$smooth.construct$model.matrix$X), sep = ".")
        rval <- cbind(rval, psamples)
      }
      if(has_sterms(x[[i]]$terms)) {
        id <- 1:length(x[[i]]$smooth.construct)
        if(has_pterms(x[[i]]$terms))
          id <- id[id != which(names(x[[i]]$smooth.construct) == "model.matrix")]
        for(ii in id) {
          xsamples <- rsamples <- NULL
          kr <- if(is.null(x[[i]]$smooth.construct[[ii]]$rand$Xr)) 0 else ncol(x[[i]]$smooth.construct[[ii]]$rand$Xr)
          kx <- if(is.null(x[[i]]$smooth.construct[[ii]]$Xf)) 0 else ncol(x[[i]]$smooth.construct[[ii]]$Xf)
          if(kx) {
            pn <- grep(paste("b", ii, i, sep = ""), snames, value = TRUE, fixed = TRUE)
            pn <- pn[!grepl(paste("tau", ii, i, sep = ""), pn)]
            xsamples <- as.matrix(samples[[j]][, snames %in% pn])
          }
          if(kr) {
            pn <- grep(paste("g", ii, i, sep = ""), snames, value = TRUE, fixed = TRUE)
            pn <- pn[!grepl(paste("taug", ii, i, sep = ""), pn)]
            rsamples <- as.matrix(samples[[j]][, snames %in% pn])
          }
          psamples <- cbind("ra" = rsamples, "fx" = xsamples)
  
          ## Retransform parameter samples.
          if(kr) {
            re_trans <- function(g) {
              g <- x[[i]]$smooth.construct[[ii]]$trans.D * g
              if(!is.null(x[[i]]$smooth.construct[[ii]]$trans.U))
                g <- x[[i]]$smooth.construct[[ii]]$trans.U %*% g
              g
            }
            psamples <- t(apply(psamples, 1, re_trans))
          }

          colnames(psamples) <- paste(nx[i], "s", x[[i]]$smooth.construct[[ii]]$label,
            paste("b", 1:ncol(x[[i]]$smooth.construct[[ii]]$X), sep = ""), sep = ".")

          taug <- paste("taug", ii, i, sep = "")
          if(taug %in% snames) {
            vsamples <- as.matrix(samples[[j]][, snames %in% taug])
            colnames(vsamples) <- paste(nx[i], "s", x[[i]]$smooth.construct[[ii]]$label,
              paste("tau2", 1:length(x[[i]]$smooth.construct[[ii]]$S), sep = ""), sep = ".")
            psamples <- cbind(psamples, vsamples)
          }

          rval <- cbind(rval, psamples)
        }
      }
    }
    samps.out[[j]] <- as.mcmc(rval)
  }

  as.mcmc.list(samps.out)
}

