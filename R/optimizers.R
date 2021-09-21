##################################################
## (1) Generic setup function for smooth terms. ##
##################################################
## For each s() term, an additional list() named "state" must be supplied,
## within the list of specifications returned by smooth.construct(). This list
## contains the following objects:
##
## - fitted.values: numeric, vector containing the current fitted values.
## - parameters: numeric, vector containing the current parameters.
##               Can also contain smoothing variances/parameters,
##               which should be named with "tau2", regression coefficients
##               should be named with "b1", "b2", "b3", ..., "b"k.
## - edf: numeric, the current degrees of freedom of the term.
## - do.optim: NULL or logical, if NULL then per default the backfitting
##             algorithm is allowed to optimize possible variance
##             parameters within one update step of a term.
## - interval: numeric, if optimization is allowed, specifies the min and max
##             of the search space for optimizing variances. This can also
##             be supplied with the xt argument of s().
## - grid: integer, the grid length for searching variance parameters, if needed.
##         Can also be supplied within the xt argument in s().
##
## 4 additional functions must be supplied using the provided backfitting functions.
##
## - get.mu(X, gamma): computes the fitted values.
## - prior(parameters): computes the log prior using the parameters.
## - edf(x, weights): computes degrees of freedom of the smooth term.
## - grad(score, parameters, ...): function that computes the gradient.
##
## NOTE: model.matrix effects are also added to the smooth term list with
##       appropriate penalty structure. The name of the object in the
##       list is "model.matrix", for later identifyng the pure model.matrix
##       modeled effects.
##
## bamlss.engine.setup() sets up the basic structure, i.e., adds
## possible model.matrix terms to the smooth term list in x, also
## adds model.matrix terms of a random effect presentation of smooth
## terms to the "model.matrix" term. It calls the generic function
## bamlss.engine.setup.smooth(), which adds additional parts to the
## state list, as this could vary for special terms. A default
## method is provided.
bamlss.engine.setup <- function(x, update = "iwls", propose = "iwlsC_gp",
  do.optim = NULL, df = NULL, parametric2smooth = TRUE, ...)
{
  if(!is.null(attr(x, "bamlss.engine.setup"))) return(x)
  
  foo <- function(x, id = NULL) {
    if(!any(c("formula", "fake.formula") %in% names(x))) {
      for(j in names(x))
        x[[j]] <- foo(x[[j]], id = c(id, j))
    } else {
      if(is.null(id)) id <- ""
      if(!is.null(dim(x$model.matrix)) & parametric2smooth) {
        if(nrow(x$model.matrix) > 0 & !is.na(mean(unlist(x$model.matrix), na.rm = TRUE))) {
          if(is.null(x$smooth.construct)) x$smooth.construct <- list()
          label <- if(is.null(colnames(x$model.matrix))) {
            paste("b", 1:ncol(x$model.matrix), sep = "", collapse = "+")
          } else paste(colnames(x$model.matrix), collapse = "+")
          x$smooth.construct <- c(list("model.matrix" = list(
            "X" = x$model.matrix,
            "S" = list(diag(0, ncol(x$model.matrix))),
            "rank" = ncol(x$model.matrix),
            "term" = label,
            "label" = label,
            "bs.dim" = ncol(x$model.matrix),
            "fixed" = TRUE,
            "is.model.matrix" = TRUE,
            "by" = "NA",
            "xt" = list("binning" = x$binning)
          )), x$smooth.construct)
          if(!is.null(attr(x$model.matrix, "binning"))) {
            x$smooth.construct[["model.matrix"]]$binning <- attr(x$model.matrix, "binning")
            x$smooth.construct[["model.matrix"]]$xt$binning <- TRUE
          }
          class(x$smooth.construct[["model.matrix"]]) <- c(class(x$smooth.construct[["model.matrix"]]),
                                                           "no.mgcv", "model.matrix")
          x$model.matrix <- NULL
        }
      }
      if(length(x$smooth.construct)) {
        for(j in seq_along(x$smooth.construct)) {
          x$smooth.construct[[j]] <- bamlss.engine.setup.smooth(x$smooth.construct[[j]], ...)
          tdf <- NULL
          if(!is.null(df)) {
            if(!is.null(names(df))) {
              if((x$smooth.construct[[j]]$label %in% names(df)))
                tdf <- df[x$smooth.construct[[j]]$label]
            } else tdf <- df[1]
          }
          if(is.null(list(...)$nodf))
            x$smooth.construct[[j]] <- assign.df(x$smooth.construct[[j]], tdf, do.part = TRUE)
          if(!is.null(x$smooth.construct[[j]]$xt[["update"]]))
            x$smooth.construct[[j]]$update <- x$smooth.construct[[j]]$xt[["update"]]
          if(is.null(x$smooth.construct[[j]]$update)) {
            if(is.character(update)) {
              if(!grepl("bfit_", update))
                update <- paste("bfit", update, sep = "_")
              update <- get(update)
            }
            if(is.null(x$smooth.construct[[j]]$is.model.matrix))
              x$smooth.construct[[j]]$update <- update
            else
              x$smooth.construct[[j]]$update <- bfit_iwls
          }
          if(is.null(x$smooth.construct[[j]]$propose)) {
            if(is.character(propose)) {
              if(!grepl("GMCMC", propose))
                propose <- paste("GMCMC", propose, sep = "_")
              propose <- get(propose)
            }
            x$smooth.construct[[j]]$propose <- propose
          }
          if(is.null(do.optim))
            x$smooth.construct[[j]]$state$do.optim <- TRUE
          else
            x$smooth.construct[[j]]$state$do.optim <- do.optim
          if(!is.null(x$smooth.construct[[j]]$rank))
            x$smooth.construct[[j]]$rank <- as.numeric(x$smooth.construct[[j]]$rank)
          if(!is.null(x$smooth.construct[[j]]$Xf)) {
            x$smooth.construct[[j]]$Xfcn <- paste(paste(paste(x$smooth.construct[[j]]$term, collapse = "."),
                                                        "Xf", sep = "."), 1:ncol(x$smooth.construct[[j]]$Xf), sep = ".")
            colnames(x$smooth.construct[[j]]$Xf) <- x$smooth.construct[[j]]$Xfcn
            if(is.null(x$smooth.construct[["model.matrix"]])) {
              label <- paste(x$smooth.construct[[j]]$Xfcn, collapse = "+")
              x$smooth.construct[["model.matrix"]] <- list(
                "X" = x$smooth.construct[[j]]$Xf,
                "S" = list(diag(0, ncol(x$Xf))),
                "rank" = ncol(x$smooth.construct[[j]]$Xf),
                "term" = label,
                "label" = label,
                "bs.dim" = ncol(x$smooth.construct[[j]]$Xf),
                "fixed" = TRUE,
                "is.model.matrix" = TRUE,
                "by" = "NA"
              )
              x$smooth.construct <- c(x$smooth.construct, "model.matrix")
            } else {
              x$smooth.construct[["model.matrix"]]$X <- cbind(x$smooth.construct[["model.matrix"]]$X, x$smooth.construct[[j]]$Xf)
              x$smooth.construct[["model.matrix"]]$S <- list(diag(0, ncol(x$smooth.construct[["model.matrix"]]$X)))
              x$smooth.construct[["model.matrix"]]$bs.dim <- list(diag(0, ncol(x$smooth.construct[["model.matrix"]]$X)))
            }
          }
        }
      }
    }
    if(length(x$smooth.construct)) {
      if("model.matrix" %in% names(x$smooth.construct)) {
        if(length(nsc <- names(x$smooth.construct)) > 1) {
          nsc <- c(nsc[nsc != "model.matrix"], "model.matrix")
          x$smooth.construct <- x$smooth.construct[nsc]
        }
      }
      if(any(is.nnet <- sapply(x$smooth.construct, function(z) inherits(z, "nnet.smooth")))) {
        nsc <- names(x$smooth.construct)
        nsc <- c(nsc[-which(is.nnet)], nsc[which(is.nnet)])
        x$smooth.construct <- x$smooth.construct[nsc]
      }
    }
    x
  }
  
  x <- foo(x)
  attr(x, "bamlss.engine.setup") <- TRUE
  
  x
}


## Generic additional setup function for smooth terms.
bamlss.engine.setup.smooth <- function(x, ...) {
  UseMethod("bamlss.engine.setup.smooth")
}

## Simple extractor function.
get.state <- function(x, what = NULL) {
  if(is.null(what)) return(x$state)
  if(what %in% c("par", "parameters")) {
    return(x$state$parameters)
  } else {
    if(what %in% c("tau2", "lambda")) {
      p <- x$state$parameters
      return(p[grep("tau2", names(p))])
    } else {
      if(what %in% "b") {
        p <- x$state$parameters
        return(p[!grepl("tau2", names(p)) & !grepl("edf", names(p)) & !grepl("lasso", names(p))])
      } else return(x$state[[what]])
    }
  }
}

get.par <- function(x, what = NULL) {
  if(is.null(what) | is.null(names(x))) return(x)
  if(what %in% c("tau2", "lambda")) {
    return(x[grep("tau2", names(x))])
  } else {
    if(what %in% "b") {
      return(x[!grepl("tau2", names(x)) & !grepl("edf", names(x)) & !grepl("lasso", names(x)) & !grepl("bw", names(x))])
    } else return(x[what])
  }
}

set.par <- function(x, replacement, what) {
  if(is.null(replacement))
    return(x)
  if(what %in% c("tau2", "lambda")) {
    x[grep("tau2", names(x))] <- replacement
  } else {
    if(what %in% "b") {
      if(as.integer(sum(!grepl("tau2", names(x)) & !grepl("edf", names(x)) & !grepl("lasso", names(x)) & !grepl("bw", names(x)))) != length(replacement)) {
        stop("here")
      }
      x[!grepl("tau2", names(x)) & !grepl("edf", names(x)) & !grepl("lasso", names(x)) & !grepl("bw", names(x))] <- replacement
    } else x[what] <- replacement
  }
  x
}

## The default method.
bamlss.engine.setup.smooth.default <- function(x, Matrix = FALSE, ...)
{
  if(inherits(x, "special"))
    return(x)
  if(!is.null(x$margin)) {
    x$xt <- c(x$xt, x$margin[[1]]$xt)
    x$xt <- x$xt[unique(names(x$xt))]
    x$fixed <- x$margin[[1]]$fixed
  }
  if(is.null(x$binning) & !is.null(x$xt[["binning"]])) {
    if(is.logical(x$xt[["binning"]])) {
      if(x$xt[["binning"]]) {
        x$binning <- match.index(x$X)
        x$binning$order <- order(x$binning$match.index)
        x$binning$sorted.index <- x$binning$match.index[x$binning$order]
        assign <- attr(x$X, "assign")
        x$X <- x$X[x$binning$nodups, , drop = FALSE]
        attr(x$X, "assign") <- assign
      }
    } else {
      x$binning <- match.index(x$X)
      x$binning$order <- order(x$binning$match.index)
      x$binning$sorted.index <- x$binning$match.index[x$binning$order]
      assign <- attr(x$X, "assign")
      x$X <- x$X[x$binning$nodups, , drop = FALSE]
      attr(x$X, "assign") <- assign
    }
  }
  if(!is.null(x$binning)) {
    if(nrow(x$X) != length(x$binning$nodups)) {
      assign <- attr(x$X, "assign")
      x$X <- x$X[x$binning$nodups, , drop = FALSE]
      attr(x$X, "assign") <- assign
    }
  }
  if(is.null(x$binning)) {
    nr <- nrow(x$X)
    x$binning <- list(
      "match.index" = 1:nr,
      "nodups" = 1:nr,
      "order" = 1:nr,
      "sorted.index" = 1:nr
    )
  }
  x$nobs <- length(x$binning$match.index)
  k <- length(x$binning$nodups)
  x$weights <- rep(0, length = k)
  x$rres <- rep(0, length = k)
  x$fit.reduced <- rep(0, length = k)
  state <- if(is.null(x$xt[["state"]])) list() else x$xt[["state"]]
  if(is.null(x$fixed))
    x$fixed <- if(!is.null(x$fx)) x$fx[1] else FALSE
  if(!x$fixed & is.null(state$interval))
    state$interval <- if(is.null(x$xt[["interval"]])) tau2interval(x) else x$xt[["interval"]]
  if(!is.null(x$xt[["pSa"]])) {
    x$S <- c(x$S, list("pSa" = x$xt[["pSa"]]))
    priors <- make.prior(x)
    x$prior <- priors$prior
    x$grad <- priors$grad
    x$hess <- priors$hess
  }
  ntau2 <- length(x$S)
  if(length(ntau2) < 1) {
    if(x$fixed) {
      x$sp <- 1e+20
      ntau2 <- 1
      x$S <- list(diag(ncol(x$X)))
    } else {
      x$sp <- NULL
    }
  }
  if(!is.null(x$xt[["sp"]])) {
    x$sp <- x$xt[["sp"]]
    for(j in seq_along(x$sp))
      if(x$sp[j] == 0) x$sp[j] <- .Machine$double.eps^0.5
      x$xt[["tau2"]] <- 1 / x$sp
  }
  if(!is.null(x$sp)) {
    if(all(is.numeric(x$sp))) {
      x$sp <- rep(x$sp, length.out = ntau2)
      for(j in seq_along(x$sp))
        if(x$sp[j] == 0) x$sp[j] <- .Machine$double.eps^0.5
        x$fxsp <- TRUE
    } else x$fxsp <- FALSE
  } else x$fxsp <- FALSE
  if(is.null(state$parameters)) {
    state$parameters <- rep(0, ncol(x$X))
    names(state$parameters) <- if(is.null(colnames(x$X))) {
      paste("b", 1:length(state$parameters), sep = "")
    } else colnames(x$X)
    if(is.null(x$is.model.matrix)) {
      if(ntau2 > 0) {
        tau2 <- if(is.null(x$sp)) {
          if(x$fixed) {
            rep(1e+20, length.out = ntau2)
          } else {
            rep(if(!is.null(x$xt[["tau2"]])) {
              x$xt[["tau2"]]
            } else {
              if(!is.null(x$xt[["lambda"]])) 1 / x$xt[["lambda"]] else 1000
            }, length.out = ntau2)
          }
        } else rep(x$sp, length.out = ntau2)
        names(tau2) <- paste("tau2", 1:ntau2, sep = "")
        state$parameters <- c(state$parameters, tau2)
      }
    }
  }
  if((ntau2 > 0) & !any(grepl("tau2", names(state$parameters))) & is.null(x$is.model.matrix)) {
    tau2 <- if(is.null(x$sp)) {
      if(x$fixed) {
        rep(1e+20, length.out = ntau2)
      } else {
        rep(if(!is.null(x$xt[["tau2"]])) {
          x$xt[["tau2"]]
        } else {
          if(!is.null(x$xt[["lambda"]])) 1 / x$xt[["lambda"]] else 100
        }, length.out = ntau2)
      }
    } else rep(x$sp, length.out = ntau2)
    names(tau2) <- paste("tau2", 1:ntau2, sep = "")
    state$parameters <- c(state$parameters, tau2)
  }
  x$a <- if(is.null(x$xt[["a"]])) 1e-04 else x$xt[["a"]]
  x$b <- if(is.null(x$xt[["b"]])) 1e-04 else x$xt[["b"]]
  if(is.null(x$edf)) {
    x$edf <- function(x) {
      tau2 <- get.state(x, "tau2")
      if(x$fixed | !length(tau2)) return(ncol(x$X))
      if(is.null(x$state$XX))
        x$state$XX <- crossprod(x$X)
      S <- 0
      for(j in seq_along(tau2))
        S <- S + 1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](get.state(x, "b")) else x$S[[j]]
      P <- matrix_inv(x$state$XX + S, index = x$sparse.setup)
      edf <- try(sum_diag(x$state$XX %*% P), silent = TRUE)
      if(inherits(edf, "try-error"))
        edf <- ncol(x$X)
      return(edf)
    }
  }
  ng <- length(get.par(state$parameters, "b"))
  x$lower <- c(rep(-Inf, ng),
               if(is.list(state$interval)) {
                 unlist(sapply(state$interval, function(x) { x[1] }))
               } else state$interval[1])
  x$upper <- c(rep(Inf, ng),
               if(is.list(state$interval)) {
                 unlist(sapply(state$interval, function(x) { x[2] }))
               } else state$interval[2])
  names(x$lower) <- names(x$upper) <- names(state$parameters)[1:length(x$upper)]
  if(!is.null(x$sp)) {
    if(length(x$sp) < 1)
      x$sp <- NULL
    if(is.logical(x$sp))
      x[["sp"]] <- NULL
  }
  state$interval <- NULL
  x$state <- state
  if(!is.null(x$xt[["do.optim"]]))
    x$state$do.optim <- x$xt[["do.optim"]]
  x$sparse.setup <- sparse.setup(x$X, S = x$S)
  x$added <- c("nobs", "weights", "rres", "state", "a", "b", "prior", "edf",
    "grad", "hess", "lower", "upper")
  
  args <- list(...)
  
  force.Matrix <- if(is.null(args$force.Matrix)) FALSE else args$force.Matrix
  if(!is.null(x$xt$force.Matrix))
    force.Matrix <- x$xt$force.Matrix
  if(!is.null(x$sparse.setup$crossprod)) {
    if((ncol(x$sparse.setup$crossprod) < ncol(x$X) * 0.5) & force.Matrix)
      Matrix <- TRUE
    if(Matrix) {
      x$X <- Matrix(x$X, sparse = TRUE)
      for(j in seq_along(x$S))
        x$S[[j]] <- Matrix(if(is.function(x$S[[j]])) x$S[[j]](c("b" = rep(0, attr(x$S[[j]], "npar")))) else x$S[[j]], sparse = TRUE)
      if(force.Matrix)
        x$update <- bfit_iwls_Matrix
      priors <- make.prior(x)
      x$prior <- priors$prior
      x$grad <- priors$grad
      x$hess <- priors$hess
    }
  }
  if(ntau2 > 0) {
    tau2 <- NULL
    if(length(x$margin)) {
      for(j in seq_along(x$margin)) {
        if(!is.null(x$margin[[j]]$xt$tau2))
          tau2 <- c(tau2, x$margin[[j]]$xt$tau2)
      }
    } else {
      if(!is.null(x$xt$tau2))
        tau2 <- x$xt$tau2
      if(!is.null(x$xt[["lambda"]])) {
        tau2 <- 1 / x$xt[["lambda"]]
      }
    }
    if(!is.null(tau2)) {
      tau2 <- rep(tau2, length.out = ntau2)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    }
  }
  pid <- !grepl("tau2", names(x$state$parameters)) & !grepl("edf", names(x$state$parameters))
  x$pid <- list("b" = which(pid), "tau2" = which(!pid))
  if(!length(x$pid$tau2))
    x$pid$tau2 <- NULL
  if(is.null(x$prior)) {
    if(!is.null(x$xt[["prior"]]))
      x$prior <- x$xt[["prior"]]
    if(is.null(x$prior) | !is.function(x$prior)) {
      priors <- make.prior(x)
      x$prior <- priors$prior
      x$grad <- priors$grad
      x$hess <- priors$hess
    }
  }

  x$fit.fun <- make.fit.fun(x)
  x$state$fitted.values <- x$fit.fun(x$X, get.par(x$state$parameters, "b"))
  x$state$edf <- x$edf(x)
  
  x
}


## Function to find tau2 interval according to the
## effective degrees of freedom
tau2interval <- function(x, lower = .Machine$double.eps^0.8, upper = 1e+10)
{
  if(length(x$S) < 2) {
    return(c(lower, upper))
  } else {
    return(rep(list(c(lower, upper)), length.out = length(x$S)))
  }
}


## Assign degrees of freedom.
assign.df <- function(x, df, do.part = FALSE, ret.tau2 = FALSE)
{
  if(inherits(x, "special"))
    return(x)
  if(!is.null(x$is.model.matrix)) {
    if(x$is.model.matrix)
      return(x)
  }
  if(is.null(x$S))
    return(x)
  tau2 <- get.par(x$state$parameters, "tau2")
  if(x$fixed | !length(tau2))
    return(x)
  if(!is.null(x$fxsp)) {
    if(x$fxsp)
      return(x)
  }
  if(!is.null(x$no.assign.df))
    return(x)
  df <- if(is.null(x$xt$df)) df else x$xt$df
  if(is.null(df)) {
    nc <- ncol(x$X)
    df <- ceiling(nc * 0.5)
  }
  if(df > ncol(x$X))
    df <- ncol(x$X)
  XX <- crossprod(x$X)
  if(length(tau2) > 1) {
    objfun <- function(tau2, ret.edf = FALSE) {
      S <- 0
      for(i in seq_along(x$S))
        S <- S + 1 / tau2[i] * (if(is.function(x$S[[i]])) x$S[[i]](c("b" = rep(0, attr(x$S[[i]], "npar")))) else x$S[[i]])
      edf <- sum_diag(XX %*% matrix_inv(XX + S, index = x$sparse.setup))
      if(ret.edf)
        return(edf)
      else
        return((df - edf)^2)
    }
    if(do.part) {
      opt <- tau2.optim(objfun, start = tau2, maxit = 1000, scale = 100,
        add = FALSE, force.stop = FALSE, eps = .Machine$double.eps^0.8)
      if(!inherits(opt, "try-error"))
        tau2 <- opt
    }
  } else {
    objfun <- function(tau2, ret.edf = FALSE) {
      edf <- sum_diag(XX %*% matrix_inv(XX + 1 / tau2 * (if(is.function(x$S[[1]])) {
        x$S[[1]](c("b" = rep(0, attr(x$S[[1]], "npar")), x$fixed.hyper))
      } else x$S[[1]]), index = x$sparse.setup))
      if(ret.edf)
        return(edf)
      else
        return((df - edf)^2)
    }
    tau2 <- tau2.optim(objfun, start = tau2, maxit = 1000, scale = 10,
      add = FALSE, force.stop = FALSE, eps = .Machine$double.eps^0.8)
    if(inherits(tau2, "try-error"))
      return(x)
  }
  if(ret.tau2)
    return(tau2)
  x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
  x$state$edf <- objfun(tau2, ret.edf = TRUE)
  return(x)
}


get.eta <- function(x, expand = TRUE)
{
  nx <- names(x)
  np <- length(nx)
  eta <- vector(mode = "list", length = np)
  names(eta) <- nx
  for(j in 1:np) {
    eta[[j]] <- 0
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      par <- x[[nx[j]]]$smooth.construct[[sj]]$state$parameters
      par <- if(!is.null(x[[nx[j]]]$smooth.construct[[sj]]$pid)) {
        par[x[[nx[j]]]$smooth.construct[[sj]]$pid$b]
      } else get.state(x[[nx[j]]]$smooth.construct[[sj]], "b")
      fit <- x[[nx[j]]]$smooth.construct[[sj]]$fit.fun(x[[nx[j]]]$smooth.construct[[sj]]$X, par, expand)
      eta[[j]] <- eta[[j]] + fit
    }
  }
  eta
}


ffdf_eval <- function(x, FUN)
{
  #  res <- NULL
  #  for(i in bamlss_chunk(x)) {
  #    res <- ffappend(res, FUN(x[i, ]))
  #  }
  #  res
  ## FIXME: ff support!
  FUN(x)
}

ffdf_eval_sh <- function(y, par, FUN)
{
  #  res <- NULL
  #  for(i in bamlss_chunk(y)) {
  #    tpar <- list()
  #    for(j in names(par))
  #      tpar[[j]] <- par[[j]][i]
  #    res <- ffappend(res, FUN(y[i, ], tpar))
  #  }
  #  res
  ## FIXME: ff support!
  FUN(y, par)
}

ff_eval <- function(x, FUN, lower = NULL, upper = NULL)
{
  #  res <- NULL
  #  for(i in bamlss_chunk(x)) {
  #    tres <- FUN(x[i])
  #    if(!is.null(lower)) {
  #      if(any(jj <- tres == lower[1]))
  #        tres[jj] <- lower[2]
  #    }
  #    if(!is.null(upper)) {
  #      if(any(jj <- tres == upper[1]))
  #        tres[jj] <- upper[2]
  #    }
  #    res <- ffappend(res, tres)
  #  }
  #  res
  ## FIXME: ff support!
  FUN(x)
}


## Initialze.
init.eta <- function(eta, y, family, nobs)
{
  if(is.null(family$initialize))
    return(eta)
  for(j in family$names) {
    if(!is.null(family$initialize[[j]])) {
      linkfun <- make.link2(family$links[j])$linkfun
      if(inherits(y, "ffdf")) {
        eta[[j]] <- ffdf_eval(y, function(x) { linkfun(family$initialize[[j]](x)) })
      } else {
        eta[[j]] <- linkfun(family$initialize[[j]](y))
        eta[[j]] <- rep(eta[[j]], length.out = nobs)
      }
    }
  }
  return(eta)
}


get.edf <- function(x, type = 1)
{
  nx <- names(x)
  np <- length(nx)
  edf <- 0
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      edf <- edf + if(type < 2) {
        x[[nx[j]]]$smooth.construct[[sj]]$edf(x[[nx[j]]]$smooth.construct[[sj]])
      } else x[[nx[j]]]$smooth.construct[[sj]]$state$edf
    }
  }
  edf
}

get.log.prior <- function(x, type = 1)
{
  nx <- names(x)
  np <- length(nx)
  lp <- 0
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      lp <- lp + if(type < 2) {
        x[[nx[j]]]$smooth.construct[[sj]]$prior(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters)
      } else x[[nx[j]]]$smooth.construct[[sj]]$state$log.prior
    }
  }
  lp
}

get.all.par <- function(x, drop = FALSE, list = TRUE)
{
  nx <- names(x)
  np <- length(nx)
  par <- list()
  for(i in nx) {
    par[[i]] <- list()
    if(!all(c("formula", "fake.formula") %in% names(x[[i]]))) {
      for(k in names(x[[i]])) {
        if(!is.null(x[[i]][[k]]$smooth.construct)) {
          par[[i]][[k]]$s <- list()
          for(j in names(x[[i]][[k]]$smooth.construct)) {
            if(j == "model.matrix") {
              par[[i]][[k]]$p <- x[[i]][[k]]$smooth.construct[[j]]$state$parameters
            } else {
              if(is.null(par[[i]][[k]]$s))
                par[[i]][[k]]$s <- list()
              par[[i]][[k]]$s[[j]] <- x[[i]][[k]]$smooth.construct[[j]]$state$parameters
              if(!is.null(edf <- x[[i]][[k]]$smooth.construct[[j]]$state$edf))
                par[[i]][[k]]$s[[j]] <- c(par[[i]][[k]]$s[[j]], "edf" = edf)
            }
          }
        }
      }
    } else {
      if(!is.null(x[[i]]$smooth.construct)) {
        for(j in names(x[[i]]$smooth.construct)) {
          if(j == "model.matrix") {
            par[[i]]$p <- x[[i]]$smooth.construct[[j]]$state$parameters
          } else {
            if(is.null(par[[i]]$s))
              par[[i]]$s <- list()
            par[[i]]$s[[j]] <- x[[i]]$smooth.construct[[j]]$state$parameters
            if(!is.null(edf <- x[[i]]$smooth.construct[[j]]$state$edf))
              par[[i]]$s[[j]] <- c(par[[i]]$s[[j]], "edf" = edf)
          }
        }
      }
    }
  }
  if(!list) {
    par <- unlist(par)
    if(drop) {
      for(j in c(".edf", ".tau2", ".alpha"))
        par <- par[!grepl(j, names(par), fixed = TRUE)]
    }
  }
  par
}


get.hessian <- function(x)
{
  npar <- names(get.all.par(x, list = FALSE, drop = TRUE))
  hessian <- list(); nh <- NULL
  for(i in names(x)) {
    for(j in names(x[[i]]$smooth.construct)) {
      pn <- if(j == "model.matrix") paste(i, "p", sep = ".") else paste(i, "s", j, sep = ".")
      if(is.null(x[[i]]$smooth.construct[[j]]$state$hessian))
        x[[i]]$smooth.construct[[j]]$state$hessian <- diag(1e-07, ncol(x[[i]]$smooth.construct[[j]]$X))
      hessian[[pn]] <- as.matrix(x[[i]]$smooth.construct[[j]]$state$hessian)
      if(is.null(colnames(hessian[[pn]]))) {
        cn <- colnames(x[[i]]$smooth.construct[[j]]$X)
        if(is.null(cn))
          cn <- paste("b", 1:ncol(x[[i]]$smooth.construct[[j]]$X), sep = "")
      } else cn <- colnames(hessian[[pn]])
      pn <- paste(pn, cn, sep = ".")
      nh <- c(nh, pn)
    }
  }
  hessian <- -1 * as.matrix(do.call("bdiag", hessian))
  rownames(hessian) <- colnames(hessian) <- nh
  hessian <- hessian[npar, npar]
  return(hessian)
}


## Formatting for printing.
fmt <- Vectorize(function(x, width = 8, digits = 2) {
  txt <- formatC(round(x, digits), format = "f", digits = digits, width = width)
  if(nchar(txt) > width) {
    txt <- strsplit(txt, "")[[1]]
    txt <- paste(txt[1:width], collapse = "", sep = "")
  }
  txt
})

opt_bfit <- bfit <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
  update = "iwls", criterion = c("AICc", "BIC", "AIC"),
  eps = .Machine$double.eps^0.25, maxit = 400,
  outer = NULL, inner = FALSE, mgcv = FALSE,
  verbose = TRUE, digits = 4, flush = TRUE, nu = TRUE, stop.nu = NULL, ...)
{
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("design construct names mismatch with family names!")
  
  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, update = update, ...)

  plot <- if(is.null(list(...)$plot)) {
    FALSE
  } else {
    list(...)$plot
  }
  
  criterion <- match.arg(criterion)
  np <- length(nx)
  
  if(!is.null(nu)) {
    if(!is.logical(nu)) {
      if(nu < 0)
        nu <- NULL
    }
  }
  
  no_ff <- !inherits(y, "ffdf")
  
  nobs <- nrow(y)
  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }

  if(!is.null(start))
    x <- set.starting.values(x, start)
  eta <- get.eta(x)

  if(!is.null(weights))
    weights <- as.data.frame(weights)
  if(!is.null(offset)) {
    offset <- as.data.frame(offset)
    for(j in nx) {
      if(!is.null(offset[[j]]))
        eta[[j]] <- eta[[j]] + offset[[j]]
    }
  } else {
    if(is.null(start))
      eta <- init.eta(eta, y, family, nobs)
  }
  
  ia <- if(flush) interactive() else FALSE

  if(is.null(outer)) {
    outer <- FALSE
    max_outer <- 1
  } else {
    max_outer <- if(outer) {
      maxit + 1
    } else {
      0
    }
  }
  if(mgcv) {
    outer <- TRUE
    inner <- TRUE
    max_outer <- maxit + 1
  }
  
  inner_bf <- if(!mgcv) {
    function(x, y, eta, family, edf, id, nu, logprior, ...) {
      eps0 <- eps + 1; iter <- 1
      while(eps0 > eps & iter < maxit) {
        eta0 <- eta
        for(sj in seq_along(x)) {
          ## Get updated parameters.
          p.state <- x[[sj]]$update(x[[sj]], family, y, eta, id, edf = edf, ...)
          
          if(!is.null(nu)) {
            lpost0 <- family$loglik(y, family$map2par(eta))## + logprior
            ##lp <- logprior - x[[sj]]$prior(x[[sj]]$state$parameters)
            
            eta2 <- eta
            eta2[[id]] <- eta2[[id]] - x[[sj]]$state$fitted.values
            
            b0 <- get.par(x[[sj]]$state$parameters, "b")
            b1 <- get.par(p.state$parameters, "b")
            
            objfun <- function(nu, diff = TRUE) {
              p.state$parameters <- set.par(p.state$parameters, nu * b1 + (1 - nu) * b0, "b")
              eta2[[id]] <- eta2[[id]] + x[[sj]]$fit.fun(x[[sj]]$X,
                                                         get.par(p.state$parameters, "b"))
              lp2 <- family$loglik(y, family$map2par(eta2)) ##+ lp + x[[sj]]$prior(p.state$parameters)
              if(diff) {
                return(-1 * (lp2 - lpost0))
              } else
                return(lp2)
            }
            
            lpost1 <- objfun(1, diff = FALSE)
            
            if(lpost1 < lpost0) {
              if(!is.numeric(nu)) {
                nuo <- optimize(f = objfun, interval = c(0, 1))$minimum
              } else {
                nuo <- nu
                while((objfun(nuo, diff = FALSE) < lpost0) & (.Machine$double.eps < nuo)) {
                  nuo <- nuo / 2
                }
              }

              if(!(objfun(nuo, diff = FALSE) < lpost0)) {
                p.state$parameters <- set.par(p.state$parameters, nuo * b1 + (1 - nuo) * b0, "b")
                p.state$fitted.values <- x[[sj]]$fit.fun(x[[sj]]$X,
                                                       get.par(p.state$parameters, "b"))
                eta2[[id]] <- eta2[[id]] + p.state$fitted.values
                lpost1 <- family$loglik(y, family$map2par(eta2)) ##+ lp + x[[sj]]$prior(p.state$parameters)
              } else {
                next
              }
            }
          }
          
          ## Compute equivalent degrees of freedom.
          edf <- edf - x[[sj]]$state$edf + p.state$edf
          
          ## Update log priors.
          logprior <- logprior - x[[sj]]$prior(x[[sj]]$state$parameters) + x[[sj]]$prior(p.state$parameters)
          
          ## Update predictor and smooth fit.
          eta[[id]] <- eta[[id]] - fitted(x[[sj]]$state) + fitted(p.state)
          
          x[[sj]]$state <- p.state
        }
        eps0 <- do.call("cbind", eta)
        eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
        if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1
        iter <- iter + 1
      }
      return(list("x" = x, "eta" = eta, "edf" = edf))
    }
  } else {
    function(x, y, eta, family, edf, id, z, hess, weights, ...) {
      X <- lapply(x, function(x) { x$X })
      S <- lapply(x, function(x) { x$S })
      nt <- nt0 <- names(X)
      nt <- rmf(nt)
      names(X) <- names(S) <- nt
      if("modelmatrix" %in% nt)
        S <- S[!(nt %in% "modelmatrix")]
      X$z <- z
      f <- paste("z", paste(c("-1", nt), collapse = " + "), sep = " ~ ")
      f <- as.formula(f)
      if(!is.null(weights))
        hess <- hess * weights
      b <- gam(f, data = X, weights = hess, paraPen = S)
      cb <- coef(b)
      ncb <- names(cb)
      tau2 <- if(length(b$sp)) 1 / b$sp else NULL
      fitted <- 0
      for(sj in seq_along(x)) {
        tn <- rmf(nt0[sj])
        par <- cb[grep(tn, ncb, fixed = TRUE)]
        tedf <- sum(b$edf[grep(tn, ncb, fixed = TRUE)])
        names(par) <- paste("b", 1:length(par), sep = "")
        if(!is.null(tau2) & (tn != "modelmatrix")) {
          ttau2 <- tau2[grep(tn, names(tau2), fixed = TRUE)]
          names(ttau2) <- paste("tau2", 1:length(ttau2), sep = "")
          lo <- x[[sj]]$lower[grep("tau2", names(x[[sj]]$lower), fixed = TRUE)]
          up <- x[[sj]]$upper[grep("tau2", names(x[[sj]]$upper), fixed = TRUE)]
          if(any(j <- ttau2 < lo))
            ttau2[j] <- lo[j]
          if(any(j <- ttau2 > up))
            ttau2[j] <- up[j]
          par <- c(par, ttau2)
        } else {
          names(par) <- colnames(x[[sj]]$X)
          par <- c(par, "tau21" = 1e+20)
        }
        x[[sj]]$state$parameters <- par
        x[[sj]]$state$fitted.values <- x[[sj]]$fit.fun(x[[sj]]$X, par)
        fitted <- fitted + x[[sj]]$state$fitted.values
        edf <- edf - x[[sj]]$state$edf + tedf
        x[[sj]]$state$edf <- tedf
        x[[sj]]$state$prior <- x[[sj]]$prior(par)
      }
      eta[[id]] <- fitted
      return(list("x" = x, "eta" = eta, "edf" = edf))
    }
  }
  
  ## Backfitting main function.
  backfit <- function(verbose = TRUE) {
    eps0 <- eps + 1; iter <- 0
    edf <- get.edf(x, type = 2)
    ll_save <- NULL
    ptm <- proc.time()
    while(eps0 > eps & iter < maxit) {
      eta0 <- eta
      ## Cycle through all parameters
      for(j in 1:np) {
        if(iter < max_outer) {
          peta <- family$map2par(eta)
          
          if(no_ff) {
            ## Compute weights.
            hess <- process.derivs(family$hess[[nx[j]]](y, peta, id = nx[j]), is.weight = TRUE)
            
            ## Score.
            score <- process.derivs(family$score[[nx[j]]](y, peta, id = nx[j]), is.weight = FALSE)
            if(length(score) != nobs) { 
              stop("something wrong in processing the family $score() function! More elements in return value of $score() than the response!")
            }
            if(length(hess) != nobs) { 
              stop("something wrong in processing the family $hess() function! More elements in return value of $hess() than the response!")
            }
          } else {
            ## Same for large files.
            hess <- ffdf_eval_sh(y, peta, FUN = function(y, par) {
              process.derivs(family$hess[[nx[j]]](y, par, id = nx[j]), is.weight = TRUE)
            })
            
            score <- ffdf_eval_sh(y, peta, FUN = function(y, par) {
              process.derivs(family$score[[nx[j]]](y, par, id = nx[j]), is.weight = FALSE)
            })
          }
          
          ## Compute working observations.
          z <- eta[[nx[j]]] + 1 / hess * score

        } else z <- hess <- score <- NULL
        
        if(iter < 2) {
          eta[[nx[j]]] <- get.eta(x)[[nx[j]]]
          if(!is.null(offset)) {
            if(!is.null(offset[[nx[j]]]))
              eta[[nx[j]]] <- eta[[nx[j]]] + offset[[nx[j]]]
          }
        }
        
        ## And all terms.
        if(inner) {
          tbf <- inner_bf(x[[nx[j]]]$smooth.construct, y, eta, family,
            edf = edf, id = nx[j], z = z, hess = hess, weights = weights[[nx[j]]],
            criterion = criterion, iteration = iter, nu = nu, score = score,
            logprior = get.log.prior(x))
          x[[nx[j]]]$smooth.construct <- tbf$x
          edf <- tbf$edf
          eta <- tbf$eta
          rm(tbf)
        } else {
          for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
            ## Get updated parameters.
            p.state <- x[[nx[j]]]$smooth.construct[[sj]]$update(x[[nx[j]]]$smooth.construct[[sj]],
              family, y, eta, nx[j], edf = edf, z = z, hess = hess, weights = weights[[nx[j]]],
              iteration = iter, criterion = criterion, score = score)
            
            ## Update predictor and smooth fit.
            if(!is.null(nu) & !inherits(x[[nx[j]]]$smooth.construct[[sj]], "nnet0.smooth")) {
              lp0 <- get.log.prior(x)
              lpost0 <- family$loglik(y, family$map2par(eta)) ##+ lp0
              ##lp <- lp0 - x[[nx[j]]]$smooth.construct[[sj]]$prior(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters)
              
              eta2 <- eta
              eta2[[nx[j]]] <- eta2[[nx[j]]] - x[[nx[j]]]$smooth.construct[[sj]]$state$fitted.values
              
              b0 <- get.par(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, "b")
              b1 <- get.par(p.state$parameters, "b")
              
              objfun <- function(nu, diff = TRUE) {
                p.state$parameters <- set.par(p.state$parameters, nu * b1 + (1 - nu) * b0, "b")
                eta2[[nx[j]]] <- eta2[[nx[j]]] + x[[nx[j]]]$smooth.construct[[sj]]$fit.fun(x[[nx[j]]]$smooth.construct[[sj]]$X,
                                                                                           get.par(p.state$parameters, "b"))
                lp2 <- family$loglik(y, family$map2par(eta2)) ##+ lp + x[[nx[j]]]$smooth.construct[[sj]]$prior(p.state$parameters)
                if(diff) {
                  return(-1 * (lp2 - lpost0))
                } else
                  return(lp2)
              }
              
              lpost1 <- objfun(1, diff = FALSE)
              
              if(lpost1 < lpost0) {
                if(!is.numeric(nu)) {
                  nuo <- optimize(f = objfun, interval = c(0, 1))$minimum
                } else {
                  nuo <- nu
                  while((objfun(nuo, diff = FALSE) < lpost0) & (.Machine$double.eps < nuo)) {
                    nuo <- nuo / 2
                  }
                }

                if(!(objfun(nuo, diff = FALSE) < lpost0)) {
                  p.state$parameters <- set.par(p.state$parameters, nuo * b1 + (1 - nuo) * b0, "b")
                  p.state$fitted.values <- x[[nx[j]]]$smooth.construct[[sj]]$fit.fun(x[[nx[j]]]$smooth.construct[[sj]]$X, get.par(p.state$parameters, "b"))
                  eta2[[nx[j]]] <- eta2[[nx[j]]] + p.state$fitted.values
                  lpost1 <- family$loglik(y, family$map2par(eta2)) ##+ lp + x[[nx[j]]]
                } else {
                  next
                }
              }
            }
            
            ## Compute equivalent degrees of freedom.
            edf <- edf - x[[nx[j]]]$smooth.construct[[sj]]$state$edf + p.state$edf
            
            eta[[nx[j]]] <- eta[[nx[j]]] - fitted(x[[nx[j]]]$smooth.construct[[sj]]$state) + fitted(p.state)
            
            x[[nx[j]]]$smooth.construct[[sj]]$state <- p.state
          }
        }
      }
      
      if(!is.null(stop.nu)) {
        if(iter > stop.nu)
          nu <- NULL
      }

      eps0 <- do.call("cbind", eta)
      eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
      if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1
      
      peta <- family$map2par(eta)
      
      IC <- get.ic(family, y, peta, edf, nobs, criterion)
      
      iter <- iter + 1

      logLik <- family$loglik(y, peta)
      
      if(verbose) {
        cat(if(ia) "\r" else if(iter > 1) "\n" else NULL)
        vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
                      " logPost ", fmt(family$loglik(y, peta) + get.log.prior(x), width = 8, digits = digits),
                      " logLik ", fmt(logLik, width = 8, digits = digits),
                      " edf ", fmt(edf, width = 6, digits = digits),
                      " eps ", fmt(eps0, width = 6, digits = digits + 2),
                      " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
        cat(vtxt)
        
        if(.Platform$OS.type != "unix" & ia) flush.console()
      }

      ll_save <- c(ll_save, logLik)
      if(iter > 2)
        slope <- ll_save[length(ll_save)] - ll_save[length(ll_save) - 1L]
      else
        slope <- NA

      if(plot) {
        plot(ll_save, xlab = "Iteration", ylab = "logLik",
          main = paste("Slope", fmt(slope, width = 6, digits = digits + 2)),
          type = "l", ylim = c(0.9 * max(ll_save), max(ll_save)))
      }
    }
    
    elapsed <- c(proc.time() - ptm)[3]
    
    IC <- get.ic(family, y, peta, edf, nobs, criterion)
    logLik <- family$loglik(y, peta)
    logPost <- as.numeric(logLik + get.log.prior(x))

    ll_save <- c(ll_save, logLik)
    if(iter > 2)
      slope <- ll_save[length(ll_save)] - ll_save[length(ll_save) - 1L]
    else
      slope <- NA
    
    if(verbose) {
      cat(if(ia) "\r" else "\n")
      vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
                    " logPost ", fmt(logPost, width = 8, digits = digits),
                    " logLik ", fmt(family$loglik(y, peta), width = 8, digits = digits),
                    " edf ", fmt(edf, width = 6, digits = digits),
                    " eps ", fmt(eps0, width = 6, digits = digits + 2),
                    " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
      cat(vtxt)
      if(.Platform$OS.type != "unix" & ia) flush.console()
      et <- if(elapsed > 60) {
        paste(formatC(format(round(elapsed / 60, 2), nsmall = 2), width = 5), "min", sep = "")
      } else paste(formatC(format(round(elapsed, 2), nsmall = 2), width = 5), "sec", sep = "")
      cat("\nelapsed time: ", et, "\n", sep = "")
    }

    if(plot) {
      plot(ll_save, xlab = "Iteration", ylab = "logLik",
        main = paste("Slope", fmt(slope, width = 6, digits = digits + 2)),
        type = "l", ylim = c(0.9 * max(ll_save), max(ll_save)))
    }
    
    if(iter == maxit)
      warning("the backfitting algorithm did not converge!")
    
    names(IC) <- criterion
    
    rval <- list("fitted.values" = eta, "parameters" = get.all.par(x), "edf" = edf,
                 "logLik" = logLik, "logPost" = logPost, "nobs" = nobs,
                 "converged" = iter < maxit, "runtime" = elapsed)
    rval[[names(IC)]] <- IC
    rval
  }
  
  backfit(verbose = verbose)
}


## Extract information criteria.
get.ic <- function(family, y, par, edf, n, type = c("AIC", "BIC", "AICc", "MP"), ...)
{
  type <- match.arg(type)
  ll <- family$loglik(y, par)
  if(is.na(edf))
    edf <- n - 1
  denom <- (n - edf - 1)
  if(is.na(denom)) {
    add <- 0
  } else {
    if(denom < 1e-10) {
      add <- 0
    } else {
      add <- (2 * edf * (edf + 1)) / denom
    }
  }
  pen <- switch(type,
                "AIC" = -2 * ll + 2 * edf,
                "BIC" = -2 * ll + edf * log(n),
                "AICc" = -2 * ll + 2 * edf + add,
                "MP" = -1 * (ll + edf)
  )
  return(pen)
}

get.ic2 <- function(logLik, edf, n, type = c("AIC", "BIC", "AICc", "MP"), ...)
{
  type <- match.arg(type)
  denom <- (n - edf - 1)
  if(denom < 1e-10) {
    add <- 0
  } else {
    add <- (2 * edf * (edf + 1)) / denom
  }
  pen <- switch(type,
                "AIC" = -2 * logLik + 2 * edf,
                "BIC" = -2 * logLik + edf * log(n),
                "AICc" = -2 * logLik + 2 * edf + add,
                "MP" = -1 * (logLik + edf)
  )
  return(pen)
}


cround <- function(x, digits = 2)
{
  return(x)
  cdigits <- Vectorize(function(x) {
    if(abs(x) >= 1)
      return(0)
    scipen <- getOption("scipen")
    on.exit(options("scipen" = scipen))
    options("scipen" = 100)
    x <- strsplit(paste(x), "")[[1]]
    x <- x[which(x == "."):length(x)][-1]
    i <- which(x != "0")
    x <- x[1:(i[1] - 1)]
    n <- length(x)
    if(n < 2) {
      if(x != "0")
        return(1)
      else return(n + 1)
    } else return(n + 1)
  })
  
  round(x, digits = cdigits(x) + digits)
}


## Naive smoothing parameter optimization.
tau2.optim <- function(f, start, ..., scale = 10, eps = .Machine$double.eps^0.5, maxit = 1, add = TRUE, force.stop = TRUE, optim = FALSE)
{
  if(optim) {
    lower <- start / scale
    upper <- start * scale + if(add) 1 else 0
    start <- optim(start, fn = f, method = "L-BFGS-B",
      lower = lower, upper = upper)$par
    return(start)
  }

  foo <- function(par, start, k) {
    start[k] <- cround(par)
    return(f(start, ...))
  }
  
  start <- cround(start)
  ic0 <- f(start, ...)
  
  iter <- 0; eps0 <- eps + 1
  while((eps0 > eps) & (iter < maxit)) {
    start0 <- start
    for(k in seq_along(start)) {
      xr <- c(start[k] / scale, start[k] * scale + if(add) 1 else 0)
      tpar <- try(optimize(foo, interval = xr, start = start, k = k, tol = eps), silent = TRUE)
      if(!inherits(tpar, "try-error")) {
        if(tpar$objective < ic0) {
          start[k] <- tpar$minimum
          ic0 <- tpar$objective
        }
      }
    }

    if((length(start) < 2) & force.stop)
      break
    
    eps0 <- mean(abs((start - start0) / start0))

    iter <- iter + 1
  }
  
  return(start)
}


## Function to create full parameter vector.
make_par <- function(x, type = 1, add.tau2 = FALSE) {
  family <- attr(x, "family")
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  np <- length(nx)
  par <- lower <- upper <- NULL
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      tpar <- x[[nx[j]]]$smooth.construct[[sj]]$state$parameters
      tlower <- x[[nx[j]]]$smooth.construct[[sj]]$lower
      tupper <- x[[nx[j]]]$smooth.construct[[sj]]$upper
      if(!add.tau2) {
        tlower <- tlower[!grepl("tau2", names(tlower))]
        tupper <- tupper[!grepl("tau2", names(tupper))]
        tpar <- tpar[!grepl("tau2", names(tpar))]
      }
      g <- get.par(tpar, "b")
      npar <- paste(paste(nx[j], "h1", x[[nx[j]]]$smooth.construct[[sj]]$label, sep = ":"), 1:length(g), sep = ".")
      if(length(tau2 <- get.par(tpar, "tau2"))) {
        npar <- c(npar, paste(nx[j], "h1", paste(x[[nx[j]]]$smooth.construct[[sj]]$label,
                                                 paste("tau2", 1:length(tau2), sep = ""), sep = "."), sep = ":"))
      }
      names(tpar) <- names(tlower) <- names(tupper) <- if(type < 2) {
        paste("p", j, ".t", sj, ".", names(tpar), sep = "")
      } else npar
      par <- c(par, tpar)
      lower <- c(lower, tlower)
      upper <- c(upper, tupper)
    }
  }
  return(list("par" = par, "lower" = lower, "upper" = upper))
}


## Backfitting updating functions.
bfit_newton <- function(x, family, y, eta, id, ...)
{
  args <- list(...)
  
  eta[[id]] <- eta[[id]] - fitted(x$state)
  
  tau2 <- if(!x$fixed) get.par(x$state$parameters, "tau2") else NULL
  
  lp <- function(g) {
    eta[[id]] <- eta[[id]] + x$fit.fun(x$X, g)
    family$loglik(y, family$map2par(eta)) + x$prior(c(g, tau2))
  }
  
  if(is.null(family$gradient[[id]])) {
    gfun <- NULL
  } else {
    gfun <- list()
    gfun[[id]] <- function(g, y, eta, x, ...) {
      gg <- family$gradient[[id]](g, y, eta, x, ...)
      if(!is.null(x$grad)) {
        gg <- gg + x$grad(score = NULL, c(g, tau2), full = FALSE)
      }
      drop(gg)
    }
  }
  
  if(is.null(family$hessian[[id]])) {
    hfun <- NULL
  } else {
    hfun <- list()
    hfun[[id]] <- function(g, y, eta, x, ...) {
      hg <- family$hessian[[id]](g, y, eta, x, ...)
      if(!is.null(x$hess)) {
        hg <- hg + x$hess(score = NULL, c(g, tau2), full = FALSE)
      }
      hg
    }
  }
  
  g <- get.par(x$state$parameters, "b")
  nu <- if(is.null(x$nu)) 0.1 else x$nu
  
  g.grad <- grad(fun = lp, theta = g, id = id, prior = NULL,
                 args = list("gradient" = gfun, "x" = x, "y" = y, "eta" = eta))
  
  g.hess <- hess(fun = lp, theta = g, id = id, prior = NULL,
                 args = list("gradient" = gfun, "hessian" = hfun, "x" = x, "y" = y, "eta" = eta))
  
  Sigma <- matrix_inv(g.hess, index = x$sparse.setup)
  
  g <- drop(g + nu * Sigma %*% g.grad)
  
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))
  x$state$hessian <- Sigma
  
  return(x$state)
}


boostm_fit0 <- function(x, grad, hess, z, nu, stop.criterion, family, y, eta, edf, id, do.optim, ...)
{
  b0 <- get.par(x$state$parameters, "b")

  xbin.fun(x$binning$sorted.index, hess, rep(0, length(grad)),
    x$weights, x$rres, x$binning$order, x$binning$uind)
  XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)

  if(x$fixed) {
    k <- length(b0)
    S <- matrix(0, k, k)
  } else {
    if(do.optim) {
      tpar <- x$state$parameters
      edf <- edf - x$state$edf

      eta[[id]] <- eta[[id]] - fitted(x$state)

      objfun <- function(tau2) {
        tpar2 <- set.par(tpar, tau2, "tau2")
        S <- 0
        for(j in seq_along(tau2))
          S <- S + (1 / tau2[j]) * if(is.function(x$S[[j]])) x$S[[j]](c(tpar2, x$fixed.hyper)) else x$S[[j]]
        xgrad <- t(x$X) %*% grad + S %*% b0
        xhess <- XWX + S
        Sigma <- matrix_inv(xhess, index = x$sparse.setup)
        b1 <- b0 + drop(nu * Sigma %*% xgrad)
        eta[[id]] <- eta[[id]] + x$fit.fun(x$X, b1)
        edf <- edf + sum_diag(XWX %*% Sigma)
        return(get.ic(family, y, family$map2par(eta), edf, length(eta[[1]]), type = stop.criterion))
      }

      tau2 <- tau2.optim(objfun, start = get.par(x$state$parameters, "tau2"), scale = 10, maxit = 10, force.stop = FALSE, add = FALSE)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    }

    tau2 <- get.par(x$state$parameters, "tau2")

    S <- 0
    for(j in seq_along(tau2))
      S <- S + (1 / tau2[j]) * if(is.function(x$S[[j]])) x$S[[j]](c(x$state$parameters, x$fixed.hyper)) else x$S[[j]]
  }

  xgrad <- t(x$X) %*% grad + S %*% b0
  xhess <- XWX + S
  Sigma <- matrix_inv(xhess, index = x$sparse.setup)
  b1 <- b0 + drop(nu * Sigma %*% xgrad)
    
  x$state$parameters <- set.par(x$state$parameters, b1, "b")
  x$state$fitted.values <- x$fit.fun(x$X, b1)
  x$state$hessian <- Sigma
  x$state$edf <- sum_diag(XWX %*% Sigma)
  
  return(x$state)
}


boostm_fit <- function(x, grad, hess, z, nu, stop.criterion, family, y, eta, edf, id, do.optim, ...)
{
  x$state$do.optim <- do.optim
  b0 <- get.par(x$state$parameters, "b")
  state <- bfit_iwls(x = x, family = family, y = y, eta = eta, id = id, criterion = stop.criterion,
    grad = grad, hess = hess, z = z, edf = edf, ...)
  b1 <- get.par(state$parameters, "b")
  state$parameters <- set.par(state$parameters, nu * b1 + (1 - nu) * b0, "b")
  state$fitted.values <- x$fit.fun(x$X, get.par(state$parameters, "b"))
  return(state)
}

#x smooth construct selbst x$X design matrix x$S als Funktion
bfit_lm <- function(x, family, y, eta, id, weights, criterion, ...)
{
  args <- list(...)
  
  peta <- family$map2par(eta)
  
  hess <- family$hess[[id]](y, peta, id = id, ...)
  
  ## Score.
  score <- family$score[[id]](y, peta, id = id, ...)
  
  ## Compute working observations.
  z <- eta[[id]] + 1 / hess * score
  
  ## Compute reduced residuals.
  e <- z - eta[[id]] + fitted(x$state)
  
  if(!is.null(weights))
    hess <- hess * weights
  if(x$fixed | x$fxsp) {
    b <- lm.wfit(x$X, e, hess)
  } else {
    tau2 <- get.par(x$state$parameters, "tau2")
    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    n <- nrow(S)
    w <- c(hess, rep(0, n))
    e <- c(e, rep(1, n))
    b <- lm.wfit(rbind(x$X, S), e, w)
  }
  
  x$state$parameters <- set.par(x$state$parameters, coef(b), "b")
  x$state$fitted.values <- x$X %*% coef(b)
  
  x$state
}


bfit_iwls <- function(x, family, y, eta, id, weights, criterion, ...)
{
  args <- list(...)

  no_ff <- !inherits(y, "ff")
  peta <- family$map2par(eta)

  nobs <- length(eta[[1L]])
  
  if(is.null(args$hess)) {
    ## Compute weights.
    if(no_ff) {
      hess <- process.derivs(family$hess[[id]](y, peta, id = id, ...), is.weight = TRUE)
    } else {
      hess <- ffdf_eval_sh(y, peta, FUN = function(y, par) {
        process.derivs(family$hess[[id]](y, par, id = id), is.weight = TRUE)
      })
    }

    if(length(hess) != nobs) { 
      stop("something wrong in processing the family $hess() function! More elements in return value of $hess() than the response!")
    }
  } else hess <- args$hess
  
  if(!is.null(weights))
    hess <- hess * weights
  
  if(is.null(args$z)) {
    ## Score.
    if(no_ff) {
      score <- process.derivs(family$score[[id]](y, peta, id = id, ...), is.weight = FALSE)
    } else {
      score <- ffdf_eval_sh(y, peta, FUN = function(y, par) {
        process.derivs(family$score[[id]](y, par, id = id), is.weight = FALSE)
      })
    }

    if(length(score) != nobs) { 
      stop("something wrong in processing the family $score() function! More elements in return value of $score() than the response!")
    }
    
    ## Compute working observations.
    z <- eta[[id]] + 1 / hess * score
  } else z <- args$z
  
  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - fitted(x$state)
  
  ## Compute reduced residuals.
  e <- z - eta[[id]]

  xbin.fun(x$binning$sorted.index, hess, e, x$weights, x$rres, x$binning$order, x$binning$uind)
  
  ## Old parameters.
  g0 <- get.state(x, "b")
  
  ## Compute mean and precision.
  XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)

  if(!x$state$do.optim | x$fixed | x$fxsp) {
    if(x$fixed) {
      P <- matrix_inv(XWX + if(!is.null(x$xt[["pS"]])) x$xt[["pS"]] else 0, index = x$sparse.setup)
    } else {
      S <- 0
      tau2 <- get.state(x, "tau2")
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](c(g0, x$fixed.hyper)) else x$S[[j]]
      P <- matrix_inv(XWX + S + if(!is.null(x$xt[["pS"]])) x$xt[["pS"]] else 0, index = x$sparse.setup)
    }
    if(is.null(x$xt[["pm"]])) {
      x$state$parameters <- set.par(x$state$parameters, drop(P %*% crossprod(x$X, x$rres)), "b")
    } else {
      pS <- if(!is.null(x$xt[["pS"]])) {
        x$xt[["pS"]]
      } else {
        if(!is.null(x$xt[["pSa"]])) {
          1 / tau2[length(tau2)] * x$xt[["pSa"]]
        } else 0
      }
      x$state$parameters <- set.par(x$state$parameters, drop(P %*% (crossprod(x$X, x$rres) + pS %*% x$xt[["pm"]])), "b")
    }
  } else {
    args <- list(...)
    edf0 <- args$edf - x$state$edf
    eta2 <- eta
    
    env <- new.env()

    objfun <- function(tau2, ...) {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](c(g0, x$fixed.hyper)) else x$S[[j]]
      P <- matrix_inv(XWX + S + if(!is.null(x$xt[["pS"]])) x$xt[["pS"]] else 0, index = x$sparse.setup)
      if(inherits(P, "try-error")) return(NA)
      if(is.null(x$xt[["pm"]])) {
        g <- drop(P %*% crossprod(x$X, x$rres))
      } else {
        pS <- if(!is.null(x$xt[["pS"]])) {
          x$xt[["pS"]]
        } else {
          if(!is.null(x$xt[["pSa"]])) {
            1 / tau2[length(tau2)] * x$xt[["pSa"]]
          } else 0
        }
        g <- drop(P %*% (crossprod(x$X, x$rres) + pS %*% x$xt[["pm"]]))
      }

      if(!is.null(x$doCmat)) {
        V <- P %*% t(x$C)
        W <- x$C %*% V
        U <- chol2inv(chol(W)) %*% t(V)
        g <- drop(g - t(U) %*% x$C %*% g)
      }

      if(any(is.na(g)) | any(g %in% c(-Inf, Inf))) g <- rep(0, length(g))
      fit <- x$fit.fun(x$X, g)

      if(!is.null(x$doCmat))
        fit <- fit - mean(fit, na.rm = TRUE)
      edf <- sum_diag(XWX %*% P)
      eta2[[id]] <- eta2[[id]] + fit
      ic <- get.ic(family, y, family$map2par(eta2), edf0 + edf, length(z), criterion, ...)
      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          par <- c(g, tau2)
          names(par) <- names(x$state$parameters)
          x$state$parameters <- par
          x$state$fitted.values <- fit
          x$state$edf <- edf
          if(!is.null(x$prior)) {
            if(is.function(x$prior))
              x$state$log.prior <- x$prior(par)
          }
          assign("state", x$state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }

    assign("ic00_val", objfun(tau2 <- get.state(x, "tau2")), envir = env)

    tau2 <- tau2.optim(objfun, start = tau2)
    
    if(!is.null(env$state))
      return(env$state)
    
    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](c(x$state$parameters, x$fixed.hyper)) else x$S[[j]]
    P <- matrix_inv(XWX + S + if(!is.null(x$xt[["pS"]])) x$xt[["pS"]] else 0, index = x$sparse.setup)
    if(is.null(x$xt[["pm"]])) {
      g <- drop(P %*% crossprod(x$X, x$rres))
    } else {
      pS <- if(!is.null(x$xt[["pS"]])) {
        x$xt[["pS"]]
      } else {
        if(!is.null(x$xt[["pSa"]])) {
          1 / tau2[length(tau2)] * x$xt[["pSa"]]
        } else 0
      }
      g <- drop(P %*% (crossprod(x$X, x$rres) + pS %*% x$xt[["pm"]]))
    }

    if(!is.null(x$doCmat)) {
      V <- P %*% t(x$C)
      W <- x$C %*% V
      U <- chol2inv(chol(W)) %*% t(V)
      g <- drop(g - t(U) %*% x$C %*% g)
    }
  }
  
  ## Compute fitted values.
  g <- get.state(x, "b")
  if(any(is.na(g)) | any(g %in% c(-Inf, Inf))) {
    x$state$parameters <- set.par(x$state$parameters, rep(0, length(get.state(x, "b"))), "b")
  }
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))

  if(!is.null(x$doCmat))
    x$state$fitted.values <- x$state$fitted.values - mean(x$state$fitted.values, na.rm = TRUE)

  x$state$edf <- sum_diag(XWX %*% P)

  if(!is.null(x$prior)) {
    if(is.function(x$prior))
      x$state$log.prior <- x$prior(x$state$parameters)
  }
  
  return(x$state)# returing!
}


bfit_iwls_Matrix <- function(x, family, y, eta, id, weights, criterion, ...)
{
  args <- list(...)
  
  peta <- family$map2par(eta)
  if(is.null(args$hess)) {
    ## Compute weights.
    hess <- family$hess[[id]](y, peta, id = id, ...)
  } else hess <- args$hess
  
  if(!is.null(weights))
    hess <- hess * weights
  
  hess <- process.derivs(hess, is.weight = TRUE)
  
  if(is.null(args$z)) {
    ## Score.
    score <- process.derivs(family$score[[id]](y, peta, id = id, ...), is.weight = FALSE)
    
    ## Compute working observations.
    z <- eta[[id]] + 1 / hess * score
  } else z <- args$z
  
  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - fitted(x$state)
  
  ## Compute reduced residuals.
  e <- z - eta[[id]]
  xbin.fun(x$binning$sorted.index, hess, e, x$weights, x$rres, x$binning$order)
  
  ## Compute mean and precision.
  XWX <- crossprod(Diagonal(x = x$weights) %*% x$X, x$X)
  Xr <- crossprod(x$X, x$rres)
  if(!x$state$do.optim | x$fixed | x$fxsp) {
    if(!x$fixed) {
      tau2 <- get.state(x, "tau2")
      S <- Matrix(0, ncol(x$X), ncol(x$X))
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      U <- chol(XWX + S)
    } else {
      U <- chol(XWX)
    }
    P <- chol2inv(U)
    b <- P %*% Xr
    x$state$parameters <- set.par(x$state$parameters, as.numeric(b), "b")
  } else {
    args <- list(...)
    edf0 <- args$edf - x$state$edf
    eta2 <- eta
    
    env <- new.env()
    
    objfun <- function(tau2, ...) {
      S <- Matrix(0, ncol(x$X), ncol(x$X))
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      U <- chol(XWX + S)
      P <- chol2inv(U)
      b <- P %*% Xr
      fit <- x$fit.fun(x$X, b)
      edf <- sum_diag(XWX %*% P)
      eta2[[id]] <- eta2[[id]] + fit
      ic <- get.ic(family, y, family$map2par(eta2), edf0 + edf, length(z), criterion, ...)
      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          par <- c(as.numeric(b), tau2)
          names(par) <- names(x$state$parameters)
          x$state$parameters <- par
          x$state$fitted.values <- fit
          x$state$edf <- edf
          if(!is.null(x$prior)) {
            if(is.function(x$prior))
              x$state$log.prior <- x$prior(par)
          }
          assign("state", x$state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }
    
    assign("ic00_val", objfun(get.state(x, "tau2")), envir = env)
    tau2 <- tau2.optim(objfun, start = get.state(x, "tau2"))
    
    if(!is.null(env$state))
      return(env$state)
    
    S <- Matrix(0, ncol(x$X), ncol(x$X))
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    U <- chol(XWX + S)
    P <- chol2inv(U)
    b <- P %*% Xr
    x$state$parameters <- set.par(x$state$parameters, as.numeric(b), "b")
    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
  }
  
  ## Compute fitted values.
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))
  x$state$edf <- sum_diag(XWX %*% P)
  if(!is.null(x$prior)) {
    if(is.function(x$prior))
      x$state$log.prior <- x$prior(x$state$parameters)
  }
  
  return(x$state)
}


bfit_glmnet <- function(x, family, y, eta, id, weights, criterion, ...)
{
  requireNamespace("glmnet")

  args <- list(...)
  peta <- family$map2par(eta)
  
  hess <- if(is.null(args$hess)) {
    hess <- process.derivs(family$hess[[id]](y, peta, id = id, ...), is.weight = TRUE)
  } else args$hess
  
  if(!is.null(weights))
    hess <- hess * weights
 
  if(is.null(args$z)) {
    score <- process.derivs(family$score[[id]](y, peta, id = id, ...), is.weight = FALSE)
    
    ## Compute working observations.
    z <- eta[[id]] + 1 / hess * score
  } else z <- args$z
  
  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - fitted(x$state)
  
  ## Compute residuals.
  e <- z - eta[[id]]

  if(is.null(x$xt$alpha))
    x$xt$alpha <- 1
  if(is.null(x$xt$nlambda))
    x$xt$nlambda <- 100
  if(is.null(x$xt$lambda.min.ratio))
    x$xt$lambda.min.ratio <- 1e-20

  b <- glmnet::glmnet(x$X, e, alpha = x$xt$alpha,
    nlambda = x$xt$nlambda, standardize = FALSE, intercept = FALSE,
    lambda.min.ratio = x$xt$lambda.min.ratio,
    weights = hess)

  tLL <- b$nulldev - deviance(b)
  k <- b$df
  n <- b$nobs
  IC <- switch(criterion,
    "AICc" = -tLL + 2*k + 2*k*(k + 1)/(n - k - 1),
    "BIC" = -tLL + k * log(n),
    "AIC" = -tLL + 2 * k
  )
  i <- which.min(IC)

  x$state$parameters <- set.par(x$state$parameters, as.numeric(b$lambda[i]), "tau2")
  cb <- as.numeric(coef(b, s = b$lambda[i])[-1])
  x$state$parameters <- set.par(x$state$parameters, cb, "b")
  x$state$fitted.values <- x$X %*% cb
  x$state$edf <- b$df[i]

  return(x$state)
}


## Updating based on optim.
bfit_optim <- function(x, family, y, eta, id, weights, criterion, ...)
{
  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - fitted(x$state)
  eta2 <- eta
  
  tpar <- x$state$parameters
  
  ## Objective for regression coefficients.
  objfun <- function(b, tau2 = NULL) {
    tpar <- set.par(tpar, b, "b")
    if(!is.null(tau2) & !x$fixed)
      tpar <- set.par(tpar, tau2, "tau2")
    fit <- x$fit.fun(x$X, b)
    eta2[[id]] <- eta[[id]] + fit
    ll <- if(is.null(weights[[id]])) {
      family$loglik(y, family$map2par(eta2))
    } else {
      sum(family$d(y, family$map2par(eta2)) * weights[[id]], na.rm = TRUE)
    }
    lp <- x$prior(tpar)
    val <- -1 * (ll + lp)
    if(!is.finite(val)) val <- NA
    val
  }
  
  ## Gradient function.
  grad <- if(!is.null(family$score[[id]]) & is.function(x$grad)) {
    function(gamma, tau2 = NULL) {
      tpar <- set.par(tpar, gamma, "b")
      if(!is.null(tau2) & !x$fixed)
        tpar <- set.par(tpar, tau2, "tau2")
      eta2[[id]] <- eta[[id]] + x$fit.fun(x$X, tpar)
      peta <- family$map2par(eta2)
      score <- drop(family$score[[id]](y, peta))
      grad <- x$grad(score, tpar, full = FALSE)
      return(drop(-1 * grad))
    }
  } else NULL

  suppressWarnings(opt <- try(optim(get.par(tpar, "b"), fn = objfun, gr = grad,
    method = "BFGS", control = list(), tau2 = get.par(tpar, "tau2"), hessian = TRUE,
    lower = if(!is.null(x$force.positive)) 1e-10 else -Inf),
    silent = TRUE))
  
  if(!inherits(opt, "try-error")) {
    tpar <- set.par(tpar, opt$par, "b")
    x$state$fitted.values <- x$fit.fun(x$X, tpar)
    x$state$parameters <- tpar
    x$state$hessian <- opt$hessian
  }
  
  return(x$state)
}


## Compute fitted.values from set of parameters.
get.eta.par <- function(par, x)
{
  nx <- names(x)
  eta <- vector(mode = "list", length = length(nx))
  names(eta) <- nx
  for(j in nx) {
    eta[[j]] <- 0.0
    for(sj in names(x[[j]]$smooth.construct)) {
      xl <- if(sj != "model.matrix") {
        paste(j, "s", x[[j]]$smooth.construct[[sj]]$label, sep = ".")
      } else {
        paste(j, "p", strsplit(x[[j]]$smooth.construct[[sj]]$label, "+", fixed = TRUE)[[1]], sep = ".")
      }
      tpar <- par[grep2(xl, names(par), fixed = TRUE)]
      x[[j]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[j]]$smooth.construct[[sj]]$state$parameters, tpar, "b")
      x[[j]]$smooth.construct[[sj]]$state$fitted.values <- x[[j]]$smooth.construct[[sj]]$fit.fun(x[[j]]$smooth.construct[[sj]]$X,
                                                                                                 get.par(tpar, "b"))
      eta[[j]] <- eta[[j]] + fitted(x[[j]]$smooth.construct[[sj]]$state)
    }
    if(!is.null(x[[j]]$model.matrix)) {
      xl <- paste(j, "p", colnames(x[[j]]$model.matrix), sep = ".")
      tpar <- par[grep(xl, names(par), fixed = TRUE)]
      eta[[j]] <- eta[[j]] + drop(x[[j]]$model.matrix %*% tpar)
    }
  }
  return(eta)
}


## The log-posterior.
log_posterior <- function(par, x, y, family, verbose = TRUE, digits = 3, scale = NULL, ienv = NULL)
{
  nx <- names(x)
  eta <- vector(mode = "list", length = length(nx))
  names(eta) <- nx
  lprior <- 0.0
  for(j in nx) {
    eta[[j]] <- 0.0
    for(sj in names(x[[j]]$smooth.construct)) {
      xl <- if(sj != "model.matrix") {
        paste(j, "s", x[[j]]$smooth.construct[[sj]]$label, sep = ".")
      } else {
        paste(j, "p", strsplit(x[[j]]$smooth.construct[[sj]]$label, "+", fixed = TRUE)[[1]], sep = ".")
      }
      tpar <- par[grep2(xl, names(par), fixed = TRUE)]
      if(x[[j]]$smooth.construct[[sj]]$by == "NA") {
        tpar <- tpar[!grepl(":", names(tpar), fixed = TRUE)]
      }
      bb <- get.par(tpar, "b")
      x[[j]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[j]]$smooth.construct[[sj]]$state$parameters, bb, "b")
      x[[j]]$smooth.construct[[sj]]$state$fitted.values <- x[[j]]$smooth.construct[[sj]]$fit.fun(x[[j]]$smooth.construct[[sj]]$X, bb)
      eta[[j]] <- eta[[j]] + fitted(x[[j]]$smooth.construct[[sj]]$state)
      if(any(grepl("tau2", names(tpar)))) {
        lprior <- lprior + x[[j]]$smooth.construct[[sj]]$prior(c(tpar, x[[j]]$smooth.construct[[sj]]$fixed.hyper))
      } else {
        lprior <- lprior + x[[j]]$smooth.construct[[sj]]$prior(c(tpar, get.state(x[[j]]$smooth.construct[[sj]], "tau2"), x[[j]]$smooth.construct[[sj]]$fixed.hyper))
      }
    }
  }
  ll <- family$loglik(y, family$map2par(eta))
  lp <- as.numeric(ll + lprior)
  
  if(verbose) {
    cat(if(interactive()) "\r" else "\n")
    vtxt <- paste("logLik ", fmt(ll, width = 8, digits = digits),
                  " logPost ", fmt(lp, width = 8, digits = digits),
                  " iteration ", formatC(ienv$bamlss_log_posterior_iteration, width = 4), sep = "")
    cat(vtxt)
    if(.Platform$OS.type != "unix" & interactive()) flush.console()
    bamlss_log_posterior_iteration <- ienv$bamlss_log_posterior_iteration + 1
    assign("bamlss_log_posterior_iteration", bamlss_log_posterior_iteration, envir = ienv)
  }
  
  if(!is.null(scale))
    lp <- lp * scale
  
  return(lp)
}


## Gradient vector of the log-posterior.
grad_posterior <- function(par, x, y, family, ...)
{
  nx <- names(x)
  eta <- vector(mode = "list", length = length(nx))
  names(eta) <- nx
  grad <- NULL
  for(j in nx) {
    eta[[j]] <- 0
    for(sj in names(x[[j]]$smooth.construct)) {
      xl <- if(sj != "model.matrix") {
        paste(j, "s", x[[j]]$smooth.construct[[sj]]$label, sep = ".")
      } else {
        paste(j, "p", strsplit(x[[j]]$smooth.construct[[sj]]$label, "+", fixed = TRUE)[[1]], sep = ".")
      }
      tpar <- par[grep2(xl, names(par), fixed = TRUE)]
      if((x[[j]]$smooth.construct[[sj]]$by == "NA") & (sj != "model.matrix")) {
        tpar <- tpar[!grepl(":", names(tpar), fixed = TRUE)]
      }
      x[[j]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[j]]$smooth.construct[[sj]]$state$parameters, tpar, "b")
      x[[j]]$smooth.construct[[sj]]$state$fitted.values <- x[[j]]$smooth.construct[[sj]]$fit.fun(x[[j]]$smooth.construct[[sj]]$X,
                                                                                                 get.par(tpar, "b"))
      eta[[j]] <- eta[[j]] + fitted(x[[j]]$smooth.construct[[sj]]$state)
    }
  }
  for(j in nx) {
    score <- family$score[[j]](y, family$map2par(eta), id = j)
    for(sj in names(x[[j]]$smooth.construct)) {
      tgrad <- x[[j]]$smooth.construct[[sj]]$grad(score, c(x[[j]]$smooth.construct[[sj]]$state$parameters, x[[j]]$smooth.construct[[sj]]$fixed.hyper), full = FALSE)
      grad <- c(grad, tgrad)
    }
  }
  return(grad)
}


## Optimizer based on optim().
opt_optim <- function(x, y, family, start = NULL, verbose = TRUE, digits = 3,
                gradient = TRUE, hessian = FALSE, eps = .Machine$double.eps^0.5, maxit = 100, ...)
{
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("design construct names mismatch with family names!")
  
  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, ...)
  
  if(!is.null(start))
    x <- set.starting.values(x, start)
  
  nobs <- nrow(y)
  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }
  
  for(i in names(x)) {
    for(j in seq_along(x[[i]]$smooth.construct)) {
      if(is.null(x[[i]]$smooth.construct[[j]]$grad)) {
        gradient <- FALSE
      } else {
        if(!all(c("score", "parameters", "full") %in% names(formals(x[[i]]$smooth.construct[[j]]$grad))))
          gradient <- FALSE
      }
    }
  }
  
  par <- get.all.par(x, list = FALSE, drop = TRUE)
  
  ienv <- NULL
  if(verbose) {
    ienv <- new.env()
    bamlss_log_posterior_iteration <- 1
    assign("bamlss_log_posterior_iteration", bamlss_log_posterior_iteration, envir = ienv)
  }
  
  if(!hessian) {
    opt <- optim(par, fn = log_posterior,
                 gr = if(!is.null(family$score) & gradient) grad_posterior else NULL,
                 x = x, y = y, family = family, method = "BFGS", verbose = verbose,
                 digits = digits, ienv = ienv, control = list(fnscale = -1, reltol = eps, maxit = maxit),
                 hessian = TRUE)
    
    if(verbose) {
      cat("\n")
      rm(ienv)
    }
    
    eta <- get.eta.par(opt$par, x)
    
    return(list("parameters" = opt$par, "fitted.values" = eta,
                "logPost" = opt$value, "logLik" = family$loglik(y, family$map2par(eta)),
                "nobs" = nobs, "hessian" = opt$hessian, "converged" = opt$convergence < 1))
  } else {
    fn <- if(is.null(family$p2d)) {
      log_posterior
    } else function(par, ...) { sum(family$p2d(par, log = TRUE), na.rm = TRUE) }
    opt <- optimHess(par, fn = fn,
                     gr = if(!is.null(family$score) & gradient & is.null(family$p2d)) grad_posterior else NULL,
                     x = x, y = y, family = family, verbose = verbose, digits = digits, ienv = ienv,
                     control = list(fnscale = -1, reltol = eps, maxit = maxit))
    
    if(verbose) {
      cat("\n")
      rm(ienv)
    }
    
    return(opt)
  }
}


## Fast computation of weights and residuals when binning.
xbin.fun <- function(ind, weights, e, xweights, xrres, oind, uind = NULL)
{
  if(inherits(ind, "ff")) {
    stop("ff support stops here!")
  } else {
    .Call("xbin_fun", as.integer(ind), as.numeric(weights), 
          as.numeric(e), as.numeric(xweights), as.numeric(xrres),
          as.integer(oind), PACKAGE = "bamlss")
  }
  invisible(NULL)
}


xcenter <- function(x)
{
  .Call("xcenter", as.numeric(x), PACKAGE = "bamlss")
}


## Modified likelihood based boosting.
opt_boostm <- boostm <- function(x, y, family, offset = NULL,
  nu = 0.1, df = 3, maxit = 400, mstop = NULL,
  verbose = TRUE, digits = 4, flush = TRUE,
  eps = .Machine$double.eps^0.25, plot = TRUE,
  initialize = TRUE, stop.criterion = "BIC",
  force.stop = !is.null(stop.criterion),
  do.optim = TRUE, always = FALSE, ...)
{
  ## FIXME: hard coded!
  offset <- weights <- NULL

  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")
  
  if(!is.null(mstop))
    maxit <- mstop

  if(!is.null(stop.criterion))
    stop.criterion <- toupper(stop.criterion)
  
  if(is.null(maxit))
    stop("please set either argument 'maxit' or 'mstop'!")

  always2 <- FALSE
  if(!is.logical(always)) {
    if(is.character(always)) {
      if(!is.na(pmatch(always, "best"))) {
        always2 <- TRUE
        always <- TRUE
      } else {
        always <- FALSE
      }
    }
  }
  
  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, df = NULL, ...)
  
  np <- length(nx)
  nobs <- nrow(y)
  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }
  
  ## Setup boosting structure, i.e, all parametric
  ## terms get an entry in $smooth.construct object.
  ## Intercepts are initalized.
  x <- boost_transform(x = x, y = y, df = df, family = family,
    maxit = maxit, eps = eps, initialize = initialize, offset = offset,
    weights = weights)
  for(i in nx) {
    for(j in names(x[[i]]$smooth.construct))
      x[[i]]$smooth.construct[[j]]$criterion <- x[[i]]$smooth.construct[[j]]$loglik
  }
  
  ## Create a list() that saves the states for
  ## all parameters and model terms.
  states <- make.state.list(x)
  
  ## Matrix of all parameters.
  parm <- make.par.list(x, iter = maxit)
  
  ## Term selector help vectors.
  select <- rep(NA, length = length(nx))
  names(select) <- nx
  loglik <- select

  ## Criterion used.
  sic <- if(is.null(stop.criterion)) "BIC" else stop.criterion
  
  ## Save criterion in list().
  crit <- ll.contrib <- make.state.list(x, type = 2)

  ## Extract actual predictor.
  eta <- get.eta(x)

  if(!is.null(offset)) {
    offset <- as.data.frame(offset)
    for(j in nx) {
      if(!is.null(offset[[j]]))
        eta[[j]] <- eta[[j]] + offset[[j]]
    }
  }
  
  ## Print stuff.
  ia <- if(flush) interactive() else FALSE
 
  ## Save edf and IC?
  edf <- save.ic <- rep(NA, maxit)
  
  ## Env for C.
  rho <- new.env()
  
  ## Start boosting.
  eps0 <- 1; iter <- if(initialize) 2 else 1
  save.ll <- NULL; stopped <- FALSE
  ll <- family$loglik(y, family$map2par(eta))
  medf <- get.edf(x, type = 2) + if(initialize) length(nx) else 0
  ic0 <- -2 * ll + medf * (if(tolower(sic) == "aic") 2 else log(nobs))
  ptm <- proc.time()
  while(iter <= maxit) {
    eta0 <- eta

    ## Actual parameters.
    peta <- family$map2par(eta)
    
    ## Cycle through all parameters and terms.
    for(i in nx) {
      ## Actual gradient.
      grad <- process.derivs(family$score[[i]](y, peta, id = i), is.weight = FALSE)

      ## Actual hessian.
      hess <- process.derivs(family$score[[i]](y, peta, id = i), is.weight = FALSE)

      ## Working response.
      z <- eta[[i]] + 1 / hess * grad
 
      for(j in names(x[[i]]$smooth.construct)) {
        if(always) {
          if(j == "(Intercept)") {
            crit[[i]][j] <- -Inf
            next
          }
        }

        ## Get update.
        states[[i]][[j]] <- if(is.null(x[[i]]$smooth.construct[[j]][["boostm.fit"]])) {
          boostm_fit(x[[i]]$smooth.construct[[j]], grad, hess, z, nu, stop.criterion,
            family, y, eta, medf, id = i, do.optim = do.optim, ...)
        } else {
          x[[i]]$smooth.construct[[j]][["boostm.fit"]](x[[i]]$smooth.construct[[j]],
            grad = grad, hess = hess, z = z, nu = nu, criterion = stop.criterion,
            family = family, y = y, eta = eta, edf = medf, id = i,
            do.optim = do.optim, iteration = iter, ...)
        }
        
        ## Get contribution.
        eta[[i]] <- eta[[i]] + fitted(states[[i]][[j]]) - fitted(x[[i]]$smooth.construct[[j]]$state)
        tll <- family$loglik(y, family$map2par(eta))
        if(is.null(stop.criterion)) {
          crit[[i]][j] <- -1 * (ll - tll)
        } else {
          tedf <- medf - x[[i]]$smooth.construct[[j]]$state$edf + states[[i]][[j]]$edf
          ic1 <- -2 * tll + tedf * (if(tolower(stop.criterion) == "aic") 2 else log(nobs))
          crit[[i]][j] <- ic1
        }
        ll.contrib[[i]][j] <- tll - ll
        eta[[i]] <- eta0[[i]]
      }
      
      ## Which one is best?
      select[i] <- which.min(crit[[i]])
    }

    i <- which.min(sapply(crit, function(x) { min(x) }))
    
    ## Which term to update.
    take <- c(nx[i], names(crit[[i]])[select[i]])

    ## Update selected term.
    eta[[take[1]]] <- eta[[take[1]]] + fitted(states[[take[1]]][[take[2]]]) - fitted(x[[take[1]]]$smooth.construct[[take[2]]]$state)

    ## Save parameters.
    parm[[take[1]]][[take[2]]][iter, ] <- get.par(states[[take[1]]][[take[2]]]$parameters, "b") - get.par(x[[take[1]]]$smooth.construct[[take[2]]]$state$parameters, "b")
    medf <- medf - x[[take[1]]]$smooth.construct[[take[2]]]$state$edf + states[[take[1]]][[take[2]]]$edf

    ## Write to x.
    x[[take[1]]]$smooth.construct[[take[2]]]$state <- states[[take[1]]][[take[2]]]
    x[[take[1]]]$smooth.construct[[take[2]]]$selected[iter] <- 1
    x[[take[1]]]$smooth.construct[[take[2]]]$loglik[iter] <- ll.contrib[[take[1]]][take[2]]
    x[[take[1]]]$smooth.construct[[take[2]]]$criterion[iter] <- -1 * crit[[take[1]]][take[2]]

    ## Intercept updating.
    if(always) {
      nxa <- if(always2) take[1] else nx

      for(ii in nxa) {
        if("(Intercept)" %in% names(x[[ii]]$smooth.construct)) {
          ## Actual gradient.
          grad <- process.derivs(family$score[[ii]](y, peta, id = ii), is.weight = FALSE)

          ## Actual hessian.
          hess <- process.derivs(family$score[[ii]](y, peta, id = ii), is.weight = FALSE)

          ## Working response.
          z <- eta[[ii]] + 1 / hess * grad

          ## Get update.
          states[[ii]][["(Intercept)"]] <- boostm_fit(x[[ii]]$smooth.construct[["(Intercept)"]],
            grad, hess, z, nu, stop.criterion, family, y, eta, medf, id = i, do.optim = do.optim, ...)

          ll <- family$loglik(y, family$map2par(eta))

          ## Update predictor.
          eta[[ii]] <- eta[[ii]] + fitted(states[[ii]][["(Intercept)"]]) - fitted(x[[ii]]$smooth.construct[["(Intercept)"]]$state)

          ## Save parameters.
          parm[[ii]][["(Intercept)"]][iter, ] <- get.par(states[[ii]][["(Intercept)"]]$parameters, "b") - get.par(x[[ii]]$smooth.construct[["(Intercept)"]]$state$parameters, "b")
          medf <- medf - x[[ii]]$smooth.construct[["(Intercept)"]]$state$edf + states[[ii]][["(Intercept)"]]$edf

          tll <- family$loglik(y, family$map2par(eta))
          ll.contrib[[ii]]["(Intercept)"] <- tll - ll

          ## Write to x.
          x[[ii]]$smooth.construct[["(Intercept)"]]$state <- states[[ii]][["(Intercept)"]]
          x[[ii]]$smooth.construct[["(Intercept)"]]$selected[iter] <- 1
          x[[ii]]$smooth.construct[["(Intercept)"]]$loglik[iter] <- ll.contrib[[ii]]["(Intercept)"]
        }
      }
    }

    edf[iter] <- medf
    
    ## Change.
    eps0 <- do.call("cbind", eta)
    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1
    
    ll <- family$loglik(y, family$map2par(eta))
    ic0 <- -2 * ll + medf * (if(tolower(sic) == "aic") 2 else log(nobs))

    save.ll <- c(save.ll, ll)
    save.ic[iter] <- ic0

    qsel <- get.qsel(x, iter)
    
    if(verbose) {
      cat(if(ia) "\r" else "\n")
      vtxt <- paste(
        paste(sic, " ", fmt(save.ic[iter], width = 8, digits = digits), " ", sep = ""),
        "logLik ", fmt(ll, width = 8, digits = digits),
        " edf ", fmt(edf[iter], width = 4, digits = digits), " ",
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)),
        " qsel ", qsel, sep = "")
      cat(vtxt)
      
      if(.Platform$OS.type != "unix" & ia) flush.console()
    }

    if(!is.null(stop.criterion)) {
      if(iter > 2) {
        if(!is.na(save.ic[iter - 1]) & force.stop) {
          if(save.ic[iter - 1] < save.ic[iter]) {
            stopped <- TRUE
            break
          }
        }
      }
    }
    
    iter <- iter + 1
  }

  elapsed <- c(proc.time() - ptm)[3]
  
  if(verbose) {
    cat("\n")
    et <- if(elapsed > 60) {
      paste(formatC(format(round(elapsed / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(elapsed, 2), nsmall = 2), width = 5), "sec", sep = "")
    cat("\n elapsed time: ", et, "\n", sep = "")
  }
  
  itr <- if(!stopped) maxit else (iter - 1)
  bsum <- make.boost_summary(x, itr, save.ll, edf, FALSE, nobs)
  bsum$criterion <- list(
    "bic" = -2 * save.ll[1:itr] + edf[1:itr] * log(nobs),
    "aic" = -2 * save.ll[1:itr] + edf[1:itr] * 2,
    "edf" = edf[1:itr]
  )
  if(plot)
    plot.boost_summary(bsum)
  
  return(list("parameters" = parm2mat(parm, itr),
    "fitted.values" = eta, "nobs" = nobs, "boost_summary" = bsum,
    "runtime" = elapsed))
}


## Gradient boosting.
opt_boost <- boost <- function(x, y, family, weights = NULL, offset = NULL,
  nu = 0.1, nu.adapt = TRUE, df = 4, maxit = 400, mstop = NULL,
  maxq = NULL, qsel.splitfactor = FALSE,
  verbose = TRUE, digits = 4, flush = TRUE,
  eps = .Machine$double.eps^0.25, nback = NULL, plot = TRUE,
  initialize = TRUE, stop.criterion = NULL, select.type = 1, force.stop = TRUE,
  hatmatrix = !is.null(stop.criterion), reverse.edf = FALSE, approx.edf = FALSE,
  always = FALSE, ...)
{
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")

  if(reverse.edf | approx.edf)
    hatmatrix <- FALSE
  
  if(!is.null(mstop))
    maxit <- mstop

  light <- list(...)$boost.light
  if(is.null(light))
    light <- FALSE
  
  if(!is.null(nback)) {
    if(is.null(maxit))
      maxit <- 10000
  }
  nu <- rep(nu, length.out = length(nx))
  names(nu) <- nx

  always2 <- always3 <- FALSE
  if(!is.logical(always)) {
    if(is.character(always)) {
      if(!is.na(pmatch(always, "best"))) {
        always2 <- TRUE
        always <- TRUE
      } else {
        if(!is.na(pmatch(always, "yes"))) {
          always3 <- TRUE
          always <- TRUE
        } else {
          always <- FALSE
        }
      }
    }
  }
  
  if(is.null(maxit))
    stop("please set either argument 'maxit' or 'mstop'!")
  
  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, df = NULL, nodf = TRUE, ...)
  
  np <- length(nx)
  nobs <- nrow(y)

  CRPS <- !is.null(list(...)$crps) | !is.null(list(...)$CRPS)

  yname <- names(y)

  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }

  if(!is.null(offset))
    initialize <- FALSE
  
  ## Setup boosting structure, i.e, all parametric
  ## terms get an entry in $smooth.construct object.
  ## Intercepts are initalized.
  x <- boost_transform(x = x, y = y, df = df, family = family,
    maxit = maxit, eps = eps, initialize = initialize, offset = offset,
    weights = weights, always3 = always3, ...)

  if(!is.null(list(...)$ret.x)) {
    if(list(...)$ret.x)
      return(x)
  }

  ## Create a list() that saves the states for
  ## all parameters and model terms.
  states <- make.state.list(x)
  
  ## Matrix of all parameters.
  parm <- make.par.list(x, iter = if(light) 1L else maxit)

  ## Term selector help vectors.
  select <- rep(NA, length = length(nx))
  names(select) <- nx
  loglik <- select
  
  ## Save rss in list().
  rss <- make.state.list(x, type = 2)
  
  ## Extract actual predictor.
  eta <- get.eta(x)

  if(!is.null(offset)) {
    offset <- as.data.frame(offset)
    for(j in nx) {
      if(!is.null(offset[[j]]))
        eta[[j]] <- eta[[j]] + offset[[j]]
    }
  }

  W <- NULL
  if(!is.null(weights)) {
    if(attr(weights, "identical"))
      W <- as.numeric(weights[, 1])
  }
  
  ## Print stuff.
  ia <- if(flush) interactive() else FALSE
 
  ## Hat matrix?
  HatMat <- list()
  edf <- Imat <- save.ic <- NULL
  if(hatmatrix) {
    for(i in nx)
      HatMat[[i]] <- diag(length(eta[[1]]))
    edf <- rep(0, maxit)
    if(!is.null(stop.criterion))
      save.ic <- rep(NA, maxit)
    Imat <- diag(nobs)
  }
  if(reverse.edf | approx.edf) {
    edf <- rep(0, maxit)
    if(!is.null(stop.criterion))
      save.ic <- rep(NA, maxit)
  }
  selectfun <- list(...)$selectfun
  selectmodel <- list(...)$selectmodel
  nthreads <- list(...)$nthreads
  if(is.null(selectmodel))
    selectmodel <- TRUE
  if(!is.null(selectfun))
    save.ic <- rep(NA, maxit)

  if(!is.null(selectfun))
    stop.criterion <- "userIC"
  
  ## Env for C.
  rho <- new.env()

  if(is.null(maxq))
    maxq <- Inf
  qsel <- 0
  
  ## Start boosting.
  eps0 <- 1; iter <- if(initialize) 2 else 1
  save.ll <- NULL
  ll <- if(is.null(W)) {
    family$loglik(y, family$map2par(eta))
  } else {
    sum(family$d(y, family$map2par(eta)) * W)
  }
  redf <- if(initialize) length(nx) else 0
  loglik <- loglik2 <- NULL
  iter_ll2 <- 0
  nu0 <- nu
  ptm <- proc.time()
  while(iter <= maxit & qsel < maxq) {
    if(iter > 2)
      loglik2 <- loglik

    eta0 <- eta
    
    ## Cycle through all parameters
    for(i in nx) {
      peta <- family$map2par(eta)
      
      ## Actual gradient.
      grad <- process.derivs(family$score[[i]](y, peta, id = i), is.weight = FALSE)

      if(length(grad) != nobs)
        stop("something wrong in processing the family $score() function! More elements in return value of $score() than the response!")
      
      ## Fit to gradient.
      for(j in names(x[[i]]$smooth.construct)) {
        if(always) {
          if(j == "(Intercept)") {
            rss[[i]][j] <- Inf
            next
          }
        }

        ## Get updated parameters.
        nu2 <- if(inherits(x[[i]]$smooth.construct[[j]], "nnet.boost")) nu[i] else nu[i]
        states[[i]][[j]] <- if(is.null(x[[i]]$smooth.construct[[j]][["boost.fit"]])) {
          if(hatmatrix) {
            boost_fit(x[[i]]$smooth.construct[[j]], grad, nu2,
              hatmatrix = hatmatrix, weights = if(!is.null(weights)) weights[, i] else NULL,
              nthreads = nthreads)
          } else {
            try(.Call("boost_fit", x[[i]]$smooth.construct[[j]], grad, nu2,
              if(!is.null(weights)) as.numeric(weights[, i]) else numeric(0), rho, PACKAGE = "bamlss"), silent = TRUE)
          }
        } else {
          x[[i]]$smooth.construct[[j]][["boost.fit"]](x = x[[i]]$smooth.construct[[j]],
            y = grad, nu = nu2, hatmatrix = hatmatrix,
            weights = if(!is.null(weights)) weights[, i] else NULL,
            rho = rho, nthreads = nthreads, always3 = always3)
        }
        
        ## Get rss.
        if(is.null(selectfun)) {
          if(is.null(stop.criterion)) {
            rss[[i]][j] <- states[[i]][[j]]$rss
          } else {
            if(select.type == 1) {
              rss[[i]][j] <- states[[i]][[j]]$rss
            } else {
              teta <- eta
              teta[[i]] <- teta[[i]] + fitted(states[[i]][[j]])
              if(is.null(W))
                tll <- family$loglik(y, family$map2par(teta))
              else
                tll <- sum(family$d(y, family$map2par(teta), log = TRUE) * W)
              if(!light) {
                if(approx.edf) {
                  tredf <- redf
                  if(!is.null(x[[i]]$smooth.construct[[j]]$is.model.matrix) | inherits(x[[i]]$smooth.construct[[j]], "nnet.boost")) {
                    if(x[[i]]$smooth.construct[[j]]$state$init.edf < 1)
                      tredf <- tredf + 1
                  } else {
                    if(inherits(x[[i]]$smooth.construct[[j]], "lasso.smooth") | inherits(x[[i]]$smooth.construct[[j]], "nnet.smooth")) {
                      if(iter < 2) {
                        aset <- if(x[[i]]$smooth.construct[[j]]$fuse) {
                          sum(abs(unique.fuse(get.par(states[[i]][[j]]$parameters, "b"))) > 1e-10)
                        } else {
                          sum(abs(get.par(states[[i]][[j]]$parameters, "b")) > 1e-10)
                        }
                        tredf <- tredf + aset
                      } else {
                        aset0 <- apply(parm[[i]][[j]][1:(iter - 1L), , drop = FALSE], 2, sum)
                        aset1 <- apply(rbind(parm[[i]][[j]][1:(iter - 1L), , drop = FALSE],
                          get.par(states[[i]][[j]]$parameters, "b")), 2, sum)
                        if(x[[i]]$smooth.construct[[j]]$fuse) {
                          aset0 <- sum(abs(unique.fuse(aset0)) > 1e-10)
                          aset1 <- sum(abs(unique.fuse(aset1)) > 1e-10)
                        } else {
                          aset0 <- sum(abs(aset0) > 1e-10)
                          aset1 <- sum(abs(aset1) > 1e-10)
                        }
                        aset <- aset1 - aset0
                        tredf <- tredf + aset
                      }
                    } else {
                      tredf <- tredf + nu[i] * x[[i]]$smooth.construct[[j]]$state$init.edf
                    }
                  }
                  rss[[i]][j] <- -2 * tll + tredf * (if(tolower(stop.criterion) == "aic") 2 else log(nobs))
                } else {
                  if(reverse.edf) {
                    states[[i]][[j]]$redf <- reverse_edf(x = x[[i]]$smooth.construct[[j]], bn = get.par(states[[i]][[j]]$parameters, "b"),
                      bmat = parm[[i]][[j]][1:iter, , drop = FALSE], nobs, grad, teta[[i]])
                    tredf <- redf + states[[i]][[j]]$redf$edf
                    rss[[i]][j] <- -2 * tll + tredf * (if(tolower(stop.criterion) == "aic") 2 else log(nobs))
                  } else {
                    ## tedf0 <- sum(diag(Imat - HatMat[[i]] %*% (Imat - states[[i]][[j]]$hat)))
                    tedf <- hatmat_trace(HatMat[[i]], states[[i]][[j]]$hat)
                    if(length(nxr <- nx[nx != i])) {
                      for(ii in nxr)
                        tedf <- tedf + hatmat_sumdiag(HatMat[[i]])
                    }
                    rss[[i]][j] <- -2 * tll + tedf * (if(tolower(stop.criterion) == "aic") 2 else log(nobs))
                  }
                }
              }
            }
          }
        } else {
          rss[[i]][j] <- selfun(iter = iter, i = i, j = j, state = states[[i]][[j]],
            parm = parm, x = x, family = family, sfun = selectfun, yname = yname, weights = weights,
            selectmodel = selectmodel)
        }
      }
      
      ## Which one is best?
      if(always & (length(rss[[i]]) < 2)) {
        if(names(rss[[i]]) == "(Intercept)")
          next
      }
      select[i] <- which.min(rss[[i]])

      if(nu.adapt) {
        tbeta <- get.par(states[[i]][[select[i]]]$parameters, "b") * 1 / nu[i]
        fv <- function(v) {
          beta <- v * tbeta
          eta[[i]] <- eta[[i]] + x[[i]]$smooth.construct[[select[i]]]$fit.fun(x[[i]]$smooth.construct[[select[i]]]$X, beta)
          family$loglik(y, family$map2par(eta))
        }
        v <- optimize(fv, interval = c(.Machine$double.eps^0.5, 1), maximum = TRUE)$maximum
        beta <- nu[i] * v * tbeta
        states[[i]][[select[i]]]$parameters <- set.par(states[[i]][[select[i]]]$parameters, beta, "b")
        states[[i]][[select[i]]]$fitted.values <- x[[i]]$smooth.construct[[select[i]]]$fit.fun(x[[i]]$smooth.construct[[select[i]]]$X, beta)
      }
      
      ## Compute likelihood contribution.
      eta[[i]] <- eta[[i]] + fitted(states[[i]][[select[i]]])
      llf <- if(is.null(W)) {
        if(CRPS & !is.null(family$crps)) {
          -1 * family$crps(y, family$map2par(eta))
        } else {
          family$loglik(y, family$map2par(eta))
        }
      } else {
        sum(family$d(y, family$map2par(eta), log = TRUE) * W)
      }
      loglik[i] <- -1 * (ll - llf)
      
      eta[[i]] <- eta0[[i]]
    }

    if(is.null(stop.criterion) & is.null(selectfun)) {
      i <- which.max(loglik)
    } else {
      i <- if(select.type == 1) {
        which.max(loglik)
      } else {
        which.min(sapply(rss, function(x) { min(x) }))
      }
    }
    
    ## Which term to update.
    take <- c(nx[i], names(rss[[i]])[select[i]])
    
    ## Update selected base learner.
    eta[[take[1]]] <- eta[[take[1]]] + states[[take[1]]][[take[2]]]$fitted.values
    
    ## Write to x.
    if(is.null(x[[take[1]]]$smooth.construct[[take[2]]][["increase"]])) {
      x[[take[1]]]$smooth.construct[[take[2]]]$state <- increase(x[[take[1]]]$smooth.construct[[take[2]]]$state,
        states[[take[1]]][[take[2]]])
    } else {
      x[[take[1]]]$smooth.construct[[take[2]]]$state <- x[[take[1]]]$smooth.construct[[take[2]]][["increase"]](x[[take[1]]]$smooth.construct[[take[2]]]$state,
        states[[take[1]]][[take[2]]])
    }
    x[[take[1]]]$smooth.construct[[take[2]]]$selected[iter] <- 1
    x[[take[1]]]$smooth.construct[[take[2]]]$loglik[iter] <- loglik[i]
    
    ## Save parameters.
    if(always3) {
      tpar <- get.par(states[[take[1]]][[take[2]]]$parameters, "b")
      x[[take[1]]]$smooth.construct[["(Intercept)"]]$selected[iter] <- 1
      ##parm[[take[1]]][["(Intercept)"]][iter, ] <- tpar[1]
      if(light) {
        parm[[take[1]]][[take[2]]] <- parm[[take[1]]][[take[2]]] + tpar[-1]
      } else {
        parm[[take[1]]][[take[2]]][iter, ] <- tpar[-1]
      }
    } else {
      if(light) {
        parm[[take[1]]][[take[2]]] <- parm[[take[1]]][[take[2]]] + get.par(states[[take[1]]][[take[2]]]$parameters, "b")
      } else {
        parm[[take[1]]][[take[2]]][iter, ] <- get.par(states[[take[1]]][[take[2]]]$parameters, "b")
      }
    }

    ## Intercept updating.
    if(always) {
      ll <- if(is.null(W)) {
        family$loglik(y, family$map2par(eta))
      } else {
        sum(family$d(y, family$map2par(eta), log = TRUE) * W)
      }
      nxa <- if(always2) take[1] else nx
      for(ii in nxa) {
        if("(Intercept)" %in% names(x[[ii]]$smooth.construct)) {
          if(always3) {
            if(ii == take[1])
              next
          }
          peta <- family$map2par(eta)
      
          ## Actual gradient.
          grad <- process.derivs(family$score[[ii]](y, peta, id = ii), is.weight = FALSE)

          ## Update.
          states[[ii]][["(Intercept)"]] <- if(hatmatrix) {
            boost_fit(x[[ii]]$smooth.construct[["(Intercept)"]], grad, nu[ii],
              hatmatrix = hatmatrix, weights = if(!is.null(weights)) weights[, ii] else NULL)
          } else {
            .Call("boost_fit", x[[ii]]$smooth.construct[["(Intercept)"]], grad, nu[ii],
              if(!is.null(weights)) as.numeric(weights[, ii]) else numeric(0), rho, PACKAGE = "bamlss")
          }
          eta[[ii]] <- eta[[ii]] + fitted(states[[ii]][["(Intercept)"]])
          x[[ii]]$smooth.construct[["(Intercept)"]]$state <- increase(x[[ii]]$smooth.construct[["(Intercept)"]]$state,
            states[[ii]][["(Intercept)"]])
          x[[ii]]$smooth.construct[["(Intercept)"]]$selected[iter] <- 1
          x[[ii]]$smooth.construct[["(Intercept)"]]$loglik[iter] <- -1 * (ll - family$loglik(y, family$map2par(eta)))
          if(light) {
            parm[[ii]][["(Intercept)"]] <- parm[[ii]][["(Intercept)"]] + get.par(states[[ii]][["(Intercept)"]]$parameters, "b")
          } else {
            parm[[ii]][["(Intercept)"]][iter, ] <- get.par(states[[ii]][["(Intercept)"]]$parameters, "b")
          }
          if(approx.edf) {
            if(x[[ii]]$smooth.construct[["(Intercept)"]]$state$init.edf < 1) {
              redf <- redf + 1
              x[[ii]]$smooth.construct[["(Intercept)"]]$state$init.edf <- 1
            }
          }
        }
      }
    }

    eps0 <- do.call("cbind", eta)
    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1
    
    peta <- family$map2par(eta)
    ll <- if(is.null(W)) {
      if(CRPS & !is.null(family$crps)) {
        -1 * family$crps(y, peta)
      } else {
        family$loglik(y, peta)
      }
    } else {
      sum(family$d(y, peta, log = TRUE) * W)
    }
    save.ll <- c(save.ll, ll)

    if(hatmatrix) {
      HatMat[[take[1]]] <- HatMat[[take[1]]] %*% (Imat - x[[take[1]]]$smooth.construct[[take[2]]]$state$hat)
      for(i in nx)
        edf[iter] <- edf[iter] + hatmat_sumdiag(HatMat[[i]])
      if(!is.null(stop.criterion)) {
        save.ic[iter] <- -2 * ll + edf[iter] * (if(tolower(stop.criterion) == "aic") 2 else log(nobs))
        if(iter > (if(initialize) 2 else 1)) {
          if(!is.na(save.ic[iter - 1]) & force.stop) {
            if(save.ic[iter - 1] < save.ic[iter]) {
              nback <- TRUE
              break
            }
          }
        }
      }
    }

    if(reverse.edf | approx.edf) {
      if(approx.edf) {
        if(!is.null(x[[take[1]]]$smooth.construct[[take[2]]]$is.model.matrix) | inherits(x[[take[1]]]$smooth.construct[[take[2]]], "nnet.boost")) {
          if(x[[take[1]]]$smooth.construct[[take[2]]]$state$init.edf < 1) {
            redf <- redf + 1
            x[[take[1]]]$smooth.construct[[take[2]]]$state$init.edf <- 1
          }
        } else {
          if(inherits(x[[take[1]]]$smooth.construct[[take[2]]], "lasso.smooth") | inherits(x[[take[1]]]$smooth.construct[[take[2]]], "nnet.smooth")) {
            if(iter < 2) {
              aset <- if(x[[take[1]]]$smooth.construct[[take[2]]]$fuse) {
                  sum(abs(unique.fuse(parm[[take[1]]][[take[2]]][if(light) 1L else iter, ])) > 1e-10)
                } else {
                  sum(abs(parm[[take[1]]][[take[2]]][if(light) 1L else iter, ]) > 1e-10)
                }
              redf <- redf + aset
            } else {
              aset0 <- apply(parm[[take[1]]][[take[2]]][if(light) 1L else 1:(iter - 1L), , drop = FALSE], 2, sum)
              aset1 <- apply(parm[[take[1]]][[take[2]]][if(light) 1L else 1:iter, , drop = FALSE], 2, sum)
              if(x[[take[1]]]$smooth.construct[[take[2]]]$fuse) {
                aset0 <- sum(abs(unique.fuse(aset0)) > 1e-10)
                aset1 <- sum(abs(unique.fuse(aset1)) > 1e-10)
              } else {
                aset0 <- sum(abs(aset0) > 1e-10)
                aset1 <- sum(abs(aset1) > 1e-10)
              }
              aset <- aset1 - aset0
              redf <- redf + aset
            }
          } else {
            redf <- redf + nu[take[1]] * x[[take[1]]]$smooth.construct[[take[2]]]$state$init.edf
          }
        }
      } else {
        if(is.null(stop.criterion))
          stop("reverse.edf not implemented!")
        redf <- redf + states[[take[1]]][[take[2]]]$redf$edf
      }
      edf[iter] <- redf
      if(!is.null(stop.criterion)) {
        save.ic[iter] <- -2 * ll + edf[iter] * (if(tolower(stop.criterion) == "aic") 2 else log(nobs))
        if(iter > ((if(initialize) 2 else 1) * 100)) {
          if(!is.na(save.ic[iter - 1]) & force.stop) {
            if(save.ic[iter - 1] < save.ic[iter]) {
              nback <- TRUE
              break
            }
          }
        }
      }
    }

    if(!is.null(selectfun)) {
      save.ic[iter] <- min(unlist(rss))
      if(force.stop & (iter > (if(initialize) 2 else 1))) {
        if(save.ic[iter - 1] < save.ic[iter]) {
          nback <- TRUE
          break
        }
      }
    }

    ## Compute number of selected base learners.
    qsel <- get.qsel(x, if(light) 1L else iter, qsel.splitfactor = qsel.splitfactor)
    
    if(verbose) {
      cat(if(ia) "\r" else "\n")
      vtxt <- paste(
        if(!is.null(stop.criterion)) paste(stop.criterion, " ", fmt(save.ic[iter], width = 8, digits = digits), " ", sep = "") else NULL,
        if(CRPS & !is.null(family$crps)) "CRPS" else "logLik ", fmt(ll, width = 8, digits = digits),
        if(hatmatrix | reverse.edf | approx.edf) paste(" edf ", fmt(edf[iter], width = 4, digits = digits), " ", sep = "") else NULL,
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)),
        " qsel ", qsel, sep = "")
      cat(vtxt)
      
      if(.Platform$OS.type != "unix" & ia) flush.console()
    }

    if((iter > 2) & all(loglik2 == loglik)) {
      warning("no more improvements in the log-likelihood, setting nu = nu * 0.9!")
      ## nu[take[1]] <- nu[take[1]] * 0.9
      iter_ll2 <- iter_ll2 + 1
    }

    if(all(nu < .Machine$double.eps^0.5) & (iter_ll2 > 10)) {
      nback <- TRUE
      warning(paste("no more improvements after", iter_ll2, "iterations in the log-likelihood, stopped!"))
      break
    }
    
    iter <- iter + 1
    
    if(!is.null(nback)) {
      if(iter > nback) {
        dll <- abs(diff(tail(save.ll, nback)))
        if(any(!is.finite(dll)) | any(is.na(dll)))
          break
        if(all(dll < eps))
          break
      }
    }
  }
  elapsed <- c(proc.time() - ptm)[3]
  
  if(verbose) {
    cat("\n")
    et <- if(elapsed > 60) {
      paste(formatC(format(round(elapsed / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(elapsed, 2), nsmall = 2), width = 5), "sec", sep = "")
    cat("\n elapsed time: ", et, "\n", sep = "")
  }
  
  bsum <- make.boost_summary(x, if(is.null(nback)) maxit else (iter - 1), save.ll, edf,
    (hatmatrix | approx.edf | reverse.edf), length(eta[[1]]))
  if(plot)
    plot.boost_summary(bsum)

  if(!is.null(selectfun)) {
    if(is.null(bsum$criterion))
      bsum$criterion <- list()
    bsum$criterion$userIC <- save.ic[1:(if(is.null(nback)) maxit else (iter - 1))]
  }
  
  return(list("parameters" = parm2mat(parm, if(light) { 1L} else { if(is.null(nback)) maxit else (iter - 1) }),
    "fitted.values" = eta, "nobs" = nobs, "boost_summary" = bsum, "runtime" = elapsed))
}


reverse_edf <- function(x, bn, bmat, nobs, y, eta, approx = TRUE)
{
  beta <- bn + apply(bmat, 2, sum)

  fit <- x$X %*% beta
  y <- y + fit

  tX <- t(x$X)
  XX <- crossprod(x$X)

  objfun <- function(tau2) {
    if(!x$fixed) {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1/tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](beta) else x$S[[j]]
    } else {
      S <- 1/tau2 * diag(1, ncol(x$X))
    }
    beta2 <- matrix_inv(XX + S, index = x$sparse.setup) %*% tX %*% y
    mean((fit - x$X %*% beta2)^2)
  }

  tau2 <- tau2.optim(objfun, start = x$boost.tau2, maxit = 100)

  if(!x$fixed) {
    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1/tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](beta) else x$S[[j]]
  } else {
    I <- diag(1, ncol(x$X))
    S <- 1/tau2 * I
  }

  P <- matrix_inv(XX + S, index = x$sparse.setup)
  edf <- sum_diag(XX %*% P)

  return(list("edf" = edf - x$state$edf, "tau2" = tau2, "fedf" = edf))
}


unique.fuse <- function(x, digits = 4) {
  unique(round(x, digits = digits))
}


if(FALSE) {
  n <- 1000
  d <- data.frame("x" = runif(n, -3, 3))
  d$y <- 1.2 + sin(d$x) + rnorm(n, sd = 0.3)
  plot(d)

  b1 <- bamlss(y ~ s(x,k=50), data = d, sampler = FALSE, optimizer = boost, stop.criterion = "AIC", reverse = TRUE)
  b1 <- bamlss(y ~ s(x,k=50), data = d, sampler = FALSE, optimizer = boost, stop.criterion = "AIC", reverse = FALSE)

  d$p1 <- predict(b1, model = "mu")
  d$p2 <- predict(b1, model = "mu")

  plot2d(p1 ~ x, data = d)
  plot2d(p2 ~ x, data = d, add = TRUE, col.lines = 4)
  plot2d(I(1.2 + sin(x)) ~ x, data = d, add = TRUE, col.lines = 2)

  b1 <- bamlss(y ~ s(x,k=50), data = d, sampler = FALSE)
}

selfun <- function(iter, i, j, state, parm, x, family, sfun, yname, weights, selectmodel = TRUE)
{
  if(is.null(selectmodel))
    selectmodel <- FALSE
  parm[[i]][[j]][iter, ] <- get.par(state$parameters, "b")
  parm <- parm2mat(parm, mstop = iter, fixed = iter)
  if(!selectmodel) {
    formula <- list()
    for(i in names(x))
      formula[[i]] <- x[[i]][c("formula", "fake.formula")]
    class(formula) <- c("bamlss.formula", "list")
    environment(formula) <- environment(formula[[1]]$formula)
    attr(formula, "response.name") <- yname
    m <- list("formula" = formula, "x" = x, "family" = family, "parameters" = parm)
    class(m) <- c("bamlss", "bamlss.frame", "list")
    return(sfun(m))
  } else {
    return(sfun(parm))
  }
}


boost_frame <- function(formula, train, test, family = "gaussian", ...)
{
  if(!all(names(test) == names(train)))
    stop("test and training data must contain the same variables!")

  bf <- bamlss.frame(formula, data = train, family = family, ...)

  for(i in names(bf$x)) {
    for(j in seq_along(bf$x[[i]]$smooth.construct)) {
      if(!inherits(bf$x[[i]]$smooth.construct[[j]], "no.mgcv") & !inherits(bf$x[[i]]$smooth.construct[[j]], "special")) {
        if(!is.null(bf$x[[i]]$smooth.construct[[j]]$is.refund)) {
          rfcall <- bf$x[[i]]$smooth.construct[[j]]$refund.call
          tfm <- eval(parse(text = rfcall), envir = test)
          tfme <- eval(tfm$call, envir = tfm$data)
          bf$x[[i]]$smooth.construct[[j]]$X <- smoothCon(tfme, data = tfm$data, n = nrow(tfm$data[[1L]]),
            knots = NULL, absorb.cons = TRUE)[[1]]$X
          rm(tfm)
          rm(tfme)
        } else {
          bf$x[[i]]$smooth.construct[[j]]$X <- PredictMat(bf$x[[i]]$smooth.construct[[j]], test)
        }
      } else {
        if(is.null(bf$x[[i]]$smooth.construct[[j]]$PredictMat)) {
          bf$x[[i]]$smooth.construct[[j]]$X <- PredictMat(bf$x[[i]]$smooth.construct[[j]], test)
        } else {
          bf$x[[i]]$smooth.construct[[j]]$X <- bf$x[[i]]$smooth.construct[[j]]$PredictMat(bf$x[[i]]$smooth.construct[[j]], test)
        }
      }
    }
    if(!is.null(bf$x[[i]]$model.matrix)) {
      sc <- attr(bf$x[[i]]$model.matrix, "scale")
      bf$x[[i]]$model.matrix <- model.matrix(drop.terms.bamlss(bf$x[[i]]$terms,
        sterms = FALSE, keep.response = FALSE, data = test), data = test)
      if(ncol(bf$x[[i]]$model.matrix) > 0) {
        if(!is.null(sc)) {
          for(name in unique(unlist(lapply(sc, names)))) {
            bf$x[[i]]$model.matrix[,name] <- (bf$x[[i]]$model.matrix[,name] - sc$center[name] ) / sc$scale[name]
          }
        }
      } else bf$x[[i]]$model.matrix <- NULL
    }
  }

  yname <- names(bf$y)
  family <- bf$family

  bf <- opt_boost(x = bf$x, y = bf$y, family = bf$family,
    weights = model.weights(bf$model.frame),
    offset = model.offset(bf$model.frame), ret.x = TRUE, initialize = FALSE, ...)

  formula <- list()
  for(i in names(bf))
    formula[[i]] <- bf[[i]][c("formula", "fake.formula")]
  class(formula) <- c("bamlss.formula", "list")
  environment(formula) <- environment(formula[[1]]$formula)
  attr(formula, "response.name") <- yname

  bf <- list("formula" = formula, "x" = bf, "family" = family)
  class(bf) <- c("boost_frame", "list")

  bf
}

predict.boost_frame <- function(object, type = c("link", "parameter"), ...)
{
  type <- match.arg(type)
  object$x <- set.starting.values(object$x, object$parameters)
  fit <- get.eta(object$x, expand = TRUE)
  if(type == "parameter")
    fit <- object$family$map2par(fit)
  return(fit)
}


## Updating the hat-matrix.
hatmat_trace <- function(H0, H1)
{
  .Call("hatmat_trace", H0, H1, PACKAGE = "bamlss")
}

hatmat_sumdiag <- function(H)
{
  .Call("hatmat_sumdiag", H, PACKAGE = "bamlss")
}


## Boost setup.
boost_transform <- function(x, y, df = NULL, family,
  weights = NULL, offset = NULL, maxit = 100,
  eps = .Machine$double.eps^0.25, initialize = TRUE,
  nu = 0.1, nu.adapt = TRUE, ...)
{
  np <- length(x)
  nx <- names(x)

  ## Initialize select indicator and intercepts.
  for(j in 1:np) {
    nid <- NULL
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      if(!is.null(df) & !inherits(x[[nx[j]]]$smooth.construct[[sj]], "randombits.smooth") & !inherits(x[[nx[j]]]$smooth.construct[[sj]], "nnet.smooth") & !inherits(x[[nx[j]]]$smooth.construct[[sj]], "nnet2.smooth")) {
        if(inherits(x[[nx[j]]]$smooth.construct[[sj]], "lasso.smooth"))
          x[[nx[j]]]$smooth.construct[[sj]]$xt$df <- df
        x[[nx[j]]]$smooth.construct[[sj]] <- assign.df(x[[nx[j]]]$smooth.construct[[sj]], df, do.part = TRUE)
      }
      if(!is.null(x[[nx[j]]]$smooth.construct[[sj]]$fxsp)) {
        if(!x[[nx[j]]]$smooth.construct[[sj]]$fxsp & !x[[nx[j]]]$smooth.construct[[sj]]$fixed) {
          x[[nx[j]]]$smooth.construct[[sj]]$old.optimize <- x[[nx[j]]]$smooth.construct[[sj]]$state$do.optim
          x[[nx[j]]]$smooth.construct[[sj]]$state$do.optim <- FALSE
          x[[nx[j]]]$smooth.construct[[sj]]$do.optim <- FALSE
        }
      }
    }
    if(has_pterms(x[[nx[j]]]$terms)) {
      ii <- which(names(x[[nx[j]]]$smooth.construct) == "model.matrix")
      model.matrix <- list()
      cn <- colnames(x[[nx[j]]]$smooth.construct[[ii]]$X)
      g0 <- get.par(x[[nx[j]]]$smooth.construct[[ii]]$state$parameters, "b")
      nm <- NULL
      assign <- attr(x[[nx[j]]]$smooth.construct[[ii]]$X, "assign")
      for(pj in 1:ncol(x[[nx[j]]]$smooth.construct[[ii]]$X)) {
        model.matrix[[pj]] <- list()
        model.matrix[[pj]]$label <- cn[pj]
        model.matrix[[pj]]$term <- cn[pj]
        model.matrix[[pj]]$X <- x[[nx[j]]]$smooth.construct[[ii]]$X[, pj, drop = FALSE]
        model.matrix[[pj]]$binning <- x[[nx[j]]]$smooth.construct[[ii]]$binning
        model.matrix[[pj]]$nobs <- x[[nx[j]]]$smooth.construct[[ii]]$nobs
        model.matrix[[pj]]$fixed <- TRUE
        model.matrix[[pj]]$fxsp <- FALSE
        model.matrix[[pj]]$weights <- x[[nx[j]]]$smooth.construct[[ii]]$weights
        model.matrix[[pj]]$rres <- x[[nx[j]]]$smooth.construct[[ii]]$rres
        model.matrix[[pj]]$fit.reduced <- x[[nx[j]]]$smooth.construct[[ii]]$fit.reduced
        model.matrix[[pj]]$fit.fun <- x[[nx[j]]]$smooth.construct[[ii]]$fit.fun
        model.matrix[[pj]]$state <- list("parameters" = g0[pj])
        model.matrix[[pj]]$state$fitted.values <- drop(model.matrix[[pj]]$X %*% g0[pj])
        if(!is.null(model.matrix[[pj]]$binning$match.index))
          model.matrix[[pj]]$state$fitted.values <- model.matrix[[pj]]$state$fitted.values[model.matrix[[pj]]$binning$match.index]
        model.matrix[[pj]]$state$edf <- 0
        model.matrix[[pj]]$state$rss <- 0
        model.matrix[[pj]]$state$do.optim <- FALSE
        model.matrix[[pj]]$is.model.matrix <- TRUE
        model.matrix[[pj]]$selected <- rep(0, length = maxit)
        model.matrix[[pj]]$sparse.setup <- sparse.setup(model.matrix[[pj]]$X, S = model.matrix[[pj]]$S)
        model.matrix[[pj]]$upper <- Inf
        model.matrix[[pj]]$lower <- -Inf
        model.matrix[[pj]]$assign <- assign[pj]
        class(model.matrix[[pj]]) <- class(x[[nx[j]]]$smooth.construct[[ii]])
      }
      names(model.matrix) <- cn
      x[[nx[j]]]$smooth.construct[[ii]] <- NULL
      x[[nx[j]]]$smooth.construct <- c(model.matrix, x[[nx[j]]]$smooth.construct)
      attr(x[[nx[j]]], "assign") <- assign
    }
  }

  always3 <- list(...)$always3
  if(is.null(always3))
    always3 <- FALSE
  
  ## Save more info.
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      if(always3 & (x[[nx[j]]]$smooth.construct[[sj]]$label != "(Intercept)")) {
        x[[nx[j]]]$smooth.construct[[sj]]$X <- cbind(1, x[[nx[j]]]$smooth.construct[[sj]]$X)
        x[[nx[j]]]$smooth.construct[[sj]]$with.itcpt <- TRUE
        x[[nx[j]]]$smooth.construct[[sj]]$state$parameters <- c("b0" = 0, x[[nx[j]]]$smooth.construct[[sj]]$state$parameters)
      }
      x[[nx[j]]]$smooth.construct[[sj]]$state$init.edf <- x[[nx[j]]]$smooth.construct[[sj]]$state$edf
      x[[nx[j]]]$smooth.construct[[sj]]$state$edf <- 0
      nc <- ncol(x[[nx[j]]]$smooth.construct[[sj]]$X)
      nr <- nrow(x[[nx[j]]]$smooth.construct[[sj]]$X)
      x[[nx[j]]]$smooth.construct[[sj]]$XWX <- matrix(0, nc, nc)
      x[[nx[j]]]$smooth.construct[[sj]]$XW <- matrix(0, nc, nr)
      x[[nx[j]]]$smooth.construct[[sj]]$selected <- rep(0, length = maxit)
      x[[nx[j]]]$smooth.construct[[sj]]$loglik <- rep(0, length = maxit)
      x[[nx[j]]]$smooth.construct[[sj]]$state$rss <- 0
      if(is.null(x[[nx[j]]]$smooth.construct[[sj]]$is.model.matrix))
        x[[nx[j]]]$smooth.construct[[sj]]$boost.tau2 <- get.par(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, "tau2")
      else
        x[[nx[j]]]$smooth.construct[[sj]]$boost.tau2 <- 1000
      if(!is.null(x[[nx[j]]]$smooth.construct[[sj]]$S))
        x[[nx[j]]]$smooth.construct[[sj]]$penaltyFunction <- as.integer(sapply(x[[nx[j]]]$smooth.construct[[sj]]$S, is.function))
      else
        x[[nx[j]]]$smooth.construct[[sj]]$penaltyFunction <- 0L
      if(inherits(x[[nx[j]]]$smooth.construct[[sj]], "nnet.smooth") | inherits(x[[nx[j]]]$smooth.construct[[sj]], "nnet2.smooth"))
        x[[nx[j]]]$smooth.construct[[sj]]$fuse <- FALSE
    }
  }
  
  if(initialize) {
    eta <- get.eta(x)
    eta <- init.eta(eta, y, family, nobs)
    if(!is.null(offset)) {
      offset <- as.data.frame(offset)
      for(j in nx) {
        if(!is.null(offset[[j]]))
          eta[[j]] <- eta[[j]] + offset[[j]]
      }
    }
    nobs <- length(eta[[1]])
    start <- unlist(lapply(eta, mean, na.rm = TRUE))

    par <- rep(0, length(nx))
    names(par) <- nx

    eps0 <- eps + 1L
    k2 <- 0
    while((eps < eps0) & (k2 < 100)) {
      eta0 <- eta
      for(i in nx) {
        ll0 <- family$loglik(y, family$map2par(eta))

        peta <- family$map2par(eta)
        hess <- process.derivs(family$hess[[i]](y, peta, id = i), is.weight = TRUE)
        score <- process.derivs(family$score[[i]](y, peta, id = i), is.weight = FALSE)
        z <- eta[[i]] + 1 / hess * score

        b0 <- par[i]

        if(!is.null(weights)) {
          if(attr(weights, "identical"))
            hess <- hess * as.numeric(weights[, 1])
        }

        par[i] <- 1 / (sum(hess) + 1e-20) * sum(hess * z, na.rm = TRUE)

        eta[[i]] <- rep(par[i], nobs)

        ll1 <- family$loglik(y, family$map2par(eta))

        if(ll1 < ll0) {
          fnu <- function(nu2) {
            b <- nu2 * par[i] + (1 - nu2) * b0
            eta[[i]] <- rep(b, nobs)
            return(family$loglik(y, family$map2par(eta)))
          }
          nu2 <- 1
          while((fnu(nu2) < ll0) & (.Machine$double.eps < nu2)) {
            nu2 <- nu2 / 2
          }
          par[i] <- nu2 * par[i] + (1 - nu2) * b0
          eta[[i]] <- rep(par[i], nobs)
        }
      }
      eps0 <- do.call("cbind", eta)
      eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)

      k2 <- k2 + 1
    }
    
    for(i in nx) {
      if(!is.null(x[[i]]$smooth.construct[["(Intercept)"]])) {
        x[[i]]$smooth.construct[["(Intercept)"]]$state$parameters[1] <- par[i]
        x[[i]]$smooth.construct[["(Intercept)"]]$state$fitted.values <- rep(par[i], length = nobs)
        x[[i]]$smooth.construct[["(Intercept)"]]$state$edf <- 1
        x[[i]]$smooth.construct[["(Intercept)"]]$state$init.edf <- 1
      }
    }
  }
  
  return(x)
}


boost_fit_nnet <- function(nu, X, N, y, ind, nthreads = NULL)
{
  if(is.null(nthreads))
    nthreads <- 1L
  .Call("boost_fit_nnet", nu, X, N, y, ind, as.integer(nthreads))
}


## Simple list() generator for
## saving states of model terms.
make.state.list <- function(x, type = 1, intercept = TRUE)
{
  elmts <- c("formula", "fake.formula")
  if(all(elmts %in% names(x))) {
    rval <- list()
    if(!is.null(x$model.matrix))
      rval$model.matrix <- NA
    if(!is.null(x$smooth.construct)) {
      for(j in names(x$smooth.construct)) {
        if(j == "(Intercept)" & intercept)
          rval[[j]] <- NA
        if(j != "(Intercept)")
          rval[[j]] <- NA
      }
    }
    if(type > 1)
      rval <- unlist(rval)
  } else {
    rval <- list()
    for(j in names(x)) {
      rval[[j]] <- make.state.list(x[[j]], type, intercept = intercept)
    }
  }
  return(rval)
}


make.par.list <- function(x, iter)
{
  elmts <- c("formula", "fake.formula")
  if(all(elmts %in% names(x))) {
    rval <- list()
    if(!is.null(x$smooth.construct)) {
      for(j in names(x$smooth.construct)) {
        rval[[j]] <- if(is.null(x$smooth.construct[[j]]$special.npar)) {
          matrix(0, nrow = iter, ncol = ncol(x$smooth.construct[[j]]$X))
        } else {
          matrix(0, nrow = iter, ncol = x$smooth.construct[[j]]$special.npar)
        }
        colnames(rval[[j]]) <- names(get.par(x$smooth.construct[[j]]$state$parameters, "b"))
        rval[[j]][1, ] <- get.par(x$smooth.construct[[j]]$state$parameters, "b")
        if(!is.null(x$smooth.construct[[j]]$with.itcpt)) {
          if(length(i <- grep("b0", colnames(rval[[j]]))))
            rval[[j]] <- rval[[j]][, -i, drop = FALSE]
        }
        if(!is.null(x$smooth.construct[[j]]$is.model.matrix))
          attr(rval[[j]], "is.model.matrix") <- TRUE
        if(inherits(x$smooth.construct[[j]], "nnet.smooth"))
          class(rval[[j]]) <- c(class(rval[[j]]), "nnet.smooth")
      }
    }
  } else {
    rval <- list()
    for(j in names(x)) {
      rval[[j]] <- make.par.list(x[[j]], iter)
    }
  }
  return(rval)
}

parm2mat <- function(x, mstop, fixed = NULL)
{
  nx <- names(x)
  for(i in seq_along(x)) {
    is.mm <- NULL
    for(j in names(x[[i]])) {
      if(!is.null(attr(x[[i]][[j]], "is.model.matrix")))
        is.mm <- c(is.mm, j)
      cn <- colnames(x[[i]][[j]])
      if(!inherits(x[[i]][[j]], "nnet.smooth")) {
        x[[i]][[j]] <- apply(x[[i]][[j]][1:mstop, , drop = FALSE], 2, cumsum)
      } else {
        x[[i]][[j]] <- x[[i]][[j]][1:mstop, , drop = FALSE]
      }
      if(!is.matrix(x[[i]][[j]]))
        x[[i]][[j]] <- matrix(x[[i]][[j]], ncol = length(cn))
      colnames(x[[i]][[j]]) <- cn
    }
    if(!is.null(is.mm)) {
      x[[i]][["p"]] <- do.call("cbind", x[[i]][is.mm])
      colnames(x[[i]][["p"]]) <- is.mm
      x[[i]][is.mm[is.mm != "p"]] <- NULL
    }
    sm <- names(x[[i]])
    sm <- sm[sm != "p"]
    if(length(sm)) {
      x[[i]][["s"]] <- x[[i]][sm]
      x[[i]][sm[sm != "s"]] <- NULL
    }
    n <- names(x[[i]])
    for(j in names(x[[i]])) {
      if(j != "s") {
        colnames(x[[i]][[j]]) <- paste(nx[i], j, colnames(x[[i]][[j]]), sep = ".")
      } else {
        for(k in names(x[[i]][[j]])) {
          colnames(x[[i]][[j]][[k]]) <- paste(nx[i], j, k, colnames(x[[i]][[j]][[k]]), sep = ".")
        }
        x[[i]][[j]] <- do.call("cbind", x[[i]][[j]])
      }
    }
    x[[i]] <- do.call("cbind", x[[i]])
  }
  x <- do.call("cbind", x)
  if(!is.null(fixed))
    x <- x[fixed, ]
  return(x)
}


## Retransform 'x' to 'bamlss.frame' structure.
boost.retransform <- function(x) {
  for(i in names(x)) {
    if(has_pterms(x[[i]]$terms)) {
      state <- list()
      X <- drop <- xscales <- NULL
      for(j in names(x[[i]]$smooth.construct)) {
        if(inherits(x[[i]]$smooth.construct[[j]], "model.matrix")) {
          drop <- c(drop, j)
          b <- get.par(x[[i]]$smooth.construct[[j]]$state$parameters, "b")
          X <- cbind(X, x[[i]]$smooth.construct[[j]]$X)
          state$parameters <- c(state$parameters, b)
        }
      }
      label <- paste(drop, collapse = "+")
      binning <- x[[i]]$smooth.construct[[drop[1]]]$binning
      state$fitted.values <- drop(X %*% state$parameters)
      x[[i]]$smooth.construct[drop] <- NULL
      x[[i]]$smooth.construct$model.matrix <- list(
        "X" = X,
        "S" = list(diag(0, ncol(X))),
        "rank" = ncol(X),
        "term" = label,
        "label" = label,
        "bs.dim" = ncol(X),
        "fixed" = TRUE,
        "is.model.matrix" = TRUE,
        "by" = "NA",
        "xt" = list("binning" = binning),
        "state" = state
      )
      x[[i]]$smooth.construct$model.matrix$fit.fun <- make.fit.fun(x[[i]]$smooth.construct$model.matrix)
    }
  }
  return(x)
}


## Boosting iwls.
boost_iwls <- function(x, hess, resids, nu)
{
  ## Initial parameters and fit.
  g0 <- get.par(x$state$parameters, "b")
  fit0 <- fitted(x$state)
  
  ## Compute reduced residuals.
  xbin.fun(x$binning$sorted.index, hess, resids, x$weights, x$rres, x$binning$order)
  
  ## Compute mean and precision.
  XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)
  if(x$fixed) {
    P <- matrix_inv(XWX, index = x$sparse.setup)
  } else {
    S <- 0
    tau2 <- get.state(x, "tau2")
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    P <- matrix_inv(XWX + S, index = x$sparse.setup)
  }
  
  ## New parameters.
  g <- nu * drop(P %*% crossprod(x$X, x$rres))
  
  ## Finalize.
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))
  
  ## Find edf.
  xbin.fun(x$binning$sorted.index, hess, resids + fit0 + fitted(x$state), x$weights, x$rres, x$binning$order)
  
  XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)
  if(x$fixed) {
    P <- matrix_inv(XWX, index = x$sparse.setup)
  } else {
    g0 <- g0 + g
    
    objfun <- function(tau2) {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      P <- matrix_inv(XWX + S, index = x$sparse.setup)
      g1 <- drop(P %*% crossprod(x$X, x$rres))
      sum((g1 - g0)^2)
    }
    
    if(length(get.state(x, "tau2")) < 2) {
      tau2 <- optimize(objfun, interval = x$state$interval)$minimum
    } else {
      i <- grep("tau2", names(x$lower))
      tau2 <- if(!is.null(x$state$true.tau2)) x$state$true.tau2 else get.state(x, "tau2")
      opt <- try(optim(tau2, fn = objfun, method = "L-BFGS-B",
                       lower = x$lower[i], upper = x$upper[i]), silent = TRUE)
      if(!inherits(opt, "try-error"))
        tau2 <- opt$par
    }
    if(inherits(tau2, "try-error"))
      stop(paste("problem in finding optimum smoothing parameter for term ", x$label, "!", sep = ""))
    
    attr(x$state$parameters, "true.tau2") <- tau2
    
    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    P <- matrix_inv(XWX + S, index = x$sparse.setup)
  }
  
  ## Assign degrees of freedom.
  x$state$edf <- sum_diag(XWX %*% P)
  attr(x$state$parameters, "edf") <- x$state$edf
  
  return(x$state)
}


## Boosting gradient fit.
boost_fit <- function(x, y, nu, hatmatrix = TRUE, weights = NULL, ...)
{
  ## process weights.
  if(is.null(weights))
    weights <- rep(1, length = length(y))

  ## Compute reduced residuals.
  xbin.fun(x$binning$sorted.index, weights, y, x$weights, x$rres, x$binning$order)
  
  ## Compute mean and precision.
  XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)

  if(x$fixed) {
    P <- matrix_inv(XWX, index = x$sparse.setup)
  } else {
    S <- 0
    tau2 <- get.state(x, "tau2")
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](x$state$parameters) else x$S[[j]]
    P <- matrix_inv(XWX + S, index = x$sparse.setup)
  }
  
  ## New parameters.
  g <- nu * drop(P %*% crossprod(x$X, x$rres))
  
  ## Finalize.
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))

if(any(is.na(x$state$fitted.values))) {
  stop("why?")
}

  x$state$rss <- sum((x$state$fitted.values - y)^2 * weights)

  if(hatmatrix)
    x$state$hat <- nu * x$X %*% P %*% t(x$X)
  
  return(x$state)
}


## Increase coefficients.
increase <- function(state0, state1)
{
  g <- get.par(state0$parameters, "b") + get.par(state1$parameters, "b")
  state0$fitted.values <- fitted(state0) + fitted(state1)
  state0$parameters <- set.par(state0$parameters, g, "b")
  state0$edf <- state1$edf
  state0$parameters <- set.par(state0$parameters, get.par(state1$parameters, "tau2"), "tau2")
  attr(state0$parameters, "true.tau2") <- attr(state1$parameters, "true.tau2")
  attr(state0$parameters, "edf") <- attr(state1$parameters, "edf")
  state0$special <- state1$special
  state0$hat <- state1$hat
  if(!is.null(state1$redf)) {
    state0$boost.tau2 <- state1$redf$tau2
    state0$edf <- state1$redf$fedf
  }
  state0
}


## Extract number of selected base learners
get.qsel <- function(x, iter, qsel.splitfactor = FALSE)
{
  rval <- 0
  for(i in names(x)) {
    assign  <- as.character(attr(x[[i]], "assign"))
    if(!length(assign))
      return(1)
    uassign <- unique(assign)
    facID   <- sapply(uassign, function(x) { sum(assign == x) > 1 })
    assign  <- uassign[facID]
    asssel  <- list()
    for(m in assign)
      asssel[[m]] <- 0

    for(j in names(x[[i]]$smooth.construct)) {
      if(inherits(x[[i]]$smooth.construct[[j]], "linear.smooth") | inherits(x[[i]]$smooth.construct[[j]], "randombits.smooth") | inherits(x[[i]]$smooth.construct[[j]], "nnet2.smooth")) {
        np <- names(x[[i]]$smooth.construct[[j]]$state$parameters)
        rval <- rval + sum(abs(x[[i]]$smooth.construct[[j]]$state$parameters[grep("b", np)]) > 1e-10)
        next
      }
      rval <- rval + 1 * (any(x[[i]]$smooth.construct[[j]]$selected[1:iter] > 0) &
                          j != "(Intercept)")

      if(!is.null(x[[i]]$smooth.construct[[j]]$assign)) {
        m <- as.character(x[[i]]$smooth.construct[[j]]$assign)
        if(m %in% assign) {
          asssel[[m]] <- asssel[[m]] +
            1 * any(x[[i]]$smooth.construct[[j]]$selected[1:iter] > 0)
        }
      }
    }
    
    if(!qsel.splitfactor) {
        for(m in assign) {
            rval <- rval - max(0, asssel[[m]] - 1)
        }
    }

    rm(asssel)
  }
  rval
}

get.maxq <- function(x)
{
  rval <- 0
  for(i in names(x)) {
    rval <- rval + length(names(x[[i]]$smooth.construct))
  }
  rval
}


## Extract summary for boosting.
make.boost_summary <- function(x, mstop, save.ic, edf, hatmatrix, nobs)
{
  nx <- names(x)
  labels <- NULL
  ll.contrib <- crit.contrib <- NULL
  bsum <- lmat <- list()
  for(i in nx) {
    rn <- NULL
    for(j in names(x[[i]]$smooth.construct)) {
      labels <- c(labels, paste(x[[i]]$smooth.construct[[j]]$label, i, sep = "."))
      rn <- c(rn, x[[i]]$smooth.construct[[j]]$label)
      bsum[[i]] <- rbind(bsum[[i]], sum(x[[i]]$smooth.construct[[j]]$selected[1:mstop]) / mstop * 100)
      lmat[[i]] <- rbind(lmat[[i]], sum(x[[i]]$smooth.construct[[j]]$loglik[1:mstop]))
      ll.contrib <- cbind(ll.contrib, cumsum(x[[i]]$smooth.construct[[j]]$loglik[1:mstop]))
      if(!is.null(x[[i]]$smooth.construct[[j]]$criterion))
        crit.contrib <- cbind(crit.contrib, cumsum(x[[i]]$smooth.construct[[j]]$criterion[1:mstop]))
    }
    if(!is.matrix(bsum[[i]])) bsum[[i]] <- matrix(bsum[[i]], nrow = 1)
    bsum[[i]] <- cbind(bsum[[i]], lmat[[i]])
    if(!is.matrix(bsum[[i]])) bsum[[i]] <- matrix(bsum[[i]], nrow = 1)
    colnames(bsum[[i]]) <- c(paste(i, "% selected"), "LogLik contrib.")
    rownames(bsum[[i]]) <- rownames(lmat[[i]]) <- rn
    bsum[[i]] <- bsum[[i]][order(bsum[[i]][, 2], decreasing = TRUE), , drop = FALSE]
  }
  colnames(ll.contrib) <- labels
  if(!is.null(crit.contrib))
    colnames(crit.contrib) <- labels
  names(bsum) <- nx
  bsum <- list("summary" = bsum, "mstop" = mstop,
    "ic" = save.ic[1:mstop], "loglik" = ll.contrib)
  if(hatmatrix) {
    bsum$criterion <- list()
    bsum$criterion$bic <- -2 * bsum$ic + edf[1:mstop] * log(nobs)
    bsum$criterion$aic <- -2 * bsum$ic + edf[1:mstop] * 2
    bsum$criterion$edf <- edf[1:mstop]
  }
  if(!is.null(crit.contrib))
    bsum$crit.contrib <- crit.contrib
  class(bsum) <- "boost_summary"
  return(bsum)
}


boost_summary <- function(object, ...)
{
  if(!is.null(object$model.stats$optimizer$boost_summary))
    print.boost_summary(object$model.stats$optimizer$boost_summary, ...)
  invisible(object$model.stats$optimizer$boost_summary)
}


## Smallish print function for boost summaries.
print.boost_summary <- function(x, summary = TRUE, plot = TRUE,
  which = c("loglik", "loglik.contrib"), intercept = TRUE,
  spar = TRUE, ...)
{
  if(inherits(x, "bamlss"))
    x <- x$model.stats$optimizer$boost_summary
  if(is.null(x))
    stop("no summary for boosted model available")
  if(summary) {
    np <- length(x$summary)
    cat("\n")
    cat("logLik. =", if(is.na(x$ic[x$mstop])) x$ic[x$mstop - 1] else x$ic[x$mstop], "-> at mstop =", x$mstop, "\n---\n")
    for(j in 1:np) {
      if(length(x$summary[[j]]) < 2) {
        print(round(x$summary[[j]], digits = 4))
      } else printCoefmat(x$summary[[j]], digits = 4)
      if(j != np)
        cat("---\n")
    }
    cat("\n")
  }
  
  if(plot) {
    if(!is.character(which)) {
      which <- c("loglik", "loglik.contrib", "parameters", "aic", "bic", "user")[as.integer(which)]
    } else {
      which <- tolower(which)
      which <- match.arg(which, c("loglik", "loglik.contrib", "parameters", "aic", "bic", "user"), several.ok = TRUE)
    }
    
    if(spar) {
      op <- par(no.readonly = TRUE)
      on.exit(par(op))
      par(mfrow = c(1, length(which)))
    }
    
    for(w in which) {
      if(w == "loglik") {
        if(spar)
          par(mar = c(5.1, 4.1, 2.1, 2.1))
        plot(x$ic, type = "l", xlab = "Iteration", ylab = "logLik", ...)
        abline(v = x$mstop, lwd = 3, col = "lightgray")
        axis(3, at = x$mstop, labels = paste("mstop =", x$mstop))
      }
      if(w == "loglik.contrib") {
        if(spar)
          par(mar = c(5.1, 4.1, 2.1, 10.1))
        if(!intercept) {
          j <- grep("(Intercept)", colnames(x$loglik), fixed = TRUE)
          x$loglik <- x$loglik[, -j]
        }
        args <- list(...)
        if(!is.null(args$name)) {
          x$loglik <- x$loglik[, grep2(args$name, colnames(x$loglik), fixed = TRUE), drop = FALSE]
        }
        xn <- sapply(strsplit(colnames(x$loglik), ".", fixed = TRUE), function(x) { x[length(x)] })
        if(is.null(cols <- args$mcol)) {
          cols <- rainbow_hcl(length(unique(xn)))
        } else {
          cols <- rep(cols, length.out = length(unique(xn)))
        }
        matplot(x$loglik, type = "l", lty = 1,
          xlab = "Iteration", ylab = "LogLik contribution", col = cols[as.factor(xn)],
          lwd = args$lwd)
        abline(v = x$mstop, lwd = 3, col = "lightgray")
        cn <- colnames(x$loglik)
        if(!is.null(args$drop)) {
          for(dn in args$drop)
            cn <- gsub(dn, "", cn, fixed = TRUE)
        }
        axis(4, at = x$loglik[nrow(x$loglik), ], labels = cn, las = 1)
        axis(3, at = x$mstop, labels = paste("mstop =", x$mstop))
      }
      if(w %in% c("aic", "bic", "user")) {
        if(!is.null(x$criterion)) {
          if(spar)
            par(mar = c(5.1, 4.1, 2.1, 2.1))
          args <- list()
          if(is.null(args$xlab))
            args$xlab <- "Iteration"
          if(is.null(args$ylab))
            args$ylab <- if(w == "user") "User IC" else toupper(w)
          plot(x$criterion[[if(w == "user") "userIC" else w]], type = "l", xlab = args$xlab, ylab = args$ylab)
          i <- which.min(x$criterion[[w]])
          abline(v = i, lwd = 3, col = "lightgray")
          if(!is.null(x$criterion$edf))
            axis(3, at = i, labels = paste("mstop = ", i, ", edf = ", round(x$criterion$edf[i], digits = 2), sep = ""))
        }
      }
    }
  }
  
  return(invisible(x))
}


plot.boost_summary <- function(x, ...)
{
  print.boost_summary(x, summary = FALSE, plot = TRUE, ...) 
}

boost_plot <- function(x, which = c("loglik", "loglik.contrib", "parameters", "aic", "bic", "user"),
  intercept = TRUE, spar = TRUE, mstop = NULL, name = NULL, drop = NULL, labels = NULL, color = NULL, ...)
{
  if(!is.character(which)) {
    which <- c("loglik", "loglik.contrib", "parameters")[as.integer(which)]
  } else {
    which <- tolower(which)
    which <- match.arg(which, several.ok = TRUE)
  }
  
  if(spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = c(1, length(which)))
  }
  
  if(is.null(mstop))
    mstop <- x$model.stats$optimizer$boost_summary$mstop
  x$model.stats$optimizer$boost_summary$mstop <- mstop
  x$model.stats$optimizer$boost_summary$ic <- x$model.stats$optimizer$boost_summary$ic[1:mstop]
  x$model.stats$optimizer$boost_summary$loglik <- x$model.stats$optimizer$boost_summary$loglik[1:mstop, , drop = FALSE]
  
  for(w in which) {
    if(w %in% c("loglik", "loglik.contrib", "aic", "bic", "user")) {
      if((w == "loglik") & spar)
        par(mar = c(5.1, 4.1, 2.1, 2.1))
      if((w == "loglik.contrib") & spar)
        par(mar = c(5.1, 4.1, 2.1, 10.1))
      plot.boost_summary(x, which = w, spar = FALSE, intercept = intercept, name = name, ...)
    }
    if(w == "parameters") {
      if(spar)
        par(mar = c(5.1, 4.1, 2.1, 10.1))
      if(!is.null(drop)) {
        x$parameters <- x$parameters[, -grep2(drop, colnames(x$parameters), fixed = TRUE), drop = FALSE]
      }
      if(!is.null(name)) {
        x$parameters <- x$parameters[, grep2(name, colnames(x$parameters), fixed = TRUE), drop = FALSE]
      }
      
      p <- x$parameters[1:mstop, , drop = FALSE]
      if(!intercept)
        p <- p[, -grep("(Intercept)", colnames(p), fixed = TRUE), drop = FALSE]
      
      xn <- sapply(strsplit(colnames(x$parameters), ".", fixed = TRUE), function(x) { x[1] })
      if(length(unique(xn)) < 2)
        xn <- sapply(strsplit(colnames(x$parameters), ".", fixed = TRUE), function(x) { x[3] })
      
      cols <- if(is.null(color)) {
        if(length(unique(xn)) < 2) "black" else rainbow_hcl(length(unique(xn)))
      } else {
        if(is.function(color)) {
          color(length(unique(xn)))
        } else {
          rep(color, length.out = length(unique(xn)))
        }
      }

      if(is.null(labels)) {
        labs <- labs0 <- colnames(p)
        plab <- p[nrow(p), ]
        o <- order(plab, decreasing = TRUE)
        labs <- labs[o]
        plab <- plab[o]
        rplab <- diff(range(plab))
        for(i in 1:(length(plab) - 1)) {
          dp <- abs(plab[i] - plab[i + 1]) / rplab
          if(length(dp) < 1)
            dp <- 0
          if(is.na(dp))
            dp <- 0
          if(dp <= 0.02) {
            labs[i + 1] <- paste(c(labs[i], labs[i + 1]), collapse = ",")
            labs[i] <- ""
          }
        }
        labs <- labs[order(o)]
        if(!is.null(name)) {
          for(j in seq_along(name))
            labs <- gsub(name[j], "", labs, fixed = TRUE)
        }
      } else labs <- rep(labels, length.out = ncol(p))
      at <- p[nrow(p), ]
      at <- at[labs != ""]
      labs <- labs[labs != ""]
      matplot(p, type = "l", lty = 1, col = cols[as.factor(xn)], xlab = "Iteration", ...)
      abline(v = mstop, lwd = 3, col = "lightgray")
      axis(4, at = at, labels = labs, las = 1)
      axis(3, at = mstop, labels = paste("mstop =", mstop))
    }
  }
}


## Assign starting values.
set.starting.values <- function(x, start)
{
  if(!is.null(start)) {
    if(is.list(start)) {
      if("parameters" %in% names(start))
        start <- start$parameters
    }
    if(is.list(start))
      start <- unlist(start)
    if(is.matrix(start)) {
      nstart <- colnames(start)
      start <- as.vector(start[nrow(start), , drop = TRUE])
      names(start) <- nstart
    }
    nstart <- names(start)
    tns <- sapply(strsplit(nstart, ".", fixed = TRUE), function(x) { x[1] })
    nx <- names(x)
    for(id in nx) {
      if(!is.null(x[[id]]$smooth.construct)) {
        if(!is.null(x[[id]]$smooth.construct$model.matrix)) {
          if(length(take <- grep(paste(id, "p", sep = "."), nstart[tns %in% id], fixed = TRUE, value = TRUE))) {
            cn <- paste(id, "p", colnames(x[[id]]$smooth.construct$model.matrix$X), sep = ".")
            i <- grep2(take, cn, fixed = TRUE)
            if(length(i)) {
              tpar <- start[take[i]]
              i <- grep2(c(".edf", ".accepted", ".alpha"), names(tpar), fixed = TRUE)
              if(length(i))
                tpar <- tpar[-i]
              names(tpar) <- gsub(paste(id, "p.", sep = "."), "", names(tpar), fixed = TRUE)
              if(any(l <- grepl("tau2", take))) {
                tau2 <- start[take[l]]
                names(tau2) <- gsub(paste(id, "p.", sep = "."), "", names(tau2), fixed = TRUE)
                tpar <- c(tpar, tau2)
              }
              if(all(names(tpar) %in% names(x[[id]]$smooth.construct$model.matrix$state$parameters))) {
                x[[id]]$smooth.construct$model.matrix$state$parameters[names(tpar)] <- tpar
                x[[id]]$smooth.construct$model.matrix$state$fitted.values <- x[[id]]$smooth.construct$model.matrix$fit.fun(x[[id]]$smooth.construct$model.matrix$X, x[[id]]$smooth.construct$model.matrix$state$parameters)
              }
            }
          }
        }
        for(j in seq_along(x[[id]]$smooth.construct)) {
          tl <- x[[id]]$smooth.construct[[j]]$label
          tl <- paste(id, "s", tl, sep = ".")
          if(inherits(x[[id]]$smooth.construct[[j]], "nnet.boost")) {
            take <- tl
          } else {
            take <- grep(paste0(tl, "."), nstart[tns %in% id], fixed = TRUE, value = TRUE)
          }
          if(is.null(x[[id]]$smooth.construct[[j]]$by))
            x[[id]]$smooth.construct[[j]]$by <- "NA"
          if(x[[id]]$smooth.construct[[j]]$by == "NA") {
            take <- take[!grepl(paste(tl, ":", sep = ""), take, fixed = TRUE)]
          }
          if(length(take)) {
            tpar <- start[take]
            i <- grep2(c(".edf", ".accepted", ".alpha"), names(tpar), fixed = TRUE)
            tpar <- if(length(i)) tpar[-i] else tpar
            names(tpar) <- gsub(paste(tl, ".", sep = ""), "", names(tpar), fixed = TRUE)
            if(inherits(x[[id]]$smooth.construct[[j]], "nnet0.smooth")) {
              spar <- tpar
            } else {
              spar <- x[[id]]$smooth.construct[[j]]$state$parameters
              if(length(get.par(tpar, "b")))
                spar <- set.par(spar, get.par(tpar, "b"), "b")
              if(any(grepl("tau2", names(tpar)))) {
                spar <- set.par(spar, get.par(tpar, "tau2"), "tau2")
              }
            }
            x[[id]]$smooth.construct[[j]]$state$parameters <- spar
            x[[id]]$smooth.construct[[j]]$state$fitted.values <- x[[id]]$smooth.construct[[j]]$fit.fun(x[[id]]$smooth.construct[[j]]$X, x[[id]]$smooth.construct[[j]]$state$parameters)
          }
        }
      }
    }
  }
  
  return(x)
}


opt_lasso <- lasso <- function(x, y, start = NULL, adaptive = TRUE,
  lower = 0.001, upper = 1000,  nlambda = 100, lambda = NULL, multiple = FALSE,
  verbose = TRUE, digits = 4, flush = TRUE,
  nu = NULL, stop.nu = NULL, ridge = .Machine$double.eps^0.5,
  zeromodel = NULL, ...)
{
  method <- list(...)$method
  if(is.null(method))
    method <- 1
  
  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, update = bfit_iwls, ...)
  
  start2 <- start

  if(lower < 1e-20)
    lower <- 1e-20
  
  lambdas <- if(is.null(lambda)) {
    exp(seq(log(upper), log(lower), length = nlambda))
  } else lambda

  lambdas <- rep(list(lambdas), length = length(x))
  names(lambdas) <- names(x)

  lambdas <- as.matrix(do.call(if(multiple) "expand.grid" else "cbind", lambdas))
  
  if(length(verbose) < 2)
    verbose <- c(verbose, FALSE)
  
  ia <- if(flush) interactive() else FALSE
  
  par <- list(); ic <- NULL
  
  ptm <- proc.time()
  
  fuse <- NULL
  
  for(i in names(x)) {
    for(j in names(x[[i]]$smooth.construct)) {
      if(inherits(x[[i]]$smooth.construct[[j]], "lasso.smooth")) {
        x[[i]]$smooth.construct[[j]]$state$do.optim <- FALSE
        x[[i]]$smooth.construct[[j]]$fxsp <- TRUE
        fuse <- c(fuse, x[[i]]$smooth.construct[[j]]$fuse)
        if(adaptive) {
          tau2 <- get.par(x[[i]]$smooth.construct[[j]]$state$parameters, "tau2")
          tau2 <- rep(1/ridge, length.out = length(tau2))
          x[[i]]$smooth.construct[[j]]$state$parameters <- set.par(x[[i]]$smooth.construct[[j]]$state$parameters, tau2, "tau2")
          x[[i]]$smooth.construct[[j]]$LAPEN <- x[[i]]$smooth.construct[[j]]$S
          x[[i]]$smooth.construct[[j]]$S <- list(diag(length(get.par(x[[i]]$smooth.construct[[j]]$state$parameters, "b"))))
        }
      }
    }
  }
  
  fuse <- if(is.null(fuse)) FALSE else any(fuse)
  
  if(!is.null(nu))
    nu <- rep(nu, length.out = 2)
  if(!is.null(stop.nu))
    stop.nu <- rep(stop.nu, length.out = 2)
  
  if(adaptive & fuse) {
    if(verbose[1] & is.null(zeromodel))
      cat("Estimating adaptive weights\n---\n")
    if(is.null(zeromodel)) {
      if(method == 1) {
        zeromodel <- opt_bfit(x = x, y = y, start = start, verbose = verbose[1], nu = nu[2], stop.nu = stop.nu[2], ...)
      } else {
        zeromodel <- opt_optim(x = x, y = y, start = start, verbose = verbose[1], ...)
      }
    }
    x <- lasso_transform(x, zeromodel, nobs = nrow(y))
  } else {
    if(!is.null(zeromodel)) {
      x <- lasso_transform(x, zeromodel, nobs = nrow(y))
    }
  }
  
  for(l in 1:nrow(lambdas)) {
    if(l > 1)
      start <- unlist(par[[l - 1]])
    tau2 <- NULL
    for(i in names(x)) {
      for(j in names(x[[i]]$smooth.construct)) {
        if(inherits(x[[i]]$smooth.construct[[j]], "lasso.smooth")) {
          tau2 <- get.par(x[[i]]$smooth.construct[[j]]$state$parameters, "tau2")
          nt <- names(tau2)
          tau2 <- rep(1 / lambdas[l, i], length.out = length(tau2))
          names(tau2) <- paste(i, "s", x[[i]]$smooth.construct[[j]]$label, nt, sep = ".")
          if(!is.null(start) & (l > 1)) {
            if(all(names(tau2) %in% names(start))) {
              start[names(tau2)] <- tau2
            } else {
              start <- c(start, tau2)
            }
          } else {
            start <- c(start, tau2)
          }
        }
      }
    }
    
    if((l < 2) & !is.null(start2)) {
      start <- c(start, start2)
      start <- start[!duplicated(names(start))]
    }
    
    if(method == 1) {
      b <- opt_bfit(x = x, y = y, start = start, verbose = verbose[2], nu = nu[2], stop.nu = stop.nu[2], ...)
    } else {
      b <- opt_optim(x = x, y = y, start = start, verbose = verbose[2], ...)
    }
    
    nic <- grep("ic", names(b), value = TRUE, ignore.case = TRUE)
    if(!length(nic)) {
      b$edf <- sum(abs(unlist(b$parameters)) > .Machine$double.eps^0.25)
      b$BIC <- -2 * b$logLik + b$edf * log(nrow(y))
    }
    nic <- grep("ic", names(b), value = TRUE, ignore.case = TRUE)
    par[[l]] <- unlist(b$parameters)
    mstats <- c(b$logLik, b$logPost, b[[nic]], b[["edf"]])
    names(mstats) <- c("logLik", "logPost", nic, "edf")
    ic <- rbind(ic, mstats)
    
    if(!is.null(list(...)$track)) {
      plot(ic[, nic] ~ c(1:l), type = "l", xlab = "Iteration", ylab = nic)
    }
    
    if(!is.null(stop.nu)) {
      if(l > stop.nu)
        nu <- NULL
    }
    
    if(verbose[1]) {
      cat(if(ia) "\r" else if(l > 1) "\n" else NULL)
      vtxt <- paste(nic, " ", fmt(b[[nic]], width = 8, digits = digits),
                    " edf ", fmt(mstats["edf"], width = 6, digits = digits),
                    " lambda ", paste(fmt(if(!multiple) lambdas[l, 1] else lambdas[l, ], width = 6, digits = digits), collapse = ","),
                    " iteration ", formatC(l, width = nchar(nlambda)), sep = "")
      cat(vtxt)
      
      if(.Platform$OS.type != "unix" & ia) flush.console()
    }
  }
  
  elapsed <- c(proc.time() - ptm)[3]
  
  if(verbose[1]) {
    et <- if(elapsed > 60) {
      paste(formatC(format(round(elapsed / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(elapsed, 2), nsmall = 2), width = 5), "sec", sep = "")
    cat("\nelapsed time: ", et, "\n", sep = "")
  }

  colnames(lambdas) <- paste("lambda", names(x), sep = ".")
  ic <- cbind(ic, "lambda" = lambdas)
  rownames(ic) <- NULL
  attr(ic, "multiple") <- multiple
  class(ic) <- c("lasso.stats", "matrix")
  
  list("parameters" = do.call("rbind", par), "lasso.stats" = ic, "nobs" = nrow(y))
}

lasso_transform <- function(x, zeromodel, nobs = NULL, ...)
{
  if(bframe <- inherits(x, "bamlss.frame")) {
    if(is.null(x$x))
      stop("no 'x' object in 'bamlss.frame'!")
    x <- x$x
  }
  for(i in names(x)) {
    for(j in names(x[[i]]$smooth.construct)) {
      if(inherits(x[[i]]$smooth.construct[[j]], "lasso.smooth")) {
        if(!is.null(x[[i]]$smooth.construct[[j]]$LAPEN)) {
          x[[i]]$smooth.construct[[j]]$S <- x[[i]]$smooth.construct[[j]]$LAPEN
          x[[i]]$smooth.construct[[j]]$LAPEN <- NULL
        }
        if(x[[i]]$smooth.construct[[j]]$fuse) {
          if(is.list(zeromodel$parameters)) {
            beta <- get.par(zeromodel$parameters[[i]]$s[[j]], "b")
          } else {
            if(is.matrix(zeromodel$parameters)) {
              beta <- grep(paste(i, ".s.", j, ".", sep = ""), colnames(zeromodel$parameters), fixed = TRUE)
              beta <- get.par(zeromodel$parameters[nrow(zeromodel$parameters), beta], "b")
            } else {
              beta <- grep(paste(i, ".s.", j, ".", sep = ""), names(zeromodel$parameters), fixed = TRUE)
              beta <- get.par(zeromodel$parameters[beta], "b")
            }
          }
          df <- x[[i]]$smooth.construct[[j]]$lasso$df
          Af <- x[[i]]$smooth.construct[[j]]$Af
          w <- rep(0, ncol(Af))
          fuse_type <- x[[i]]$smooth.construct[[j]]$fuse_type
          k <- ncol(x[[i]]$smooth.construct[[j]]$X)
          if(is.null(nobs))
            nobs <- nrow(x[[i]]$smooth.construct[[j]]$X)
          if(x[[i]]$smooth.construct[[j]]$xt$gfx) {
            w <- NULL
            for(ff in 1:ncol(Af))
              w <- c(w, 1/abs(t(Af[, ff]) %*% beta))
          } else {
            nref <- nobs - sum(df)
            for(ff in 1:ncol(Af)) {
              ok <- which(Af[, ff] != 0)
              w[ff] <- if(fuse_type == "nominal") {
                if(length(ok) < 2) {
                  2 / (k + 1) * sqrt((df[ok[1]] + nref) / nobs)
                } else {
                  2 / (k + 1) * sqrt((df[ok[1]] + df[ok[2]]) / nobs)
                }
              } else {
                if(length(ok) < 2) {
                  sqrt((df[ok[1]] + nref) / nobs)
                } else {
                  sqrt((df[ok[1]] + df[ok[2]]) / nobs)
                }
              }
              w[ff] <- w[ff] * 1 / abs(t(Af[, ff]) %*% beta)
            }
          }
          names(w) <- paste("lasso", 1:length(w), sep = "")
          w[!is.finite(w)] <- 1e10
          x[[i]]$smooth.construct[[j]]$fixed.hyper <- w
        } else {
          w <- get.par(zeromodel$parameters[[i]]$s[[j]], "b")
          names(w) <- paste("lasso", 1:length(w), sep = "")
          w[!is.finite(w)] <- 1e10
          x[[i]]$smooth.construct[[j]]$fixed.hyper <- w
        }
      }
    }
  }
  
  if(bframe) {
    return(list("x" = x))
  } else {
    return(x)
  }
}

print.lasso.stats <- function(x, digits = 4, ...)
{
  ls <- attr(lasso_stop(x), "stats")
  ic <- grep("ic", names(ls), ignore.case = TRUE, value = TRUE)
  cat(ic, "=", ls[ic], "-> at lambda =", ls[grep("lambda", names(ls))], "\n")
  ls <- ls[!grepl("lambda", names(ls))]
  ls <- paste(names(ls), "=", round(ls, digits = digits), collapse = " ")
  cat(ls, "\n---\n")
  return(invisible(NULL))
}


lasso_coef <- function(x, ...) {
  cx <- coef.bamlss(x, ...)
  ncx <- if(!is.null(dim(cx))) colnames(cx) else names(cx)
  if(is.null(x$x))
    x$x <- smooth.construct(x)
  for(i in names(x$x)) {
    for(j in names(x$x[[i]]$smooth.construct)) {
      if(inherits(x$x[[i]]$smooth.construct[[j]], "lasso.smooth")) {
        for(jj in names(x$x[[i]]$smooth.construct[[j]]$lasso$trans)) {
          cid <- paste(i, ".s.", x$x[[i]]$smooth.construct[[j]]$label, ".",
                       x$x[[i]]$smooth.construct[[j]]$lasso$trans[[jj]]$colnames, sep = "")
          if(is.null(x$x[[i]]$smooth.construct[[j]]$lasso$trans[[jj]]$blockscale)) {
            if(is.null(dim(cx))) {
              cx[cid] <- cx[cid] / x$x[[i]]$smooth.construct[[j]]$lasso$trans[[jj]]$scale
            } else {
              cx[, cid] <- cx[, cid, drop = FALSE] / x$x[[i]]$smooth.construct[[j]]$lasso$trans[[jj]]$scale
            }
          } else {
            if(is.null(dim(cx))) {
              cx[cid] <- solve(x$x[[i]]$smooth.construct[[j]]$lasso$trans[[jj]]$blockscale, cx[cid])
            } else {
              for(ii in 1:nrow(cx)) {
                cx[ii, cid] <- solve(x$x[[i]]$smooth.construct[[j]]$lasso$trans[[jj]]$blockscale, cx[ii, cid])
              }
            }
          }
        }
      }
    }
  }
  cx
}


lasso_plot <- function(x, which = c("criterion", "parameters"), spar = TRUE, model = NULL, name = NULL,
  mstop = NULL, retrans = FALSE, color = NULL, show.lambda = TRUE, labels = NULL,
  digits = 2, ...)
{
  if(is.null(model))
    model <- x$family$names
  model <- x$family$names[pmatch(model, x$family$names)]
  if(any(is.na(model))) {
    model <- model[!is.na(model)]
    if(!length(model))
      stop("argument model is spcified wrong")
    else
      warning("argument model is spcified wrong")
  }
  if(!is.character(which)) {
    which <- c("criterion", "parameters")[as.integer(which)]
  } else {
    which <- tolower(which)
    which <- match.arg(which, several.ok = TRUE)
  }
  if(is.null(mstop))
    mstop <- 1:nrow(x$parameters)
  if(retrans)
    x$parameters <- lasso_coef(x)
  npar <- colnames(x$parameters)
  for(j in c("Intercept", ".edf", ".lambda", ".tau"))
    npar <- npar[!grepl(j, npar, fixed = TRUE)]
  x$parameters <- x$parameters[, npar, drop = FALSE]
  ic <- x$model.stats$optimizer$lasso.stats
  multiple <- attr(ic, "multiple")
  log_lambda <- log(ic[, grep("lambda", colnames(ic)), drop = FALSE])
  nic <- grep("ic", colnames(ic), value = TRUE, ignore.case = TRUE)
  if(spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    n <- if("criterion" %in% which) {
      if(multiple) length(model) else 1
    } else 0
    if("parameters" %in% which)
      n <- n + length(model)
    par(mfrow = n2mfrow(n), mar = c(5.1, 5.1, 4.1, 1.1))
  }
  at <- pretty(1:nrow(ic))
  at[at == 0] <- 1
  
  if("criterion" %in% which) {
    if(!multiple) {
      plot(ic[, nic], type = "l",
        xlab = expression(log(lambda[, 1])), ylab = nic, axes = FALSE, lwd = list(...)$lwd)
      at <- pretty(mstop)
      at[at == 0] <- 1
      axis(1, at = at, labels = as.numeric(fmt(log_lambda[, 1][mstop][at], digits)))
      axis(2)
      if(show.lambda) {
        i <- which.min(ic[, nic])
        abline(v = i, col = "lightgray", lwd = 2, lty = 2)
        val <- round(ic[i, grep("lambda", colnames(ic))[1]], 4)
        axis(3, at = i, labels = substitute(paste(lambda, '=', val)))
      }
      if(!is.null(main <- list(...)$main))
        mtext(main, side = 3, line = 2.5, cex = 1.2, font = 2)
      box()
    } else {
      main <- list(...)$main
      if(is.null(main))
        main <- model
      main <- rep(main, length.out = length(model))
      k <- 1
      for(m in model) {
        imin <- which.min(ic[, nic])
        lambda_min <- ic[imin, grep("lambda", colnames(ic))]
        tlambda <- names(lambda_min)
        tlambda <- tlambda[!grepl(m, tlambda)]
        take <- NULL
        for(j in tlambda)
          take <- cbind(take, ic[, j] == lambda_min[j])
        take <- apply(take, 1, all)
        tic <- ic[take, nic]
        plot(tic, type = "l", xlab = expression(log(lambda[, 1])), ylab = nic, axes = FALSE, lwd = list(...)$lwd)
        at <- pretty(1:length(tic))
        at[at == 0] <- 1
        axis(1, at = at, labels = as.numeric(fmt(log_lambda[take, paste("lambda", m, sep = ".")][at], digits)))
        axis(2)
        if(show.lambda) {
          i <- which.min(tic)
          abline(v = i, col = "lightgray", lwd = 2, lty = 2)
          val <- lambda_min[paste("lambda", m, sep = ".")]
          axis(3, at = i, labels = substitute(paste(lambda, '=', val)))
        }
        box()
        if(!is.expression(main[k])) {
          if(main[k] != "")
            mtext(main[k], side = 3, line = 2.5, cex = 1.2, font = 2)
        } else {
          mtext(main[k], side = 3, line = 2.5, cex = 1.2, font = 2)
        }
        k <- k + 1
      }
    }
  }
  
  if("parameters" %in% which) {
    main <- list(...)$main
    if(is.null(main))
      main <- model
    main <- rep(main, length.out = length(model))
    imin <- which.min(ic[, nic])
    lambda_min <- ic[imin, grep("lambda", colnames(ic))]
    k <- 1
    for(m in model) {
      if(spar)
        par(mar = c(5.1, 5.1, 4.1, 10.1))

      tpar <- x$parameters[, grep(paste(m, ".", sep = ""), colnames(x$parameters), fixed = TRUE), drop = FALSE]

      if(multiple) {
        tlambda <- names(lambda_min)
        tlambda <- tlambda[!grepl(m, tlambda)]
        take <- NULL
        for(j in tlambda)
          take <- cbind(take, ic[, j] == lambda_min[j])
        take <- apply(take, 1, all)
        tpar <- tpar[take, , drop = FALSE]
      } else {
        take <- 1:nrow(tpar)
      }
      if(!is.null(name))
        tpar <- tpar[, grep2(name, colnames(tpar), fixed = TRUE), drop = FALSE]
      if(max(mstop) > nrow(tpar))
        mstop <- nrow(tpar)
      tpar <- tpar[if(length(mstop) < 2) 1:mstop else mstop, , drop = FALSE]
      xn <- sapply(strsplit(colnames(tpar), ".", fixed = TRUE), function(x) { x[1] })
      if(length(unique(xn)) < 2)
        xn <- sapply(strsplit(colnames(tpar), ".", fixed = TRUE), function(x) { x[3] })
    
      cols <- if(is.null(color)) {
        if(length(unique(xn)) < 2) "black" else rainbow_hcl(length(unique(xn)))
      } else {
        if(is.function(color)) {
          color(length(unique(xn)))
        } else {
          rep(color, length.out = length(unique(xn)))
        }
      }
      add <- if(is.null(list(...)$add)) FALSE else list(...)$add
      nolabels <- if(is.null(list(...)$nolabels)) FALSE else list(...)$nolabels
      matplot(tpar, type = "l", lty = 1, col = cols[as.factor(xn)],
        xlab = expression(log(lambda)), ylab = expression(beta[j]), axes = FALSE, add = add,
        lwd = list(...)$lwd)
      if(!nolabels) {
        if(is.null(labels)) {
          labs <- labs0 <- colnames(tpar)
          plab <- tpar[nrow(tpar), ]
          o <- order(plab, decreasing = TRUE)
          labs <- labs[o]
          plab <- plab[o]
          rplab <- diff(range(plab))
          for(i in 1:(length(plab) - 1)) {
            dp <- abs(plab[i] - plab[i + 1]) / rplab
            if(dp <= 0.02) {
              labs[i + 1] <- paste(c(labs[i], labs[i + 1]), collapse = ",")
              labs[i] <- ""
            }
          }
          labs <- labs[order(o)]
          if(!is.null(name)) {
            for(j in seq_along(name))
              labs <- gsub(name[j], "", labs, fixed = TRUE)
          }
        } else labs <- rep(labels, length.out = ncol(tpar))
        at <- tpar[nrow(tpar), ]
        at <- at[labs != ""]
        labs <- labs[labs != ""]
        axis(4, at = at, labels = labs, las = 1, cex.axis = list(...)$labcex)
      }
      at <- pretty(1:nrow(tpar))
      at[at == 0] <- 1
      axis(1, at = at, labels = as.numeric(fmt(log_lambda[take, paste("lambda", m, sep = ".")][at], digits)))
      axis(2)
      if(show.lambda) {
        i <- which.min(ic[take, nic])
        abline(v = i, col = "lightgray", lwd = 1, lty = 2)
        cat(colnames(tpar)[abs(tpar[i, ]) > 0.009], "\n")
        val <- round(lambda_min[paste("lambda", m, sep = ".")], digits)
        if(multiple) {
          lval <- parse(text = paste('paste(lambda[', m, '], "=", ', val, ')', sep = ''))
          axis(3, at = i, labels = lval)
        } else {
          axis(3, at = i, labels = substitute(paste(lambda, '=', val)))
        }
      }
      box()
      if(!is.expression(main[k])) {
        if(main[k] != "")
          mtext(main[k], side = 3, line = 2.5, cex = 1.2, font = 2)
      } else {
        mtext(main[k], side = 3, line = 2.5, cex = 1.2, font = 2)
      }
      k <- k + 1
    }
  }
  
  return(invisible(NULL))
}


lasso_stop <- function(x)
{
  if(!inherits(x, "lasso.stats"))
    x <- x$model.stats$optimizer$lasso.stats
  nic <- grep("ic", colnames(x), value = TRUE, ignore.case = TRUE)
  i <- which.min(x[, nic])
  attr(i, "stats") <- x[i, ]
  i
}


## Deep learning bamlss.
dl.bamlss <- function(object,
  optimizer = "adam", epochs = 30, batch_size = NULL,
  nlayers = 2, units = 100, activation = "sigmoid", l1 = NULL, l2 = NULL,
  verbose = TRUE, ...)
{
  stopifnot(requireNamespace("keras"))
  stopifnot(requireNamespace("tensorflow"))

  if(!inherits(object, "bamlss")) {
    object <- bamlss.frame(object, ...)
  }

  y <- object$y
  nobs <- nrow(y)
  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }

  nx <- names(object$formula)
  X <- list()
  for(i in nx) {
    X[[i]] <- object$x[[i]]$model.matrix
    if(!is.null(object$x[[i]]$smooth.construct)) {
      for(j in seq_along(object$x[[i]]$smooth.construct))
        X[[i]] <- cbind(X[[i]], object$x[[i]]$smooth.construct[[j]]$X)
    }
  }

  family <- family(object)

  if(is.null(family$keras$nloglik))
    stop("no keras negative loglik() function is specified!")

  nll <- family$keras$nloglik

  if(is.null(optimizer))
    optimizer <- keras::optimizer_rmsprop(lr = 0.0001)
  if(is.character(optimizer)) {
      optimizer <- match.arg(optimizer, c("adam", "sgd", "rmsprop",
        "adagrad", "adadelta", "adamax", "nadam"))
  }

  models <- inputs <- outputs <- list()

  units <- rep(units, length.out = nlayers)
  activation <- rep(activation, length.out = nlayers)

  if(!is.null(l1))
    l1 <- rep(l1, length.out = nlayers)
  if(!is.null(l2))
    l2 <- rep(l2, length.out = nlayers)

  pen <- NULL
  if(!is.null(l1) & is.null(l2))
    pen <- paste0('regularizer_l1(', l1, ')')
  if(is.null(l1) & !is.null(l2))
    pen <- paste0('regularizer_l2(', l2, ')')
  if(!is.null(l1) & !is.null(l2))
    pen <- paste0('regularizer_l1_l2(', l1,, ', ',  l2, ')')

  for(j in 1:length(X)) {
    inputs[[j]] <- keras::layer_input(shape = ncol(X[[j]]))

    itcpt <- FALSE
    if(ncol(X[[j]]) < 2) {
      if(colnames(X[[i]]) == "(Intercept)")
        itcpt <- TRUE
    }

    if(itcpt) {
      opts <- c('outputs[[j]] <- inputs[[j]] %>%',
        'layer_dense(units = 1, use_bias = FALSE)')
    } else {
      opts <- c('outputs[[j]] <- inputs[[j]] %>%',
        paste0('layer_dense(units = ', units,
          if(!is.null(pen)) paste0(', kernel_regularizer = ', pen) else NULL,
          ', activation = "', activation, '") %>%'),
        'layer_dense(units = 1)')
    }

    opts <- paste(opts, collapse = " ")
    eval(parse(text = opts))
  }

  final_output <- keras::layer_concatenate(outputs)

  model <- keras::keras_model(inputs, final_output)

  model <- keras::compile(model,
    loss = nll, 
    optimizer = optimizer
  )

  names(X) <- NULL
  if(length(X) > 1) {
    Y <- matrix(y, ncol = 1)
    for(j in 1:(length(nx) - 1))
      Y <- cbind(Y, 1)
  } else {
    Y <- matrix(y, ncol = 1)
  }
  if(length(X) < 2)
    X <- X[[1L]]

  ptm <- proc.time()

  history <- keras::fit(model,
    x = X, 
    y = Y, 
    epochs = epochs, batch_size = batch_size,
    verbose = as.integer(verbose)
  )

  elapsed <- c(proc.time() - ptm)[3]
  
  if(verbose) {
    cat("\n")
    et <- if(elapsed > 60) {
      paste(formatC(format(round(elapsed / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(elapsed, 2), nsmall = 2), width = 5), "sec", sep = "")
    cat("\n elapsed time: ", et, "\n", sep = "")
  }

  object$dnn <- model
  object$fitted.values <- as.data.frame(predict(model, X))
  colnames(object$fitted.values) <- nx
  object$elapsed <- elapsed
  object$history <- history

  class(object) <- "dl.bamlss"

  return(object)
}


## Extractor functions.
fitted.dl.bamlss <- function(object, ...) { object$fitted.values }
family.dl.bamlss <- function(object) { object$family }
residuals.dl.bamlss <- function(object, ...) { residuals.bamlss(object, ...) }
plot.dl.bamlss <- function(x, ...) { plot(x$history, ...) }


## Predict function.
predict.dl.bamlss <- function(object, newdata, model = NULL,
  type = c("link", "parameter"), drop = TRUE, ...)
{
  ## If data have been scaled (scale.d = TRUE)
  if(!missing(newdata) & ! is.null(attr(object$model.frame, 'scale')) ) {
    sc <- attr(object$model.frame, 'scale')
    for( name in unique(unlist(lapply(sc,names))) ) {
      newdata[,name] <- (newdata[,name] - sc$center[name] ) / sc$scale[name]
    }
  }
  if(missing(newdata))
    newdata <- NULL
  if(is.null(newdata)) {
    newdata <- model.frame.bamlss.frame(object)
  } else {
    if(is.character(newdata)) {
      if(file.exists(newdata <- path.expand(newdata)))
        newdata <- read.table(newdata, header = TRUE, ...)
      else stop("cannot find newdata")
    }
    if(is.matrix(newdata) || is.list(newdata))
      newdata <- as.data.frame(newdata)
  }

  type <- match.arg(type)

  nx <- names(object$formula)
  X <- list()
  for(i in seq_along(nx)) {
    tfi <- drop.terms.bamlss(object$x[[i]]$terms,
      sterms = FALSE, keep.response = FALSE, data = newdata, specials = NULL)
    X[[i]] <- model.matrix(tfi, data = newdata)
    if(!is.null(object$x[[i]]$smooth.construct)) {
      for(j in seq_along(object$x[[i]]$smooth.construct)) {
        X[[i]] <- cbind(X[[i]], PredictMat(object$x[[i]]$smooth.construct[[j]], newdata))
      }
    }
  }

  if(length(X) < 2)
    X <- X[[1L]]

  pred <- as.data.frame(predict(object$dnn, X))
  colnames(pred) <- nx

  if(!is.null(model)) {
    if(is.character(model))
      nx <- nx[grep(model, nx)[1]]
    else
      nx <- nx[as.integer(model)]
  }

  pred <- pred[nx]

  if(type == "parameter") {
    links <- object$family$links[nx]
    for(j in seq_along(links)) {
      if(links[j] != "identity") {
        linkinv <- make.link2(links[j])$linkinv
        pred[[j]] <- linkinv(pred[[j]])
      }
    }
  }

  if((length(pred) < 2) & drop)
    pred <- pred[[1L]]

  return(pred)
}


## Most likely transformations.
opt_mlt <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
  criterion = c("AICc", "BIC", "AIC"),
  eps = .Machine$double.eps^0.25, maxit = 400,
  verbose = TRUE, digits = 4, flush = TRUE, nu = NULL, stop.nu = NULL, ...)
{
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("design construct names mismatch with family names!")
  
  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, ...)

  nobs <- nrow(y)
  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }

  Fy <- ecdf(y)
  Fy <- Fy(y)
  Fy[Fy > 0.9999] <- 0.999
  Fy[Fy < 0.0001] <- 0.001
  Yhat <- family$distr$q(Fy)

  opt <- opt_bfit(x = x, y = data.frame("y" = Yhat),
    family = complete.bamlss.family(Gaussian_bamlss()),
    eps = eps, maxit = maxit, nu = nu, update = bfit_optim())
  
  criterion <- match.arg(criterion)
  np <- length(nx)
  
  if(!is.null(nu)) {
    if(nu < 0)
      nu <- NULL
  }
  
  if(!is.null(start))
    x <- set.starting.values(x, start)
  else
    x <- set.starting.values(x, opt$parameters)
  eta <- get.eta(x)
  
  if(!is.null(weights))
    weights <- as.data.frame(weights)
  if(!is.null(offset)) {
    offset <- as.data.frame(offset)
    for(j in nx) {
      if(!is.null(offset[[j]]))
        eta[[j]] <- eta[[j]] + offset[[j]]
    }
  } else {
    if(is.null(start))
      eta <- init.eta(eta, y, family, nobs)
  }
  
  ia <- if(flush) interactive() else FALSE

  eps0 <- eps + 1; iter <- 0
  edf <- get.edf(x, type = 2)
  ptm <- proc.time()
  while(eps0 > eps & iter < maxit) {
    eta0 <- eta
    ## Cycle through all terms.
    for(sj in seq_along(x$mu$smooth.construct)) {
      ## Get updated parameters.
      p.state <- mlt_update(x$mu$smooth.construct[[sj]],
        family$distr, y, eta, edf = edf, weights = weights$mu,
        iteration = iter, criterion = criterion)
    }
  }

  stop("here!")
}

mlt_update <- function(x, distr, y, eta, edf, weights, iteration, criterion)
{
  beta <- get.par(x$state$parameters, "b")
  score <- t(x$X) %*% (distr$dd(eta$mu) / distr$d(eta$mu))
  if(inherits(x, "mlt.smooth"))
    score <- score + t(x$dX) %*% as.numeric((1 / (x$dX %*% beta + 1e-10)))
  w <- distr$ddd(eta$mu) / distr$d(eta$mu) - (distr$dd(eta$mu) / distr$d(eta$mu))^2
  hess <- crossprod(x$X * w, x$X)
  if(inherits(x, "mlt.smooth"))
    hess <- hess - crossprod(x$dX * as.numeric(1 / (x$dX %*% beta + 1e-10)^2), x$dX)

print(beta)

  beta <- beta + hess %*% score

print(beta)

  stop("yess!\n")
}


boost.net <- function(formula, maxit = 1000, nu = 1, nodes = 10, df = 4,
  lambda = NULL, dropout = NULL, flush = TRUE, initialize = TRUE,
  eps = .Machine$double.eps^0.25, verbose = TRUE, digits = 4,
  activation = "sigmoid", 
  r = list("sigmoid" = 0.01, "gauss" = 0.01, "sin" = 0.01, "cos" = 0.01),
  s = list("sigmoid" = 10000, "gauss" = 20, "sin" = 20, "cos" = 20),
  select = FALSE, ...)
{
  bf <- bamlss.frame(formula, ...)
  y <- bf$y
  has_offset <- any(grepl("(offset)", names(bf$model.frame), fixed = TRUE))

  nx <- names(bf$x)
  np <- length(nx)
  nobs <- nrow(y)
  nu <- rep(nu, length.out = np)
  names(nu) <- nx

  nodes <- rep(nodes, length.out = np)
  names(nodes) <- nx

  if(!is.null(lambda)) {
    lambda <- rep(lambda, length.out = np)
    names(lambda) <- nx
  }

  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }

  if(!is.null(dropout)) {
    dropout <- rep(dropout, length.out = np)
    if(any(dropout > 1) | any(dropout < 0))
      stop("argument dropout must be between [0,1]!")
    names(dropout) <- nx
  }

  ntake <- rep(NA, length.out = np)
  names(ntake) <- nx
  Xn <- s01 <- list()
  beta <- taken <- list()
  for(i in nx) {
    k <- ncol(bf$x[[i]]$model.matrix)
    if(!is.null(dropout))
      ntake[i] <- ceiling(k * (1 - dropout[i]))
    else
      ntake[i] <- k - 1L
    if(k > 1) {
      w <- list()
      for(j in activation)
        w[[j]] <- n.weights(nodes[i], k = ntake[i], type = j)
      w <- unlist(w)
      if(!is.null(dropout))
        taken[[i]] <- matrix(0L, nrow = maxit, ncol = ntake[i])
      beta[[i]] <- matrix(0, nrow = maxit, ncol = k + (nodes[i] * length(activation)) + length(w))
      colnames(beta[[i]]) <- c(colnames(bf$x[[i]]$model.matrix),
        paste0("b", 1:(nodes[i] * length(activation))), names(w))
      Xn[[i]] <- bf$x[[i]]$model.matrix
      for(j in 2:k) {
        xmin <- min(Xn[[i]][, j], 2, na.rm = TRUE)
        xmax <- max(Xn[[i]][, j], 2, na.rm = TRUE)
        if((xmax - xmin) < sqrt(.Machine$double.eps)) {
          xmin <- 0
          xmax <- 1
        }
        Xn[[i]][, j] <- (Xn[[i]][, j] - xmin) / (xmax - xmin)
        s01[[i]]$xmin <- c(s01[[i]]$xmin, xmin)
        s01[[i]]$xmax <- c(s01[[i]]$xmax, xmax)
      }
    } else {
      beta[[i]] <- matrix(0, nrow = maxit, ncol = 1)
      colnames(beta[[i]]) <- "(Intercept)"
      nodes[i] <- -1
    }
  }

  if(initialize & !has_offset) {
    objfun <- function(par) {
      eta <- list()
      for(i in seq_along(nx))
        eta[[nx[i]]] <- rep(par[i], length = nobs)
      ll <- bf$family$loglik(y, bf$family$map2par(eta))
      return(ll)
    }
    
    gradfun <- function(par) {
      eta <- list()
      for(i in seq_along(nx))
        eta[[nx[i]]] <- rep(par[i], length = nobs)
      peta <- bf$family$map2par(eta)
      grad <- par
      for(j in nx) {
        score <- process.derivs(bf$family$score[[j]](y, peta, id = j), is.weight = FALSE)
        grad[i] <- mean(score)
      }
      return(grad)
    }

    start <- init.eta(get.eta(bf$x), y, bf$family, nobs)

    start <- unlist(lapply(start, mean, na.rm = TRUE))    
    opt <- optim(start, fn = objfun, gr = gradfun, method = "BFGS", control = list(fnscale = -1))

    eta <- list()
    for(i in nx) {
      beta[[i]][1, "(Intercept)"] <- as.numeric(opt$par[i])
      eta[[i]] <- rep(as.numeric(opt$par[i]), length = nobs)
    }
  } else {
    eta <- get.eta(bf$x)
  }

  if(has_offset) {
    for(i in nx) {
      eta[[i]] <- eta[[i]] + bf$model.frame[["(offset)"]][, i]
    }
  }

  logLik <- rep(0, maxit)
  logLik[1] <- bf$family$loglik(y, bf$family$map2par(eta))
  ia <- if(flush) interactive() else FALSE
  iter <- 2
  eps0 <- eps + 1
  ll_contrib <- rep(NA, np)
  names(ll_contrib) <- nx
  ll_contrib_save <- list()
  for(i in nx)
    ll_contrib_save[[i]] <- rep(0, maxit)
  par <- bpar <- Z <- list()
  tau2o <- rep(0.1, np)
  names(tau2o) <- nx
  edfn <- rep(NA, np)
  names(edfn) <- nx

  ptm <- proc.time()

  while(iter <= maxit & eps0 > eps) {
    eta0 <- eta
    ll0 <- bf$family$loglik(y, bf$family$map2par(eta))
    for(i in nx) {
      peta <- bf$family$map2par(eta)
      grad <- process.derivs(bf$family$score[[i]](y, peta, id = i), is.weight = FALSE)
      if(nodes[i] > 0) {
        Z[[i]] <- NULL
        w <- list()
        k <- ncol(bf$x[[i]]$model.matrix)
        for(j in activation) {
          if(is.null(dropout)) {
            w[[j]] <- n.weights(nodes[i], k = ntake[i],
              rint = r[[j]], sint = s[[j]], type = j,
              x = Xn[[i]][sample(1:nobs, size = nodes[i], replace = FALSE), -1, drop = FALSE])
            Z[[i]] <- cbind(Z[[i]], nnet2Zmat(Xn[[i]], w[[j]], j))
          } else {
            taken[[i]][iter, ] <- sample(2:k, size = ntake[i], replace = FALSE)
            w[[j]] <- n.weights(nodes[i], k = ntake[i],
              rint = r[[j]], sint = s[[j]], type = j,
              x = Xn[[i]][sample(1:nobs, size = nodes[i], replace = FALSE), taken[[i]][iter, ], drop = FALSE])
            Z[[i]] <- cbind(Z[[i]], nnet2Zmat(Xn[[i]][, c(1, taken[[i]][iter, ]), drop = FALSE], w[[j]], j))
          }
        }
        Z[[i]] <- cbind(bf$x[[i]]$model.matrix, Z[[i]])
        S <- diag(c(rep(0, k), rep(1, ncol(Z[[i]]) - k)))
        ZZ <- crossprod(Z[[i]])

        if(is.null(lambda)) {
          fn <- function(tau2) {
            Si <- 1 / tau2 * S
            P <- matrix_inv(ZZ + Si)
            b <- drop(P %*% crossprod(Z[[i]], grad))
            fit <- Z[[i]] %*% b
            edf <- sum_diag(ZZ %*% P) - k
            ic <- if(is.null(df)) {
              sum((grad - fit)^2) + 2 * edf
            } else {
              (df - edf)^2
            }
            return(ic)
          }

          tau2o[i] <- tau2.optim(fn, tau2o[i], maxit = 1e+04, force.stop = FALSE)
        } else {
          tau2o[i] <- 1/lambda[i]
        }

        S <- 1 / tau2o[i] * S
        P <- matrix_inv(ZZ + S)
        b <- nu[i] * drop(P %*% crossprod(Z[[i]], grad))
        par[[i]] <- c(b, unlist(w))
        bpar[[i]] <- b
        eta[[i]] <- eta[[i]] + Z[[i]] %*% b
        edfn[i] <- sum_diag(ZZ %*% P) - k          
      } else {
        mgrad <- nu[i] * mean(grad)
        eta[[i]] <- eta[[i]] + mgrad
        par[[i]] <- bpar[[i]] <- mgrad
        Z[[i]] <- matrix(1, nrow = length(grad), ncol = 1)
      }
      ll1 <- bf$family$loglik(y, bf$family$map2par(eta))
      if(ll1 < ll0) {
        nu[i] <- nu[i] * 0.9
        next
      }
      ll_contrib[i] <- ll1 - ll0
      if(select) {
        eta[[i]] <- eta0[[i]]
      } else {
        ll_contrib_save[[i]][iter] <- ll_contrib[i]
        beta[[i]][iter, ] <- par[[i]]
      }
    }

    if(select) {
      i <- nx[which.max(ll_contrib)]
      beta[[i]][iter, ] <- par[[i]]
      eta[[i]] <- eta[[i]] + Z[[i]] %*% bpar[[i]]
      ll_contrib_save[[i]][iter] <- ll_contrib[i]
    }

    eps0 <- do.call("cbind", eta)
    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

    ll <- bf$family$loglik(y, bf$family$map2par(eta))
    logLik[iter] <- ll

    iter <- iter + 1
    if(verbose) {
      cat(if(ia) "\r" else if(iter > 1) "\n" else NULL)
      vtxt <- paste(
        "logLik ", fmt(ll, width = 8, digits = digits),
        " edf ", paste(paste(nx, fmt(edfn, digits = 2, width = 4)), collapse = " "),
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter - 1L, width = nchar(maxit)), sep = "")
      cat(vtxt)
        
      if(.Platform$OS.type != "unix" & ia) flush.console()
    }
  }

  elapsed <- c(proc.time() - ptm)[3]

  if(verbose) {
    cat("\n")
    et <- if(elapsed > 60) {
      paste(formatC(format(round(elapsed / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(elapsed, 2), nsmall = 2), width = 5), "sec", sep = "")
    cat("elapsed time: ", et, "\n", sep = "")
  }

  scale <- list()
  for(i in nx) {
    beta[[i]] <- beta[[i]][1L:(iter - 1L), , drop = FALSE]
    ll_contrib_save[[i]]<- cumsum(ll_contrib_save[[i]][1:(iter - 1L)])
    scale[[i]] <- attr(bf$x[[i]]$model.matrix, "scale")
    if(!is.null(dropout)) {
      taken[[i]] <- taken[[i]][1L:(iter - 1L), , drop = FALSE]
      taken[[i]][1L, ] <- taken[[i]][2L, ]
    }
  }

  rval <- list(
    "parameters" = beta,
    "fitted.values" = eta,
    "loglik" = data.frame("loglik" = logLik[1L:(iter - 1L)]),
    "family" = bf$family,
    "formula" = bf$formula,
    "nodes" = nodes,
    "elapsed" = elapsed,
    "activation" = activation,
    "scale" = scale,
    "s01" = s01,
    "taken" = taken,
    "ntake" = ntake,
    "dropout" = dropout
  )
  rval$loglik[["contrib"]] <- do.call("cbind", ll_contrib_save)
  rval$call <- match.call()

  class(rval) <- "boost.net"

  return(rval)
}


predict.boost.net <- function(object, newdata, model = NULL,
  verbose = FALSE, cores = 1, mstop = NULL, matrix = FALSE, ...)
{
  nx <- object$family$names
  formula <- object$formula
  for(i in nx) {
    formula[[i]]$formula <- delete.response(formula[[i]]$formula)
    formula[[i]]$fake.formula <- delete.response(formula[[i]]$fake.formula)
  }
  bf <- bamlss.frame(formula, data = newdata, family = object$family)
  Xn <- list()
  for(i in nx) {
    if(!is.null(object$scale[[i]])) {
      for(j in 1:ncol(bf$x[[i]]$model.matrix)) {
        bf$x[[i]]$model.matrix[, j] <- (bf$x[[i]]$model.matrix[, j] - object$scale[[i]]$center[j]) / object$scale[[i]]$scale[j]
      }
    }
    if(!is.null(object$s01[[i]])) {
      Xn[[i]] <- bf$x[[i]]$model.matrix
      for(j in 1:length(object$s01[[i]]$xmin)) {
        Xn[[i]][, j + 1L] <- (Xn[[i]][, j + 1L] - object$s01[[i]]$xmin[j]) / (object$s01[[i]]$xmax[j] - object$s01[[i]]$xmin[j])
      }
    }
  }
  activation <- object$activation
  nodes <- object$nodes
  if(is.null(model))
    model <- nx
  for(j in seq_along(model))
    model[j] <- grep(model[j], nx, fixed = TRUE, value = TRUE)
  fit <- list()
  for(j in model) {
    fit[[j]] <- 0.0
    k <- ncol(bf$x[[j]]$model.matrix)
    if(!is.null(object$dropout))
      ind <- as.factor(sort(rep(rep(1:nodes[j]), object$ntake[i] + 1L)))
    else
      ind <- as.factor(sort(rep(rep(1:nodes[j]), k)))
    nr <- nrow(object$parameters[[j]])
    if(!is.null(mstop))
      nr <- min(c(nr, mstop))
    if(cores < 2) {
      for(i in 1:nr) {
        if(verbose)
          cat(i, "/", sep = "")
        if(nodes[j] > 0) {
          b <- object$parameters[[j]][i, 1:(k + nodes[j] * length(activation))]
          w <- object$parameters[[j]][i, -c(1:(k + nodes[j] * length(activation)))]
          Z <- NULL
          for(a in activation) {
            wa <- split(w[grep(a, names(w))], ind)
            if(is.null(object$dropout)) {
              Z <- cbind(Z, nnet2Zmat(Xn[[j]], wa, a))
            } else {
              Z <- cbind(Z, nnet2Zmat(Xn[[j]][, c(1, object$taken[[j]][i, ]), drop = FALSE], wa, a))
            }
          }
          Z <- cbind(bf$x[[j]]$model.matrix, Z)
          fit[[j]] <- fit[[j]] + drop(Z %*% b)
        } else {
          fit[[j]] <- fit[[j]] + object$parameters[[j]][i, "(Intercept)"]
        }
      }
    } else {
      jind <- split(1:nr, as.factor(sort(rep(1:cores, length.out = nr))))
      parallel_fun <- function(cid) {
        if(verbose)
          cat(j, ": started core", cid, "\n", sep = "")
        fit2 <- 0
        for(i in jind[[cid]]) {
          if(nodes[j] > 0) {
            b <- object$parameters[[j]][i, 1:(k + nodes[j] * length(activation))]
            w <- object$parameters[[j]][i, -c(1:(k + nodes[j] * length(activation)))]
            Z <- NULL
            for(a in activation) {
              wa <- split(w[grep(a, names(w))], ind)
              if(is.null(object$dropout)) {
                Z <- cbind(Z, nnet2Zmat(Xn[[j]], wa, a))
              } else {
                Z <- cbind(Z, nnet2Zmat(Xn[[j]][, c(1, object$taken[[j]][i, ]), drop = FALSE], wa, a))
              }
            }
            Z <- cbind(bf$x[[j]]$model.matrix, Z)
            fit2 <- fit2 + drop(Z %*% b)
          } else {
            fit2 <- fit2 + object$parameters[[j]][i, "(Intercept)"]
          }
        }
        if(verbose)
          cat(j, ": finished core", cid, "\n", sep = "")
        fit2
      }
      fit[[j]] <- parallel::mclapply(1:cores, parallel_fun, mc.cores = cores)
      fit[[j]] <- do.call("cbind", fit[[j]])
      fit[[j]] <- rowSums(fit[[j]])
    }
  }
  if(verbose)
    cat("\n")
  if(length(fit) < 2) {
    fit <- fit[[1L]]
  } else {
    fit <- as.data.frame(fit)
  }
  return(fit)
}


################################################################################
####                     STOCHASTIC GRADIENT DESCENT                        ####
################################################################################

####  ## sgd fitter
####  sgdfit <- function(x, y, gammaFun = function(i) 1/i, shuffle = TRUE,
####                     CFun = function(beta) diag(length(beta)),
####                     start = rep(0, ncol(x)), i.state = 0, link = function(x) x) {
####  
####      N <- length(y)
####      
####      ## shuffle observations
####      shuffle <- if(shuffle) sample(1L:N) else 1L:N
####      
####      ## Explicit SVG
####      beta     <- start
####      betaXVec <- matrix(0, nrow = N, ncol = length(beta))
####      
####      for (i in seq.int(N)) {
####         mu   <- drop(link(beta %*% x[shuffle[i],]))
####         grad <- (y[shuffle[i]] - mu)
####         beta <- beta + gammaFun(i + i.state) * grad * drop(x[shuffle[i],] %*% CFun(beta))
####         betaXVec[i,] <- beta
####      }
####  
####      rval <- list()
####      rval$shuffle <- shuffle
####      rval$coef    <- beta
####      rval$y       <- y
####      rval$x       <- x
####      rval$i.state <- i
####      rval$diagnostics <- list("betaMat" = betaXVec)
####      class(rval) <- "sgdfit"
####  
####      rval
####  }
####  

## Implicit SGD
opt_isgd <- function(x, y, family, weights = NULL, offset = NULL,
                 gammaFun = function(i) 1/(1+i), shuffle = TRUE,
                 CFun = function(beta) diag(length(beta)),
                 start = NULL, i.state = 0) {

    ## constants
    nx <- family$names
    if(!all(nx %in% names(x)))
        stop("parameter names mismatch with family names!")

    N  <- nrow(y)
    y  <- as.matrix(y[[1]])

    ## shuffle observations
    shuffle <- if(shuffle) sample(1L:N) else 1L:N
   
    ## grep design matrices
    X <- sgd_grep_X(x)
    m <- sapply(X, ncol)      ## number of columns in each design matrix
    
    ## make a list where each elements contains the indices for selecting the
    ##   coefficients corresponding to a distributional parameter.
    rng <- list()
    for(j in 1:length(m)) {
        rng[[j]] <- seq(c(1, cumsum(m)[-length(m)] + 1)[j], cumsum(m)[j])
    }
    names(rng) <- names(m)

    ## Implicit SVG
    beta        <- if(is.null(start)) rep(0, sum(m)) else start
    names(beta) <- do.call("c", lapply(X, colnames))
    betaXVec    <- matrix(0, nrow = N, ncol = length(beta))
    colnames(betaXVec) <- names(beta)

    ## grad and link functions
    gfun <- family$score
    lfun <- lapply(family$links, make.link)
    
    zetaVec <- list()
    for(nxi in nx) zetaVec[[nxi]] <- numeric(N)
    ptm <- proc.time()
    for(i in seq.int(N)) {
        cat(sprintf("   * no. obs %i\r", i))

        ## evaluate gammaFun for current iteration
        gamma <- gammaFun(i + i.state) 

        ## predictor
        eta <- list()
        for(nxi in nx) {
            eta[[nxi]] <- drop(beta[rng[[nxi]]] %*% X[[nxi]][shuffle[i],])
        }

        for(nxi in nx) {
            ## find zeta: see slide 110 (Ioannis Big Data Course)
            XCX <- c(X[[nxi]][shuffle[i],, drop = FALSE] %*%
                   CFun(beta[rng[[nxi]]]) %*%
                   t(X[[nxi]][shuffle[i],, drop = FALSE]))
           
            zeta_fun <- make_zeta_fun(y = y[shuffle[i], , drop = FALSE],
                                      eta = eta, XCX = XCX, gfun = gfun,
                                      lfun = lfun, gamma = gamma, parname = nxi)
            upper <- .1
            lower <- -upper 
            root     <- tryCatch(uniroot(zeta_fun, c(lower, upper))$root,
                                 error = function(e) e)
            ## if the first try fails, the interval is enlarged 3 times,
            ##    if no root is found zeta/root is set to 0.
            ierror <- 0
            while(inherits(root, "error")) {
                ierror <- ierror + 1
                if(ierror > 3) {
                    root <- 0
                } else {
                    lower <- lower * 10
                    upper <- upper * 10
                    root  <- tryCatch(uniroot(zeta_fun, c(lower, upper))$root,
                                      error = function(e) e)
                }
            }
            zetaVec[[nxi]][i] <- root

            ## update beta, eta
            beta[rng[[nxi]]] <- beta[rng[[nxi]]] + c(root) * c(X[[nxi]][shuffle[i],] %*%
                                CFun(beta[rng[[nxi]]]))
            eta[[nxi]] <- eta[[nxi]] + root * XCX
        }

        ## keep betapath
        betaXVec[i,] <- beta
    }
    elapsed <- c(proc.time() - ptm)[3]
    cat(sprintf("\n   * runtime = %.3f\n", elapsed))

    rval <- list()
    rval$parameters <- betaXVec
    
    ## fitted values
    rval$fitted.values <- eta

    ## summary
    sgdsum <- list()
    sgdsum$shuffle <- shuffle
    sgdsum$coef    <- beta
    sgdsum$path    <- betaXVec
    sgdsum$y       <- y
    sgdsum$x       <- x
    sgdsum$i.state <- i
    sgdsum$nobs    <- N
    sgdsum$runtime <- elapsed
    sgdsum$zeta    <- zetaVec

    class(sgdsum) <- "sgd.summary"

    rval$sgd.summary <- sgdsum

    rval
}


print.sgd.summary <- function(x, ...) {
    print(x$coef)
    invisible(x)
}

plot.sgd.summary <- function(x, ...) {
    k <- length(x$beta)

    ## coef paths
    matplot(x$path, type = "l", col = colorspace::rainbow_hcl(k), lty = 1)
###    if(!is.null(betaref))
###        abline(h = betaref, col = colorspace::rainbow_hcl(3), lty = 3)

    invisible(x)
}

### helper functions
make_zeta_fun <- function(y, eta, XCX, gfun, lfun, gamma, parname) {

    rfun <- function(zeta) {
        eta[[parname]] <- eta[[parname]] + zeta * XCX

        par <- list()
        for(nxi in names(eta)) { par[[nxi]] <- lfun[[nxi]]$linkinv(eta[[nxi]]) }

        rval <- gamma * gfun[[parname]](y, par) - zeta

        rval
    }

    rfun
}

sgd_grep_X <- function(x) {
    
    X <- list()
    for(nxi in names(x)) {
        X[[nxi]] <- x[[nxi]]$model.matrix
        colnames(X[[nxi]]) <- paste(nxi, "p", colnames(X[[nxi]]), sep = ".")
        for(sci in names(x[[nxi]]$smooth.construct)) {
            xx <- x[[nxi]]$smooth.construct[[sci]]$X
            colnames(xx) <- paste(nxi, "s", sci, 1L:ncol(xx), sep = ".")
            X[[nxi]] <- cbind(X[[nxi]], xx)
        }
    }

    return(X)
}

#opt_sgd.ff <- function(x, y, family, weights = NULL, offset = NULL,
#  gammaFun = function(i) 1/(1+i),
#  shuffle = TRUE, start = NULL, i.state = 0,
#  batch = 1L)
#{
#  nx <- family$names
#  if(!all(nx %in% names(x)))
#    stop("parameter names mismatch with family names!")

#  N  <- nrow(y)
#  y  <- y[[1]]

#  ## Shuffle observations.
#  shuffle_id <- NULL
#  for(i in bamlss_chunk(y)) {
#    ind <- i[1]:i[2]
#    shuffle_id <- ffbase::ffappend(shuffle_id, if(shuffle) sample(ind) else ind)
#  }

#  if(!is.null(start))
#    start <- unlist(start)

#  beta <- list()
#  for(i in nx) {
#    beta[[i]] <- list()
#    if(!is.null(x[[i]]$model.matrix)) {
#      if(!is.null(start)) {
#        beta[[i]][["p"]] <- start[paste0(i, ".p.", colnames(x[[i]]$model.matrix))]
#      } else {
#        beta[[i]][["p"]] <- rep(0, ncol(x[[i]]$model.matrix))
#      }
#      names(beta[[i]][["p"]]) <- colnames(x[[i]]$model.matrix)
#    }
#    if(!is.null(x[[i]]$smooth.construct)) {
#      for(j in names(x[[i]]$smooth.construct)) {
#        ncX <- ncol(x[[i]]$smooth.construct[[j]]$X)
#        if(is.null(start)) {
#          beta[[i]][[paste0("s.", j)]] <- rep(0, ncX)
#        } else {
#          beta[[i]][[paste0("s.", j)]] <- start[paste0(i, ".s.", j, ".b", 1:ncX)]
#        }
#        names(beta[[i]][[paste0("s.", j)]]) <- paste0("b", 1:ncX)
#      }
#    }
#  }

#  ## Init eta.
#  k <- batch
#  eta <- list()
#  for(i in nx) {
#    eta[[i]] <- 0
#    if(!is.null(x[[i]]$model.matrix))
#      eta[[i]] <- eta[[i]] + sum(beta[[i]][["p"]] * x[[i]]$model.matrix[shuffle_id[1:k], ])
#    if(!is.null(x[[i]]$smooth.construct)) {
#      for(j in names(x[[i]]$smooth.construct)) {
#        eta[[i]] <- eta[[i]] + sum(beta[[i]][[paste0("s.", j)]] * x[[i]]$smooth.construct[[j]]$X[shuffle_id[1:k], ])
#      }
#    }
#  }

#  iter <- 1L

#  ptm <- proc.time()
#  while(k <= N) {
#    cat(sprintf("   * no. obs %i\r", k))

#    take <- (k - batch + 1L):k

#    ## Evaluate gammaFun for current iteration.
#    gamma <- gammaFun(iter + i.state)

#    ## Extract response.
#    yn <- y[shuffle_id[take]]

#    for(i in nx) {
#      eta[[i]] <- 0
#      if(!is.null(x[[i]]$model.matrix))
#        eta[[i]] <- eta[[i]] + sum(beta[[i]][["p"]] * x[[i]]$model.matrix[shuffle_id[take], ])
#      if(!is.null(x[[i]]$smooth.construct)) {
#        for(j in names(x[[i]]$smooth.construct)) {
#          eta[[i]] <- eta[[i]] + sum(beta[[i]][[paste0("s.", j)]] * x[[i]]$smooth.construct[[j]]$X[shuffle_id[take], ])
#        }
#      }

#      ## Linear part.
#      if(!is.null(x[[i]]$model.matrix)) {
#        Xn <- x[[i]]$model.matrix[shuffle_id[take], , drop = FALSE]

#        rn <- gamma * family$score[[i]](yn, family$map2par(eta))

#        foo <- function(zeta) {
#          eta[[i]] <- eta[[i]] + drop(Xn %*% (t(Xn) %*% zeta))
#          rval <- gamma * family$score[[i]](yn, family$map2par(eta)) - zeta
#          rval
#        }

#        zeta <- multiroot(foo, start = rn)
#        zeta <- zeta$root

#        beta[[i]][["p"]] <- beta[[i]][["p"]] + drop(t(Xn) %*% zeta)

#        eta[[i]] <- drop(x[[i]]$model.matrix[shuffle_id[take], , drop = FALSE] %*% beta[[i]][["p"]])
#        if(!is.null(x[[i]]$smooth.construct)) {
#          for(j in names(x[[i]]$smooth.construct)) {
#            eta[[i]] <- eta[[i]] + drop(x[[i]]$smooth.construct[[j]]$X[shuffle_id[take], , drop = FALSE] %*% beta[[i]][[paste0("s.", j)]])
#          }
#        }
#      }

#      ## Nonlinear.
#      if(!is.null(x[[i]]$smooth.construct)) {
#        for(j in names(x[[i]]$smooth.construct)) {
#          Xn <- x[[i]]$smooth.construct[[j]]$X[shuffle_id[take], , drop = FALSE]

#          rn <- gamma * family$score[[i]](yn, family$map2par(eta))

#          foo <- function(zeta) {
#            eta[[i]] <- eta[[i]] + drop(Xn %*% (t(Xn) %*% zeta))
#            rval <- gamma * family$score[[i]](yn, family$map2par(eta)) - zeta
#            rval
#          }

#          zeta <- multiroot(foo, start = rn)
#          zeta <- zeta$root

#          ## Update beta, eta.
#          beta[[i]][[paste0("s.", j)]] <- beta[[i]][[paste0("s.", j)]] + drop(t(Xn) %*% zeta)

#          eta[[i]] <- 0
#          if(!is.null(x[[i]]$model.matrix))
#            eta[[i]] <- eta[[i]] + drop(x[[i]]$model.matrix[shuffle_id[take], , drop = FALSE] %*% beta[[i]][["p"]])
#          for(jj in names(x[[i]]$smooth.construct)) {
#            eta[[i]] <- eta[[i]] + drop(x[[i]]$smooth.construct[[jj]]$X[shuffle_id[take], , drop = FALSE] %*% beta[[i]][[paste0("s.", jj)]])
#          }
#        }
#      }
#    }

#    k <- k + batch
#    iter <- iter + 1L
#  }

#  elapsed <- c(proc.time() - ptm)[3]
#  cat(sprintf("\n   * runtime = %.3f\n", elapsed))

#  rval <- list()
#  rval$parameters <- unlist(beta)
#  rval$fitted.values <- eta
#  rval$shuffle <- shuffle
#  rval$runtime <- elapsed

#  rval
#}


opt_bbfit <- bbfit <- function(x, y, family, shuffle = TRUE, start = NULL, offset = NULL,
  epochs = 1, nbatch = 10, verbose = TRUE, ...)
{
  ## Paper: https://openreview.net/pdf?id=ryQu7f-RZ
  ##        https://www.ijcai.org/proceedings/2018/0753.pdf
  aic <- list(...)$aic
  loglik <- list(...)$loglik
  if(is.null(loglik))
    loglik <- FALSE
  if(loglik)
    aic <- FALSE
  if(is.null(aic))
    aic <- FALSE

  eps_loglik <- list(...)$eps_loglik
  if(is.null(eps_loglik))
    eps_loglik <- 0.001

  select <- list(...)$select
  if(is.null(select))
    select <- FALSE

  lasso <- list(...)$lasso
  if(is.null(lasso))
    lasso <- FALSE

  OL <- list(...)$OL
  if(is.null(OL))
    OL <- FALSE
  if(OL)
    lasso <- TRUE

  initialize <- list(...)$initialize
  if(is.null(initialize))
    initialize <- TRUE

  K <- list(...)$K
  if(is.null(K))
    K <- 2

  always <- list(...)$always
  if(is.null(always))
    always <- FALSE

  slice <- list(...)$slice
  if(is.null(slice))
    slice <- FALSE

  nu <- if(is.null(list(...)$nu)) 0.05 else list(...)$nu

  sslice <- NULL
  if(!is.logical(slice)) {
    sslice <- slice
    eps_loglik <- -Inf
    always <- TRUE
    nu <- 1
    slice <- FALSE
  }

  if(slice) {
    eps_loglik <- -Inf
    always <- TRUE
    nu <- 1
  }

  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")

  N  <- nrow(y)

  batch_ids <- list(...)$batch_ids

  if(!is.null(batch_ids)) {
    if(!is.list(batch_ids)) {
      if(length(batch_ids) == 2L) {
        yind <- 1:N
        nb <- floor(batch_ids[1])
        ni <- batch_ids[2]
        batch_ids <- vector(mode = "list", length = ni)
        for(i in 1:ni)
          batch_ids[[i]] <- sample(yind, size = nb, replace = FALSE)
        rm(yind)
      }
    }
  }

  random <- all(nbatch < 1) & all(nbatch > 0)
  batch_select <- srandom <- samp_ids <- FALSE
  if(is.null(batch_ids)) {
    if(!random) {
      batch <- floor(seq.int(1, N, length.out = nbatch + 1L)[-1])
      batch[length(batch)] <- N
      batch <- as.list(batch)
      start0 <- 1L
      for(i in 1:length(batch)) {
        batch[[i]] <- c(start0, batch[[i]])
        start0 <- batch[[i]][-1L] + 1L
      }
    } else {
      if(length(nbatch) < 2L) {
        batch <- floor(N * nbatch)
        batch <- list(c(1, batch), c(batch + 1L, N))
      } else {
        batch <- list(nbatch)
        srandom <- TRUE
        samp_ids <- 1:N
      }
    }
  } else {
    if(is.factor(batch_ids)) {
      batch <- split(1:N, batch_ids)
      batch <- lapply(batch, range)
    } else {
      if(is.list(batch_ids)) {
        batch <- batch_ids
        rm(batch_ids)
      }
    }
    if(!is.list(batch))
      stop("argument batch_ids specified wrong!")
    nbatch <- length(batch)
    batch_select <- TRUE
  }

  if(!is.null(dim(y))) {
    if(ncol(y) < 2)
      y  <- y[[1]]
  }

  noff <- !inherits(y, "ff") #

  if(!is.null(start))
    start <- unlist(start)

  beta <- eta <- etas <- tau2 <- ll_contrib <- medf <- parm <- LLC <- ionly <- list()
  for(i in nx) {
    beta[[i]] <- list()
    tau2[[i]] <- list()
    medf[[i]] <- list()
    parm[[i]] <- list()
    LLC[[i]] <- list()
    ll_contrib[[i]] <- list()
    eta[[i]] <- etas[[i]] <- 0
    if(!is.null(x[[i]]$model.matrix)) {
      ll_contrib[[i]][["p"]] <- medf[[i]][["p.edf"]] <- NA
      LLC[[i]][["p"]] <- 0
      parm[[i]][["p"]] <- matrix(nrow = 0, ncol = ncol(x[[i]]$model.matrix))
      colnames(parm[[i]][["p"]]) <- colnames(x[[i]]$model.matrix)
      if(ncol(x[[i]]$model.matrix) < 2) {
        if(colnames(x[[i]]$model.matrix) == "(Intercept)")
          ionly[[i]] <- TRUE
      } else {
        ionly[[i]] <- FALSE
      }
      beta[[i]][["p"]] <- rep(0, ncol(x[[i]]$model.matrix))
      names(beta[[i]][["p"]]) <- colnames(x[[i]]$model.matrix)
      if(!is.null(start)) {
        start2 <- start[paste0(i, ".p.", colnames(x[[i]]$model.matrix))]
        start2 <- start2[!is.na(start2)]
        names(start2) <- gsub(paste0(i, ".p."), "", names(start2))
        beta[[i]][["p"]][names(start2)] <- start2
      } else {
        if(!is.null(family$initialize) & is.null(offset) & initialize) {
          if(noff) {
            shuffle_id <- sample(seq_len(N))
          } else {
            shuffle_id <- NULL
            for(ii in bamlss_chunk(y)) {
              shuffle_id <- ffbase::ffappend(shuffle_id, if(shuffle) sample(ii) else ii)
            }
          }
          if(!srandom) {
            take <- if(length(batch[[1L]]) > 2) batch[[1L]] else batch[[1L]][1L]:batch[[1L]][2L]
          } else {
            take <- sample(samp_ids, floor(batch[[1L]][1L] * N))
          }
          if(is.null(dim(y))) {
            yn <- y[shuffle_id[take]]
          } else {
            yn <- y[shuffle_id[take], , drop = FALSE]
          }
          if(i %in% names(family$initialize)) {
            if(any(is.na(yn))) {
              stop("NA values in response, please check!")
            }
            yinit <- make.link2(family$links[i])$linkfun(family$initialize[[i]](yn))
            beta[[i]][["p"]]["(Intercept)"] <- mean(yinit, na.rm = TRUE)
          }
        }
      }
      names(beta[[i]][["p"]]) <- colnames(x[[i]]$model.matrix)
    }
    if(!is.null(x[[i]]$smooth.construct)) {
      for(j in names(x[[i]]$smooth.construct)) {
        if(!is.null(x[[i]]$smooth.construct[[j]]$orig.class))
          class(x[[i]]$smooth.construct[[j]]) <- x[[i]]$smooth.construct[[j]]$orig.class
        ll_contrib[[i]][[paste0("s.", j)]] <- medf[[i]][[paste0("s.", j, ".edf")]] <- -1
        LLC[[i]][[paste0("s.", j)]] <- 0
        ncX <- ncol(x[[i]]$smooth.construct[[j]]$X)
        if(is.null(x[[i]]$smooth.construct[[j]]$xt$center))
          x[[i]]$smooth.construct[[j]]$xt$center <- TRUE
        if(inherits(x[[i]]$smooth.construct[[j]], "nnet0.smooth")) {
          tpar <- x[[i]]$smooth.construct[[j]]$state$parameters
          tpar <- tpar[!grepl("tau2", names(tpar))]
          ncX <- length(tpar)
        }
        if(OL) {
          x[[i]]$smooth.construct[[j]]$S <- list()
        }
        ncS <- length(x[[i]]$smooth.construct[[j]]$S) + if(lasso) 1L else 0L
        parm[[i]][[paste0("s.", j)]] <- matrix(nrow = 0L, ncol = ncX + ncS + 1L)
        if(inherits(x[[i]]$smooth.construct[[j]], "nnet0.smooth")) {
          tpar <- x[[i]]$smooth.construct[[j]]$state$parameters
          tpar <- tpar[!grepl("tau2", names(tpar))]
          colnames(parm[[i]][[paste0("s.", j)]]) <- c(names(tpar), paste0("tau2", 1:ncS), "edf")
        } else {
          colnames(parm[[i]][[paste0("s.", j)]]) <- c(paste0("b", 1:ncX), paste0("tau2", 1:ncS), "edf")
        }
        
        if(lasso) {
          lS <- length(x[[i]]$smooth.construct[[j]]$S)
          x[[i]]$smooth.construct[[j]]$S[[lS + 1]] <- function(parameters, ...) {
            b <- get.par(parameters, "b")
            A <- 1 / sqrt(b^2 + 1e-05)
            A <- if(length(A) < 2) matrix(A, 1, 1) else diag(A)
            A
          }
          attr(x[[i]]$smooth.construct[[j]]$S[[lS + 1]], "npar") <- ncX
        }
        if(!inherits(x[[i]]$smooth.construct[[j]], "nnet0.smooth")) {
          if(noff) {
            shuffle_id <- sample(1:N)
          } else {
            shuffle_id <- NULL
            for(ii in bamlss_chunk(y)) {
              shuffle_id <- ffbase::ffappend(shuffle_id, if(shuffle) sample(ii) else ii)
            }
          }
          if(!srandom) {
            if(length(batch[[1L]]) < 3)
              take <- batch[[1L]][1L]:batch[[1L]][2L]
            else
              take <- batch[[1L]]
          } else {
            take <- sample(samp_ids, floor(batch[[1L]][1L] * N))
          }
          Xn <- x[[i]]$smooth.construct[[j]]$X[shuffle_id[take], , drop = FALSE]
          XX <- crossprod(Xn)
          objfun1 <- function(tau2, retedf = FALSE) {
            S <- 0
            for(l in seq_along(x[[i]]$smooth.construct[[j]]$S)) {
              S <- S + 1 / tau2[l] * if(is.function(x[[i]]$smooth.construct[[j]]$S[[l]])) {
                x[[i]]$smooth.construct[[j]]$S[[l]](c("b" = rep(0, ncol(XX))))
              } else {
                x[[i]]$smooth.construct[[j]]$S[[l]]
              }
            }
            edf <- tryCatch(sum_diag(XX %*% matrix_inv(XX + S)), error = function(e) { 2 * ncX })
            if(retedf) {
              return(edf)
            } else {
              tedf <- if(ncX <= 40) {
                4
              } else {
                10
              }
              return((tedf - edf)^2)
            }
          }
          tau2[[i]][[j]] <- rep(0.1, length(x[[i]]$smooth.construct[[j]]$S))
          opt <- try(tau2.optim(objfun1, start = tau2[[i]][[j]],
            scale = 10000, optim = TRUE), silent = TRUE)
          if(!inherits(opt, "try-error")) {
            cat("  .. df", i, "term", x[[i]]$smooth.construct[[j]]$label,
              objfun1(opt, retedf = TRUE), "\n")
            tau2[[i]][[j]] <- opt
          } else {
            cat("  .. df", i, "term", x[[i]]$smooth.construct[[j]]$label, "-1\n")
            tau2[[i]][[j]] <- rep(1, length(x[[i]]$smooth.construct[[j]]$S))
          }
        } else {
          tau2[[i]][[j]] <- rep(100, length(x[[i]]$smooth.construct[[j]]$S))
        }
        if(is.null(start)) {
          beta[[i]][[paste0("s.", j)]] <- rep(0, ncX)
          if(inherits(x[[i]]$smooth.construct[[j]], "nnet0.smooth")) {
            npar <- x[[i]]$smooth.construct[[j]]$state$parameters
            npar <- npar[!grepl("tau2", names(npar))]
            beta[[i]][[paste0("s.", j)]] <- npar
          }
        } else {
          if(inherits(x[[i]]$smooth.construct[[j]], "nnet0.smooth")) {
            start2 <- start[grep(paste0(i, ".s.", j, "."), names(start), fixed = TRUE)]
            start2 <- start2[!grepl("tau2", names(start2))]
          } else {
            start2 <- start[paste0(i, ".s.", j, ".b", 1:ncX)]
          }
          if(any(is.na(start2)))
            stop("dimensions do not match, check starting values!")
          beta[[i]][[paste0("s.", j)]] <- if(all(is.na(start2))) rep(0, ncX) else start2
          if(inherits(x[[i]]$smooth.construct[[j]], "nnet0.smooth")) {
            npar <- x[[i]]$smooth.construct[[j]]$state$parameters
            npar <- npar[!grepl("tau2", names(npar))]
            names(beta[[i]][[paste0("s.", j)]]) <- names(npar)
          }
        }

        if(!inherits(x[[i]]$smooth.construct[[j]], "nnet0.smooth")) {
          names(beta[[i]][[paste0("s.", j)]]) <- paste0("b", 1:ncX)
        }

        x[[i]]$smooth.construct[[j]]$xt[["prior"]] <- "ig"
        x[[i]]$smooth.construct[[j]]$xt[["a"]] <- 0.0001
        x[[i]]$smooth.construct[[j]]$xt[["b"]] <- 10

        priors <- make.prior(x[[i]]$smooth.construct[[j]])
        x[[i]]$smooth.construct[[j]]$prior <- priors$prior
        x[[i]]$smooth.construct[[j]]$grad <- priors$grad
        x[[i]]$smooth.construct[[j]]$hess <- priors$hess
      }
    }
  }
  tbeta <- if(select) beta else NA
  tau2f <- rep(1, length(nx))
  names(tau2f) <- nx

  iter <- 1L

  ptm <- proc.time()
  for(ej in 1:epochs) {
    if(verbose & (epochs > 1))
      cat("starting epoch", ej, "\n")

    ## Shuffle observations.
    if(!batch_select) {
      if(noff) {
        shuffle_id <- sample(1:N)
      } else {
        shuffle_id <- NULL
        for(ii in bamlss_chunk(y)) {
          shuffle_id <- ffbase::ffappend(shuffle_id, if(shuffle) sample(ii) else ii)
        }
      }
    } else {
      shuffle_id <- 1:N
    }

    edf <- NA

    bind <- if(!random & !srandom) {
      seq_along(batch)
    } else 1L

    for(bid in bind) {
      if(!srandom) {
        if(length(batch[[bid]]) > 2) {
          take <- batch[[bid]]
          take2 <- if(bid < 2) {
            batch[[bid + 1L]]
          } else {
            batch[[bid - 1L]]
          }
        } else {
          take <- batch[[bid]][1L]:batch[[bid]][2L]
          take2 <- if(bid < 2) {
            batch[[bid + 1L]][1L]:batch[[bid + 1L]][2L]
          } else {
            batch[[bid - 1L]][1L]:batch[[bid - 1L]][2L]
          }
        }
      } else {
        take <- sample(samp_ids, floor(batch[[bid]][1L] * N))
        take2 <- sample(samp_ids, floor(batch[[bid]][2L] * N))
      }

      ## Extract responses.
      if(is.null(dim(y))) {
        yn <- y[shuffle_id[take]]
        yt <- y[shuffle_id[take2]]
      } else {
        yn <- y[shuffle_id[take], , drop = FALSE]
        yt <- y[shuffle_id[take2], , drop = FALSE]
      }

      for(i in nx) {
        eta[[i]] <- etas[[i]] <- 0
        if(!is.null(x[[i]]$model.matrix)) {
          eta[[i]] <- eta[[i]] + drop(x[[i]]$model.matrix[shuffle_id[take], , drop = FALSE] %*% beta[[i]][["p"]])
          etas[[i]] <- etas[[i]] + drop(x[[i]]$model.matrix[shuffle_id[take2], , drop = FALSE] %*% beta[[i]][["p"]])
        }
        if(!is.null(x[[i]]$smooth.construct)) {
          for(j in names(x[[i]]$smooth.construct)) {
            if(inherits(x[[i]]$smooth.construct[[j]], "nnet0.smooth")) {
              eta[[i]] <- eta[[i]] + x[[i]]$smooth.construct[[j]]$fit.fun(x[[i]]$smooth.construct[[j]]$X[shuffle_id[take], , drop = FALSE],
                beta[[i]][[paste0("s.", j)]])
              etas[[i]] <- etas[[i]] + x[[i]]$smooth.construct[[j]]$fit.fun(x[[i]]$smooth.construct[[j]]$X[shuffle_id[take2], , drop = FALSE],
                beta[[i]][[paste0("s.", j)]])
            } else {
              eta[[i]] <- eta[[i]] + xcenter(x[[i]]$smooth.construct[[j]]$X[shuffle_id[take], , drop = FALSE] %*% beta[[i]][[paste0("s.", j)]])
              etas[[i]] <- etas[[i]] + xcenter(x[[i]]$smooth.construct[[j]]$X[shuffle_id[take2], , drop = FALSE] %*% beta[[i]][[paste0("s.", j)]])
            }
          }
        }
        if(!is.null(offset)) {
          if(i %in% colnames(offset)) {
            eta[[i]] <- eta[[i]] + offset[shuffle_id[take], i]
            etas[[i]] <- etas[[i]] + offset[shuffle_id[take2], i]
          }
        }
      }

      eta00 <- eta
      edf <- 0

      for(i in nx) {
        ## Linear part.
        if(!is.null(x[[i]]$model.matrix)) {
          Xn <- x[[i]]$model.matrix[shuffle_id[take], , drop = FALSE]
          Xt <- x[[i]]$model.matrix[shuffle_id[take2], , drop = FALSE]

          peta <- family$map2par(eta)
          petas <- family$map2par(etas)

          ll0 <- family$loglik(yt, petas)

          score <- process.derivs(family$score[[i]](yn, peta, id = i), is.weight = FALSE)
          hess <- process.derivs(family$hess[[i]](yn, peta, id = i), is.weight = TRUE)

          scores <- process.derivs(family$score[[i]](yt, petas, id = i), is.weight = FALSE)
          hesss <- process.derivs(family$hess[[i]](yt, petas, id = i), is.weight = TRUE)

          b0 <- beta[[i]][["p"]]

          eta_0 <- eta[[i]]
          etas_0 <- etas[[i]]

          z <- eta[[i]] + 1/hess * score
          zs <- etas[[i]] + 1/hesss * scores

          eta[[i]] <- eta[[i]] - drop(Xn %*% b0)
          e <- z - eta[[i]]

          XWX <- crossprod(Xn * hess, Xn)
          I <- diag(1, ncol(XWX))
          I[1, 1] <- 1e-10

          if(!ionly[[i]]) {
            if(ncol(I) > 1)
              I[1, 1] <- 0
          }

          etas[[i]] <- etas[[i]] - drop(Xt %*% b0)

          objfun2 <- function(tau2, retLL = FALSE, step = FALSE) {
            P <- matrix_inv(XWX + 1/tau2 * I)
            if(step) {
              b <- b0 + nu * (drop(P %*% crossprod(Xn * hess, e)) - b0)
            } else {
              b <- drop(P %*% crossprod(Xn * hess, e))
            }
            etas[[i]] <- etas[[i]] + drop(Xt %*% b)
            if(retLL) {
              return(family$loglik(yt, family$map2par(etas)))
            }
            if(aic | loglik) {
              if(aic) {
                iedf <- sum_diag(XWX %*% P)
                ll <- -2 * family$loglik(yt, family$map2par(etas)) + K * iedf
              } else {
                ll <- -1 * family$loglik(yt, family$map2par(etas))
              }
            } else {
              ll <- mean((zs - etas[[i]])^2, na.rm = TRUE)
            }
            return(ll)
          }
          if(ionly[[i]]) {
            tau2fe <- 1e+10
          } else {
            tau2fe <- try(tau2.optim(objfun2, tau2f[i], optim = TRUE), silent = TRUE)
          }

          ll_contrib[[i]][["p"]] <- NA

          if(!inherits(tau2fe, "try-error")) {
            ll1 <- objfun2(tau2fe, retLL = TRUE, step = TRUE)
            epsll <- (ll1 - ll0)/abs(ll0)
            if(is.na(epsll)) {
              ll1 <- ll0 <- 1
              epsll <- -1
            }
            accept <- epsll >= -0.5
            if((((ll1 > ll0) & (epsll > eps_loglik)) | always) & accept) {
              tau2f[i] <- tau2fe
              P <- matrix_inv(XWX + 1/tau2f[i] * I)
              if(select) {
                tbeta[[i]][["p"]] <- b0 + nu * (drop(P %*% crossprod(Xn * hess, e)) - b0)
              } else {
                beta[[i]][["p"]] <- b0 + nu * (drop(P %*% crossprod(Xn * hess, e)) - b0)
              }
              tedf <- sum_diag(XWX %*% P)
              edf <- edf + tedf
              ll_contrib[[i]][["p"]] <- ll1 - ll0
              medf[[i]][["p.edf"]] <- c(medf[[i]][["p.edf"]], tedf)
            }
          }

          if(!select) {
            eta[[i]] <- eta[[i]] + drop(Xn %*% beta[[i]][["p"]])
            etas[[i]] <- etas[[i]] + drop(Xt %*% beta[[i]][["p"]])
          } else {
            eta[[i]] <- eta_0
            etas[[i]] <- etas_0
          }
        }

        ## Nonlinear.
        if(!is.null(x[[i]]$smooth.construct)) {
          for(j in names(x[[i]]$smooth.construct)) {
            Xn <- x[[i]]$smooth.construct[[j]]$X[shuffle_id[take], , drop = FALSE]
            Xt <- x[[i]]$smooth.construct[[j]]$X[shuffle_id[take2], , drop = FALSE]

            b0 <- beta[[i]][[paste0("s.", j)]]

            if(inherits(x[[i]]$smooth.construct[[j]], "nnet0.smooth")) {
              Xn <- x[[i]]$smooth.construct[[j]]$getZ(Xn, b0)
              Xt <- x[[i]]$smooth.construct[[j]]$getZ(Xt, b0)
              b0 <- b0[1:ncol(Xn)]
            }

            eta_0 <- eta[[i]]
            etas_0 <- etas[[i]]

            peta <- family$map2par(eta)
            petas <- family$map2par(etas)

            ll0 <- family$loglik(yt, petas)

            score <- process.derivs(family$score[[i]](yn, peta, id = i), is.weight = FALSE)
            hess <- process.derivs(family$hess[[i]](yn, peta, id = i), is.weight = TRUE)

            scores <- process.derivs(family$score[[i]](yt, petas, id = i), is.weight = FALSE)
            hesss <- process.derivs(family$hess[[i]](yt, petas, id = i), is.weight = TRUE)

            z <- eta[[i]] + 1/hess * score
            zs <- etas[[i]] + 1/hesss * scores

            if(x[[i]]$smooth.construct[[j]]$xt$center) {
              eta[[i]] <- eta[[i]] - xcenter(Xn %*% b0)
            } else {
              eta[[i]] <- eta[[i]] - drop(Xn %*% b0)
            }

            e <- z - eta[[i]]

            if(x[[i]]$smooth.construct[[j]]$xt$center) {
              etas[[i]] <- etas[[i]] - xcenter(Xt %*% b0)
            } else {
              etas[[i]] <- etas[[i]] - drop(Xt %*% b0)
            }

            wts <- NULL
            if(inherits(x[[i]]$smooth.construct[[j]], "nnet0.smooth")) {
              wts <- unlist(x[[i]]$smooth.construct[[j]]$sample_weights(
                x = x[[i]]$smooth.construct[[j]]$X[shuffle_id[take], -1, drop = FALSE],
                y = e, weights = hess, wts = beta[[i]][[paste0("s.", j)]])
              )
              Xn <- x[[i]]$smooth.construct[[j]]$getZ(x[[i]]$smooth.construct[[j]]$X[shuffle_id[take], , drop = FALSE], wts)
              Xt <- x[[i]]$smooth.construct[[j]]$getZ(x[[i]]$smooth.construct[[j]]$X[shuffle_id[take2], , drop = FALSE], wts)
            }

            XWX <- crossprod(Xn * hess, Xn)

            objfun3 <- function(tau2, retLL = FALSE, step = FALSE) {
              S <- 0
              for(l in 1:length(tau2)) {
                S <- S + 1/tau2[l] * if(is.function(x[[i]]$smooth.construct[[j]]$S[[l]])) {
                  x[[i]]$smooth.construct[[j]]$S[[l]](c(b0, x[[i]]$smooth.construct[[j]]$fixed.hyper))
                } else {
                  x[[i]]$smooth.construct[[j]]$S[[l]]
                }
              }
              P <- matrix_inv(XWX + S)
              if(step) {
                b <- b0 + nu * (drop(P %*% crossprod(Xn * hess, e)) - b0)
              } else {
                b <- drop(P %*% crossprod(Xn * hess, e))
              }

              if(x[[i]]$smooth.construct[[j]]$xt$center) {
                etas[[i]] <- etas[[i]] + xcenter(Xt %*% b)
              } else {
                etas[[i]] <- etas[[i]] + drop(Xt %*% b)
              }
              if(retLL) {
                return(family$loglik(yt, family$map2par(etas)))
              }
              if(aic | loglik) {
                if(aic) {
                  iedf <- sum_diag(XWX %*% P)
                  ll <- -2 * family$loglik(yt, family$map2par(etas)) + K * iedf
                } else {
                  names(b) <- paste0("b", 1:length(b))
                  names(tau2) <- paste0("tau2", 1:length(tau2))
                  ll <- -1 * (family$loglik(yt, family$map2par(etas)) + x[[i]]$smooth.construct[[j]]$prior(c(b, tau2)))
                }
              } else {
                ll <- mean((zs - etas[[i]])^2, na.rm = TRUE)
              }
              return(ll)
            }

            if(!is.null(sslice)) {
              if(iter > sslice)
                slice <- TRUE
            }

            if(!slice) {
              tau2s <- try(tau2.optim(objfun3, tau2[[i]][[j]], optim = TRUE), silent = TRUE)
            } else {
              theta <- c(b0, "tau2" = tau2[[i]][[j]])
              ii <- grep("tau2", names(theta))
              logP <- function(g, x, ll, ...) {
                -1 * objfun3(get.par(g, "tau2"))
              }
              sok <- TRUE
              for(jj in ii) {
                theta <- try(uni.slice(theta, x[[i]]$smooth.construct[[j]], family, NULL,
                  NULL, i, jj, logPost = logP, lower = 0, ll = ll0), silent = TRUE)
                if(inherits(theta, "try-error")) {
                  sok <- FALSE
                  break
                }
              }
              if(sok) {
                tau2s <- as.numeric(get.par(theta, "tau2"))
              } else {
                tau2s <- NA
                class(tau2s) <- "try-error"
              }
            }
            ll_contrib[[i]][[paste0("s.", j)]] <- NA
            accept <- TRUE
            if(!inherits(tau2s, "try-error")) {
              ll1 <- objfun3(tau2s, retLL = TRUE, step = TRUE)
              epsll <- (ll1 - ll0)/abs(ll0)
              if(is.na(epsll)) {
                ll1 <- ll0 <- 1
                epsll <- -1
              }
#              if(!slice) {
#                accept <- TRUE
#              } else {
#                epsll < 0.5
#              }
#              if(!always) {
                accept <- epsll >= -0.5
#              } else {
#                accept <- TRUE
#              }
              if((((ll1 > ll0) & (epsll > eps_loglik)) | always) & accept) {
                tau2[[i]][[j]] <- tau2s
                S <- 0
                for(l in 1:length(tau2[[i]][[j]])) {
                  S <- S + 1/tau2[[i]][[j]][l] * if(is.function(x[[i]]$smooth.construct[[j]]$S[[l]])) {
                    x[[i]]$smooth.construct[[j]]$S[[l]](c(b0, x[[i]]$smooth.construct[[j]]$fixed.hyper))
                  } else {
                    x[[i]]$smooth.construct[[j]]$S[[l]]
                  }
                }

                P <- matrix_inv(XWX + S)

                if(select) {
                  tbeta[[i]][[paste0("s.", j)]] <- b0 + nu * (drop(P %*% crossprod(Xn * hess, e)) -  b0)
                  if(!is.null(wts)) {
                    names(tbeta[[i]][[paste0("s.", j)]]) <- paste0("bb", 1:length(tbeta[[i]][[paste0("s.", j)]]))
                    tbeta[[i]][[paste0("s.", j)]] <- c(tbeta[[i]][[paste0("s.", j)]], wts)
                  }
                } else {
                  beta[[i]][[paste0("s.", j)]] <- b0 + nu * (drop(P %*% crossprod(Xn * hess, e)) - b0)
                  if(!is.null(wts)) {
                    names(beta[[i]][[paste0("s.", j)]]) <- paste0("bb", 1:length(beta[[i]][[paste0("s.", j)]]))
                    beta[[i]][[paste0("s.", j)]] <- c(beta[[i]][[paste0("s.", j)]], wts)
                  }
                }
                tedf <- sum_diag(XWX %*% P)
                edf <- edf + tedf
                ll_contrib[[i]][[paste0("s.", j)]] <- ll1 - ll0
                medf[[i]][[paste0("s.", j, ".edf")]] <- c(medf[[i]][[paste0("s.", j, ".edf")]], tedf)
              }
            } else {
              warning(paste0("check distribution of term ", j, "!"))
            }
            if(!select & accept) {
              if(inherits(x[[i]]$smooth.construct[[j]], "nnet0.smooth")) {
                nid <- 1:x[[i]]$smooth.construct[[j]]$nodes
                if(x[[i]]$smooth.construct[[j]]$xt$center) {
                  eta[[i]] <- eta[[i]] + xcenter(Xn %*% beta[[i]][[paste0("s.", j)]][nid])
                  etas[[i]] <- etas[[i]] + xcenter(Xt %*% beta[[i]][[paste0("s.", j)]][nid])
                } else {
                  eta[[i]] <- eta[[i]] + drop(Xn %*% beta[[i]][[paste0("s.", j)]][nid])
                  etas[[i]] <- etas[[i]] + drop(Xt %*% beta[[i]][[paste0("s.", j)]][nid])
                }
#fit <- Xn %*% beta[[i]][[paste0("s.", j)]][nid]
#Z <- x[[i]]$smooth.construct[[j]]$X[shuffle_id[take], , drop = FALSE]
#plot(Z[, 2], e)
#plot2d(fit ~ Z[,2], add = TRUE)
              } else {
                if(x[[i]]$smooth.construct[[j]]$xt$center) {
                  eta[[i]] <- eta[[i]] + xcenter(Xn %*% beta[[i]][[paste0("s.", j)]])
                  etas[[i]] <- etas[[i]] + xcenter(Xt %*% beta[[i]][[paste0("s.", j)]])
                } else {
                  eta[[i]] <- eta[[i]] + drop(Xn %*% beta[[i]][[paste0("s.", j)]])
                  etas[[i]] <- etas[[i]] + drop(Xt %*% beta[[i]][[paste0("s.", j)]])
                }
#fit <- Xn %*% beta[[i]][[paste0("s.", j)]]
#Z <- d$x2[shuffle_id[take]]
#plot(Z, e)
#plot2d(fit ~ Z, add = TRUE)
              }
            }
            if(!accept | select) {
              eta[[i]] <- eta_0
              etas[[i]] <- etas_0
            }
          }
        }
      }

      if(select) {
        llc <- unlist(ll_contrib)
        if(!all(is.na(llc))) {
          llval <- max(llc, na.rm = TRUE)
          llc <- names(llc)[which.max(llc)]
          llc <- strsplit(llc, ".", fixed = TRUE)[[1]]
          llc <- c(llc[1], paste0(llc[-1], collapse = "."))
          beta[[llc[1]]][[llc[2]]] <- b0 <- tbeta[[llc[1]]][[llc[2]]]
          csm <- TRUE
          if(llc[2] != "p") {
            llc2 <- gsub("s.", "", llc[2], fixed = TRUE)
            Xn <- x[[llc[1]]]$smooth.construct[[llc2]]$X[shuffle_id[take], , drop = FALSE]
            Xt <- x[[llc[1]]]$smooth.construct[[llc2]]$X[shuffle_id[take2], , drop = FALSE]
            if(inherits(x[[llc[1]]]$smooth.construct[[llc2]], "nnet0.smooth")) {
              Xn <- x[[llc[1]]]$smooth.construct[[llc2]]$getZ(Xn, b0)
              Xt <- x[[llc[1]]]$smooth.construct[[llc2]]$getZ(Xt, b0)
              b0 <- b0[1:ncol(Xn)]
            }
            csm <- x[[llc[1]]]$smooth.construct[[llc2]]$xt$center
          } else {
            csm <- FALSE
            Xn <- x[[llc[1]]]$model.matrix[shuffle_id[take], , drop = FALSE]
            Xt <- x[[llc[1]]]$model.matrix[shuffle_id[take2], , drop = FALSE]
          }
          ll_iter <- attr(LLC[[llc[1]]][[llc[2]]], "iteration")
          ll_iter <- c(ll_iter, iter)
          LLC[[llc[1]]][[llc[2]]] <- c(LLC[[llc[1]]][[llc[2]]], llval)
          attr(LLC[[llc[1]]][[llc[2]]], "iteration") <- ll_iter
          if(csm) {
            eta[[llc[1]]] <- eta[[llc[1]]] + xcenter(Xn %*% b0)
            etas[[llc[1]]] <- etas[[llc[1]]] + xcenter(Xt %*% b0)
          } else {
            eta[[llc[1]]] <- eta[[llc[1]]] + drop(Xn %*% b0)
            etas[[llc[1]]] <- etas[[llc[1]]] + drop(Xt %*% b0)
          }
        }
      }

      for(i in nx) {
        for(j in names(parm[[i]])) {
          jj <- paste0(strsplit(j, ".", fixed = TRUE)[[1]][-1], collapse = ".")
          tedf <- medf[[i]][[paste0(j, ".edf")]]
          tpar <-  if(j != "p") {
            c(beta[[i]][[j]], tau2[[i]][[jj]], tedf[length(tedf)])
          } else {
            beta[[i]][[j]]
          }
          names(tpar) <- NULL
          parm[[i]][[j]] <- rbind(parm[[i]][[j]], tpar)
        }
      }

      eta00 <- do.call("cbind", eta00)
      eta01 <- do.call("cbind", eta)

      if(iter < 2L)
        eta00[abs(eta00) < 1e-20] <- 1e-20

      eps <- mean(abs((eta01 - eta00) / eta00), na.rm = TRUE)

      if(verbose) {
        edf <- abs(edf)
        btxt <- if(srandom) {
          NA
        } else {
          if(length(batch[[bid]]) > 2) {
            length(batch[[bid]]) * iter
          } else {
            batch[[bid]][2L]
          }
        }
        if(iter < 2) {
          cat(sprintf("   * iter %i, no. obs %i, edf %f\r", iter, btxt, round(edf, 4)))
        } else {
          cat(sprintf("   * iter %i, no. obs %i, eps %f, edf %f\r", iter, btxt, round(eps, 4), round(edf, 2)))
        }
      }

      iter <- iter + 1L
    }

    if(verbose)
      cat("\n")
  }

  elapsed <- c(proc.time() - ptm)[3]

  if(verbose) {
    cat("\n")
    et <- if(elapsed > 60) {
      paste(formatC(format(round(elapsed / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(elapsed, 2), nsmall = 2), width = 5), "sec", sep = "")
    cat("elapsed time: ", et, "\n", sep = "")
  }

  for(i in nx) {
    for(j in seq_along(medf[[i]])) {
      medf[[i]][[j]] <-  if(all(is.na(medf[[i]][[j]]))) {
        0
      } else median(medf[[i]][[j]], na.rm = TRUE)
    }
    for(j in names(parm[[i]])) {
      colnames(parm[[i]][[j]]) <- paste0(i, ".", j, ".", colnames(parm[[i]][[j]]))
    }
    parm[[i]] <- do.call("cbind", parm[[i]])
  }
  parm <- do.call("cbind", parm)
  rownames(parm) <- NULL
  ##if(nrow(parm) > 1L)
  ##  parm <- parm[-1L, , drop = FALSE]

  rval <- list()
  rval$parameters <- c(unlist(beta), unlist(medf))
  rval$fitted.values <- eta
  rval$shuffle <- shuffle
  rval$runtime <- elapsed
  rval$edf <- edf
  rval$nbatch <- nbatch
  rval$parpaths <- parm
  rval$epochs <- epochs
  rval$n.iter <- iter
  if(select) {
    rval$llcontrib <- LLC
  }

  rval
}

contribplot <- function(x, ...) {
  if(is.null(ll <- x$model.stats$optimizer$llcontrib))
    stop("nothing to plot")
  iter <- x$model.stats$optimizer$n.iter - 1L
  iter2 <- NULL
  sf <- list()
  for(i in names(ll)) {
    sf[[i]] <- list()
    for(j in names(ll[[i]])) {
      if(!is.null(ll[[i]][[j]])) {
        ii <- attr(ll[[i]][[j]], "iteration")
        sf[[i]][[j]] <- length(ii)
        iter2 <- c(iter2, length(ii))
        llv <- rep(0, iter)
        llv[ii] <- ll[[i]][[j]][-1]
        llv <- cumsum(llv)
        ll[[i]][[j]] <- c(0, llv)
      } else {
        ll[[i]][[j]] <- rep(0, iter)
        sf[[i]][[j]] <- 0
      }
    }
    sf[[i]] <- do.call("rbind", sf[[i]])
    sf[[i]] <- sf[[i]][order(sf[[i]][, 1], decreasing = TRUE), , drop = FALSE]
    colnames(sf[[i]]) <- "Sel. freq."
    ll[[i]] <- do.call("cbind", ll[[i]])
    colnames(ll[[i]]) <- paste0(i, ".", colnames(ll[[i]]))
  }
  iter2 <- sum(iter2)
  for(i in names(ll)) {
    sf[[i]] <- sf[[i]]/iter2
    cat(i, "\n", sep = "")
    printCoefmat(sf[[i]])
    cat("\n")
  }
  ll <- do.call("cbind", ll)
  print.boost_summary(list("loglik" = ll, "mstop" = iter),
    summary = FALSE, which = "loglik.contrib", ...)
  invisible(list("loglik" = ll, "selfreqs" = sf))
}


opt_bbfitp <- bbfitp <- function(x, y, family, mc.cores = 1, ...)
{
  seeds <- ceiling(runif(mc.cores, 1, 1000000))
  parallel_fun <- function(i) {
    set.seed(seeds[i])
    opt_bbfit(x, y, family, ...)
  }
  b <- parallel::mclapply(1:mc.cores, parallel_fun, mc.cores = mc.cores)
  rval <- list()
  rval$samples <- lapply(b, function(x) {
    if(inherits(x, "try-error")) {
      writeLines(x)
      return(x)
    } else {
      return(as.mcmc(x$parpaths))
    }
  })
  is_err <- sapply(rval$samples, is.character)
  if(all(is_err))
    stop("something went wrong in bbfitp()!")
  if(any(is_err))
    warning("one core reports an error.")
  b <- b[!is_err]
  rval$samples <- as.mcmc.list(rval$samples[!is_err])
  rval$parameters <- colMeans(do.call("rbind", lapply(b, function(x) x$parpaths)))
  rval$nbatch <- b[[1]]$nbatch
  rval$runtime <- mean(sapply(b, function(x) x$runtime))
  rval$epochs <- b[[1]]$epochs
  rval
}


bbfit_plot <- function(x, name = NULL, ...)
{
  x <- x$model.stats$optimizer$parpaths
  if(is.null(x)) {
    warning("there is nothing to plot")
    return(invisible(NULL))
  }
  if(!is.null(name)) {
    for(i in name) {
      x <- x[, grep(i, colnames(x), fixed = TRUE), drop = FALSE]
    }
  }
  cn <- colnames(x)
  cn2 <- strsplit(cn, ".", fixed = TRUE)
  cn2 <- lapply(cn2, function(x) { paste0(x[-length(x)], collapse = ".") })
  cn2 <- as.factor(unlist(cn2))
  cat(levels(cn2), "\n")
  col <- rainbow_hcl(nlevels(cn2))[cn2]
  matplot(x, type = "l", lty = 1, xlab = "Iteration", ylab = "Coefficients", col = col, ...)
  return(invisible(x))
}

