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
          x$smooth.construct[[j]] <- assign.df(x$smooth.construct[[j]], tdf)
          if(!is.null(x$smooth.construct[[j]]$xt[["update"]]))
            x$smooth.construct[[j]]$update <- x$smooth.construct[[j]]$xt[["update"]]
          if(is.null(x$smooth.construct[[j]]$update)) {
            if(is.character(update)) {
              if(!grepl("bfit_", update))
                update <- paste("bfit", update, sep = "_")
              update <- get(update)
            }
            x$smooth.construct[[j]]$update <- update
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
    if(what %in% c("tau2", "tau", "lambda")) {
      p <- x$state$parameters
      return(p[grep("tau", names(p))])
    } else {
      if(what %in% "b") {
        p <- x$state$parameters
        return(p[!grepl("tau", names(p)) & !grepl("edf", names(p)) & !grepl("lasso", names(p))])
      } else return(x$state[[what]])
    }
  }
}

get.par <- function(x, what = NULL) {
  if(is.null(what) | is.null(names(x))) return(x)
  if(what %in% c("tau2", "tau", "lambda")) {
    return(x[grep("tau", names(x))])
  } else {
    if(what %in% "b") {
      return(x[!grepl("tau", names(x)) & !grepl("edf", names(x)) & !grepl("lasso", names(x))])
    } else return(x[what])
  }
}

set.par <- function(x, replacement, what) {
  if(is.null(replacement))
    return(x)
  if(what %in% c("tau2", "tau", "lambda")) {
    x[grep("tau", names(x))] <- replacement
  } else {
    if(what %in% "b") {
      if(as.integer(sum(!grepl("tau", names(x)) & !grepl("edf", names(x)) & !grepl("lasso", names(x)))) != length(replacement)) {
        stop("here")
      }
      x[!grepl("tau", names(x)) & !grepl("edf", names(x)) & !grepl("lasso", names(x))] <- replacement
    } else x[what] <- replacement
  }
  x
}

## The default method.
bamlss.engine.setup.smooth.default <- function(x, spam = FALSE, Matrix = FALSE, ...)
{
  if(inherits(x, "special"))
    return(x)
  if(is.null(x$binning) & !is.null(x$xt[["binning"]])) {
    x$binning <- match.index(x$X)
    x$binning$order <- order(x$binning$match.index)
    x$binning$sorted.index <- x$binning$match.index[x$binning$order]
    x$X <- x$X[x$binning$nodups, , drop = FALSE]
  }
  if(!is.null(x$binning)) {
    if(nrow(x$X) != length(x$binning$nodups))
      x$X <- x$X[x$binning$nodups, , drop = FALSE]
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
      edf <- sum_diag(x$state$XX %*% P)
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

  force.spam <- if(is.null(args$force.spam)) FALSE else args$force.spam
  if(!is.null(x$sparse.setup$crossprod)) {
    if((ncol(x$sparse.setup$crossprod) < ncol(x$X) * 0.5) & force.spam)
      spam <- TRUE
    if(spam) {
      x$X <- as.spam(x$X)
      xx <- crossprod.spam(x$X)
      for(j in seq_along(x$S)) {
        x$S[[j]] <- as.spam( if(is.function(x$S[[j]])) x$S[[j]](c("b" = rep(0, attr(x$S[[j]], "npar")))) else x$S[[j]] )
        xx <- xx + x$S[[j]]
      }
      x$sparse.setup$spam.cholFactor <- chol.spam(xx)
      if(force.spam) {
        x$update <- bfit_iwls_spam
      }
    }
  }

  force.Matrix <- if(is.null(args$force.Matrix)) FALSE else args$force.Matrix
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
      if(!is.null(x$xt$lambda))
        tau2 <- 1 / x$xt$lambda
    }
    if(!is.null(tau2)) {
      tau2 <- rep(tau2, length.out = ntau2)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    }
  }
  pid <- !grepl("tau", names(x$state$parameters)) & !grepl("edf", names(x$state$parameters))
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

  if(is.null(x$all_diagonal)) {
    XX <- crossprod(x$X)
    XX_is_diagonal <- all(XX[!diag(nrow(XX))] == 0)

    b <- runif(length(get.par(x$state$parameters, "b")))
    tau2 <- get.par(x$state$parameters, "tau2")
    S <- 0
    for(j in seq_along(tau2))
       S <- S + 1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](c("b" = b, "tau2" = tau2)) else x$S[[j]]
    S_is_diagonal <- all(S[!diag(nrow(S))] == 0)
    x$all_diagonal <- XX_is_diagonal & (if(x$fixed) TRUE else S_is_diagonal)
  }

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
assign.df <- function(x, df)
{
  if(inherits(x, "special"))
    return(x)
  tau2 <- get.par(x$state$parameters, "tau2")
  if(x$fixed | !length(tau2))
    return(x)
  if(!is.null(x$fxsp)) {
    if(x$fxsp)
      return(x)
  }
  df <- if(is.null(x$xt$df)) df else x$xt$df
  if(is.null(df)) {
    nc <- ncol(x$X)
    df <- ceiling(nc * 0.5)
  }
  if(df > ncol(x$X))
    df <- ncol(x$X)
  if(df < 1)
    df <- 1
  int <- c(.Machine$double.eps^0.25, 1e+10)
  if(inherits(x$X, "spam")) {
    XX <- crossprod.spam(x$X)
  } else {
    XX <- crossprod(x$X)
  }
  if(length(tau2) > 1) {
    if(FALSE) {
      df.part <- df / length(tau2)
      for(j in seq_along(tau2)) {
        objfun <- function(val) {
          tau2[j] <- val
          S <- 0
          for(i in seq_along(x$S))
            S <- S + 1 / tau2[i] * (if(is.function(x$S[[i]])) x$S[[i]](c("b" = rep(0, attr(x$S[[i]], "npar")))) else x$S[[i]])
          edf <- sum_diag(XX %*% matrix_inv(XX + S, index = x$sparse.setup))
          return((df - edf)^2)
        }
        opt <- try(optimize(objfun, int)$minimum, silent = TRUE)
        if(!inherits(opt, "try-error"))
          tau2[j] <- opt
      }
      tau2 <- rep(1000, length(tau2))
    }
  } else {
    objfun <- function(tau2) {
      edf <- sum_diag(XX %*% matrix_inv(XX + 1 / tau2 * (if(is.function(x$S[[1]])) {
          x$S[[1]](c("b" = rep(0, attr(x$S[[1]], "npar"))))
        } else x$S[[1]]), index = x$sparse.setup))
      return((df - edf)^2)
    }
    tau2 <- try(optimize(objfun, int)$minimum, silent = TRUE)
    if(inherits(tau2, "try-error"))
      return(x)
  }
  x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
  x$state$edf <- df
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
#  for(i in bit::chunk(x)) {
#    res <- ffappend(res, FUN(x[i, ]))
#  }
#  res
## FIXME: ff support!
  FUN(x)
}

ffdf_eval_sh <- function(y, par, FUN)
{
#  res <- NULL
#  for(i in bit::chunk(y)) {
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
#  for(i in bit::chunk(x)) {
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
fmt <- function(x, width = 8, digits = 2) {
  txt <- formatC(round(x, digits), format = "f", digits = digits , width = width)
  if(nchar(txt) > width) {
    txt <- strsplit(txt, "")[[1]]
    txt <- paste(txt[1:width], collapse = "", sep = "")
  }
  txt
}

bfit <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
  update = "iwls", criterion = c("AICc", "BIC", "AIC"),
  eps = .Machine$double.eps^0.25, maxit = 400,
  outer = FALSE, inner = FALSE, mgcv = FALSE,
  verbose = TRUE, digits = 4, flush = TRUE, nu = NULL, stop.nu = NULL, ...)
{
  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("design construct names mismatch with family names!")

  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, update = update, ...)

  criterion <- match.arg(criterion)
  np <- length(nx)

  if(!is.null(nu)) {
    if(nu < 0)
      nu <- NULL
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
  if(is.null(start))
    eta <- init.eta(eta, y, family, nobs)
 
  if(!is.null(weights))
    weights <- as.data.frame(weights)
  if(!is.null(offset)) {
    offset <- as.data.frame(offset)
    for(j in nx) {
      if(!is.null(offset[[j]]))
        eta[[j]] <- eta[[j]] + offset[[j]]
    }
  }

  ia <- if(flush) interactive() else FALSE

  if(mgcv) {
    outer <- TRUE
    inner <- TRUE
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
            lpost0 <- family$loglik(y, family$map2par(eta)) + logprior
            lp <- logprior - x[[sj]]$prior(x[[sj]]$state$parameters)

            eta2 <- eta
            eta2[[id]] <- eta2[[id]] - x[[sj]]$state$fitted.values

            b0 <- get.par(x[[sj]]$state$parameters, "b")
            b1 <- get.par(p.state$parameters, "b")

            objfun <- function(nu, diff = TRUE) {
              p.state$parameters <- set.par(p.state$parameters, nu * b1 + (1 - nu) * b0, "b")
              eta2[[id]] <- eta2[[id]] + x[[sj]]$fit.fun(x[[sj]]$X,
                get.par(p.state$parameters, "b"))
              lp2 <- family$loglik(y, family$map2par(eta2)) + lp + x[[sj]]$prior(p.state$parameters)
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

              p.state$parameters <- set.par(p.state$parameters, nuo * b1 + (1 - nuo) * b0, "b")
              p.state$fitted.values <- x[[sj]]$fit.fun(x[[sj]]$X,
                get.par(p.state$parameters, "b"))
              eta2[[id]] <- eta2[[id]] + p.state$fitted.values
              lpost1 <- family$loglik(y, family$map2par(eta2)) + lp + x[[sj]]$prior(p.state$parameters)
              if(lpost1 < lpost0) {
                warning(paste("logPost is decreasing updating term: ", id, ", ",
                  x[[sj]]$label, "; diff: ", lpost1 - lpost0, sep = ""))
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
    ptm <- proc.time()
    while(eps0 > eps & iter < maxit) {
      eta0 <- eta
      ## Cycle through all parameters
      for(j in 1:np) {
        if(outer | iter < 2) {
          peta <- family$map2par(eta)

          if(no_ff) {
            ## Compute weights.
            hess <- process.derivs(family$hess[[nx[j]]](y, peta, id = nx[j]), is.weight = TRUE)

            ## Score.
            score <- process.derivs(family$score[[nx[j]]](y, peta, id = nx[j]), is.weight = FALSE)
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

        if(iter < 2)
          eta[[nx[j]]] <- get.eta(x)[[nx[j]]]

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
            if(!is.null(nu)) {
              lp0 <- get.log.prior(x)
              lpost0 <- family$loglik(y, family$map2par(eta)) + lp0
              lp <- lp0 - x[[nx[j]]]$smooth.construct[[sj]]$prior(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters)

              eta2 <- eta
              eta2[[nx[j]]] <- eta2[[nx[j]]] - x[[nx[j]]]$smooth.construct[[sj]]$state$fitted.values

              b0 <- get.par(x[[nx[j]]]$smooth.construct[[sj]]$state$parameters, "b")
              b1 <- get.par(p.state$parameters, "b")

              objfun <- function(nu, diff = TRUE) {
                p.state$parameters <- set.par(p.state$parameters, nu * b1 + (1 - nu) * b0, "b")
                eta2[[nx[j]]] <- eta2[[nx[j]]] + x[[nx[j]]]$smooth.construct[[sj]]$fit.fun(x[[nx[j]]]$smooth.construct[[sj]]$X,
                  get.par(p.state$parameters, "b"))
                lp2 <- family$loglik(y, family$map2par(eta2)) + lp + x[[nx[j]]]$smooth.construct[[sj]]$prior(p.state$parameters)
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

                p.state$parameters <- set.par(p.state$parameters, nuo * b1 + (1 - nuo) * b0, "b")
                p.state$fitted.values <- x[[nx[j]]]$smooth.construct[[sj]]$fit.fun(x[[nx[j]]]$smooth.construct[[sj]]$X,
                  get.par(p.state$parameters, "b"))
                eta2[[nx[j]]] <- eta2[[nx[j]]] + p.state$fitted.values
                lpost1 <- family$loglik(y, family$map2par(eta2)) + lp + x[[nx[j]]]$smooth.construct[[sj]]$prior(p.state$parameters)
                if(lpost1 < lpost0) {
                  warning(paste("logPost is decreasing updating term: ", nx[j], ", ",
                    x[[nx[j]]]$smooth.construct[[sj]]$label, "; diff: ", lpost1 - lpost0, sep = ""))
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

      if(verbose) {
        cat(if(ia) "\r" else if(iter > 1) "\n" else NULL)
        vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
          " logPost ", fmt(family$loglik(y, peta) + get.log.prior(x), width = 8, digits = digits),
          " logLik ", fmt(family$loglik(y, peta), width = 8, digits = digits),
          " edf ", fmt(edf, width = 6, digits = digits),
          " eps ", fmt(eps0, width = 6, digits = digits + 2),
          " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
        cat(vtxt)

        if(.Platform$OS.type != "unix" & ia) flush.console()
      }
    }

    elapsed <- c(proc.time() - ptm)[3]

    IC <- get.ic(family, y, peta, edf, nobs, criterion)
    logLik <- family$loglik(y, peta)
    logPost <- as.numeric(logLik + get.log.prior(x))

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
  pen <- switch(type,
    "AIC" = -2 * ll + 2 * edf,
    "BIC" = -2 * ll + edf * log(n),
    "AICc" = -2 * ll + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1),
    "MP" = -1 * (ll + edf)
  )
  return(pen)
}

get.ic2 <- function(logLik, edf, n, type = c("AIC", "BIC", "AICc", "MP"), ...)
{
  type <- match.arg(type)
  pen <- switch(type,
    "AIC" = -2 * logLik + 2 * edf,
    "BIC" = -2 * logLik + edf * log(n),
    "AICc" = -2 * logLik + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1),
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
tau2.optim <- function(f, start, ..., scale = 10, eps = 0.0001, maxit = 1)
{
  foo <- function(par, start, k) {
    start[k] <- cround(par)
    return(f(start, ...))
  }

  start <- cround(start)
  ic0 <- f(start)

  iter <- 0; eps0 <- eps + 1
  while((eps0 > eps) & (iter < maxit)) {
    start0 <- start
    for(k in seq_along(start)) {
      xr <- c(start[k] / scale, start[k] * scale + 1)
      tpar <- try(optimize(foo, interval = xr, start = start, k = k), silent = TRUE)
      if(!inherits(tpar, "try-error")) {
        if(tpar$objective < ic0) {
          start[k] <- cround(tpar$minimum)
          ic0 <- tpar$objective
        }
      }
    }
    if(length(start) < 2)
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

  if(is.null(args$hess)) {
    ## Compute weights.
    if(no_ff) {
      hess <- process.derivs(family$hess[[id]](y, peta, id = id, ...), is.weight = TRUE)
    } else {
      hess <- ffdf_eval_sh(y, peta, FUN = function(y, par) {
        process.derivs(family$hess[[id]](y, par, id = id), is.weight = TRUE)
      })
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
      P <- matrix_inv(XWX, index = x$sparse.setup, all_diagonal = x$all_diagonal)
    } else {
      S <- 0
      tau2 <- get.state(x, "tau2")
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](c(g0, x$fixed.hyper)) else x$S[[j]]
      P <- matrix_inv(XWX + S, index = x$sparse.setup, all_diagonal = x$all_diagonal)
    }
    x$state$parameters <- set.par(x$state$parameters, drop(P %*% crossprod(x$X, x$rres)), "b")
  } else {
    args <- list(...)
    edf0 <- args$edf - x$state$edf
    eta2 <- eta

    env <- new.env()

    objfun <- function(tau2, ...) {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](c(g0, x$fixed.hyper)) else x$S[[j]]
      P <- matrix_inv(XWX + S, index = x$sparse.setup, all_diagonal = x$all_diagonal)
      if(inherits(P, "try-error")) return(NA)
      g <- drop(P %*% crossprod(x$X, x$rres))
      if(any(is.na(g)) | any(g %in% c(-Inf, Inf))) g <- rep(0, length(g))
      fit <- x$fit.fun(x$X, g)
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

    assign("ic00_val", objfun(get.state(x, "tau2")), envir = env)

    tau2 <- tau2.optim(objfun, start = get.state(x, "tau2"))

    if(!is.null(env$state))
      return(env$state)

    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](c(x$state$parameters, x$fixed.hyper)) else x$S[[j]]
    P <- matrix_inv(XWX + S, index = x$sparse.setup, all_diagonal = x$all_diagonal)
    x$state$parameters <- set.par(x$state$parameters, drop(P %*% crossprod(x$X, x$rres)), "b")
  }

  ## Compute fitted values.
  g <- get.state(x, "b")
  if(any(is.na(g)) | any(g %in% c(-Inf, Inf))) {
    x$state$parameters <- set.par(x$state$parameters, rep(0, length(get.state(x, "b"))), "b")
  }
  x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))
  x$state$edf <- sum_diag(XWX %*% P)

  if(!is.null(x$prior)) {
    if(is.function(x$prior))
      x$state$log.prior <- x$prior(x$state$parameters)
  }

  return(x$state)
}


bfit_iwls_spam <- function(x, family, y, eta, id, weights, criterion, ...)
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
  XWX <- (t(diag.spam(x$weights) %*% x$X)) %*% x$X
  Xr <- crossprod.spam(x$X, x$rres)

  if(!x$state$do.optim | x$fixed | x$fxsp) {
    if(!x$fixed) {
      tau2 <- get.state(x, "tau2")
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      U <- update.spam.chol.NgPeyton(x$sparse.setup$spam.cholFactor, XWX + S)
    } else {
      U <- update.spam.chol.NgPeyton(x$sparse.setup$spam.cholFactor, XWX)
    }
    P <- chol2inv.spam(U)
    b <- P %*% Xr
    x$state$parameters <- set.par(x$state$parameters, b, "b")
  } else {
    args <- list(...)
    edf0 <- args$edf - x$state$edf
    eta2 <- eta

    env <- new.env()

    objfun <- function(tau2, ...) {
      S <- 0
      for(j in seq_along(x$S))
        S <- S + 1 / tau2[j] * x$S[[j]]
      U <- update.spam.chol.NgPeyton(x$sparse.setup$spam.cholFactor, XWX + S)
      P <- chol2inv.spam(U)
      b <- P %*% Xr
      fit <- x$fit.fun(x$X, b)
      edf <- sum_diag(XWX %*% P)
      eta2[[id]] <- eta2[[id]] + fit
      ic <- get.ic(family, y, family$map2par(eta2), edf0 + edf, length(z), criterion, ...)
      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          par <- c(b, tau2)
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

    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S[[j]]
    U <- update.spam.chol.NgPeyton(x$sparse.setup$spam.cholFactor, XWX + S)
    P <- chol2inv.spam(U)
    b <- P %*% Xr
    x$state$parameters <- set.par(x$state$parameters, b, "b")
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
    method = "BFGS", control = list(), tau2 = get.par(tpar, "tau2"), hessian = TRUE),
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
      x[[j]]$smooth.construct[[sj]]$state$parameters <- set.par(x[[j]]$smooth.construct[[sj]]$state$parameters, tpar, "b")
      x[[j]]$smooth.construct[[sj]]$state$fitted.values <- x[[j]]$smooth.construct[[sj]]$fit.fun(x[[j]]$smooth.construct[[sj]]$X,
        get.par(tpar, "b"))
      eta[[j]] <- eta[[j]] + fitted(x[[j]]$smooth.construct[[sj]]$state)
      lprior <- lprior + x[[j]]$smooth.construct[[sj]]$prior(c(tpar, get.state(x[[j]]$smooth.construct[[sj]], "tau2"), x[[j]]$smooth.construct[[sj]]$fixed.hyper))
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


## Gradient vecor of the log-posterior.
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
opt <- function(x, y, family, start = NULL, verbose = TRUE, digits = 3,
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
      as.integer(oind))
  }
  invisible(NULL)
}


## Likelihood based boosting.
#boost_logLik <- function(x, y, family, weights = NULL, offset = NULL,
#  criterion = c("AICc", "BIC", "AIC"),
#  nu = 1, df = 4, maxit = 100, mstop = NULL, best = TRUE,
#  verbose = TRUE, digits = 4,
#  eps = .Machine$double.eps^0.25, plot = TRUE, ...)
#{
#  if(!is.null(mstop))
#    maxit <- mstop

#  if(is.null(attr(x, "bamlss.engine.setup")))
#    x <- bamlss.engine.setup(x, ...)

#  nx <- family$names
#  if(!all(nx %in% names(x)))
#    stop("parameter names mismatch with family names!")
#  criterion <- match.arg(criterion)

#  np <- length(nx)
#  y <- y[[1]]
#  nobs <- if(is.null(dim(y))) length(y) else nrow(y)

#  ## Setup boosting structure, i.e, all parametric
#  ## terms get an entry in $smooth.construct object.
#  ## Intercepts are initalized.
#  x <- boost.transform(x, y, df, family, weights, offset, maxit, eps, ...)

#  ## Create a list() that saves the states for
#  ## all parameters and model terms.
#  states <- make.state.list(x)

#  ## Term selector help vectors.
#  select <- unlist(states)
#  term.names <- names(select)

#  ## Extract actual predictor.
#  eta <- get.eta(x)

#  ## Initial parameters.
#  parameters <- get.all.par(x)

#  ## Start boosting.
#  eps0 <- 1; iter <- 1
#  save.ic <- save.ll <- NULL
#  ll <- family$loglik(y, family$map2par(eta))
#  while(iter <= maxit) {
#    eta0 <- eta

#    ## Cycle through all parameters
#    for(i in nx) {
#      peta <- family$map2par(eta)

#      ## Compute weights.
#      hess <- family$hess[[i]](y, peta, id = i)

#      ## Score.
#      score <- family$score[[i]](y, peta, id = i)

#      ## Compute working observations.
#      z <- eta[[i]] + 1 / hess * score

#      ## Residuals.
#      resids <- z - eta[[i]]

#      for(j in names(x[[i]]$smooth.construct)) {
#        ## Get updated parameters.
#        states[[i]][[j]] <- boost_iwls(x[[i]]$smooth.construct[[j]], hess, resids, nu)

#        ## Compute likelihood contribution.
#        eta[[i]] <- eta[[i]] + fitted(states[[i]][[j]])
#        select[paste(i, j, sep = ".")] <- -1 * (ll - family$loglik(y, family$map2par(eta)))
#        eta[[i]] <- eta0[[i]]
#      }
#    }

#    ## Which term to update.
#    take <- strsplit(term.names[which.max(select)], ".", fixed = TRUE)[[1]]

#    ## Update selected base learner.
#    eta[[take[1]]] <- eta[[take[1]]] + fitted(states[[take[1]]][[take[2]]])

#    ## Write to x.
#    x[[take[1]]]$smooth.construct[[take[2]]]$state <- increase(x[[take[1]]]$smooth.construct[[take[2]]]$state, states[[take[1]]][[take[2]]])
#    x[[take[1]]]$smooth.construct[[take[2]]]$selected[iter] <- 1
#    x[[take[1]]]$smooth.construct[[take[2]]]$loglik[iter] <- max(select)

#    edf <- get.edf(x, type = 2)

#    eps0 <- do.call("cbind", eta)
#    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
#    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

#    peta <- family$map2par(eta)
#    IC <- get.ic(family, y, peta, edf, nobs, criterion)
#    if(!is.null(save.ic) & best) {
#      if(all(IC < save.ic))
#        parameters <- get.all.par(x)
#    }
#    ll <- family$loglik(y, peta)

#    save.ic <- c(save.ic, IC)
#    save.ll <- c(save.ll, ll)

#    if(verbose) {
#      cat(if(interactive()) "\r" else "\n")
#      vtxt <- paste(criterion, " ", fmt(IC, width = 8, digits = digits),
#        " logLik ", fmt(ll, width = 8, digits = digits),
#        " edf ", fmt(edf, width = 6, digits = digits),
#        " eps ", fmt(eps0, width = 6, digits = digits + 2),
#        " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
#      cat(vtxt)

#      if(.Platform$OS.type != "unix" & interactive()) flush.console()
#    }

#    iter <- iter + 1
#  }

#  if(verbose) cat("\n")

#  mstop <- which.min(save.ic)

#  ## Overwrite parameter state.
#  if(best) {
#    for(i in nx) {
#      rn <- NULL
#      for(j in names(x[[i]]$smooth.construct)) {
#        x[[i]]$smooth.construct[[j]]$state$parameters <- parameters[[i]]$s[[j]]
#        g <- get.par(x[[i]]$smooth.construct[[j]]$state$parameters, "b")
#        if(!is.null(tau2 <- attr(parameters[[i]]$s[[j]], "true.tau2"))) {
#          x[[i]]$smooth.construct[[j]]$state$parameters <- set.par(x[[i]]$smooth.construct[[j]]$state$parameters,
#            tau2, "tau2")
#        }
#        x[[i]]$smooth.construct[[j]]$state$edf <- attr(parameters[[i]]$s[[j]], "edf")
#        if(is.null(x[[i]]$smooth.construct[[j]]$state$edf))
#          x[[i]]$smooth.construct[[j]]$state$edf <- 0
#      }
#    }
#  }

#  if(verbose) {
#    cat("---\n", criterion, "=", save.ic[mstop], "-> at mstop =", mstop, "\n---\n")
#  }

#  bsum <- make.boost.summary(x, mstop, criterion, save.ic)
#  x <- boost.retransform(x)

#  list("parameters" = get.all.par(x), "fitted.values" = get.eta(x), "boost.summary" = bsum)
#}


## Gradient boosting.
boost <- function(x, y, family,
  nu = 0.1, df = 4, maxit = 400, mstop = NULL,
  verbose = TRUE, digits = 4, flush = TRUE,
  eps = .Machine$double.eps^0.25, nback = NULL, plot = TRUE,
  initialize = TRUE, ...)
{
  ## FIXME: hard coded.
  weights <- offset <- NULL

  nx <- family$names
  if(!all(nx %in% names(x)))
    stop("parameter names mismatch with family names!")

  if(!is.null(mstop))
    maxit <- mstop

  if(!is.null(nback)) {
    if(is.null(maxit))
      maxit <- 10000
  }

  if(is.null(maxit))
    stop("please set either argument 'maxit' or 'mstop'!")

  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, df = df, ...)

  np <- length(nx)
  nobs <- nrow(y)
  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }

  ## Setup boosting structure, i.e, all parametric
  ## terms get an entry in $smooth.construct object.
  ## Intercepts are initalized.
  x <- boost.transform(x = x, y = y, df = NULL, family = family,
    maxit = maxit, eps = eps, initialize = initialize, ...)

  ## Create a list() that saves the states for
  ## all parameters and model terms.
  states <- make.state.list(x)

  ## Matrix of all parameters.
  parm <- make.par.list(x, iter = maxit)

  ## Term selector help vectors.
  select <- rep(NA, length = length(nx))
  names(select) <- nx
  loglik <- select

  ## Save rss in list().
  rss <- make.state.list(x, type = 2)

  ## Extract actual predictor.
  eta <- get.eta(x)

  ## Print stuff.
  ia <- if(flush) interactive() else FALSE

  ## Env for C.
  rho <- new.env()

  ## Start boosting.
  eps0 <- 1; iter <- if(initialize) 2 else 1
  save.ll <- NULL
  ll <- family$loglik(y, family$map2par(eta))
  ptm <- proc.time()
  while(iter <= maxit) {
    eta0 <- eta

    ## Cycle through all parameters
    for(i in nx) {
      peta <- family$map2par(eta)

      ## Actual gradient.
      grad <- process.derivs(family$score[[i]](y, peta, id = i), is.weight = FALSE)

      ## Fit to gradient.
      for(j in names(x[[i]]$smooth.construct)) {
        ## Get updated parameters.
        states[[i]][[j]] <- if(is.null(x[[i]]$smooth.construct[[j]][["boost.fit"]])) {
          .Call("boost_fit", x[[i]]$smooth.construct[[j]], grad, nu, rho, PACKAGE = "bamlss")
        } else {
          x[[i]]$smooth.construct[[j]][["boost.fit"]](x[[i]]$smooth.construct[[j]], grad, nu, rho)
        }

        ## Get rss.
        rss[[i]][j] <- states[[i]][[j]]$rss
      }

      ## Which one is best?
      select[i] <- which.min(rss[[i]])

      ## Compute likelihood contribution.
      eta[[i]] <- eta[[i]] + fitted(states[[i]][[select[i]]])
      llf <- family$loglik(y, family$map2par(eta))
      loglik[i] <- -1 * (ll - llf)

      if(loglik[i] > 40000) {
        print(as.character(c(ll, llf)))
        stop()
      }

      eta[[i]] <- eta0[[i]]
    }

    i <- which.max(loglik)

    ## Which term to update.
    take <- c(nx[i], names(rss[[i]])[select[i]])

    ## Update selected base learner.
    eta[[take[1]]] <- eta[[take[1]]] + states[[take[1]]][[take[2]]]$fitted.values

    ## Write to x.
    x[[take[1]]]$smooth.construct[[take[2]]]$state <- increase(x[[take[1]]]$smooth.construct[[take[2]]]$state, states[[take[1]]][[take[2]]])
    x[[take[1]]]$smooth.construct[[take[2]]]$selected[iter] <- 1
    x[[take[1]]]$smooth.construct[[take[2]]]$loglik[iter] <- loglik[i]

    ## Save parameters.
    parm[[take[1]]][[take[2]]][iter, ] <- get.par(states[[take[1]]][[take[2]]]$parameters, "b")
    eps0 <- do.call("cbind", eta)
    eps0 <- mean(abs((eps0 - do.call("cbind", eta0)) / eps0), na.rm = TRUE)
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps + 1

    peta <- family$map2par(eta)
    ll <- family$loglik(y, peta)

    save.ll <- c(save.ll, ll)

    if(verbose) {
      cat(if(ia) "\r" else "\n")
      vtxt <- paste(
        " logLik ", fmt(ll, width = 8, digits = digits),
        " eps ", fmt(eps0, width = 6, digits = digits + 2),
        " iteration ", formatC(iter, width = nchar(maxit)), sep = "")
      cat(vtxt)

      if(.Platform$OS.type != "unix" & ia) flush.console()
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

  bsum <- make.boost.summary(x, if(is.null(nback)) maxit else (iter - 1), save.ll)
  if(plot)
    plot.boost.summary(bsum)

  return(list("parameters" = parm2mat(parm, if(is.null(nback)) maxit else (iter - 1)),
    "fitted.values" = get.eta(x), "nobs" = nobs, "boost.summary" = bsum, "runtime" = elapsed))
}


## Boost setup.
boost.transform <- function(x, y, df = NULL, family,
  weights = NULL, offset = NULL, maxit = 100,
  eps = .Machine$double.eps^0.25, initialize = TRUE, ...)
{
  np <- length(x)
  nx <- names(x)

  ## Initialize select indicator and intercepts.
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      if(!is.null(df))
        x[[nx[j]]]$smooth.construct[[sj]] <- assign.df(x[[nx[j]]]$smooth.construct[[sj]], df)
      if(!x[[nx[j]]]$smooth.construct[[sj]]$fxsp & !x[[nx[j]]]$smooth.construct[[sj]]$fixed) {
        x[[nx[j]]]$smooth.construct[[sj]]$old.optimize <- x[[nx[j]]]$smooth.construct[[sj]]$state$do.optim
        x[[nx[j]]]$smooth.construct[[sj]]$state$do.optim <- FALSE
        x[[nx[j]]]$smooth.construct[[sj]]$do.optim <- FALSE
      }
    }
    if(has_pterms(x[[nx[j]]]$terms)) {
      ii <- which(names(x[[nx[j]]]$smooth.construct) == "model.matrix")
      model.matrix <- list()
      cn <- colnames(x[[nx[j]]]$smooth.construct[[ii]]$X)
      g0 <- get.par(x[[nx[j]]]$smooth.construct[[ii]]$state$parameters, "b")
      nm <- NULL
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
        class(model.matrix[[pj]]) <- class(x[[nx[j]]]$smooth.construct[[ii]])
      }
      names(model.matrix) <- cn
      x[[nx[j]]]$smooth.construct[[ii]] <- NULL
      x[[nx[j]]]$smooth.construct <- c(model.matrix, x[[nx[j]]]$smooth.construct)
    }
  }

  ## Save more info.
  for(j in 1:np) {
    for(sj in seq_along(x[[nx[j]]]$smooth.construct)) {
      x[[nx[j]]]$smooth.construct[[sj]]$selected <- rep(0, length = maxit)
      x[[nx[j]]]$smooth.construct[[sj]]$loglik <- rep(0, length = maxit)
      x[[nx[j]]]$smooth.construct[[sj]]$state$edf <- 0
      x[[nx[j]]]$smooth.construct[[sj]]$state$rss <- 0
      nc <- ncol(x[[nx[j]]]$smooth.construct[[sj]]$X)
      nr <- nrow(x[[nx[j]]]$smooth.construct[[sj]]$X)
      x[[nx[j]]]$smooth.construct[[sj]]$XWX <- matrix(0, nc, nc)
      x[[nx[j]]]$smooth.construct[[sj]]$XW <- matrix(0, nc, nr)
      if(!is.null(x[[nx[j]]]$smooth.construct[[sj]]$S))
        x[[nx[j]]]$smooth.construct[[sj]]$penaltyFunction <- as.integer(sapply(x[[nx[j]]]$smooth.construct[[sj]]$S, is.function))
      else
        x[[nx[j]]]$smooth.construct[[sj]]$penaltyFunction <- 0L
    }
  }

  if(initialize) {
    eta <- get.eta(x)
    eta <- init.eta(eta, y, family, nobs)
    nobs <- length(eta[[1]])
    start <- unlist(lapply(eta, mean, na.rm = TRUE))

    objfun <- function(par) {
      eta <- list()
      for(i in seq_along(nx))
        eta[[nx[i]]] <- rep(par[i], length = nobs)
      ll <- family$loglik(y, family$map2par(eta))
      return(ll)
    }

    gradfun <- function(par) {
      eta <- list()
      for(i in seq_along(nx))
        eta[[nx[i]]] <- rep(par[i], length = nobs)
      peta <- family$map2par(eta)
      grad <- par
      for(j in nx) {
        score <- process.derivs(family$score[[j]](y, peta, id = j), is.weight = FALSE)
        grad[i] <- mean(score)
      }
      return(grad)
    }

    opt <- optim(start, fn = objfun, gr = gradfun, method = "BFGS", control = list(fnscale = -1))

    for(i in nx) {
      if(!is.null(x[[i]]$smooth.construct[["(Intercept)"]])) {
        x[[i]]$smooth.construct[["(Intercept)"]]$state$parameters[1] <- opt$par[i]
        x[[i]]$smooth.construct[["(Intercept)"]]$state$fitted.values <- rep(opt$par[i], length = nobs)
      }
    }
  }

  return(x)
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
        rval[[j]] <- matrix(0, nrow = iter, ncol = ncol(x$smooth.construct[[j]]$X))
        colnames(rval[[j]]) <- names(get.par(x$smooth.construct[[j]]$state$parameters, "b"))
        rval[[j]][1, ] <- get.par(x$smooth.construct[[j]]$state$parameters, "b")
        if(!is.null(x$smooth.construct[[j]]$is.model.matrix))
          attr(rval[[j]], "is.model.matrix") <- TRUE
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

parm2mat <- function(x, mstop)
{
  nx <- names(x)
  for(i in seq_along(x)) {
    is.mm <- NULL
    for(j in names(x[[i]])) {
      if(!is.null(attr(x[[i]][[j]], "is.model.matrix")))
        is.mm <- c(is.mm, j)
      cn <- colnames(x[[i]][[j]])
      x[[i]][[j]] <- apply(x[[i]][[j]][1:mstop, , drop = FALSE], 2, cumsum)
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
  return(x)
}


#	for(i in seq_along(coefs))
#	{
#		curname<-names(coefs)[i]
#		cpos<-match(curname, colnames(dfr))
#		usedScale<-attr(dfr[[cpos]], "scaled:scale")
#		usedCenter<-attr(dfr[[cpos]], "scaled:center")
#		if(! is.null(usedScale))
#		{
#			catwif(verbosity > 0, "Scaling back for variable", curname)
#			catwif(verbosity >1, "usedScale structure")
#			if(verbosity > 1) str(usedScale)
#			catwif(verbosity >1, "usedCenter structure")
#			if(verbosity > 1) str(usedCenter)
#			oldcoef<-coefs[i]
#			itc<-itc - ((oldcoef * usedCenter)/usedScale)
#			coefs[i]<-oldcoef / usedScale
#		}
#	}
#	coefs<-c(itc, coefs)
#	names(coefs)[1]<-itcname
#	return(coefs)


#          if(!is.null(x[[i]]$smooth.construct[[j]]$boost.scale)) {
#            mx <- x[[i]]$smooth.construct[[j]]$boost.scale$mean
#            sdx <- x[[i]]$smooth.construct[[j]]$boost.scale$sd
#            x[[i]]$smooth.construct[[j]]$X <- x[[i]]$smooth.construct[[j]]$X * sdx + mx
#            if(!is.null(intercept))
#              intercept <- intercept - ((b * mx) / sdx)
#            b <- b / sdx
#          }


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
boost_fit <- function(x, y, nu)
{
  ## Compute reduced residuals.
  xbin.fun(x$binning$sorted.index, rep(1, length = length(y)), y, x$weights, x$rres, x$binning$order)

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
  x$state$rss <- sum((x$state$fitted.values - y)^2)

  return(x$state)
}


## Increase coefficients.
increase <- function(state0, state1)
{
  g <- get.par(state0$parameters, "b") + get.par(state1$parameters, "b")
  state0$fitted.values <- fitted(state0) + fitted(state1)
  state0$parameters <- set.par(state0$parameters, g, "b")
  state0$edf <- state1$edf
  attr(state0$parameters, "true.tau2") <- attr(state1$parameters, "true.tau2")
  attr(state0$parameters, "edf") <- attr(state1$parameters, "edf")
  state0
}


## Extract summary for boosting.
make.boost.summary <- function(x, mstop, save.ic)
{
  nx <- names(x)
  labels <- NULL
  ll.contrib <- NULL
  bsum <- lmat <- list()
  for(i in nx) {
    rn <- NULL
    for(j in names(x[[i]]$smooth.construct)) {
      labels <- c(labels, paste(x[[i]]$smooth.construct[[j]]$label, i, sep = "."))
      rn <- c(rn, x[[i]]$smooth.construct[[j]]$label)
      bsum[[i]] <- rbind(bsum[[i]], sum(x[[i]]$smooth.construct[[j]]$selected[1:mstop]) / mstop * 100)
      lmat[[i]] <- rbind(lmat[[i]], sum(x[[i]]$smooth.construct[[j]]$loglik[1:mstop]))
      ll.contrib <- cbind(ll.contrib, cumsum(x[[i]]$smooth.construct[[j]]$loglik[1:mstop]))
    }
    if(!is.matrix(bsum[[i]])) bsum[[i]] <- matrix(bsum[[i]], nrow = 1)
    bsum[[i]] <- cbind(bsum[[i]], lmat[[i]])
    if(!is.matrix(bsum[[i]])) bsum[[i]] <- matrix(bsum[[i]], nrow = 1)
    colnames(bsum[[i]]) <- c(paste(i, "% selected"), "LogLik contrib.")
    rownames(bsum[[i]]) <- rownames(lmat[[i]]) <- rn
    bsum[[i]] <- bsum[[i]][order(bsum[[i]][, 2], decreasing = TRUE), , drop = FALSE]
  }
  colnames(ll.contrib) <- labels
  names(bsum) <- nx
  bsum <- list("summary" = bsum, "mstop" = mstop,
    "ic" = save.ic[1:mstop], "loglik" = ll.contrib)
  class(bsum) <- "boost.summary"
  return(bsum)
}


boost.summary <- function(object, ...)
{
  if(!is.null(object$model.stats$optimizer$boost.summary))
    print.boost.summary(object$model.stats$optimizer$boost.summary, ...)
  invisible(object$model.stats$optimizer$boost.summary)
}


## Smallish print function for boost summaries.
print.boost.summary <- function(x, summary = TRUE, plot = TRUE,
  which = c("loglik", "loglik.contrib"), intercept = TRUE,
  spar = TRUE, ...)
{
  if(inherits(x, "bamlss"))
    x <- x$model.stats$optimizer$boost.summary
  if(is.null(x))
    stop("no summary for boosted model available")
  if(summary) {
    np <- length(x$summary)
    cat("\n")
    cat("logLik. =", x$ic[x$mstop], "-> at mstop =", x$mstop, "\n---\n")
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
        xn <- sapply(strsplit(colnames(x$loglik), ".", fixed = TRUE), function(x) { x[length(x)] })
        cols <- rainbow_hcl(length(unique(xn)))
        matplot(x$loglik, type = "l", lty = 1,
          xlab = "Iteration", ylab = "LogLik contribution", col = cols[as.factor(xn)], ...)
        abline(v = x$mstop, lwd = 3, col = "lightgray")
        axis(4, at = x$loglik[nrow(x$loglik), ], labels = colnames(x$loglik), las = 1)
        axis(3, at = x$mstop, labels = paste("mstop =", x$mstop))
      }
    }
  }

  return(invisible(x))
}


plot.boost.summary <- function(x, ...)
{
  print.boost.summary(x, summary = FALSE, plot = TRUE, ...) 
}

boost.plot <- function(x, which = c("loglik", "loglik.contrib", "parameters"),
  intercept = TRUE, spar = TRUE, mstop = NULL, name = NULL, labels = NULL, color = NULL, ...)
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
    mstop <- x$model.stats$optimizer$boost.summary$mstop
  x$model.stats$optimizer$boost.summary$mstop <- mstop
  x$model.stats$optimizer$boost.summary$ic <- x$model.stats$optimizer$boost.summary$ic[1:mstop]
  x$model.stats$optimizer$boost.summary$loglik <- x$model.stats$optimizer$boost.summary$loglik[1:mstop, , drop = FALSE]

  for(w in which) {
    if(w %in% c("loglik", "loglik.contrib")) {
      if((w == "loglik") & spar)
        par(mar = c(5.1, 4.1, 2.1, 2.1))
      if((w == "loglik.contrib") & spar)
        par(mar = c(5.1, 4.1, 2.1, 10.1))
      plot.boost.summary(x, which = w, spar = FALSE, intercept = intercept, ...)
    }
    if(w == "parameters") {
      if(spar)
        par(mar = c(5.1, 4.1, 2.1, 10.1))
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
        labs <- colnames(x$parameters)
        if(!is.null(name)) {
          for(j in seq_along(name))
            labs <- gsub(name[j], "", labs, fixed = TRUE)
        }
      } else labs <- rep(labels, length.out = ncol(x$parameters))

      matplot(p, type = "l", lty = 1, col = cols[as.factor(xn)], xlab = "Iteration",
        ylab = "Value", ...)
      abline(v = mstop, lwd = 3, col = "lightgray")
      axis(4, at = p[nrow(p), ], labels = labs, las = 1)
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
              x[[id]]$smooth.construct$model.matrix$state$parameters <- tpar
              x[[id]]$smooth.construct$model.matrix$state$fitted.values <- x[[id]]$smooth.construct$model.matrix$fit.fun(x[[id]]$smooth.construct$model.matrix$X, x[[id]]$smooth.construct$model.matrix$state$parameters)
            }
          }
        }
        for(j in seq_along(x[[id]]$smooth.construct)) {
          tl <- x[[id]]$smooth.construct[[j]]$label
          take <- grep(tl <- paste(id, "s", tl, sep = "."),
            nstart[tns %in% id], fixed = TRUE, value = TRUE)
          if(x[[id]]$smooth.construct[[j]]$by == "NA") {
            take <- take[!grepl(paste(tl, ":", sep = ""), take, fixed = TRUE)]
          }
          if(length(take)) {
            tpar <- start[take]
            i <- grep2(c(".edf", ".accepted", ".alpha"), names(tpar), fixed = TRUE)
            tpar <- if(length(i)) tpar[-i] else tpar
            names(tpar) <- gsub(paste(tl, ".", sep = ""), "", names(tpar), fixed = TRUE)
            spar <- x[[id]]$smooth.construct[[j]]$state$parameters
            if(length(get.par(tpar, "b")))
              spar <- set.par(spar, get.par(tpar, "b"), "b")
            if(any(grepl("tau2", names(tpar)))) {
              spar <- set.par(spar, get.par(tpar, "tau2"), "tau2")
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


lasso <- function(x, y, start = NULL, adaptive = TRUE,
  lower = 0.001, upper = 1000,  nlambda = 100, lambda = NULL,
  verbose = TRUE, digits = 4, flush = TRUE,
  nu = NULL, stop.nu = NULL, ridge = .Machine$double.eps^0.25,
  zeromodel = NULL, ...)
{
  method <- list(...)$method
  if(is.null(method))
    method <- 1

  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, update = bfit_iwls, ...)

  start2 <- start

  lambdas <- if(is.null(lambda)) {
    exp(seq(log(upper), log(lower), length = nlambda))
  } else lambda

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
        zeromodel <- bfit(x = x, y = y, start = start, verbose = verbose[1], nu = nu[2], stop.nu = stop.nu[2], ...)
      } else {
        zeromodel <- opt(x = x, y = y, start = start, verbose = verbose[1], ...)
      }
    }
    x <- lasso.transform(x, zeromodel, nobs = nrow(y))
  }

  for(l in seq_along(lambdas)) {
    if(l > 1)
      start <- unlist(par[[l - 1]])
    tau2 <- NULL
    for(i in names(x)) {
      for(j in names(x[[i]]$smooth.construct)) {
        if(inherits(x[[i]]$smooth.construct[[j]], "lasso.smooth")) {
          tau2 <- get.par(x[[i]]$smooth.construct[[j]]$state$parameters, "tau2")
          nt <- names(tau2)
          tau2 <- rep(1 / lambdas[l], length.out = length(tau2))
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
      b <- bfit(x = x, y = y, start = start, verbose = verbose[2], nu = nu[2], stop.nu = stop.nu[2], ...)
    } else {
      b <- opt(x = x, y = y, start = start, verbose = verbose[2], ...)
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
        " lambda ", fmt(lambdas[l], width = 6, digits = digits),
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

  ic <- cbind(ic, "lambda" = lambdas)
  rownames(ic) <- NULL
  class(ic) <- c("lasso.stats", "matrix")

  list("parameters" = do.call("rbind", par), "lasso.stats" = ic, "nobs" = nrow(y))
}

lasso.transform <- function(x, zeromodel, nobs = NULL, ...)
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
          names(w) <- paste("lasso", 1:length(w), sep = "")
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
  ls <- attr(lasso.stop(x), "stats")
  ic <- grep("ic", names(ls), ignore.case = TRUE, value = TRUE)
  cat(ic, "=", ls[ic], "-> at lambda =", ls["lambda"], "\n")
  ls <- ls[!(names(ls) %in% c("lambda", ic))]
  ls <- paste(names(ls), "=", round(ls, digits = digits), collapse = " ")
  cat(ls, "\n---\n")
  return(invisible(NULL))
}


lasso.coef <- function(x, ...) {
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


lasso.plot <- function(x, which = c("criterion", "parameters"), spar = TRUE, name = NULL,
  mstop = NULL, retrans = FALSE, color = NULL, show.lambda = TRUE, labels = NULL,
  digits = 2, ...)
{
  if(!is.character(which)) {
    which <- c("criterion", "parameters")[as.integer(which)]
  } else {
    which <- tolower(which)
    which <- match.arg(which, several.ok = TRUE)
  }
  if(is.null(mstop))
    mstop <- 1:nrow(x$parameters)
  if(retrans)
    x$parameters <- lasso.coef(x, mstop = mstop)
  par <- x$parameters[mstop, , drop = FALSE]
  npar <- colnames(par)
  for(j in c("Intercept", ".edf", ".lambda", ".tau"))
    npar <- npar[!grepl(j, npar, fixed = TRUE)]
  x$parameters <- x$parameters[mstop, npar, drop = FALSE]
  ic <- x$model.stats$optimizer$lasso.stats
  log_lambda <- log(ic[, "lambda"])
  nic <- grep("ic", colnames(ic), value = TRUE, ignore.case = TRUE)

  if(spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = c(1, length(which)), mar = c(5.1, 5.1, 4.1, 1.1))
  }

  at <- pretty(1:nrow(ic))
  at[at == 0] <- 1

  fmt2 <- Vectorize(fmt)

  if("criterion" %in% which) {
    plot(ic[mstop, nic], type = "l",
      xlab = expression(log(lambda)), ylab = nic, axes = FALSE)
    at <- pretty(mstop)
    at[at == 0] <- 1
    axis(1, at = at, labels = fmt2(log_lambda[mstop][at], digits))
    axis(2)
    if(show.lambda) {
      i <- which.min(ic[, nic])
      abline(v = i, col = "lightgray", lwd = 2, lty = 2)
      val <- round(ic[i, "lambda"], 4)
      axis(3, at = i, labels = substitute(paste(lambda, '=', val)))
    }
    box()
  }

  if("parameters" %in% which) {
    if(spar)
      par(mar = c(5.1, 5.1, 4.1, 10.1))

    if(!is.null(name)) {
      x$parameters <- x$parameters[, grep2(name, colnames(x$parameters), fixed = TRUE), drop = FALSE]
    }
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

    matplot(x$parameters, type = "l", lty = 1, col = cols[as.factor(xn)],
      xlab = expression(log(lambda)), ylab = expression(beta[j]), axes = FALSE, ...)
    if(is.null(labels)) {
      labs <- colnames(x$parameters)
      if(!is.null(name)) {
        for(j in seq_along(name))
          labs <- gsub(name[j], "", labs, fixed = TRUE)
      }
    } else labs <- rep(labels, length.out = ncol(x$parameters))
    axis(4, at = x$parameters[nrow(x$parameters), ],
      labels = labs, las = 1)
    at <- pretty(mstop)
    at[at == 0] <- 1
    axis(1, at = at, labels = fmt2(log_lambda[mstop][at], digits))
    axis(2)
    if(show.lambda) {
      i <- which.min(ic[, nic])
      abline(v = i, col = "lightgray", lwd = 2, lty = 2)
      val <- round(ic[i, "lambda"], 4)
      axis(3, at = i, labels = substitute(paste(lambda, '=', val)))
    }
    box()
  }

  return(invisible(NULL))
}

lasso.stop <- function(x)
{
  if(!inherits(x, "lasso.stats"))
    x <- x$model.stats$optimizer$lasso.stats
  nic <- grep("ic", colnames(x), value = TRUE, ignore.case = TRUE)
  i <- which.min(x[, nic])
  attr(i, "stats") <- x[i, ]
  i
}

