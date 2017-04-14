## Create a 'bamlss.frame'.
bamlss.frame <- function(formula, data = NULL, family = "gaussian",
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit,
  contrasts = NULL, knots = NULL, specials = NULL, reference = NULL,
  model.matrix = TRUE, smooth.construct = TRUE, ytype = c("matrix", "vector", "integer"),
  scale.x = FALSE, scale.d = FALSE, ...)
{
  ## Parse family object.
  family <- bamlss.family(family, ...)

  ## Parse formula.
  if(!inherits(formula, "bamlss.formula"))
    formula <- bamlss.formula(formula, family, specials, env = parent.frame())
  if(!is.null(attr(formula, "orig.formula")))
    formula <- attr(formula, "orig.formula")

  ## Setup return object.
  bf <- list()
  bf$call <- match.call()

  ## Create the model frame.
  bf$model.frame <- bamlss.model.frame(formula, data, family, weights,
    subset, offset, na.action, specials, contrasts)

  if(!inherits(bf$model.frame, "ffdf")) {
    ## Type of y.
    ytype <- match.arg(ytype)

    ## Process categorical responses and assign 'y'.
    cf <- bamlss.formula.cat(formula, family, bf$model.frame, reference)
    if(!is.null(cf) & (ytype != "integer")) {
      rn <- response.name(formula, hierarchical = FALSE, na.rm = TRUE)
      hrn <- response.name(formula, hierarchical = TRUE, na.rm = TRUE)
      orig.formula <- formula
      formula <- cf$formula
      if(ytype == "matrix") {
        if(is.factor(bf$model.frame[[rn[1]]])) {
          f <- as.formula(paste("~ -1 +", rn[1]), env = NULL)
          bf$y <- bf$model.frame[rn]
          bf$y[rn[1]] <- model.matrix(f, data = bf$model.frame)
          colnames(bf$y[[rn[1]]]) <- rmf(gsub(rn[1], "", colnames(bf$y[[rn[1]]])))
          bf$y[[rn[1]]] <- bf$y[[rn[1]]][, rmf(c(names(formula), cf$reference))]
        } else {
          bf$y <- bf$model.frame[rn]
        }
      } else {
        bf$y <- bf$model.frame[rn]
        if(is.factor(bf$model.frame[[rn[1]]])) {
          if(ytype == "integer")
            bf$y[[1]] <- as.integer(bf$y[[1]]) - if(nlevels(bf$model.frame[[rn[1]]]) < 3) 1L else 0L
        }
      }
      if(length(hrn) > length(rn)) {
        ynot <- hrn[!(hrn %in% rn)]
        bf$y <- cbind(bf$y, bf$model.frame[ynot])
      }
      attr(bf$y, "reference") <- cf$reference
      family$names <- names(formula)
      family$links <- rep(family$links, length.out = length(formula))
      names(family$links) <- names(formula)
      attr(formula, "orig.formula") <- orig.formula
    } else {
      rn <- response.name(formula, hierarchical = FALSE, keep.functions = TRUE)
      rn <- rn[rn %in% names(bf$model.frame)]
      bf$y <- bf$model.frame[rn]
      for(j in rn) {
        if(is.factor(bf$y[[j]]) & (ytype == "matrix")) {
          f <- as.formula(paste("~ -1 +", j), env = NULL)
          bf$y[j] <- model.matrix(f, data = bf$model.frame)
        }
        if(is.factor(bf$y[[j]]) & (ytype == "integer")) {
          bf$y[[j]] <- as.integer(bf$y[[j]]) - if(nlevels(bf$y[[j]]) < 3) 1L else 0L
        }
      }
    }
  } else {
    rn <- response.name(formula, hierarchical = FALSE, keep.functions = TRUE)
    rn <- rn[rn %in% names(bf$model.frame)]
    bf$y <- bf$model.frame[rn]
  }

  bf$formula <- formula
  attr(bf$formula, "response.name") <- rn

  ## Add the terms object.
  bf$terms <- terms.bamlss.formula(formula, data = data, drop = FALSE, specials = specials, ...)

  ## Process possible score and hess functions.
  if(!is.null(score <- family$score)) {
    if(is.function(score))
      score <- list(score)
    family$score <- rep(score, length.out = length(formula))
    names(family$score) <- names(formula)
  }
  if(!is.null(hess <- family$hess)) {
    if(is.function(hess))
      hess <- list(hess)
    family$hess <- rep(hess, length.out = length(formula))
    names(family$hess) <- names(formula)
  }

  ## Add more functions to family object.
  bf$family <- complete.bamlss.family(family)

  if(inherits(bf$model.frame, "data.frame") & scale.d)
    bf$model.frame <- scale_model.frame(bf$model.frame, not = rn)

  ## Assign the 'x' master object.
  bf$x <- design.construct(bf$terms, data = bf$model.frame, knots = knots,
    model.matrix = model.matrix, smooth.construct = smooth.construct, model = NULL,
    scale.x = scale.x, specials = specials, ...)

  bf$knots <- knots

  ## Assign class and return.
  class(bf) <- c("bamlss.frame", "list")

  return(bf)
}


## Simple print method for 'bamlss.frame'
print.bamlss.frame <- function(x, ...)
{
  cat("'bamlss.frame' structure:", "\n")  
  nx <- c("call", "model.frame", "formula", "family", "terms", "x", "y", "knots")
  nx <- c(nx, names(x)[!(names(x) %in% nx)])
  for(i in nx) {
    if(!is.null(x[[i]])) {
      cat("  ..$", i, "\n")
      if(i == "x") {
        for(j in names(x[[i]])) {
          cat("  .. ..$", j, "\n")
          if(!all(c("formula", "fake.formula") %in% names(x[[i]][[j]]))) {
            for(k in names(x[[i]][[j]])) {
              cat("  .. .. ..$", k, "\n")
              for(d in names(x[[i]][[j]][[k]])) {
               cat("  .. .. .. ..$", d, "\n")
              }
            }
          } else {
            for(k in names(x[[i]][[j]]))
              cat("  .. .. ..$", k, "\n")
          }
        }
      }
      if(i == "y") {
        for(j in names(x[[i]])) {
          cat("  .. ..$", j, "\n")
        }
      }
    }
  }
  invisible(NULL)
}


## ff version for indexing.
match.index.ff <- function(x)
{
#  nodups <- ffwhich(x, !duplicated(x))
#  ind <- ffdfmatch(x, x[nodups, , drop = FALSE])
#  ord <- fforder(ind)
#  sindex <- ind[ord]
#  
#  return(list("match.index" = ind, "nodups" = nodups, "order" = ord, "sorted.index" = sindex, "uind" = ind[nodups]))
## FIXME: ff support!
  match.index(x)
}


## Compute the 'bamlss.frame' 'x' master object.
design.construct <- function(formula, data = NULL, knots = NULL,
  model.matrix = TRUE, smooth.construct = TRUE, binning = FALSE,
  before = TRUE, gam.side = NULL, model = NULL, drop = NULL,
  scale.x = TRUE, absorb.cons = NULL, sparse.cons = 0, specials = NULL, ...)
{
  if(!model.matrix & !smooth.construct)
    return(NULL)

  if(is.null(gam.side))
    gam.side <- if(binning) FALSE else TRUE

  if(inherits(formula, "bamlss.frame")) {
    data <- if(is.null(data)) model.frame(formula) else data
    formula <- formula(formula)
  }
  if(!inherits(formula, "bamlss.terms")) {
    if(!inherits(formula, "bamlss.formula"))
      formula <- bamlss.formula(formula, ...)
    if(inherits(formula, "bamlss.formula"))
      formula <- terms.bamlss.formula(formula, data = data, ...)
  }
  formula <- formula.bamlss.terms(formula)
  if(is.null(data))
    stop("data needs to be supplied!")

  no_ff <- !inherits(data, "ffdf")

  if(!is.character(data) & no_ff) {
    if(!inherits(data, "data.frame"))
      data <- as.data.frame(data)
  }
  if(is.character(data)) {
    ## data <- read.table.ffdf(file = data,
    ##   na.strings = "", header = TRUE, sep = ",")
    ## FIXME: ff data.frames!
    data <- read.table(file = data, header = TRUE, ...)
  }
  if(inherits(data, "ffdf")) {
    before <- TRUE
    gam.side <- FALSE
    if(is.null(binning))
      binning <- TRUE
  }
  if(!is.null(model))
    formula <- model.terms(formula, model)
  if(!binning)
    binning <- NULL

  assign.design <- function(obj, dups = NULL)
  {
    if(!is.null(dups) & no_ff) {
      if(any(dups)) {
        mi <- match.index(data[, all.vars(obj$fake.formula), drop = FALSE])
        obj[names(mi)] <- mi
        data <- subset(data, !dups)
      }
    }
    obj$binning <- binning
    if(!all(c("formula", "fake.formula") %in% names(obj)))
      return(obj)
    if(model.matrix) {
      if(!inherits(data, "ffdf")) {
        obj$model.matrix <- model.matrix(drop.terms.bamlss(obj$terms,
          sterms = FALSE, keep.response = FALSE, data = data, specials = specials), data = data)
        if(ncol(obj$model.matrix) > 0) {
          if(scale.x)
            obj$model.matrix <- scale.model.matrix(obj$model.matrix)
        } else obj$model.matrix <- NULL
      } else {
        mm_terms <- drop.terms.bamlss(obj$terms,
          sterms = FALSE, keep.response = FALSE, data = NULL, specials = specials)
        mm_intercept <- attr(mm_terms, "intercept") > 0
        mm_vars <- all.vars.formula(mm_terms)
        if(mm_intercept) {
          ## FIXME: ff support!
          ## obj$model.matrix <- ffdf("Intercept" = ff(1, length = nrow(data)))
          obj$model.matrix <- cbind("Intercept" = rep(1, length = nrow(data)))
        }
        if(!is.null(mm_vars)) {
          if(!all(mm_vars %in% colnames(data)))
            stop("variables missing in supplied data!")
          if(mm_intercept) {
            for(v in mm_vars)
              obj$model.matrix[[v]] <- data[[v]]
          } else {
            obj$model.matrix <- data[mm_vars]
          }
          if(is.numeric(binning)) {
            for(v in mm_vars) {
              obj$model.matrix[[v]] <- round(obj$model.matrix[[v]], binning)
            }
          }
        }
        if(!is.null(obj$model.matrix)) {
          bind <- match.index.ff(obj$model.matrix)
          obj$model.matrix <- as.matrix(as.data.frame(obj$model.matrix[bind$nodups, , drop = FALSE]))
          if(is.null(colnames(obj$model.matrix)) & ncol(obj$model.matrix) < 2) {
            if(mm_intercept)
              colnames(obj$model.matrix) <- "(Intercept)"
            if(!is.null(mm_vars))
              colnames(obj$model.matrix) <- mm_vars
          }
          attr(obj$model.matrix, "binning") <- bind
        }
      }
    }
    if(smooth.construct) {
      tx <- drop.terms.bamlss(obj$terms,
        pterms = FALSE, keep.response = FALSE, data = data, specials = specials)
      sid <- unlist(attr(tx, "specials"))
      if(!length(sid))
        sid <- NULL
      smt <- NULL
      if(!is.null(sid)) {
        sterms <- sterm_labels <- attr(tx, "term.labels")[sid]
        sterms <- lapply(sterms, function(x) { eval(parse(text = x)) })
        nst <- NULL
        for(j in seq_along(sterms)) {
          sl <- sterms[[j]]$label
          if(is.null(sl))
            sl <- sterm_labels[j]
          nst <- c(nst, sl)
        }
        names(sterms) <- nst
        for(tsm in sterms) {
          if(is.null(tsm$xt))
            tsm$xt <- list()
          if(is.null(tsm$xt$binning))
            tsm$xt$binning <- binning
          if(!is.null(tsm$xt$binning)) {
            if(!is.logical(tsm$xt$binning)) {
              for(tsmt in tsm$term) {
                if(!inherits(data, "ffdf")) {
                  if(!is.factor(data[[tsmt]]))
                    data[[tsmt]] <- round(data[[tsmt]], digits = tsm$xt$binning)
                } else {
                  if(is.numeric(binning))
                    data[[tsmt]] <- round(data[[tsmt]], digits = tsm$xt$binning)
                }
              }
            }
          }
        }
        no.mgcv <- NULL
        smooth <- list()
        for(tsm in sterms) {
          special <- FALSE
          if(!is.null(tsm$special))
            special <- tsm$special
          if(!special) {
            if(is.null(tsm$xt))
              tsm$xt <- list()
            if(is.null(tsm$xt$binning))
              tsm$xt$binning <- binning
            acons <- TRUE
            if(inherits(tsm, "tensor.smooth.spec")) {
              if(!is.null(tsm$margin[[1]]$xt$center))
                acons <- tsm$margin[[1]]$xt$center
            } else {
              if(!is.null(tsm$xt$center))
                acons <- tsm$xt$center
            }
            tsm$xt$center <- acons
            tsm$xt$before <- before
            if(!is.null(tsm$xt$binning)) {
              term.names <- c(tsm$term, if(tsm$by != "NA") tsm$by else NULL)
              if(!inherits(data, "ffdf")) {
                tsm$binning <- match.index(data[, term.names, drop = FALSE])
                tsm$binning$order <- order(tsm$binning$match.index)
                tsm$binning$sorted.index <- tsm$binning$match.index[tsm$binning$order]
              } else {
                tsm$binning <- match.index.ff(data[term.names])
                attr(tsm, "ff") <- TRUE
              }
              if(!inherits(data, "ffdf")) {
                smt <- smoothCon(tsm, if(before) data[tsm$binning$nodups, term.names, drop = FALSE] else data,
                  knots, absorb.cons = if(is.null(absorb.cons)) acons else absorb.cons, sparse.cons = sparse.cons)
                smooth <- c(smooth, smt)
              } else {
                xdata <- as.data.frame(data[tsm$binning$nodups, term.names, drop = FALSE])
                if(!inherits(xdata, "data.frame")) {
                  xdata <- data.frame(xdata)
                  names(xdata) <- term.names
                }
                smt <- smoothCon(tsm, xdata,
                  knots, absorb.cons = if(is.null(absorb.cons)) acons else absorb.cons,
                  sparse.cons = sparse.cons)
                smooth <- c(smooth, smt)
              }
            } else {
              smt <- smoothCon(tsm, data, knots,
                absorb.cons = if(is.null(absorb.cons)) acons else absorb.cons,
                sparse.cons = sparse.cons)
              smooth <- c(smooth, smt)
            }
          } else {
            if(is.null(tsm$by))
              tsm$by <- "NA"
            if(inherits(tsm, "mrf.smooth.spec")) {
              if(!is.null(tsm$xt$map)) {
                vl <- levels(data[[tsm$term]])
                mapn <- names(tsm$xt$map)
                if(!all(mapn %in% vl))
                  levels(data[[tsm$term]]) <- c(vl, mapn[!(mapn %in% vl)])
                tsm$xt$polys <- as.list(tsm$xt$map)
              }
            }
            if((tsm$by != "NA") & is.factor(data[[tsm$by]])) {
              fm <- model.matrix(as.formula(paste("~ -1 +", tsm$by)), data = data)
              tlab <- tsm$label
              byvar <- tsm$by
              for(jj in 1:ncol(fm)) {
                tsm$by <- colnames(fm)[jj]
                tsm$label <- gsub(byvar, colnames(fm)[jj], tlab, fixed = TRUE)
                data[[colnames(fm)[jj]]] <- fm[, jj]
                smt2 <- smooth.construct(tsm, data, knots)
                if(inherits(tsm, "no.mgcv") | inherits(smt2, "no.mgcv")) {
                  no.mgcv <- c(no.mgcv, list(smt2))
                } else {
                  class(smt2) <- c(class(smt2), "mgcv.smooth")
                  smt <- list(smt2)
                  smooth <- c(smooth, smt)
                }
              }
            } else {
              smt2 <- smooth.construct(tsm, data, knots)
              if(inherits(tsm, "no.mgcv") | inherits(smt2, "no.mgcv")) {
                no.mgcv <- c(no.mgcv, if(!inherits(smt2, "smooth.list")) list(smt2) else smt2)
              } else {
                class(smt2) <- c(class(smt2), "mgcv.smooth")
                smt <- if(!inherits(smt2, "smooth.list")) list(smt2) else smt2
                smooth <- c(smooth, smt)
              }
            }
          }
        }
        if(length(smooth) > 0) {
          if(gam.side) {
            if(is.null(obj$model.matrix)) {
              Xp <- model.matrix(drop.terms.bamlss(obj$terms,
                sterms = FALSE, keep.response = FALSE, data = data, specials = specials), data = data)
              smooth <- try(gam.side(smooth, Xp, tol = .Machine$double.eps^.5), silent = TRUE)
            } else {
              smooth <- try(gam.side(smooth, obj$model.matrix, tol = .Machine$double.eps^.5), silent = TRUE)
            }
            if(inherits(smooth, "try-error")) {
              cat("---\n", smooth, "---\n")
              if(binning)
                stop("gam.side() produces an error when binning, try to set before = FALSE or set gam.side = FALSE!")
              else
                stop("gam.side() produces an error, try to set gam.side = FALSE!")
            }
          }
          sme <- NULL
          if(smooth.construct)
            sme <- expand.t2.smooths(smooth)
          if(is.null(sme)) {
            original.smooth <- NULL
          } else {
            original.smooth <- smooth
            smooth <- sme
            rm(sme)
          }
        }
        for(j in seq_along(smooth))
          smooth[[j]][["X.dim"]] <- ncol(smooth[[j]]$X)
        if(!is.null(no.mgcv))
          smooth <- c(smooth, no.mgcv)
        if(length(smooth))
          obj$smooth.construct <- smooth
      }
    }
    if(!is.null(obj$smooth.construct)) {
      sl <- NULL
      for(j in seq_along(obj$smooth.construct)) {
        slj <- obj$smooth.construct[[j]]$label
        if(!is.null(obj$smooth.construct[[j]]$by)) {
          if(obj$smooth.construct[[j]]$by != "NA") {
            if(grepl(pat <- paste("):", obj$smooth.construct[[j]]$by, sep = ""), slj, fixed = TRUE)) {
              slj <- gsub(pat, paste(",by=", obj$smooth.construct[[j]]$by, "):", sep = ""), slj, fixed = TRUE)
              slj <- strsplit(slj, "", fixed = TRUE)[[1]]
              if(slj[length(slj)] == ":")
                slj <- slj[-length(slj)]
              slj <- paste(slj, collapse = "")
              obj$smooth.construct[[j]]$label <- slj
            }
          }
        } else obj$smooth.construct[[j]]$by <- "NA"
        sl <- c(sl, slj)
      }
      if(length(unique(sl)) < length(sl)) {
        sld <- sl[duplicated(sl)]
        for(j in seq_along(sld)) {
          for(jj in which(sl == sld[j])) {
            clj <- class(obj$smooth.construct[[jj]])
            clj <- strsplit(clj, ".", fixed = TRUE)[[1]][1]
            if(clj == "random")
              clj <- "re"
            sl[jj] <- paste(sl[jj], clj, sep = ":")
          }
        }
      }
      names(obj$smooth.construct) <- sl
    }
    if(!is.null(drop)) {
      take <- c("model.matrix", "smooth.construct")[c(model.matrix, smooth.construct)]
      obj[!(names(obj) %in% take)] <- NULL
    }

    obj
  }

  if(!all(c("formula", "fake.formula") %in% names(formula))) {
    for(j in seq_along(formula)) {
      if(!all(c("formula", "fake.formula") %in% names(formula[[j]]))) {
        for(i in seq_along(formula[[j]])) {
          formula[[j]][[i]] <- assign.design(formula[[j]][[i]],
            if(i > 1) duplicated(data[, all.vars(formula[[j]][[i]]$fake.formula), drop = FALSE]) else NULL)
        }
      } else formula[[j]] <- assign.design(formula[[j]])
    }
  } else formula <- assign.design(formula)

  if((!all(c("formula", "fake.formula") %in% names(formula))) & smooth.construct) {
    for(i in seq_along(formula)) {
      if(!all(c("formula", "fake.formula") %in% names(formula[[i]]))) {
        for(j in seq_along(formula[[i]])) {
          if(!is.null(formula[[i]][[j]]$smooth.construct)) {
            for(k in seq_along(formula[[i]][[j]]$smooth.construct)) {
              if(is.null(formula[[i]][[j]]$smooth.construct[[k]]$fit.fun))
                formula[[i]][[j]]$smooth.construct[[k]]$fit.fun <- make.fit.fun(formula[[i]][[j]]$smooth.construct[[k]])
              if(is.null(formula[[i]][[j]]$smooth.construct[[k]]$prior)) {
                priors <- make.prior(formula[[i]][[j]]$smooth.construct[[k]])
                formula[[i]][[j]]$smooth.construct[[k]]$prior <- priors$prior
                formula[[i]][[j]]$smooth.construct[[k]]$grad <- priors$grad
                formula[[i]][[j]]$smooth.construct[[k]]$hesss <- priors$hess
              }
            }
          }
        }
      } else {
        if(!is.null(formula[[i]]$smooth.construct)) {
          for(j in seq_along(formula[[i]]$smooth.construct)) {
            if(is.null(formula[[i]]$smooth.construct[[j]]$fixed))
              formula[[i]]$smooth.construct[[j]]$fixed <- FALSE
            if(length(formula[[i]]$smooth.construct[[j]]$S)) {
              for(sj in seq_along(formula[[i]]$smooth.construct[[j]]$S)) {
                if(!is.list(formula[[i]]$smooth.construct[[j]]$S[[sj]]) & !is.function(formula[[i]]$smooth.construct[[j]]$S[[sj]])) {
                  nc <- ncol(formula[[i]]$smooth.construct[[j]]$S[[sj]])
                  formula[[i]]$smooth.construct[[j]]$S[[sj]] <- formula[[i]]$smooth.construct[[j]]$S[[sj]] + diag(1e-05, nc, nc)
                }
              }
            }
            if(is.null(formula[[i]]$smooth.construct[[j]]$fit.fun))
              formula[[i]]$smooth.construct[[j]]$fit.fun <- make.fit.fun(formula[[i]]$smooth.construct[[j]])
            if(is.null(formula[[i]]$smooth.construct[[j]]$prior)) {
              priors <- make.prior(formula[[i]]$smooth.construct[[j]])
              formula[[i]]$smooth.construct[[j]]$prior <- priors$prior
              formula[[i]]$smooth.construct[[j]]$grad <- priors$grad
              formula[[i]]$smooth.construct[[j]]$hess <- priors$hess
            }
          }
        }
      }
    }
  } else {
    if(!is.null(formula$smooth.construct)) {
      for(j in seq_along(formula$smooth.construct)) {
        if(is.null(formula$smooth.construct[[j]]$fixed))
          formula$smooth.construct[[j]]$fixed <- FALSE
        if(length(formula[[i]]$smooth.construct[[j]]$S)) {
          for(sj in seq_along(formula$smooth.construct[[j]]$S)) {
            if(!is.list(formula$smooth.construct[[j]]$S[[sj]]) & !is.function(formula$smooth.construct[[j]]$S[[sj]])) {
              nc <- ncol(formula$smooth.construct[[j]]$S[[sj]])
              formula$smooth.construct[[j]]$S[[sj]] <- formula$smooth.construct[[j]]$S[[sj]] + diag(1e-05, nc, nc)
            }
          }
        }
        if(is.null(formula$smooth.construct[[j]]$fit.fun))
          formula$smooth.construct[[j]]$fit.fun <- make.fit.fun(formula$smooth.construct[[j]])
        if(is.null(formula$smooth.construct[[j]]$prior)) {
          priors <- make.prior(formula$smooth.construct[[j]])
          formula$smooth.construct[[j]]$prior <- priors$prior
          formula$smooth.construct[[j]]$grad <- priors$grad
          formula$smooth.construct[[j]]$hess <- priors$hess
        }
      }
    }
  }

  attr(formula, "specials") <- NULL
  attr(formula, ".Environment") <- NULL
  class(formula) <- "list"
  if(!is.null(drop)) {
    if(drop & (length(formula) < 2))
      formula <- formula[[1]]
  }

  return(formula)
}


## Functions for sparse matrices.
sparse.matrix.index <- function(x, ...)
{
  if(is.null(dim(x)))
    return(NULL)
  index <- apply(x, 1, function(x) {
    which(x != 0)
  })
  if(length(index) < 1)
    return(NULL)
  if(is.list(index)) {
    n <- max(sapply(index, length))
    index <- lapply(index, function(x) {
      if((nx <- length(x)) < n)
        x <- c(x, rep(-1L, length = n - nx))
      x
    })
    index <- do.call("rbind", index)
  } else {
    index <- if(is.null(dim(index))) {
      matrix(index, ncol = 1)
    } else t(index)
  }
  storage.mode(index) <- "integer"
  index
}


## Bandwidth minimization permutation.
sparse.matrix.ordering <- function(x, ...)
{
  x <- as.spam(x)
  i <- try(ordering(chol.spam(x)), silent = TRUE)
  if(inherits(i, "try-error"))
    i <- 1:nrow(x)
  return(i)
}


## Setup sparse indeces for various algorithms.
sparse.setup <- function(x, S = NULL, ...)
{
  symmetric <- nrow(x) == ncol(x)
  index.matrix <- sparse.matrix.index(x, ...)
  if(!symmetric)
    x <- crossprod(x)
  if(!is.null(S)) {
    if(!is.list(S))
      S <- list(S)
    for(j in seq_along(S)) {
      x <- x + if(length(S[[j]]) < 1) 0 else { if(is.function(S[[j]])) S[[j]](c("b" = rep(0, attr(S[[j]], "npar")))) else S[[j]] }
    }
  }
  index.crossprod <- if(!symmetric) sparse.matrix.index(x, ...) else NULL
  setup <- list(
    "matrix" = index.matrix,
    "crossprod" = index.crossprod
  )
  if(!is.null(index.crossprod)) {
    idf <- as.factor(apply(setup$crossprod, 1, paste, collapse = ","))
    if((nlevels(idf) > 1) & (nlevels(idf) < nrow(setup$crossprod))) {
      setup$block.index <- split(as.integer(1:nrow(setup$crossprod)), idf)
      setup$is.diagonal <- all(sapply(setup$block.index, length) == 1)
    }
  }
  return(setup)
}


## Sparse cholesky decomposition,
## returns the lower triangle.
sparse.chol <- function(x, index, ...)
{
  if(all(dim(x) < 2))
    return(sqrt(x))

  imat <- index[["matrix"]]
  p <- index[["ordering"]]

  # imat: index matrix of a[p,p]
  # ??? check for positive definiteness?
  n <- nrow(x)
  l <- matrix(0, nrow = n, ncol = n)
  
  # First column simplified (no elements to sum up)
  l[1, 1] <- (x[p[1], p[1]])^0.5
  for(i in imat[1,][imat[1,]>1]) {
    l[i, 1] <- x[p[i], p[1]] / l[1, 1]
  }
  c <- 1
  for(j in p[2:(n-1)]) {
    c <- c + 1
    l[c, c] <- (x[j, j] - sum(l[c, 1:(c - 1)]^2))^0.5
    # use only non-zero entries in lower subdiagonal
    for(i in imat[c,][imat[c,] > c]) {   
      l[i, c] <- (x[p[i], p[c]] -
                    sum(l[i, 1:(c - 1)] * l[c, 1:(c - 1)])) / l[c, c]
    }
  }
  # last column simplified: no subdiagonal - maybe still leave in loop?
  l[n, n] <- (x[p[n], p[n]] - sum(l[n, 1:(n - 1)]^2))^0.5
  j <- c(1:n)[p]
  return(l[j,j])
}


## Sparse forward substitution.
## L %*% x = bn
## with bn = t(P) %*% b
sparse.forwardsolve <- function(l, x, index, ...)
{
  if(all(dim(l) < 2))
    return(x / l)
  imat <- index[["matrix"]]
  p <- index[["ordering"]]
  n <- ncol(l)
  Pt <- diag(n)[p,]
  xn <- Pt %*% x
  y <- matrix(rep(NA, n), ncol=1)
  y[1] <- xn[1]/l[1, 1]
  for(i in 2:n){
    y[i] <- xn[i]/l[i, i]  
    if(max(imat[i,]) > 0){
      y[i] <- y[i] - sum(l[i, imat[i,][imat[i,] > 0] ] * y[imat[i,][imat[i,] > 0]])/l[i, i]
    }
  }
  return(y)
}


## Sparse backward substitution.
## t(L) %*% xn = x,
## with x = P %*% xn.
sparse.backsolve <- function(r, x, index = NULL, ...)
{
  if(all(dim(x) < 2))
    return(r / x)
  imat <- index[["matrix"]]
  p <- index[["ordering"]]
  n <- ncol(r)
  P <- diag(n)[,p]
  xn <- rep(NA, n)
  xn[n] <- x[n]/r[n,n]
  for(i in (n-1):1){
    xn[i] <- x[i]/r[i,i]
    if(max(imat[i,]) > 0){
      xn[i] <- xn[i] - sum(r[imat[i,][imat[i,] > 0],i ] * xn[imat[i,][imat[i,] > 0]])/r[i, i]
    }
  }
  x <- P %*% xn
  return(x)
}


## Sparse matrix solve.
sparse.solve <- function(a, b, index, ...)
{
  if(all(dim(a) < 2))
    return(b / a)
  id <- if(!("crossprod" %in% names(index))) "matrix" else "crossprod"
  L <- sparse.chol(a, index = list("matrix" = index[[id]], "ordering" = index[["ordering"]]), ...)
  y <- sparse.forwardsolve(L, b, index = list("matrix" = index$forward, "ordering" = index[["ordering"]]), ...)
  z <- sparse.backsolve(L, y, list("matrix" = index$backward, "ordering" = index[["ordering"]]), ...)
  return(z)
}


## Computation of fitted values with index matrices.
sparse.matrix.fit.fun <- function(X, b, index = NULL)
{
  if(!is.null(index)) {
    if(nrow(index) != nrow(X))
      return(drop(X %*% b))
  }
  fit <- if(inherits(X, "dgCMatrix") | is.null(index) | inherits(X, "Matrix")) {
    drop(X %*% b)
  } else .Call("sparse_matrix_fit_fun", X, b, index, PACKAGE = "bamlss")
  return(fit)
}


## The model term fitting function.
make.fit.fun <- function(x, type = 1)
{
  ff <- function(X, b, expand = TRUE, no.sparse.setup = FALSE) {
    if(!is.null(names(b))) {
      b <- if(!is.null(x$pid)) b[x$pid$b] else get.par(b, "b")
    }
    if(inherits(X, "spam") | inherits(X, "Matrix")) {
      f <- as.matrix(X %*% b)
    } else {
      what <- if(type < 2) "matrix" else "grid.matrix"
      f <- if(is.null(x$sparse.setup[[what]]) | no.sparse.setup) {
        drop(X %*% b)
      } else sparse.matrix.fit.fun(X, b, x$sparse.setup[[what]])
    }
    if(!is.null(x$binning$match.index) & expand) {
      if(inherits(x$binning$match.index, "ff")) {
        ## f <- as.ff(f)
        ## FIXME: ff support!
        return(f[x$binning$match.index])
      }
      f <- f[x$binning$match.index]
    }
    if(!is.null(x$xt$force.center))
      f <- f - mean(f, na.rm = TRUE)
    return(as.numeric(f))
  }
  return(ff)
}


## The prior function.
make.prior <- function(x, sigma = 0.1)
{
  prior <- NULL
  if(!is.null(x$xt$prior)) {
    prior <- x$xt$prior
    if(is.character(x$xt$prior)) {
      prior <- tolower(prior)
      if(!(prior %in% c("ig", "hc", "sd", "hn", "hn.lasso")))
        stop(paste('smoothing variance prior "', prior, '" not supported!', sep = ''))
    }
  } else {
    prior <- "ig"
  }

  if(!is.null(x$margin)) {
    if(!is.null(x$margin[[1]]$xt)) {
      xt <- x$margin[[1]]$xt
      if(is.null(names(xt)) & (length(xt) == 1)) {
        if(length(xt[[1]]) > 0)
          xt <- xt[[1]]
      }
      x$xt <- c(x$xt, xt)
    }
  }

  if(!is.function(prior)) {
    rval <- list()

    a <- if(is.null(x$xt[["a"]])) {
      if(is.null(x[["a"]])) 1e-04 else x[["a"]]
    } else x$xt[["a"]]
    b <- if(is.null(x$xt[["b"]])) {
      if(is.null(x[["b"]])) 1e-04 else x[["b"]]
    } else x$xt[["b"]]
    theta <- if(is.null(x$xt[["theta"]])) {
      x[["theta"]]
    } else x$xt[["theta"]]

    if(is.null(theta)) {
      theta <- switch(prior,
        "sd" = 0.00877812,
        "hc" = 0.01034553,
        "hn" = 0.1457644
      )
    }

    fixed <- if(is.null(x$fixed)) FALSE else x$fixed

    igs <- log((b^a)) - log(gamma(a))
    var_prior_fun <- switch(prior,
      "ig" = function(tau2) { igs + (-a - 1) * log(tau2) - b / tau2 },
      "hc" = function(tau2) { -log(1 + tau2 / (theta^2)) - 0.5 * log(tau2) - log(theta^2) },
      "sd" = function(tau2) { -0.5 * log(tau2) + 0.5 * log(theta) - (tau2 / theta)^(0.5) },
      "hn" = function(tau2) { -0.5 * log(tau2) - tau2 / (2 * theta^2) },
      "hn.lasso0" = function(tau2) { -0.2257913 - log(sigma) - tau2^2/(2 * sigma^2) },
      "hn.lasso" = function(tau2) {
        theta <- sqrt(pi) / (sigma * sqrt(2))
        log(2 * theta / pi) - (tau2^2 * theta^2) / pi
      }
    )
    
    rval$prior <- function(parameters) {
      if(is.null(x$pid)) {
        if(!is.null(names(parameters))) {
          gamma <- get.par(parameters, "b")
          tau2 <- get.par(parameters, "tau2")
        } else {
          gamma <- parameters
          tau2 <- numeric(0)
        }
      } else {
        gamma <- parameters[x$pid$b]
        tau2 <-  parameters[x$pid$tau2]
      }
      if(fixed | !length(tau2)) {
        lp <- sum(dnorm(gamma, sd = 1000, log = TRUE))
      } else {
        if(length(tau2) < 2) {
          K <- if(is.function(x$S[[1]])) x$S[[1]](c(parameters, x$fixed.hyper)) else x$S[[1]]
          if(is.null(x$rank))
            x$rank <- qr(K)$rank
          lp <- -log(tau2) * x$rank / 2 + drop(-0.5 / tau2 * t(gamma) %*% K %*% gamma) + var_prior_fun(tau2)
        } else {
          ld <- 0
          P <- if(inherits(x$X, "Matrix")) Matrix(0, ncol(x$X), ncol(x$X)) else 0
          for(j in seq_along(tau2)) {
            P <- P + 1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](c(parameters, x$fixed.hyper)) else x$S[[j]]
            ld <- ld + var_prior_fun(tau2[j])
          }
          ##lp <- dmvnorm(gamma, sigma = matrix_inv(P), log = TRUE) + ld
          dP <- determinant(P, logarithm = TRUE)
          dP <- as.numeric(dP$modulus) * as.numeric(dP$sign)
          lp <- 0.5 * dP - 0.5 * (t(gamma) %*% P %*% gamma) + ld
        }
      }
      return(as.numeric(lp))
    }

    attr(rval$prior, "var_prior") <- prior

    rval$grad <- function(score = NULL, parameters, full = TRUE) {
      gamma <- get.par(parameters, "b")
      tau2 <-  get.par(parameters, "tau2")
      grad2 <- NULL
      if(x$fixed | !length(tau2)) {
        grad <- rep(0, length(gamma))
      } else {
        grad <- 0; grad2 <- NULL
        for(j in seq_along(tau2)) {
          tauS <- -1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](c(parameters, x$fixed.hyper)) else x$S[[j]]
          grad <- grad + tauS %*% gamma
          if(full & !is.null(tau2[j])) {
            grad2 <- c(grad2, drop(-x$rank[j] / (2 * tau2[j]) - 1 / (2 * tau2[j]^2) * (if(is.function(x$S[[j]])) x$S[[j]](c(parameters, x$fixed.hyper)) else x$S[[j]]) %*% gamma + (-x$a - 1) / tau2[j] + x$b / (tau2[j]^2)))
            x$X <- cbind(x$X, 0)
          }
          grad <- drop(grad)
        }
      }
      if(!is.null(score)) {
        grad <- if(!is.null(x$xbin.ind)) {
          drop(crossprod(x$X[x$xbin.ind, , drop = FALSE], score)) + c(grad, grad2)
        } else drop(crossprod(x$X, score)) + c(grad, grad2)
      } else grad <- c(grad, grad2)
      return(grad)
    }

    rval$hess <- function(score = NULL, parameters, full = FALSE) {
      tau2 <- get.par(parameters, "tau2")
      if(x$fixed | !length(tau2)) {
        k <- length(get.par(parameters, "b"))
        hx <- matrix(0, k, k)
      } else {
        hx <- 0
        for(j in seq_along(tau2)) {
          hx <- hx + (1 / tau2[j]) * if(is.function(x$S[[j]])) x$S[[j]](c(parameters, x$fixed.hyper)) else x$S[[j]]
        }
      }
      return(hx)
    }

    return(rval)
  } else {
    return(prior(x))
  }
}


## Fast block diagonal crossproduct with weights.
do.XWX <- function(x, w, index = NULL)
{
  if(is.null(index) | inherits(x, "dgCMatrix")) {
    rval <- crossprod(x / w, x)
  } else {
    if(is.null(dim(index)))
      index <- matrix(index, ncol = 1)
    rval <- .Call("do_XWX", x, w, index, PACKAGE = "bamlss")
  }
  rval
}


## Get the model.frame.
model.frame.bamlss <- model.frame.bamlss.frame <- function(formula, ...) 
{
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
  mf <- if(length(nargs) || is.null(formula$model.frame)) {
    fcall <- formula$call
    fcall[[1L]] <- quote(bamlss.model.frame)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$formula)
    if(is.null(env))
      env <- parent.frame()
    ft <- eval(fcall[["formula"]], env)
    if(!is.null(attr(ft, "orig.formula"))) {
      fcall["formula"] <- parse(text = paste("attr(", fcall["formula"], ", 'orig.formula')", sep = ""))
    }
    nf <- names(fcall)
    nf <- nf[!(nf %in% names(formals(bamlss.model.frame)))]
    nf <- nf[nchar(nf) > 0]
    fcall[nf] <- NULL
    fcall["drop.unused.levels"] <- FALSE
    if(is.null(fcall["family"]))
      fcall["family"] <- parse(text = "gaussian_bamlss()")
    eval(fcall, env)
  } else formula$model.frame
  mf
}


## Search for parts in models, optionally extract.
model.search <- function(x, what, model = NULL, part = c("x", "formula", "terms"),
  extract = FALSE, drop = FALSE)
{
  if(!inherits(x, "bamlss.formula") & !inherits(x, "bamlss.frame"))
    stop("x must be a 'bamlss.formula' or 'bamlss.frame' object!")
  part <- match.arg(part)
  if(is.null(x[[part]]))
    return(FALSE)
  x <- model.terms(x, model = model, part = part)
  elmts <- c("formula", "fake.formula")
  nx <- names(x)
  rval <- list()
  for(i in nx) {
    if(!all(elmts %in% names(x[[i]]))) {
      rval[[i]] <- list()
      for(j in names(x[[i]])) {
        rval[[i]][[j]] <- if(is.null(x[[i]][[j]][[what]])) FALSE else TRUE
        if(extract & rval[[i]][[j]])
          rval[[i]][[j]] <- x[[i]][[j]][[what]]
      }
    } else {
      rval[[i]] <- if(is.null(x[[i]][[what]])) FALSE else TRUE
      if(extract & rval[[i]])
        rval[[i]] <- x[[i]][[what]]
    }
  }
  if(!extract) {
    rval <- unlist(rval)
  } else {
    if(drop & (length(rval) < 2))
      rval <- rval[[1]]
  }
  rval
}


## Wrapper for design construct extraction.
extract.design.construct <- function(object, data = NULL,
  knots = NULL, model = NULL, drop = TRUE, what = c("model.matrix", "smooth.construct"),
  specials = NULL, ...)
{
  if(!inherits(object, "bamlss.frame") & !inherits(object, "bamlss.formula") & !inherits(object, "bamlss.terms"))
    stop("object must be a 'bamlss.frame', 'bamlss.formula' or 'bamlss.terms' object!")
  what <- match.arg(what)
  model.matrix <- what == "model.matrix"
  smooth.construct <- what == "smooth.construct"
  if(inherits(object, "bamlss.frame")) {
    if(!is.null(data)) {
      object$model.frame <- NULL
      object <- design.construct(object, data = data, knots = knots,
        model.matrix = model.matrix, smooth.construct = smooth.construct,
        model = model, drop = drop, specials = specials, ...)
    } else {
      if(!all(model.search(object, what, model, part = "x"))) {
        object <- design.construct(object, model.matrix = model.matrix,
          smooth.construct = smooth.construct, model = model, drop = TRUE,
          specials = specials, ...)
      } else {
        object <- model.search(object, what, model, extract = TRUE, drop = drop, part = "x")
      }
    }
  } else {
    if(is.null(data))
      stop("argument data is missing!")
    object <- design.construct(object, data = data, knots = knots,
      model.matrix = model.matrix, smooth.construct = smooth.construct, model = model,
      drop = drop, specials = specials, ...)
  }
  if(!is.null(drop)) {
    if(length(object) & drop & (length(object) < 2))
      object <- object[[1]]
  }
  if(!length(object))
    return(NULL)
  mostattributes(object) <- NULL
  attr(object, "orig.formula") <- NULL
  if(what == "model.matrix") {
    if(is.list(object)) {
      for(j in seq_along(object)) {
        if(is.list(object[[j]])) {
          if((length(object[[j]]) < 2) & (names(object[[j]]) == "model.matrix")) {
            object[[j]] <- object[[j]]$model.matrix
          }
        }
      }
    }
  }
  return(object)
}


## Model matrix extractor.
model.matrix.bamlss.frame <- model.matrix.bamlss.formula <- model.matrix.bamlss.terms <- function(object, data = NULL, model = NULL, drop = TRUE, scale.x = FALSE, ...)
{
  extract.design.construct(object, data = data,
    knots = NULL, model = model, drop = drop, what = "model.matrix",
    scale.x = scale.x)
}


## Extract smooth constructs.
smooth.construct <- function(object, data, knots, ...)
{
  UseMethod("smooth.construct")
}

smooth.construct.bamlss.frame <- smooth.construct.bamlss.formula <- smooth.construct.bamlss.terms <- function(object, data = NULL, knots = NULL, model = NULL, drop = TRUE, ...)
{
  extract.design.construct(object, data = data,
    knots = knots, model = model, drop = drop, what = "smooth.construct")
}


## Extract/initialize parameters.
parameters <- function(x, model = NULL, start = NULL, fill = c(0, 0.0001),
  list = FALSE, simple.list = FALSE, extract = FALSE, ...)
{
  if(inherits(x, "bamlss") | extract) {
    if(!is.null(x$parameters)) {
      if(is.null(model)) {
        if(list) {
          return(x$parameters)
        } else {
          if(inherits(x$parameters, "data.frame") | inherits(x$parameters, "matrix")) {
            args <- list(...)
            mstop <- if(is.null(args$mstop)) nrow(x$parameters) else args$mstop
            return(x$parameters[mstop, ])
          } else return(unlist(x$parameters))
        }
      } else {
        if(is.list(x$parameters)) {
          if(list) return(x$parameters[model]) else return(unlist(x$parameters[model]))
        } else {
          if(!is.character(model))
            model <- names(x$terms)[model]
          rp <- grep(paste(model, ".", sep = ""), names(x$parameters), fixed = TRUE, value = TRUE)
          if(inherits(x$parameters, "data.frame") | inherits(x$parameters, "matrix")) {
            rp <- grep(paste(model, ".", sep = ""), colnames(x$parameters), fixed = TRUE, value = TRUE)
            args <- list(...)
            mstop <- if(is.null(args$mstop)) nrow(x$parameters) else args$mstop
            return(x$parameters[mstop, rp])
          } else return(x$parameters[rp])
        }
      }
    }
  }
  if(inherits(x, "bamlss.frame")) {
    if(is.null(x$x)) {
      x <- design.construct(x, data = x$model.frame,
        knots = x$knots, model.matrix = TRUE, smooth.construct = TRUE, model = NULL)
    } else x <- x$x
  }
  fill <- rep(fill, length.out = 2)
  if(!is.null(start)) {
    if(is.list(start))
      start <- unlist(start)
  }
  par <- list()
  for(i in names(x)) {
    par[[i]] <- list()
    if(!all(c("formula", "fake.formula") %in% names(x[[i]]))) {
      for(j in names(x[[i]])) {
        par[[i]][[j]] <- list()
        if(!is.null(x[[i]][[j]]$model.matrix)) {
          nc <- ncol(x[[i]][[j]]$model.matrix)
          if(simple.list) {
            par[[i]][[j]]$p <- fill[1]
          } else {
            par[[i]][[j]]$p <- rep(fill[1], length = nc)
            if(is.null(cn <- colnames(x[[i]][[j]]$model.matrix)))
              cn <- paste("b", 1:nc, sep = "")
            names(par[[i]][[j]]$p) <- cn
            if(!is.null(start)) {
              if(length(ii <- grep(paste(i, j, "p", sep = "."), names(start), fixed = TRUE))) {
                spar <- start[ii]
                spn <- names(spar)
                cn2 <- paste(i, j, "p", cn, sep = ".")
                take <- which(spn %in% cn2)
                if(length(take)) {
                  par[[i]][[j]]$p[which(cn2 %in% spn)] <- spar[take]
                }
              }
            }
          }
        }
        if(!is.null(x[[i]][[j]]$smooth.construct)) {
          par[[i]][[j]]$s <- list()
          for(k in names(x[[i]][[j]]$smooth.construct)) {
            if(simple.list) {
              par[[i]][[j]]$s[[k]] <- fill[1]
            } else {
              if(!is.null(x[[i]][[j]]$smooth.construct[[k]]$rand)) {
                tpar1 <- rep(fill[1], ncol(x[[i]][[j]]$smooth.construct[[k]]$rand$Xr))
                tpar2 <- rep(fill[1], ncol(x[[i]][[j]]$smooth.construct[[k]]$Xf))
                cn1 <- colnames(x[[i]][[j]]$smooth.construct[[k]]$rand$Xr)
                cn2 <- colnames(x[[i]][[j]]$smooth.construct[[k]]$Xf)
                if(is.null(cn1))
                  cn1 <- paste("b", 1:length(tpar1), ".re", sep = "")
                if(is.null(cn2))
                  cn2 <- paste("b", 1:length(tpar2), ".fx", sep = "")
                names(tpar1) <- cn1
                names(tpar2) <- cn2
                tpar <- c(tpar1, tpar2)
              } else {
                nfill <- if(is.null(x[[i]][[j]]$smooth.construct[[k]]$special.npar)) {
                  ncol(x[[i]][[j]]$smooth.construct[[k]]$X)
                } else x[[i]][[j]]$smooth.construct[[k]]$special.npar
                tpar <- rep(fill[1], nfill)
                cn <- colnames(x[[i]][[j]]$smooth.construct[[k]]$X)
                if(is.null(cn))
                  cn <- paste("b", 1:length(tpar), sep = "")
                names(tpar) <- cn
              }
              if(length(x[[i]][[j]]$smooth.construct[[k]]$S)) {
                tpar3 <- NULL
                for(kk in seq_along(x[[i]][[j]]$smooth.construct[[k]]$S)) {
                  tpar3 <- c(tpar3, fill[2])
                }
                names(tpar3) <- paste("tau2", 1:length(tpar3), sep = "")
                tpar <- c(tpar, tpar3)
              }
              par[[i]][[j]]$s[[k]] <- tpar
              if(!is.null(start)) {
                if(length(ii <- grep(paste(i, j, "s", k, sep = "."), names(start), fixed = TRUE))) {
                  spar <- start[ii]
                  cn <- names(par[[i]][[j]]$s[[k]])
                  if(length(tau2 <- grep("tau2", names(spar)))) {
                    tau2 <- spar[tau2]
                    if(length(jj <- grep("tau2", cn, fixed = TRUE))) {
                      tau2 <- rep(tau2, length.out = length(jj))
                      par[[i]][[j]]$s[[k]][jj] <- tau2
                    }
                  }
                  if(any(b <- !grepl("tau2", names(spar)))) {
                    b <- spar[b]
                    if(any(jj <- !grepl("tau2", cn, fixed = TRUE))) {
                      b <- rep(b, length.out = sum(jj))
                      par[[i]][[j]]$s[[k]][jj] <- b
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else {
      if(!is.null(x[[i]]$model.matrix)) {
        if(ncol(x[[i]]$model.matrix) > 0) {
          if(simple.list) {
            par[[i]]$p <- fill[1]
          } else {
            nc <- ncol(x[[i]]$model.matrix)
            par[[i]]$p <- rep(fill[1], length = nc)
            if(is.null(cn <- colnames(x[[i]]$model.matrix)))
              cn <- paste("b", 1:nc, sep = "")
            names(par[[i]]$p) <- cn
            if(!is.null(start)) {
              if(length(ii <- grep(paste(i, "p", sep = "."), names(start), fixed = TRUE))) {
                spar <- start[ii]
                spn <- names(spar)
                cn2 <- paste(i, "p", cn, sep = ".")
                take <- which(spn %in% cn2)
                if(length(take)) {
                  par[[i]]$p[which(cn2 %in% spn)] <- spar[take]
                }
              }
            }
          }
        }
      }
      if(!is.null(x[[i]]$smooth.construct)) {
        par[[i]]$s <- list()
        for(k in names(x[[i]]$smooth.construct)) {
          re.effect <- !is.null(x[[i]]$smooth.construct[[k]]$rand)
          if(re.effect) {
            if(is.logical(x[[i]]$smooth.construct[[k]]$rand))
              re.effect <- FALSE
          }
          if(re.effect) {
            tpar1 <- rep(fill[1], ncol(x[[i]]$smooth.construct[[k]]$rand$Xr))
            tpar2 <- rep(fill[1], ncol(x[[i]]$smooth.construct[[k]]$Xf))
            cn1 <- colnames(x[[i]]$smooth.construct[[k]]$rand$Xr)
            cn2 <- colnames(x[[i]]$smooth.construct[[k]]$Xf)
            if(is.null(cn1))
              cn1 <- paste("b", 1:length(tpar1), ".re", sep = "")
            if(is.null(cn2))
              cn2 <- paste("b", 1:length(tpar2), ".fx", sep = "")
            names(tpar1) <- cn1
            names(tpar2) <- cn2
            tpar <- c(tpar1, tpar2)
          } else {
            nfill <- if(is.null(x[[i]]$smooth.construct[[k]]$special.npar)) {
              ncol(x[[i]]$smooth.construct[[k]]$X)
            } else x[[i]]$smooth.construct[[k]]$special.npar
            tpar <- rep(fill[1], nfill)
            cn <- colnames(x[[i]]$smooth.construct[[k]]$X)
            if(is.null(cn))
              cn <- paste("b", 1:length(tpar), sep = "")
            if(!simple.list)
              names(tpar) <- cn
          }
          if(length(x[[i]]$smooth.construct[[k]]$S)) {
            tpar3 <- NULL
            for(kk in seq_along(x[[i]]$smooth.construct[[k]]$S)) {
              tpar3 <- c(tpar3, fill[2])
            }
            if(!simple.list)
              names(tpar3) <- paste("tau2", 1:length(tpar3), sep = "")
            tpar <- c(tpar, tpar3)
          }
          if(!is.null(x[[i]]$smooth.construct[[k]]$special.mpar)) {
            tpar <- x[[i]]$smooth.construct[[k]]$special.mpar()
          }
          par[[i]]$s[[k]] <- tpar
          if(!is.null(start)) {
            if(length(ii <- grep(paste(i, "s", k, sep = "."), names(start), fixed = TRUE))) {
              spar <- start[ii]
              cn <- names(par[[i]]$s[[k]])
              if(length(tau2 <- grep("tau2", names(spar)))) {
                tau2 <- spar[tau2]
                if(length(jj <- grep("tau2", cn, fixed = TRUE))) {
                  tau2 <- rep(tau2, length.out = length(jj))
                  par[[i]]$s[[k]][jj] <- tau2
                }
              }
              if(any(b <- !grepl("tau2", names(spar)))) {
                b <- spar[b]
                if(any(jj <- !grepl("tau2", cn, fixed = TRUE))) {
                  b <- rep(b, length.out = sum(jj))
                  par[[i]]$s[[k]][jj] <- b
                }
              }
            }
          }
        }
      }
    }
  }
  if(!is.null(model))
    par <- par[model]
  if(!list)
    par <- unlist(par)
  return(par)
}


par2list <- function(x)
{
  xl <- list()
  nx <- names(x)
  npl <- strsplit(nx, ".", fixed = TRUE)
  np <- unique(sapply(npl, function(x) { x[1] }))
  for(j in np) {
    xl[[j]] <- list()
    tmp <- grep(paste(j, ".", sep = ""), nx, value = TRUE)
    tmp2 <- unique(sapply(strsplit(tmp, ".", fixed = TRUE), function(x) { x[2] }))
    for(jj in tmp2) {
      tmp3 <- grep(paste(j, jj, sep = "."), nx, fixed = TRUE, value = TRUE)
      if(jj == "p") {
        xl[[j]][[jj]] <- x[tmp3]
        names(xl[[j]][[jj]]) <- gsub(paste(j, ".", jj, ".", sep = ""), "", names(xl[[j]][[jj]]), fixed = TRUE)
      } else {
        tmp4 <- grep(paste(j, jj, sep = "."), nx, fixed = TRUE, value = TRUE)
        tmp4 <- unique(sapply(strsplit(tmp4, ".", fixed = TRUE), function(x) { x[3] }))
        xl[[j]][[jj]] <- list()
        for(jjj in tmp4) {
          xl[[j]][[jj]][[jjj]] <- x[grep(paste(j, ".", jj, ".", jjj, ".", sep = ""), nx, fixed = TRUE)]
          names(xl[[j]][[jj]][[jjj]]) <- gsub(paste(j, ".", jj, ".", jjj, ".", sep = ""), "", names(xl[[j]][[jj]][[jjj]]), fixed = TRUE)
        }
      }
    }
  }
  xl
}


## Main bamlss().
bamlss <- function(formula, family = "gaussian", data = NULL, start = NULL, knots = NULL,
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit, contrasts = NULL,
  reference = NULL, transform = NULL, optimizer = NULL, sampler = NULL, samplestats = NULL, results = NULL,
  cores = NULL, sleep = NULL, combine = TRUE, model = TRUE, x = TRUE, light = FALSE, ...)
{
  ## Search for functions in family object.
  family <- bamlss.family(family, ...)
  if(!is.null(family$transform) & is.null(transform))
    transform <- family$transform
  if(!is.null(family$optimizer) & is.null(optimizer))
    optimizer <- family$optimizer
  if(!is.null(family$sampler) & is.null(sampler))
    sampler <- family$sampler
  if(!is.null(family$samplestats) & is.null(samplestats))
    samplestats <- family$samplestats
  if(!is.null(family$results) & is.null(results))
    results <- family$results

  ## Switch for light variant.
  if(light) {
    results <- FALSE
    samplestats <- FALSE
  }

  ## Setup all processing functions.
  foo <- list("transform" = transform, "optimizer" = optimizer,
    "sampler" = sampler, "samplestats" = samplestats, "results" = results)
  nf <- names(foo)
  default_fun <- c("no.transform", "bfit", "GMCMC", "samplestats", "results.bamlss.default")
  functions <- list()
  for(j in 1:length(foo)) {
    if(is.null(foo[[nf[j]]])) {
      foo[[nf[j]]] <- if(default_fun[j] != "no.transform") {
        get(default_fun[j], envir = asNamespace("bamlss"))
      } else FALSE
    }
    if(is.list(foo[[nf[j]]])) {
      args <- foo[[nf[j]]]
      fun <- default_fun[j]
      functions[[nf[j]]] <- function(x, ...) {
        args <- c(args, list(...))
        args$x <- x
        do.call(fun, args)
      }
    } else functions[[nf[j]]] <- foo[[nf[j]]]
    if(!is.function(functions[[nf[j]]])) {
      if(!is.logical(functions[[nf[j]]])) {
        stop(paste("argument", nf[j], "is not a function!"))
      } else {
        if(functions[[nf[j]]]) {
          functions[[nf[j]]] <- get(default_fun[j])
        }
      }
    }
  }

  ## Create the 'bamlss.frame'.
  bf <- match.call(expand.dots = TRUE)
  bf[c("transform", "optimizer", "sampler", "samplestats",
    "results", "cores", "sleep", "combine", "model", "x")] <- NULL
  bf[[1]] <- as.name("bamlss.frame")
  bf <- eval(bf, envir = parent.frame())

  ## Transform.
  if(is.function(functions$transform)) {
    tbf <- functions$transform(bf, ...)
    bf[names(tbf)] <- tbf
    rm(tbf)
  }

  ## Start optimizer.
  if(is.function(functions$optimizer)) {
    opt <- functions$optimizer(x = bf$x, y = bf$y, family = bf$family,
      start = start, weights = model.weights(bf$model.frame),
      offset = model.offset(bf$model.frame), ...)
    if(!is.list(opt)) {
      if(inherits(opt, "numeric")) {
        opt <- list("parameters" = drop(opt))
      } else stop("the optimizer should return the parameters as named numeric vector!")
    }
    if(is.null(opt$parameters))
      stop("the optimizer must return an element $parameters!")
    if(inherits(opt$parameters, "data.frame") | inherits(opt$parameters, "matrix")) {
      if(is.null(colnames(opt$parameters)))
        stop("the returned parameters must be a named numeric data.frame or matrix!")
    } else {
      if(is.null(names(opt$parameters)))
        stop("the returned parameters must be a named numeric vector!")
    }
    bf$parameters <- opt$parameters
    if(!is.null(opt$fitted.values))
      bf$fitted.values <- opt$fitted.values
    if(!is.null(opt$hessian))
      bf$hessian <- opt$hessian
    ne <- names(opt)
    bf$model.stats <- opt[ne[!(ne %in% c("parameters", "fitted.values", "hessian"))]]
    rm(opt)
  }

  ## Start sampling.
  if(is.function(functions$sampler)) {
    if(is.null(cores)) {
      bf$samples <- functions$sampler(x = bf$x, y = bf$y, family = bf$family,
        weights = model.weights(bf$model.frame),
        offset = model.offset(bf$model.frame),
        start = if(is.null(bf$parameters)) start else unlist(bf$parameters),
        hessian = bf$hessian, ...)
    } else {
      parallel_fun <- function(j) {
        if(j > 1 & !is.null(sleep)) Sys.sleep(sleep)
        functions$sampler(x = bf$x, y = bf$y, family = bf$family,
          weights = model.weights(bf$model.frame),
          offset = model.offset(bf$model.frame),
          start = if(is.null(bf$parameters)) start else unlist(bf$parameters),
          hessian = bf$hessian, ...)
      }
      bf$samples <- parallel::mclapply(1:cores, parallel_fun, mc.cores = cores)
    }
    if(!inherits(bf$samples, "mcmc")) {
      if(!is.null(bf$samples)) {
        if(is.list(bf$samples)) {
          bf$samples <- as.mcmc.list(lapply(bf$samples, as.mcmc))
        } else {
          bf$samples <- as.mcmc(bf$samples)
        }
      }
    }

    ## Process samples.
    bf$samples <- process.chains(bf$samples, combine)

    ## Optionally, compute more model stats from samples.
    if(is.function(functions$samplestats)) {
      ms <- functions$samplestats(samples = bf$samples, x = bf$x, y = bf$y, family = bf$family, ...)
      if(is.null(bf$model.stats)) {
        bf$model.stats <- list("sampler" = list())
        bf$model.stats$sampler[names(ms)] <- ms
      } else {
        bf$model.stats <- list("optimizer" = bf$model.stats, "sampler" = list())
        bf$model.stats$sampler[names(ms)] <- ms
      }
    } else {
      if(!is.null(bf$model.stats))
        bf$model.stats <- list("optimizer" = bf$model.stats)
    }
  } else {
    if(!is.null(bf$model.stats))
      bf$model.stats <- list("optimizer" = bf$model.stats)
  }

  ## Compute results.
  if(is.function(functions$results))
    bf$results <- try(functions$results(bf, bamlss = TRUE,  ...))

  ## Save the model frame?
  if(!model | light)
    bf$model.frame <- NULL

  ## Save 'x' master object?
  if(!x)
    bf$x <- NULL
  if(light) {
    if(!is.null(bf$x)) {
      for(j in names(bf$x)) {
        if(length(bf$x[[j]]$smooth.construct)) {
          for(i in seq_along(bf$x[[j]]$smooth.construct)) {
            bf$x[[j]]$smooth.construct[[i]][c("X", "S", "Xr", "Xf", "binning", "prior", "grad", "hess")] <- NULL
            bf$x[[j]]$smooth.construct[[i]][["xt"]][["binning"]] <- NULL
            if(!is.null(bf$x[[j]]$smooth.construct[[i]][["margin"]])) {
              for(jj in seq_along(bf$x[[j]]$smooth.construct[[i]][["margin"]])) {
                bf$x[[j]]$smooth.construct[[i]][["margin"]][[jj]][c("X", "S", "Xr", "Xf", "binning", "prior", "grad", "hess")] <- NULL
                bf$x[[j]]$smooth.construct[[i]][["margin"]][[jj]][["xt"]][["binning"]] <- NULL
              }
            }
          }
        }
      }
    }
    bf$y <- NULL
    bf$fitted.values <- NULL
    attr(bf$formula, ".Environment") <- NULL
    attr(bf$terms, ".Environment") <- NULL
  }

  bf$call <- match.call()
  class(bf) <- c("bamlss", "bamlss.frame", "list")
  attr(bf, "functions") <- functions

  bf
}


## Basic engine setup transformer.
bamlss.setup <- function(x, ...)
{
  list("x" = bamlss.engine.setup(x$x, ...))
}

## family extractor.
family.bamlss <- family.bamlss.frame <- function(object, ...)
{
  return(object$family)
}


## Extract all parameter names.
get.all.parnames <- function(x, rename.p = TRUE)
{
  pn <- names(parameters(if(inherits(x, "bamlss.frame")) x$x else x, list = FALSE))
  if(rename.p)
    pn <- gsub("s.model.matrix", "p.model.matrix", pn, fixed = TRUE)
  pn
}


## Model stats based on samples.
samplestats <- function(samples, x = NULL, y = NULL, family = NULL, logLik = FALSE, ...)
{
  if(inherits(samples, "bamlss")) {
    if(is.null(samples$samples))
      stop("no samples in 'bamlss' object!")
    x <- if(is.null(samples$x)) smooth.construct(samples) else samples$x
    y <- samples$y
    family <- samples$family
    samples <- samples$samples
  }
  what <- c("logLik", "logPost", "DIC", "pd")
  if(inherits(samples, "mcmc.list"))
    samples <- process.chains(samples)
  if(is.null(samples)) return(NULL)
  samples <- as.matrix(samples)
  sn <- colnames(samples)
  stats <- NULL
  taken <- what[what %in% sn]
  if(length(taken)) {
    what <- what[!(what %in% taken)]
    stats <- samples[, taken, drop = FALSE]
    if(logLik) {
      if("loglik" %in% tolower(taken))
        return(stats[, tolower(taken) == "loglik"])
    }
    stats <- as.list(apply(stats, 2, mean, na.rm = TRUE))
  }
  if(is.data.frame(y)) {
    if(ncol(y) < 2)
      y <- y[[1]]
  }
  if(length(what) | logLik) {
    pn <- get.all.parnames(x, rename.p = TRUE)
    pn <- gsub(".p.model.matrix.", ".p.", pn, fixed = TRUE)
    pn <- pn[pn %in% colnames(samples)]
    if(length(pn) & ("DIC" %in% what)) {
      samples <- samples[, pn, drop = FALSE]
      if(is.null(family$p2d) & is.null(family$p2logLik)) {
        nx <- names(x)
        par <- rep(list(0), length = length(x))
        names(par) <- nx
        mpar <- par
        for(i in nx)
          par[[i]] <- make.link2(family$links[i])$linkinv(.fitted.bamlss(i, x[[i]], samples))
        msamples <- matrix(apply(samples, 2, mean, na.rm = TRUE), nrow = 1)
        colnames(msamples) <- pn
        for(i in nx)
          mpar[[i]] <- make.link2(family$links[i])$linkinv(.fitted.bamlss(i, x[[i]], msamples))
        tpar <- mpar
        dev <- ll <- rep(NA, ncol(par[[1]]))
        for(j in 1:ncol(par[[1]])) {
          for(i in nx)
            tpar[[i]] <- par[[i]][, j]
          llt <- try(family$loglik(y, tpar), silent = TRUE)
          if(!inherits(llt, "try-error")) {
            ll[j] <- llt
            dev[j] <- -2 * ll[j]
          }
        }
        if(logLik)
          return(ll)
        ll <- try(family$loglik(y, mpar), silent = TRUE)
      } else {
        mpar <- apply(samples, 2, mean, na.rm = TRUE)
        ll <- try(apply(samples, 1, function(x) {
          names(x) <- colnames(samples)
          if(is.null(family$p2logLik)) {
            sum(family$p2d(x, log = TRUE), na.rm = TRUE)
          } else family$p2logLik(x)
        }), silent = TRUE)
        if(inherits(ll, "try-error")) {
          warning("no DIC, cannot evaluate the $p2d() function in 'bamlss' family object!")
        } else {
          if(logLik)
            return(ll)
          dev <- -2 * ll
          ll <- try(if(is.null(family$p2logLik)) {
              sum(family$p2d(mpar, log = TRUE), na.rm = TRUE)
            } else family$p2logLik(mpar), silent = TRUE)
        }
      }
      if(!inherits(ll, "try-error") & !all(!is.finite(dev))) {
        mdev <- -2 * ll
        pd <- mean(dev, na.rm = TRUE) - mdev
        DIC <- mdev + 2 * pd
        if(is.null(stats))
          stats <- list()
        stats$DIC <- DIC
        stats$pd <- pd
      } else {
        warning("no DIC, cannot evaluate the $loglik() function in 'bamlss' family object!")
      }
    }
  }
  
  return(stats)
}


#### -----------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------
#### -----------------------------------------------------------------------------------------------
## Could be interesting: http://people.duke.edu/~neelo003/r/
##                       http://www.life.illinois.edu/dietze/Lectures2012/
#########################
## (2) Engine stacker. ##
#########################
#stacker <- function(x, optimizer = bfit0, sampler = samplerJAGS, ...)
#{
#  if(is.function(optimizer) | is.character(optimizer))
#    optimizer <- list(optimizer)
#  if(is.integer(sampler) | is.numeric(sampler)) {
#    n.samples <- as.integer(sampler)
#    sampler <- function(x, ...) { null.sampler(x, n.samples = n.samples) }
#  }
#  if(is.null(sampler))
#    sampler <- null.sampler
#  if(is.function(sampler) | is.character(sampler))
#    sampler <- list(sampler)
#  if(length(optimizer)) {
#    for(j in optimizer) {
#      if(is.character(j)) j <- eval(parse(text = j))
#      if(!is.function(j)) stop("the optimizer must be a function!")
#      x <- j(x, ...)
#    }
#  }
#  if(length(sampler)) {
#    for(j in sampler) {
#      if(is.character(j)) j <- eval(parse(text = j))
#      if(!is.function(j)) stop("the sampler must be a function!")
#      x <- j(x, ...)
#    }
#  }

#  x
#}


"[.bamlss" <- function(x, ...) {
  rval <- NextMethod("[", ...)
  mostattributes(rval) <- attributes(x)
  rval
}


## Create the model.frame.
bamlss.model.frame <- function(formula, data, family = gaussian_bamlss(),
  weights = NULL, subset = NULL, offset = NULL, na.action = na.omit,
  specials = NULL, contrasts.arg = NULL, drop.unused.levels = TRUE, ...)
{
  if(!missing(data)) {
    if(is.character(data)) {
      if(!file.exists(data))
        stop("data path is not existing!")
      ## data <- read.table.ffdf(file = data,
      ##   na.strings = "", header = TRUE, sep = ",")
      ## return(data)
      ## FIXME: ff data.frames!
      data <- read.table(file = data, header = TRUE, ...)
    }
  } else data <- NULL
  if(inherits(formula, "bamlss.frame") | inherits(formula, "bamlss")) {
    if(!is.null(formula$model.frame))
      return(formula$model.frame)
    fcall <- formula$call
    fcall[[1L]] <- quote(bamlss.model.frame)
    env <- environment(formula$formula)
    if(is.null(env))
      env <- parent.frame()
    return(eval(fcall, env))
  } else {
    family <- bamlss.family(family, ...)
    formula <- bamlss.formula(formula, family, specials, env = parent.frame())
  }

  if(is.null(na.action))
    na.action <- get(getOption("na.action"))
  if(missing(data))
    data <- environment(formula)

  if(!is.data.frame(data))
    data <- as.data.frame(data)

  ## Make fake "Formula" object.
  fF <- make_fFormula(formula)
  attr(fF, ".Environment") <- environment(formula)

  ## Resulting terms object.
  mterms <- terms(formula(fF), data = data)

  ## Set up the model.frame.
  data <- list(formula = fF, data = data, subset = subset,
    na.action = na.action, drop.unused.levels = drop.unused.levels, ...)
  data <- do.call("model.frame", data)
  rownames(data) <- NULL

  ## Code from stats model.matrix()
  contr.funs <- as.character(getOption("contrasts"))
  namD <- names(data)
  for(i in namD) {
    if(is.character(data[[i]])) 
      data[[i]] <- factor(data[[i]])
  }
  isF <- vapply(data, function(x) is.factor(x) || is.logical(x), NA)
  isF[attr(mterms, "response")] <- FALSE
  isOF <- vapply(data, is.ordered, NA)
  for(nn in namD[isF]) {
    if(is.null(attr(data[[nn]], "contrasts"))) {
      contrasts(data[[nn]]) <- contr.funs[1 + isOF[nn]]
    }
  }
  if(!is.null(contrasts.arg) && is.list(contrasts.arg)) {
    if(is.null(namC <- names(contrasts.arg))) 
      stop("invalid 'contrasts' argument")
    for(nn in namC) {
      if(is.na(ni <- match(nn, namD))) 
        warning(gettextf("variable '%s' is absent, its contrast will be ignored", nn), domain = NA)
      else {
        ca <- contrasts.arg[[nn]]
        if(is.matrix(ca)) {
          contrasts(data[[ni]], ncol(ca)) <- ca
        } else {
          contrasts(data[[ni]]) <- contrasts.arg[[nn]]
        }
      }
    }
  }

  ## Process weights and offset.
  if(!is.null(weights)) {
    if(!is.list(weights)) {
      weights <- rep(list(weights), length = length(family$names))
      names(weights) <- family$names
    }
    weights <- do.call("cbind", weights)
    colnames(weights) <- names(formula)[1:ncol(weights)]
    if(!is.null(subset)) {
      weights <- if(!is.logical(subset)) {
        weights[subset, , drop = FALSE]
      } else subset(weights, subset)
    }
    if(nrow(weights) < 2)
      weights <- do.call(rbind, replicate(nrow(data), weights, simplify = FALSE))
    for(j in 1:ncol(weights))
      weights[weights[, j] == 0, j] <- .Machine$double.eps
    data[["(weights)"]] <- weights
  }
  if(!is.null(offset)) {
    if(!is.list(offset)) {
      offset <- rep(list(offset), length = length(family$names))
      names(offset) <- family$names
    }
    offset <- do.call("cbind", offset)
    colnames(offset) <- names(formula)[1:ncol(offset)]
    if(!is.null(subset)) {
      offset <- if(!is.logical(subset)) {
        offset[subset, , drop = FALSE]
      } else subset(offset, subset)
    }
    if(nrow(offset) < 2)
      offset <- do.call(rbind, replicate(nrow(data), offset, simplify = FALSE))
    data[["(offset)"]] <- offset
  }

  ## Remove inf values.
  data <- rm_infinite(data)

  ## Assign terms object.
  attr(data, "terms") <- mterms

  ## Check response.
  if(!is.null(family$valid.response)) {
    family$valid.response(model.response(data))
  }

  data
}


## Remove Inf values from data.
rm_infinite <- function(x) {
  if(is.null(dim(x))) return(x)
  if(ncol(x) > 0) {
    for(j in 1:ncol(x)) {
      if(any(class(x[, j]) %in% c("numeric", "integer"))) {
        if(any(!is.finite(x[, j]))) {
          warning("infinite values in data, removing these observations in model.frame!")
          x <- x[is.finite(x[, j]), ]
        }
      }
    }
  }
  x
}


## Parse families and get correct family object, depending on type.
bamlss.family <- function(family, type = "bamlss", ...)
{
  family <- if(is.function(family)) family() else {
    if(is.character(family)) {
      if(!is.null(type)) {
        if(!grepl("gF(", family, fixed = TRUE) & !grepl("gF2(", family, fixed = TRUE))
          if(!grepl(type, family))
            family <- paste(family, type, sep = "_")
      }
      family <- eval(parse(text = family[1]))
      if(is.function(family))
        family()
      else family
    } else family
  }
  if(inherits(family, "gamlss.family"))
    family <- tF(family, ...)
  if(!inherits(family, "family.bamlss")) {
    if(!is.character(family)) {
      if(is.null(family$family)) stop("family is specified wrong, no family name available!")
      family <- family$family
    }
    txt <- paste(family, type, sep = if(!is.null(type)) "_" else "")
    txt <- gsub("bamlss.bamlss", "bamlss", txt, fixed = TRUE)
    family <- eval(parse(text = txt[1]))
    family <- family()
  }
  if(is.null(family)) family <- list()

  family
}


complete.bamlss.family <- function(family)
{
  if(is.null(names(family$links)))
    names(family$links) <- family$names

  linkinv <- linkfun <- list()
  for(j in family$names) {
    link <- make.link2(family$links[j])
    linkinv[[j]] <- link$linkinv
    linkfun[[j]] <- link$linkfun
  }

  if(is.null(family$map2par)) {
    family$map2par <- function(eta) {
      if(inherits(eta[[1L]], "ff")) {
        for(j in family$names)
          eta[[j]] <- ff_eval(eta[[j]], FUN = function(x) { linkinv[[j]](x) },
            lower = c(-Inf, -10), upper = c(Inf, 10))
      } else {
        for(j in family$names) {        
          eta[[j]] <- linkinv[[j]](eta[[j]])
          eta[[j]][is.na(eta[[j]])] <- 0
          if(any(jj <- eta[[j]] == Inf))
            eta[[j]][jj] <- 10
          if(any(jj <- eta[[j]] == -Inf))
            eta[[j]][jj] <- -10
        }
      }
      return(eta)
    }
  }

  if(is.null(family$mu)) {
    family$mu <- function(par) { make.link2(family$links[1])$linkinv(par[[1]]) }
  }

  if(is.null(family$loglik)) {
    if(!is.null(family$d))
      family$loglik <- function(y, par, ...) { sum(family$d(y, par, log = TRUE), na.rm = TRUE) }
  }

  err01 <- .Machine$double.eps^(1/3)
  err02 <- err01 * 2
  err11 <- .Machine$double.eps^(1/4)
  err12 <- err11 * 2

  if(is.null(family$score))
    family$score <- list()
  for(i in family$names) {
    if(is.null(family$score[[i]])) {
      fun <- c(
        "function(y, par, ...) {",
        paste("  eta <- linkfun[['", i, "']](par[['", i, "']]);", sep = ""),
        paste("  par[['", i, "']] <- linkinv[['", i, "']](eta + err01);", sep = ""),
        "  d1 <- family$d(y, par, log = TRUE);",
        paste("  par[['", i, "']] <- linkinv[['", i, "']](eta - err01);", sep = ""),
        "  d2 <- family$d(y, par, log = TRUE);",
        "  return((d1 - d2) / err02)",
        "}"
      )
      family$score[[i]] <- eval(parse(text = paste(fun, collapse = "")))
      attr(family$score[[i]], "dnum") <- TRUE
    }
  }

  if(is.null(family$hess))
    family$hess <- list()
  for(i in family$names) {
    if(is.null(family$hess[[i]])) {
      fun <- if(!is.null(attr(family$score[[i]], "dnum"))) {
        c(
          "function(y, par, ...) {",
          paste("  eta <- linkfun[['", i, "']](par[['", i, "']]);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta + err11);", sep = ""),
          paste("  d1 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta - err11);", sep = ""),
          paste("  d2 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          "  return(-1 * (d1 - d2) / err12)",
          "}"
        )
      } else {
        c(
          "function(y, par, ...) {",
          paste("  eta <- linkfun[['", i, "']](par[['", i, "']]);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta + err01);", sep = ""),
          paste("  d1 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta - err01);", sep = ""),
          paste("  d2 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          "  return(-1 * (d1 - d2) / err02)",
          "}"
        )
      }
      family$hess[[i]] <- eval(parse(text = paste(fun, collapse = "")))
    }
  }

  return(family)
}


## Formula to list().
as.list.Formula <- function(x)
{
  if(!inherits(x, "Formula"))
    x <- as.Formula(x)
  env <- environment(x)
  lhs <- attr(x, "lhs")
  rhs <- attr(x, "rhs")
  nl <- length(lhs)
  nr <- length(rhs)
  if(nl < nr)
    lhs <- c(lhs, rep(list(NA), length = nr - nl))
  if(nr < nl)
    rhs <- c(rhs, rep(list(1), length = nl - nr))
  x <- mapply(c, lhs, rhs, SIMPLIFY = FALSE)
  formula <- list()
  for(i in seq_along(x)) {
    check <- inherits(x[[i]][[1]], "call") | inherits(x[[i]][[1]], "name")
    f <- if(check) {
      as.call(c(as.symbol("~"), x[[i]]))
    } else as.call(c(as.symbol("~"), x[[i]][[2]]))
    formula[[i]] <- eval(f, envir = env)
    attr(formula[[i]], ".Environment") <- NULL
  }
  environment(formula) <- env
  formula
}


## Special formula parser, can deal with multi parameter models
## and hierarchical structures.
bamlss.formula <- function(formula, family = NULL, specials = NULL, env = NULL, ...)
{
  if(is.null(specials))
    specials <- c("s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "la", "n")
  if(inherits(formula, "bamlss.formula"))
    return(formula)
  if(!is.list(formula)) {
    if(!inherits(formula, "Formula"))
      formula <- as.Formula(formula)
    if(inherits(formula, "Formula"))
      formula <- as.list.Formula(formula)
  }
  if(!is.null(family))
    family <- bamlss.family(family, ...)
  if(!is.list(formula)) formula <- list(formula)
  if(!length(formula)) stop("formula is specified wrong!")
  
  if(is.null(env))
    env <- get_formula_envir(formula)

  complete_formula <- function(formula) {
    if(!is.null(family)) {
      if(length(formula) < length(family$names))
        formula <- c(formula, rep(list(), length = length(family$names) - length(formula)))
    }
    fn <- NULL
    for(j in seq_along(formula)) {
      ft <- if(!inherits(formula[[j]], "formula")) formula[[j]][[1]] else formula[[j]]
      if(!is.null(ft)) {
        yok <- attr(terms(formula(as.Formula(ft), rhs = FALSE)), "response") > 0
        fn <- c(fn, if(yok) all.vars(ft)[1] else NULL)
      }
    }
    fn[fn %in% c("1", "-1")] <- NA

    nas <- which(is.na(fn))
    if(!is.null(family)) {
      if(length(nas))
        fn[nas] <- family$names[nas]
      if(is.null(family$names))
        family$names <- NA
      if(!all(is.na(family$names[1:length(fn)])))
        fn <- family$names[1:length(fn)]
      else
        family$names[1:length(fn)] <- fn
    } else fn[nas] <- paste("par", 1:length(fn[nas]), sep = ".")

    if(is.null(family)) {
      if(length(fn) < length(formula)) {
        k <- length(formula) - length(fn)
        if(k > 1)
          fn <- c(fn, paste("?par", 1:k, sep = ""))
        else
          fn <- c(fn, "?par")
      }
    }
    names(formula) <- fn
    if(!is.null(family)) {
      if(any(i <- is.na(names(formula))))
        names(formula)[i] <- family$names[i]
      for(j in family$names) {
        if(is.null(formula[[j]])) {
          formula[[j]] <- as.formula(paste(j, "1", sep = " ~ "), env = NULL)
        }
      }
    }
    for(j in seq_along(formula)) {
      if(!inherits(formula[[j]], "formula")) {
        if(is.null(names(formula[[j]])))
          names(formula[[j]]) <- paste("h", 1:length(formula[[j]]), sep = "")
      } else {
        attr(formula[[j]], ".Environment") <- NULL
      }
    }
    formula
  }
  formula <- formula_and(formula)
  formula <- formula_at(formula)
  formula <- complete_formula(formula_hierarchical(formula))
  formula <- formula_extend(formula, family, specials)

  environment(formula) <- env
  class(formula) <- c("bamlss.formula", "list")

  formula
}


## Formula environment.
get_formula_envir <- function(formula)
{
  env <- environment(formula)
  if(is.null(env)) {
    get_env <- function(x) {
      if(inherits(x, "list")) {
        env <- NULL
        for(j in x)
          env <- c(env, get_env(j))
        return(env)
      } else return(environment(x))
    }
    env <- get_env(formula)
  }
  if(is.null(env)) env <- .GlobalEnv
  if(is.list(env))
    env <- env[[1]]
  return(env)
}


## Process categorical responses.
bamlss.formula.cat <- function(formula, family, data, reference)
{
  env <- environment(formula)

  rn <- y <- NULL
  for(j in seq_along(formula)) {
    ft <- if(!inherits(formula[[j]]$formula, "formula")) {
      formula[[j]][[1]]$formula
    } else formula[[j]]$formula
    yok <- attr(terms(formula(as.Formula(ft), rhs = FALSE)), "response") > 0
    if(yok)
      rn <- c(rn, all.vars(ft)[1])
  }
  rn2 <- rn[rn %in% names(data)]
  cat <- !is.null(family$cat) & (length(rn2) > 1)
  if(is.factor(data[[rn2[1]]]) | cat) {
    if((nlevels(data[[rn2[1]]]) > 2) | cat) {
      if(!cat | is.factor(data[[rn2[1]]])) {
        ft <- as.formula(paste("~ -1 +", rn2[1]), env = NULL)
        y <- model.matrix(ft, data = data)
        colnames(y) <- rmf(gsub(rn2[1], "", colnames(y), fixed = TRUE))
        if(is.null(reference)) {
          ty <- table(data[[rn2[1]]])
          reference <- c(names(ty)[ty == max(ty)])[1]
        } else {
          ld <- levels(data[[rn2[1]]])
          reference <- ld[match(gsub(rn2[1], "", reference), ld)]
        }
        if(is.na(reference))
          stop(paste("cannot find reference category within response levels!"))
        reference <- rmf(reference)
        ylevels <- rmf(levels(data[[rn2[1]]]))
        ylevels <- ylevels[ylevels != reference]
        y <- y[, colnames(y) %in% ylevels, drop = FALSE]
      } else {
        y <- data[, rn2]
        ylevels <- colnames(y)
        reference <- ""
      }
      if(length(formula) < ncol(y)) {
        formula <- c(formula, rep(formula, length = ncol(y) - length(formula)))
      }
      if(!(names(formula)[[1]] %in% colnames(y))) {
        names(formula)[[1]] <- colnames(y)[1]
        ft <- if(!inherits(formula[[1]]$formula, "formula")) {
          formula[[1]][[1]]$formula
        } else formula[[1]]$formula
        env <- environment(ft)
        ft <- update(ft, as.formula(paste(colnames(y)[1], ".", sep = "~"), env = NULL))
        environment(ft) <- env
        if(!inherits(formula[[1]]$formula, "formula")) {
          formula[[1]][[1]]$formula <- ft
          formula[[1]][[1]]$response <- colnames(y)[1]
        } else {
          formula[[1]]$formula <- ft
          formula[[1]]$response <- colnames(y)[1]
        }
        ft <- if(!inherits(formula[[1]]$formula, "formula")) {
          formula[[1]][[1]]$fake.formula
        } else formula[[1]]$fake.formula
        ft <- update(ft, as.formula(paste(colnames(y)[1], ".", sep = "~"), env = NULL))
        if(!inherits(formula[[1]]$formula, "formula")) {
          formula[[1]][[1]]$fake.formula <- ft
        } else {
          formula[[1]]$fake.formula <- ft
        }
      }
      if(length(i <- !(names(formula) %in% ylevels))) {
        k <- 1
        ynot <- ylevels[!(ylevels %in% names(formula))]
        for(j in which(i)) {
          names(formula)[[j]] <- ynot[k]
          ft <- if(!all(c("formula", "fake.formula") %in% names(formula[[j]]))) {
            formula[[j]][[1]]$formula
          } else formula[[j]]$formula
          env <- environment(ft)
          ft <- update(ft, as.formula(paste(ynot[k], "~ ."), env = NULL))
          environment(ft) <- env
          if(!inherits(formula[[j]]$formula, "formula")) {
            formula[[j]][[1]]$formula <- ft
            formula[[j]][[1]]$response <- ynot[k]
          } else {
            formula[[j]]$formula <- ft
            formula[[j]]$response <- ynot[k]
          }
          attr(formula[[j]], "name") <- ynot[k]
          k <- k + 1
        }
      }
    }
  }

  rval <- if(!is.null(y)) {
    class(formula) <- "bamlss.formula"
    environment(formula) <- env
    list("formula" = formula, "ylevels" = ylevels, "reference" = reference)
  } else NULL
  rval
}


## Make "Formula" object from fake.formulas.
make_fFormula <- function(formula)
{
  fF <- NULL
  for(j in seq_along(formula)) {
    if(!all(c("formula", "fake.formula") %in% names(formula[[j]]))) {
      for(i in seq_along(formula[[j]]))
        fF <- c(fF, formula[[j]][[i]]$fake.formula)
    } else {
      fF <- c(fF, formula[[j]]$fake.formula)
    }
  }
  fF <- do.call("as.Formula", fF)
  fF
}


all.vars.formula <- function(formula, lhs = TRUE, rhs = TRUE, specials = NULL, intercept = FALSE, type = 1)
{
  env <- environment(formula)
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "la", "n"))
  tf <- terms(formula, specials = specials, keep.order = TRUE)
  ## sid <- unlist(attr(tf, "specials")) - attr(tf, "response")
  tl <- attr(tf, "term.labels")
  sid <- NULL
  for(j in specials)
    sid <- c(sid, grep2(paste(j, "(", sep = ""), tl, fixed = TRUE))
  if(length(sid))
    sid <- sort(unique(sid))
  vars <- NULL
  if(rhs) {
    if(length(sid)) {
      vars <- tl[-sid]
      if(!length(vars))
        vars <- NULL
      for(j in tl[sid]) {
        tcall <- parse(text = j)[[1]]
        tcall[c("k","fx","bs","m","xt","id","sp","pc","d","mp","np")] <- NULL
        tcall <- eval(tcall)
        vars <- c(vars, tcall$term)
        if(!is.null(tcall$by)) {
          if(tcall$by != "NA")
            vars <- c(vars, tcall$by)
        }
      }
    } else {
      vars <- tl
    }
  }
  if(lhs & (attr(tf, "response") > 0))
    vars <- c(vars, response.name(formula, keep.functions = TRUE))
  if(intercept & (attr(tf, "intercept") > 0))
    vars <- c("1", vars)
  vars <- vars[vars != "."]
  if(length(vars) < 1)
    vars <- NULL
  if(any(i <- grep(":", vars, fixed = TRUE))) {
    dv <- unlist(strsplit(vars[i], ":", fixed = TRUE))
    vars <- c(vars[-i], dv)
  }
  vars <- unique(vars)
  if(type == 1)
    vars <- all.vars(as.formula(paste("~", paste(vars, collapse = "+")), env = NULL))
  unique(vars)
}


## From nlme.
splitFormula <- function(form, sep = "+") 
{
  if(inherits(form, "formula") || mode(form) == "call" && 
    form[[1]] == as.name("~")) 
    return(splitFormula(form[[length(form)]], sep = sep))
  if(mode(form) == "call" && form[[1]] == as.name(sep)) 
    return(do.call("c", lapply(as.list(form[-1]), splitFormula, 
      sep = sep)))
  if(mode(form) == "(") 
    return(splitFormula(form[[2]], sep = sep))
  if(length(form) < 1) 
    return(NULL)
  list(stats::asOneSidedFormula(form))
}


terms.formula2 <- function(formula, specials, keep.order = TRUE, ...)
{
  fs <- splitFormula(formula, sep = c("+", "-"))
  tl <- rep("", length = length(fs))
  adds <- NULL
  for(j in seq_along(fs)) {
    tlj <- attr(terms(fs[[j]], specials = specials), "term.labels")
    if(length(tlj))
      tl[j] <- tlj[1]
    if(length(tlj) > 1)
      adds <- c(adds, tlj[-1])
  }
  tl <- c(tl, adds)
  tl <- tl[tl != ""]
  if(!length(tl))
    tl <- ""
  attr(formula, "term.labels") <- tl
  formula
}


all.labels.formula <- function(formula, specials = NULL, full.names = FALSE)
{
  env <- environment(formula)
  specials <- unique(c("s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "la", "n", specials))
  tf <- terms.formula2(formula, specials = specials, keep.order = FALSE)
  ## sid <- unlist(attr(tf, "specials")) - attr(tf, "response")
  tl <- attr(tf, "term.labels")
  sid <- NULL
  for(j in specials)
    sid <- c(sid, grep2(paste(j, "(", sep = ""), tl, fixed = TRUE))
  if(length(sid))
    sid <- sort(unique(sid))
  labs <- NULL

  if(length(sid)) {
    labs <- tl[-sid]
    if(full.names & length(labs))
      labs <- paste("p", labs, sep = ".")
    if(!length(labs))
      labs <- NULL
    else
      tl[-sid] <- labs
    tl[sid] <- gsub(" ", "", tl[sid])
    for(j in sid) {
      tcall <- parse(text = tl[j])[[1]]
      tcall[c("k","fx","bs","m","xt","sp","pc","d","mp","np")] <- NULL
      tcall <- eval(tcall)
      tl[j] <- gsub(" ", "", tcall$label)
      if(!is.null(tcall$by)) {
        if(tcall$by != "NA") {
          if(!grepl("by=", tl[j], fixed = TRUE)) {
            tlt <- strsplit(tl[j], "")[[1]]
            tlt <- paste(tlt[1:(length(tlt) - 1)], collapse = "")
            tl[j] <- paste(tlt, ",by=", tcall$by, ")", sep = "")
          }
        }
      }
    }
    if(full.names)
      tl[sid] <- paste("s", tl[sid], sep = ".")
    labs <- tl
  } else labs <- if(full.names) paste("p", tl, sep = ".") else tl
  unique(labs)
}


fake.formula <- function(formula, lhs = TRUE, rhs = TRUE, specials = NULL)
{
  if(all(!lhs & !rhs))
    return(0 ~ 0)
  if(all(rhs)) {
    f <- paste(all.vars.formula(formula, lhs = FALSE, rhs = TRUE, specials, intercept = TRUE, type = 2), collapse = "+")
    if(f == "") f <- "-1"
  }
  if(all(lhs))
    f <- paste(all.vars.formula(formula, lhs = TRUE, rhs = FALSE, type = 2), "~", if(!is.null(f)) f else 0)
  else
    f <- paste("~", if(!is.null(f)) f else 0)
  f <- as.formula(f, env = NULL)
  f
}


## Extend formula by a fake formula with all variables
## to compute a model.frame, create smooth objects.
formula_extend <- function(formula, family, specials = NULL)
{
  if(is.list(formula)) {
    for(j in seq_along(formula))
      formula[[j]] <- formula_extend(formula[[j]], family, specials)
    return(formula)
  } else {
    rn <- response.name(formula)
    ff <- fake.formula(formula, lhs = !(rn %in% family$names), specials = specials)
    return(list("formula" = formula, "fake.formula" = ff))
  }
}


## Get response name.
response.name <- function(formula, hierarchical = TRUE, keep.functions = FALSE, na.rm = FALSE)
{
  rn <- NA
  if(inherits(formula, "bamlss.frame")) {
    if(!is.null(formula$formula)) {
      if(!is.null(attr(formula$formula, "response.name")))
        return(attr(formula$formula, "response.name"))
    }
    formula <- terms(model.frame(formula))
  }
  if(!is.null(attr(formula, "terms")))
    formula <- attr(formula, "terms")
  if(inherits(formula, "formula")) {
    f <- as.Formula(formula)
    f <- formula(f, lhs = TRUE, rhs = FALSE)
    if(keep.functions) {
      cf <- as.character(formula)
      rn <- if(length(cf) < 3) character(0) else cf[2]
    } else {
      rn <- all.vars(f)
    }
  } else {
    if(inherits(formula, "list")) {
      rn <- NULL
      for(i in seq_along(formula)) {
        if(is.null(formula[[i]]$formula) & inherits(formula[[i]], "list")) {
          for(j in seq_along(formula[[i]])) {
            if(!hierarchical & (j > 1)) {
              next
            } else {
              tf <- if(is.null(formula[[i]][[j]]$formula)) {
                formula[[i]][[j]]
              } else formula[[i]][[j]]$formula
              rn <- c(rn, response.name(tf, keep.functions = keep.functions))
            }
          }
        } else {
          rn <- c(rn , response.name(formula[[i]]$formula, keep.functions = keep.functions))
        }
      }
    }
  }
  if(!length(rn))
    rn <- NA
  if(na.rm)
    rn <- rn[!is.na(rn)]
  rn
}

## Search and process "&"
formula_and <- function(formula)
{
  if(nol <- !is.list(formula))
    formula <- list(formula)
  for(j in seq_along(formula)) {
    if(!inherits(formula[[j]], "formula")) {
      formula[[j]] <- formula_and(formula[[j]])
    } else {
      ft <- deparse(formula[[j]])
      if(any(grep("&", ft, fixed = TRUE))) {
        formula[[j]] <- as.list(strsplit(ft, "&", fixed = TRUE)[[1]])
        for(i in seq_along(formula[[j]])) {
          if(!any(grepl("~", formula[[j]][[i]], fixed = TRUE)))
            formula[[j]][[i]] <- paste("~", formula[[j]][[i]])
          formula[[j]][[i]] <- as.formula(formula[[j]][[i]], env = NULL)
        }
      }
    }
  }
  if(nol) {
    names(formula) <- response.name(formula[[1]][[1]])
    if(inherits(formula[[1]], "formula"))
      formula <- formula[[1]]
  }
  formula
}

## Search and process "@"
formula_at <- function(formula)
{
  if(nol <- !is.list(formula))
    formula <- list(formula)
  for(j in seq_along(formula)) {
    if(!inherits(formula[[j]], "formula")) {
      formula[[j]] <- formula_at(formula[[j]])
    } else {
      ft <- deparse(formula[[j]])
      if(any(grep("@", ft, fixed = TRUE))) {
        formula[[j]] <- strsplit(ft, "@", fixed = TRUE)[[1]]
        control <- formula[[j]][-1]
        formula[[j]] <- as.formula(formula[[j]][1], env = NULL)
        control <- gsub(":", "=", control)
        if(any(grepl("+", control)))
          control <- strsplit(control, "+", fixed = TRUE)[[1]]
        control <- gsub("^ +", "", control)
        control <- gsub(" +$", "", control)
        attr(formula[[j]], "control") <- gsub("using=", "using ", control, fixed = TRUE)
      }
    }
  }
  if(nol) {
    names(formula) <- response.name(formula[[1]][[1]])
    if(inherits(formula[[1]], "formula"))
      formula <- formula[[1]]
  }
  formula
}

formula_rm_at <- function(formula)
{
  ctr <- attr(formula, "control")
  attr(formula, "control") <- NULL
  if(isf <- !is.character(formula)) {
    formula <- deparse(formula)
  }
  if(any(grepl("@", formula)))
    formula <- strsplit(formula, "@")[[1]][1]
  if(!isf) {
    formula <- as.formula(formula, env = NULL)
    formula <- deparse(formula)
  }
  if(isf) {
    formula <- as.formula(formula, env = NULL)
  }
  if(!is.null(ctr)) attr(formula, "control") <- ctr
  formula
}

## Hierarchical formulae.
formula_hcheck <- function(formula)
{
  if(!is.list(formula))
    return(formula)
  check <- vector(mode = "list", length = length(formula))
  for(j in seq_along(formula)) {
      for(i in seq_along(formula)) {
        if(j != i) {
          fi <- if(!is.list(formula[[i]])) list(formula[[i]]) else formula[[i]]
          rnj <- response.name(formula[[j]], keep.functions = TRUE)
          for(jj in seq_along(fi)) {
            av <- all.vars(fi[[jj]])
            rn <- response.name(fi[[jj]])
            if(!any(is.na(rn))) {
              if(length(rn[is.na(rn)]))
                av <- av[av != rn[is.na(rn)]]
            }
            if(!has_dot(fi[[jj]])) {
              if(attr(terms(fi[[jj]]), "intercept") < 1) {
                av <- c(av, "-1")
              }
            }
            if(any(av %in% rnj)) {
              check[[j]] <- c(check[[j]], i)
            }
          }
        }
      }
  }
  check
}

formula_insert <- function(from, to, formula)
{
  nf <- names(formula)
  hm <- sapply(to, max)
  o <- order(hm, decreasing = TRUE)
  from <- from[o]
  to <- to[o]
  for(j in seq_along(from)) {
    for(i in seq_along(to[[j]])) {
      formula[[to[[j]][i]]] <- c(formula[[to[[j]][i]]], formula[[from[j]]])
    }
  }
  formula <- formula[take <- !(1:length(formula) %in% from)]
  names(formula) <- nf[take]
  formula
}

formula_hierarchical <- function(formula)
{
  if(!is.list(formula))
    return(formula)
  j <- formula_hcheck(formula)
  while(any(!sapply(j, is.null))) {
    i <- which(!sapply(j, is.null))
    formula <- formula_insert(i, j[i], formula)
    j <- formula_hcheck(formula)
  }
  formula
}


## Transform smooth terms to mixed model representation.
randomize <- function(x)
{
  if(bframe <- inherits(x, "bamlss.frame")) {
    if(is.null(x$x))
      stop("no 'x' object to randomize in 'bamlss.frame'!")
    x <- x$x
  }

  rand_fun <- function(x)
  {
    if(m <- length(x$smooth.construct)) {
      for(j in 1:m) {
        if(!inherits(x$smooth.construct[[j]], "no.mgcv")) {
          if(is.null(x$smooth.construct[[j]]$rand) & is.null(x$smooth.construct[[j]]$Xf)) {
            vnames <- x$smooth.construct[[j]]$term
            if(x$smooth.construct[[j]]$by != "NA")
              vnames <- c(vnames, x$smooth.construct[[j]]$by)
            tmp <- smooth2random(x$smooth.construct[[j]], vnames = vnames, type = 2)
            if(is.null(x$smooth.construct[[j]]$xt$nolin))
              x$smooth.construct[[j]]$Xf <- tmp$Xf
#          if(inherits(x$smooth.construct[[j]], "random.effect")) {
#            tmp$rand$Xr[tmp$rand$Xr > 0] <- 1
#            tmp$rand$Xr <- scale(tmp$rand$Xr)
#            tmp$trans.D <- rep(1, ncol(tmp$rand$Xr))
#            tmp$trans.U <- diag(1, ncol(tmp$rand$Xr))
#          }
            x$smooth.construct[[j]]$rand <- tmp$rand
            x$smooth.construct[[j]]$trans.D <- tmp$trans.D
            x$smooth.construct[[j]]$trans.U <- tmp$trans.U
            if(!is.null(x$smooth.construct[[j]]$state$parameters)) {
              b2 <- get.par(x$smooth.construct[[j]]$state$parameters, "b")
              if(!is.null(x$smooth.construct[[j]]$trans.U))
                b2 <- solve(x$smooth.construct[[j]]$trans.U) %*% b2
              b2 <- drop(b2 / x$smooth.construct[[j]]$trans.D)
              x$smooth.construct[[j]]$state$parameters <- set.par(x$smooth.construct[[j]]$state$parameters, b2, "b")
            }
          }
        }
      }
    }
    x
  }

  elmts <- c("formula", "fake.formula")
  for(j in seq_along(x)) {
    if(!all(elmts %in% names(x[[j]]))) {
      for(i in seq_along(x[[j]]))
        x[[j]][[i]] <- rand_fun(x[[j]][[i]])
    } else x[[j]] <- rand_fun(x[[j]])
  }
 
  if(bframe) {
    return(list("x" = x))
  } else {
    return(x)
  }
}


## Combine sample chains.
process.chains <- function(x, combine = TRUE, drop = FALSE)
{
  if(is.null(x)) return(NULL)
  if(!is.list(x))
    x <- list(x)
  n <- sapply(x, nrow)
  if((length(unique(n)) > 1) & combine) {
    x <- lapply(x, as.matrix)
    x <- list(as.mcmc(do.call("rbind", x)))
  }
  model.specs <- attr(x[[1]], "model.specs")
  if(inherits(x[[1]], "mcmc.list")) {
    x <- as.mcmc.list(do.call("c", x))
  } else {
    stopifnot(inherits(x[[1]], "mcmc"))
    x <- as.mcmc.list(x)
  }
  if(combine) {
    x <- do.call("rbind", x)
    x <- as.mcmc.list(list(as.mcmc(x)))
  }
  if(drop & (length(x) < 2))
    x <- x[[1]]
  attr(x, "model.specs") <- model.specs
  return(x)
}


## Combine method for "bamlss" objects.
c.bamlss <- function(...)
{
  objects <- list(...)
  x <- NULL
  for(i in 1L:length(objects))
    x <- c(x, objects[i])
  Call <- match.call()
  names(x) <- as.character(Call[-1L])
  class(x) <- c("bamlss", "cbamlss")

  return(x)
}


## Fast computation of quantiles.
quick_quantiles <- function(X, samples)
{
  rval <- .Call("quick_quantiles", X, samples, PACKAGE = "bamlss")
  rval <- as.data.frame(rval)
  names(rval) <- c("2.5%", "50%", "97.5%")
  rval
}

fitted_matrix <- function(X, samples)
{
  if(ncol(X) != ncol(samples))
    stop("dimensions of design matrix and samples do not match!")
  fit <- .Call("fitted_matrix", X, as.matrix(samples), PACKAGE = "bamlss")
}


## Function to compute statistics from samples of a model term.
compute_s.effect <- function(x, get.X, fit.fun, psamples,
  FUN = NULL, snames, data, grid = -1, rug = TRUE)
{
  nt <- length(x$term)
  if(nt > 2)
    return(NULL)

  if(x$by != "NA") grid <- NA
  if(!is.na(grid)) {
    if(grid < 0) {
      grid <- if(nt > 1) {
        NA
      } else 100
    }
  }

  tterms <- NULL
  for(l in nt:1) {
    tterm <- x$term[l]
    for(char in c("(", ")", "[", "]")) {
      tterm <- gsub(char, ".", tterm, fixed = TRUE)
    }
    if(inherits(data[[tterm]], "ts"))
      data[[tterm]] <- as.numeric(data[[tterm]])
    tterms <- c(tterms, tterm)
  }

  for(char in c("(", ")", "[", "]")) {
    colnames(data) <- gsub(char, ".", colnames(data), fixed = TRUE)
    x$by <- gsub(char, ".", x$by, fixed = TRUE)
  }

  ## Data for rug plotting.
  rugp <- if(nt < 2 & rug) data[[x$term]] else NULL

  ## New x values for which effect should
  ## be calculated, n = 100.
  if(!is.na(grid)) {
    if(!is.factor(data[[tterms[1]]]) & !any(grepl("mrf", class(x))) &
      !any(grepl("re.", class(x), fixed = TRUE)) & !any(grepl("random", class(x)))) {
      xsmall <- TRUE
      nd <- list()
      for(j in tterms) {
        xr <- range(data[[j]], na.rm = TRUE)
        nd[[j]] <- seq(xr[1], xr[2], length = grid)
      }
      nd <- expand.grid(nd)
      grid <- nrow(nd)
      if(x$by != "NA") { ## FIXME: check by variables!
        if(!is.factor(data[[x$by]])) {
          xr <- range(data[[x$by]], na.rm = TRUE)
          nd[[x$by]] <- seq(xr[1], xr[2], length = grid)
        } else nd[[x$by]] <- rep(data[[x$by]], length.out = grid)
      }
      nd <- as.data.frame(nd)
      if(nt == 2L) {
        pid <- chull(as.matrix(data[, tterms]))
        pol <- data[c(pid, pid[1]), ]
        pip <- point.in.polygon(nd[, 1], nd[, 2], pol[, 1], pol[, 2])
        nd[pip < 1, ] <- NA
        nd <- na.omit(nd)
      }
      data0 <- data
      data <- nd
    } else xsmall <- FALSE
  } else {
    data0 <- data[, c(tterms, if(x$by != "NA") x$by else NULL), drop = FALSE]
    if(nt < 2) {
      if(x$by != "NA") {
        if(!is.factor(data[[x$by]]))
          data <- unique(data0)
      } else data <- unique(data0)
    }
    xsmall <- if((nrow(data) != nrow(data0)) & (nt < 2)) TRUE else FALSE
  }
  if(is.null(x$special)) {
    X <- get.X(data)
  } else {
    if(x$special) {
      X <- get.X(data[, c(tterms, if(x$by != "NA") x$by else NULL), drop = FALSE])
    } else
      X <- get.X(data)
  }

  ## Compute samples of fitted values.
  if((inherits(x, "mgcv.smooth") | inherits(x, "deriv.smooth")) & (nrow(psamples) > 39L) & is.null(FUN)) {
    smf <- quick_quantiles(X, psamples)
  } else {
    if(is.null(FUN)) {
      FUN <- c95
    }
    if(nt < 2) {
      fsamples <- apply(psamples, 1, function(g) {
        f <- fit.fun(X, g, expand = FALSE, no.sparse.setup = (nrow(psamples) < 2))
        f
      })
      smf <- t(apply(fsamples, 1, FUN))
    } else {
      smf <- 0
      for(i in 1:nrow(psamples)) {
        smf <- smf + drop(fit.fun(X, psamples[i, ], expand = FALSE))
      }
      smf <- as.matrix(smf / nrow(psamples), ncol = 1)
      colnames(smf) <- "50%"
    }
  }

  cnames <- colnames(smf)
  smf <- as.data.frame(smf)
  for(l in 1:nt) {
    smf <- cbind(data[[tterms[l]]], smf)
  }
  names(smf) <- c(x$term, cnames)

  if(is.null(FUN)) {
    FUN <- c95
  }

  if(x$by != "NA") { ## FIXME: hard coded fix for plotting varying coefficient terms!
    if(!is.factor(data[[x$by]])) {
      X <- X / data[[x$by]]

      fsamples <- apply(psamples, 1, function(g) {
        fit.fun(X, g, expand = FALSE)
      })

      smf <- t(apply(fsamples, 1, FUN = FUN))

      cnames <- colnames(smf)
      smf <- as.data.frame(smf)
      for(l in 1:nt) {
        smf <- cbind(data[[tterms[l]]], smf)
      }
      names(smf) <- c(x$term, cnames)
    }
  }

  by.drop <- NULL
  if(x$by != "NA" & !is.null(x$by.level)) {
    by.drop <- (if(xsmall) data0[[x$by]] else data[[x$by]]) == x$by.level
    if(!xsmall)
      smf <- smf[by.drop, ]
  }

  ## Assign class and attributes.
  smf <- unique(smf)
  class(smf) <- c(class(x), "data.frame")
  x[!(names(x) %in% c("term", "label", "bs.dim", "dim"))] <- NULL
  attr(x, "qrc") <- NULL
  attr(smf, "specs") <- x
  class(attr(smf, "specs")) <- class(x)
  attr(smf, "x") <- if(xsmall & nt < 2) data0[, tterms] else data[, tterms]
  attr(smf, "by.drop") <- by.drop
  attr(smf, "rug") <- rugp

  return(smf)
}


## Function to add partial residuals based on weights() and score() function.
add.partial <- function(x, samples = FALSE, nsamps = 100)
{
  if(!inherits(x, "bamlss"))
    stop("x must be a 'bamlss' object!")

  nx <- names(x$terms)
  family <- x$family

  if(!is.null(family$hess) & !is.null(family$score)) {
    y <- model.response(model.frame(x))
    eta <- fitted.bamlss(x, samples = samples, nsamps = nsamps)
    if(is.null(x$model.frame))
      x$model.frame <- model.frame(x)
    for(j in seq_along(x$terms)) {
      if(!is.null(x$terms[[j]]$effects)) {
        peta <- family$map2par(eta)
        weights <- family$hess[[nx[j]]](y, peta, id = nx[j])
        score <- family$score[[nx[j]]](y, peta, id = nx[j])
        z <- eta[[nx[j]]] + 1 / weights * score
        ne <- names(x$terms[[j]]$effects)
        for(sj in seq_along(ne)) {
          f <- predict.bamlss(x, model = nx[j], term = ne[sj], nsamps = nsamps)
          term <- attr(x$terms[[j]]$effects[[ne[sj]]], "specs")$term
          e <- z - eta[[nx[j]]] + f
          if(is.null(attr(x$terms[[j]]$effects[[ne[sj]]], "specs")$xt$center)) {
            e <- e - mean(e)
          } else {
            if(attr(x$terms[[j]]$effects[[ne[sj]]], "specs")$xt$center)
              e <- e - mean(e)
          }
          e <- data.frame(x$model.frame[, term], e)
          names(e) <- c(term, "partial.resids")
          attr(x$terms[[j]]$effects[[ne[sj]]], "partial.resids") <- e
        }
      }
    }
  } else {
    stop("cannot compute partial residuals, no score() and hess() function in family object!")
  }
  x
}


## Helper function for prediction mean an 95% credible interval.
c95 <- function(x)
{
  qx <- quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  return(c(qx[1], "Mean" = mean(x, na.rm = TRUE), qx[2]))
}


## A prediction method for "bamlss" objects.
## Prediction can also be based on multiple chains.
predict.bamlss <- function(object, newdata, model = NULL, term = NULL, match.names = TRUE,
  intercept = TRUE, type = c("link", "parameter"), FUN = function(x) { mean(x, na.rm = TRUE) },
  trans = NULL, what = c("samples", "parameters"), nsamps = NULL, verbose = FALSE, drop = TRUE,
  cores = NULL, chunks = 1, ...)
{
  ## If data have been scaled (scale.d=TRUE)
  if ( ! missing(newdata) & ! is.null(attr(object$model.frame,'scale')) ) {
    sc <- attr(object$model.frame,'scale')
    for ( name in unique(unlist(lapply(sc,names))) ) {
      newdata[,name] <- (newdata[,name] - sc$center[name] ) / sc$scale[name]
    }
  }
  FUN2 <- function(x, ...) FUN(x)
  if(missing(newdata))
    newdata <- NULL
  family <- object$family
  if(!is.null(family$predict)) {
    if(is.function(family$predict)) {
      return(family$predict(object = object, newdata = newdata, model = model, term = term,
        intercept = intercept, type = type, FUN = FUN2, trans = trans, what = what, nsamps = nsamps,
        verbose = verbose, drop = drop, cores = cores, chunks = chunks, ...))
    }
  }
  if(is.null(object$x))
    object$x <- smooth.construct(object)
  if(is.null(newdata)) {
    newdata <- model.frame(object)
  } else {
    if(is.character(newdata)) {
      if(file.exists(newdata <- path.expand(newdata)))
        newdata <- read.table(newdata, header = TRUE, ...)
      else stop("cannot find newdata")
    }
    if(is.matrix(newdata) || is.list(newdata))
      newdata <- as.data.frame(newdata)
    ## FIXME: ??? newdata <- model.frame.bamlss.frame(object, data = newdata)
  }
  if(!is.null(attr(object, "fixed.names")))
    names(newdata) <- rmf(names(newdata))
  nn <- names(newdata)
  nn <- all.vars(as.formula(paste("~", paste(nn, collapse = "+")), env = NULL))
  rn <- response.name(object, keep.functions = TRUE)
  nn <- nn[!(nn %in% rn)]
  tl <- term.labels2(object, model = model, intercept = intercept, type = 2)
  nx <- names(tl)
  if(!is.null(term)) {
    enames <- vector(mode = "list", length = length(nx))
    for(j in term) {
      for(i in seq_along(tl)) {
        if(!is.character(j)) {
          if(j > 0 | j < length(tl[[i]]))
            enames[[i]] <- c(enames[[i]], tl[[i]][j])
        } else {
          if(grepl("intercept", tolower(j), fixed = TRUE)) {
            if("(Intercept)" %in% tl[[i]])
              enames[[i]] <- c(enames[[i]], "(Intercept)")
          } else {
            k <- if(match.names) grep(j, tl[[i]], fixed = TRUE) else which(tl[[i]] == j)
            if(length(k))
              enames[[i]] <- c(enames[[i]], tl[[i]][k])
          }
        }
      }
    }
    names(enames) <- nx
  } else enames <- tl
  if(intercept) {
    intcpt <- unlist(lapply(enames, function(x) { any(grepl("intercept", tolower(x))) }))
    if(any(!intcpt)) {
      for(i in seq_along(intcpt)) {
        if(!intcpt[i])
          enames[[i]] <- c(enames[[i]], "(Intercept)")
      }
    }
  }
  enames <- lapply(lapply(enames, unique), function(x) {
    x <- x[!is.na(x)]
    return(if(length(x) < 1) NULL else x)
  })
  if(all(is.null(unlist(enames))))
    stop("argument term is specified wrong!")
  ff <- as.formula(paste("~", paste(unique(unlist(enames)), collapse = "+")), env = NULL)
  vars <- all.vars.formula(ff)
  if(!all(vars[vars != "Intercept"] %in% nn))
    stop("cannot compute prediction, variables missing in newdata!")
  type <- match.arg(type)
  what <- match.arg(what)
  if(!is.null(object$samples) & what == "samples") {
    samps <- samples(object, model = model, ...)
    if(!is.null(nsamps)) {
      i <- seq(1, nrow(samps), length = nsamps)
      samps <- samps[i, , drop = FALSE]
    }
  } else {
    if(is.null(object$parameters))
      stop("cannot find any parameters!")
    samps <- parameters(object, model = model, list = FALSE, extract = TRUE, ...)
    cn <- names(samps)
    samps <- matrix(samps, nrow = 1)
    colnames(samps) <- cn
    samps <- as.mcmc(samps)
  }

  ## Remove samples not needed for predictions!
  cn <- colnames(samps)
  drop2 <- grep2(c(".tau2", ".alpha", ".edf", ".accepted", ".dic", ".loglik", ".logpost"),
    tolower(cn), fixed = TRUE)
  if(length(drop2))
    cn <- cn[-drop2]
  samps <- samps[, cn, drop = FALSE]

  env <- environment(object$formula)

  enames <- lapply(enames, function(x) {
    if(is.null(x)) return(NULL)
    f <- as.formula(paste("~", paste(x, collapse = "+")), env = NULL)
    all.labels.formula(f, full.names = TRUE)
  })

  if(!is.null(list(...)$get.bamlss.predict.setup)) {
    return(list("samps" = samps, "enames" = enames, "intercept" = intercept,
      "FUN" = FUN2, "trans" = trans, "type" = type, "nsamps" = nsamps, "env" = env))
  }

  pred_fun <- function(pred, id) {
    if(type != "link") {
      links <- family$links[nx]
      if(length(links) > 0) {
        if(links[id] != "identity") {
          linkinv <- make.link2(links[id])$linkinv
          pred <- linkinv(pred)
        }
      } else {
        warning(paste("could not compute predictions on the scale of parameter",
          ", predictions on the scale of the linear predictor are returned!", sep = ""))
      }
    }
    if(!is.null(trans)) {
      if(!is.list(trans)) {
        trans <- rep(list(trans), length = length(nx))
        names(trans) <- nx
      }
      if(!is.null(trans[[id]])) {
        if(!is.function(trans[[id]]))
          stop("argument trans must be a list of transformer functions!")
        pred <- trans[[id]](pred)
      }
    }
    pred <- apply(pred, 1, FUN2, ...)
    if(!is.null(dim(pred)))
      pred <- t(pred)
    return(pred)
  }

  ia <- interactive()

  if(is.null(cores)) {
    pred <- list()
    if(chunks > 1) {
      chunk_id <- sort(rep(1:chunks, length.out = nrow(newdata)))
      newdata <- split(newdata, chunk_id)
      chunks <- length(newdata)
      for(i in nx) {
        pred[[i]] <- NULL
        for(j in 1:chunks) {
          if(verbose) {
            cat(if(ia) "\r" else "\n")
            cat("predicting chunk", j, "of", chunks, "...")
            if(.Platform$OS.type != "unix" & ia) flush.console()
          }
          if(j < 2) {
            pred[[i]] <- pred_fun(.predict.bamlss(i, object$x[[i]], samps,
              enames[[i]], intercept, nsamps, newdata[[j]]), id = i)
          } else {
            if(is.null(dim(pred[[i]]))) {
              pred[[i]] <- c(pred[[i]], pred_fun(.predict.bamlss(i, object$x[[i]], samps,
                enames[[i]], intercept, nsamps, newdata[[j]]), id = i))
            } else {
              pred[[i]] <- rbind(pred[[i]], pred_fun(.predict.bamlss(i, object$x[[i]], samps,
                enames[[i]], intercept, nsamps, newdata[[j]]), id = i))
            }
          }
        }
        if(verbose) cat("\n")
      }
    } else {
      for(i in nx) {
        pred[[i]] <- pred_fun(.predict.bamlss(i, object$x[[i]], samps,
          enames[[i]], intercept, nsamps, newdata), id = i)
      }
    }
  } else {
    parallel_fun <- function(k) {
      pred <- list()
      if(chunks > 1) {
        chunk_id <- sort(rep(1:chunks, length.out = nrow(newdata[[k]])))
        nd <- split(newdata[[k]], chunk_id)
        chunks <- length(nd)
        for(i in nx) {
          for(j in 1:chunks) {
            if(verbose)
              cat("\npredicting chunk", j, "of", chunks, "...")

            if(j < 2) {
              pred[[i]] <- pred_fun(.predict.bamlss(i, object$x[[i]], samps,
                enames[[i]], intercept, nsamps, nd[[j]]), id = i)
            } else {
              if(is.null(dim(pred[[i]]))) {
                pred[[i]] <- c(pred[[i]], pred_fun(.predict.bamlss(i, object$x[[i]], samps,
                  enames[[i]], intercept, nsamps, nd[[j]]), id = i))
              } else {
                pred[[i]] <- rbind(pred[[i]], pred_fun(.predict.bamlss(i, object$x[[i]], samps,
                  enames[[i]], intercept, nsamps, nd[[j]]), id = i))
              }
            }
          }
        }
        if(verbose) cat("\n")
      } else {
        for(i in nx) {
          pred[[i]] <- pred_fun(.predict.bamlss(i, object$x[[i]], samps,
            enames[[i]], intercept, nsamps, newdata[[k]]), id = i)
        }
      }
      return(pred)
    }

    core_id <- sort(rep(1:cores, length.out = nrow(newdata)))
    newdata <- split(newdata, core_id)
    cores <- length(newdata)

    pred <- parallel::mclapply(1:cores, parallel_fun, mc.cores = cores)

    if(cores < 2) {
      pred <- pred[[1]]
    } else {
      for(i in nx) {
        for(j in 2:cores) {
          pred[[1]][[i]] <- if(is.matrix(pred[[1]][[i]])) {
            rbind(pred[[1]][[i]], pred[[j]][[i]])
          } else c(pred[[1]][[i]], pred[[j]][[i]])
          pred[[j]][[i]] <- NA
        }
      }
      pred <- pred[[1]]
    }
  }

  if((length(pred) < 2) & drop)
    pred <- pred[[1]]

  return(pred)
}


.predict.bamlss <- function(id, x, samps, enames, intercept, nsamps, data)
{
  if("smooth.construct" %in% names(x))
    x <- x$smooth.construct
  snames <- colnames(samps)
  enames <- gsub("p.Intercept", "p.(Intercept)", enames, fixed = TRUE)
  has_intercept <- any(grepl(paste(id, "p", "(Intercept)", sep = "."), snames, fixed = TRUE))
  if(!has_intercept) {
    if(any(grepl(".p.model.matrix.(Intercept)", snames, fixed = TRUE)))
      has_intercept <- TRUE
    if(any(grepl(".p.model.matrix.Intercept", snames, fixed = TRUE)))
      has_intercept <- TRUE
  }
  if(intercept & has_intercept)
    enames <- c("p.(Intercept)", enames)
  enames <- unique(enames)
  if(!has_intercept) {
    if(length(i <- grep("intercept", enames, ignore.case = TRUE)))
      enames <- enames[-i]  
  }
  ec <- sapply(enames, function(x) {
    paste(strsplit(x, "", fixed = TRUE)[[1]][1:2], collapse = "")
  })
  enames2 <- sapply(enames, function(x) {
    paste(strsplit(x, "", fixed = TRUE)[[1]][-c(1:2)], collapse = "")
  })
  eta <- matrix(0, nrow = nrow(data), ncol = nrow(samps))
  if(length(i <- grep("p.", ec))) {
    for(j in enames2[i]) {
      if(j != "(Intercept)") {
        f <- as.formula(paste("~", if(has_intercept) "1" else "-1", "+", j))
        X <- model.matrix(f, data = data)
        if(has_intercept)
          X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
        if(grepl(":", j)) {
          jt <- strsplit(j, ":", fixed = TRUE)[[1]]
          ok <- NULL
          for(jjt in jt)
            ok <- cbind(ok, grepl(jjt, snames, fixed = TRUE), grepl(paste(id, "p", sep = "."), snames, fixed = TRUE) | grepl(paste(id, "p.model.matrix", sep = "."), snames, fixed = TRUE))
          ok <- apply(ok, 1, all)
          sn <- snames[ok]
        } else {
          sn <- snames[grep2(paste(id, "p", j, sep = "."), snames, fixed = TRUE)]
          sn2 <- paste(sn, ".", sep = "")
          sn <- sn[grepl(paste(id, ".p.", j, ".", sep = ""), sn2, fixed = TRUE)]
          sn <- sn[!grepl(":", sn, fixed = TRUE)]
        }
        if(!length(sn)) {
          sn <- snames[grep2(paste(id, "p.model.matrix", j, sep = "."), snames, fixed = TRUE)]
        }
        if(ncol(X) > length(sn)) {
          sn2 <- gsub(paste(id, "p.", sep = "."), "", sn, fixed = TRUE)
          X <- X[, sn2, drop = FALSE]
        }
        eta <- eta + fitted_matrix(X, samps[, sn, drop = FALSE])
      } else {
        sn <- snames[grep2(paste(id, "p", j, sep = "."), snames, fixed = TRUE)]
        if(!length(sn))
          sn <- snames[grep2(paste(id, "p.model.matrix", j, sep = "."), snames, fixed = TRUE)]
        if(length(sn))
          eta <- eta + fitted_matrix(matrix(1, nrow = nrow(data), ncol = 1), samps[, sn, drop = FALSE])
      }
    }
  }
  if(length(i <- grep("s.", ec))) {
    for(j in enames2[i]) {
      for(jj in grep(j, names(x), fixed = TRUE, value = TRUE)) {
        sn <- snames[grep2(paste(id, "s", jj, sep = "."), snames, fixed = TRUE)]
        if(!inherits(x[[jj]], "no.mgcv") & !inherits(x[[jj]], "special")) {
          X <- PredictMat(x[[jj]], data)
          eta <- eta + fitted_matrix(X, samps[, sn, drop = FALSE])
        } else {
          if(is.null(x[[jj]]$PredictMat)) {
            X <- PredictMat(x[[jj]], data)
          } else {
            X <- x[[jj]]$PredictMat(x[[jj]], data)
          }
          fit <- apply(samps[, sn, drop = FALSE], 1, function(b) {
            x[[jj]]$fit.fun(X, b)
          })
          eta <- eta + fit
        }
      }
    }
  }
  eta
}


.fitted.bamlss <- function(id, x, samps)
{
  snames <- colnames(samps)
  if(length(i <- grep2(c(".alpha", ".edf", ".tau2", ".accepted"), snames, fixed = TRUE)))
    snames <- snames[-i]
  eta <- 0
  if(!is.null(x$model.matrix)) {
    if(ncol(x$model.matrix) > 0) {
      sn <- paste(id, "p", if(is.null(colnames(x$model.matrix))) {
        paste("b", 1:ncol(x$model.matrix))
      } else colnames(x$model.matrix), sep = ".")
      eta <- eta + fitted_matrix(x$model.matrix, samps[, sn, drop = FALSE])
    }
  }
  if(!is.null(x$smooth.construct)) {
    for(j in names(x$smooth.construct)) {
      sn <- grep(paste(id, if(j != "model.matrix") "s" else "p", j, sep = "."), snames,
        fixed = TRUE, value = TRUE)
      if(!length(sn)) {
        if(j == "model.matrix")
          sn <- grep(paste(id, ".p.", sep = ""), snames, fixed = TRUE, value = TRUE)
      }
      if(!length(sn))
        stop(paste('no fitted matrix for "', id, '", "', j, '"!', sep = ""))
      if(j != "model.matrix") {
        if(!inherits(x$smooth.construct[[j]], "no.mgcv") & !inherits(x$smooth.construct[[j]], "special")) {
          fit <- fitted_matrix(x$smooth.construct[[j]]$X, samps[, sn, drop = FALSE])
          if(!is.null(x$smooth.construct[[j]]$binning$match.index))
            fit <- fit[x$smooth.construct[[j]]$binning$match.index, , drop = FALSE]
          eta <- eta + fit
        } else {
          fit <- apply(samps[, sn, drop = FALSE], 1, function(b) {
            x$smooth.construct[[j]]$fit.fun(x$smooth.construct[[j]]$X, b)
          })
          eta <- eta + fit
        }
      } else {
        fit <- fitted_matrix(x$smooth.construct[[j]]$X, samps[, sn, drop = FALSE])
        if(!is.null(x$smooth.construct[[j]]$binning$match.index))
          fit <- fit[x$smooth.construct[[j]]$binning$match.index, , drop = FALSE]
        eta <- eta + fit
      }
    }
  }
  eta
}


## Setup function for handling "special" model terms.
s2 <- function(...)
{
  rval <- s(...)
  rval$special <- TRUE
  rval$label <- gsub("s(", "s2(", rval$label, fixed = TRUE)
  rval
}


## Setup function for random scaling terms.
rsc <- function(..., by = NA)
{
  by <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  if(by != "NA") {
    if(!grepl("~", by, fixed = TRUE)) {
      if(by == ".") 
        stop("by=. not allowed")
      by <- paste("~", by)
    }
  }
  rval <- s(...)
  rval$by.formula <- if(by != "NA") as.formula(by) else NULL
  rval$class <- class(rval)
  rval$special <- TRUE
  class(rval) <- "rsc.smooth.spec"
  rval
}


## Smooth constructor function for random scaling terms.
smooth.construct.rsc.smooth.spec <- function(object, data, knots) {
  class(object) <- object$class
  acons <- TRUE
  if(!is.null(object$xt$center))
    acons <- object$xt$center
  rval <- smoothCon(object, data, knots, absorb.cons = acons)
  rval <- rval[[1]]
  rval$class <- class(rval)
  if(!is.null(object$by.formula)) {
    ft <- terms(object$by.formula, keep.order = TRUE)
    vars <- attr(ft, "term.labels")
    if(length(vars)) {
      rs.by <- list()
      for(j in vars) {
        rs.by[[j]] <- data[[j]]
        if(!is.factor(rs.by[[j]])) stop("random scaling by variables must be factors!")
      }
      n <- length(vars)
      g <- paste("g[", 1:n, "]", sep = "", collapse = " + ")
      fun <- paste("function(X, g) { ", "(", if(attr(ft, "intercept")) "1 + ",
        g, ") * (X %*% ", if(n > 1) paste("g[-(1:", n, ")]", sep = "") else "g[-1]", ") }", sep = "")
      rval$fit.fun <- eval(parse(text = fun))
      rval$rs.by <- rs.by
      rval$by.vars <- vars
      rval$by.formula <- object$by.formula
      rval$one <- attr(ft, "intercept")
    }
  } else {
    rval$fit.fun <- function(X, g, ...) {
      X %*% as.numeric(g)
    }
  }

  class(rval) <- "rsc.smooth"
  rval
}


## Rational smooths constructor.
rs <- function(formula, link = "log", ...)
{
  formula <- deparse(substitute(formula), backtick = TRUE, width.cutoff = 500)
  formula <- gsub("[[:space:]]", "", formula)
  if(!grepl("~", strsplit(formula, "")[[1]][1]))
    formula <- paste("~", formula, sep = "")
  formula <- as.Formula(formula)
  nd <- TRUE
  if(length(formula)[2] < 2L) {
    formula <- as.Formula(formula(formula), formula(formula))
    nd <- FALSE
  }
  formula <- formula(formula, lhs = 0, drop = FALSE)
  fn <- formula(formula, rhs = 1)
  fd <- formula(formula, rhs = 2)

  kn <- sum(unlist(attr(terms(fn, specials = c("s", "te", "t2", "ti")), "specials")))
  kd <- sum(unlist(attr(terms(fn, specials = c("s", "te", "t2", "ti")), "specials")))

  fnl <- all.labels.formula(fn)
  fdl <- all.labels.formula(fd)
  fnl <- paste(fnl, collapse = "+")
  fdl <- paste(fdl, collapse = "+")
  if(fdl == "")
    nd <- FALSE
  label <- if(nd) paste("rs(", fnl, "|", fdl, ")", sep = "") else paste("rs(", fnl, ")", sep = "")
  vn <- all.vars.formula(fn)
  vd <- all.vars.formula(fd)
  formula <- bamlss.formula(list(fn, fd))
  names(formula) <- c("numerator", "denominator")
  xt <- list(...)
  if(is.null(xt$df))
    xt$df <- 5
  if(is.null(xt$update.nu))
    xt$update.nu <- FALSE
  if(is.null(xt$nu))
    xt$nu <- 0.1
  if(is.null(xt$do.optim))
    xt$do.optim <- TRUE

  rval <- list(
    "formula" = formula,
    "term" = unique(c(vn, vd)),
    "label" = label,
    "special" = TRUE,
    "link" = link,
    "xt" = xt,
    "by" = "NA",
    "nspecials" = kn + kd
  )
  rval$dim <- length(rval$term)

  class(rval) <- "rs.smooth.spec"
  rval
}

smooth.construct.rs.smooth.spec <- function(object, data, knots) 
{
  object$linkfun <- make.link2("identity")$linkfun
  object$linkinv <- make.link2("identity")$linkinv
  object$scale_linkfun <- make.link2(object$link)$linkfun
  object$scale_linkinv <- make.link2(object$link)$linkinv
  object$scale_mu.eta <- make.link2(object$link)$mu.eta
  object$scale_mu.eta2 <- make.link2(object$link)$mu.eta2
  object$dev.resids <- gaussian()$dev.resids
  object$aic <- gaussian()$aic
  object$mu.eta <- gaussian()$mu.eta
  object$variance <- gaussian()$variance
  object$dispersion <- function(wresiduals, wweights) {
    sum(wresiduals^2, na.rm = TRUE) / sum(wweights, na.rm = TRUE)
  }

  center <- if(!is.null(object$xt$center)) {
    object$xt$center
  } else TRUE

  object$X <- bamlss.engine.setup(design.construct(object$formula, data = data, knots = knots),
    df = rep(object$xt$df, object$nspecials))

  if(!is.null(object$X[[2]]$smooth.construct$model.matrix)) {
    cn <- colnames(object$X[[2]]$smooth.construct$model.matrix$X)
    if("(Intercept)" %in% cn)
      object$X[[2]]$smooth.construct$model.matrix$X <- object$X[[2]]$smooth.construct$model.matrix$X[, cn != "(Intercept)", drop = FALSE]
    if(ncol(object$X[[2]]$smooth.construct$model.matrix$X) < 1) {
      object$X[[2]]$smooth.construct$model.matrix <- NULL
      object$X[[2]]$terms <- drop.terms.bamlss(object$X[[2]]$terms, pterms = FALSE, keep.intercept = FALSE)
    } else {
      object$X[[2]]$smooth.construct$model.matrix$term <- gsub("(Intercept)+", "",
        object$X[[2]]$smooth.construct$model.matrix$term, fixed = TRUE)
      object$X[[2]]$smooth.construct$model.matrix$state$parameters <- object$X[[2]]$smooth.construct$model.matrix$state$parameters[-1]
      attr(object$X[[2]]$terms, "intercept") <- 0
    }
  }

  parameters <- NULL
  npar <- edf <- 0
  object$xmat <- object$zmat <- NULL
  for(j in 1:2) {
    for(sj in seq_along(object$X[[j]]$smooth.construct)) {
      pn <- if(j < 2) "n" else "d"
      names(object$X[[j]]$smooth.construct[[sj]]$state$parameters) <- paste(paste(pn, sj, sep = ""),
        names(object$X[[j]]$smooth.construct[[sj]]$state$parameters), sep = ".")
      tpar <- object$X[[j]]$smooth.construct[[sj]]$state$parameters
      parameters <- c(parameters, object$X[[j]]$smooth.construct[[sj]]$state$parameters)
      npar <- npar + ncol(object$X[[j]]$smooth.construct[[sj]]$X)
      edf <- edf + object$X[[j]]$smooth.construct[[sj]]$state$edf
      if(j < 2)
        object$xmat <- cbind(object$xmat, object$X[[j]]$smooth.construct[[sj]]$X)
      else
        object$zmat <- cbind(object$zmat, object$X[[j]]$smooth.construct[[sj]]$X)
    }
  }

  object$prior <- function(parameters) {
    lp <- 0
    for(j in 1:2) {
      for(sj in seq_along(object$X[[j]]$smooth.construct)) {
        id <- paste(if(j < 2) "n" else "d", sj, sep = "")
        tpar <- parameters[grep(id, names(parameters), fixed = TRUE)]
        lp <- lp + object$X[[j]]$smooth.construct[[sj]]$prior(tpar)
      }
    }
    return(lp)
  }

  object$grad <- function(parameters) {
    lg <- NULL
    for(j in 1:2) {
      for(sj in seq_along(object$X[[j]]$smooth.construct)) {
        id <- paste(if(j < 2) "n" else "d", sj, sep = "")
        tpar <- parameters[grep(id, names(parameters), fixed = TRUE)]
        lg <- c(lg, object$X[[j]]$smooth.construct[[sj]]$grad(score = NULL, tpar, full = FALSE))
      }
    }
    return(lg)
  }

  object$hess <- function(parameters) {
    lh <- list(); k <- 1
    for(j in 1:2) {
      for(sj in seq_along(object$X[[j]]$smooth.construct)) {
        id <- paste(if(j < 2) "n" else "d", sj, sep = "")
        tpar <- parameters[grep(id, names(parameters), fixed = TRUE)]
        lh[[k]] <- object$X[[j]]$smooth.construct[[sj]]$hess(score = NULL, tpar, full = FALSE)
        k <- k + 1
      }
    }
    lh <- as.matrix(do.call("bdiag", lh))
    return(lh)
  }

  rs_intcpt <- 1 - .Machine$double.eps

  object$fit.fun <- function(X, b, ..., nocenter = FALSE, mu = FALSE, scale = FALSE) {
    if(!is.null(names(b)))
      b <- get.par(b, "b")
    fn <- fd <- 0
    k1 <- 1
    for(sj in seq_along(X[[1]]$smooth.construct)) {
      k2 <- k1 + ncol(X[[1]]$smooth.construct[[sj]]$X) - 1
      fn <- fn + X[[1]]$smooth.construct[[sj]]$X %*% b[k1:k2]
      k1 <- k2 + 1
    }
    if(mu) return(drop(fn))
    for(sj in seq_along(X[[2]]$smooth.construct)) {
      k2 <- k1 + ncol(X[[2]]$smooth.construct[[sj]]$X) - 1
      fd <- fd + X[[2]]$smooth.construct[[sj]]$X %*% b[k1:k2]
      k1 <- k2 + 1
    }
    if(scale) return(drop(fd))
    fd <- object$scale_linkinv(drop(object$scale_linkfun(rs_intcpt) + fd))
    f <- fn / fd
    if(center & !nocenter)
      f <- f - mean(f)
    return(drop(f))
  }

  object$update <- function(x, family, y, eta, id, weights, criterion, ...)
  {
    peta <- family$map2par(eta)
    score <- drop(family$score[[id]](y, peta))
    hess <- drop(family$hess[[id]](y, peta))

    gradfun <- function(b) {
      eta_mu <- x$fit.fun(x$X, b, mu = TRUE)
      eta_scale <- x$scale_linkfun(rs_intcpt) + x$fit.fun(x$X, b, scale = TRUE)
      scale <- x$scale_linkinv(eta_scale)

      w_mu <- 1 / scale
      w_scale <- eta_mu * (-1 / (scale^2)) * x$scale_mu.eta(eta_scale)

      grad <- c(colSums(x$xmat * score * w_mu), if(!is.null(x$zmat)) colSums(x$zmat * score * w_scale) else NULL)

      grad
    }

    hessfun <- function(b) {
      eta_mu <- x$fit.fun(x$X, b, mu = TRUE)
      eta_scale <- x$scale_linkfun(rs_intcpt) + x$fit.fun(x$X, b, scale = TRUE)
      scale <- x$scale_linkinv(eta_scale)

      w_mu <- 1 / scale
      w_scale <- eta_mu * (-1 / (scale^2)) * x$scale_mu.eta(eta_scale)

      Hd <- crossprod(x$xmat * (w_mu^2 * hess), x$xmat)
      Hn <- if(!is.null(x$zmat)) crossprod(x$zmat * (w_scale^2 * hess), x$zmat) else NULL

      return(if(!is.null(Hn)) as.matrix(do.call("bdiag", list(Hd, Hn))) else as.matrix(Hd))
    }

    b0 <- get.state(x, "b")
    par0 <- x$state$parameters

    hess <- hessfun(b0)
    grad <- gradfun(b0)

    eta[[id]] <- eta[[id]] - fitted(x$state)

    if(length(start <- get.par(x$state$parameters, "tau2")) & x$xt$do.optim) {
      env <- new.env()
      args <- list(...)
      edf0 <- args$edf - x$state$edf
      k <- ncol(x$xmat)

      objfun1 <- function(tau2) {
        par1 <- set.par(par0, tau2, "tau2")
        grad <- -1 * (grad + x$grad(par1))
        Sigma <- matrix_inv(hess + x$hess(par1))
        Hs <- Sigma %*% grad
        if(x$xt$update.nu) {
          objfun_nu1 <- function(nu) {
            b1 <- drop(b0 - nu * Hs)
            par2 <- set.par(par1, b1, "b")
            eta[[id]] <- eta[[id]] + x$fit.fun(x$X, par2)
            logLik <- family$loglik(y, family$map2par(eta))
            logPost <- logLik + x$prior(par2)
            return(-1 * logPost)
          }
          nu <- optimize(f = objfun_nu1, interval = c(0, 1))$minimum
        } else {
          nu <- x$xt$nu
        }
        b1 <- drop(b0 - nu * Hs)
        par2 <- set.par(par1, b1, "b")
        fit <- x$fit.fun(x$X, par2)
        eta[[id]] <- eta[[id]] + fit
        logLik <- family$loglik(y, family$map2par(eta))
        edf1 <- sum_diag(hess[1:k, 1:k] %*% Sigma[1:k, 1:k])
        edf2 <- sum_diag(hess[-(1:k), -(1:k)] %*% Sigma[-(1:k), -(1:k)])
        edf3 <- edf1 + edf2 - 1
        edf <- edf0 + edf3
        ic <- get.ic2(logLik, edf, length(eta[[id]]), criterion)
        if(!is.null(env$ic_val)) {
          if((ic < env$ic_val) & (ic < env$ic00_val)) {
            opt_state <- list("parameters" = par2,
              "fitted.values" = fit, "edf" = edf3, "nu" = nu)
            assign("state", opt_state, envir = env)
            assign("ic_val", ic, envir = env)
          }
        } else assign("ic_val", ic, envir = env)
        return(ic)
      }

      assign("ic00_val", objfun1(get.state(x, "tau2")), envir = env)

      tau2 <- tau2.optim(objfun1, start = start)

      if(!is.null(env$state))
        return(env$state)

      par0 <- set.par(par0, tau2, "tau2")
    }

    hess <- hess + x$hess(par0)
    Sigma <- matrix_inv(hess)
    grad <- -1 * (grad + x$grad(par0))

    if(x$xt$update.nu) {
      objfun_nu2 <- function(nu) {
        b1 <- b0 - nu * Sigma %*% grad
        par0 <- set.par(par0, b1, "b")
        eta[[id]] <- eta[[id]] + x$fit.fun(x$X, par0)
        logLik <- family$loglik(y, family$map2par(eta))
        logPost <- logLik + x$prior(par0)
        return(-1 * logPost)
      }
      nu <- optimize(f = objfun_nu2, interval = c(0, 1))$minimum
    } else {
      nu <- x$xt$nu
    }
    b <- b0 - nu * grad %*% Sigma

    x$state$parameters <- set.par(x$state$parameters, b, "b")
    x$state$fitted.values <- x$fit.fun(x$X, x$state$parameters)

    return(x$state)
  }

  object$update99 <- function(x, family, y, eta, id, weights, criterion, ...)
  {
    args <- list(...)

    for(i in 1:2) {
      for(j in seq_along(x$X[[i]]$smooth.construct)) {
        peta <- family$map2par(eta)

        hess <- family$hess[[id]](y, peta, id = id, ...)

        if(!is.null(weights))
          hess <- hess * weights

        score <- family$score[[id]](y, peta, id = id, ...)

        eta_scale <- x$scale_linkfun(rs_intcpt) + x$fit.fun(x$X, x$state$parameters, scale = TRUE)
        scale <- x$scale_linkinv(eta_scale)
        if(i < 2) {
          hess <- hess * (1 / scale)^2
          score <- score * (1 / scale)
        } else {
          eta_mu <- x$fit.fun(x$X, x$state$parameters, mu = TRUE)
          hess <- hess * (eta_mu * (-1 / (scale^2)) * x$scale_mu.eta(eta_scale))^2
          score <- score * (eta_mu * (-1 / (scale^2)) * x$scale_mu.eta(eta_scale))
        }

        ## Compute working observations.
        z <- eta[[id]] + 1 / hess * score

        eta[[id]] <- eta[[id]] - fitted(x$state)
        e <- z - eta[[id]]

        XWX <- crossprod(x$X[[i]]$smooth.construct[[j]]$X * hess, x$X[[i]]$smooth.construct[[j]]$X)

        idj <- paste(if(i < 2) "n" else "d", j, sep = "")

        if(x$X[[i]]$smooth.construct[[j]]$fixed) {
          P <- matrix_inv(XWX, index = x$X[[i]]$smooth.construct[[j]]$sparse.setup)
        } else {
          S <- 0
          ij <- grep(paste(idj, ".tau2", sep = ""), names(x$state$parameters), fixed = TRUE)
          tau2 <- x$state$parameters[ij]
          for(jj in seq_along(x$X[[i]]$smooth.construct[[j]]$S))
            S <- S + 1 / tau2[jj] * x$X[[i]]$smooth.construct[[j]]$S[[jj]]
          P <- matrix_inv(XWX + S, index = x$X[[i]]$smooth.construct[[j]]$sparse.setup)
        }
        b <- drop(P %*% crossprod(x$X[[i]]$smooth.construct[[j]]$X * hess, e))
        ij <- grep(paste(idj, ".b", sep = ""), names(x$state$parameters), fixed = TRUE)
        x$state$parameters[ij] <- b
        x$state$fitted.values <- x$fit.fun(x$X, x$state$parameters)
        eta[[id]] <- eta[[id]] + fitted(x$state)
      }
    }

    return(x$state)
  }

  object$propose <- function(family, theta, id, eta, y, data, weights = NULL, ...)
  {
    theta <- theta[[id[1]]][[id[2]]]

    if(is.null(attr(theta, "fitted.values")))
      attr(theta, "fitted.values") <- data$fit.fun(data$X, theta)

    gradfun <- function(b, score) {
      eta_mu <- data$fit.fun(data$X, b, mu = TRUE)
      eta_scale <- data$scale_linkfun(rs_intcpt) + data$fit.fun(data$X, b, scale = TRUE)
      scale <- data$scale_linkinv(eta_scale)

      w_mu <- 1 / scale
      w_scale <- eta_mu * (-1 / (scale^2)) * data$scale_mu.eta(eta_scale)

      grad <- c(colSums(data$xmat * score * w_mu), if(!is.null(data$zmat)) colSums(data$zmat * score * w_scale) else NULL)

      grad
    }

    hessfun <- function(b, score, hess) {
      eta_mu <- data$fit.fun(data$X, b, mu = TRUE)
      eta_scale <- data$scale_linkfun(rs_intcpt) + data$fit.fun(data$X, b, scale = TRUE)
      scale <- data$scale_linkinv(eta_scale)

      w_mu <- 1 / scale
      w_scale <- eta_mu * (-1 / (scale^2)) * data$scale_mu.eta(eta_scale)

      Hd <- crossprod(data$xmat * (w_mu^2 * hess), data$xmat)
      Hn <- if(!is.null(data$zmat)) crossprod(data$zmat * (w_scale^2 * hess), data$zmat) else NULL

      return(if(!is.null(Hn)) as.matrix(do.call("bdiag", list(Hd, Hn))) else as.matrix(Hd))
    }

    peta <- family$map2par(eta)

    score <- process.derivs(family$score[[id[1]]](y, peta, id = id[1]))
    hess <- process.derivs(family$hess[[id[1]]](y, peta, id = id[1]))

    pibeta <- family$loglik(y, peta)
    p1 <- data$prior(theta)

    hess0 <- hessfun(theta, score, hess)
    Sigma <- matrix_inv(hess0 + data$hess(theta))
    xgrad <- -1 * (gradfun(theta, score) + data$grad(theta))

    if(all(is.na(Sigma)) | all(is.na(xgrad)))
      return(list("parameters" = theta, "alpha" = -Inf, "extra" = c("edf" = NA)))

    edf <- sum_diag(hess0 %*% Sigma) - 1

    ## Old position.
    g0 <- get.par(theta, "b")

    ## Get new position.
    mu <- drop(g0 - Sigma %*% xgrad)

    ## Sample new parameters.
    g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma))
    names(g) <- names(g0)
    theta2 <- set.par(theta, g, "b")

    ## Compute log priors.
    p2 <- data$prior(theta2)
    qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)

    ## Map predictor to parameter scale.
    fit <- data$fit.fun(data$X, theta2)
    eta[[id[1]]] <- eta[[id[1]]] - attr(theta, "fitted.values") + fit

    peta <- family$map2par(eta)
    score2 <- process.derivs(family$score[[id[1]]](y, peta, id = id[1]))
    hess2 <- process.derivs(family$hess[[id[1]]](y, peta, id = id[1]))

    ## Compute new log likelihood.
    pibetaprop <- family$loglik(y, peta)
    Sigma2 <- matrix_inv(hessfun(theta2, score2, hess2) + data$hess(theta2))
    xgrad2 <- -1 * (gradfun(theta2, score2) + data$grad(theta2))

    if(all(is.na(Sigma2)) | all(is.na(xgrad2)))
      return(list("parameters" = theta, "alpha" = -Inf, "extra" = c("edf" = NA)))

    mu2 <- drop(g - Sigma2 %*% xgrad2)
    qbeta <- dmvnorm(g0, mean = mu2, sigma = Sigma2, log = TRUE)

    ## Sample variance parameter.
    i <- grep("tau2", names(theta2))
    if(length(i)) {
      for(j in i) {
        theta2 <- uni.slice(theta2, data, NULL, NULL,
          NULL, id[1], j, logPost = gmcmc_logPost, lower = 0, ll = pibetaprop)
      }
    }

    ## Compute acceptance probablity.
    alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))

    ## New theta.
    attr(theta2, "fitted.values") <- fit

    return(list("parameters" = theta2, "alpha" = alpha, "extra" = c("edf" = edf)))
  }

  object$state <- list(
    "parameters" = parameters,
    "fitted.values" = rep(0, nrow(object$X[[1]]$smooth.construct[[1]]$X)),
    "edf" = edf
  )
  object$PredictMat <- function(object, data) {
    Predict.matrix.rs.smooth(object, data)
  }
  object$special.npar <- length(get.par(object$state$parameters, "b"))
  object$special.mpar <- function(...) { object$state$parameters }
  object$fixed <- FALSE
  object$fxsp <- FALSE
  object$S <- list(matrix(0, 1, 1))

  class(object) <- c("rs.smooth", "no.mgcv", "special")
  object
}

Predict.matrix.rs.smooth <- function(object, data)
{
  data <- as.data.frame(data)
  Xl <- list()
  for(j in 1:2) {
    Xl[[j]] <- list("smooth.construct" = list())
    for(sj in names(object$X[[j]]$smooth.construct)) {
      if(sj == "model.matrix") {
        f <- drop.terms.bamlss(object$formula[[j]], keep.response = FALSE, sterms = FALSE)
        Xl[[j]]$smooth.construct[[sj]] <- list("X" = model.matrix(f, data = data))
      } else {
        Xl[[j]]$smooth.construct[[sj]] <- list("X" = PredictMat(object$X[[j]]$smooth.construct[[sj]], data))
      }
    }
  }
  Xl
}

rs.plot <- function(x, model = NULL, term = NULL,
  what = c("numerator", "denominator"), type = "link", ...)
{
  tl <- term.labels2(x, model = model, pterms = FALSE, intercept = FALSE, type = 2, list = FALSE)
  tl <- grep("rs(", tl, fixed = TRUE, value = TRUE)
  if(length(tl) < 1)
    return(invisible(NULL))
  term <- if(!is.null(term)) {
    grep(term, tl, fixed = TRUE, value = TRUE)
  } else tl

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  if(!is.null(x$samples)) {
    samps <- samples(x, model = model)
  } else {
    if(is.null(x$parameters))
      stop("cannot find any parameters!")
    samps <- parameters(x, model = model, list = FALSE, extract = TRUE)
    cn <- names(samps)
    samps <- matrix(samps, nrow = 1)
    colnames(samps) <- cn
    samps <- as.mcmc(samps)
  }

  pl <- list()
  k <- 1
  for(j in seq_along(term)) {
    pn <- names(term)[j]
    tlj <- names(x$x[[pn]]$smooth.construct)
    i <- grep(term[j], tlj, fixed = TRUE)
    for(w in what) {
      nw <- nw0 <- names(x$x[[pn]]$smooth.construct[[i]]$X[[w]]$smooth.construct)
      nw <- nw[nw != "model.matrix"]
      if(length(nw) > 0) {
        for(ii in nw) {
          get.X <- function(data) {
            PredictMat(x$x[[pn]]$smooth.construct[[i]]$X[[w]]$smooth.construct[[ii]], data)
          }
          iii <- which(nw0 %in% ii)
          iii <- paste(if(w == "denominator") "d" else "n", iii, sep = "")
          iii <- paste(pn, ".s.", tlj[i], ".", iii, sep = "")
          sn <- colnames(samps)
          sn <- grep(iii, sn, fixed = TRUE, value = TRUE)
          if((w == "denominator") & (type != "link")) {
            FUN <- function(z) {
              c95(x$x[[pn]]$smooth.construct[[i]]$scale_linkinv(z))
            }
          } else {
            FUN <- NULL
          }

          pl[[k]] <- compute_s.effect(x$x[[pn]]$smooth.construct[[i]]$X[[w]]$smooth.construct[[ii]], get.X,
            x$x[[pn]]$smooth.construct[[i]]$X[[w]]$smooth.construct[[ii]]$fit.fun, samps[, sn, drop = FALSE],
            FUN = FUN, sn, model.frame(x), grid = -1, rug = TRUE)
          attr(pl[[k]], "specs")$label <- paste(attr(pl[[k]], "specs")$label, ".", w, sep = "")
          k <- k + 1
        }
      }
    }
  }

  par(mfrow = n2mfrow(length(pl)))
  for(j in seq_along(pl)) {
    plot.bamlss.effect(pl[[j]], ...)
  }
  invisible(pl)
}


## Lasso smooth constructor.
la <- function(formula, type = c("single", "multiple"), ...)
{
  formula <- deparse(substitute(formula), backtick = TRUE, width.cutoff = 500)
  formula <- gsub("[[:space:]]", "", formula)
  label <- NULL
  if(!grepl("+", formula, fixed = TRUE) & !grepl("-", formula, fixed = TRUE)) {
    if(formula %in% ls(envir = .GlobalEnv)) {
      label <- paste("la(", formula, ")", sep = "")
      f0 <- formula
      formula <- get(formula, envir = .GlobalEnv)
      if(!inherits(formula, "formula"))
        formula <- f0
    }
  }
  if(is.character(formula)) {
    if(!grepl("~", strsplit(formula, "")[[1]][1]))
      formula <- paste("~", formula, sep = "")
    formula <- as.formula(formula)
  }
  formula <- as.formula(formula)
  if(!any(grepl("+", formula, fixed = TRUE)) & !any(grepl("-", formula, fixed = TRUE)) & is.null(label))
    label <- paste("la(", paste(all.vars.formula(as.formula(formula)), collapse = "+"), ")", sep = "")
  vars <- unique(all.vars.formula(formula))
  rval <- list(
    "formula" = formula,
    "term" = vars,
    "label" = if(is.null(label)) paste("la(~", paste(vars, collapse = "+"), ")", sep = "") else label,
    "type" = match.arg(type),
    "by" = "NA"
  )
  rval$dim <- length(rval$term)
  rval$special <- TRUE
  rval$xt <- list(...)
  class(rval) <- "la.smooth.spec"
  rval
}


blockstand <- function(x, n)
{
  cn <- colnames(x)
  decomp <- qr(x)
  if(decomp$rank < ncol(x))
    stop("block standardization cannot be computed, matrix is not of full rank!")
  scale <- qr.R(decomp) * 1 / sqrt(n)
  x <- qr.Q(decomp) * sqrt(n)
  attr(x, "blockscale") <- scale
  colnames(x) <- cn
  x
}


smooth.construct.la.smooth.spec <- function(object, data, knots, ...)
{
  ridge <- if(is.null(object$xt[["ridge"]])) FALSE else object$xt[["ridge"]]
  fuse <- if(is.null(object$xt[["fuse"]])) FALSE else object$xt[["fuse"]]
  standardize <- if(is.null(object$xt[["standardize"]])) FALSE else object$xt[["standardize"]]
  fuse_type <- "nominal"
  if(is.logical(fuse)) {
    if(fuse)
      fuse <- "nominal"
  }
  if(!is.logical(fuse)) {
    if(is.character(fuse)) {
      fuse_type <- match.arg(fuse, c("nominal", "ordered"))
    } else {
      fuse_type <- switch(as.integer(fuse),
        "1" = "nominal",
        "2" = "ordered"
      )
    }
    fuse <- TRUE
  }
  object$fuse <- fuse
  object$fuse_type <- fuse_type
  object$standardize <- standardize
  
  data <- as.data.frame(data)
  nobs <- nrow(data)
  tl <- term.labels2(terms(object$formula), intercept = FALSE, list = FALSE)
  if(any(grepl("la(", tl, fixed = TRUE)))
    tl <- object$term
  object$X <- df <- group <- list()
  object$lasso <- list("trans" = list())
  k <- 1
  for(j in tl) {
    object$X[[j]] <- as.matrix(model.matrix(as.formula(paste("~", j)), data = data))
    if(length(i <- grep("Intercept", colnames(object$X[[j]]))))
      object$X[[j]] <- object$X[[j]][, -i, drop = FALSE]
    is_f <- is.factor(data[[j]])
    if(is_f) {
      group[[j]] <- k:(k + ncol(object$X[[j]]) - 1)
    } else {
      group[[j]] <- NA
    }
    k <- k + ncol(object$X[[j]])
    if(grepl(":", j, fixed = TRUE)) {
      j2 <- strsplit(j, ":")[[1]]
      is_f <- any(sapply(j2, function(i) is.factor(data[[i]])))
    }
    if(grepl("*", j, fixed = TRUE)) {
      j2 <- strsplit(j, "*")[[1]]
      is_f <- any(sapply(j2, function(i) is.factor(data[[i]])))
    }
    if(!fuse | standardize) {
      if(!is_f) {
        object$X[[j]] <- scale(object$X[[j]])
        object$lasso$trans[[j]] <- list(
          "center" = attr(object$X[[j]], "scaled:center"),
          "scale" = attr(object$X[[j]], "scaled:scale")
        )
      } else {
        object$X[[j]] <- blockstand(object$X[[j]], n = nobs)
        object$lasso$trans[[j]] <- list("blockscale" = attr(object$X[[j]], "blockscale"))
      }
    }
    if(is.null(colnames(object$X[[j]])))
      colnames(object$X[[j]]) <- paste("X", 1:ncol(object$X[[j]]), sep = "")
    if(!fuse | standardize) {
      df[[j]] <- sqrt(rep(ncol(object$X[[j]]), ncol(object$X[[j]])))
    } else {
      object$lasso$trans[[j]] <- list(
        "center" = 0.0,
        "scale" = 1.0
      )
      if(is.factor(data[[j]])) {
        df[[j]] <- colSums(object$X[[j]])
      } else {
        df[[j]] <- rep(nobs, ncol(object$X[[j]]))
        names(df[[j]]) <- colnames(object$X[[j]])
      }
    }
    object$lasso$trans[[j]]$colnames <- colnames(object$X[[j]])
  }
  df <- unlist(df)
  object$lasso$df <- df
  object$X <- do.call("cbind", object$X)
  object$S <- list()
  const <- object$xt$const
  if(is.null(const))
    const <- 1e-05
  if(!fuse) {
    if(object$type == "single") {
      object$S[[1]] <- function(parameters) {
        b <- get.par(parameters, "b")
        A <- df / sqrt(b^2 + const)
        for(j in seq_along(group)) {
          if(all(!is.na(group[[j]]))) {
            A[group[[j]]] <- df[group[[j]]] / rep(sqrt(sum(b[group[[j]]]^2) + const), length(group[[j]]))
          }
        }
        ## FIXME: adaptive weights: A <- A * MLpen ## 1 / abs(beta)
        A <- if(length(A) < 2) matrix(A, 1, 1) else diag(A)
        A
      }
      attr(object$S[[1]], "npar") <- ncol(object$X)
    } else {
      A <- list()
      for(j in 1:ncol(object$X)) {
        f <- c('function(parameters) {',
          '  b <- get.par(parameters, "b")',
          '  A <- diag(0, length(b))',
          paste('  A[', j, ',', j, '] <- df[', j, '] / sqrt(b[', j, ']^2 + const)', sep = ''),
          '  A',
          '}')
        A[[j]] <- eval(parse(text = paste(f, collapse = "\n")))
        attr(A[[j]], "npar") <- ncol(object$X)
      }
      object$S <- A
    }
  }
  if(fuse) {
    k <- ncol(object$X)
    if(fuse_type == "nominal") {
      Af <- matrix(0, ncol = choose(k, 2), nrow = k)
      combis <- combn(k, 2)
      for(ff in 1:ncol(combis)){
        Af[combis[1, ff], ff] <- 1
        Af[combis[2, ff], ff] <- -1
      }
      Af <- cbind(diag(k), Af)
    } else {
      Af <- diff(diag(k + 1))
      Af[1, 1] <- 1
      Af <- Af[, -ncol(Af), drop = FALSE]
    }
    beta <- object$xt$beta
    w <- rep(0, length = ncol(Af))
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
      if(!is.null(beta))
        w[ff] <- w[ff] * 1 / abs(t(Af[, ff]) %*% beta)
    }
    object$Af <- Af
    object$S[[ls <- length(object$S) + 1]] <- function(parameters, fixed.hyper = NULL) {
      b <- get.par(parameters, "b")
      if(!is.null(fixed.hyper)) {
        w <- fixed.hyper
      } else {
        if(length(i <- grep("lasso", names(parameters))))
          w <- parameters[i]
      }
      S <- 0
      for(k in 1:ncol(Af)) {
        tAf <- t(Af[, k])
        d <- drop(tAf %*% b)
        S <- S + w[k] / sqrt(d^2 + const) * Af[, k] %*% tAf
      }
      S
    }
    attr(object$S[[ls]], "npar") <- ncol(object$X)
  }

  object$xt[["prior"]] <- "ig"
  object$xt[["a"]] <- 10
  object$xt[["b"]] <- 1e-5
  object$fixed <- if(is.null(object$xt[["fx"]])) FALSE else object$xt[["fx"]]
  priors <- make.prior(object)
  object$prior <- priors$prior
  object$grad <- priors$grad
  object$hess <- priors$hess
  if(is.null(object$xt$lambda)) {
    object$xt$lambda <- if(is.null(object$xt[["sp"]])) 0.0001 else object$xt[["sp"]]
  } else {
    if(!is.null(object$xt[["sp"]]))
      object$xt$lambda <- object$xt[["sp"]]
  }
  object$xt$do.optim <- TRUE
  object$lasso$const <- const

  if(ridge)
    object$S <- list(diag(.Machine$double.eps, ncol(object$X)))

  object$xt[["binning"]] <- TRUE
  if(is.null(object$xt[["df"]]))
    object$xt[["df"]] <- if(!ridge) ceiling(ncol(object$X) * 0.9) else ceiling(ncol(object$X) * 0.3)
  object$ctype <- switch(object$type,
    "single" = 0,
    "multiple" = 1
  )
  object$C <- matrix(nrow = 0, ncol = ncol(object$X))
  object$side.constrain <- FALSE

  if(object$fixed)
    object$S <- NULL

  class(object) <- "lasso.smooth"

  object
}

Predict.matrix.lasso.smooth <- function(object, data)
{
  data <- as.data.frame(data)
  tl <- term.labels2(terms(object$formula), intercept = FALSE, list = FALSE)
  X <- list()
  for(j in tl) {
    X[[j]] <- as.matrix(model.matrix(as.formula(paste("~", j)), data = data))
    if(length(i <- grep("Intercept", colnames(X[[j]]))))
      X[[j]] <- X[[j]][, -i, drop = FALSE]
    is_f <- is.factor(data[[j]])
    if(grepl(":", j, fixed = TRUE)) {
      j2 <- strsplit(j, ":")[[1]]
      is_f <- any(sapply(j2, function(i) is.factor(data[[i]])))
    }
    if(grepl("*", j, fixed = TRUE)) {
      j2 <- strsplit(j, "*")[[1]]
      is_f <- any(sapply(j2, function(i) is.factor(data[[i]])))
    }
    if(is_f & is.null(object$lasso$trans[[j]]$blockscale))
      is_f <- FALSE
    if(!is_f) {
      X[[j]] <- (X[[j]] - object$lasso$trans[[j]]$center) / object$lasso$trans[[j]]$scale
    } else {
      X[[j]] <- X[[j]] %*% object$lasso$trans[[j]]$blockscale
    }
  }
  return(do.call("cbind", X))
}


## Neural networks.
n <- function(..., k = 10)
{
  ret <- la(..., k = k)
  ret$label <- gsub("la(", "n(", ret$label, fixed = TRUE)
  if(!is.null(node <- ret$xt$node)) {
    lab <- strsplit(ret$label, "")[[1]]
    lab <- paste(c(lab[-length(lab)], paste(',node="', node, '")', sep = '')), collapse = "", sep = "")
    ret$label <- lab
  }
  class(ret) <- "nnet.smooth.spec"
  ret
}

smooth.construct.nnet.smooth.spec <- function(object, data, knots, ...)
{
  split <- if(is.null(object$xt$split)) FALSE else object$xt$split
  object <- smooth.construct.la.smooth.spec(object, data, knots)
  object[!(names(object) %in% c("formula", "term", "label", "dim", "X", "xt", "lasso"))] <- NULL
  nodes <- if(split) 1 else object$xt$k
  object$X <- cbind(1, object$X)
  colnames(object$X) <- NULL
  nc <- ncol(object$X) + 1
  npar <- nodes * nc

  nid <- split(1:npar, factor(sort(rep(1:nodes, times = nc))))

  object$fit.fun <- function(X, b, ...) {
    if(!is.null(names(b)))
      b <- get.par(b, "b")
    fit <- 0
    for(j in seq_along(nid)) {
      f <- X %*% b[nid[[j]]][-1]
      fit <- fit + b[nid[[j]]][1] / (1 + exp(-f))
    }
    if(!split)
      fit <- fit - mean(fit, na.rm = TRUE)
    return(fit)
  }

  object$fixed <- FALSE
  object$state$parameters <- rnorm(npar, sd = 0.5)
  for(j in nid)
    object$state$parameters[j[1]] <- rnorm(1, sd = 1e-10)
  names(object$state$parameters) <- paste("b", 1:npar, sep = "")
  object$state$parameters <- c(object$state$parameters, "tau21" = 1000, "tau22" = 1000)
  object$state$fitted.values <- object$fit.fun(object$X, object$state$parameters)
  object$special.npar <- npar
  object$prior <- function(b) { sum(dnorm(get.par(b, "b"), sd = 1000, log = TRUE)) }
  object$nnodes <- nodes

  X <- object$X

  getU <- function(b) {
    if(!is.null(names(b)))
      b <- get.par(b, "b")
    k <- ncol(X)
    U <- matrix(0, nrow(X), length(b))
    i <- 1
    for(j in seq_along(nid)) {
      u <- X %*% b[nid[[j]]][-1]
      o <- 1 / (1 + exp(-u))
      U[, i] <- o
      i <- i + 1
      o2 <- b[nid[[j]]][1] * o * (1 - o)
      U[, i] <- o2
      i <- i + 1
      for(jj in 2:k) {
        U[, i] <- o2 * X[, jj]
        i <- i + 1
      }
    }
    return(U)
  }

  object$S <- list()
  object$S[[1]] <- diag(npar)
  object$S[[2]] <- function(parameters) {
    UU <- diag(crossprod(getU(parameters)))
    diag(UU)
  }
  prior <- make.prior(object)
  object[names(prior)] <- prior
  object$sparse.setup <- FALSE

  update_nn <- function(x, family, y, eta, id, weights, criterion, ...)
  {
    args <- list(...)
  
    peta <- family$map2par(eta)
  
    if(is.null(args$hess)) {
      hess <- process.derivs(family$hess[[id]](y, peta, id = id, ...), is.weight = TRUE)
    } else hess <- args$hess
  
    if(!is.null(weights))
      hess <- hess * weights
  
    if(is.null(args$z)) {
      score <- process.derivs(family$score[[id]](y, peta, id = id, ...), is.weight = FALSE)
      z <- eta[[id]] + 1 / hess * score
    } else z <- args$z
  
    e <- z - eta[[id]]
    eta[[id]] <- eta[[id]] - fitted(x$state)

    b0 <- get.state(x, "b")
    nb <- names(b0)

    U <- getU(b0)

    tUW <- t(U / hess)
    UWU <- tUW %*% U

    edf0 <- args$edf - x$state$edf

    objfun <- function(tau2) {
      H <- matrix_inv(UWU + 1/tau2[1] * x$S[[1]] + 1/tau2[2] * x$S[[2]](b0))
      b0 <- drop(b0 + H %*% tUW %*% e)
      names(b0) <- nb
      edf <- sum_diag(UWU %*% H)
      eta[[id]] <- eta[[id]] + x$fit.fun(x$X, b0)
      ic <- get.ic(family, y, family$map2par(eta), edf0 + edf, length(z), criterion)
      ic
    }

    tau2 <- tau2.optim(objfun, start = get.state(x, "tau2"))

    H <- matrix_inv(UWU + 1/tau2[1] * x$S[[1]] + 1/tau2[2] * x$S[[2]](b0))
    b0 <- drop(b0 + H %*% tUW %*% e)
    names(b0) <- nb

    x$state$parameters <- set.par(x$state$parameters, b0, "b")
    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    x$state$fitted.values <- x$fit.fun(x$X, b0)
    x$state$edf <- sum_diag(UWU %*% H)

    return(x$state)
  }

#  propose_nn <- function(family, theta, id, eta, y, data, weights = NULL, ...) {
#  }

  tau2 <- get.par(object$state$parameters, "tau2")

  UU <- crossprod(getU(object$state$parameters))
  df <- min(c(length(object$state$parameters), 4))

  S1 <- diag(ncol(UU))
  S2 <- diag(diag(UU))

  objfun <- function(tau2) {
    H <- matrix_inv(UU + 1/tau2[1] * S1 + 1/tau2[2] * S2)
    edf <- sum_diag(UU %*% H)
    (df - edf)^2
  }

  opt <- optim(tau2, objfun, method = "L-BFGS-B", lower = rep(.Machine$double.eps^0.5, 2), upper = rep(Inf, 2))
  tau2 <- opt$par
  H <- matrix_inv(UU + 1/tau2[1] * diag(ncol(UU)) + 1/tau2[2] * diag(diag(UU)))
  edf <- sum_diag(UU %*% H)
  object$state$parameters <- set.par(object$state$parameters, tau2, "tau2")
  object$state$edf <- edf
  object$state$special <- object$state$parameters

  object$update <- update_nn

  object$boost.fit <- function(x, y, nu, ...)
  {
    b0 <- get.par(x$state$special, "b")
    tau2 <- get.par(x$state$parameters, "tau2")
    nb <- names(b0)

    U <- getU(b0)
    UU <- crossprod(U, U)
    H <- matrix_inv(UU + 1/1000000 * x$S[[1]] + 1/1000000 * x$S[[2]](b0))
    b1 <- nu * drop(b0 + H %*% t(U) %*% (y - x$fit.fun(x$X, b0)))

    names(b1) <- nb
    fit <- x$fit.fun(x$X, b1)

    x$state$parameters <- set.par(x$state$parameters, b1, "b")
    x$state$fitted.values <- fit
    x$state$rss <- sum((y - fit)^2)
    x$state$special <- x$state$parameters

    return(x$state)
  }

  class(object) <- c("nnet.smooth", "no.mgcv", "special")

  if(split) {
    nodes <- object$xt$k
    object <- rep(list(object), nodes)
    for(j in seq_along(object)) {
      b <- rnorm(length(get.par(object[[j]]$state$parameters, "b")), sd = 0.5)
      b[1] <- rnorm(1, sd = 1e-10)
      names(b) <- paste("b", 1:length(b), sep = "")
      object[[j]]$state$parameters <- set.par(object[[j]]$state$parameters, b, "b")
      object[[j]]$state$fitted.values <- object[[j]]$fit.fun(object[[j]]$X, object[[j]]$state$parameters)
      object[[j]]$state$special <- object[[j]]$state$parameters
      object[[j]]$term_id <- object[[j]]$label
      lab <- strsplit(object[[j]]$label, "")[[1]]
      lab <- paste(c(lab[-length(lab)], paste(',node="', j, '")', sep = '')), collapse = "", sep = "")
      object[[j]]$label <- lab
    }
    class(object) <- c("nnet.smooth", "no.mgcv", "special", "smooth.list")
  }

  object
}

Predict.matrix.nnet.smooth <- function(object, data)
{
  X <- cbind(1, Predict.matrix.lasso.smooth(object, data))
  colnames(X) <- NULL
  X
}

fit_nn <- function(X, y, weights = NULL, start = NULL,
  nnodes = NULL, eps = .Machine$double.eps^0.25)
{
  if(is.null(nnodes))
    nnodes <- 10

  if(is.vector(X))
    X <- matrix(X, ncol = 1)
  X <- cbind(1, X)
  k <- ncol(X)
  nr <- nrow(X)

  if(is.null(weights))
    weights <- rep(1, length(y))

  nc <- k + 1
  npar <- nc * nnodes
  nid <- split(1:npar, factor(sort(rep(1:nnodes, times = nc))))
  b <- rnorm(npar, sd = 0.5)
  for(j in nid)
    b[j[1]] <- rnorm(1, sd = 1e-10)
  names(b) <- paste("b", 1:length(b), sep = "")

  getU <- function(b) {
    U <- matrix(0, nr, npar)
    i <- 1
    for(j in seq_along(nid)) {
      u <- X %*% b[nid[[j]]][-1]
      o <- 1 / (1 + exp(-u))
      U[, i] <- o
      i <- i + 1
      o2 <- b[nid[[j]]][1] * o * (1 - o)
      U[, i] <- o2
      i <- i + 1
      for(jj in 2:k) {
        U[, i] <- o2 * X[, jj]
        i <- i + 1
      }
    }
    return(U)
  }

  fit.fun <- function(b) {
    fit <- 0
    for(j in seq_along(nid)) {
      f <- X %*% b[nid[[j]]][-1]
      fit <- fit + b[nid[[j]]][1] / (1 + exp(-f))
    }
    return(fit - mean(fit))
  }

  I <- diag(npar)
  f0 <- fit.fun(b)

  for(i in 1:1) {
    U <- getU(b)
    tUW <- t(U / weights)
    UWU <- tUW %*% U
    H <- matrix_inv(UWU + 1/10000000000 * I + 1/10000000000 * diag(diag(UWU)))
    b <- drop(b + H %*% tUW %*% (y - fit.fun(b)))
  }
  
  plot(y ~ x)
  plot2d(fit.fun(b) ~ x, add = TRUE, col.lines = 2)
  plot2d(f0 ~ x, add = TRUE, col.lines = 4)
}


if(FALSE) {
  set.seed(123)
  nobs <- 500
  x <- runif(nobs, -3, 3)
  y <- sin(x) + rnorm(nobs, sd = 0.1)
  plot(x, y)

  b <- fit_nn(x, y)
}


## Penalized harmonic smooth.
smooth.construct.ha.smooth.spec <- function(object, data, knots) 
{
  x <- data[[object$term]]

  freq <- if(is.null(object$xt$frequency)) as.integer(max(x, na.rm = TRUE))
  stopifnot(freq > 1 && identical(all.equal(freq, round(freq)), TRUE))

  if(length(object$p.order) < 2) {
    if(is.na(object$p.order))
      object$p.order <- c(2, 2)
    else
      object$p.order <- c(object$p.order, 2)
  }
  object$p.order[is.na(object$p.order)] <- 2

  order <- object$p.order[1]
  order <- min(freq, order)
  x <- x / freq
  X <- outer(2 * pi * x, 1:order)
  X <- cbind(apply(X, 2, cos), apply(X, 2, sin))
  colnames(X) <- if(order == 1) {
    c("cos", "sin")
  } else {
    c(paste("cos", 1:order, sep = ""), paste("sin", 1:order, sep = ""))
  }
  if((2 * order) == freq) X <- X[, -(2 * order)]
  object$X <- X

  gsin1 <- function(x) { cos(2 * pi * order * x) * 2 *pi *order }
  gsin2 <- function(x) { 4 * pi^2 * order^2 * -sin(2 * pi * order * x) }
  gcos1 <- function(x) { -sin(2 * pi * order * x) * 2 * pi * order }
  gcos2 <- function(x) { -4 * pi^2 * order^2 * cos(2 * pi * order * x) }

  if(!object$fixed) {
#    S <- outer(2 * pi * x, 1:order)
#    S <- if(object$p.order[2] < 2) {
#      cbind(apply(S, 2, gcos1), apply(S, 2, gsin1))
#    } else cbind(apply(S, 2, gcos2), apply(S, 2, gsin2))
#    object$S <- list(diag(rep(order:1, 2)))
    K <- t(diff(diag(order))) %*% diff(diag(order))
    K <- rbind(cbind(K, matrix(0, order, order)), cbind(matrix(0, order, order), K))
    object$S <- list(K)
  } else object$S <- list(diag(0, ncol(X)))

  object$frequency <- freq
  object$bs.dim <- ncol(X)
  object$rank <- qr(object$S[[1]])$rank
  object$null.space.dim <- ncol(X)
  object$C <- matrix(nrow = 0, ncol = ncol(X))
#  object$no.rescale <- 1
#  object$side.constrain <- FALSE
  class(object) <- "harmon.smooth"
  object
}


Predict.matrix.harmon.smooth <- function(object, data, knots)
{
  x <- data[[object$term]]
  x <- x / object$frequency
  order <- object$p.order[1]
  X <- outer(2 * pi * x, 1:order)
  X <- cbind(apply(X, 2, cos), apply(X, 2, sin))
  colnames(X) <- if (order == 1) {
    c("cos", "sin")
  } else {
    c(paste("cos", 1:order, sep = ""), paste("sin", 1:order, sep = ""))
  }
  if((2 * order) == object$frequency) X <- X[, -(2 * order)]
  X
}


## Kriging smooth constructor.
## Evaluate a kriging
## design and penalty matrix.
krDesign1D <- function(z, knots = NULL, rho = NULL,
  phi = NULL, v = NULL, c = NULL, ...)
{
  rho <- if(is.null(rho)) {
    geoR::matern
  } else rho
  knots <- if(is.null(knots)) sort(unique(z)) else knots
  v <- if(is.null(v)) 2.5 else v
  c <- if(is.null(c)) {
    optim(1, geoR::matern, phi = 1, kappa = v, method = "L-BFGS-B", lower = 1e-10)$par
  } else c
  phi <- if(is.null(phi)) max(abs(diff(range(knots)))) / c else phi
  B <- NULL
  K <- as.matrix(dist(knots, diag = TRUE, upper = TRUE))
  for(j in seq_along(knots)) {
    h <- abs(z - knots[j])
    B <- cbind(B, rho(h, phi, v))
    K[, j] <- rho(K[, j], phi, v)
  }
  return(list("B" = B, "K" = K, "phi" = phi, "v" = v, "c" = c, "knots" = knots))
}

krDesign2D <- function(z1, z2, knots = 10, rho = NULL,
  phi = NULL, v = NULL, c = NULL, psi = NULL, delta = 1,
  isotropic = TRUE, ...)
{
  rho <- if(is.null(rho)) {
    geoR::matern
  } else rho
  if(is.null(psi)) psi <- 1
  if(is.null(delta)) delta <- 1
  if(is.null(isotropic)) isotropic <- TRUE
  if(is.null(knots)) knots <- min(c(10, nrow(unique(cbind(z1, z2)))), na.rm = TRUE)
  knots <- if(length(knots) < 2) {
    if(knots == length(z1)) {
      unique(cbind(z1, z2))
    } else {
      fields::cover.design(R = unique(cbind(z1, z2)), nd = knots)
    }
  } else knots
  v <- if(is.null(v)) 2.5 else v
  c <- if(is.null(c)) {
    optim(1, rho, phi = 1, kappa = v,
      method = "L-BFGS-B", lower = 1e-10)$par
  } else c
  z <- cbind(z1, z2)
  if(class(knots) == "spatial.design")
    knots <- knots[, 1:2]
  if(!is.matrix(knots))
    knots <- matrix(knots, ncol = 2)
  nk <- nrow(knots)
  phi <- if(is.null(phi)) {
    max(abs(diff(range(knots)))) / c
  } else phi
  if(phi == 0)
    phi <- max(abs(fields::rdist(z1, z2))) / c
  K <- rho(fields::rdist(knots, knots), phi, v)
  if(isotropic) {
    B <- NULL
    for(j in 1:nk) {
      kn <- matrix(knots[j, ], nrow = 1, ncol = 2)
      h <- fields::rdist(z, kn)
      B <- cbind(B, rho(h, phi, v))
    }
  } else {
    B <- matrix(0, nrow(z), nk)
    R <- matrix(c(cos(psi), -1 * sin(psi),
      sin(psi), cos(psi)), 2, 2)
    D <- matrix(c(delta^(-1), 0, 0, 1), 2, 2)
    for(i in 1:nrow(z)) {
      for(j in 1:nk) {
        kn <- matrix(knots[j, ], nrow = 1, ncol = 2)
        h <- as.numeric(z[i, ] - kn)
        h <- drop(sqrt(t(h) %*% t(R) %*% D %*% R %*% h))
        B[i, j] <- rho(h, phi, v)
      }
    }
  }

  return(list("B" = B, "K" = K, "knots" = knots,
    "phi" = phi, "v" = v, "c" = c, "psi" = psi,
    "delta" = delta))
}


## Kriging smooth constructor functions.
smooth.construct.kr.smooth.spec <- function(object, data, knots, ...)
{
  if(object$dim > 2) stop("more than 2 covariates not supported using kriging terms!")
  if(object$bs.dim < 0) object$bs.dim <- 10
  if(object$dim < 2) {
    k <- knots[[object$term]]
    x <- data[[object$term]]
    if(is.null(k))
      k <- seq(min(x), max(x), length = object$bs.dim)
    D <- krDesign1D(x, knots = k, rho = object$xt$rho,
      phi = object$xt$phi, v = object$xt$v, c = object$xt$c)
  } else {
    knots <- if(is.null(object$xt$knots)) object$bs.dim else object$xt$knots
    D <- krDesign2D(data[[object$term[1]]], data[[object$term[2]]],
      knots = knots,
      phi = object$xt$phi, v = object$xt$v, c = object$xt$c,
      psi = object$xt$psi, delta = object$xt$delta,
      isotropic = object$xt$isotropic)
  }

  X <- D$B
  object$X <- X
  object$S <- list(D$K)
  object$rank <- qr(D$K)$rank
  object$knots <- D$knots
  object$null.space.dim <- ncol(D$K)
 
  class(object) <- "kriging.smooth"
  object
}

## Predict function for the new kriging smooth.
Predict.matrix.kriging.smooth <- function(object, data)
{
  if(object$dim < 2) {
    X <- krDesign1D(data[[object$term]], knots = object$knots, rho = object$xt$rho,
      phi = object$xt$phi, v = object$xt$v, c = object$xt$c)$B
  } else {
    X <- krDesign2D(data[[object$term[1]]], data[[object$term[2]]],
      knots = object$knots,
      phi = object$xt$phi, v = object$xt$v, c = object$xt$c,
      psi = object$xt$psi, delta = object$xt$delta,
      isotropic = object$xt$isotropic)$B
  }
  X
}


## Space-time random effect constructor functions.
smooth.construct.str.smooth.spec <- function(object, data, knots)
{
  if(object$dim < 3) stop("need at least 3 variables!")
  if(object$bs.dim < 0) object$bs.dim <- 10

  knots <- if(is.null(object$xt$knots)) object$bs.dim else object$xt$knots

  trend <- data[object$term[3:object$dim]]
  if(length(trend) > 1)
    trend <- trend[[1]] + scale2(trend[[2]], 0, 1)
  trend <- as.vector(trend)

  co0 <- cbind(data[[object$term[1]]], data[[object$term[2]]])
  coid <- match.index(co0)
  co1 <- co0[coid$nodups, ]
  mid <- c(1:nrow(co1))[coid$match.index]
  
  D <- krDesign2D(co1[, 1], co1[, 2],
    knots = knots,
    phi = object$xt$phi, v = object$xt$v, c = object$xt$c,
    psi = object$xt$psi, delta = object$xt$delta,
    isotropic = object$xt$isotropic)

  object$X <- D$B %*% chol2inv(chol(D$K))
  b <- list()
  time <- sort(unique(trend))
  b <- rep(list(rep(0, length = length(time))), length = ncol(D$B))
  b <- do.call("cbind", b)
  rownames(b) <- paste("t", time, sep = "")
  colnames(b) <- paste("k", 1:ncol(b), sep = "")

  object$fit.fun <- function(X, b, ...) {
    fit <- apply(b$b, 1, function(g) {
      X %*% g
    })
    fit <- as.numeric(fit)
print(fit)
    fit
  }

  object$prior <- function(parameters) {
    b <- parameters$b
    tau <- parameters$tau
    print(b)
stop()
  }

  object$update <- bfit_optim

  object$knots <- D$knots
  object$state <- list("parameters" = list("b" = b, "tau2" = c(0.001, 0.001)),
    "fitted.values" = rep(0, length(trend)))
 
  class(object) <- c("strandom.smooth", "no.mgcv", "special")
  object
}


Predict.matrix.strandom.smooth <- function(object, data, knots) 
{
  D <- krDesign2D(data[[object$term[1]]], data[[object$term[2]]],
    knots = object$knots,
    phi = object$xt$phi, v = object$xt$v, c = object$xt$c,
    psi = object$xt$psi, delta = object$xt$delta,
    isotropic = object$xt$isotropic)
  return(D$X %*% chol2inv(chol(D$K)))
}


## Smooth constructor for lag function.
## (C) Viola Obermeier; Flexible distributed lags for modelling earthquake data,
##                      DOI: 10.1111/rssc.12077.
smooth.construct.fdl.smooth.spec <- function(object, data, knots)
{
  ## Modify object so that it's fitted as a p-spline signal regression term.
  object$bs <- "ps"
  object <- smooth.construct.ps.smooth.spec(object, data, knots)

  if(!is.null(object$xt$fullrankpen) && object$xt$fullrankpen){
    ## Add ridge penalty to first <order of B-spline>+1 (=m+2) basis functions.
    ## With same variance as difference penalty: penalty = lambda * coef' (DiffPen + RidgePen) coef.
    object$S[[1]][cbind(1:(object$m[1]+2), 1:(object$m[1]+2))] <- object$S[[1]][cbind(1:(object$m[1]+2), 1:(object$m[1]+2))] + 1
    object$rank <- min(object$bs.dim, object$rank + object$m[1]+2)
  }
  if(!is.null(object$xt$ridge) && object$xt$ridge){
    ## Add ridge penalty to first <order of B-spline>+1 (=m+2) basis functions
    ## penalty = coef' (lambda_1*DiffPen + lambda_2*RidgePen) coef.
    object$S[[2]] <- matrix(0, object$bs.dim, object$bs.dim)
    object$S[[2]][cbind(1:(object$m[1]+2), 1:(object$m[1]+2))] <- 1
    object$rank <- c(object$rank, object$m[1]+2)
  }
  if(!is.null(object$xt$constrain) && object$xt$constrain){
    ## Optionally one can constrain the last lag coefficient to be zero,
    ## not recommended as we favor a soft, data-driven shrinkage to a hard constraint!
    ## Constrain to end in zero (i.e (X%*%coefficients)[1] == 0).
    ## --> Constraint matric C = X[1,]
    object$C <- matrix(object$X[1,],nrow=1)
    object$C <- structure(object$C, always.apply=TRUE)
  }

  return(object)
}

## gam1 <- gam(y ~ 1 + s(lags, K=15, by=X, bs="fdl",
##   xt=list(ridge=TRUE), data=simul, family="poisson", method="REML")


traceplot2 <- function(theta, n.plot=100, ylab = "", ...) {
  cuq <- Vectorize(function(n, x) {
    as.numeric(quantile(x[1:n],c(.025,.5,.975)))
  }, vectorize.args = "n")
  n.rep <- length(theta)
  plot(1:n.rep, theta, col = "lightgrey", xlab = "iter",
    ylab = ylab, type = "l", ...)
  iter <- round(seq(1, n.rep, length = n.plot + 1)[-1])
  tq <- cuq(iter,theta)
  lines(iter, tq[2,])
  lines(iter, tq[1,], lty = 2)
  lines(iter, tq[3,], lty = 2)
}


## Plotting method for "bamlss" objects.
plot.bamlss <- function(x, model = NULL, term = NULL, which = "effects",
  parameters = FALSE, ask = dev.interactive(), spar = TRUE, ...)
{
  if(spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
  }

  ## What should be plotted?
  which.match <- c("effects", "samples", "hist-resid", "qq-resid",
    "scatter-resid", "max-acf", "param-samples", "boost.summary", "results",
    "max-acf")
  if(!is.character(which)) {
    if(any(which > 8L))
      which <- which[which <= 8L]
    which <- which.match[which]
  } else which <- which.match[pmatch(tolower(which), which.match)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")
  if(length(which) < 2) {
    if(which == "results")
      which <- "effects"
  }
  which <- which[which != "results"]

  ok <- any(c("hist-resid", "qq-resid") %in% which)

  if(length(which) > 1 | ok) {
    which <- which[which %in% c("hist-resid", "qq-resid")]
    res <- residuals.bamlss(x, ...)
    plot(res, which = which, spar = spar, ...)
  } else {
    if(which %in% c("samples", "max-acf")) {
      par <- if(parameters) {
        if(is.null(x$parameters)) NULL else unlist(x$parameters)
      } else NULL
      samps <- samples(x, model = model, term = term, drop = TRUE, combine = TRUE, ...)
      snames <- colnames(samps)
      snames <- snames[!grepl(".p.edf", snames, fixed = TRUE) & !grepl(".accepted", snames, fixed = TRUE)]
      snames <- snames[!grepl("DIC", snames, fixed = TRUE) & !grepl("pd", snames, fixed = TRUE)]
      snames <- snames[!grepl(".model.matrix.edf", snames, fixed = TRUE)]
      samps <- samps[, snames, drop = FALSE]
      if(which == "samples") {
        np <- ncol(samps)
        par(mfrow = if(np <= 4) c(np, 2) else c(4, 2))
        devAskNewPage(ask)
        tx <- as.vector(time(samps))
        for(j in 1:np) {
          traceplot2(samps[, j, drop = FALSE], main = paste("Trace of", snames[j]))
          lines(lowess(tx, samps[, j]), col = "red")
          nu <- length(unique(samps[, j, drop = FALSE]))
          acf(if(nu < 2) jitter(samps[, j, drop = FALSE]) else samps[, j, drop = FALSE], main = paste("ACF of", snames[j]), ...)
        }
      } else {
        snames <- snames[!grepl(".edf", snames, fixed = TRUE)]
        snames <- snames[!grepl(".alpha", snames, fixed = TRUE)]
        snames <- snames[!grepl("logLik", snames, fixed = TRUE)]
        samps <- samps[, snames, drop = FALSE]
        macf <- apply(samps, 2, function(x) { acf(x, plot = FALSE, ...) })
        acfx <- macf[[1]]
        acfx$acf <- array(apply(do.call("rbind", lapply(macf, function(x) { x$acf })), 2, max), dim = c(length(acfx$acf), 1L, 1L))
        args <- list(...)
        if(is.null(args$main))
          args$main <- "Maximum ACF of samples"
        if(is.null(args$xlab))
          args$xlab <- "Lag"
        if(is.null(args$ylab))
          args$ylab <- "ACF"
        getS3method("plot", class = "acf")(acfx, main = args$main, xlab = args$xlab,
          ylab = args$ylab, xlim = args$xlim, ylim = args$ylim)
      }
    }

    if(which == "effects") {
      if(is.null(x$results)) {
        plot(results.bamlss.default(x), model = model, term = term, ask = ask, ...)
      } else plot(x$results, model = model, term = term, spar = spar, ask = ask, ...)
    }

    if(which == "boost.summary") {
      if(!is.null(x$boost.summary))
        plot(x$boost.summary, ...)
    }
  }

  return(invisible(NULL))
}


plot.bamlss.results <- function(x, model = NULL, term = NULL,
  ask = dev.interactive(), scale = 1, spar = TRUE, ...)
{
  args <- list(...)
  cx <- class(x)

  if(spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
  }

  if(!is.null(model)) {
    if(!is.character(model)) {
      if(any(model < 0 | model > length(x)))
        stop("model specified wrong!")
      model <- names(x)[model]
    } else {
      i <- grep2(model, names(x))
      if(!length(i))
        stop("model specified wrong!")
      model <- names(x)[i]
    }
    x <- x[model]
  }

  if(FALSE) {
    ## What should be plotted?
    which.match <- which <- "effects"
    if(!is.character(which)) {
      if(any(which > 8L))
        which <- which[which <= 8L]
      which <- which.match[which]
    } else which <- which.match[pmatch(tolower(which), which.match)]
    if(length(which) > length(which.match) || !any(which %in% which.match))
      stop("argument which is specified wrong!")

    args2 <- args
    args2$object <- x
    res0 <- do.call("residuals.bamlss", delete.args("residuals.bamlss", args2, not = "mstop"))
    ny <- if(is.null(dim(res0))) 1 else ncol(res0)
    if(spar) {
      if(!ask) {
        par(mfrow = n2mfrow(length(which) * ny))
      } else par(ask = ask)
    }
    if(any(which %in% c("scatter-resid", "scale-resid"))) {
      fit0 <- fitted.bamlss(x, type = "parameter", samples = TRUE,
        model = if(ny < 2) 1 else NULL, nsamps = args$nsamps)
    }
    rtype <- args$type
    if(is.null(rtype)) rtype <- "quantile"
    if(rtype == "quantile2") rtype <- "quantile"
    if(rtype == "ordinary2") rtype <- "ordinary"
    for(j in 1:ny) {
      res <- if(ny > 1) res0[, j] else res0
      dropi <- !(res %in% c(Inf, -Inf)) & !is.na(res)
      res <- res[dropi]
      if(any(which %in% c("scatter-resid", "scale-resid"))) {
        fit <- if(ny < 2) {
          if(is.list(fit0)) fit0[[1]] else fit0
        } else fit0[[j]]
      }
      for(w in which) {
        args2 <- args
        if(w == "hist-resid") {
          rdens <- density(res)
          rh <- hist(res, plot = FALSE)
          args2$ylim <- c(0, max(c(rh$density, rdens$y)))
          args2$freq <- FALSE
          args2$x <- res
          args2 <- delete.args("hist.default", args2, package = "graphics")
          if(is.null(args$xlab)) args2$xlab <- paste(if(rtype == "quantile") {
              "Quantile"
            } else "Ordinary", "residuals")
          if(is.null(args$ylab)) args2$ylab <- "Density"
          if(is.null(args$main)) {
            args2$main <- "Histogramm and density"
            if(ny > 1)
              args2$main <- paste(names(res0)[j], args2$main, sep = ": ")
          }
          ok <- try(do.call(get("hist.default"), args2))
          if(!inherits(ok, "try-error"))
            lines(rdens)
          box()
        }
        if(w == "qq-resid") {
          args2$y <- if(rtype == "quantile") (res) else (res - mean(res)) / sd(res)
          args2 <- delete.args("qqnorm.default", args2, package = "stats", not = c("col", "pch"))
          if(is.null(args$main)) {
            args2$main <- "Normal Q-Q Plot"
            if(ny > 1)
              args2$main <- paste(names(res0)[j], args2$main, sep = ": ")
          }
          ok <- try(do.call(qqnorm, args2))
          if(!inherits(ok, "try-error"))
  		      if(rtype == "quantile") abline(0,1) else qqline(args2$y)
        }
        if(w == "scatter-resid") {
          args2$x <- fit[dropi]
          args2$y <- res
          args2 <- delete.args("scatter.smooth", args2, package = "stats", not = c("col", "pch"))
          if(is.null(args$xlab)) args2$xlab <- "Fitted values"
          if(is.null(args$xlab)) args2$ylab <- paste(if(rtype == "quantile") {
              "Quantile"
            } else "Ordinary", "residuals")
          if(is.null(args$xlab)) {
            args2$main <- "Fitted values vs. residuals"
            if(ny > 1)
              args2$main <- paste(names(res0)[j], args2$main, sep = ": ")
          }
          ok <- try(do.call(scatter.smooth, args2))
          if(!inherits(ok, "try-error"))
            abline(h = 0, lty = 2)
        }
        if(w == "scale-resid") {
          args2$x <- fit[dropi]
          args2$y <- sqrt(abs((res - mean(res)) / sd(res)))
          args2 <- delete.args("scatter.smooth", args2, package = "stats", not = c("col", "pch"))
          if(is.null(args$xlab)) args2$xlab <- "Fitted values"
          if(is.null(args$ylab)) args2$ylab <- expression(sqrt(abs("Standardized residuals")))
          if(is.null(args$main)) {
            args2$main <- "Scale-location"
            if(ny > 1)
              args2$main <- paste(names(res0)[j], args2$main, sep = ": ")
          }
          try(do.call(scatter.smooth, args2))
        }
      }
    }
  } else {
    ## Get number of plots.
    get_k_n <- function(x) {
      kn <- c(0, length(x))
      ne <- pterms <- list()
      for(i in 1:kn[2]) {
        if(!any(c("s.effects", "p.effects") %in% names(x[[i]]))) {
          kn <- kn + get_k_n(x[[i]])
        } else {
          ne[[i]] <- if(!is.null(names(x[[i]]$s.effects))) names(x[[i]]$s.effects) else NA
          if(is.null(term))
            pterms[[i]] <- 1:length(ne[[i]])
          else {
            if(is.character(term)) {
              tterm <- NULL
              for(j in term)
                tterm <- c(tterm, grep(j, ne[[i]], fixed = TRUE))
              pterms[[i]] <- if(length(tterm)) tterm else NA
            } else pterms[[i]] <- term[term <= length(ne[[i]])]
          }
          if(!is.null(x[[i]]$s.effects) & length(x[[i]]$s.effects)) {
            kn[1] <- kn[1] + length(na.omit(pterms[[i]]))
          }
        }
      }
      kn
    }

    if(any(c("s.effects", "p.effects") %in% names(x)))
      x <- list(x)

    kn <- get_k_n(x)

    if(kn[1] < 1) on.exit(warning("no terms to plot!"), add = TRUE)

    if(spar & (kn[1] > 0)) {
      if(!ask) {
        if("cbamlss" %in% cx) {
          par(mfrow = c(length(x), kn[1] / length(x)))
        } else par(mfrow = n2mfrow(kn[1]))
      } else par(ask = ask)
    }

    mmain <- if(any(c("h1", "Chain_1") %in% (nx <- names(x)))) TRUE else FALSE
    main <- args$main
    if((is.null(args$main) & mmain) | !is.null(args$mmain)) {
      main <- if(!is.null(args$main)) paste(args$main, nx, sep = "-") else nx
      args$mmain <- TRUE
    }
    if(!is.null(main)) main <- rep(main, length.out = length(x))
    for(i in seq_along(x)) {
      args[c("x", "term", "ask", "scale")] <- list(x[i], term, ask, scale)
      args$main <- if(!is.null(main)) main[i] else NULL
      if(!any(c("s.effects", "p.effects") %in% names(x[[i]]))) {
        do.call("plot.bamlss", args)
      } else {
        args$mmain <- NULL
        do.call(".plot.bamlss.results", args)
      }
    }
  }

  invisible(NULL)
}

.plot.bamlss.results <- function(x, model = NULL, term = NULL,
  ask = FALSE, scale = 1, spar = TRUE, ...)
{
  n <- length(x)
  args <- list(...)

  ## Effect plotting.
  k <- 0; ylim <- NULL
  ylim <- args$ylim
  args$residuals <- if(is.null(args$residuals)) FALSE else args$residuals
  if(!is.null(args$ylim))
    scale <- 0
  ne <- pterms <- list()
  for(i in 1:n) {
    ne[[i]] <- if(!is.null(names(x[[i]]$s.effects))) names(x[[i]]$s.effects) else NA
    if(is.null(term))
      pterms[[i]] <- 1:length(ne[[i]])
    else {
      if(is.character(term)) {
        tterm <- NULL
        for(j in term)
          tterm <- c(tterm, grep(j, ne[[i]], fixed = TRUE))
        pterms[[i]] <- if(length(tterm)) tterm else NA
      } else pterms[[i]] <- term[term <= length(ne[[i]])]
    }
  }
  for(i in 1:n) {
    if(!is.null(x[[i]]$s.effects) & length(na.omit(pterms[[i]])) & length(x[[i]]$s.effects)) {
      k <- k + length(na.omit(pterms[[i]]))
      if(scale > 0) {
        term <- term[1:length(x[[i]]$s.effects)]
        for(e in pterms[[i]]) {
          et <- x[[i]]$s.effects[[e]]
          de <- attr(et, "specs")$dim + 1
          ylim <- c(ylim, range(et[, de:ncol(et)], na.rm = TRUE))
          if(args$residuals) {
            if(!is.null(attr(et, "partial.resids"))) {
              res <- attr(et, "partial.resids")
              ylim <- c(ylim, range(res[, de:ncol(res)], na.rm = TRUE))
            }
          }
        }
      }
    }
  }
  if(k < 1) return(NULL)
  if(scale > 0)
    ylim <- range(ylim, na.rm = TRUE)
  args$residuals <- NULL
  for(i in 1:n) {
    if(!is.null(x[[i]]$s.effects) & length(x[[i]]$s.effects)) {
      for(e in pterms[[i]]) {
        lim <- c("ylim", "zlim")[(attr(x[[i]]$s.effects[[e]], "specs")$dim > 1) * 1 + 1]
        setlim <- FALSE
        if(!is.null(ylim) & is.null(args[[lim]])) {
          args[[lim]] <- ylim
          setlim <- TRUE
        }
        args$x <- x[[i]]$s.effects[[e]]
        do.call("plot.bamlss.effect", args)
        if(setlim) args[[lim]] <- NULL
      }
    }
  }

  return(invisible(NULL))
}


## Generic plotting method for model terms.
plot.bamlss.effect <- function(x, ...) {
  UseMethod("plot.bamlss.effect")
}


## Default model term plotting method.
plot.bamlss.effect.default <- function(x, ...) {
  args <- list(...)

  names(x) <- gsub("Mean", "50%", names(x), fixed = TRUE)

  if(attr(x, "specs")$dim > 1 & inherits(x, "rs.smooth")) {
    if(identical(x[, 1], x[, 2])) {
      cn <- colnames(x)[-2]
      xattr <- attributes(x)
      xattr$specs$dim <- 1
      x <- x[, -2, drop = FALSE]
      xattr$names <- colnames(x) <- cn
      cn <- colnames(xattr$partial.resids)[-2]
      xattr$partial.resids <- xattr$partial.resids[, -2, drop = FALSE]
      colnames(xattr$partial.resids) <- cn
      mostattributes(x) <- xattr
    }
  }

  if(length(terms <- attr(x, "specs")$term) > 2) {
    plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "")
    text(0.5, 0.5, paste("use predict() to plot ", attr(x, "specs")$label, "!", sep = ""))
    box()
    warning(paste("use predict() to plot ", attr(x, "specs")$label, "!", sep = ""))
    return(NULL)
  }

  args$x <- x

  lim <- c("ylim", "zlim")[(attr(x, "specs")$dim > 1) * 1 + 1]
  limNULL <- FALSE
  if(is.null(args[[lim]])) {
    limNULL <- TRUE
    if(all(c("2.5%", "97.5%") %in% names(x)))
      args[[lim]] <- range(x[, c("2.5%", "97.5%")], na.rm = TRUE)
    else {
      if(all("50%" %in% names(x))) {
        args[[lim]] <- range(x[, "50%"], na.rm = TRUE)
      }
    }
    if(!is.null(args$residuals)) {
      if(args$residuals & !is.null(attr(x, "partial.resids")))
        args[[lim]] <- range(c(args[[lim]], attr(x, "partial.resids")[, -1]), na.rm = TRUE)
    }
  }
  if((length(unique(args[[lim]])) < 2) & lim == "zlim") {
    add <- if(args[[lim]][1] == 0) 0.01 else 0.01 * abs(args[[lim]][1])
    args[[lim]] <- c(args[[lim]][1] - add, args[[lim]][1] + add)
  }
  if(!is.null(args$shift))
    args[[lim]] <- args[[lim]] + args$shift
  if((attr(x, "specs")$dim > 1) & inherits(x, "mrf.smooth"))
    attr(x, "specs")$dim <- 1
  if(attr(x, "specs")$dim < 2) {
    if(is.null(args$fill.select))
      args$fill.select <- c(0, 1, 0, 1)
    if(is.null(args$lty) & is.null(args$map))
      args$lty <- c(2, 1, 2)
    if(is.null(args$col.lines))
      args$col.lines <- c(NA, "black", NA)
    if(inherits(x, "random.effect") | inherits(x, "re.smooth.spec") |
       inherits(x, "mrf.smooth.spec") | inherits(x, "mrf.smooth") | is.factor(x[[1]])) {
      if(if(!is.null(args$density)) args$density else FALSE) {
        args$density <- NULL
        if(is.null(args$main))
          args$main <- attr(x, "specs")$label
        args$x <- density(x[, "50%"], na.rm = TRUE)
        if(!limNULL)
          args$xlim <- args$ylim
        do.call("plot", delete.args(plot.density2, args, c("main", "xlim")))
      } else {
        if(!is.null(args$map)) {
          if(inherits(args$map, "bnd") | inherits(args$map, "list"))
            args$map <- list2sp(args$map)
          args$x <- data.frame("x" = as.numeric(x[, grepl("50%", colnames(x), fixed = TRUE)]),
            "ID" = as.character(x[, 1]), stringsAsFactors = FALSE)
          idvar <- NULL
          for(j in names(args$map@data)) {
            if(any(args$x$ID %in% as.character(args$map@data[[j]])))
              idvar <- j
          }
          if(!is.null(idvar))
            names(args$x)[2] <- idvar
          args$id <- if(!is.null(idvar)) idvar else as.character(x[, 1])
          args$xlim <- args$ylim <- NULL
          do.call("plotmap", delete.args("plotmap", args,
            not = c("border", "lwd", "lty", names(formals("colorlegend")), "main", "names", "names_id")))
        } else {
          if(is.null(args$ylab))
            args$ylab <- attr(x, "specs")$label
            args$xlab <- attr(x, "specs")$term
          do.call("plotblock", delete.args("plotblock", args,
            c("xlim", "ylim", "pch", "main", "xlab", "ylab", "lwd", "axes", "add", "scheme")))
        }
      }
    } else {
      do.call("plot2d", delete.args("plot2d", args,
        c("xlim", "ylim", "pch", "main", "xlab", "ylab", "lwd", "axes", "add")))
    }
  } else {
    if(is.null(args$c.select))
      args$c.select <- grep("50%", colnames(x), fixed = TRUE)
    if(!is.null(args$slice)) {
      do.call("sliceplot", delete.args("sliceplot", args,
        c("xlim", "ylim", "zlim", "main", "xlab", "ylab", "col", "lwd", "lty")))
    } else {
      if(inherits(x, "random.effect")) {
        do.call("bamlss_random_plot", args)
      } else {
        specs <- attr(x, "specs")
        isf <- sapply(args$x[, specs$term], is.factor)
        if(any(isf)) {
          args$ylim <- args$zlim
          do.call("bamlss_factor2d_plot", args)
        } else {
          do.call("plot3d", delete.args("plot3d", args,
            c("xlim", "ylim", "zlim", "pch", "main", "xlab", "ylab", "ticktype",
            "zlab", "phi", "theta", "r", "d", "scale", "range", "lrange", "pos", "image.map",
            "symmetric", "border", "lwd")))
        }
      }
    }
  }
}


bamlss_random_plot <- function(x, ...)
{
  term <- attr(x, "specs")$term
  cn <- colnames(x)
  isf <- sapply(x[, term], is.factor)
  plot(x[, "50%"] ~ x[, term[!isf]], type = "n", xlab = term[!isf], ylab = attr(x, "specs")$label)
  id <- x[, term[isf]]
  col <- rainbow_hcl(length(unique(id)))
  ii <- 1
  for(j in unique(id)) {
    d <- subset(x, x[, term[isf]] == j)
    i <- order(d[, term[!isf]])
    lines(d[i, "50%"] ~ d[i, term[!isf]], col = col[ii])
    ii <- ii + 1
  }
  return(invisible(NULL))   
}

bamlss_factor2d_plot <- function(x, ids = NULL, add = FALSE, rug = FALSE, ...)
{
  args <- list(...)
  y <- args$response
  specs <- attr(x, "specs")
  if(is.null(specs)) {
    specs <- list("term" = colnames(x)[1:2],
      label = paste("f(", colnames(x)[1], ",", colnames(x)[2], ")", sep = ""))
  }
  isf <- sapply(x[, specs$term], is.factor)
  xd <- x[, specs$term]
  fx <- unlist(x[, grepl("50", colnames(x), fixed = TRUE)])
  isf <- isf[1:length(specs$term)]
  xlab <- if(is.null(args$xlab)) colnames(xd)[!isf] else args$xlab
  ylab <- if(is.null(args$ylab)) specs$label else args$ylab
  id <- xd[, isf]
  xd <- xd[, !isf]
  if(!is.null(ids)) {
    if(!is.character(ids))
      ids <- levels(id)[as.integer(ids)]
    i <- id %in% ids
    id <- droplevels(id[i])
    xd <- xd[i]
    fx <- fx[i]
    if(!is.null(y))
      y <- y[i]
  }
  xlim <- if(is.null(args$xlim)) range(xd) else args$xlim
  ylim <- if(is.null(args$ylim)) range(fx) else args$ylim
  if(!add) {
    plot(1, 1, type = "n",
      xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, main = args$main)
  }
  col <- if(is.null(args$col)) rainbow_hcl(nlevels(id)) else args$col
  if(is.function(col))
     col <- col(nlevels(id))
  lwd <- if(is.null(args$lwd)) 1 else args$lwd
  lty <- if(is.null(args$lty)) 1 else args$lty
  col <- rep(col, length.out = nlevels(id))
  lwd <- rep(lwd, length.out = nlevels(id))
  lty <- rep(lty, length.out = nlevels(id))
  i <- 1
  for(j in levels(id)) {
    fid <- fx[id == j]
    tid <- xd[id == j]
    o <- order(tid)
    lines(fid[o] ~ tid[o], col = col[i], lwd = lwd[i], lty = lty[i])
    if(!is.null(y))
      points(tid, y[id == j], col = col[i], cex = args$cex, pch = args$pch)
    i <- i + 1
  }
  if(rug) {
    jitter <- if(is.null(args$jitter)) TRUE else args$jitter
    if(jitter)
      xd <- jitter(xd)
    rug(xd, col = args$rug.col)
  }
  return(invisible(NULL))   
}


## Other helping functions.
delete.args <- function(fun = NULL, args = NULL, not = NULL, package = NULL)
{
  if(is.character(fun) & !is.null(package))
    fun <- eval(parse(text = paste(package, paste(rep(":", 3), collapse = ""), fun, sep = "")))
  nf <- names(formals(fun))
  na <- names(args)
  for(elmt in na)
    if(!elmt %in% nf) {
      if(!is.null(not)) {
        if(!elmt %in% not)
          args[elmt] <- NULL
      } else args[elmt] <- NULL
    }

  return(args)
}

delete.NULLs <- function(x.list) 
{
  x.list[unlist(lapply(x.list, length) != 0)]
}


## Model summary functions.
summary.bamlss <- function(object, model = NULL, FUN = NULL, parameters = TRUE, ...)
{
  if(!is.null(object$results)) {
    sfun <- try(get(paste("summary", class(object$results), sep = ".")), silent = TRUE)
    if(!inherits(sfun, "try-error"))
      return(sfun(object, model = model, FUN = FUN, parameters = parameters, ...))
  }
  rval <- list()
  rval$call <- object$call
  rval$family <- object$family
  rval$formula <- object$formula
  if(is.null(FUN)) {
    FUN <- function(x) {
      c("Mean" = mean(x, na.rm = TRUE),
         quantile(x, probs = c(0.025, 0.5, 0.975)))
    }
  }
  rval$model.matrix <- .coef.bamlss(object, model = model, FUN = FUN,
     sterms = FALSE, full.names = FALSE, list = TRUE, parameters = parameters, ...)
  rval$model.matrix <- lapply(rval$model.matrix, function(x) {
    if(!is.matrix(x)) {
      rn <- names(x)
      x <- matrix(x, ncol = 1)
      rownames(x) <- rn
      colnames(x) <- ""
      x
    }
    x
  })
  rval$smooth.construct <- .coef.bamlss(object, model = model, FUN = FUN,
     sterms = TRUE, full.names = FALSE, list = TRUE, parameters = parameters, hyper.parameters = TRUE,
     summary = TRUE, ...)
  rval$model.stats <- object$model.stats
  class(rval) <- "summary.bamlss"
  rval
}

print.summary.bamlss <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("---\n")
  print(x$family, full = FALSE)
  cat("*---\n")
  for(i in names(x$formula)) {
    print.bamlss.formula(x$formula[i])
    if(!is.null(x$model.matrix[[i]])) {
      cat("-\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$model.matrix[[i]], digits = digits)
    }
    if(!is.null(x$smooth.construct) & length(x$smooth.construct)) {
      if(!is.null(x$smooth.construct[[i]])) {
        cat("-\n")
        cat("Smooth terms:\n")
        printCoefmat(x$smooth.construct[[i]], digits = digits)
      }
    }
    cat("---\n")
  }
  if(!is.null(x$model.stats)) {
    if(!length(x$model.stats$sampler))
      x$model.stats$sampler <- NULL
    if(!is.null(x$model.stats$sampler)) {
      cat("Sampler summary:\n-\n")
      k <- 1; ok <- FALSE
      for(j in sort(names(x$model.stats$sampler))) {
        if(length(x$model.stats$sampler[[j]]) < 2) {
          ok <- TRUE
          cat(if(k > 1) " " else "", j, " = ", round(x$model.stats$sampler[[j]], digits), sep = "")
          k <- k + 1
          if(k == 4) {
            k <- 1
            cat("\n")
            ok <- FALSE
          }
        }
      }
      if(ok) {
        cat("\n---\n")
      } else {
        if(!is.null(x$model.stats$optimizer))
          cat("---\n")
      }
    }
    if(!is.null(x$model.stats$optimizer)) {
      cat("Optimizer summary:\n-\n")
      k <- 1
      nmo <- sort(names(x$model.stats$optimizer))
      cl <- sapply(x$model.stats$optimizer, class)
      nmo <- nmo[cl != "matrix"]
      if(length(nmo)) {
        for(j in nmo) {
          if(length(x$model.stats$optimizer[[j]]) < 2) {
            ok <- TRUE
            cat(if(k > 1) " " else "", j, " = ", round(x$model.stats$optimizer[[j]], digits), sep = "")
            k <- k + 1
            if(k == 4) {
              k <- 1
              cat("\n")
              ok <- FALSE
            }
          } else {
            print(x$model.stats$optimizer[[j]], ...)
            ok <- FALSE
          }
        }
        if(ok) cat("\n---\n")
      }
    }
  }
  cat("\n")
  return(invisible(x))
}


## Simple "bamlss" print method.
print.bamlss <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  print(family(x), full = FALSE)
  cat("*---\n")
  print(formula(x))
  if(any(c("logLik", "logPost", "IC", "edf") %in% names(x))) {
    cat("*---\n")
    sep <- ""
    if(!is.null(x$logLik)) {
      cat("logLik =", fmt(x$logLik, width = digits, digits = digits))
      sep <- " "
    }
    if(!is.null(x$logPost)) {
      cat(sep, "logPost =", fmt(x$logPost, width = digits, digits = digits))
      sep <- " "
    }
    if(!is.null(x$IC)) {
      if(!is.null(names(x$IC))) {
        cat(sep, names(x$IC), "=", fmt(x$IC, width = digits, digits = digits))
        sep <- " "
      }
    }
    if(!is.null(x$edf))
      cat(sep, "edf =", fmt(x$edf, width = 4, digits = digits))
    cat("\n")
  }
  return(invisible(NULL))
}


## More extractor functions.
DIC.bamlss <- function(object, ..., samples = TRUE, nsamps = NULL)
{
  object <- c(object, ...)
  rval <- NULL

  if(!samples) {
    for(i in 1:length(object)) {
      xs <- summary(object[[i]])
      n <- attr(xs, "n")
      if(n < 2)
        xs <- list(xs)
      rval <- rbind(rval, data.frame(
        "DIC" = if(is.null(xs[[n]]$DIC)) NA else xs[[n]]$DIC,
        "pd" = if(is.null(xs[[n]]$DIC)) NA else xs[[n]]$pd
      ))
    }
  } else {
    for(i in 1:length(object)) {
      family <- attr(object[[i]], "family")
      if(is.null(family$d)) stop("no d() function available in model family object!")
      y <- model.response2(object[[i]])
      d0 <- -2 * sum(family$d(y, family$map2par(fitted.bamlss(object[[i]], type = "link")), log = TRUE), na.rm = TRUE)
      eta <- fitted.bamlss(object[[i]], type = "link",
        samples = TRUE, nsamps = nsamps,
        FUN = function(x) { x })
      iter <- ncol(eta[[1]])
      d1 <- NULL
      for(j in 1:iter) {
        teta <- NULL
        for(ii in 1:length(eta))
          teta <- cbind(teta, eta[[ii]][, j])
        teta <- as.data.frame(teta)
        names(teta) <- names(eta)
        d1 <- c(d1, -2 * sum(family$d(y, family$map2par(teta), log = TRUE), na.rm = TRUE))
      }
      md1 <- mean(d1)
      pd <- md1 - d0
      dic <- md1 + pd
      rval <- rbind(rval, data.frame(
        "DIC" = dic,
        "pd" = pd
      ))
    }
  }

  Call <- match.call()
  row.names(rval) <- if(nrow(rval) > 1) as.character(Call[-1L]) else ""
  rval
}


logLik.bamlss <- function(object, ..., optimizer = FALSE, samples = FALSE)
{
  Call <- match.call()
  Call <- Call[!(names(Call) %in% c("optimizer", "samples"))]
  mn <- as.character(Call)[-1L]
  object <- list(object, ...)
  mstop <- object$mstop
  if(any(names(object) != "")) {
    i <- names(object) == ""
    object <- object[i]
    mn <- mn[i]
  }
  object <- object[mn != "mstop"]
  ll <- edf <- nobs <- NULL
  if(samples)
    ll <- list()
  for(j in seq_along(object)) {
    if(samples) {
      if(is.null(object[[j]]$samples)) {
        warning(paste("no samples available for object ", mn[j], ", cannot compute logLik!", sep = ""))
        ll[[j]] <- as.mcmc(NA)
      } else {
        ll[[j]] <- mcmc(samplestats(object[[j]], logLik = TRUE),
          start = start(object[[j]]$samples), end = end(object[[j]]$samples),
          thin = thin(object[[j]]$samples))
        samps <- process.chains(object[[j]]$samples, combine = TRUE, drop = TRUE)
        sn <- sapply(strsplit(colnames(samps), ".", fixed = TRUE), function(x) { x[length(x)] })
        if(any(i <- (sn == "edf"))) {
          edf <- mcmc(apply(samps[, i, drop = FALSE], 1, sum, na.rm = TRUE),
            start = start(object[[j]]$samples), end = end(object[[j]]$samples),
            thin = thin(object[[j]]$samples))
          ll[[j]] <- mcmc(cbind("logLik" = ll[[j]], "edf" = edf),
            start = start(object[[j]]$samples), end = end(object[[j]]$samples),
            thin = thin(object[[j]]$samples))
        }
      }
    } else {
      ms <- ms0 <- object[[j]]$model.stats
      ms <- if(is.null(ms$sampler) | optimizer) {
        if(is.null(ms$optimizer) & !optimizer) {
          samplestats(object[[j]])
        } else ms$optimizer
      } else ms$sampler
      if(is.null(ms) & !optimizer)
        ms <- samplestats(object[[j]])
      if(is.null(ms)) {
        warning(paste("no logLik available for model ", mn[j], "!", sep = ""))
      } else {
        if(!("logLik" %in% names(ms))) {
          if(!is.null(ms0$optimizer)) {
            if("logLik" %in% names(ms0$optimizer)) {
              ms <- ms0$optimizer
            }
          }
        }        
        if(!("logLik" %in% names(ms))) {
          dfun <- object[[j]]$family$d
          pred <- predict(object[[j]], type = "parameter", mstop = mstop)
          ms <- list("logLik" = sum(dfun(object[[j]]$y[[1]], pred, log = TRUE), na.rm = TRUE))
        }
        if(!("logLik" %in% names(ms))) {
          warning(paste("no logLik available for model ", mn[j], "!", sep = ""))
        }
        ll <- c(ll, ms$logLik)
        edf <- c(edf, if(is.null(ms$edf)) NA else ms$edf)
        nobs <- c(nobs, if(is.null(ms$nobs)) nrow(object[[j]]$y) else ms$nobs)
      }
    }
  }
  if(!is.null(edf)) {
    if(all(is.na(edf)))
      edf <- NULL
  }
  if(!is.null(ll)) {
    if(samples) {
      names(ll) <- mn[1:length(ll)]
      if(length(ll) > 1)
        rval <- as.mcmc.list(ll)
      else
        rval <- ll[[1]]
    } else {
      rval <- cbind("logLik" = ll, "edf" = edf, "nobs" = nobs)
      row.names(rval) <- if(nrow(rval) > 1) mn[1:nrow(rval)] else ""
    }
  } else rval <- NULL
  rval
}


## Extract model formulas.
formula.bamlss.frame <- formula.bamlss <- function(x, model = NULL, ...)
{
  f <- model.terms(x$formula, model)
  class(f) <- "bamlss.formula"
  return(f)
}

formula.bamlss.terms <- function(x, model, ...)
{
  if(!inherits(x, "list") & !inherits(x, "bamlss.formula")) {
    x <- list(x)
    names(x) <- "formula.1"
  }
  f <- list()
  for(i in names(x)) {
    f[[i]] <- list()
    if(!inherits(x[[i]], "terms")) {
      for(j in names(x[[i]])) {
        f[[i]][[j]] <- list()
        f[[i]][[j]]$formula <- x[[i]][[j]]
        env <- environment(x[[i]][[j]])
        attributes(f[[i]][[j]]$formula) <- NULL
        environment(f[[i]][[j]]$formula) <- env
        vars <- all.vars(x[[i]][[j]])
        response <- response.name(x[[i]][[j]])
        if(all(is.na(response)))
          response <- NULL
        if(!is.null(response)) {
          response <- NULL
          vars <- vars[-1]
        }
        f[[i]][[j]]$fake.formula <- as.formula(paste(response, "~1", if(length(vars)) "+" else NULL,
          paste(vars, collapse = "+")), env = environment(x[[i]][[j]]))
        f[[i]][[j]]$terms <- x[[i]][[j]]
      }
    } else {
      f[[i]]$formula <- x[[i]]
      env <- environment(x[[i]])
      attributes(f[[i]]$formula) <- NULL
      environment(f[[i]]$formula) <- env
      vars <- all.vars(x[[i]])
      response <- response.name(x[[i]])
      if(all(is.na(response)))
        response <- NULL
      if(!is.null(response)) {
        response <- NULL
        vars <- vars[-1]
      }
      f[[i]]$fake.formula <- as.formula(paste(response, "~1", if(length(vars)) "+" else NULL,
        paste(vars, collapse = "+")), env = environment(x[[i]]))
      f[[i]]$terms <- x[[i]]
    }
  }
  class(f) <- c("bamlss.formula", "list")
  environment(f) <- environment(x)
  return(f)
}

print.bamlss.formula <- function(x, ...) {
  if(!inherits(x, "list") & !inherits(x, "bamlss.formula")) {
    print(x)
  } else {
    nx <- names(x)
    if(is.null(nx))
      nx <- as.character(1:length(x))
    for(i in seq_along(x)) {
      cat("Formula ", nx[i], ":\n---\n", sep = "")
      if(inherits(x[[i]], "list") & "h1" %in% names(x[[i]])) {
        for(j in seq_along(x[[i]])) {
          cat("h", j, ": ", sep = "")
          attr(x[[i]][[j]], "name") <- NULL
          attr(x[[i]][[j]]$formula, ".Environment") <- NULL
          print(x[[i]][[j]]$formula, showEnv = FALSE)
        }
      } else {
        attr(x[[i]], "name") <- NULL
        attr(x[[i]]$formula, "name") <- NULL
        attr(x[[i]]$formula, ".Environment") <- NULL
        if("formula" %in% names(x[[i]])) print(x[[i]]$formula, showEnv = FALSE) else print(x[[i]])
      }
      if(i < length(x))
      cat("\n")
    }
  }
  invisible(NULL)
}


## Drop terms from "bamlss.terms'.
drop.terms.bamlss <- function(f, pterms = TRUE, sterms = TRUE,
  specials = NULL, keep.response = TRUE, keep.intercept = TRUE, data = NULL)
{
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "la", "n"))
  if(!inherits(f, "formula")) {
    if(!is.null(f$terms)) {
      f <- f$terms
    } else {
      if(!is.null(f$formula))
        f <- f$formula
    }
  }
  tx <- if(!inherits(f, "terms")) {
    terms.formula(f, specials = specials, keep.order = TRUE, data = data)
  } else f
  specials <- unique(c(names(attr(tx, "specials")), specials))
  tl <- attr(tx, "term.labels")
  sid <- NULL
  for(j in specials)
    sid <- c(sid, grep2(paste(j, "(", sep = ""), tl, fixed = TRUE))
  if(length(sid))
    sid <- sort(unique(sid))
  sub <- attr(tx, "response")
  if(length(sid)) {
    st <- tl[sid]
    pt <- tl[-sid]
  } else {
    st <- character(0)
    pt <- tl
  }
  if(!sterms & length(st)) {
    st <- paste("-", st, collapse = "")
    st <- as.formula(paste(". ~ .", st), env = NULL)
    tx <- terms.formula(update(tx, st), specials = specials, keep.order = TRUE, data = data)
  }
  if(!pterms & length(pt)) {
    tl <- attr(tx, "term.labels")
    sid <- NULL
    for(j in specials)
      sid <- c(sid, grep2(paste(j, "(", sep = ""), tl, fixed = TRUE))
    if(length(sid))
      sid <- sort(unique(sid))
    if(length(sid)) {
      st <- tl[sid]
      pt <- tl[-sid]
    } else {
      st <- character(0)
      pt <- tl
    }
    pt <- paste("-", pt, collapse = "")
    pt <- as.formula(paste(". ~ .", pt), env = NULL)
    tx <- terms.formula(update(tx, pt), specials = specials, keep.order = TRUE, data = data)
  }
  class(tx) <- c("formula", "terms")
  environment(tx) <- environment(f)
  if(!keep.response)
    tx <- delete.response(tx)
  if(!keep.intercept) {
    if(attr(tx, "intercept") > 0)
      tx <- terms.formula(update(tx, . ~ -1 + .), specials = specials, keep.order = TRUE, data = data)
  }
  tx
}

has_dot <- function(formula) {
  inherits(try(terms(formula), silent = TRUE), "try-error")
}

terms.bamlss <- terms.bamlss.frame <- terms.bamlss.formula <- function(x, specials = NULL,
  data = NULL, model = NULL, pterms = TRUE, sterms = TRUE, drop = TRUE, ...)
{
  if(inherits(x, "bamlss.frame"))
    x <- formula(x)
  if(!inherits(x, "bamlss.formula"))
    x <- bamlss.formula(x, ...)
  env <- environment(x)
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "la", "n"))
  elmts <- c("formula", "fake.formula")
  if(!any(names(x) %in% elmts) & !inherits(x, "formula")) {
    if(!is.null(model)) {
      if(is.character(model)) {
        if(all(is.na(pmatch(model[1], names(x)))))
          stop("argument model is specified wrong!")
      } else {
        if(max(model[1]) > length(x) || is.na(model[1]) || min(model[1]) < 1) 
          stop("argument model is specified wrong!")
      }
      if(length(model) > 1)
        model <- model[1:2]
      if(length(model) < 2) {
        x <- x[model]
      } else {
        x <- x[[model[1]]]
        if(is.character(model)) {
          if(all(is.na(pmatch(model[2], names(x)))))
            stop("argument model is specified wrong!")
        } else {
          if(max(model[2]) > length(x) || is.na(model[2]) || min(model[2]) < 1) 
            stop("argument model is specified wrong!")
        }
        x <- x[model[2]]
      }
    }
  } else x <- list(x)

  rval <- list()
  if(is.null(nx <- names(x))) {
    nx <- paste("formula", 1:length(x), sep = ".")
    names(x) <- nx
  }
  for(i in seq_along(nx)) {
    if(!any(names(x[[nx[i]]]) %in% elmts) & !inherits(x[[nx[i]]], "formula")) {
      rval[[nx[i]]] <- list()
      nx2 <- names(x[[nx[i]]])
      for(j in seq_along(nx2)) {
        rval[[nx[i]]][[nx2[j]]] <- drop.terms.bamlss(x[[nx[i]]][[nx2[j]]],
          pterms = pterms, sterms = sterms, specials = specials, data = data)
      }
    } else {
      rval[[nx[i]]] <- drop.terms.bamlss(x[[nx[i]]], pterms = pterms,
        sterms = sterms, specials = specials, data = data)
    }
  }

  if(drop & (length(rval) < 2)) {
    rval <- rval[[1]]
  } else {
    class(rval) <- c("bamlss.terms", "list")
  }
  environment(rval) <- env

  rval
}


## Model terms extractor function for formulas and 'bamlss.frame'.
model.terms <- function(x, model = NULL, part = c("x", "formula", "terms"))
{
  if(!inherits(x, "bamlss.formula")) {
    if(inherits(x, "bamlss.frame")) {
      part <- match.arg(part)
      if(is.null(x[[part]]))
        stop(paste("cannot find object", part, "in 'bamlss.frame' object!"))
      x <- x[[part]]
    } else stop(paste("cannot extract parts from object of class '", class(x), "'!", sep = ""))
  }
  if(is.null(model))
    return(x)
  cx <- class(x)
  env <- environment(x)
  elmts <- c("formula", "fake.formula")
  if(!any(names(x) %in% elmts)) {
    if(is.character(model)) {
      if(all(is.na(pmatch(model[1], names(x)))))
        stop("argument model is specified wrong!")
    } else {
      if(max(model[1]) > length(x) || is.na(model[1]) || min(model[1]) < 1) 
        stop("argument model is specified wrong!")
    }
    if(length(model) > 1)
      model <- model[1:2]
    if(length(model) < 2) {
      x <- x[model]
    } else {
      x <- x[[model[1]]]
      if(is.character(model)) {
        if(all(is.na(pmatch(model[2], names(x)))))
          stop("argument model is specified wrong!")
      } else {
        if(max(model[2]) > length(x) || is.na(model[2]) || min(model[2]) < 1) 
          stop("argument model is specified wrong!")
      }
      x <- x[model[2]]
    }
  } else x <- list(x)
  class(x) <- cx
  environment(x) <- env
  return(x)
}


## Some simple check functions for 'term' objects.
has_intercept <- function(x)
{
  if(inherits(x, "formula"))
    x <- terms(x)
  if(!inherits(x, "terms"))
    stop("x must be a 'terms' object!")
  return(attr(x, "intercept") > 0)
}

has_response <- function(x)
{
  if(inherits(x, "formula"))
    x <- terms(x)
  if(!inherits(x, "terms"))
    stop("x must be a 'terms' object!")
  return(attr(x, "response") > 0)
}

has_sterms <- function(x, specials = NULL)
{
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "la", "n"))
  if(inherits(x, "formula"))
    x <- terms(x, specials = specials)
  if(!inherits(x, "terms"))
    stop("x must be a 'terms' object!")
  return(length(unlist(attr(x, "specials"))) > 0)
}

has_pterms <- function(x, specials = NULL)
{
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "la", "n"))
  if(inherits(x, "formula"))
    x <- terms(x, specials = specials)
  if(!inherits(x, "terms"))
    stop("x must be a 'terms' object!")
  x <- drop.terms.bamlss(x, pterms = TRUE, sterms = FALSE, specials = specials, keep.response = FALSE)
  fc <- length(attr(x, "factors")) > 0
  ic <- attr(x, "intercept") > 0
  return(fc | ic)
}

get_pterms_labels <- function(x, specials = NULL)
{
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "la", "n"))
  tl <- if(has_pterms(x, specials)) {
    x <- drop.terms.bamlss(x, pterms = TRUE, sterms = FALSE,
      keep.response = FALSE, specials = specials)
    c(attr(x, "term.labels"), if(attr(x, "intercept") > 0) "(Intercept)" else NULL)
  } else character(0)
  tl
}

get_sterms_labels <- function(x, specials = NULL)
{
  env <- environment(x)
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "la", "n"))
  if(has_sterms(x, specials)) {
    x <- drop.terms.bamlss(x, pterms = FALSE, sterms = TRUE,
      keep.response = FALSE, specials = specials)
    tl <- all.labels.formula(x)
  } else tl <- character(0)
  tl
}


## Process results with samples and bamlss.frame.
results.bamlss.default <- function(x, what = c("samples", "parameters"), grid = -1, nsamps = NULL,
  burnin = NULL, thin = NULL, ...)
{
  if(!inherits(x, "bamlss.frame") & !inherits(x, "bamlss"))
    stop("x must be a 'bamlss' object!")
  if(is.null(x$samples) & is.null(x$parameters)) {
    warning("nothing to do!")
    return(NULL)
  }

  if(is.null(x$x))
    stop("cannot compute results, 'x' object is missing, see design.construct()!")

  what <- match.arg(what)
  if(!is.null(x$samples) & what == "samples") {
    if(!is.null(list(...)$bamlss)) {
      burnin = NULL; thin <- NULL
    }
    samps <- samples(x, burnin = burnin, thin = thin)
    if(!is.null(nsamps)) {
      i <- seq(1, nrow(samps), length = nsamps)
      samps <- samps[i, , drop = FALSE]
    }
  } else {
    if(is.null(x$parameters)) {
      warning("nothing to do!")
      return(NULL)
    }
    samps <- parameters(x, extract = TRUE, list = FALSE)
    cn <- names(samps)
    samps <- matrix(samps, nrow = 1)
    colnames(samps) <- cn
    samps <- as.mcmc(samps)
  }

  family <- x$family
  snames <- colnames(samps)
  mf <- model.frame(x)

  make_results <- function(obj, id = NULL)
  {
    DIC <- pd <- NA
    if(any(grepl("deviance", snames))) {
      DIC <- as.numeric(samps[, grepl("deviance", snames)])
      pd <- var(DIC, na.rm = TRUE) / 2
      DIC <- mean(DIC, na.rm = TRUE)
    }
    if(any(grepl("logLik", snames))) {
      DIC <- -2 * as.numeric(samps[, grepl("logLik", snames)])
      pd <- var(DIC, na.rm = TRUE) / 2
      DIC <- mean(DIC, na.rm = TRUE)
    }
    IC <- c("DIC" = DIC, "pd" = pd)

    ## Compute model term effects.
    p.effects <- s.effects <- s.effects.resmat <- NULL

    ## Parametric effects.
    if(has_pterms(obj$terms)) {
      tl <- get_pterms_labels(obj$terms)
      sn <- paste(id, "p", tl, sep = ".")
      i <- grep2(sn, snames, fixed = TRUE)
      if(length(i)) {
        psamples <- as.matrix(samps[, snames[i], drop = FALSE])
        nas <- apply(psamples, 1, function(x) { any(is.na(x)) } )
        psamples <- psamples[!nas, , drop = FALSE]
        qu <- t(apply(psamples, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
        sd <- drop(apply(psamples, 2, sd, na.rm = TRUE))
        me <- drop(apply(psamples, 2, mean, na.rm = TRUE))
        p.effects <- cbind(me, sd, qu)
        rownames(p.effects) <- gsub(paste(id, "p.", sep = "."), "", snames[i], fixed = TRUE)
        colnames(p.effects) <- c("Mean", "Sd", "2.5%", "50%", "97.5%")
      }
    }

    ## Smooth effects.
    if(has_sterms(obj$terms)) {
      tl <- names(obj$smooth.construct)
      tl2 <- get_sterms_labels(obj$terms)
      if(length(ib <- grep("by=", tl2, fixed = TRUE))) {
        tl2[ib] <- gsub(")", "", tl2[ib], fixed = TRUE)
      }
      if(length(nn <- grep("n(", tl2, fixed = TRUE))) {
        tl2[nn] <- sapply(strsplit(tl2[nn], ""), function(x) {
          paste(x[-length(x)], collapse = "")
        })
      }
      tl <- tl[grep2(tl2, tl, fixed = TRUE)]
      sn <- paste(id, "s", tl, sep = ".")
      i <- grep2(sn, snames, fixed = TRUE)
      if(length(i)) {
        for(j in tl) {
          sn <- paste(id, "s", j, sep = ".")
          psamples <- as.matrix(samps[, snames[grep2(sn, snames, fixed = TRUE)], drop = FALSE])
          nas <- apply(psamples, 1, function(x) { any(is.na(x)) } )
          psamples <- psamples[!nas, , drop = FALSE]
       
          ## FIXME: retransform!
          if(!is.null(obj$smooth.construct[[j]]$Xf) & FALSE) {
            stop("no randomized terms supported yet!")
            kx <- ncol(obj$smooth.construct[[j]]$Xf)
            if(kx) {
              pn <- paste(paste(id, ":h1:linear.",
                paste(paste(obj$smooth.construct[[j]]$term, collapse = "."), "Xf", sep = "."), sep = ""),
                1:kx, sep = ".")
              xsamps <- matrix(samples[[j]][, snames %in% pn], ncol = kx)
              psamples <- cbind("ra" = psamples, "fx" = xsamps)
              re_trans <- function(g) {
                g <- obj$smooth.construct[[j]]$trans.D * g
                if(!is.null(obj$smooth.construct[[j]]$trans.U))
                  g <- obj$smooth.construct[[j]]$trans.U %*% g
                g
              }
              psamples <- t(apply(psamples, 1, re_trans))
            }
          }

          ## Prediction matrix.
          get.X <- function(x) { ## FIXME: time(x)
            for(char in c("(", ")", "[", "]")) {
              obj$smooth.construct[[j]]$term <- gsub(char, ".", obj$smooth.construct[[j]]$term, fixed = TRUE)
              obj$smooth.construct[[j]]$by <- gsub(char, ".", obj$smooth.construct[[j]]$by, fixed = TRUE)
            }
            if(is.null(obj$smooth.construct[[j]]$PredictMat)) {
              X <- PredictMat(obj$smooth.construct[[j]], x)
            } else {
              X <- obj$smooth.construct[[j]]$PredictMat(obj$smooth.construct[[j]], x)
            }
            X
          }

          ## Compute effect.
          if(!is.list(s.effects))
            s.effects <- list()
          if(length(s.effects)) {
            if(obj$smooth.construct[[j]]$label %in% names(s.effects)) {
              ct <- gsub(".smooth.spec", "", class(obj$smooth.construct[[j]]))[1]
              if(ct == "random.effect") ct <- "re"
              obj$smooth.construct[[j]]$label <- paste(obj$smooth.construct[[j]]$label, ct, sep = ".")
            }
          }
          if(is.null(obj$smooth.construct[[j]]$fit.fun)) {
            obj$smooth.construct[[j]]$fit.fun <- function(X, b, ...) {
              drop(X %*% b)
            }
          }

          if(is.null(obj$smooth.construct[[j]][["X"]])) {
            b <- paste(id, "s", j, paste("b", 1:obj$smooth.construct[[j]][["X.dim"]], sep = ""), sep = ".")
          } else {
            b <- paste(id, "s", j,
              if(is.null(colnames(obj$smooth.construct[[j]]$X))) {
                if(!inherits(obj$smooth.construct[[j]], "special")) {
                  paste("b", 1:ncol(obj$smooth.construct[[j]]$X), sep = "")
                } else {
                  npar  <- if(inherits(obj$smooth.construct[[j]], "rs.smooth")) {
                    names(get.par(obj$smooth.construct[[j]]$state$parameters, "b"))
                  } else {
                    npar <- if(!is.null(obj$smooth.construct[[j]]$state$parameters)) {
                      length(get.state(obj$smooth.construct[[j]], "b"))
                    } else {
                      ncol(obj$smooth.construct[[j]]$X)
                    }
                    paste("b", 1:npar, sep = "")
                  }
                }
              } else colnames(obj$smooth.construct[[j]]$X), sep = ".")
          }

          tn <- c(obj$smooth.construct[[j]]$term, if(obj$smooth.construct[[j]]$by != "NA") {
            obj$smooth.construct[[j]]$by
          } else NULL)

          if(!all(ii <- tn %in% names(mf))) {
            ii <- tn[which(!ii)]
            take <- NULL  ## FIXME: by dummies!
          }

          if(!any(b %in% colnames(psamples))) {
            b <- grep(paste(id, "s", j, "", sep = "."), colnames(psamples), fixed = TRUE, value = TRUE)
            if(length(drop <- grep2(c("tau2", "edf", "alpha", "hyper"), colnames(psamples), fixed = TRUE)))
              b <- b[-drop]
          }

          s.effects[[obj$smooth.construct[[j]]$label]] <- compute_s.effect(obj$smooth.construct[[j]],
            get.X = get.X, fit.fun = obj$smooth.construct[[j]]$fit.fun, psamples = psamples[, b, drop = FALSE],
            FUN = NULL, snames = snames, data = mf[, tn, drop = FALSE], grid = grid)
        }
      }
    }

    rval <- list(
      "model" = list("formula" = obj$formula,
        "DIC" = DIC, "pd" = pd, "N" = nrow(mf)),
      "p.effects" = p.effects, "s.effects" = s.effects
    )

    class(rval) <- "bamlss.results"
    return(rval)
  }

  rval <- list()
  nx <- names(x$x)
  for(j in nx) {
    rval[[j]] <- make_results(x$x[[j]], id = j)
    if(!is.null(rval[[j]]$s.effects)) {
      for(i in seq_along(rval[[j]]$s.effects)) {
        specs <- attr(rval[[j]]$s.effects[[i]], "specs")
        specs$label <- paste(specs$label, j, sep = ".")
        attr(rval[[j]]$s.effects[[i]], "specs") <- specs
      }
    }
  }

  class(rval) <- "bamlss.results"
  return(rval)
}


## Fitted values/terms extraction
fitted.bamlss <- function(object, model = NULL, term = NULL,
  type = c("link", "parameter"), samples = TRUE, FUN = c95,
  nsamps = NULL, ...)
{
  type <- match.arg(type)
  if(!samples & !is.null(object$fitted.values)) {
    if(!is.null(term))
      stop("term specific fitted values must be computed with 'samples = TRUE'!")
    return(if(is.null(model)) {
        object$fitted.values
      } else {
        if(length(model) < 2) {
          object$fitted.values[[model]]
        } else object$fitted.values[model]})
  } else {
    return(predict.bamlss(object, model = model, term = term,
      type = type, FUN = FUN, nsamps = nsamps, ...))
  }
}

## Functions for model samples
grep2 <- function(pattern, x, ...) {
  i <- NULL
  for(p in pattern)
    i <- c(i, grep(p, x, ...))
  sort(unique(i))
}

samples <- function(object, ...)
{
  UseMethod("samples")
}

samples.bamlss <- samples.bamlss.frame <- function(object, model = NULL, term = NULL, combine = TRUE, drop = TRUE,
  burnin = NULL, thin = NULL, ...)
{
  if(!inherits(object, "bamlss") & !inherits(object, "bamlss.frame"))
    stop("object is not a 'bamlss' object!")
  if(is.null(object$samples))
    stop("no samples to extract!")
  tx <- terms(object, drop = FALSE)
  x <- object$samples
  x <- process.chains(x, combine, drop = FALSE)

  snames <- colnames(x[[1]])
  nx <- names(tx)

  if(!is.null(model)) {
    model <- model[1]
    i <- if(is.character(model)) {
      pmatch(model, nx)
    } else {
      if(length(model) > length(tx)) NA else model
    }
    if(is.na(i))
      stop("cannot find model!")
    j <- grep(paste(nx[i], ".", sep = ""), snames, fixed = TRUE, value = TRUE)
    for(k in seq_along(x)) {
      x[[k]] <- x[[k]][, j]
    }
    tx <- tx[nx[i]]
    snames <- colnames(x[[1]])
  }

  if(!is.null(term)) {
    term <- term[1]
    if(!is.character(term)) {
      if(term < 1)
        term <- "(Intercept)"
    }
    rval <- vector(mode = "list", length = length(x))
    nx <- names(tx)
    for(i in seq_along(tx)) {
      tl <- all.labels.formula(tx[[i]], full.names = TRUE)
      if(attr(tx[[i]], "intercept") > 0)
        tl <- c(tl, "(Intercept)")
      if(is.character(term)) {
        j <- grep(term, tl, fixed = TRUE)
        if(!length(j))
          j <- NA
      } else {
        j <- if(length(term) > length(tl)) NA else term
      }
      if(is.na(j))
        next
      jj <- grep(tl[j], snames, fixed = TRUE, value = TRUE)
      for(ii in jj) {
        for(k in seq_along(x))
          rval[[k]] <- cbind(rval[[k]], x[[k]][, ii, drop = FALSE])
      }
    }
    for(k in seq_along(x))
      rval[[k]] <- as.mcmc(rval[[k]], start = start(x[[k]]), end = end(x[[k]]))
    x <- as.mcmc.list(rval)
  }

  if(!is.null(burnin)) {
    for(i in seq_along(x)) {
      x[[i]] <- mcmc(x[[i]][burnin:nrow(x[[i]]), , drop = FALSE], start = burnin)
    }
  }
  if(!is.null(thin)) {
    iterthin <- as.integer(seq(1, nrow(x[[1]]), by = thin))
    for(i in seq_along(x)) {
      x[[i]] <- mcmc(x[[i]][iterthin, , drop = FALSE],
        start = if(!is.null(burnin)) burnin else 1, thin = thin)
    }
  }

  if(drop & (length(x) < 2))
    x <- x[[1]]

  return(x)
}


## Continue sampling.
continue <- function(object, cores = NULL, combine = TRUE,
  sleep = NULL, results = TRUE, ...)
{
  if(is.null(object$samples))
    stop("no samples to continue from!")
  start <- drop(tail(process.chains(object$samples, combine = TRUE, drop = TRUE), 0))
  i <- grep2(c(".edf", ".alpha", ".accepted", "logLik", "DIC"), names(start), fixed = TRUE)
  start <- start[-i]

  sampler <- attr(object, "functions")$sampler
  results <- if(results) attr(object, "functions")$results else FALSE

  if(is.null(cores)) {
    samples <- sampler(x = object$x, y = object$y, family = object$family,
      weights = model.weights(object$model.frame),
      offset = model.offset(object$model.frame),
      start = start, hessian = object$hessian, ...)
  } else {
    parallel_fun <- function(j) {
      if(j > 1 & !is.null(sleep)) Sys.sleep(sleep)
      sampler(x = object$x, y = object$y, family = object$family,
        weights = model.weights(object$model.frame),
        offset = model.offset(object$model.frame), start = start,
        hessian = object$hessian, ...)
    }
    samples <- parallel::mclapply(1:cores, parallel_fun, mc.cores = cores)
  }
  if(!inherits(samples, "mcmc")) {
    if(is.list(samples)) {
      samples <- as.mcmc.list(lapply(samples, as.mcmc))
    } else {
      samples <- as.mcmc(samples)
    }
  }

  ## Process samples.
  samples <- process.chains(samples, TRUE)
  object$samples <- process.chains(c(object$samples, samples), combine = combine)

  ## Compute results.
  if(is.function(results))
    object$results <- try(results(object, bamlss = TRUE,  ...))

  return(object)
}


## Credible intervals of coefficients.
confint.bamlss <- function(object, parm, level = 0.95, model = NULL,
  pterms = TRUE, sterms = FALSE, full.names = FALSE, hyper.parameters = FALSE, ...)
{
  args <- list(...)
  if(!is.null(args$term))
    parm <- args$term
  if(missing(parm))
    parm <- NULL
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  FUN <- function(x) {
    quantile(x, probs = probs, na.rm = TRUE)
  }
  return(.coef.bamlss(object, model = model, term = parm,
    FUN = FUN, parameters = FALSE, pterms = pterms, sterms = sterms,
    full.names = full.names, hyper.parameters = hyper.parameters, ...))
}


## Extract model coefficients.
coef.bamlss <- function(object, model = NULL, term = NULL,
  FUN = NULL, parameters = NULL, pterms = TRUE, sterms = TRUE,
  hyper.parameters = TRUE, list = FALSE, full.names = TRUE, rescale = FALSE, ...)
{
  .coef.bamlss(object, model = model, term = term,
    FUN = FUN, parameters = parameters, pterms = pterms, sterms = sterms,
    s.variances = TRUE, hyper.parameters = hyper.parameters,
    summary = FALSE, list = list, full.names = full.names, rescale = rescale, ...)
}

.coef.bamlss <- function(object, model = NULL, term = NULL,
  FUN = NULL, parameters = NULL, pterms = TRUE, sterms = TRUE,
  s.variances = FALSE, hyper.parameters = FALSE, summary = FALSE,
  list = FALSE, full.names = TRUE, rescale = FALSE, ...)
{
  if(is.null(object$samples) & is.null(object$parameters))
    stop("no coefficients to extract!")
  if(is.null(parameters))
    parameters <- is.null(object$samples)
  if(hyper.parameters) {
    pterms <- if(s.variances) TRUE else FALSE
    if(summary) {
      drop <- c(".accepted", "logLik", "logPost", "AIC", "BIC", "DIC", "pd")
    } else {
      drop <- c(".accepted", "logLik", "logPost", "AIC", "BIC", "DIC", "pd", ".edf")
    }
    if(is.null(FUN)) {
      FUN <- function(x) {
        c("Mean" = mean(x, na.rm = TRUE),
          quantile(x, probs = c(0.025, 0.5, 0.975)))
      }
    }
  } else {
    drop <- c(".tau2", ".lambda", ".edf", ".accepted",
      ".alpha", "logLik", "logPost", "AIC", "BIC", "DIC", "pd")
    if(is.null(FUN))
      FUN <- function(x) { mean(x, na.rm = TRUE) }
  }
  if(!pterms)
    drop <- c(drop, ".p.")
  if(!sterms)
    drop <- c(drop, ".s.") 
  par <- samps <- NULL
  rval <- list()
  if(!is.null(object$samples)) {
    rval$samples <- samples(object, model = model, term = term, ...)
    tdrop <- grep2(drop, colnames(rval$samples), fixed = TRUE)
    if(length(tdrop))
      rval$samples <- rval$samples[, -tdrop, drop = FALSE]
    if(hyper.parameters & summary) {
      ttake <- grep2(c(".tau2", ".lambda", ".edf", ".alpha"), colnames(rval$samples), fixed = TRUE)
      if(length(ttake)) {
        rval$samples <- rval$samples[, ttake, drop = FALSE]
      } else rval$samples <- numeric(0)
    }
    if(length(rval$samples)) {
      rval$samples <- apply(rval$samples, 2, function(x) { FUN(na.omit(x), ...) })
      rval$samples <- if(!is.null(dim(rval$samples))) {
        t(rval$samples)
      } else {
        as.matrix(rval$samples, ncol = 1)
      }
      if(is.null(colnames(rval$samples))) {
        fn <- deparse(substitute(FUN), backtick = TRUE, width.cutoff = 500)
        colnames(rval$samples) <- rep(fn, length = ncol(rval$samples))
      }
    }
  }
  if(!is.null(object$parameters) & parameters) {
    rval$parameters <- parameters(object, list = FALSE, ...)
    if(length(di <- grep2(drop, names(rval$parameters), fixed = TRUE)))
      rval$parameters <- rval$parameters[-di]
    if(summary)
      rval$parameters <- rval$parameters[grep2(c(".tau2", ".edf"), names(rval$parameters), fixed = TRUE)]
    rval$parameters <- as.matrix(rval$parameters, ncol = 1)
    if(!is.null(rval$samples) & length(rval$samples)) {
      pc <- NULL
      rns <- gsub(".model.matrix", "", rownames(rval$samples), fixed = TRUE)
      for(j in rns) {
        if(j %in% rownames(rval$parameters)) {
          pc <- rbind(pc, rval$parameters[j, , drop = FALSE])
        } else {
          tpc <- matrix(NA, nrow = 1, ncol = ncol(rval$parameters))
          rownames(tpc) <- j
          pc <- rbind(pc, tpc)
        }
      }
      rval$parameters <- pc
    }
    if((ncol(rval$parameters) > 1) & !is.null(list(...)$mstop))
      return(rval$parameters)
    colnames(rval$parameters) <- "parameters"
  }
  if(!length(rval)) return(NULL)
  rval <- if(length(rval) < 2) {
    as.matrix(rval[[1]], ncol = 1)
  } else {
    do.call("cbind", rval)
  }
  if(!length(rval)) return(numeric(0))
  nx <- sapply(strsplit(rownames(rval), ".", fixed = TRUE), function(x) { x[1] })
  if(list) {
    rval2 <- list()
    for(i in unique(nx)) {
      rval2[[i]] <- rval[nx == i, , drop = FALSE]
      rownames(rval2[[i]]) <- gsub("model.matrix.", "", rownames(rval2[[i]]), fixed = TRUE)
      if(!full.names) {
        rownames(rval2[[i]]) <- gsub(paste(i, "p.", sep = "."), "", rownames(rval2[[i]]), fixed = TRUE)
        rownames(rval2[[i]]) <- gsub(paste(i, "s.", sep = "."), "", rownames(rval2[[i]]), fixed = TRUE)
      }
      if(ncol(rval2[[i]]) < 2 & !summary) {
        rn <- rownames(rval2[[i]])
        rval2[[i]] <- rval2[[i]][, 1]
        names(rval2[[i]]) <- rn
      }
    }
    rval <- rval2
  } else {
    rownames(rval) <- gsub("model.matrix.", "", rownames(rval), fixed = TRUE)
    if(!full.names) {
      rownames(rval) <- gsub("p.", "", rownames(rval), fixed = TRUE)
      rownames(rval) <- gsub("s.", "", rownames(rval), fixed = TRUE)
      if(!is.null(model) & (length(model) < 2)) {
        for(i in nx)
          rownames(rval) <- gsub(paste(i, ".", sep = ""), "", rownames(rval), fixed = TRUE)
      }
    } 
    if(ncol(rval) < 2 & !summary) {
      rn <- rownames(rval)
      rval <- rval[, 1]
      names(rval) <- rn
    }
  }
  ## If data have been scaled (scale.d=TRUE)
  if ( ! is.null(attr(object$model.frame,'scale')) & rescale) {
    ## Get scaling
    sc <- attr(object$model.frame,'scale')
    for ( par in names(object$terms) ) {
      for ( nam in names(sc$scale) ) {
         # Descaling coefficients
         idx <- which(grepl(sprintf("%s.p.%s",par,nam),names(rval)))
         if ( length(idx) > 0 )
            rval[idx] <- rval[idx] / sc$scale[nam]
         # Descaling intercepts
         idx <- which(grepl(sprintf("%s.p.\\(Intercept\\)",par),names(rval)))
         if ( length(idx) > 0 & sprintf("%s.p.%s",par,nam) %in% names(rval) )
            rval[idx] <- rval[idx] - sc$center[nam] * rval[sprintf('%s.p.%s',par,nam)]
      }
    }
  }
  rval
}


## Get all terms names used.
term.labels <- function(x, model = NULL, pterms = TRUE, sterms = TRUE,
  intercept = TRUE, list = TRUE, ...)
{
  if(inherits(x, "bamlss") | inherits(x, "bamlss.frame")) {
    x <- terms(x)
  } else {
    if(!inherits(x, "bamlss.terms")) {
      if(inherits(x, "terms")) {
        x <- list("p" = x)
      } else stop("x must be a 'terms' or 'bamlss.terms' object!")
    }
  }

  nx <- names(x)

  if(is.null(model)) {
    model <- nx
  } else {
    if(!is.character(model)) {
      if(max(model) > length(nx) | min(model) < 1)
        stop("model is specified wrong")
    } else {
      for(j in seq_along(model)) {
        mm <- pmatch(model[j], nx)
        if(is.na(mm))
          stop("model is specified wrong")
        model[j] <- nx[mm]
      }
    }
  }

  x <- x[model]
  nx <- names(x)
  rval <- vector(mode = "list", length = length(nx))
  for(j in seq_along(x)) {
    rval[[j]] <- list()
    txj <- drop.terms.bamlss(x[[j]], pterms = pterms, sterms = sterms, keep.response = FALSE)
    tl <- attr(txj, "term.labels")
    specials <- unlist(attr(txj, "specials"))
    if(length(specials)) {
      sub <- if(attr(txj, "response") > 0) 1 else 0
      rval[[j]]$p <- tl[-1 * c(specials - sub)]
      rval[[j]]$s <- tl[specials - sub]
    } else {
      rval[[j]]$p <- tl
    }
    if(intercept & (attr(x[[j]], "intercept") > 0)) {
      rval[[j]]$p <- if(!length(tl)) "(Intercept)" else c("(Intercept)", rval[[j]]$p)
    }
  }
  names(rval) <- nx
  if(!list)
    rval <- unlist(rval)

  rval
}

term.labels2 <- function(x, model = NULL, pterms = TRUE, sterms = TRUE,
  intercept = TRUE, list = TRUE, type = 1, rm.by = TRUE, ...)
{
  stl <- NULL
  is.bamlss <- FALSE
  if(inherits(x, "bamlss") | inherits(x, "bamlss.frame")) {
    stl <- lapply(x$x, function(x) {
      if(!is.null(x$smooth.construct)) {
        nst <- names(x$smooth.construct)
        nst <- nst[nst != "model.matrix"]
        if(length(nst)) return(nst) else return(NULL)
      } else return(NULL)
    })
    is.bamlss <- TRUE
    x <- terms(x, drop = FALSE)
  } else {
    if(!inherits(x, "bamlss.terms")) {
      if(inherits(x, "terms")) {
        x <- list("p" = x)
      } else stop("x must be a 'terms' or 'bamlss.terms' object!")
    }
  }

  nx <- names(x)
  if(is.null(model)) {
    model <- nx
  } else {
    if(!is.character(model)) {
      if(max(model) > length(nx) | min(model) < 1)
        stop("model is specified wrong")
    } else {
      for(j in seq_along(model)) {
        mm <- pmatch(model[j], nx)
        if(is.na(mm))
          stop("model is specified wrong")
        model[j] <- nx[mm]
      }
    }
  }

  x <- x[model]
  if(!is.null(stl))
    stl <- stl[model]
  nx <- names(x)
  rval <- vector(mode = "list", length = length(nx))
  for(j in seq_along(x)) {
    txj <- drop.terms.bamlss(x[[j]], pterms = TRUE, sterms = !is.bamlss, keep.response = FALSE)
    if(type < 2) {
      rval[[j]] <- attr(txj, "term.labels")
    } else {
      rval[[j]] <- all.labels.formula(txj)
    }
    if(intercept & (attr(txj, "intercept") > 0))
      rval[[j]] <- c(rval[[j]], "(Intercept)")
    if(is.bamlss)
      rval[[j]] <- c(rval[[j]], stl[[j]])
  }
  names(rval) <- nx
  if(rm.by) {
    for(j in seq_along(rval)) {
      if(any(by <- (grepl("by=", rval[[j]], fixed = TRUE) & grepl("):", rval[[j]], fixed = TRUE)))) {
        for(i in which(by)) {
          rval[[j]][i] <- paste(strsplit(rval[[j]][i], "):", fixed = TRUE)[[1]][1], ")", sep = "")
        }
      }
      rval[[j]] <- unique(rval[[j]])
      rval[[j]] <- rval[[j]][rval[[j]] != ""]
    }
  }
  if(!list) {
    rval2 <- NULL
    for(j in seq_along(nx)) {
      names(rval[[j]]) <- rep(nx[j], length = length(rval[[j]]))
      rval2 <- c(rval2, rval[[j]])
    }
    rval <- rval2
  }
  rval
}


## Scores for model comparison.
score <- function(x, limits = NULL, FUN = function(x) { mean(x, na.rm = TRUE) },
  type = c("mean", "samples"), kfitted = TRUE, nsamps = NULL, ...)
{
  stopifnot(inherits(x, "bamlss"))
  family <- attr(x, "family")
  stopifnot(!is.null(family$d))
  type <- match.arg(type)
  y <- model.response2(x)
  n <- if(is.null(dim(y))) length(y) else nrow(y)
  maxy <- max(y, na.rm = TRUE)

  if(is.null(family$nscore)) {
    nscore <- function(eta) {
      integrand <- function(x) {
        int <- family$d(x, family$map2par(eta))^2
        int[int == Inf | int == -Inf] <- 0
        int
      }
      rval <- if(is.null(limits)) {
          try(integrate(integrand, lower = -Inf, upper = Inf), silent = TRUE)
        } else try(integrate(integrand, lower = limits[1], upper = limits[2]), silent = TRUE)
      if(inherits(rval, "try-error")) {
        rval <- try(integrate(integrand, lower = min(y, na.rm = TRUE),
          upper = max(y, na.rm = TRUE)))
      }
      rval <- if(inherits(rval, "try-error")) NA else rval$value
      rval
    }
  } else {
    nscore <- function(eta) {
	    integrand <- function(x) {
        family$d(x, family$map2par(eta))^2
	    }
	    rval <- sum(integrand(seq(0, maxy)))
	    rval
    }

	  nscore2 <- function(y, eta) {
	    integrand <- function(x) {
         -sum(((x == y) * 1 - family$d(x, family$map2par(eta)))^2)
	    }
	    rval <- (integrand(seq(0, maxy)))
	    rval
    }
  }

  scorefun <- function(eta) {
    norm <- rep(0, n)
    for(i in 1:n) {
      ni <- try(nscore(eta[i, , drop = FALSE]), silent = TRUE)
      if(inherits(ni, "try-error")) ni <- NA
      norm[i] <- ni
    }
    pp <- family$d(y, family$map2par(eta))
    pp[pp == Inf | pp == -Inf] <- 0
    loglik <- log(pp)
    if(is.null(family$nscore)) {
      quadratic <- 2 * pp - norm
    } else {
      quadratic <- rep(0, n)
      for(i in 1:n) {
        ni <- try(nscore2(y[i], eta[i, , drop = FALSE]), silent = TRUE)
        if(inherits(ni, "try-error")) ni <- NA
        quadratic[i] <- ni
      }
    }
    spherical <- pp / sqrt(norm)

    return(data.frame(
      "log" = FUN(loglik),
      "quadratic" = FUN(quadratic),
      "spherical" = FUN(spherical)
    ))
  }

  if(type == "mean") {
    eta <- if(kfitted) {
      kfitted(x, nsamps = nsamps,
        FUN = function(x) { mean(x, na.rm = TRUE) }, ...)
    } else fitted(x, samples = if(!is.null(h_response(x))) TRUE else FALSE)
    if(!inherits(eta, "list")) {
      eta <- list(eta)
      names(eta) <- family$names[1]
    }
    eta <- as.data.frame(eta)
    res <- unlist(scorefun(eta))
  } else {
    nx <- names(x)
    eta <- if(kfitted) {
      kfitted(x, FUN = function(x) { x }, nsamps = nsamps, ...)
    } else fitted(x, samples = TRUE, FUN = function(x) { x }, nsamps = nsamps)
    if(!inherits(eta, "list")) {
      eta <- list(eta)
      names(eta) <- family$names[1]
    }
    for(j in nx) {
      colnames(eta[[j]]) <- paste("i",
        formatC(1:ncol(eta[[j]]), width = nchar(ncol(eta[[1]])), flag = "0"),
        sep = ".")
    }
    nc <- ncol(eta[[1]])
    eta <- as.data.frame(eta)
    res <- list()
    for(i in 1:nc) {
      eta2 <- eta[, grep(ni <- paste(".i.",
        formatC(i, width = nchar(nc), flag = "0"), sep = ""),
        names(eta)), drop = FALSE]
      names(eta2) <- gsub(ni, "", names(eta2))
      res[[i]] <- scorefun(eta2)
    }
    res <- do.call("rbind", res)
  }

  res
}


## Compute fitted values with dropping data.
kfitted <- function(x, k = 5, weighted = FALSE, random = FALSE,
  engine = NULL, verbose = TRUE, FUN = mean, nsamps = NULL, ...)
{
  if(!inherits(x, "bamlss")) stop('argument x is not a "bamlss" object!')
  if(is.null(engine))
    engine <- attr(x, "engine")
  if(is.null(engine)) stop("please choose an engine!")
  mf <- model.frame(x)
  i <- rep(1:k, length.out = nrow(mf))
  if(random)
    i <- sample(i)
  k <- sort(unique(i))
  f <- formula(x)
  family <- family(x)
  ny <- length(unique(attr(mf, "response.name")))
  rval <- NULL
  jj <- 1
  for(j in k) {
    if(verbose) cat("subset:", jj, "\n")
    drop <- mf[i == j, ]
    if(!weighted) {
      take <- mf[i != j, ]
      bcv <- bamlss(f, data = take, family = family,
        engine = engine, verbose = verbose, ...)
    } else {
      w <- 1 * (i != j)
      bcv <- bamlss(f, data = mf, family = family,
        engine = engine, verbose = verbose, weights = w, ...)
    }
    if(!is.null(attr(mf, "orig.names")))
      names(drop) <- rmf(names(drop))
    fit <- fitted.bamlss(bcv, newdata = drop, samples = TRUE, FUN = FUN, nsamps = nsamps)
    if(!inherits(fit, "list")) {
      fit <- list(fit)
      names(fit) <- family$names
    }
    if(is.null(rval)) {
      rval <- list()
      for(ii in names(fit)) {
        rval[[ii]] <- matrix(NA, nrow = nrow(mf),
          ncol = if(is.null(dim(fit[[ii]]))) 1 else ncol(fit[[ii]]))
      }
    }
    for(ii in names(fit)) {
      rval[[ii]][i == j, ] <- fit[[ii]]
    }
    jj <- jj + 1
  }

  for(ii in names(fit)) {
    rval[[ii]] <- if(ncol(rval[[ii]]) > 1) {
      as.data.frame(rval[[ii]])
    } else drop(rval[[ii]])
  }

  if(length(rval) < 2)
    rval <- rval[[1]]

  rval
}


## Modified p() and d() functions.
create.dp <- function(family)
{
  if(is.null(names(family$links)))
    names(family$links) <- family$names
  links <- list()
  for(j in names(family$links))
    links[[j]] <- make.link2(family$links[j])$linkinv
  d <- function(y, eta, ...) {
    for(j in names(eta))
      eta[[j]] <- links[[j]](eta[[j]])
    family$d(y, eta, ...)
  }
  p <- function(y, eta, ...) {
    for(j in names(eta))
      eta[[j]] <- links[[j]](eta[[j]])
    family$p(y, eta, ...)
  }
  return(list("d" = d, "p" = p))
}


## Extract model residuals.
residuals.bamlss <- function(object, type = c("quantile", "response"), nsamps = NULL, ...)
{
  family <- family(object)

  if(!is.null(family$residuals)) {
    res <- family$residuals(object, type = type, nsamps = nsamps, ...)
    if(length(class(res)) < 2) {
      if(inherits(res, "numeric"))
        class(res) <- c("bamlss.residuals", class(res))
    }
  } else {
    type <- match.arg(type)

    if(is.null(object$y))
      stop("response variable is missing, cannot compute residuals!")

    nobs <- nrow(object$y)
    y <- if(is.data.frame(object$y)) {
      if(ncol(object$y) < 2) {
        object$y[[1]]
      } else object$y
    } else {
      object$y
    }

    par <- predict(object, nsamps = nsamps, drop = FALSE, ...)
    for(j in family$names)
      par[[j]] <- make.link2(family$links[j])$linkinv(par[[j]])

    if(type == "quantile") {
      if(is.null(family$p)) {
        type <- "response"
        warning(paste("no $p() function in family '", family$family,
          "', cannot compute quantile residuals, computing response resdiuals instead!", sep = ""))
      } else {
        res <- qnorm(family$p(y, par))
        attr(res, "type") <- "Quantile"
      }
    }

    if(type == "response") {
      mu <- if(is.null(family$mu)) {
        function(par, ...) { par[[1]] }
      } else family$mu
      res <- y - mu(par)
      attr(res, "type") <- "Response"
    }
   
    class(res) <- c("bamlss.residuals", class(res))
  }

  return(res)
}


## Residuals plotting functions.
plot.bamlss.residuals <- function(x, which = c("hist-resid", "qq-resid"), spar = TRUE, ...)
{
  ## What should be plotted?
  which.match <- c("hist-resid", "qq-resid")
  if(!is.character(which)) {
    if(any(which > 2L))
      which <- which[which <= 2L]
    which <- which.match[which]
  } else which <- which.match[pmatch(tolower(which), which.match)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")

  if(is.null(dim(x)))
    x <- matrix(x, ncol = 1)
  nc <- ncol(x)
  cn <- colnames(x)

  if(spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = n2mfrow(length(which) * nc))
  }

  type <- attr(x, "type")
  x <- x[is.finite(x)]

  for(j in 1:nc) {
    for(w in which) {
      args <- list(...)
      if(w == "hist-resid") {
        rdens <- density(x)
        rh <- hist(x, plot = FALSE)
        args$ylim <- c(0, max(c(rh$density, rdens$y)))
        args$freq <- FALSE
        args$x <- x
        args <- delete.args("hist.default", args, package = "graphics")
        if(is.null(args$xlab))
          args$xlab <- if(is.null(type)) "Residuals" else paste(type, "residuals")
        if(is.null(args$ylab))
          args$ylab <- "Density"
        if(is.null(args$main)) 
          args$main <- paste("Histogramm and density", if(!is.null(cn[j])) paste(":", cn[j]) else NULL)
        ok <- try(do.call("hist", args))
        if(!inherits(ok, "try-error"))
          lines(rdens)
        box()
      }
      if(w == "qq-resid") {
        args$y <- x
        args$x <- NULL
        args <- delete.args("qqnorm.default", args, package = "stats", not = c("col", "pch"))
        if(is.null(args$main))
          args$main <- paste("Normal Q-Q Plot", if(!is.null(cn[j])) paste(":", cn[j]) else NULL)
        ok <- try(do.call(qqnorm, args))
        if(!inherits(ok, "try-error"))
          qqline(args$y) ## abline(0,1)
      }
    }
  }

  return(invisible(NULL))
}


## Extract the model response.
model.response2 <- function(data, hierarchical = FALSE, ...)
{
  if(!inherits(data, "data.frame")) {
    f <- if(inherits(data, "bamlss")) formula(data) else NULL
    data <- model.frame(data)
    if(!is.null(f)) {
      if("h1" %in% names(f)) {
        rn <- all.vars(f$h1)[1]
        attr(data, "response.name") <- rn
      } else {
        rn <- NULL
        for(j in seq_along(f)) {
          if(is.list(f[[j]])) {
            if("h1" %in% names(f[[j]]))
              rn <- c(rn, all.vars(f[[j]]$h1)[1])
          }
        }
        rn <- rn[rn %in% names(data)]
        if(length(rn))
          attr(data, "response.name") <- rn
      }
    }
  }
  rn <- attr(data, "response.name")
  y <- if(is.null(rn)) {
    model.response(data, ...)
  } else data[, unique(rn), ...]
  y
}

## find hierarchical responses
h_response <- function(x)
{
  rval <- NULL
  if(!all(c("model", "fitted.values") %in% names(x))) {
    for(j in seq_along(x))
      rval <- c(rval, h_response(x[[j]]))
  } else {
    if(!is.null(x$model$hlevel)) {
      if(x$model$hlevel > 1)
        rval <- response.name(x$model$formula)
    }
  }
  rval
}


blockMatrixDiagonal<-function(...){  
  matrixList<-list(...)
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
 
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  finalMatrix
}

## Create the inverse of a matrix.
matrix_inv <- function(x, index = NULL, force = FALSE)
{
  if(!is.null(index$block.index)) {
    if(!is.matrix(x)) x <- as.matrix(x)
    return(.Call("block_inverse", x, index$block.index, index$is.diagonal))
  }
  if(inherits(x, "Matrix")) {
    if(!is.null(index$crossprod)) {
      if(ncol(index$crossprod) < 2) {
        return(Diagonal(x = 1 / diag(x)))
      } else {
        return(chol2inv(chol(x)))
      }
    }
  }
  if(!is.null(index$crossprod)) {
    if(ncol(index$crossprod) < ncol(x)) {
      if(ncol(index$crossprod) < 2) {
        return(diag(1 / diag(x)))
      } else {
        if(FALSE) {
          ju <- unique(index$crossprod[, 1])
          if(length(ju) < nrow(x)) {
            inv <- list()
            for(i in ju) {
              take <- index$crossprod[, 1] == i
              inv[[as.character(i)]] <- solve(x[take, take, drop = FALSE])
            }
            return(as.matrix(do.call("bdiag", inv)))
          }
        }
      }
    }
  }
  if(length(x) < 2)
    return(1 / x)
  rn <- rownames(x)
  cn <- colnames(x)
  if(!is.null(index) & is(x, "spam")) {
    p <- update.spam.chol.NgPeyton(index$spam.cholFactor, x)
  } else {
    p <- try(chol(x), silent = TRUE)
  }
  p <- if(inherits(p, "try-error")) {
    try(solve(x), silent = TRUE)
  } else {
    try(if(inherits(x, "spam")) chol2inv.spam(p) else chol2inv(p), silent = TRUE)
  }
  if(inherits(p, "try-error")) {
    diag(x) <- jitter(diag(x), amount = 1e-5)
    p <- try(solve(x), silent = TRUE)
  }
  if(inherits(p, "try-error") & force) {
    p <- diag(ncol(x))
  }
  rownames(p) <- rn
  colnames(p) <- cn
  return(p)
}


## Compute matching index for duplicates in data.
match.index <- function(x)
{
  if(!is.vector(x)) {
    if(!inherits(x, "matrix") & !inherits(x, "data.frame"))
      stop("x must be a matrix or a data.frame!")
    x <- if(inherits(x, "matrix")) {
      apply(x, 1, paste, sep = "\r", collapse = ";")
    } else do.call("paste", c(x, sep = "\r"))
  }
  nodups <- which(!duplicated(x))
  ind <- match(x, x[nodups])
  return(list("match.index" = ind, "nodups" = nodups))
}


XinY <-
    function(x, y, by = intersect(names(x), names(y)), by.x = by, by.y = by,
             notin = FALSE, incomparables = NULL,
             ...)
{
    fix.by <- function(by, df)
    {
        ## fix up 'by' to be a valid set of cols by number: 0 is row.names
        if(is.null(by)) by <- numeric(0L)
        by <- as.vector(by)
        nc <- ncol(df)
        if(is.character(by))
            by <- match(by, c("row.names", names(df))) - 1L
        else if(is.numeric(by)) {
            if(any(by < 0L) || any(by > nc))
                stop("'by' must match numbers of columns")
        } else if(is.logical(by)) {
            if(length(by) != nc) stop("'by' must match number of columns")
            by <- seq_along(by)[by]
        } else stop("'by' must specify column(s) as numbers, names or logical")
        if(any(is.na(by))) stop("'by' must specify valid column(s)")
        unique(by)
    }

    nx <- nrow(x <- as.data.frame(x)); ny <- nrow(y <- as.data.frame(y))
    by.x <- fix.by(by.x, x)
    by.y <- fix.by(by.y, y)
    if((l.b <- length(by.x)) != length(by.y))
        stop("'by.x' and 'by.y' specify different numbers of columns")
    if(l.b == 0L) {
        ## was: stop("no columns to match on")
        ## returns x
        return(x)
    }
    else {
        if(any(by.x == 0L)) {
            x <- cbind(Row.names = I(row.names(x)), x)
            by.x <- by.x + 1L
        }
        if(any(by.y == 0L)) {
            y <- cbind(Row.names = I(row.names(y)), y)
            by.y <- by.y + 1L
        }
        ## create keys from 'by' columns:
        if(l.b == 1L) {                  # (be faster)
            bx <- x[, by.x]; if(is.factor(bx)) bx <- as.character(bx)
            by <- y[, by.y]; if(is.factor(by)) by <- as.character(by)
        } else {
            ## Do these together for consistency in as.character.
            ## Use same set of names.
            bx <- x[, by.x, drop=FALSE]; by <- y[, by.y, drop=FALSE]
            names(bx) <- names(by) <- paste("V", seq_len(ncol(bx)), sep="")
            bz <- do.call("paste", c(rbind(bx, by), sep = "\r"))
            bx <- bz[seq_len(nx)]
            by <- bz[nx + seq_len(ny)]
        }
        comm <- match(bx, by, 0L)
        if (notin) {
            res <- x[comm == 0,]
        } else {
            res <- x[comm > 0,]
        }
    }
    ## avoid a copy
    ## row.names(res) <- NULL
    attr(res, "row.names") <- .set_row_names(nrow(res))
    res
}


XnotinY <-
    function(x, y, by = intersect(names(x), names(y)), by.x = by, by.y = by,
             notin = TRUE, incomparables = NULL,
             ...)
{
    XinY(x,y,by,by.x,by.y,notin,incomparables)
}


## Small helper function to scale the model.matrix.
scale.model.matrix <- function(x)
{
  if(!is.matrix(x))
    x <- as.matrix(x)
  cn <- colnames(x)
  center <- as.numeric(colMeans(x, na.rm = TRUE))
  scale <- as.numeric(apply(x, 2, sd, na.rm = TRUE))
  if(length(i <- grep("(Intercept)", cn, fixed = TRUE))) {
    center[i] <- 0.0
    scale[i] <- 1.0
  }
  x <- .Call("scale_matrix", x, center, scale, PACKAGE = "bamlss")
  attr(x, "scale") <- list("center" = center, "scale" = scale)
  x
}

## Small helper function to scale the model.matrix.
scale_model.frame <- function(x, not = "")
{
  if(!inherits(x, "data.frame"))
    x <- as.data.frame(x)
  cn <- colnames(x)
  cn2 <- scales <- centers <- NULL
  for(j in cn) {
    if(!is.factor(x[[j]]) & !(j %in% not)) {
      cx <- mean(as.numeric(x[[j]]), na.rm = TRUE)
      sx <- sd(as.numeric(x[[j]]), na.rm = TRUE)
      x[[j]] <- (x[[j]] - cx) / sx
      cn2 <- c(cn2, j)
      scales <- c(scales, sx)
      centers <- c(centers, cx)
    }
  }
  names(centers) <- cn2
  names(scales) <- cn2
  attr(x, "scale") <- list("center" = centers, "scale" = scales)
  x
}


## Sum of diagonal elements.
sum_diag <- function(x)
{
  if(inherits(x, "spam"))
    return(sum(diag.spam(x), na.rm = TRUE))
  if(inherits(x, "Matrix"))
    return(sum(diag(x), na.rm = TRUE))
  if(is.null(dx <- dim(x)))
    stop("x must be a matrix!")
  if(dx[1] != dx[2])
    stop("x must be symmetric!")
  .Call("sum_diag", x, dx[1], PACKAGE = "bamlss")
}

sum_diag2 <- function(x, y)
{
  if(inherits(x, "spam"))
    x <- as.matrix(x)
  if(inherits(y, "spam"))
    y <- as.matrix(y)
  if(is.null(dx <- dim(x)))
    stop("x must be a matrix!")
  if(is.null(dy <- dim(y)))
    stop("y must be a matrix!")
  if(dx[1] != dx[2])
    stop("x must be symmetric!")
  if(dy[1] != dy[2])
    stop("y must be symmetric!")
  .Call("sum_diag2", x, y, PACKAGE = "bamlss")
}

AIC.bamlss <- function(..., k = 2, optimizer = FALSE, samples = FALSE, FUN = mean)
{
  val <- logLik.bamlss(..., optimizer = optimizer, samples = samples)
  if(!is.list(val)) {
    val <- list(val)
  }
  val <- lapply(val, function(x) {
    if(!("edf" %in% colnames(x)))
      stop("cannot compute AIC, edf are missing!")
    x <- as.data.frame(cbind("AIC" = -2 * x[,"logLik"] + k * x[,"edf"], "edf" = x[,"edf"]))
    x
  })
  if(samples) {
    val <- lapply(val, function(x) {
      xn <- colnames(x)
      x <- apply(x, 2, FUN)
      names(x) <- xn
      x
    })
  }
  Call <- match.call()
  Call$k <- NULL
  Call$optimizer <- NULL
  Call$samples <- NULL
  if(!samples) {
    val <- val[[1]]
    row.names(val) <- if(nrow(val) > 1) Call[-1L] else ""
  } else {
    names(val) <- rep(as.character(Call[-1L]), length.out = length(val))
    val <- do.call("rbind", val)
  }
  val
}

BIC.bamlss <- function(..., k = 2, optimizer = FALSE, samples = FALSE, FUN = mean)
{
  nobs <- sapply(list(...), function(x) { if(is.null(nrow(x$y[[1]]))) nrow(model.frame(x)) else nrow(x$y[[1]]) })
  val <- logLik.bamlss(..., optimizer = optimizer, samples = samples)
  if(!is.list(val)) {
    val <- list(val)
  }
  val <- lapply(1:length(val), function(i) {
    if(!("edf" %in% colnames(val[[i]])))
      stop("cannot compute BIC, edf are missing!")
    x <- as.data.frame(cbind("BIC" = -2 * val[[i]][,"logLik"] + val[[i]][,"edf"] * log(nobs[i]), "edf" = val[[i]][,"edf"]))
    x
  })
  if(samples) {
    val <- lapply(val, function(x) {
      xn <- colnames(x)
      x <- apply(x, 2, FUN)
      names(x) <- xn
      x
    })
  }
  Call <- match.call()
  Call$k <- NULL
  Call$optimizer <- NULL
  Call$samples <- NULL
  if(!samples) {
    val <- val[[1]]
    row.names(val) <- if(nrow(val) > 1) Call[-1L] else ""
  } else {
    names(val) <- rep(as.character(Call[-1L]), length.out = length(val))
    val <- do.call("rbind", val)
  }
  val
}

#.First.lib <- function(lib, pkg)
#{
#  library.dynam("bamlss", pkg, lib)
#}


## TS-Decomp.
stg <- function(x, interp = FALSE, k = -1, ...)
{
  if(interp) {
    x <- zoo::na.approx(x, rule = 2)
  }

  xf <- stats::frequency(x)
  xc <- stats::cycle(x)
  mf <- na.omit(data.frame("x" = as.numeric(x), "trend" = 1:NROW(x), "season" = as.integer(xc),
    "xlag" = c(NA, as.numeric(x)[-length(x)])))

  k <- rep(k, length.out = 3)
  if(k[1] < 0)
    k[1] <- floor(length(unique(mf$trend)) * 0.1)
  if(k[2] < 0)
    k[2] <- floor(length(unique(mf$season)) * 0.9)
  if(k[3] < 0)
    k[3] <- 10

  b <- gam(x ~ s(trend,by=xlag,k=k[1]) + ti(trend,k=k[1],bs="cr") + ti(season,k=k[2],bs="cc") +
    ti(trend,season,bs=c("cr","cc"),k=k[3]), data = mf, method = "REML")

  plot(b, pages = 1)

  p <- predict(b, type = "terms")

  rval <- cbind("raw" = mf[["x"]], "fitted" = fitted(b),
    "trend" = p[, "ti(trend)"] + coef(b)[1],
    "seasonal" = p[, "ti(season)"] + p[, "ti(trend,season)"],
    "lag1" = p[, "s(trend):xlag"] / mf$xlag,
    "remainder" = stats::residuals(b))

  rval <- stats::ts(rval, start = start(x), frequency = xf)

  rval
}


if(FALSE) {
  x <- runif(3000, -3, 3)
  f <- bamlss:::simfun("pick")
  y <- sin(x) + rnorm(3000, sd = scale2(f(scale2(x, 0, 1)), 0.01, 0.3))
  plot(x, y)
  b <- bamlss(list(y ~ n(x), ~ n(x)), sampler = FALSE)
}
