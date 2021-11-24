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
      if(!is.null(family$bayesx)) {
        family$bayesx <- rep(family$bayesx, length(family$names))
        names(family$bayesx) <- family$names
      }
      family$links <- rep(family$links, length.out = length(formula))
      names(family$links) <- names(formula)
      attr(formula, "orig.formula") <- orig.formula
    } else {
      rn <- response.name(formula, hierarchical = FALSE, keep.functions = TRUE)
      if(any(i <- grepl("|", rn, fixed = TRUE))) {
        rnt <- rn[i]
        rn <- rn[!i]
        rnt <- gsub(" ", "", unlist(strsplit(rnt, "|", fixed = TRUE)))
        rn <- c(rn, rnt)
      }
      rn <- rn[rn %in% names(bf$model.frame)]
      bf$y <- bf$model.frame[rn]
      if(is.null(family$nocat)) {
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
#      if(family$family == "dirichlet") {
#        names(formula) <- colnames(bf$y)
#        family$bayesx <- rep(list(c("dirichlet", "alpha")), length(formula))
#        family$names <- names(family$bayesx) <- names(formula)
#        family$links <- rep(family$links, length(formula))
#        names(family$links) <- family$names
#        attr(family$bayesx, "nrcat") <- length(formula)
#      }
    }
  } else {
    overwrite <- list(...)$overwrite
    ff_name <- list(...)$ff_name
    if(is.null(ff_name))
      ff_name <- "ff_data_bamlss"
    if(is.null(overwrite))
      overwrite <- TRUE
    if(file.exists(ff_name) & overwrite) {
      unlink(ff_name, recursive = TRUE, force = TRUE)
    }
    if(!file.exists(ff_name)) {
      cat(paste0("  .. creating directory '", ff_name, "' for storing matrices. Note, the directory is may not deleted and matrices can be used for another model. Use delete = TRUE in the bamlss call. Before starting a new model you can set overwrite = TRUE to overwrite existing data.\n"))
      dir.create(ff_name)
    }
    rn <- response.name(formula, hierarchical = FALSE, keep.functions = FALSE)
    if(!any(rn %in% names(bf$model.frame))) {
      rn <- grep2(paste0(rn, "."), names(bf$model.frame), fixed = TRUE, value = TRUE)
    } else {
      rn <- rn[rn %in% names(bf$model.frame)]
    }
    bf$y <- bf$model.frame[rn]
    bf$ff_name <- ff_name
  }

  bf$formula <- formula
  attr(bf$formula, "response.name") <- rn

  ## Add the terms object.
  bf$terms <- terms.bamlss.formula(formula, data = data, drop = FALSE, specials = specials, ...)

  ## Process possible score and hess functions.
  if(!is.null(score <- family$score)) {
    if(is.function(score)) {
      score <- list(score)
      family$score <- rep(score, length.out = length(formula))
      names(family$score) <- names(formula)
    }
  }
  if(!is.null(hess <- family$hess)) {
    if(is.function(hess)) {
      hess <- list(hess)
      family$hess <- rep(hess, length.out = length(formula))
      names(family$hess) <- names(formula)
    }
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

  ## Delete ff directory?
  bf$delete <- list(...)$delete
  if(is.null(bf$delete))
    bf$delete <- TRUE

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


bamlss_chunk <- function(x) {
  N <- if(is.null(dim(x))) {
    length(x)
  } else {
    nrow(x)
  }
  if(N > 100000L) {
    n <- max(c(100000, 0.02 * N))
    bn <- floor(N/n)
    if(bn < 2)
      bn <- 10
    xc <- cut(seq_len(N), breaks = bn, include.lowest = TRUE)
    chunks <- split(seq_len(N), xc)
  } else {
    chunks <- list(seq_len(N))
  }
  chunks
}


## Compute the 'bamlss.frame' 'x' master object.
design.construct <- function(formula, data = NULL, knots = NULL,
  model.matrix = TRUE, smooth.construct = TRUE, binning = FALSE,
  before = TRUE, gam.side = NULL, model = NULL, drop = NULL,
  scale.x = TRUE, absorb.cons = NULL, sparse.cons = 0, specials = NULL, ...)
{
  if(!model.matrix & !smooth.construct)
    return(NULL)

  doCmat <- list(...)$Cmat
  if(is.null(doCmat))
    doCmat <- FALSE

  if(is.null(gam.side))
    gam.side <- FALSE

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
  if(!no_ff) {
    stopifnot(requireNamespace("bit"))
    stopifnot(requireNamespace("ff"))
    stopifnot(requireNamespace("ffbase"))
  }

  ff_name <- list(...)$ff_name
  if(is.null(ff_name))
    ff_name <- "ff_data_bamlss"

  nthres <- list(...)$nthres
  if(is.null(nthres))
    nthres <- 1

  if(!is.character(data) & no_ff) {
    if(!inherits(data, "data.frame"))
      data <- as.data.frame(data)
  }
  if(is.character(data)) {
    data <- ff::read.table.ffdf(file = data,
      na.strings = "", header = TRUE, sep = ",")
  }
  if(inherits(data, "ffdf")) {
    before <- TRUE
    gam.side <- FALSE
    if(is.null(binning))
      binning <- FALSE
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
        drop_terms_attr <- function(x) {
          attr(x, "terms") <- NULL
          return(x)
        }
        obj$model.matrix <- try(model.matrix(drop.terms.bamlss(obj$terms,
          sterms = FALSE, keep.response = FALSE, data = data, specials = specials), data = drop_terms_attr(data)), silent = TRUE)
        if(inherits(obj$model.matrix, "try-error")) {
          lmt <- drop.terms.bamlss(obj$terms,
            sterms = FALSE, keep.response = FALSE, data = data, specials = specials)
          environment(lmt) <- .GlobalEnv
          obj$model.matrix <- model.matrix(lmt, data = drop_terms_attr(data))
        }
        if(ncol(obj$model.matrix) > 0) {
          if(scale.x)
            obj$model.matrix <- scale.model.matrix(obj$model.matrix)
        } else obj$model.matrix <- NULL
      } else {
        mm_terms <- drop.terms.bamlss(obj$terms,
          sterms = FALSE, keep.response = FALSE, data = NULL, specials = specials)
        obj$model.matrix <- NULL
        ff_mm <- function(x) {
          X <- model.matrix(mm_terms, data = x)
          cn <- colnames(X)
          X <- ff::ff(X, dim = dim(X), dimorder = c(2, 1))
          colnames(X) <- cn
          return(X)
        }
        mm_test <- model.matrix(mm_terms, data = data[1:10, , drop = FALSE])
        if(ncol(mm_test) > 0) {
          nobs <- nrow(data)
          obj$model.matrix <- ff::ff(0.0,
            length = nobs * ncol(mm_test),
            dim = c(nobs, ncol(mm_test)))
          k <- 1
          np <- 0
          for(ic in bamlss_chunk(data)) {
            obj$model.matrix[ic, ] <- model.matrix(mm_terms, data = data[ic, ])
            np <- np + length(ic)
            if(k > 1)
              cat("\r")
            cat("  .. ..", paste0(formatC(np / nobs * 100, width = 7), "%"))
            k <- k + 1
          }
          colnames(obj$model.matrix) <- colnames(mm_test)
        }
      }
    }
    if(smooth.construct) {
      tx <- drop.terms.bamlss(obj$terms,
        pterms = FALSE, keep.response = FALSE, data = data, specials = specials)
      sid <- unlist(attr(tx, "specials"))
      smt <- NULL
      fterms <- NULL
      if(any(fj <- names(sid) %in% c("lf", "af",  "lf.vd", "re", "peer", "fpc"))) {
        fterms <- attr(tx, "term.labels")[sid[fj]]
        sid <- sid[!fj]
      }
      if(!length(sid))
        sid <- NULL
      if(!is.null(sid) | !is.null(fterms)) {
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
            if(inherits(tsm, "tensor.smooth.spec")) {
              if(!is.null(tsm$margin)) {
                tsm$xt <- tsm$margin[[1]]$xt
                if(is.list(tsm$xt[[1]]))
                  tsm$xt <- tsm$xt[[1]]
              }
            }
            if(is.null(tsm$xt))
              tsm$xt <- list()
            if(is.null(tsm$xt$binning))
              tsm$xt$binning <- binning
            acons <- TRUE
            if(!is.null(tsm$xt$center))
              acons <- tsm$xt$center
            tsm$xt$center <- acons
            tsm$xt$before <- before
            if(!is.null(tsm$xt$binning)) {
              term.names <- c(tsm$term, if(tsm$by != "NA") tsm$by else NULL)
              if(!inherits(data, "ffdf")) {
                tsm$binning <- match.index(data[, term.names, drop = FALSE])
                tsm$binning$order <- order(tsm$binning$match.index)
                tsm$binning$sorted.index <- tsm$binning$match.index[tsm$binning$order]
              } else {
                stop('binning not allowd using data of class "ffdf"!')
              }
              if(!inherits(data, "ffdf")) {
                if(!doCmat) {
                  smt <- smoothCon(tsm, if(before) data[tsm$binning$nodups, term.names, drop = FALSE] else data,
                    knots, absorb.cons = if(is.null(absorb.cons)) acons else absorb.cons, sparse.cons = sparse.cons, scale.penalty=TRUE)
                } else {
                  smt <- smooth.construct(tsm, if(before) data[tsm$binning$nodups, term.names, drop = FALSE] else data, knots)
                  smt$C <- Cmat(smt)
                  smt$doCmat <- TRUE
                  smt <- list(smt)
                }
                smooth <- c(smooth, smt)
              } else {
                smt <- smooth.construct_ff(tsm, data,
                  knots, absorb.cons = if(is.null(absorb.cons)) acons else absorb.cons,
                  ff_name = ff_name, nthres = nthres)
                smooth <- c(smooth, list(smt))
              }
            } else {
              if(inherits(data, "ffdf")) {
                smt <- smooth.construct_ff(tsm, data,
                  knots, absorb.cons = if(is.null(absorb.cons)) acons else absorb.cons,
                  ff_name = ff_name, nthres = nthres)
                smt <- list(smt)
              } else {
                smt <- smoothCon(tsm, data, knots,
                  absorb.cons = if(is.null(absorb.cons)) acons else absorb.cons,
                  sparse.cons = sparse.cons, scale.penalty=TRUE)
              }
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
              if(is.null(tsm$xt$nrep)) {
                if(!inherits(data, "ffdf")) {
                  smt2 <- smooth.construct(tsm, data, knots)
                } else {
                  smt2 <- smooth.construct_ff(tsm, data, knots, ff_name = ff_name, nthres = nthres)
                }
              } else {
                smt2 <- list()
                for(jnr in 1:tsm$xt$nrep) {
                  if(!inherits(data, "ffdf")) {
                    smt2[[jnr]] <- smooth.construct(tsm, data, knots)
                  } else {
                    smt2[[jnr]] <- smooth.construct_ff(tsm, data, knots, ff_name = ff_name, nthres = nthres)
                  }
                }
                class(smt2) <- c("no.mgcv", "smooth.list")
              }
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
        if(!is.null(fterms)) {
          for(j in seq_along(fterms)) {
            nenv <- new.env()
            for(nd in names(data))
              assign(nd, data[[nd]], envir = nenv)
            vars <- all.vars(as.formula(paste("~", fterms[j])))
            if(length(vars) > 1) {
              for(vj in vars[-1]) {
                tmp <- get(vj, envir = .GlobalEnv)
                assign(vj, tmp, envir = nenv)
              }
            }
            tfm <- eval(parse(text = fterms[j]), envir = nenv)
            rm(nenv)
            tfme <- eval(tfm$call, envir = tfm$data)
            smt <- smoothCon(tfme, data = tfm$data, n = nrow(tfm$data[[1L]]),
              knots = knots, absorb.cons = TRUE,scale.penalty=TRUE)
            lab <- all.labels.formula(as.formula(paste("~", fterms[j])))
            for(jj in seq_along(smt)) {
              smt[[jj]]$model.frame <- tfm$data
              smt[[jj]]$orig.label <- smt[[jj]]$label
              smt[[jj]]$label <- lab
              smt[[jj]]$is.refund <- TRUE
              smt[[jj]]$refund.call <- fterms[j]
            }
            smooth <- c(smooth, smt)
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
        k <- 1
        for(j in seq_along(sld)) {
          for(jj in which(sl == sld[j])) {
            clj <- class(obj$smooth.construct[[jj]])
            clj <- strsplit(clj, ".", fixed = TRUE)[[1]][1]
            if(clj == "random")
              clj <- "re"
            if(clj == "nnet")
              clj <- ""
            sl[jj] <- paste(sl[jj], clj, sep = ":")
            sl[jj] <- gsub(paste0("):", clj), paste0(",id='", clj, k, "')"), sl[jj], fixed = TRUE)
            obj$smooth.construct[[jj]]$label <- sl[jj]
            k <- k + 1
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
                formula[[i]][[j]]$smooth.construct[[k]]$hess <- priors$hess
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


## Package ff smooth constructors.
smooth.construct_ff <- function(object, data, knots, ...)
{
  UseMethod("smooth.construct_ff")
}


ff_nrow <- function(x, value)
{
  d <- dim(x)
  if(is.null(d) || length(d)!=2)
    stop("not a two-dimensional array")
  dim(x) <- c(as.integer(value), d[[2]])
  x
}

ff_ncol <- function(x, value)
{
  d <- dim(x)
  if (is.null(d) || length(d)!=2)
    stop("not a two-dimensional array")
  dim(x) <- c(d[[1]], as.integer(value))
  x
}


#ff_matrix_append <- function(x, dat, recode = TRUE, adjustvmode = TRUE, ...) 
#{
#  stopifnot(requireNamespace("bit"))
#  stopifnot(requireNamespace("ff"))
#  stopifnot(requireNamespace("ffbase"))

#  w <- getOption("warn")
#  options("warn" = -1)
#  if(is.null(x))
#    return(dat)
#  n <- nrow(dat)
#  nff <- nrow(x)
#  cn <- colnames(x)
#  ## x <- ff_nrow(x, nff + n)
#  nrow(x) <- nff + n
#  if(!identical(colnames(x), colnames(dat))) {
#    warning("column names are not identical")
#  }
#  if(ncol(x) != ncol(dat)) {
#    stop("Number of columns does not match")
#  }
#  i <- hi(nff + 1, nff + n)
#  colnames(x) <- NULL
#  colnames(dat) <- NULL
#  x[i, ] <- dat[,]
#  colnames(x) <- cn
#  options("warn" = w)
#  x
#}

#ffdf_2_ff_matrix <- function(x, ...) 
#{
#  result <- ff::ff(NA, dim = dim(x), vmode = names(maxffmode(vmode(x)))[1], ...)
#  dimnames(result) <- dimnames(x)
#  for(i in chunk(x)) {
#    Log$chunk(i)
#    result[i, ] <- as.matrix(x[i, ])
#  }
#  result
#}

smooth.construct_ff.default <- function(object, data, knots, ff_name, nthres = NULL, ...)
{
  object$xt$center <- TRUE
  object$xt$nocenter <- FALSE
  terms <- object$term
  if(object$by != "NA") {
    object$label <- strsplit(object$label, "")[[1]]
    object$label <- paste0(object$label[-length(object$label)], collapse = "")
    object$label <- paste0(object$label, ",by=", object$by, ")")
    object$xt$center <- FALSE
    object$xt$nocenter <- TRUE
    terms <- c(terms, object$by)
  }
  nd <- list()
  cat("  .. ff processing term", object$label, "\n")
  xfile <- rmf(object$label)
  xfile <- file.path(ff_name, xfile)
  is_f <- sapply(data, is.factor)
  if(is.null(nthres))
    nthres <- 1
  if(nrow(data) > nthres) {
    if((length(terms) > 1) & !any(is_f) & FALSE) {
      ud <- nrow(unique(data[, terms]))
      km <- kmeans(data[, terms], min(c(1000, floor(0.9 * ud))))
      uc <- unique(km$cluster)
      nd <- matrix(NA, length(uc), length(terms))
      for(i in seq_along(uc)) {
        nd[i, ] <- as.numeric(data[sample(which(km$cluster == uc[i]), size = 1), terms])
      }
      nd <- as.data.frame(nd)
      names(nd) <- terms
    ##nd <- data[sample(1:nrow(data), size = 1000L), ]
    } else {
      for(j in terms) {
        if(!is.factor(data[[j]][1:2])) {
          ux <- ffbase::unique.ff(data[[j]])
          uxn <- length(ux)
          if(uxn > 2) {
            uxl <- if(uxn < 1000L) uxn - 1L else 1000L
            xq <- ffbase::quantile.ff(data[[j]], probs = seq(0, 1, length = uxl), na.rm = TRUE)
            names(xq) <- NULL
            if(length(unique(xq)) < 100) {
              xq <- rep(ux[], length.out = 1000L)
            }
          } else {
            xq <- rep(ux[], length.out = 1000L)
          }
          if(length(xq) == 1000L) {
            nd[[j]] <- sample(xq)
          } else {
            nd[[j]] <- sample(rep(xq, length.out = 1000L))
          }
        } else {
          nd[[j]] <- sample(rep(unique(data[[j]]), length.out = 1000L))
        }
      }
    }
    nd <- as.data.frame(nd)
  }
  object <- smoothCon(object, data = if(nrow(data) > nthres) nd else as.data.frame(data),
    knots = knots, absorb.cons = TRUE)[[1L]]
  rm(nd)
  nobs <- nrow(data)
  if(file.exists(paste0(xfile, ".rds"))) {
    object[["X"]] <- readRDS(paste0(xfile, ".rds"))
    bit::physical(object[["X"]])$filename <- paste0(xfile, ".ff")
  } else {
    object[["X"]] <- ff::ff(0.0,
      length = nrow(data) * ncol(object[["X"]]),
      dim = c(nrow(data), ncol(object[["X"]])),
      ## vmode = names(maxffmode(vmode(object[["X"]])))[1],
      filename = paste0(xfile, ".ff"))
    sX <- function(x) {
      if(is.null(object$PredictMat)) {
        X <- PredictMat(object, data = x)
      } else {
        X <- object$PredictMat(object, data = x)
      }
      return(X)
    }
    k <- 1
    np <- 0
    for(ic in bamlss_chunk(data)) {
      object[["X"]][ic, ] <- sX(data[ic, ])
      np <- np + length(ic)
      if(k > 1)
        cat("\r")
      cat("  .. ..", paste0(formatC(np / nobs * 100, width = 7), "%"))
      k <- k + 1
    }
    saveRDS(object[["X"]], file = paste0(xfile, ".rds"))
    cat("\n")
  }
  if(!inherits(object, "nnet0.smooth") & FALSE) {
    csum <- 0
    for(ic in bamlss_chunk(object[["X"]])) {
      csum <- csum + colSums(object[["X"]][ic, ])
    }
    QR <- qr(matrix(csum, ncol = 1L))
    object[["Z"]] <- qr.Q(QR, complete = TRUE)[, -1]
    tX <- try(ffmatrixmult(object[["X"]], object[["Z"]]), silent = TRUE)
    if(!inherits(tX, "try-error")) {
      object[["X"]] <- tX
      for(j in seq_along(object[["S"]])) {
        if(!is.function(object[["S"]][[j]])) {
          object[["S"]][[j]] <- crossprod(object[["Z"]], object[["S"]][[j]]) %*% object[["Z"]]
        }
      }
    } else {
      stop(paste("could not process term", object$label))
    }
  }
  object$orig.class <- class(object)
  class(object) <- "ff_smooth.smooth.spec"
  return(object)
}

Predict.matrix.ff_smooth.smooth.spec <- function(object, data)
{
  class(object) <- object$orig.class
  if(is.null(object$PredictMat)) {
    X <- Predict.matrix(object, data)
    if(!is.null(object[["Z"]]))
      X <- X %*% object[["Z"]]
  } else {
    X <- object$PredictMat(object, data)
  }
  return(X)
}

## Copy from bootSVD.
ffmatrixmult <- function(x,y=NULL,xt=FALSE,yt=FALSE,ram.output=FALSE, override.big.error=FALSE,...) {	
	{i1<-NULL; i2<- NULL} #To avoid errors in R CMD check

  stopifnot(requireNamespace("ff"))
  stopifnot(requireNamespace("ffbase"))

	dimx<-dim(x)
	if(!is.null(y)) dimy<-dim(y)
	if(is.null(y))  dimy<-dimx

	p <- max(c(dimx,dimy))
	n <- max(min(dimx),min(dimy))
	outDim <-
	inDim  <- rep(NA,2)
	outDim[1] <- dimx[xt+1]
	outDim[2] <- dimy[2-yt]
	inDim[1] <- dimx[2-xt]
	inDim[2] <- dimy[yt+1]
	if(inDim[1]!=inDim[2]) stop('non-conformable arguments')

	if(all(outDim>n) & (!override.big.error)) stop('Returned value is at risk of being extremely large. Both dimensions of output will be fairly large.')

	if(xt & yt) stop('For ff matrix algebra, set only one of xt or yt to TRUE')


	if(all(outDim==n) | (!'ff'%in% c(class(x),class(y)))|ram.output){ 
		out <- matrix(0,outDim[1],outDim[2])
	}else{
		out <- ff::ff(0,dim=outDim,...)
	}

	if(all(outDim==n)){
		if( (xt) &(!yt)) ff::ffapply({
			out<-out+crossprod(x[i1:i2,], y[i1:i2,])
			},X=x,MARGIN=1)
		if((!xt) & (yt)) ff::ffapply({
			out<-out+tcrossprod(x[,i1:i2], y[,i1:i2])
			},X=x,MARGIN=2)
		if((!xt) &(!yt)) ff::ffapply({
			out<-out+x[,i1:i2]%*% y[i1:i2,]
			},X=x,MARGIN=2)
	}
	if(outDim[1]>outDim[2] | (outDim[1]==p & outDim[2]==p)){
		if( (xt) & (!yt)) ff::ffapply({
			out[i1:i2,]<-crossprod(x[,i1:i2], y)
			},X=x,MARGIN=2)
		if((!xt) &  (yt)) ff::ffapply({
			out[i1:i2,]<-tcrossprod(x[i1:i2,], y)
			},X=x,MARGIN=1)
		if((!xt) & (!yt)) ff::ffapply({
			out[i1:i2,]<-x[i1:i2,]%*% y
			},X=x,MARGIN=1)
	}
	if(outDim[1]< outDim[2]){ 
		if( (xt) & (!yt))  ff::ffapply({
			out[,i1:i2]<-crossprod(x, y[,i1:i2])
			},X=y,MARGIN=2)
		if((!xt) &  (yt))  ff::ffapply({
			out[,i1:i2]<-tcrossprod(x, y[i1:i2,])
			},X=y,MARGIN=1)
		if((!xt) & (!yt))  ff::ffapply({
			out[,i1:i2]<- x %*% y[,i1:i2]
			},X=y,MARGIN=2)
		#Here, if y=NULL, we would've already gotten an error
	}
	return(out)
}

#chunk_mat <- function (x, RECORDBYTES = sum(ff::.rambytes[ff::vmode(x)]),
#  BATCHBYTES = getOption("ffbatchbytes"), ...) 
#{
#    n <- nrow(x)
#    if (n) {
#        l <- list(...)
#        if (is.null(l$from)) 
#            l$from <- 1L
#        if (is.null(l$to)) 
#            l$to <- n
#        if (is.null(l$by) && is.null(l$len)) {
#            b <- BATCHBYTES%/%RECORDBYTES
#            if (b == 0L) {
#                b <- 1L
#                warning("single record does not fit into BATCHBYTES")
#            }
#            l$by <- b
#        }
#        l$maxindex <- n
#        ret <- do.call(bit::chunk.default, l)
#    }
#    else {
#        ret <- list()
#    }
#    ret
#}

#smooth.construct_ff.ps.smooth.spec <- function(object, data, knots, ...)
#{
#  xr <- ffbase::range.ff(data[[object$term]])
#  print(xr)
#  stop()
#}


## Functions for sparse matrices.
sparse.matrix.index <- function(x, ...)
{
  if(is.null(dim(x)))
    return(NULL)
  if(inherits(x, "ff"))
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

## Setup sparse indices for various algorithms.
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
      x <- x + if(length(S[[j]]) < 1) {
        0
      } else {
        if(is.function(S[[j]])) {
          S[[j]](c("b" = rep(0, attr(S[[j]], "npar"))))
        } else {
          S[[j]]
        }
      }
    }
  }
  index.crossprod <- if(!symmetric) sparse.matrix.index(x, ...) else NULL
  setup <- list(
    "matrix" = index.matrix,
    "crossprod" = index.crossprod
  )
  if(!is.null(index.crossprod)) {
    # make block.index only if coefficients do not overlap
    tmp <- setup$crossprod[!duplicated(setup$crossprod), , drop = FALSE]
    l <- nrow(tmp)
    if(any(unique(tmp[duplicated(tmp, MARGIN = 0)]) > 0)){
      return(setup)
    } else {
      if((l > 1) & (l <= nrow(setup$crossprod))) {
        setup$block.index <- split(tmp, 1:l)
        setup$block.index <- lapply(1:l, function(i) setup$block.index[[i]][setup$block.index[[i]] > 0])
        setup$is.diagonal <- all(sapply(setup$block.index, length) == 1)
      }
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
    if(inherits(X, "Matrix")) {
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
  attr(ff, ".internal") <- TRUE
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
      if(!(prior %in% c("ig", "hc", "sd", "hn", "hn.lasso", "u")))
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
        "hn" = 0.1457644,
        "u" = 0.2723532
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
      },
      "u" = function(tau2) {
        1 - tau2 - exp(tau2 * 3 / theta - 3) / (1 + exp(tau2 * 3 / theta - 3))
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
          if(!is.null(x$xt[["pS"]]))
            K <- K + x$xt[["pS"]][[1]]
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
      if(!is.null(x$xt[["pm"]])) {
        pS <- if(!is.null(x$xt[["pS"]])) {
          x$xt[["pS"]]
        } else {
          if(!is.null(x$xt[["pSa"]])) {
            1 / tau2[length(tau2)] * x$xt[["pSa"]]
          } else 0
        }
        dP2 <- determinant(pS, logarithm = TRUE)
        dP2 <- as.numeric(dP2$modulus) * as.numeric(dP2$sign)
        lp2 <- 0.5 * dP2 - 0.5 * (t(gamma - x$xt[["pm"]]) %*% pS %*% (gamma - x$xt[["pm"]]))
        lp <- lp + lp2
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
        grad <- if(!is.null(x$binning)) {
          drop(crossprod(x$X[x$binning$match.index, , drop = FALSE], score)) + c(grad, grad2)
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
    formula$formula <- as.formula(formula$formula)
    env <- environment(formula$formula)
    if(is.null(env))
      env <- parent.frame()
    fcall$formula <- formula$formula
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
              if(is.null(dim(x[[i]]$smooth.construct[[k]]$X)))
                x[[i]]$smooth.construct[[k]]$X.dim
              else
                ncol(x[[i]]$smooth.construct[[k]]$X)
            } else x[[i]]$smooth.construct[[k]]$special.npar
            if(inherits(x[[i]]$smooth.construct[[k]], "nnet0.smooth")) {
              nfill <- x[[i]]$smooth.construct[[k]]$nodes * ncol(x[[i]]$smooth.construct[[k]]$X) +
                x[[i]]$smooth.construct[[k]]$nodes
            }
            tpar <- rep(fill[1], nfill)
            if(inherits(x[[i]]$smooth.construct[[k]], "nnet0.smooth")) {
              cn <- c(paste0("bb", 1:x[[i]]$smooth.construct[[k]]$nodes),
                names(unlist(x[[i]]$smooth.construct[[k]]$n.weights)))
            } else {
              cn <- colnames(x[[i]]$smooth.construct[[k]]$X)
            }
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
  default_fun <- c("no.transform", "opt_bfit", "sam_GMCMC", "samplestats", "results.bamlss.default")
  functions <- list()
  for(j in 1:length(foo)) {
    if(is.null(foo[[nf[j]]])) {
      foo[[nf[j]]] <- if(default_fun[j] != "no.transform") {
        get(default_fun[j], envir = .GlobalEnv)
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
      if(!is.logical(functions[[nf[j]]]) & !is.null(functions[[nf[j]]])) {
        stop(paste("argument", nf[j], "is not a function!"))
      } else {
        if(is.null(functions[[nf[j]]])) {
          functions[[nf[j]]] <- get(default_fun[j], envir = .GlobalEnv)
        } else {
          if(functions[[nf[j]]]) {
            functions[[nf[j]]] <- get(default_fun[j], envir = .GlobalEnv)
          }
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
      offset = bf$model.frame[["(offset)"]], ...)
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
    if(!is.null(opt$samples))
      bf$samples <- opt$samples
    ne <- names(opt)
    bf$model.stats <- opt[ne[!(ne %in% c("parameters", "fitted.values", "hessian", "samples"))]]
    rm(opt)
  }

  ## Start sampling.
  if(is.function(functions$sampler)) {
    ptm <- proc.time()
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
    elapsed <- c(proc.time() - ptm)[3]
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
    bf$model.stats$sampler$runtime <- elapsed
  } else {
    if(!is.null(bf$model.stats))
      bf$model.stats <- list("optimizer" = bf$model.stats)
  }

  #!# adjust model.frame for joint model with nonlinear
  args <- list(...)
  if(!is.null(args$nonlinear)) {
    if(args$nonlinear) {
      if(!is.null(bf$fitted.values$mu)){
        ## use estimated trajectories for setting up plot results
        bf$model.frame$mu <- bf$fitted.values$mu
      } else {
        ## use observed trajectories for setting up plot results
        bf$model.frame$mu <- bf$y[[1]][, "obs"]
      }
    }
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
  if(light)
    bf <- light_bamlss(bf)

  ## Remove ff directory?
  if(bf$delete) {
    if(dir.exists("ff_data_bamlss"))
      unlink("ff_data_bamlss", recursive = TRUE, force = TRUE)
  }


  bf$call <- match.call()
  class(bf) <- c("bamlss", "bamlss.frame", "list")
  attr(bf, "functions") <- functions

  bf
}

ff0 <- function(X, b, ...) { X %*% b }

light_bamlss <- function(object)
{
  if(!is.null(object$x)) {
    for(j in names(object$x)) {
      if(length(object$x[[j]]$smooth.construct)) {
        for(i in seq_along(object$x[[j]]$smooth.construct)) {
          object$x[[j]]$smooth.construct[[i]][c("X", "S", "Xr", "Xf", "binning", "prior", "grad", "hess", "boost.fit", "update", "propose")] <- NULL
          object$x[[j]]$smooth.construct[[i]][["xt"]][["binning"]] <- NULL
          environment(object$x[[j]]$formula) <- emptyenv()
          environment(object$x[[j]]$terms) <- emptyenv()
          environment(object$x[[j]]$fake.formula) <- emptyenv()
          if(!is.null(object$x[[j]]$smooth.construct[[i]][["margin"]])) {
            for(jj in seq_along(object$x[[j]]$smooth.construct[[i]][["margin"]])) {
              object$x[[j]]$smooth.construct[[i]][["margin"]][[jj]][c("X", "S", "Xr", "Xf", "binning", "prior", "grad", "hess", "boost.fit", "update", "propose")] <- NULL
              object$x[[j]]$smooth.construct[[i]][["margin"]][[jj]][["xt"]][["binning"]] <- NULL
            }
          }
          if(!is.null(object$x[[j]]$smooth.construct[[i]][["fit.fun"]])) {
            if(!is.null(attr(object$x[[j]]$smooth.construct[[i]][["fit.fun"]], ".internal")))
              object$x[[j]]$smooth.construct[[i]][["fit.fun"]] <- ff0
          }
          for(ff in names(object$x[[j]]$smooth.construct[[i]])) {
            if(ff != "fit.fun") {
              if(is.function(object$x[[j]]$smooth.construct[[i]][[ff]])) {
                environment(object$x[[j]]$smooth.construct[[i]][[ff]]) <- emptyenv()
              }
            }
          }
        }
      }
    }
  }
  object$y <- NULL
  object$fitted.values <- NULL
  object$formula <- as.character(object$formula)
  object$terms <- as.character(object$terms)
  return(object)
}

as.character.bamlss.formula <- as.character.bamlss.terms <- function(x, ...)
{
  if(inherits(x, "bamlss.formula.character"))
    return(x)
  for(i in seq_along(x)) {
    if(!is.null(names(x[[i]]))) {
      if(all(c("formula", "fake.formula") %in% names(x[[i]]))) {
        x[[i]]$formula <- deparse(x[[i]]$formula)
        if(length(x[[i]]$formula) > 1)
          x[[i]]$formula <- paste(x[[i]]$formula, collapse = " ")
        x[[i]]$fake.formula <- deparse(x[[i]]$fake.formula)
        if(length(x[[i]]$fake.formula) > 1)
          x[[i]]$fake.formula <- paste(x[[i]]$fake.formula, collapse = " ")
      } else {
        x[[i]] <- as.character.bamlss.formula(x[[i]])
      }
    }
  }
  attr(x, ".Environment") <- capture.output(attr(x, ".Environment"))
  class(x) <- c("bamlss.formula.character", "list")
  return(x)
}

formula.bamlss.formula.character <- function(x, ...)
{
  if(inherits(x, "bamlss.formula"))
    return(x)
  env <- attr(x, ".Environment")
  env <- if(grepl("globalenv", env, ignore.case = TRUE)) .GlobalEnv else NULL
  for(i in seq_along(x)) {
    if(!is.null(names(x[[i]]))) {
      if(all(c("formula", "fake.formula") %in% names(x[[i]]))) {
        if(length(x[[i]]$formula) > 1)
          x[[i]]$formula <- paste(x[[i]]$formula, collapse = " ")
        x[[i]]$formula <- as.formula(x[[i]]$formula, env = env)
        if(length(x[[i]]$fake.formula) > 1)
          x[[i]]$fake.formula <- paste(x[[i]]$fake.formula, collapse = " ")
        x[[i]]$fake.formula <- as.formula(x[[i]]$fake.formula)
      } else {
        x[[i]] <- formula.bamlss.formula.character(x[[i]], ...)
      }
    }
  }
  class(x) <- c("bamlss.formula", "list")
  attr(x, ".Environment") <- env
  return(x)
}

formula.bamlss.formula <- function(x, ...)
{
  if(inherits(x, "bamlss.formula"))
    return(x)
  if(inherits(x, "bamlss.formula.character"))
    return(formula(x))
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
    if(inherits(samples$formula, "bamlss.formula.character")) {
      samples$x <- design.construct(samples)
      samples$formula <- as.formula(samples$formula)
      samples$y <- model.frame(samples)[, attr(samples$formula, "response.name")]
    }
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
        for(i in nx) {
          par[[i]] <- .fitted.bamlss(i, x[[i]], samples)
          par[[i]] <- make.link2(family$links[i])$linkinv(par[[i]])
        }
        msamples <- matrix(apply(samples, 2, mean, na.rm = TRUE), nrow = 1)
        colnames(msamples) <- pn
        for(i in nx)
          mpar[[i]] <- make.link2(family$links[i])$linkinv(.fitted.bamlss(i, x[[i]], msamples))
        tpar <- mpar
        dev <- ll <- rep(NA, ncol(par[[1]]))
        for(j in 1:ncol(par[[1]])) {
          for(i in nx) {
            if(!is.null(ncol(par[[i]])))
              tpar[[i]] <- par[[i]][, j]
            else
              tpar[[i]] <- par[[i]]
          }
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
#stacker <- function(x, optimizer = opt_bfit, sampler = samplerJAGS, ...)
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
      data_ff <- ff::read.table.ffdf(file = data,
        na.strings = "", header = TRUE, sep = ",")
      if(!is.null(weights))
        data_ff[["(weights)"]] <- ff::as.ff(weights)
      if(!is.null(offset)) {
        noff <- names(offset)
        offset <- do.call("cbind", offset)
        if(is.null(dim(offset)))
          offset <- matrix(offset, ncol = 1)
        data_ff[["(offset)"]] <- ff::as.ff(offset, dim = dim(offset))
        colnames(data_ff[["(offset)"]]) <- noff
      }
      return(data_ff)
    }
    if(inherits(data, "ffdf"))
      return(data)
  } else data <- NULL

  if(inherits(formula, "bamlss.frame") | inherits(formula, "bamlss")) {
    if(!is.null(formula$model.frame))
      return(formula$model.frame)
    fcall <- formula$call
    fcall[[1L]] <- quote(bamlss.model.frame)
    formula$formula <- as.formula(formula$formula)
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
    wident <- FALSE
    if(!is.list(weights)) {
      weights <- rep(list(weights), length = length(family$names))
      names(weights) <- family$names
      wident <- TRUE
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
    attr(weights, "identical") <- wident
    data[["(weights)"]] <- weights
  }
  if(!is.null(offset)) {
    if(!is.list(offset)) {
      offset <- list(offset)
      names(offset) <- family$names[1]
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
    if(!is.null(attr(data, "na.action")))
      offset <- offset[-attr(data, "na.action"), , drop = FALSE]
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
    txt <- paste(tolower(family), type, sep = if(!is.null(type)) "_" else "")
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
        for(j in family$names) {
          eta[[j]] <- ff_eval(eta[[j]], FUN = function(x) { linkinv[[j]](x) },
            lower = c(-Inf, -10), upper = c(Inf, 10))
        }
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
    if(!is.null(family$d)) {
      family$loglik <- function(y, par, ...) {
        logdens <- family$d(y, par, log = TRUE)
        if(any(i <- !is.finite(logdens))) {
          logdens[i] <- -100
        }
        return(sum(logdens, na.rm = TRUE))
      }
    }
  }

  err01 <- .Machine$double.eps^(1/3)
  err02 <- err01 * 2
  err11 <- .Machine$double.eps^(1/4)
  err12 <- err11 * 2

  if(is.null(family$score) & !is.null(family$d))
    family$score <- list()
  for(i in family$names) {
    if(is.null(family$score[[i]]) & !is.null(family$d)) {
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

  if(is.null(family$hess) & !is.null(family$d))
    family$hess <- list()
  for(i in family$names) {
    if(is.null(family$hess[[i]]) & !is.null(family$d)) {
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
    specials <- c("s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "tx4", "la", "n", "h", "lf", "af", "lf.vd", "re", "peer", "fpc", "lin", "rb", "tree")
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
    if(any(j <- is.na(names(formula)))) {
      formula <- formula[!j]
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
  if((is.factor(data[[rn2[1]]]) | cat) & is.null(family$nocat)) {
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
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "tx4", "la", "n", "h", "lf", "af", "lf.vd", "re", "peer", "fpc", "lin", "rb", "tree"))
  tf <- terms(formula, specials = specials, keep.order = TRUE)
  ## sid <- unlist(attr(tf, "specials")) - attr(tf, "response")
  tl <- attr(tf, "term.labels")
  sid <- NULL
  for(j in specials) {
    i <- grep2(paste(j, "(", sep = ""), tl, fixed = TRUE)
    if(length(i)) {
      for(ii in i) {
        s1 <- strsplit(tl[ii], "")[[1]]
        s2 <- strsplit(paste(j, "(", sep = ""), "")[[1]]
        s1 <- paste(s1[1:length(s2)], collapse = "")
        s2 <- paste(s2, collapse = "")
        if(s1 == s2)
          sid <- c(sid, ii)
      }
    }
  }
  if(length(sid))
    sid <- sort(unique(sid))
  vars <- NULL
  if(rhs) {
    if(length(sid)) {
      vars <- tl[-sid]
      if(!length(vars))
        vars <- NULL
      for(j in tl[sid]) {
        if(any(grep2(paste0(c("af",  "lf.vd", "re", "peer", "fpc"), "("), j, fixed = TRUE)) ) {
          vars <- c(vars, all.vars(as.formula(paste("~", j)))[1])
        } else {
          tcall <- parse(text = j)[[1]]
          tcall[c("k","fx","bs","m","xt","id","sp","pc","d","mp","np","knots")] <- NULL
          tcall <- eval(tcall)
          vars <- c(vars, tcall$term)
          if(!is.null(tcall$by)) {
            if(tcall$by != "NA")
              vars <- c(vars, tcall$by)
          }
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
  if(is.null(sid)) {
    if(type == 2 & !lhs)
      type <- 1
  }
  if((type == 1) & !is.null(vars))
    vars <- all.vars(as.formula(paste("~", paste(vars, collapse = "+")), env = NULL))
  vars <- unique(vars)
  vars
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
  specials <- unique(c("s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "tx4", "la", "n", "h", "lf", "af", "lf.vd", "re", "peer", "fpc", "lin", "rb", "tree", specials))
  tf <- terms.formula2(formula, specials = specials, keep.order = FALSE)
  ## sid <- unlist(attr(tf, "specials")) - attr(tf, "response")
  tl <- attr(tf, "term.labels")
  sid <- NULL
  for(j in specials) {
    i <- grep2(paste(j, "(", sep = ""), tl, fixed = TRUE)
    if(length(i)) {
      for(ii in i) {
        s1 <- strsplit(tl[ii], "")[[1]]
        s2 <- strsplit(paste(j, "(", sep = ""), "")[[1]]
        s1 <- paste(s1[1:length(s2)], collapse = "")
        s2 <- paste(s2, collapse = "")
        if(s1 == s2)
          sid <- c(sid, ii)
      }
    }
  }
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
      if(length(i <- grep2(paste0(c("af",  "lf.vd", "re", "peer", "fpc"), "("), tl[j], fixed = TRUE)) ) {
        tl[j] <- paste0(c("af",  "lf.vd", "re", "peer", "fpc")[i], "(",
          all.vars(as.formula(paste("~", tl[j])))[1], ")")
      } else {
        tcall <- parse(text = tl[j])[[1]]
        id <- tcall$id
        tcall[c("k","fx","bs","m","xt","sp","pc","d","mp","np")] <- NULL
        tcall <- eval(tcall)
        if(is.null(tcall$label)) {
          if(!is.null(tcall$call))
            tcall <- eval(tcall$call)
        }
        if(!is.null(tcall$label))
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
        if(!is.null(id)) {
          if(grepl("mrf", id) | grepl("re", id)) {
            tlt <- paste0(",id='", id, "')")
            tl[j] <- gsub(")", tlt, tl[j], fixed = TRUE)
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
    rn <- response.name(formula, keep.functions = TRUE)
    ff <- fake.formula(formula, lhs = !(rn %in% family$names), specials = specials)
    return(list("formula" = formula, "fake.formula" = ff))
  }
}


## Get response name.
response.name <- function(formula, hierarchical = TRUE, keep.functions = FALSE, na.rm = FALSE)
{
  rn <- NA
  if(inherits(formula, "bamlss.frame")) {
    formula$formula <- as.formula(formula$formula)
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
      if(any(grepl("|", rn, fixed = TRUE)))
        rn <- all.vars(as.Formula(f))
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

response_name <- function(object, ...)
{
  return(response.name(object, ...))
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
trans_random <- randomize <- function(x)
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


trans_AR1 <- AR1 <- function(rho = 0.1) {
  function(x, ...) {
    if(bframe <- inherits(x, "bamlss.frame")) {
      if(is.null(x$x))
        stop("no 'x' object to randomize in 'bamlss.frame'!")
    }

    n <- if(is.null(dim(x$y))) length(x$y) else nrow(x$y)

    v <- rep(1, n)
    w <- rep(rho, n)
    v[1] <- sqrt(1 - rho^2)
    w[1] <- 0

    ind <- 1:(n - 1)

    trans_fun <- function(x)
    {
      if(m <- length(x$smooth.construct)) {
        for(j in 1:m) {
          x$smooth.construct[[j]]$X <- v * x$smooth.construct[[j]]$X -
            w * rbind(0, x$smooth.construct[[j]]$X[ind, , drop = FALSE])
        }
      }
      if(!is.null(x$model.matrix)) {
        x$model.matrix <- v * x$model.matrix - w * rbind(0, x$model.matrix[ind, , drop = FALSE])
      }
      x
    }

    elmts <- c("formula", "fake.formula")
    for(j in seq_along(x$x)) {
      if(!all(elmts %in% names(x$x[[j]]))) {
        for(i in seq_along(x$x[[j]]))
          x$x[[j]][[i]] <- trans_fun(x$x[[j]][[i]])
      } else x$x[[j]] <- trans_fun(x$x[[j]])
    }

    if(is.null(dim(x$y))) {
      x$y <- v * x$y - w * c(0, x$y[ind])
    } else {
      x$y <- v * x$y - w * rbind(0, x$y[ind, , drop = FALSE])
    }

    return(x)
  }
}


## Combine sample chains.
process.chains <- function(x, combine = TRUE, drop = FALSE, burnin = NULL, thin = NULL)
{
  if(is.null(x)) return(NULL)
  if(!is.list(x))
    x <- list(x)
  n <- sapply(x, nrow)
  if(!is.null(burnin) | !is.null(thin)) {
    for(i in seq_along(x)) {
      if(!is.null(burnin)) {
        if(burnin >= nrow(x[[i]]))
          stop("burnin >= nsamples!")
        x[[i]] <- as.mcmc(x[[i]][-c(1:burnin), , drop = FALSE])
      }
      if(!is.null(thin)) {
        n.iter <- nrow(x[[i]])
        iterthin <- as.integer(seq(1, n.iter, by = thin))
        x[[i]] <- as.mcmc(x[[i]][iterthin, , drop = FALSE])
      }
    }
  }
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
  names(rval) <- c("2.5%", "50%", "97.5%", "Mean")
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

  for(j in seq_along(data)) {
    if(is.matrix(data[[j]])) {
      if(ncol(data[[j]]) < 2)
        data[[j]] <- if(is.numeric(data[[j]])) as.numeric(data[[j]]) else drop(data[[j]])
    }
  }
  if(is.list(data))
    data <- as.data.frame(data)

  if(nt > 2) {
    message(paste0(".. not computing effect plot for term ", x$label, ", use predict() instead!"))
    return(NULL)
  }

  if((nt == 2) & (x$by != "NA")) {
    message(paste0(".. not computing effect plot for term ", x$label, ", use predict() instead!"))
    return(NULL)
  }

  if(x$by != "NA") grid <- NA
  if(!is.na(grid)) {
    if(grid < 0) {
      grid <- if(nt > 2) {
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

  if(!is.list(data)) {
    for(char in c("(", ")", "[", "]")) {
      colnames(data) <- gsub(char, ".", colnames(data), fixed = TRUE)
      x$by <- gsub(char, ".", x$by, fixed = TRUE)
    }
  }

  ## Data for rug plotting.
  rugp <- if(nt < 2 & rug) data[[x$term]] else NULL

  ## New x values for which effect should
  ## be calculated, n = 100.
  any_f <- any(sapply(data[tterms], is.factor))
  if(!is.na(grid)) {
    if(!any_f & !any(grepl("mrf", class(x))) &
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
        if(any(pip > 0)) {
          nd[pip < 1, ] <- NA
          nd <- na.omit(nd)
        }
      }
      data0 <- data
      data <- nd
    } else xsmall <- FALSE
  } else {
    if(is.data.frame(data)) {
      data0 <- data[, c(tterms, if(x$by != "NA") x$by else NULL), drop = FALSE]
      if(nt < 2) {
        if(x$by != "NA") {
          if(!is.factor(data[[x$by]]))
            data <- unique(data0)
        } else data <- unique(data0)
      }
      xsmall <- if((nrow(data) != nrow(data0)) & (nt < 2)) TRUE else FALSE
    } else xsmall <- FALSE
  }
  if(nrow(data) > 1e+05) {
    message(paste0(".. not computing effect plot for term ", x$label, ",\n.. .. too many observations, use predict() instead!"))
    return(NULL)
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
    smf <- quick_quantiles(X, psamples)[, c("2.5%", "Mean", "97.5%")]
  } else {
    if(is.null(FUN)) {
      FUN <- c95
    }
    if(nt < 2) {
      fsamples <- try(apply(psamples, 1, function(g) {
        f <- fit.fun(X, g, expand = FALSE, no.sparse.setup = (nrow(psamples) < 2))
        f
      }), silent = TRUE)
      if(inherits(fsamples, "try-error")) {
        fsamples <- try(apply(psamples, 1, function(g) {
          f <- X %*% g
          f
        }), silent = TRUE)
      }
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

  if(nrow(smf) < 2L) {
    smf <- t(smf)
    smf <- cbind(smf, smf, smf)
    colnames(smf) <- c("2.5%", "50%", "97.5%")
  }

  cnames <- colnames(smf)
  smf <- as.data.frame(smf)
  for(l in 1:nt) {
    if(is.matrix(data[[tterms[l]]])) {
      if(ncol(data[[tterms[l]]]) < 2)
        smf <- cbind(as.numeric(data[[tterms[l]]]), smf)
      else
        smf <- cbind(data[[tterms[l]]], smf)
    } else {
      smf <- cbind(drop(data[[tterms[l]]]), smf)
    }
  }
  names(smf) <- c(x$term, cnames)

  if(is.null(FUN)) {
    FUN <- c95
  }

  if(x$by != "NA") { ## FIXME: hard coded fix for plotting varying coefficient terms!
    if(!is.factor(data[[x$by]])) {
      X <- X / as.numeric(data[[x$by]])

      fsamples <- apply(psamples, 1, function(g) {
        fit.fun(X, g, expand = FALSE)
      })

      smf <- t(apply(fsamples, 1, FUN = FUN))

      cnames <- colnames(smf)
      smf <- as.data.frame(smf)
      for(l in 1:nt) {
        if(is.matrix(data[[tterms[l]]])) {
          smf <- cbind(as.numeric(data[[tterms[l]]]), smf)
        } else {
          smf <- cbind(data[[tterms[l]]], smf)
        }
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
  if(!any_f)
    smf <- unique(smf)

  class(smf) <- c(class(x), "data.frame")

  # this code does not work with environments/r6 classes:
  # x[!(names(x) %in% c("term", "label", "bs.dim", "dim"))] <- NULL

  # but this code does:
  specs <- list(term = x$term, label = x$label, bs.dim = x$bs.dim, dim = x$dim)
  mostattributes(specs) <- attributes(x)
  names(specs)[1:4] <- c("term", "label", "bs.dim", "dim")
  x <- specs

  attr(x, "qrc") <- NULL
  attr(smf, "specs") <- x
  class(attr(smf, "specs")) <- class(x)
  if(!is.list(data)) {
    attr(smf, "x") <- if(xsmall & nt < 2) data0[, tterms] else data[, tterms]
  } else {
    for(t in tterms)
      data[[t]] <- if(is.matrix(data[[t]])) as.numeric(data[[t]]) else drop(data[[t]])
    attr(smf, "x") <- as.data.frame(data)[, tterms]
  }
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
  family <- object$family

  object$formula <- as.formula(object$formula)
  if(any(i <- is.na(names(object$formula)))) {
    rn <- attr(object$formula, "response.name")
    object$formula <- object$formula[!i]
    class(object$formula) <- c("bamlss.formula", "list")
    if(!is.null(rn))
      attr(object$formula, "response.name") <- rn
  }

  if(!missing(newdata)) {
    if(!is.null(newdata)) {
      if(nrow(newdata) < 1)
        stop("newdata is empty with nrow < 1!")
    }
  }

  ## If data have been scaled (scale.d=TRUE)
  if(!missing(newdata) & ! is.null(attr(object$model.frame, 'scale')) ) {
    sc <- attr(object$model.frame, 'scale')
    for ( name in unique(unlist(lapply(sc,names))) ) {
      newdata[,name] <- (newdata[,name] - sc$center[name] ) / sc$scale[name]
    }
  }
  FUN2 <- function(x, ...) FUN(x)
  if(missing(newdata))
    newdata <- NULL
  if(!is.null(family$predict)) {
    if(is.function(family$predict)) {
      return(family$predict(object = object, newdata = newdata, model = model, term = term,
        intercept = intercept, type = type, FUN = FUN2, trans = trans, what = what, nsamps = nsamps,
        verbose = verbose, drop = drop, cores = cores, chunks = chunks, ...))
    }
  }
  if(is.null(object$x)) {
    object$x <- smooth.construct(object)
  }
  if(any(i <- is.na(names(object$x))))
    object$x <- object$x[!i]
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
  if(inherits(newdata, "ffdf"))
    newdata <- as.data.frame(newdata)
  nn <- names(newdata)
  rn_nn <- rownames(newdata)
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
            if(length(k)) {
              if(length(k) > 1) {
                enames[[i]] <- c(enames[[i]], tl[[i]][k])
              } else {
                enames[[i]] <- c(enames[[i]], tl[[i]][k])
              }
            }
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
  uenames <- unique(unlist(enames))
  uenames <- gsub("splines::", "", uenames, fixed = TRUE)
  ff <- as.formula(paste("~", paste(uenames, collapse = "+")), env = NULL)
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
    if(!is.null(dim(object$parameters)) & !is.null(object$model.stats$optimizer$boost_summary) & is.null(list(...)$mstop)) {
      samps <- object$parameters
      FUN <- function(x) { x }
      FUN2 <- function(x, ...) FUN(x)
    } else {
      samps <- parameters(object, model = model, list = FALSE, extract = TRUE, ...)
      if(is.null(dim(samps))) {
        cn <- names(samps)
        samps <- matrix(samps, nrow = 1)
        colnames(samps) <- cn
      }
      mstop <- list(...)$mstop
      if(!is.null(mstop)) {
        if(!is.null(object$model.stats$optimizer$boost_summary) & (length(mstop) > 1)) {
          FUN <- function(x) { x }
          FUN2 <- function(x, ...) FUN(x)
        }
      }
    }
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

  na.action <- NULL
  if(any(is.na(newdata))) {
    warning("NA values in newdata, calling na.omit()!")
    newdata <- na.omit(newdata)
    na.action <- attr(newdata, "na.action")
  }

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

  for(j in seq_along(pred)) {
    if(!is.null(dim(pred[[j]]))) {
      if(ncol(pred[[j]]) < 2L) {
        pred[[j]] <- as.vector(pred[[j]])
        if(!is.null(rn_nn)) {
          names(pred[[j]]) <- rn_nn
        } else {
          names(pred[[j]]) <- NULL
        }
      } else {
        cn <- colnames(pred[[j]])
        pred[[j]] <- as.data.frame(pred[[j]])
        names(pred[[j]]) <- cn
        if(is.null(rn_nn)) {
          rownames(pred[[j]]) <- NULL
        } else {
          rownames(pred[[j]]) <- rn_nn
        }
      }
    }
  }

  if((length(pred) < 2) & drop)
    pred <- pred[[1]]

  if(!is.null(na.action))
    attr(pred, "na.action") <- na.action

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
        if(is.character(data[[j]]))
          data[[j]] <- as.factor(data[[j]])
        if(is.factor(data[[j]])) {
          if(nlevels(data[[j]]) < 2) {
            xlev <- levels(data[[j]])
            if(!any(grepl(paste0(j, xlev), snames, fixed = TRUE))) {
              xlev2 <- grep(j, snames, fixed = TRUE, value = TRUE)
              xlev2 <- gsub(paste0(id, ".p.", j), "", xlev2, fixed = TRUE)
              levels(data[[j]]) <- c(xlev, xlev2)
              data[[j]] <- relevel(data[[j]], ref = xlev)
            } else {
              levels(data[[j]]) <- c(xlev, "NA")
              data[[j]] <- relevel(data[[j]], ref = "NA")
            }
          }
        }
        X <- try(model.matrix(f, data = data), silent = TRUE)
        if(inherits(X, "try-error")) {
          attr(data, "terms") <- NULL
          X <- model.matrix(f, data = data)
        }
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
          if(is.factor(data[[j]])) {
            xlev <- colnames(X)
            sn <- sn[grep2(paste(id, ".p.", xlev, ".", sep = ""), sn2, fixed = TRUE)]
          } else {
            sn <- sn[grepl(paste(id, ".p.", j, ".", sep = ""), sn2, fixed = TRUE)]
          }
          sn <- sn[!grepl(":", sn, fixed = TRUE)]
        }
        # to be checked if collision with other models
        if(!length(sn)) {
          sn <- snames[grep2(paste(id, "p.model.matrix", j, sep = "."), snames, fixed = TRUE)]
          if(!length(sn)) {
            sn <- snames[grep2(paste(id, "p", j, sep = "."), snames, fixed = TRUE)]
          }
        }
        if(ncol(X) > length(sn)) {
          sn2 <- gsub(paste(id, "p.", sep = "."), "", sn, fixed = TRUE)
          if(any(grepl("model.matrix.", sn2))) sn2 <- gsub("model.matrix.", "", sn2, fixed = TRUE)
          X <- X[, sn2, drop = FALSE]
        }
        eta <- eta + fitted_matrix(X, samps[, sn, drop = FALSE])
      } else {
        ## sn <- snames[grep2(paste(id, "p", j, sep = "."), snames, fixed = TRUE)]
        sn <- paste(id, "p", j, sep = ".")
        if(!any(sn %in% snames))
          sn <- snames[grep2(paste(id, "p.model.matrix", j, sep = "."), snames, fixed = TRUE)]
        if(length(sn))
          eta <- eta + fitted_matrix(matrix(1, nrow = nrow(data), ncol = 1), samps[, sn, drop = FALSE])
      }
    }
  }

  grep3 <- function(x, y) {
    m <- grep(x, y, fixed = TRUE, value = TRUE)
    if(length(m) > 1) {
      x <- paste0(strsplit(x, "(", fixed = TRUE)[[1]][1], "(")
      y2 <- sapply(strsplit(y, "(", fixed = TRUE), function(x) { paste0(x[1], "(") })
      m <- y[y2 %in% x]
    }
    m
  }
  if(length(i <- grep("s.", ec))) {
    for(j in enames2[i]) {
      for(jj in grep3(j, names(x))) {
        sn <- snames[grep2(paste0(paste(id, "s", jj, sep = "."), "."), snames, fixed = TRUE)]
        if(inherits(x[[jj]], "ff_smooth.smooth.spec")) {
          if(!is.null(x[[jj]]$orig.class)) {
            if("nnet0.smooth" %in% x[[jj]]$orig.class)
              class(x[[jj]]) <- x[[jj]]$orig.class
          }
        }
        if(!inherits(x[[jj]], "no.mgcv") & !inherits(x[[jj]], "special")) {
          if(is.null(x[[jj]]$mono))
            x[[jj]]$mono <- 0
          if(!is.null(x[[jj]]$margin)) {
            for(mjj in seq_along(x[[jj]]$margin)) {
              if(is.null(x[[jj]]$margin[[mjj]]$mono))
                x[[jj]]$margin[[mjj]]$mono <- 0
            }
          }
          random <- FALSE
          if(!is.null(x[[jj]]$is.refund)) {
            rfcall <- x[[jj]]$refund.call
            tfm <- eval(parse(text = rfcall), envir = data)
            tfme <- eval(tfm$call, envir = tfm$data)
            X <- smoothCon(tfme, data = tfm$data, n = nrow(tfm$data[[1L]]),
              knots = NULL, absorb.cons = TRUE,scale.penalty=TRUE)[[1]]$X
            rm(tfm)
            rm(tfme)
          } else {
            X <- PredictMat(x[[jj]], data)
            random <- if(!is.null(x[[jj]]$margin)) {
              any(sapply(x[[jj]]$margin, function(z) { inherits(z, "random.effect") }))
            } else inherits(x[[jj]], "random.effect")
          }
          if(any(sn %in% colnames(samps))) {
            if(random) {
              if(ncol(X) == ncol(samps[, sn, drop = FALSE]))
                eta <- eta + fitted_matrix(X, samps[, sn, drop = FALSE])
            } else {
              eta <- eta + fitted_matrix(X, samps[, sn, drop = FALSE])
            }
          } else {
            warning(paste("model term", j, "not in samples, prediction is set to 0!"))
          }
        } else {
          if(is.null(x[[jj]]$PredictMat)) {
            X <- PredictMat(x[[jj]], data)
          } else {
            if(inherits(x[[jj]], "nnet0.smooth")) {
              X <- Predict.matrix.nnet0.smooth(x[[jj]], data)
              ## X <- x[[jj]]$getZ(X, samps[1L, sn, drop = FALSE])
            } else {
              X <- x[[jj]]$PredictMat(x[[jj]], data)
            }
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
      ptxt <- paste(id, if(j != "model.matrix") "s" else "p", j, sep = ".")
      if(j != "model.matrix")
        ptxt <- paste0(ptxt, ".")
      sn <- grep(ptxt, snames, fixed = TRUE, value = TRUE)
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
  rval <- smoothCon(object, data, knots, absorb.cons = acons,scale.penalty=TRUE)
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

    ## Compute acceptance probability.
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
  env <- try(environment(formula), silent = TRUE)
  if(inherits(env, "try-error")) {
    if(grepl("not found", env)) {
      formula <- deparse(substitute(formula), backtick = TRUE, width.cutoff = 500)
      formula2 <- try(dynGet(formula), silent = TRUE)
      if(!inherits(formula2, "try-error")) {
        if(length(formula2) < 2L)
          formula <- formula2
      }
    }
  }
  formula <- deparse(substitute(formula), backtick = TRUE, width.cutoff = 500)
  formula <- gsub("[[:space:]]", "", formula)
  formula <- gsub('"', '', formula, fixed = TRUE)
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
  if(!any(grepl("+", formula, fixed = TRUE)) & !any(grepl("-", formula, fixed = TRUE)) & is.null(label)) {
    if(any(grepl(":", formula, fixed = TRUE))) {
      label <- paste0("la(", paste(as.character(formula), collapse = ""), ")")
    } else {
      label <- paste("la(", paste(all.vars.formula(as.formula(formula)), collapse = "+"), ")", sep = "")
    }
  }
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
  enet <- if(is.null(object$xt[["enet"]])) FALSE else object$xt[["enet"]]
  fuse <- if(is.null(object$xt[["fuse"]])) FALSE else object$xt[["fuse"]]
  standardize <- if(is.null(object$xt[["standardize"]])) FALSE else object$xt[["standardize"]]
  standardize01 <- if(is.null(object$xt[["standardize01"]])) FALSE else object$xt[["standardize01"]]
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
  object$standardize <- if(!fuse) TRUE else standardize
  object$standardize01 <- standardize01
  
  contr <- object$xt$contrast.arg
  if(is.null(contr))
    contr <- "contr.treatment"
  
  data <- as.data.frame(data)
  nobs <- nrow(data)
  tl <- term.labels2(terms(object$formula), intercept = FALSE, list = FALSE)
  if(any(grepl("la(", tl, fixed = TRUE)))
    tl <- object$term
  object$X <- df <- group <- list()
  object$lasso <- list("trans" = list())
  k <- 1
  for(j in tl) {
    if(is.factor(data[[j]])) {
      contr.list <- list()
      contr.list[[j]] <- contr
    } else contr.list <- NULL
    object$X[[j]] <- as.matrix(model.matrix(as.formula(paste("~", j)), data = data,
      contrasts = contr.list))
    if(length(i <- grep("Intercept", colnames(object$X[[j]]))))
      object$X[[j]] <- object$X[[j]][, -i, drop = FALSE]
    is_f <- is.factor(data[[j]])
    if(is_f) {
      group[[j]] <- k:(k + ncol(object$X[[j]]) - 1)
    } else {
      group[[j]] <- NA
    }
    k <- k + ncol(object$X[[j]])
    j2 <- NULL
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
        if(standardize) {
          if(standardize01) {
            xmin <- apply(object$X[[j]], 2, min, na.rm = TRUE)
            xmax <- apply(object$X[[j]], 2, max, na.rm = TRUE)
            if((xmax - xmin) < sqrt(.Machine$double.eps)) {
              xmin <- 0
              xmax <- 1
            }
            for(jj in 1:ncol(object$X[[j]])) {
              object$X[[j]][, jj] <- (object$X[[j]][, jj] - xmin[jj]) / (xmax[jj] - xmin[jj])
              if(!is.null(object$xt$m1p1))
                object$X[[j]][, jj] <- object$X[[j]][, jj] * 2 - 1
            }
            object$lasso$trans[[j]] <- list("xmin" = xmin, "xmax" = xmax)
          } else {
            object$X[[j]] <- scale(object$X[[j]])
            object$lasso$trans[[j]] <- list(
              "center" = attr(object$X[[j]], "scaled:center"),
              "scale" = attr(object$X[[j]], "scaled:scale")
            )
          }
        }
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
        df[[j]] <- colSums(object$X[[j]] == 1)
      } else {
        if(is.null(j2)) {
          df[[j]] <- rep(nobs, ncol(object$X[[j]]))
        } else {
          df[[j]] <- colSums(object$X[[j]] == 1)
        }
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
  if(!fuse & !ridge) {
    if(object$type == "single") {
      object$S[[1]] <- function(parameters, fixed.hyper = NULL) {
        b <- get.par(parameters, "b")
        w <- rep(1, length(b))
        if(!is.null(fixed.hyper)) {
          w <- fixed.hyper
        } else {
          if(length(i <- grep("lasso", names(parameters)))) {
            w <- parameters[i]
          }
        }
        A <- df / sqrt(b^2 + const)
        for(j in seq_along(group)) {
          if(all(!is.na(group[[j]]))) {
            A[group[[j]]] <- df[group[[j]]] / rep(sqrt(sum(b[group[[j]]]^2) + const), length(group[[j]]))
          }
        }
        ## FIXME: adaptive weights: A <- A * MLpen ## 1 / abs(beta)
        A <- A * 1 / abs(w)
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
  object$xt$gfx <- if(is.null(object$xt$gfx)) FALSE else object$xt$gfx
  if(fuse & !ridge) {
    k <- ncol(object$X)
    if(fuse_type == "nominal") {
      if(any(grep(":", colnames(object$X), fixed = TRUE)) & object$xt$gfx) {
        cn <- colnames(object$X)
        cn <- strsplit(cn, ":")
        cn <- do.call("rbind", cn)
        ucn1 <- unique(cn[, 1L])
        ucn2 <- unique(cn[, 2L])
        combis <- combn(length(ucn2), 2)
        Af <- NULL ## matrix(0, ncol = ncol(combis), nrow = k)
        nd <- nrow(cn)
        for(ff in 1:ncol(combis)) {
          for(ii in 1:length(ucn1)) {
            tAf <- rep(0, k)
            tAf[cn[, 2L] == ucn2[combis[1, ff]] & cn[, 1L] == ucn1[ii]] <- 1
            tAf[cn[, 2L] == ucn2[combis[2, ff]] & cn[, 1L] == ucn1[ii]] <- -1
            Af <- cbind(Af, tAf)
          }
        }
        Af <- cbind(diag(nrow(Af)), Af)
      } else {
        Af <- matrix(0, ncol = choose(k, 2), nrow = k)
        combis <- combn(k, 2)
        for(ff in 1:ncol(combis)){
          Af[combis[1, ff], ff] <- 1
          Af[combis[2, ff], ff] <- -1
        }
        Af <- cbind(diag(k), Af)
      }
    } else {
      Af <- diff(diag(k + 1))
      Af[1, 1] <- 1
      Af <- Af[, -ncol(Af), drop = FALSE]
    }
    beta <- object$xt$beta
    w <- rep(0, length = ncol(Af))
    nref <- nobs - sum(df)
    if(!object$xt$gfx) {
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
    } else {
      w <- rep(1, length = ncol(Af))
    }
    names(w) <- paste0("w", 1:length(w))
    object$Af <- Af
    object$S[[ls <- length(object$S) + 1]] <- function(parameters, fixed.hyper = NULL) {
      b <- get.par(parameters, "b")
      if(!is.null(fixed.hyper)) {
        w <- fixed.hyper
      } else {
        if(length(i <- grep("lasso", names(parameters)))) {
          w <- parameters[i]
        }
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

  if(ridge)
    object$S <- list(diag(1, ncol(object$X)))

  if(enet)
    object$S[[length(object$S) + 1L]] <- diag(1, ncol(object$X))

  object$xt[["prior"]] <- "ig"
  object$xt[["a"]] <- 1e-10
  object$xt[["b"]] <- 1e+04
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

  if(is.null(object$xt[["binning"]]))
    object$xt[["binning"]] <- TRUE
  if(is.null(object$xt[["df"]]))
    object$xt[["df"]] <- if(!ridge) ceiling(ncol(object$X) * 0.9) else ceiling(ncol(object$X) * 0.3)
  object$ctype <- switch(object$type,
    "single" = 0,
    "multiple" = 1
  )
  object$C <- matrix(nrow = 0, ncol = ncol(object$X))
  object$side.constrain <- FALSE
  object$boost.fit <- boost_fit

  if(object$fixed)
    object$S <- NULL

  if(!ridge) {
    object$boost.fit <- function(x, y, nu, hatmatrix = TRUE, ...) {
      ## Compute reduced residuals.
      xbin.fun(x$binning$sorted.index, rep(1, length = length(y)), y, x$weights, x$rres, x$binning$order)
  
      ## Compute mean and precision.
      XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)
      if(x$fixed) {
        P <- matrix_inv(XWX, index = x$sparse.setup)
      } else {
        S <- 0
        tau2 <- 0.01 ## get.state(x, "tau2")
        for(j in seq_along(x$S))
          S <- S + 1 / tau2[j] * if(is.function(x$S[[j]])) x$S[[j]](c(x$state$parameters, x$fixed.hyper)) else x$S[[j]]
        P <- matrix_inv(XWX + S, index = x$sparse.setup)
      }
  
      ## New parameters.
      g <- nu * drop(P %*% crossprod(x$X, x$rres))
  
      ## Finalize.
      x$state$parameters <- set.par(x$state$parameters, g, "b")
      x$state$parameters <- set.par(x$state$parameters, 1e-05, "tau2")
      x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))
      x$state$rss <- sum((x$state$fitted.values - y)^2)

      if(hatmatrix)
        x$state$hat <- nu * x$X %*% P %*% t(x$X)
  
      return(x$state)
    }
  }

  if(is.null(object$xt$alpha))
    object$xt$alpha <- 1
  if(is.null(object$xt$nlambda))
    object$xt$nlambda <- 100
  if(is.null(object$xt$lambda.min.ratio))
    object$xt$lambda.min.ratio <- 1e-20
  if(is.null(object$xt$update)) {
    object$xt$update <- bfit_iwls
  } else {
    if(is.character(object$xt$update)) {
      if(object$xt$update == "glmnet")
        object$xt$update <- bfit_glmnet
    }
  }

  class(object) <- "lasso.smooth"

  object
}

Predict.matrix.lasso.smooth <- function(object, data)
{
  data <- as.data.frame(data)
  tl <- term.labels2(terms(object$formula), intercept = FALSE, list = FALSE)
  contr <- object$xt$contrast.arg
  if(is.null(contr))
    contr <- "contr.treatment"
  X <- list()
  for(j in tl) {
    if(is.factor(data[[j]])) {
      contr.list <- list()
      contr.list[[j]] <- contr
    } else contr.list <- NULL
    X[[j]] <- as.matrix(model.matrix(as.formula(paste("~", j)), data = data, contrasts = contr.list))
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
      if(object[["standardize"]]) {
        if(object[["standardize01"]]) {
          xmin <- object$lasso$trans[[j]]$xmin
          xmax <- object$lasso$trans[[j]]$xmax
          for(jj in 1:ncol(X[[j]]))
            X[[j]][, jj] <- (X[[j]][, jj] - xmin[jj]) / (xmax[jj] - xmin[jj])
          if(!is.null(object$xt$m1p1))
            X[[j]][, jj] <- X[[j]][, jj] * 2 - 1
        } else {
          if(!is.null(object$lasso$trans[[j]]$center))
            X[[j]] <- (X[[j]] - object$lasso$trans[[j]]$center) / object$lasso$trans[[j]]$scale
        }
      }
    } else {
      X[[j]] <- X[[j]] %*% object$lasso$trans[[j]]$blockscale
    }
  }
  return(do.call("cbind", X))
}


## Linear effects constructor.
lin <- function(...)
{
  ret <- la(...)
  ret$label <- gsub("la(", "lin(", ret$label, fixed = TRUE)
  class(ret) <- "linear.smooth.spec"
  ret
}

smooth.construct.linear.smooth.spec <- function(object, data, knots, ...)
{
  object$X <- model.matrix(object$formula, data = as.data.frame(data))[, -1, drop = FALSE]
  colnames(object$X) <- paste0("b.", colnames(object$X))

  center <- scale <- rep(NA, ncol(object$X))
  for(j in 1:ncol(object$X)) {
    if(length(unique(object$X[, j])) > 3) {
      center[j] <- mean(object$X[, j])
      scale[j] <- sd(object$X[, j])
      object$X[, j] <- (object$X[, j] - center[j]) / scale[j]
    }
  }

  object$scale <- list("center" = center, "scale" = scale)

  ridge <- if(is.null(object$xt$ridge)) TRUE else object$xt$ridge
  neigh <- if(is.null(object$xt$neigh)) FALSE else object$xt$neigh
  if(ridge | neigh) {
    if(ridge)
      object$S <- list(diag(ncol(object$X)))
    else
      object$S <- list(crossprod(diff(diag(ncol(object$X)))))
    object$xt$df <- if(is.null(object$xt$df)) min(c(floor(length(object$scale$center) / 2), 2)) else object$xt$df
  } else {
    object$fixed <- TRUE
    object$fxsp <- TRUE
  }

  object$N <- apply(object$X, 2, function(x) {
    return((1/crossprod(x)) %*% t(x))
  })

  object$boost.fit <- function(x, y, nu, hatmatrix = FALSE, weights = NULL, nthreads = 1, ...) {
    ## process weights.
    if(!is.null(weights))
      stop("weights is not supported!")

    bf <- boost_fit_nnet(nu, x$X, x$N, y, x$binning$match.index, nthreads = nthreads)

    j <- which.min(bf$rss)
    g2 <- rep(0, length(bf$g))
    g2[j] <- bf$g[j]
  
    ## Finalize.
    x$state$parameters <- set.par(x$state$parameters, g2, "b")
    x$state$fitted.values <- bf$fit[, j]

    x$state$rss <- bf$rss[j]

    if(hatmatrix) {
      x$state$hat <- nu * x$X[, j] %*% (1/crossprod(x$X[, j])) %*% t(x$X[, j])
    }
  
    return(x$state)
  }

  class(object) <- c("linear.smooth", "mgcv.smooth")

  return(object)
}

Predict.matrix.linear.smooth <- function(object, data)
{
  X <- model.matrix(object$formula, data = as.data.frame(data))[, -1, drop = FALSE]
  colnames(X) <- paste0("b.", colnames(X))
  for(j in 1:ncol(X)) {
    if(!is.na(object$scale$center[j]))
      X[, j] <- (X[, j] - object$scale$center[j]) / object$scale$scale[j]
  }
  return(X)
}


#smooth.construct.ctree.smooth.spec <- function(object, data, knots, ...)
#{
#  object$X <- as.data.frame(data)[, object$term, drop = FALSE]
#  object$formula <- as.formula(paste0("y~", paste0(object$term, collapse = "+")))
#  if(is.null(object$xt$maxdepth))
#    object$xt$maxdepth <- 3
#  object$state <- list()
#  object$state$fitted.values <- rep(0, nrow(object$X))
#  object$state$parameters <- rep(NA, object$xt$maxdepth^2)
#  object$state$rules <- rep("", object$xt$maxdepth^2)
#  object$boost.fit <- function(x, y, nu, hatmatrix = FALSE, weights = NULL, ...) {
#    ## process weights.
#    if(!is.null(weights))
#      stop("weights is not supported!")

#    b <- ctree(object$formula, data = object$X, control = ctree_control(maxdepth = x$xt$maxdepth))
#  
#    ## Finalize.
#    rules <- partykit:::.list.rules.party(b)
#    x$state$fitted.values <- nu * as.numeric(fitted(b)[, "(response)"])
#    par <- unique(x$state$fitted.values)
#    x$state$parameters[1:length(rules)] <- par
#    x$state$rules[1:length(rules)] <- rules[names(par)]

#    x$state$rss <- bf$rss[j]

#    if(hatmatrix) {
#      x$state$hat <- nu * x$X[, j] %*% (1/crossprod(x$X[, j])) %*% t(x$X[, j])
#    }
#  
#    return(x$state)
#  }

#  class(object) <- c("linear.smooth", "mgcv.smooth")

#  return(object)
#}


## Neural networks.
n <- function(..., k = 10, type = 2)
{
  ret <- la(..., k = k)
  ret$label <- gsub("la(", "n(", ret$label, fixed = TRUE)
  if(!is.null(ret$xt$id)) {
    lab <- strsplit(ret$label, "")[[1]]
    lab <- paste(c(lab[-length(lab)], paste(",id='", ret$xt$id, "')", sep = "")), collapse = "", sep = "")
    ret$label <- lab
  }
  if(is.null(ret$xt$tp))
    ret$xt$tp <- FALSE
  if(type != 1) {
    if(type > 2) {
      ret$special <- TRUE
      class(ret) <- "nnet4.smooth.spec"
    } else {
      class(ret) <- "nnet2.smooth.spec"
    }
  } else {
    class(ret) <- "nnet0.smooth.spec"
  }
  ret
}


t.weights <- function(x, y, n = 20, k = 100, dropout = NULL, weights = NULL, wts = NULL, maxit = 100, ...)
{
  warn <- getOption("warn")
  options("warn" = -1)
  on.exit(options("warn" = warn))
  nw <- ncol(x)
  w0 <- build_net_w(cbind(1, x), y, k = k, n = n, plot = FALSE, weights = weights, wts = wts, maxit = maxit)
  w <- list()
  for(i in 1:ncol(w0)) {
    w[[i]] <- w0[, i]
    names(w[[i]]) <- paste0("bw", i, "_w", 0:nw)
  }
  return(w)
}


gZ <- function(x, w) {
  if(is.data.frame(x))
    x <- as.matrix(x)
  if(!is.matrix(x))
    x <- matrix(x, ncol = 1)
  x <- cbind(1, as.matrix(x))
  Z <- matrix(0, nrow = nrow(x), ncol = length(w))
  for(j in 1:length(w)) {
    Z[, j] <- 1 / (1 + exp(-(x %*% w[[j]])))
  }
  Z
}


n.weights <- function(nodes, k, r = NULL, s = NULL, type = c("sigmoid", "gauss", "softplus", "cos", "sin"), x = NULL, ...)
{
  if(!is.null(y <- list(...)$y) & !is.null(x) & !is.null(list(...)$wm)) {
    wts <- t.weights(x, y, k = nodes, n = list(...)$tntake, dropout = list(...)$dropout, weights = list(...)$weights,
      wts = list(...)$wts, maxit = list(...)$maxit)
    return(wts)
  }
  type <- match.arg(type)
  dropout <- list(...)$dropout
  if(!is.null(dropout)) {
    if(is.logical(dropout)) {
      dropout <- if(dropout) 0.1 else NULL
    } else {
      dropout <- as.numeric(dropout)
    }
  }
  if(inherits(nodes, "bamlss")) {
    if(!is.null(nodes$parameters)) {
      rval <- list()
      for(j in nodes$family$name) {
        cb <- coef.bamlss(nodes, model = j, pterms = FALSE, hyper.parameters = FALSE, ...)
        if(any(i <- grepl(paste0(j, ".s.n("), names(cb), fixed = TRUE) | grepl(paste0(j, ".s.rb("), names(cb), fixed = TRUE))) {
          cb <- cb[i]
          cb1 <- cb[grep(".b", names(cb), fixed = TRUE)]
          terms <- unique(paste0(sapply(strsplit(names(cb1), ").b", fixed = TRUE), function(x) { x[1] }), ")"))
          if(length(terms) > 1) {
            rval[[j]] <- list()
            for(tj in terms) {
              cb2 <- cb1[grep(tj, names(cb1), fixed = TRUE)]
              id <- as.integer(sapply(strsplit(names(cb2), ").b", fixed = TRUE), function(x) { x[2] }))
              id <- id[abs(cb2) > 1e-10]
              rval[[j]][[tj]] <- id
#              weights <- vector(mode = "list", length = length(id))
#              for(i in seq_along(id))
#                weights[[i]] <- as.numeric(cb[grep(paste0(tj, ".bw", id[i], "_w"), names(cb), fixed = TRUE)])
#              attr(rval[[j]][[tj]], "weights") <- weights
            }
          } else {
            id <- as.integer(sapply(strsplit(names(cb1), ").b", fixed = TRUE), function(x) { x[2] }))
            id <- id[abs(cb1) > 1e-10]
            rval[[j]] <- list()
            i <- strsplit(terms, ".", fixed = TRUE)[[1]][3]
            ww <- nodes$x[[j]]$smooth.construct[[i]]$n.weights
            for(jj in names(ww)) {
              rval[[j]][[jj]] <- ww[[jj]][id]
            }
          }
        }
      }
      class(rval) <- c("list", "n.weights")
      return(rval)
    } else return(NULL)
  }
  if(type == "relu")
    type <- "softplus"
  k <- k + 1
  if(is.null(r) & is.null(s)) {
    rint <- list(...)$rint
    sint <- list(...)$sint
    if(type == "sigmoid") {
      if(is.null(rint))
        rint <- 0.05
      if(is.null(sint))
        sint <- 5
    }
    if(type == "gauss") {
      if(is.null(rint))
        rint <- 0.05
      if(is.null(sint))
        sint <- 5
    }
    if(type == "softplus") {
      if(is.null(rint))
        rint <- 0.05
      if(is.null(sint))
        sint <- 5
    }
    if(type == "cos" | type == "sin") {
      if(is.null(rint))
        rint <- 0.05
      if(is.null(sint))
        sint <- 5
    }
    sint <- sort(rep(sint, length.out = 2))
    rint <- sort(rep(rint, length.out = 2))
  } else {
    sint <- sort(rep(s, length.out = 2))
    rint <- sort(rep(r, length.out = 2))
  }
  nodes2 <- floor(nodes/2)
#  r <- runif(nodes, rint[1], rint[2])
#  s <- runif(nodes, sint[1], sint[2])
  r <- sort(rep(c(rep(rint[1], nodes2), rep(rint[2], nodes2)), length.out = nodes))
  s <- sort(rep(c(rep(sint[1], nodes2), rep(sint[2], nodes2)), length.out = nodes))
  if(type == "sigmoid") {
    if(any(r >= 0.5))
      r[r >= 0.5] <- 0.495
  }
  if(type == "gauss") {
    if(any(r >= 1))
      r[r >= 1] <- 0.99
  }
  if(type == "softplus") {
    if(any(r >= log(2)))
      r[r >= 1] <- log(2) - 0.001
  }
  if(type == "cos" | type == "sin") {
    if(any(r > 1))
      r[r > 1] <- 1
    if(any(r < -1))
      r[r < -1] <- -1
  } else {
    if(any(r < 0))
      r[r < 0] <- 0.01
  }
  if(any(s <= 1))
    s[s <= 1] <- 1.01
  r <- rep(r, length.out = nodes)
  s <- rep(s, length.out = nodes)
  if(length(nodes) < 2) {
    weights <- lapply(1:nodes, function(i) {
      sw <- switch(type,
        "sigmoid" = runif(1, log((1 - r[i])/r[i]), s[i] * log((1 - r[i])/r[i])),
        "gauss" = runif(1, sqrt(-log(r[i])), s[i] * sqrt(-log(r[i]))),
        "softplus" = runif(1, -log(exp(r[i]) - 1), s[i] * -log(exp(r[i]) - 1)),
        "cos" = runif(1, acos(r[i]), s[i] * acos(r[i])),
        "sin" = runif(1, acos(r[i]), s[i] * acos(r[i]))
      )
      w <- runif(k - 1, -1, 1)
      if(!is.null(dropout)) {
        j <- sample(c(FALSE, TRUE), size = length(w), replace = TRUE, prob = c(1 - dropout, dropout))
        w[j] <- 0
      }
      w <- w * sw / sum(w)
      if(length(w) < 2)
        w <- w * sample(c(-1, 1), size = 1)
      b <- -1 * (t(w) %*% (if(is.null(x)) runif(k - 1, 0, 1) else x[i, ]))
      w <- c(b, w)
      names(w) <- paste0("bw", i, "_w", 0:(k - 1))
      w
    })
  } else {
    weights <- list()
    for(i in 1:length(nodes)) {
      weights[[i]] <- lapply(1:nodes[i], function(ii) {
        if(i < 2) {
          w <- runif(k, -1, 1)
        } else {
          w <- runif(nodes[i - 1] + 1, -1, 1)
        }
        w[1] <- runif(1, 0, 1)
        names(w) <- paste0("bw", i, "_w", 0:(length(w) - 1))
        w
      })
    }
  }
  attr(weights, "type") <- type
  return(weights)
}

nnet2Zmat <- function(X, weights, afun)
{
  nc <- length(weights[[1]])
  Z <- list()
  k <- 1
  for(i in names(weights)) {
    afun <- switch(i,
      "relu" = function(x) {
        x[x < 0] <- 0
        x
      },
      "sigmoid" = function(x) {
        1 / (1 + exp2(-x))
      },
      "tanh" = tanh,
      "sin" = sin,
      "cos" = sin,
      "gauss" = function(x) { exp2(-x^2) },
      "identity" = function(x) { x },
      "softplus" = function(x) { log(1 + exp2(x)) }
    )
    for(j in 1:nc) {
      Z[[k]] <- afun(X %*% weights[[i]][[j]])
      k <- k + 1
    }
  }
  return(do.call("cbind", Z))
}

smooth.construct.nnet0.smooth.spec <- function(object, data, knots, ...)
{
  if(is.null(object$formula)) {
    object$formula <- as.formula(paste("~", paste(object$term, collapse = "+")))
    object$dim <- length(object$term)
    object$by <- "NA"
    object$type <- "single"
    object$xt$k <- object$bs.dim
  }

  oc <- object$xt$oc
  if(is.null(oc))
    oc <- FALSE

  object$xt[["standardize"]] <- object[["standardize01"]] <- object$xt[["standardize01"]] <- TRUE
  object <- smooth.construct.la.smooth.spec(object, data, knots)
  object[!(names(object) %in% c("formula", "term", "label", "dim", "X", "xt", "lasso"))] <- NULL
  nodes <- object$xt$k
  if(!is.null(object$xt$weights))
    nodes <- length(object$xt$weights)

  if(length(nodes) < 2) {
    if(nodes < 0)
      nodes <- 10
  }

  object$X <- cbind(1, object$X)

  if(is.null(object$xt$weights)) {
    object$xt$afun <- if(is.null(object$xt$afun)) "sigmoid" else object$xt$afun
    type <- object$xt$afun
  } else {
    type <- attr(object$xt$weights, "type")
    object$xt$afun <- type
  }

  if(is.null(object$xt$rint))
    object$xt$rint <- 0.1
  if(is.null(object$xt$sint))
    object$xt$sint <- c(10, 100)

  if(!is.list(object$xt$rint)) {
    object$xt$rint <- rep(list(object$xt$rint), length.out = length(type))
    names(object$xt$rint) <- type
  }
  if(!is.list(object$xt$sint)) {
    object$xt$sint <- rep(list(object$xt$sint), length.out = length(type))
    names(object$xt$sint) <- type
  }

  if(is.null(object$xt$weights)) {
    nobs <- nrow(object$X)
    object$xt[["tx"]] <- object$X[sample(1:nobs, size = nodes, replace = if(nodes >= nobs) TRUE else FALSE), -1, drop = FALSE]
    if(is.null(object$xt[["nobs"]]))
      object$xt[["nobs"]] <- 10
    if(is.null(object$xt[["maxit"]]))
      object$xt[["maxit"]] <- 100
    object$sample_weights <- function(x = NULL, y = NULL, weights = NULL, wts = NULL) {
      if(!is.null(y)) {
        if(length(unique(y)) < 50)
          y <- NULL
      }
      n.weights(nodes, ncol(object$X) - 1L, rint = object$xt$rint[[1]],
        sint = object$xt$sint[[1]], type = type[1],
        x = x, dropout = object$xt[["dropout"]], y = y, tntake = object$xt[["nobs"]],
        wm = object$xt[["wm"]], weights = weights, wts = wts, maxit = object$xt[["maxit"]])
    }
    object$n.weights <- object$sample_weights(object$xt[["tx"]])
  } else {
    if(length(object$xt$weights) != nodes)
      stop("not enough weights supplied!")
    object$n.weights <- object$xt$weights
  }

  object$S <- list()

  df <- nodes
  const <- object$xt$const
  if(is.null(const))
    const <- 1e-05
  ## pt <- object$xt$pt
  object$xt$pt <- pt <- "ridge"
  if(is.null(pt))
    pt <- "ridge"

  alpha <- object$xt$alpha
  if(is.null(alpha))
    alpha <- 0.5

  k <- 1
  if("lasso" %in% pt) {
    object$S[[1]] <- function(parameters, ...) {
      b <- parameters[1:nodes]
      A <- df / sqrt(b^2 + const)
      A <- if(length(A) < 2) matrix(A, 1, 1) else diag(A)
      A
    }
    attr(object$S[[k]], "npar") <- ncol(object$X)
    k <- k + 1
  }
  if(("enet" %in% pt) | ("elasticnet" %in% pt)) {
    object$S[[1]] <- function(parameters, ...) {
      b <- parameters[1:nodes]
      A <- df / sqrt(b^2 + const)
      A <- if(length(A) < 2) matrix(A, 1, 1) else diag(A)
      A <- alpha * A + diag(1 - alpha, ncol(A))
      A
    }
    attr(object$S[[k]], "npar") <- ncol(object$X)
    k <- k + 1
  }
  if("ridge" %in% pt) {
    A <- diag(1, nodes)
    object$S[[k]] <- A
  }

  object$xt$center <- if(is.null(object$xt$center)) FALSE else object$xt$center
  object$by <- "NA"
  object$null.space.dim <- 0
  object$bs.dim <- ncol(object$X)

  object$rank <- df
  object$xt$binning <- FALSE

  object$update <- object$xt$update <- nnet0_update
  object$boost.fit <- function(...) stop("Boost not yet implemented for n()!")
  object$propose <- GMCMC_iwls ## nnet0_propose

  object$activ_fun <- switch(type[1],
    "relu" = function(x) {
      x[x < 0] <- 0
      return(x - mean(x))
    },
    "sigmoid" = function(x) {
      f <- 1 / (1 + exp2(-x))
      return(f - mean(f))
    },
    "tanh" = function(x) {
      f <- tanh(x)
      return(f - mean(f))
    },
    "sin" = function(x) {
      f <- sin(x)
      return(f - mean(f))
    },
    "cos" = function(x) {
      f <- cos(x)
      return(f - mean(f))
    },
    "gauss" = function(x) {
      f <- exp2(-x^2)
      return(f - mean(f))
    },
    "softplus" = function(x) {
      f <- log(1 + exp2(x))
      return(f - mean(f))
    }
  )

  object$activ_grad <- switch(type[1],
    "relu" = function(x) {
      x[x > 0] <- 1
      x[x <= 0] <- 0
      return(x)
    },
    "sigmoid" = function(x) {
      1 / (1 + exp2(-x)) * (1 - 1 / (1 + exp2(-x)))
    },
    "tanh" = function(x) { 1 - tanh(x)^2 },
    "sin" = cos,
    "cos" = sin,
    "gauss" = function(x) -(exp2(-x^2) * (2 * x)),
    "softplus" = function(x) exp2(x)/(1 + exp2(x))
  )

  object$activ_hess <- switch(type[1],
    "relu" = function(x) {
      return(rep(0, length(x)))
    },
    "sigmoid" = function(x) {
      exp2(-x)/(1 + exp2(-x))^2 * (1 - 1/(1 + exp2(-x))) - 1/(1 + exp2(-x)) * (exp2(-x)/(1 + exp2(-x))^2)
    },
    "tanh" = function(x) { -(2 * (1/cosh(x)^2 * tanh(x))) },
    "sin" = sin,
    "cos" = cos,
    "gauss" = function(x) -(exp2(-x^2) * 2 - exp2(-x^2) * (2 * x) * (2 * x)),
    "softplus" = function(x) exp2(x)/(1 + exp2(x)) - exp2(x) * exp2(x)/(1 + exp2(x))^2
  )

  nc <- ncol(object$X) - 1L

  if(oc) {
    ## Null( cbind(m, Null(universe)) )
    object$oc <- list()
    for(j in object$term) {
      if(is.numeric(data[[j]])) {
        sm <- eval(parse(text = paste0("ti(", j, ",k=20)")))
        object$oc[[j]] <- smoothCon(sm, data, knots)[[1]]
      }
    }
    OC <- list()
    for(j in object$term) {
      if(is.numeric(data[[j]]))
        OC[[j]] <- PredictMat(object$oc[[j]], data)
    }
    OC <- do.call("cbind", OC)
    OC <- tcrossprod(qr.Q(qr(OC)))
    attr(object$X, "oc") <- OC
    ##object$X <- object$X - object$smC %*% object$X
  }

  object$fit.fun <- function(X, b, ...) {
    nb <- names(b)
    if(is.null(nb) | (ncol(X) == nc)) {
      fit <- drop(X %*% b)
    } else {
      nb <- strsplit(nb, ".", fixed = TRUE)
      nb <- sapply(nb, function(x) { x[length(x)] })
      names(b) <- nb
      fit <- 0
      for(j in 1:nodes) {
        z <- drop(X %*% b[paste0("bw", j, "_w", 0:nc)])
        fit <- fit + b[paste0("bb", j)] * object$activ_fun(z)
      }
    }
    fit <- fit - mean(fit)
    return(fit)
  }

  object$getZ <- function(X, b, ...) {
    Z <- matrix(0, nrow(X), nodes)
    for(j in 1:nodes) {
      z <- drop(X %*% b[paste0("bw", j, "_w", 0:nc)])
      Z[, j] <- object$activ_fun(z)
    }
    if(!is.null(attr(X, "oc"))) {
      Z <- Z - attr(X, "oc") %*% Z
    }
    return(Z)
  }

  ##attr(object$fit.fun, ".internal") <- TRUE

  object$state <- list()
  tau2 <- rep(0.1, length(object$S))
  names(tau2) <- paste0("tau2", 1:length(object$S))
  object$state$parameters <- c(rep(0, nodes))
  names(object$state$parameters) <- paste0("bb", 1:nodes)
  object$state$fitted.values <- rep(0, nrow(object$X))
  object$nodes <- nodes
  class(object) <- c(class(object), "special")

  bw <- runif(length(unlist(object$n.weights)), -1e-05, 1e-05)
  bw <- rep(1, length(unlist(object$n.weights)))
  names(bw) <- names(unlist(object$n.weights))
  object$state$parameters <- c(object$state$parameters, bw, tau2)
  object$state$edf <- 0
  object$xt$prior <- "hc"
  object[["a"]] <- 1e-04
  object[["b"]] <- 1e-04
  object$rank <- qr(object$S[[1]])$rank
  pf <- make.prior(object)
  object$prior <- pf$prior
  object$grad <- pf$grad
  object$hess <- pf$hess
  object$PredictMat <- Predict.matrix.nnet0.smooth
  object$fxsp <- object$xt$fx
  if(is.null(object$fxsp))
    object$fxsp <- FALSE
  object$nobs <- nrow(object$X)
  object$weights <- rep(0, length = object$nobs)
  object$rres <- rep(0, length = object$nobs)
  object$fit.reduced <- rep(0, length = object$nobs)
  object$binning <- match.index(matrix(1:nrow(object$X), ncol = 1))
  object$binning$order <- order(object$binning$match.index)
  object$binning$sorted.index <- object$binning$match.index[object$binning$order]
  if(is.null(object$xt$decay))
    object$xt$decay <- 0.0001

  class(object) <- c("nnet0.smooth", "no.mgcv", "special")

  object
}


nnet0_update <- function(x, family, y, eta, id, weights, criterion, ...)
{
  args <- list(...)
  nobs <- length(eta[[1L]])

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

  par <- x$state$parameters
  Z <- x$getZ(x$X, par)

  nc <- ncol(x$X)

  lls <- family$loglik(y, family$map2par(eta))

  eta[[id]] <- eta[[id]] - fitted(x$state)
  e <- z - eta[[id]]

  I <- diag(c(0, rep(x$xt$decay, ncol(x$X)))) 

  gradfun <- function(w, i, j, fit) {
    Xw <- drop(x$X %*% w[-1L])
    Z[, j] <- x$activ_fun(Xw)
    fit <- fit + Z[, j] * w[1L]
    fit <- fit - mean(fit)
    eta[[id]] <- eta[[id]] + fit
    score <- family$score[[id]](y, family$map2par(eta))
    gr <- score * cbind(Z[, j], w[1L] * x$activ_grad(Xw) * x$X)
    gr <- colSums(gr) - drop(I %*% w)
    return(-gr)
  }

#  hessfun <- function(w, i, j, fit) {
#    Xw <- drop(x$X %*% w[-1L])
#    Z[, j] <- x$activ_fun(Xw)
#    fit <- fit + Z[, j] * w[1L]
#    eta[[id]] <- eta[[id]] + fit
#    score <- family$score[[id]](y, family$map2par(eta))
#    hess <- -1 * family$hess[[id]](y, family$map2par(eta))

#    h0 <- sum(Z[,j]^2 * hess)

#    dh <- w[1L]^2 * x$activ_grad(Xw)^2 * hess +  w[1L] * x$activ_hess(Xw) * score
#    h1 <- crossprod(x$X * dh, x$X)

#    h2 <- x$X * w[1L] * Z[, j] * x$activ_grad(Xw) * hess +
#      x$X * x$activ_grad(Xw) * score
#    h2 <- colSums(h2)

#    h3 <- rbind(c(h0, h2), cbind(h2, h1)) + I

#    return(-h3)    
#  }

  objfun <- function(w, i, j, fit) {
    Z[, j] <- x$activ_fun(drop(x$X %*% w[-1L]))
    fit <- fit + Z[, j] * w[1L]
    fit <- fit - mean(fit)
    eta[[id]] <- eta[[id]] + fit
    ll <- family$loglik(y, family$map2par(eta)) - t(w) %*% I %*% w

    ## attr(ll, "gradient") <- gradfun(w, i, j, fit)
    ## attr(ll, "hessian") <- hessfun(w, i, j, fit)

    return(-ll)
  }

  tau2 <- get.par(par, "tau2")
#  ZWZ <- crossprod(Z * hess, Z)
#  P <- matrix_inv(ZWZ + 1/tau2 * x$S[[1]])
#  par[1:x$nodes] <- drop(P %*% crossprod(Z * hess, e))

  fit <- drop(Z %*% par[1:x$nodes])
  fit <- fit - mean(fit)

  eta2 <- eta

  for(j in 1:x$nodes) {
    i <- paste0("bw", j, "_w", 0:(nc - 1))

    eta2[[id]] <- eta[[id]] + fit
    ll0 <- family$loglik(y, family$map2par(eta2))
    fit <- fit - Z[, j] * par[j]

    opt <- try(optim(c(par[j], par[i]), fn = objfun, gr = gradfun,
      method = "BFGS", i = i, j = j, fit = fit), silent = TRUE)

#    H <- matrix_inv(hessfun(c(par[j], par[i]), i = i, j = j, fit = fit))
#    S <- gradfun(c(par[j], par[i]), i = i, j = j, fit = fit)
#    step <- try(drop(H %*% S), silent = TRUE)

#    if(inherits(step, "try-error")) {
#      fit <- fit + Z[, j] * par[j]
#      next
#    }

#    g <- c(par[j], par[i]) - 0.000001 * S

#    opt <- list(
#      "value" = objfun(g, i = i, j = j, fit = fit),
#      "par" = g
#    )

#    opt <- try(nlm(f = objfun, p = c(par[j], par[i]), i = i, j = j, fit = fit,
#      check.analyticals = FALSE, hessian = TRUE), silent = TRUE)

#    opt <- list("par" = opt$estimate, "value" = opt$minimum)

#w <- opt$par
#print(gradfun(w, i, j, fit))
#print(numericDeriv(quote(objfun(w, i, j, fit)), "w"))
#stop()

    if(!inherits(opt, "try-error")) {
      if(((-1 * opt$value) > ll0) & (-1 * opt$value > lls)) {
        par[i] <- opt$par[-1L]
        par[j] <- opt$par[1L]
        Z[, j] <- x$activ_fun(drop(x$X %*% opt$par[-1L]))
      }
    }

    fit <- fit + Z[, j] * par[j]
    fit <- fit - mean(fit)
  }

  if(!is.null(attr(x$X, "oc")))
    Z <- Z - attr(x$X, "oc") %*% Z

  ZWZ <- crossprod(Z * hess, Z)
  edf0 <- args$edf - x$state$edf

  objfun2 <- function(tau2) {
    P <- matrix_inv(ZWZ + 1/tau2 * x$S[[1]])
    edf <- sum_diag(ZWZ %*% P) + nc * x$nodes
    g <- P %*% crossprod(Z * hess, e)
    fit <- drop(Z %*% g)
    fit <- fit - mean(fit)
    eta[[id]] <- eta[[id]] + drop(Z %*% g)
    ic <- get.ic(family, y, family$map2par(eta), edf0 + edf, nobs, type = "BIC")
    return(ic)
  }

  ic0 <- objfun2(tau2)

  tau22 <- tau2.optim(objfun2, start = tau2)

  ic1 <- objfun2(tau2)

  if(ic1 <= ic0)
    tau2 <- tau22

  P <- matrix_inv(ZWZ + 1/tau2 * x$S[[1]])
  par[1:x$nodes] <- drop(P %*% crossprod(Z * hess, e))
  par["tau21"] <- tau2

  fit <- drop(Z %*% par[1:x$nodes])
  fit <- fit - mean(fit)
  x$state$fitted.values <- fit
  x$state$parameters <- par
  x$state$edf <- sum_diag(ZWZ %*% P)
  x$state$log.prior <- x$prior(par)

#  eta[[id]] <- eta[[id]] + fit
#  lls2 <- family$loglik(y, family$map2par(eta))

#  cat("\n---\n")
#  print(c(lls, lls2))

  return(x$state)
}

nnet0_propose <- function(family, theta, id, eta, y, data, weights = NULL, offset = NULL, ...)
{
  theta <- theta[[id[1]]][[id[2]]]
  id <- id[1]

  if(!is.null(offset)) {
    for(j in names(offset))
      eta[[j]] <- eta[[j]] + offset[[j]]
  }

  if(is.null(attr(theta, "fitted.values")))
    attr(theta, "fitted.values") <- data$fit.fun(data$X, theta)

  tau2 <- 1/data$xt$decay
  I <- diag(c(1e-10, rep(1/tau2, ncol(data$X))))

  gradfun <- function(w, i, j, eta, id) {
    Xw <- drop(data$X %*% w[-1L])
    Zj <- data$activ_fun(Xw)
    eta[[id]] <- eta[[id]] + Zj * w[1L]
    score <- family$score[[id]](y, family$map2par(eta))
    gr <- score * cbind(Zj, w[1L] * data$activ_grad(Xw) * data$X)
    gr <- colSums(gr) - drop(I %*% w)
    return(gr)
  }

  hessfun <- function(w, i, j, eta, id) {
    Xw <- drop(data$X %*% w[-1L])
    Zj <- data$activ_fun(Xw)
    eta[[id]] <- eta[[id]] + Zj * w[1L]
    score <- family$score[[id]](y, family$map2par(eta))
    hess <- family$hess[[id]](y, family$map2par(eta))

    h0 <- sum(Zj^2 * hess)

    dh <- w[1L]^2 * data$activ_grad(Xw)^2 * hess +  w[1L] * data$activ_hess(Xw) * score
    h1 <- crossprod(data$X * dh, data$X)

    h2 <- data$X * w[1L] * Zj * data$activ_grad(Xw) * hess +
      data$X * data$activ_grad(Xw) * score
    h2 <- colSums(h2)

    h3 <- rbind(c(h0, h2), cbind(h2, h1)) + I

    return(-h3)    
  }

  nc <- ncol(data$X)

  pf <- function(gamma) {
    tau2 <- 1/data$xt$decay
    a <- b <- 1e-04
    igs <- log((b^a)) - log(gamma(a))
    lp <- -log(tau2) * (ncol(data$X) + 1)/2 + drop(-0.5/tau2 * 
      t(gamma) %*% I %*% gamma) + (igs + (-a - 1) * log(tau2) - b/tau2)
    lp
  }

  fit <- data$fit.fun(data$X, theta)

  for(j in 1:data$nodes) {
    i <- paste0("bw", j, "_w", 0:(nc - 1))

    pibeta <- family$loglik(y, family$map2par(eta))
    p1 <- pf(c(theta[j], theta[i]))

    fj0 <- data$activ_fun(data$X %*% theta[i]) * theta[j]

#plot(y ~ x, data = d)
#plot2d(eta[[id]] ~ d$x, add = TRUE, col.lines = 4, lwd = 5)

    eta[[id]] <- eta[[id]] - fj0

    g <- gradfun(c(theta[j], theta[i]), i = i, j = j, eta = eta, id = id)
    h <- hessfun(c(theta[j], theta[i]), i = i, j = j, eta = eta, id = id)

    S <- matrix_inv(-1 * h)

    m <- drop(c(theta[j], theta[i]) + S %*% g)

    theta2 <- drop(rmvnorm(n = 1, mean = m, sigma = S))

    p2 <- pf(theta2)
    qbetaprop <- dmvnorm(matrix(theta2, nrow = 1), mean = m, sigma = S, log = TRUE)

    g2 <- gradfun(theta2, i = i, j = j, eta = eta, id = id)
    h2 <- hessfun(theta2, i = i, j = j, eta = eta, id = id)

    S2 <- matrix_inv(-1 * h2)

    m2 <- drop(theta2 + S2 %*% g2)

    Xw <- drop(data$X %*% theta2[-1L])
    fj1 <- data$activ_fun(Xw) * theta2[1L]
    eta[[id]] <- eta[[id]] + fj1

#plot2d(eta[[id]] ~ d$x, add = TRUE, col.lines = 3, lwd = 3)
#Sys.sleep(1)

    pibetaprop <- family$loglik(y, family$map2par(eta))

    qbeta <- dmvnorm(matrix(c(theta[j], theta[i]), nrow = 1), mean = m2, sigma = S2, log = TRUE)

    alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
    if(!is.finite(alpha))
      alpha <- NA
    accepted <- if(is.na(alpha)) FALSE else log(runif(1)) <= alpha

#print(accepted)
#print(pibeta)
#print(pibetaprop)
#cat("**")
#print(qbeta)
#print(qbetaprop)
#cat("---\n")

    if(accepted) {
      theta[j] <- theta2[1L]
      theta[i] <- theta2[-1L]
      fit <- fit - fj0 + fj1
    } else {
      eta[[id]] <- (eta[[id]] - fj1) + fj0
    }
  }

  fit <- fit - mean(fit)
  attr(theta, "fitted.values") <- fit

  return(list("parameters" = theta, "alpha" = log(1),
    "extra" = c("edf" = data$nodes * ncol(data$X))))
}

nnet00_update <- function(x, family, y, eta, id, weights, criterion, ...)
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

  ll0 <- family$loglik(y, family$map2par(eta))

  ## Compute partial predictor.
  eta[[id]] <- eta[[id]] - fitted(x$state)

  ## Compute reduced residuals.
  e <- z - eta[[id]]

  par <- x$state$parameters

  objfun <- function(w, i, j, ...) {
    par[i] <- w
    Z <- x$getZ(x$X, par)
    ZWZ <- crossprod(Z * hess, Z)
    P <- matrix_inv(ZWZ)
    beta <- drop(P %*% crossprod(Z * hess, e))
    eta[[id]] <- eta[[id]] + drop(Z %*% beta)
    ll <- -1 * family$loglik(y, family$map2par(eta))

    score <- family$score[[id]](y, family$map2par(eta))
    gr <- score * beta[j] * x$activ_grad(drop(x$X %*% w)) * x$X
    attr(ll, "gradient") <- -colSums(gr)

#    hess <- family$hess[[id]](y, family$map2par(eta))
#    Xw <- drop(x$X %*% w)
#    h <- matrix(0, nc, nc)
#    for(l in 1:nobs) {
#      for(ii in 1:nc) {
#        for(jj in 1:nc) {
#          if(ii <= jj) {
#            h[ii, jj] <- h[ii, jj] + x$X[l, ii] * x$X[l, jj] * x$activ_hess(Xw[l]) * score[l] +
#              x$X[l, ii] * x$X[l, jj] * x$activ_grad(Xw[l])^2 * hess[l]
#            h[jj, ii] <- h[ii, jj]
#          }
#        }
#      }
#    }
#    #attr(ll, "hessian") <- h

#    h2 <- t(x$X) %*% diag(x$activ_hess(Xw) * score + x$activ_grad(Xw[l])^2 * hess) %*% x$X

    return(ll)
  }

  nc <- ncol(x$X)

  gradfun <- function(w, i, j, ...) {
    par[i] <- w
    Z <- x$getZ(x$X, par)
    ZWZ <- crossprod(Z * hess, Z)
    P <- matrix_inv(ZWZ)
    beta <- drop(P %*% crossprod(Z * hess, e))
    eta[[id]] <- eta[[id]] + drop(Z %*% beta)
    score <- family$score[[id]](y, family$map2par(eta))
    gr <- score * beta[j] * x$activ_grad(drop(x$X %*% w)) * x$X
    return(-colSums(gr))
  }

  hessfun <- function(w, i, j) {
    par[i] <- w
    Z <- x$getZ(x$X, par)
    ZWZ <- crossprod(Z * hess, Z)
    P <- matrix_inv(ZWZ)
    beta <- drop(P %*% crossprod(Z * hess, e))
    eta[[id]] <- eta[[id]] + drop(Z %*% beta)
    score <- family$score[[id]](y, family$map2par(eta))
    hess <- family$hess[[id]](y, family$map2par(eta))

    Xw <- drop(x$X %*% w)
    h <- matrix(0, nc, nc)
    for(l in 1:nobs) {
      for(ii in 1:nc) {
        for(jj in 1:nc) {
          if(ii <= jj) {
            h[ii, jj] <- h[ii, jj] + x$X[l, ii] * x$X[l, jj] * x$activ_hess(Xw[l]) * score[l] +
              x$X[l, ii] * x$X[l, jj] * x$activ_grad(Xw[l])^2 * hess[l]
            h[jj, ii] <- h[ii, jj]
          }
        }
      }
    }
    return(h)    
  }

  for(j in 1:x$nodes) {
    i <- paste0("bw", j, "_w", 0:(nc - 1))

#w <- par[i]
#a <- numericDeriv(quote(objfun(w, i = i)), "w")
#b <- gradfun(w, i = i, j = j)
#print(attr(a, "gradient"))
#print(b)
#stop()

    opt <- optim(par[i], fn = objfun, gr = gradfun,
      method = "L-BFGS-B", i = i, j = j)

#    opt <- try(nlm(f = objfun, p = par[i], i = i, j = j,
#      check.analyticals = TRUE, hessian = TRUE), silent = TRUE)

#oo <- optimHess(opt$estimate, objfun, gr = gradfun, i = i, j = j)

#print(oo)

#print(opt)
#cat("---\n")
##print(gradfun(opt$estimate, i, j))
#print(hessfun(opt$estimate, i, j))
#stop()

    if(inherits(opt, "try-error"))
      next

    opt <- list("par" = opt$estimate, "value" = opt$minimum)

    if((-1* opt$value) > ll0) {
      par[i] <- opt$par
      Z <- x$getZ(x$X, par)
      ZWZ <- crossprod(Z * hess, Z)
      P <- matrix_inv(ZWZ)
      beta <- drop(P %*% crossprod(Z * hess, e))
      par[grep("bb", names(par))] <- beta
    }
  }

  Z <- x$getZ(x$X, par)
  ZWZ <- crossprod(Z * hess, Z)
  P <- matrix_inv(ZWZ)
  beta <- drop(P %*% crossprod(Z * hess, e))
  fit <- drop(Z %*% beta)
  fit <- fit - mean(fit)

  par[grep("bb", names(par))] <- beta

  x$state$fitted.values <- fit
  x$state$parameters <- par
  x$state$edf <- sum_diag(ZWZ %*% P)
  x$state$log.prior <- sum(dnorm(beta, sd = 1000, log = TRUE))

  return(x$state)
}


smooth.construct.nnet2.smooth.spec <- function(object, data, knots, ...)
{
  if(is.null(object$formula)) {
    object$formula <- as.formula(paste("~", paste(object$term, collapse = "+")))
    object$dim <- length(object$term)
    object$by <- "NA"
    object$type <- "single"
    object$xt$fx <- FALSE
    object$xt$k <- object$bs.dim
  }
  object$xt[["standardize"]] <- object[["standardize01"]] <- object$xt[["standardize01"]] <-  TRUE
  object <- smooth.construct.la.smooth.spec(object, data, knots)
  object[!(names(object) %in% c("formula", "term", "label", "dim", "X", "xt", "lasso"))] <- NULL
  nodes <- object$xt$k
  if(!is.null(object$xt$weights))
    nodes <- length(object$xt$weights)
  npen <- if(is.null(object$xt$npen)) 1 else object$xt$npen
  dotake <- FALSE

  if(length(nodes) < 2) {
    if(nodes < 0)
      nodes <- 10
  }

  object$X <- cbind(1, object$X)

  if(is.null(object$xt$weights)) {
    object$xt$afun <- if(is.null(object$xt$afun)) "sigmoid" else object$xt$afun
    type <- object$xt$afun
  } else {
    type <- attr(object$xt$weights, "type")
    object$xt$afun <- type
  }

#    object$Zmat <- function(X, weights, afun) {
#      Z <- list()
#      n <- nrow(X)
#      for(i in 1:length(weights)) {
#        Z[[i]] <- matrix(0, nrow = n, ncol = length(weights[[i]]))
#        for(j in 1:length(weights[[i]])) {
#          if(i < 2) {
#            Z[[i]][, j] <- afun(X %*% weights[[i]][[j]])
#          } else {
#            Z[[i]][, j] <- afun(cbind(1, Z[[i - 1]]) %*% weights[[i]][[j]])
#          }
#        }
#      }
#      return(Z[[length(Z)]])
#    }

  if(is.null(object$xt$rint))
    object$xt$rint <- c(0.01, 0.2)
  if(is.null(object$xt$sint))
    object$xt$sint <- c(1.01, 10)

  if(!is.list(object$xt$rint)) {
    object$xt$rint <- rep(list(object$xt$rint), length.out = length(type))
    names(object$xt$rint) <- type
  }
  if(!is.list(object$xt$sint)) {
    object$xt$sint <- rep(list(object$xt$sint), length.out = length(type))
    names(object$xt$sint) <- type
  }

  if(is.null(object$xt$weights)) {
    nobs <- nrow(object$X)
    object$xt[["tx"]] <- object$X[sample(1:nobs, size = nodes, replace = if(nodes >= nobs) TRUE else FALSE), -1, drop = FALSE]
    object$n.weights <- list()
    for(j in type) {
      object$n.weights[[j]] <- n.weights(nodes, ncol(object$X) - 1L, rint = object$xt$rint[[j]],
        sint = object$xt$sint[[j]], type = j,
        x = object$xt[["tx"]], dropout = object$xt[["dropout"]])
    }
  } else {
    if(length(object$xt$weights) != nodes)
      stop("not enough weights supplied!")
    object$n.weights <- object$xt$weights
  }

  object$X <- nnet2Zmat(object$X, object$n.weights, object$xt$afun)

  object$Xkeep <- which(apply(object$X, 2, function(x) { abs(diff(range(x))) })  > 1e-05)
  object$X <- object$X[, object$Xkeep, drop = FALSE]
  ## object$Xsubset <- subset_features(object$X)
  ## object$X <- object$X[, object$Xsubset, drop = FALSE]

  ## Orthogonal complement.
  oc <- object$xt$orthc
  if(is.null(oc))
    oc <- TRUE
  if(oc) {
    object$xt$nocenter <- TRUE
    object$xt$center <- FALSE
    object$sm <- list()
    for(j in object$term) {
      if(!is.factor(data[[j]])) {
        k <- min(c(10, length(unique(data[[j]])) - 1))
        if(k > 7) {
          sj <- eval(parse(text = paste0("s(", j,")")))
          object$sm[[j]] <- smoothCon(sj, data, knots = NULL,
            absorb.cons = TRUE, sparse.cons = 0, scale.penalty = TRUE)[[1]]
        }
      }
    }
    sm <- do.call("cbind", lapply(object$sm, function(x) { x$X }))
#    for(j in seq_along(object$sm))
#      object$sm[[j]]$X <- matrix(0, nrow = 0, ncol = ncol(object$sm[[j]]$X))

#  qrL <- qr(L)
#  Q <- qr.Q(qrL)
#  XtXinvXt <- tcrossprod(Q)
#  Sorth <- S - XtXinvXt%*%S

    smL <- qr(sm)
    object$smC <- tcrossprod(qr.Q(smL))
    object$X <- object$X - object$smC %*% object$X
    ##object$X <- object$smC %*% object$X
  }
  if(is.null(object$xt$nocenter)) {
    object$QR <- qr.Q(qr(crossprod(object$X,
      rep(1, length = nrow(object$X)))), complete = TRUE)[, -1]
    object$X <- object$X %*% object$QR
  }

  if(ncol(object$X) < 1)
    stop("please check your n() specifications, no columns in the design matrix!")

  object$S <- list()

  df <- ncol(object$X)
  const <- object$xt$const
  if(is.null(const))
    const <- 1e-05
  pt <- object$xt$pt
  if(is.null(pt))
    pt <- "ridge"

  alpha <- object$xt$alpha
  if(is.null(alpha))
    alpha <- 0.5

  k <- 1
  if("lasso" %in% pt) {
    object$S[[1]] <- function(parameters, ...) {
      b <- get.par(parameters, "b")
      if(is.null(object$xt$nocenter))
        b <- c(0, b)
      A <- df / sqrt(b^2 + const)
      A <- if(length(A) < 2) matrix(A, 1, 1) else diag(A)
      if(is.null(object$xt$nocenter))
        A <- crossprod(object$QR, A) %*% object$QR
      A
    }
    attr(object$S[[k]], "npar") <- ncol(object$X)
    k <- k + 1
  }
  if(("enet" %in% pt) | ("elasticnet" %in% pt)) {
    object$S[[1]] <- function(parameters, ...) {
      b <- get.par(parameters, "b")
      if(is.null(object$xt$nocenter))
        b <- c(0, b)
      A <- df / sqrt(b^2 + const)
      A <- if(length(A) < 2) matrix(A, 1, 1) else diag(A)
      A <- alpha * A + diag(1 - alpha, ncol(A))
      if(is.null(object$xt$nocenter))
        A <- crossprod(object$QR, A) %*% object$QR
      A
    }
    attr(object$S[[k]], "npar") <- ncol(object$X)
    k <- k + 1
  }
  if("ridge" %in% pt) {
    if(is.null(object$xt$nocenter)) {
      A <- diag(1, ncol(object$X) + 1L)
      A <- crossprod(object$QR, A) %*% object$QR
    } else {
      A <- diag(1, ncol(object$X))
    }
    object$S[[k]] <- A
  }

  object$xt$center <- if(is.null(object$xt$center)) FALSE else object$xt$center
  object$by <- "NA"
  object$null.space.dim <- 0
  object$bs.dim <- ncol(object$X)

  object$rank <- df#qr(object$S[[1]](runif(df)))$rank

  object$xt$prior <- "hc"
  #object$fx <- object$xt$fx <- object$xt$fxsp <- object$fxsp <- FALSE

  if(is.null(object$xt$alpha))
    object$xt$alpha <- 1
  if(is.null(object$xt$nlambda))
    object$xt$nlambda <- 100
  if(is.null(object$xt$lambda.min.ratio))
    object$xt$lambda.min.ratio <- 1e-20
  if(is.null(object$xt$update))
    object$xt$update <- bfit_iwls
  if(!is.function(object$xt$update)) {
    if(object$xt$update == "lasso") {
      object$no.assign.df <- TRUE
      object$update <- if(!dotake) bfit_glmnet else bfit_iwls
      object$xt$update <- NULL
    }
  }

  if(is.null(object$xt$K))
    object$xt$K <- 1
  if(is.null(object$xt$single))
    object$xt$single <- FALSE
  if(is.null(object$xt$ndf))
    object$xt$ndf <- 4
  else
    object$xt$ndf <- ceiling(object$xt$ndf)
  if(object$xt$ndf < 2)
    object$xt$ndf <- NULL
  if(is.null(object$xt$frac))
    object$xt$frac <- 0.99

  object$xt[["df"]] <- object$xt[["ndf"]]

  if(object$xt$single) {
    ncX <- ncol(object$X)
    nrX <- nrow(object$X)
    if(is.null(object$xt$ndf) | TRUE) {
      object$N <- apply(object$X, 2, function(x) {
        return((1/crossprod(x)) %*% t(x))
      })
    }
    object$boost.fit <- function(x, y, nu, hatmatrix = FALSE, weights = NULL, nthreads = 1, ...) {
      ## process weights.
      if(!is.null(weights))
        stop("weights are not supported for n()!")

      g2 <- rep(0, ncX)

      if(!is.null(x$xt$ndf) & FALSE) {
#        i <- sample(1:nrX, size = ceiling(nrX * x$xt$frac), replace = FALSE)
#        X2 <- x$X[i, ]
#        y2 <- y[i]
#        bf <- forward_reg(X2, y2, n = x$xt$ndf, nu = nu)
#        if(!is.null(bf)) {
#          g2[bf$take] <- nu/x$xt$K * bf$coefficients[-1]
#          x$state$fitted.values <- drop(x$X %*% g2)
#          x$state$rss <- sum((y - x$state$fitted.values)^2)
#        } else {
#          x$state$fitted.values <- rep(0, length(y))
#          x$state$rss <- sum(y^2)
#        }
        j <- sample(1:ncol(x$X), size = x$xt$ndf)
        Z <- x$X[, j, drop = FALSE]
#        pen <- if(is.null(x$xt$pen)) {
#          1e-05
#        } else {
#          x$xt$pen
#        }
#        K <- diag(pen, ncol(Z))
#        b <- nu/x$xt$K * matrix_inv(crossprod(Z) + K) %*% t(Z) %*% y
#        b <- glmnet(Z, y, alpha = 0.5, family = "gaussian")
#        b <- nu/x$xt$K * b$beta[, ncol(b$beta)]
#        g2[j] <- b
#        x$state$fitted.values <- drop(Z %*% b)
#        x$state$rss <- sum((y - x$state$fitted.values)^2)
##print(sum_diag(crossprod(Z) %*% matrix_inv(crossprod(Z) + K)))
      } else {
        bf <- boost_fit_nnet(nu, x$X, x$N, y, x$binning$match.index, nthreads = nthreads)
        j <- which.min(bf$rss)
        g2[j] <- bf$g[j]
        x$state$fitted.values <- bf$fit[, j]
        x$state$rss <- sum((y - x$state$fitted.values)^2)
      }

      names(g2) <- paste0("b", 1:ncX)
  
      ## Finalize.
      x$state$parameters <- set.par(x$state$parameters, g2, "b")

      if(hatmatrix) {
        stop("not supported for n()!")
      }
  
      return(x$state)
    }
#    if(is.null(object$xt$nn.control))
#      object$xt$nn.control <- list(k = 50, size = min(c(floor(ncol(object$X) / 2), 50)))

#    object$boost.fit <- function(x, y, nu, hatmatrix = FALSE, weights = NULL, nthreads = 1, ...) {
#      ## process weights.
#      if(!is.null(weights))
#        stop("weights are not supported for n()!")

#      b <- nn.fit(x$X, y, k = x$xt$nn.control$k, size = x$xt$nn.control$size)
#  
#      ## Finalize.
#      x$state$parameters <- nu * set.par(x$state$parameters, coef(b)[-1], "b")
#      x$state$fitted.values <- drop(x$X %*% get.par(x$state$parameters, "b"))
#      x$state$rss <- sum((y - x$state$fitted.values)^2)

#      if(hatmatrix) {
#        stop("not supported for n()!")
#      }
#  
#      return(x$state)
#    }
  }

  object$fit.fun <- ff0
  attr(object$fit.fun, ".internal") <- TRUE

#plot2d(object$X ~ data$times, main = type)
#stop()
#Sys.sleep(5)
##stop()

  object$xt$binning <- FALSE

  if(!object$xt$single) {
    object$xt[["lambda"]] <- object$xt[["sp"]] <- object[["lambda"]] <- object[["sp"]] <- object$xt[["tau2"]] <- NULL
    object$xt$lambda.min.ratio <- NULL
    class(object) <- c("nnet3.smooth", "mgcv.smooth")
  } else {
    object[["lambda"]] <- object$xt[["lambda"]] <- NULL
    class(object) <- c("nnet2.smooth", "mgcv.smooth", "lasso.smooth")
  }
  object$fixed <- FALSE
  if(!is.null(object$xt$fx)) {
    if(object$xt$fx) {
      object$sp <- object$xt$sp <- 1e-10
      object$S <- list(diag(1, ncol(object$X)))
    }
  }

  object
}


nn.fit <- function(x, y, k = 50, size = min(c(floor(ncol(x) / 2), 50))) {
  my <- mean(y)
  y <- y - my
  coefs <- matrix(0, k, ncol(x))
  id <- 1:ncol(x)
  for(i in 1:k) {
    f <- 0
    for(j in sample(id, size = size, replace = FALSE)) {
      e <- y - f
      m <- lm.fit(x[, j, drop = FALSE], e)
      f <- f + m$fitted.values
      coefs[i, j] <- m$coefficients
    }
  }
  rval <- list(
    "coefficients" = c(my, apply(coefs, 2, mean))
  )
  rval$fitted.values <- drop(cbind(1, x) %*% rval$coefficients)
  return(rval)
}


subset_features <- function(x, eps = 0.01)
{
  cols <- 1:ncol(x)
  drop <- NULL
  err <- eps - 1
  while(err < eps) {
    rss <- NULL
    for(j in cols) {
      xs <- x[, if(is.null(drop)) -j else c(-drop, -j), drop = FALSE]
      P <- xs %*% matrix_inv(crossprod(xs) + diag(0.0001, ncol(xs))) %*% t(xs)
      rss <- c(rss, sum((x - P %*% x)^2))
    }
    drop <- c(drop, cols[which.min(rss)])
    cols <- cols[!(cols %in% drop)]
    err <- min(rss)
  }
  c(1:ncol(x))[-drop]
}


forward_reg <- function(x, y, n = 4, ...)
{
  k <- 0
  cols <- 1:ncol(x)
  take <- NULL
  while(k < n) {
    rss <- NULL
    for(j in cols)
      rss <- c(rss, sum(lm.fit(cbind(1, x[, c(take, j), drop = FALSE]), y)$residuals^2) + 2 * length(c(take, j)))
    take <- c(take, cols[which.min(rss)])
    cols <- cols[!(cols %in% take)]
    k <- k + 1
  }
  return(c(list("take" = take), lm.fit(cbind(1, x[, take, drop = FALSE]), y)))
}


forward_reg2 <- function(x, y, n = 4, nu, maxit = 100, ...)
{
  r <- y - mean(y)
  k <- 0
  g <- rep(0, ncol(x))
  l <- 0
  while(k < n) {
    cxr <- drop(cor(x, r))
    j <- which.max(abs(cxr))
    delta <- nu * sign(cxr[j])
    g[j] <- g[j] + delta
    r <- r - delta * x[, j]
    k <- sum(g != 0)
    l <- l + 1
    if(l > maxit)
      break
  }
  j <- which(g != 0)
  if(length(j))
    return(c(list("take" = j), lm.fit(x[, j, drop = FALSE], y)))
  else
    return(NULL)
}


#tanh2 <- function(x) {
#  y <- rep(0, length(x))
#  y[x > 1.92033] <- 0.96016
#  i <- x > 0 & x <= 1.92033
#  y[i] <- 0.96016 - 0.26037 * (x[i] - 1.92033)^2
#  i <- x > -1.92033 & x < 0
#  y[i] <- 0.26037 * (x[i] + 1.92033)^2 - 0.96016
#  y[x <= -1.92033] <- -0.96016
#  y
#}

#curve(tanh2, -3, 3)
#curve(tanh, -3, 3, col = 2, add = TRUE)


Predict.matrix.nnet0.smooth <- function(object, data)
{
  object[["standardize"]] <- standardize <- if(is.null(object$xt[["standardize"]])) TRUE else object$xt[["standardize"]]
  object[["standardize01"]] <- if(standardize) TRUE else FALSE
  X <- cbind(1, Predict.matrix.lasso.smooth(object, data))

  if(!is.null(object$oc)) {
    OC <- list()
    for(j in object$term) {
      if(is.numeric(data[[j]]))
        OC[[j]] <- PredictMat(object$oc[[j]], data)
    }
    OC <- do.call("cbind", OC)
    OC <- tcrossprod(qr.Q(qr(OC)))
    attr(X, "oc") <- OC
  }

  return(X)
}


Predict.matrix.nnet2.smooth <- function(object, data)
{
  object[["standardize"]] <- standardize <- if(is.null(object$xt[["standardize"]])) TRUE else object$xt[["standardize"]]
  object[["standardize01"]] <- if(standardize) TRUE else FALSE
  X <- nnet2Zmat(cbind(1, Predict.matrix.lasso.smooth(object, data)), object$n.weights, object$xt$afun)
  X <- X[, object$Xkeep, drop = FALSE]
  if(!is.null(object$Xsubset))
    X <- X[, object$Xsubset, drop = FALSE]
  if(!is.null(object$smC)) {
    smX <- list()
    for(j in seq_along(object[["sm"]])) {
      smX[[j]] <- PredictMat(object[["sm"]][[j]], as.data.frame(data))
    }
    smX <- do.call("cbind", smX)
    smL <- qr(smX)
    smC <- tcrossprod(qr.Q(smL))
    X <- X - smC %*% X
    ## X <- smC %*% X
  }
  if(is.null(object$xt$nocenter)) {
#    for(j in 1:length(object$xt$cmeans))
#      X[, j] <- X[, j] - object$xt$cmeans[j]
    X <- X %*% object$QR
  }
  return(X)
}

Predict.matrix.nnet3.smooth <- Predict.matrix.nnet2.smooth


if(FALSE) {
  n <- 600

  d <- data.frame(
    "x1" = runif(n, -3, 3),
    "x2" = runif(n,-3, 3)
  )
  d$y <- d$x1 + d$x2^2 + sin(d$x1)*cos(d$x2) + rnorm(n, sd = 0.3)

  i <- sample(1:2, size = n, replace = TRUE)

  dtrain <- d[i == 1, ]
  dtest <- d[i == 2, ]

  b0 <- bamlss(y ~ s(x1)+s(x2), data= dtrain)
  b1 <- bamlss(y ~ s(x1)+s(x2)+n(~x1+x2,k=100,orthc=F,rint=0.1,sint=1000), data = dtest)

  p0 <- predict(b0, newdata = dtest, model = "mu")
  p1 <- predict(b1, newdata = dtest, model = "mu")

  plot(dtest$y ~ p0)
  points(p1, dtest$y, col = 2)
  abline(0, 1)

  mse <- c(
    "GAM" = mean((dtest$y - p0)^2),
    "GAM+NET" = mean((dtest$y - p1)^2)
  )

  print(mse)
}


nnet.fit <- function(X, y, nodes = 20, ..., random = FALSE, w = NULL, lambda = 0.001,
  optim = FALSE, maxit = 100, nu = 0.1)
{
  nc <- ncol(X)
  if(is.null(w))
    w <- n.weights(nodes, k = nc, type = "sigmoid", x = X, ...)
  w0 <- w
  par <- unlist(w)
  nw <- names(par)

  X <- cbind(1, X)

  Z <- matrix(0, nrow = nrow(X), ncol = nodes)
  for(j in 1:nodes)
    Z[, j] <- 1 / (1 + exp(-1 * drop(X %*% par[paste0("bw", j, "_w", 0:nc)])))
  par <- drop(matrix_inv(crossprod(Z) + diag(lambda, nodes)) %*% t(Z) %*% y)

  par <- runif(nodes, -0.1, 0.1)
  names(par) <- paste0("bb", 1:nodes)
  par <- c(par, unlist(w))

  ffn <- function(X, par) {
    fit <- 0
    for(j in 1:nodes) {
      z <- drop(X %*% par[paste0("bw", j, "_w", 0:nc)])
      fit <- fit + par[paste0("bb", j)] / (1 + exp(-z))
    }
    return(fit)
  }

  ffm <- function(X, par) {
    nm <- matrix(0, nrow = nrow(X), ncol = nodes)
    for(j in 1:nodes) {
      z <- drop(X %*% par[paste0("bw", j, "_w", 0:nc)])
      nm[, j] <-  1 / (1 + exp(-z))
    }
    return(nm)
  }

  if(random) {
    rval <- list(
      "fitted.values" = ffn(X, par),
      "coefficients" = par,
      "nodes" = nodes,
      "converged" = TRUE
    )
    class(rval) <- "nnet.fit"
    return(rval)
  }

  gradfun <- function(par, X, y) {
    Z <- ffm(X, par)
    beta <- chol2inv(chol(crossprod(Z) + diag(nodes) * lambda)) %*% t(Z) %*% y
    fit <- drop(Z %*% beta)
    s1 <- -(2 * (y - fit))
    gr2 <- matrix(0, nrow = nrow(X), ncol = nodes * (nc + 1))
    k <- 1
    for(j in 1:nodes) {
      s2 <- beta[j] * Z[, j] * (1 - Z[, j])
      for(i in 1:(nc + 1)) {
        if(i < 2) {
          gr2[, k] <- s1 * s2
        } else {
          gr2[, k] <- s1 * s2 * X[, i]
        }
        k <- k + 1
      }
    }
    return(colSums(gr2))
  }

  objfun <- function(par, X, y) {
    Z <- ffm(X, par)
    beta <- chol2inv(chol(crossprod(Z) + diag(nodes) * lambda)) %*% t(Z) %*% y
    fit <- Z %*% beta
    return(sum((y - fit)^2))
  }

  i <- grep("bw", names(par))

  if(optim) {
    opt <- optim(par = par[i], fn = objfun, gr = gradfun,
      method = "BFGS", X = X, y = y)
    par[i] <- opt$par
    Z <- ffm(X, opt$par)
    beta <- chol2inv(chol(crossprod(Z) + diag(nodes) * lambda)) %*% t(Z) %*% y
    fit <- drop(Z %*% beta)
    par[-i] <- beta
  } else {
    iter <- 1
    while(iter < maxit) {
      Z <- ffm(X, par)
      beta <- chol2inv(chol(crossprod(Z) + diag(nodes) * lambda)) %*% t(Z) %*% y
      fit <- Z %*% beta
      par[-i] <- beta
      grad <- gradfun(par, X, y)
      par[i] <- par[i] - nu * grad
      err <- sum((y - fit)^2)
      cat("iter", iter, "error", round(err, 4), "\n")
      iter <- iter + 1
    }
  }

  rval <- list(
    "fitted.values" = fit,
    "coefficients" = par[-i],
    "weights" = par[i],
    "initial" = unlist(w0),
    "nodes" = nodes,
    "converged" = if(optim) opt$convergence == 1L else NA
  )
  class(rval) <- "nnet.fit"

  return(rval)
}

predict.nnet.fit <- function(object, newX, ...)
{
  nc <- ncol(newX)
  newX <- cbind(1, newX)
  fit <- 0
  for(j in 1:object$nodes) {
    z <- newX %*% object$coefficients[paste0("bw", j, "_w", 0:nc)]
    fit <- fit + object$coefficients[paste0("bb", j)] * 1 / (1 + exp(-z))
  }
  return(fit)
}


nnet.fit.fun <- function(X, b, ...) {
  fit <- 0
  nc <- ncol(X)
  X <- cbind(1, X)
  nodes <- sum(grepl("bb", names(b)))
  if(any(grepl(".", nb <- names(b), fixed = TRUE))) {
    nb <- strsplit(nb, ".", fixed = TRUE)
    nb <- sapply(nb, function(x) { x[length(x)] })
    names(b) <- nb
  }
  for(j in 1:nodes) {
    z <- drop(X %*% b[paste0("bw", j, "_w", 0:nc)])
    fit <- fit + b[paste0("bb", j)] * 1 / (1 + exp(-z))
  }
  return(fit - mean(fit))
}

boost.fit.nnet <- function(x, y, nu, hatmatrix = FALSE, weights = NULL, ...) {
  if(!is.null(weights))
    stop("weights is not supported!")

  n <- length(y)
  i <- sample(1:n, size = ceiling(x$xt$frac * n), replace = TRUE)

  b <- nnet.fit(x$X[i, , drop = FALSE], y[i], x$xt$k, rint = x$xt$rint, sint = x$xt$sint,
    w = if(x$xt$passw) x$state$parameters[grep("bw", names(x$state$parameters))] else NULL)
  b$coefficients[1:x$xt$k] <- nu * b$coefficients[1:x$xt$k]

  x$state$parameters <- set.par(x$state$parameters, b$coefficients, "b")
  x$state$fitted.values <- predict(b, newX = x$X)
  x$state$fitted.values <- x$state$fitted.values - mean(x$state$fitted.values)
  x$state$rss <- sum((y - x$state$fitted.values)^2)

  if(hatmatrix) {
    stop("hatmatrix is not supported yet!")
  }
  
  return(x$state)
}

smooth.construct.nnet.smooth.spec <- function(object, data, knots, ...)
{
  if(is.null(object$formula)) {
    object$formula <- as.formula(paste("~", paste(object$term, collapse = "+")))
    object$dim <- length(object$term)
    object$by <- "NA"
    object$type <- "single"
    object$xt$fx <- FALSE
  }
  tp <- if(is.null(object$xt$tp)) TRUE else object$xt$tp
  object[["standardize01"]] <- object$xt[["standardize01"]] <- if(tp) FALSE else TRUE
  if(is.null(object$xt[["standardize"]]))
    object$xt[["standardize"]] <- TRUE
  if(!is.null(object$xt[["standardize"]])) {
    if(!object$xt[["standardize"]]) {
      object$xt[["standardize01"]] <- FALSE
    }
  }
  object <- smooth.construct.la.smooth.spec(object, data, knots)
  object[!(names(object) %in% c("formula", "term", "label", "dim", "X", "xt", "lasso"))] <- NULL

  object$state <- list()
  object$state$parameters <- rep(0, object$xt$k)
  names(object$state$parameters) <- paste0("bb", 1:object$xt$k)
  w <- unlist(n.weights(object$xt$k, k = ncol(object$X), type = "sigmoid",
    rint = object$xt$rint, sint = object$xt$sint, x = object$X))
  object$state$parameters <- c(object$state$parameters, w)
  object$state$fitted.values <- rep(0, nrow(object$X))

  object$special.npar <- length(grep("b", names(object$state$parameters)))
  object$binning <- list("match.index" = 1:nrow(object$X))
  if(is.null(object$xt$passw))
    object$xt$passw <- TRUE
  if(is.null(object$xt$frac))
    object$xt$frac <- 0.8

  object$update <- function(...) { stop("no nnet updating function for bfit() implemented yet!") }
  object$propose <- function(...) { stop("no nnet proposal function for GMCMC() implemented yet!") }
  object$fit.fun <- nnet.fit.fun
  object$boost.fit <- boost.fit.nnet

  class(object) <- c("nnet.smooth", "no.mgcv", "special")

  object
}


Predict.matrix.nnet.smooth <- function(object, data)
{
  object[["standardize"]] <- standardize <- if(is.null(object$xt[["standardize"]])) TRUE else object$xt[["standardize"]]
  object[["standardize01"]] <- if(standardize) TRUE else FALSE
  return(Predict.matrix.lasso.smooth(object, data))
}

#smooth.construct.pnn.smooth.spec <- function(object, data, knots, ...)
#{
#  if(is.null(object$formula)) {
#    object$formula <- as.formula(paste("~", paste(object$term, collapse = "+")))
#    object$dim <- length(object$term)
#    object$by <- "NA"
#    object$type <- "single"
#    object$xt$fx <- FALSE
#  }
#  object[["standardize01"]] <- object$xt[["standardize01"]] <- TRUE
#  if(is.null(object$xt[["standardize"]]))
#    object$xt[["standardize"]] <- TRUE
#  if(!is.null(object$xt[["standardize"]])) {
#    if(!object$xt[["standardize"]]) {
#      object$xt[["standardize01"]] <- FALSE
#    }
#  }
#  object$xt$m1p1 <- TRUE
#  object <- smooth.construct.la.smooth.spec(object, data, knots)
#  object[!(names(object) %in% c("formula", "term", "label", "dim", "X", "xt", "lasso"))] <- NULL
#  if(is.null(object$xt$degree))
#    object$xt$degree <- ncol(object$X) * 20
#  object$X <- as.matrix(polyreg::getPoly(object$X, object$xt$degree)$xdata)
#plot2d(object$X ~ dtrain$x)
#  object$S <- list(diag(1, ncol(object$X)))
#  object$update <- function(...) { stop("no nnet updating function for bfit() implemented yet!") }
#  object$propose <- function(...) { stop("no nnet proposal function for GMCMC() implemented yet!") }

#  object$N <- apply(object$X, 2, function(x) {
#    return((1/crossprod(x)) %*% t(x))
#  })

#  object$boost.fit <- function(x, y, nu, hatmatrix = FALSE, weights = NULL, nthreads = 1, ...) {
#    ## process weights.
#    if(!is.null(weights))
#      stop("weights are not supported for n()!")

#    ncX <- ncol(x$X)
#    g2 <- rep(0, ncX)

#    bf <- boost_fit_nnet(nu, x$X, x$N, y, x$binning$match.index, nthreads = nthreads)
#    j <- which.min(bf$rss)
#    g2[j] <- bf$g[j]
#    x$state$fitted.values <- bf$fit[, j]
#    x$state$rss <- bf$rss[j]

#    names(g2) <- paste0("b", 1:ncX)
#  
#    ## Finalize.
#    x$state$parameters <- set.par(x$state$parameters, g2, "b")

#    if(hatmatrix) {
#      stop("not supported for n()!")
#    }
#  
#    return(x$state)
#  }

#  object$by <- "NA"

#  class(object) <- c("pnn.smooth", "mgcv.smooth")

#  object
#}

#Predict.matrix.pnn.smooth <- function(object, data)
#{
#  object[["standardize"]] <- standardize <- if(is.null(object$xt[["standardize"]])) TRUE else object$xt[["standardize"]]
#  object[["standardize01"]] <- TRUE
#  object$xt$m1p1 <- TRUE
#  X <- Predict.matrix.lasso.smooth(object, data)
#  X <- as.matrix(polyreg::getPoly(X, object$xt$degree)$xdata)
#  return(X)
#}


predictn <- function(object, newdata, model = NULL, mstop = NULL, type = c("link", "parameter"))
{
  type <- match.arg(type)
  family <- object$family
  tl <- term.labels2(object, type = 2, intercept = FALSE)
  if(!is.null(model))
    tl <- tl[model]
  p <- vector(mode = "list", length = length(tl))
  names(p) <- names(tl)
  for(i in names(tl)) {
    fit <- 0
    if(any(j <- (grepl("n(", tl[[i]], fixed = TRUE)) & !grepl("lin(", tl[[i]], fixed = TRUE))) {
      for(tj in tl[[i]][j]) {
        if("nnet.smooth" %in% class(object$x[[i]]$smooth.construct[[tj]])) {
          pt <- predict(object, newdata = newdata, model = i,
            term = tj, FUN = function(x) { x }, intercept = FALSE)
          ## pt <- .predict_nn1(object$x[[i]]$smooth.construct[[tj]], object$parameters)
          fit <- fit + t(apply(pt, 1, cumsum))
        } else {
          fit <- fit + predict(object, newdata = newdata, model = i, term = tj, intercept = FALSE)
        }
      }
    }
    tj <- if(length(tl[[i]][!j])) tl[[i]][!j] else NULL
    fit <- fit + predict(object, newdata = newdata, model = i, term = tj, intercept = TRUE)
    p[[i]] <- if(is.null(mstop)) fit else fit[, if(is.list(mstop)) mstop$mstop else mstop]
    if(type != "link") {
      link <- family$links[i]
      if(length(link) > 0) {
        if(link != "identity") {
          linkinv <- make.link2(link)$linkinv
          p[[i]] <- linkinv(p[[i]])
        }
      } else {
        warning(paste("could not compute predictions on the scale of parameter",
          ", predictions on the scale of the linear predictor are returned!", sep = ""))
      }
    }
  }
  if(length(p) < 2)
    p <- p[[1]]
  return(p)
}


#smooth.construct.nn.smooth.spec <- function(object, data, knots, ...)
#{
#  form <- as.formula(paste("~", paste(object$term, collapse = "+")))
#  term <- object$term
#  if(is.null(object$xt$tp))
#    object$xt$tp <- TRUE
#  object <- n(form, k = object$bs.dim, tp = object$xt$tp)
#  object$label <- paste0("s(",  paste(term, collapse = ","), ")")
#  object$formula <- form
#  object <- smooth.construct.nnet.smooth.spec(object, data, knots, ...)
#  object$plot.me <- TRUE
#  object$dim <- length(term)
#  object$fixed <- FALSE
#  object$term <- term
#  object
#}


## Random bits.
rb <- function(..., k = 50)
{
  ret <- la(..., k = k)
  ret$label <- gsub("la(", "rb(", ret$label, fixed = TRUE)
  if(!is.null(ret$xt$id)) {
    lab <- strsplit(ret$label, "")[[1]]
    lab <- paste(c(lab[-length(lab)], paste(",id='", ret$xt$id, "')", sep = "")), collapse = "", sep = "")
    ret$label <- lab
  }
  if(is.null(ret$xt$tp))
    ret$xt$tp <- FALSE
  class(ret) <- "randombits.smooth.spec"
  ret
}

BitsMat <- function(X, w, thres = TRUE) {
  B <- matrix(0, nrow = nrow(X), ncol = length(w))
  for(i in 1:length(w)) {
    z <- X[, attr(w[[i]], "id"), drop = FALSE] %*% w[[i]]
    if(thres)
      attr(w[[i]], "thres") <- z[attr(w[[i]], "thres")]
    ##B[, i] <- c(1, -1)[(z >= attr(w[[i]], "thres")) + 1L]
    B[, i] <- 1 * (z >= attr(w[[i]], "thres"))
  }
  if(thres) {
    return(list("X" = B, "weights" = w))
  } else {
    return(B)
  }
}

smooth.construct.randombits.smooth.spec <- function(object, data, knots, ...)
{
  object$X <- model.matrix(object$formula, data = as.data.frame(data))[, -1, drop = FALSE]
  colnames(object$X) <- paste0("b.", colnames(object$X))
  center <- scale <- rep(NA, ncol(object$X))
  for(j in 1:ncol(object$X)) {
    if(length(unique(object$X[, j])) > 3) {
      center[j] <- mean(object$X[, j])
      scale[j] <- sd(object$X[, j])
      object$X[, j] <- (object$X[, j] - center[j]) / scale[j]
    }
  }
  object$scale <- list("center" = center, "scale" = scale)
  k <- object$xt$k
  if(!is.null(object$xt$weights))
    k <- length(object$xt$weights)
  object$xt$ntake <- if(is.null(object$xt$ntake)) {
    if(ncol(object$X) < 5)
      ncol(object$X)
    else
      ceiling(ncol(object$X) / 3)
  } else object$xt$ntake
  thres <- is.null(object$xt$weights)
  if(thres) {
    object$sample_weights <- function(...) {
      weights <- vector(mode = "list", length = k)
      smp <- if(object$xt$ntake < 5) FALSE else TRUE
      for(i in 1:k) {
        weights[[i]] <- rnorm(object$xt$ntake)
        if(!smp) {
          weights[[i]] <- rnorm(object$xt$ntake)
          attr(weights[[i]], "id") <- sample(1:ncol(object$X), size = object$xt$ntake, replace = FALSE)
        } else {
          nid <- sample(3:object$xt$ntake, size = 1)
          weights[[i]] <- rnorm(nid)
          attr(weights[[i]], "id") <- sample(1:ncol(object$X), size = nid, replace = FALSE)
        }
        attr(weights[[i]], "thres") <- sample(1:nrow(object$X), size = 1L)
      }
      return(weights)
    }
    object$xt$weights <- object$sample_weights()
  }
  tXw <- BitsMat(object$X, object$xt$weights, thres = thres)
  if(thres) {
    object$X <- tXw$X
    object$xt$weights <- tXw$weights
  } else {
    object$X <- tXw
  }
  object$xt$cmeans <- colMeans(object$X)
  object$X <- object$X - rep(object$xt$cmeans, rep.int(nrow(object$X), ncol(object$X)))
  d <- apply(apply(object$X, 2, range), 2, diff)
  object$Xkeep <- which((d > 0.00001) & !duplicated(t(object$X)))
  object$X <- object$X[, object$Xkeep, drop = FALSE]

  if(ncol(object$X) < 1)
    stop("please check your rb() specifications, no columns in the design matrix!")

  if(is.null(object$xt$weights))
    object$xt$weights <- tXw$weights
  rm(tXw)

  df <- ncol(object$X)
  const <- object$xt$const
  if(is.null(const))
    const <- 1e-05

  object$S <- list()
  pt <- object$xt$pt
  if(is.null(pt))
    pt <- "ridge"
  pt <- tolower(pt)
  if("ridge" %in% pt) {
    object$S[[1]] <- diag(ncol(object$X))
    attr(object$S[[1]], "npar") <- ncol(object$X)
  }
  if("lasso" %in% pt) {
    np <- length(object$S) + 1L
    object$S[[np]] <- function(parameters, ...) {
      b <- get.par(parameters, "b")
      A <- df / sqrt(b^2 + const)
      A <- if(length(A) < 2) matrix(A, 1, 1) else diag(A)
      A
    }
    attr(object$S[[np]], "npar") <- ncol(object$X)
  }

  object$rank <- ncol(object$X)

  object$xt$prior <- "ig"
  object$fx <- object$xt$fx <- FALSE
  object$xt$df <- 4
  object$by <- "NA"
  object$null.space.dim <- 0
  object$bs.dim <- ncol(object$X)
  # object$rank <- qr(object$S[[1]](runif(df)))$rank
  object$N <- apply(object$X, 2, function(x) {
    return((1/crossprod(x)) %*% t(x))
  })
  if(is.null(object$xt$K))
    object$xt$K <- 1
  object$boost.fit <- function(x, y, nu, hatmatrix = FALSE, weights = NULL, nthreads = 1, ...) {
    ## process weights.
    if(!is.null(weights))
      stop("weights is not supported!")

    bf <- boost_fit_nnet(nu/x$xt$K, x$X, x$N, y, x$binning$match.index, nthreads = nthreads)

    j <- which.min(bf$rss)
    g2 <- rep(0, length(bf$g))
    g2[j] <- bf$g[j]
  
    ## Finalize.
    x$state$parameters <- set.par(x$state$parameters, g2, "b")
    x$state$fitted.values <- bf$fit[, j]

    x$state$rss <- bf$rss[j]

    if(hatmatrix) {
      x$state$hat <- nu/x$xt$K * x$X[, j] %*% (1/crossprod(x$X[, j])) %*% t(x$X[, j])
    }
  
    return(x$state)
  }
  class(object) <- "randombits.smooth"
  object
}

make_weights <- function(object, data, dropout = 0.2) {
  object$xt$dropout <- dropout
  rval <- smooth.construct(object, data, NULL)
  if(length(i <- grep("weights", names(rval))))
    return(rval[[i]])
  else
    return(rval$xt[[grep("weights", names(rval$xt))]])
}

Predict.matrix.randombits.smooth <- function(object, data)
{
  X <- model.matrix(object$formula, data = as.data.frame(data))[, -1, drop = FALSE]
  colnames(X) <- paste0("b.", colnames(X))
  for(j in 1:ncol(X)) {
    if(!is.na(object$scale$center[j]))
      X[, j] <- (X[, j] - object$scale$center[j]) / object$scale$scale[j]
  }
  X <- BitsMat(X, object$xt$weights, thres = FALSE)
  X <- X - rep(object$xt$cmeans, rep.int(nrow(X), ncol(X)))
  X <- X[, object$Xkeep, drop = FALSE]
  X
}


if(FALSE) {
  set.seed(123)

  nobs <- 2000
  d <- data.frame(
    "x1" = runif(nobs, 0, pi),
    "lon" = runif(nobs, -3, 3),
    "lat" = runif(nobs, -3, 3),
    "z1" = runif(nobs, 0, 1),
    "z2" = runif(nobs, 0, 1),
    "z3" = runif(nobs, 0, 1),
    "z4" = runif(nobs, 0, 1),
    "z5" = runif(nobs, 0, 1)
  )

  d$eta0 <- with(d, 10 * sin(x1) + 10 * cos(lon)*cos(lat))
  d$eta1 <- with(d, 10 * sin(pi * z1 * z2) + 20 * (z3 - 0.5)^2 + 10 * z4 + 5 * z5)
  d$eta1 <- d$eta1 - mean(d$eta1)

  sigma <- sqrt(var(d$eta0 + d$eta1) / 5)

  d$y <- d$eta0 + d$eta1 + rnorm(nobs, sd = sigma)

  f <- y ~ s(x1,bs="cc") + s(lon,lat,k=30) + n(~z1+z2+z3+z4+z5,k=20,split=TRUE)

  b <- bamlss(f, data = d, sampler = FALSE, optimizer = boost, nu = 0.1, maxit = 2000)

  p <- boost.nnet.predict(b, model = "mu")

  plot(p ~ d$eta1)
  abline(a = 0, b = 1)


  d <- GAMart()
  b <- bamlss(num ~ n(x1) + n(x2) + n(x3), data = d, sampler = FALSE, optimizer = boost, nu = 0.01)

  p1 <- boost.nnet.predict(b, model = "mu", term = "n(x1)")
  p2 <- boost.nnet.predict(b, model = "mu", term = "n(x2)")
  p3 <- boost.nnet.predict(b, model = "mu", term = "n(x3)")

  par(mfrow = c(1, 3))
  plot2d(p1 ~ d$x1)
  plot2d(p2 ~ d$x2)
  plot2d(p3 ~ d$x3)
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
  if(inherits(knots, c("spatial.design", "spatialDesign")))
    knots <- knots[, 1:2]
  if(is.null(dim(knots)))
    knots <- matrix(as.numeric(knots), ncol = 2)
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
      phi = object$xt$phi, v = object$xt[["v"]], c = object$xt[["c"]])
  } else { # caution: xt$c partially matches xt$center
    knots <- if(is.null(object$xt$knots)) object$bs.dim else object$xt$knots
    D <- krDesign2D(data[[object$term[1]]], data[[object$term[2]]],
      knots = knots,
      phi = object$xt$phi, v = object$xt[["v"]], c = object$xt[["c"]],
      psi = object$xt$psi, delta = object$xt$delta,
      isotropic = object$xt$isotropic)
  }

  X <- D$B
  object$X <- X
  object$S <- list(D$K)
  object$rank <- qr(D$K)$rank
  object$knots <- D$knots
  object$null.space.dim <- ncol(D$K)

  object$xt$phi <- D$phi
  object$xt[["v"]] <- D[["v"]]

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
    as.numeric(quantile(x[1:n], c(.025, .5, .975), na.rm = TRUE))
  }, vectorize.args = "n")
  n.rep <- length(theta)
  plot(1:n.rep, theta, col = "lightgrey", xlab = "Iterations",
    ylab = ylab, type = "l", ...)
  iter <- round(seq(1, n.rep, length = n.plot + 1)[-1])
  tq <- cuq(iter, theta)
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

  if(!is.null(list(...)$pages))
    ask <- !(list(...)$pages == 1)

  if(prod(par("mfcol")) > 1L) {
    spar <- FALSE
    ask <- FALSE
  }

  ## What should be plotted?
  which.match <- c("effects", "samples", "hist-resid", "qq-resid", "wp",
    "scatter-resid", "max-acf", "param-samples", "boost_summary", "results",
    "max-acf")
  if(!is.character(which)) {
    if(any(which > 9L))
      which <- which[which <= 9L]
    which <- which.match[which]
  } else which <- which.match[pmatch(tolower(which), which.match)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")
  if(length(which) < 2) {
    if(which == "results")
      which <- "effects"
  }
  which <- which[which != "results"]

  ok <- any(c("hist-resid", "qq-resid", "wp") %in% which)

  x$formula <- as.formula(x$formula)

  if(length(which) > 1 | ok) {
    which2 <- which[which %in% c("hist-resid", "qq-resid", "wp")]
    if(length(which2)) {
      c95 <- list(...)$c95
      if(is.null(c95))
        c95 <- FALSE
      FUN <- list(...)$FUN
      if(c95)
        FUN <- identity
      if(is.null(FUN))
        FUN <- function(x) { mean(x, na.rm = TRUE) } 
      res <- residuals.bamlss(x, FUN = FUN, ...)
      plot(res, which = which2, spar = spar, ...)
    } else {
      for(w in which) {
        plot.bamlss(x, model = model, term = term, which = w,
          parameters = parameters, ask = ask, spar = spar, ...)
      }
    }
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
        if(spar)
          par(mfrow = if(np <= 4) c(np, 2) else c(4, 2))
        devAskNewPage(ask)
        tx <- as.vector(time(samps))
        for(j in 1:np) {
          traceplot2(samps[, j, drop = FALSE], main = "")
          mtext(paste("Trace of", snames[j]), side = 3, line = 1, font = 2)
          lines(lowess(tx, samps[, j]), col = 2)
          ##lines(fitted(gam(samps[, j] ~ s(tx,bs="ps"), method = "REML")), col = 2)
          nu <- length(unique(samps[, j, drop = FALSE]))
          acf(if(nu < 2) jitter(samps[, j, drop = FALSE]) else samps[, j, drop = FALSE], main = "", ..., na.action = na.pass)
          mtext(paste("ACF of", snames[j]), side = 3, line = 1, font = 2)
        }
      } else {
        snames <- snames[!grepl(".edf", snames, fixed = TRUE)]
        snames <- snames[!grepl(".alpha", snames, fixed = TRUE)]
        snames <- snames[!grepl("logLik", snames, fixed = TRUE)]
        samps <- samps[, snames, drop = FALSE]
        macf <- apply(samps, 2, function(x) { acf(x, plot = FALSE, ...) })
        acfx <- macf[[1]]
        acfx$acf <- array(apply(do.call("rbind", lapply(macf, function(x) { x$acf })), 2, max, na.rm = TRUE), dim = c(length(acfx$acf), 1L, 1L))
        args <- list(...)
        if(is.null(args$main))
          args$main <- ""
        if(is.null(args$xlab))
          args$xlab <- "Lag"
        if(is.null(args$ylab))
          args$ylab <- "ACF"
        getS3method("plot", class = "acf")(acfx, main = args$main, xlab = args$xlab,
          ylab = args$ylab, xlim = args$xlim, ylim = args$ylim)
        mtext("Maximum ACF of samples", side = 3, line = 1, font = 2)
      }
    }

    if(which == "effects") {
      if(is.null(x$results)) {
        xres <- results.bamlss.default(x)
        any_s <- any(sapply(names(xres), function(i) { !is.null(xres[[i]]$s.effects) } ))
        if(!any_s) {
          plot.bamlss(x, which = c("hist-resid", "qq-resid"), ...)
        } else {
          plot(xres, model = model, term = term, ask = ask, spar = spar, ...)
        }
      } else {
        any_s <- any(sapply(names(x$results), function(i) { !is.null(x$results[[i]]$s.effects) } ))
        if(!any_s) {
          plot.bamlss(x, which = c("hist-resid", "qq-resid"), ...)
        } else {
          plot(x$results, model = model, term = term, spar = spar, ask = ask, ...)
        }
      }
    }

    if(which == "boost_summary") {
      if(!is.null(x$boost_summary))
        plot(x$boost_summary, ...)
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
          rdens <- density(res, na.rm = TRUE)
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
            args2$main <- "Normal Q-Q plot"
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
    if(all(c("2.5%", "97.5%") %in% names(x)) & (attr(x, "specs")$dim < 2)) {
      args[[lim]] <- range(x[, c("2.5%", "97.5%")], na.rm = TRUE)
    } else {
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
            not = c("border", "lwd", "lty", names(formals("colorlegend")), "main", "names", "names_id", "shift")))
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
  mw <- any(grepl("mean", tolower(colnames(x)), fixed = TRUE))
  if(mw) {
    fx <- unlist(x[, grepl("mean", tolower(colnames(x)), fixed = TRUE)])
  } else {
    fx <- unlist(x[, grepl("50", colnames(x), fixed = TRUE)])
  }
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
  object$formula <- as.formula(object$formula)
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
     sterms = FALSE, full.names = FALSE, list = TRUE, parameters = parameters,
     summary = TRUE, mm = TRUE, ...)
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
      alpha <- NULL
      if("alpha" %in% rownames(x$model.matrix[[i]])) {
        j <- which(rownames(x$model.matrix[[i]]) == "alpha")
        alpha <- x$model.matrix[[i]][j, , drop = FALSE]
        x$model.matrix[[i]] <- x$model.matrix[[i]][-j, , drop = FALSE]
        if("parameters" %in% colnames(alpha)) {
          j <- which(colnames(alpha) == "parameters")
          alpha <- alpha[, -j, drop = FALSE]
        }
      }
      printCoefmat(x$model.matrix[[i]], digits = digits)
      if(!is.null(alpha)) {
        cat("-\nAcceptance probability:\n")
        printCoefmat(alpha, digits = digits)
      }
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
          cat(if(k > 1) " " else "", j, " = ", if(is.numeric(x$model.stats$sampler[[j]])) {
              round(x$model.stats$sampler[[j]], digits)
            } else x$model.stats$sampler[[j]], sep = "")
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
      nmo <- names(x$model.stats$optimizer)
      cl <- sapply(x$model.stats$optimizer, class)
      nmo <- nmo[cl != "matrix"]
      if(length(nmo)) {
        for(j in sort(nmo)) {
          if(length(x$model.stats$optimizer[[j]]) < 2) {
            if(is.numeric(x$model.stats$optimizer[[j]])) {
              ok <- TRUE
              cat(if(k > 1) " " else "", j, " = ", round(x$model.stats$optimizer[[j]], digits), sep = "")
              k <- k + 1
              if(k == 4) {
                k <- 1
                cat("\n")
                ok <- FALSE
              }
            }
          } else {
            if(!is.list(x$model.stats$optimizer[[j]]) & (j != "parpaths")) {
              print(x$model.stats$optimizer[[j]], ...)
              ok <- FALSE
            }
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
DIC.bamlss <- function(object, ..., samples = TRUE, nsamps = NULL, newdata = NULL)
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
      family <- family(object[[i]])
      if(is.null(family$d))
        stop("no d() function available in model family object!", drop = FALSE)
      samps <- process.chains(object[[i]]$samples)
      if(is.null(samps))
        stop("samples are missing in object!")
      y <- if(is.null(newdata)) {
        model.response2(object[[i]])
      } else {
        newdata[, response.name(object[[i]])]
      }
      msamps <- matrix(apply(samps[[1]], 2, mean, na.rm = TRUE), nrow = 1)
      colnames(msamps) <- colnames(samps[[1]])
      msamps <- as.mcmc.list(list(mcmc(msamps, start = 1, end = 1, thin = 1)))
      object[[i]]$samples <- msamps
      mpar <- predict.bamlss(object[[i]], newdata = newdata,
        type = "parameter", FUN = mean, nsamps = nsamps, drop = FALSE)
      mdev <- -2 * sum(family$d(y, mpar, log = TRUE), na.rm = TRUE)
      object[[i]]$samples <- samps
      rm("samps")
      par <- predict.bamlss(object[[i]], newdata = newdata,
        type = "parameter", FUN = function(x) { x }, nsamps = nsamps, drop = FALSE)
      iter <- if(is.list(par)) ncol(par[[1]]) else ncol(par)
      dev <- rep(NA, iter)
      tpar <- mpar
      for(i in 1:iter) {
        for(j in seq_along(tpar))
          tpar[[j]] <- par[[j]][, i]
        dev[i] <- -2 * sum(family$d(y, tpar, log = TRUE), na.rm = TRUE)
      }
      pd <- mean(dev, na.rm = TRUE) - mdev
      DIC <- mdev + 2 * pd
      rval <- rbind(rval, data.frame(
        "DIC" = DIC,
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
  nd <- list(...)$newdata
  if(!is.null(nd))
    samples <- FALSE
  if(samples)
    ll <- list()
  for(j in seq_along(object)) {
    if(!is.null(nd)) {
      par <- predict(object[[j]], newdata = nd, type = "parameter")
      rn <- response_name(object[[j]])
      y <- eval(parse(text = rn), envir = nd)
      edfo <- summary(object[[j]])$model.stats$optimizer$edf
      edfs <- summary(object[[j]])$model.stats$sampler$pd
      ll <- c(ll, "logLik" = sum(family(object[[j]])$d(y, par, log = TRUE), na.rm = TRUE))
      edf <- c(edf, if(is.null(edfs)) edfo else edfs)
      next
    }
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
          ll[[j]] <- mcmc(cbind("logLik" = ll[[j]], "df" = edf),
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
      if(("boost_summary" %in% names(ms)) & is.null(mstop)) {
        llb <- ms$boost_summary$ic
        edfb <- ms$boost_summary$criterion$edf
        ll <- c(ll, llb)
        edf <- c(edf, edfb)
      }
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
        if(!is.null(ms$edf)) {
          edf <- c(edf, ms$edf)
        } else {
          if(!is.null(ms$pd)) {
            edf <- c(edf, ms$pd)
          } else {
            edf <- c(edf, NA)
          }
        }
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
      rval <- cbind("logLik" = ll, "df" = edf, "nobs" = nobs)
      if(length(mn) == nrow(rval))
        row.names(rval) <- if(nrow(rval) > 1) mn[1:nrow(rval)] else ""
    }
  } else rval <- NULL
  if(is.null(dim(rval))) {
    rn <- names(rval)
    rval <- matrix(rval, nrow = 1)
    colnames(rval) <- rn
  }
  if(nrow(rval) < 2) {
    ll <- rval[, "logLik"]
    for(j in colnames(rval)[-1]) {
      jn <- gsub("edf", "df", j)
      attr(ll, jn) <- rval[, j]
    }
    rval <- ll
    class(rval) <- "logLik"
  }
  rval
}


## Extract model formulas.
formula.bamlss.frame <- formula.bamlss <- function(x, model = NULL, ...)
{
  f <- model.terms(as.formula(x$formula), model)
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
        response <- response.name(x[[i]][[j]], keep.functions = TRUE)
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
      response <- response.name(x[[i]], keep.functions = TRUE)
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
          if(is.character(x[[i]][[j]]$formula)) {
            cat(x[[i]][[j]]$formula, "\n")
          } else print(x[[i]][[j]]$formula, showEnv = FALSE)
        }
      } else {
        attr(x[[i]], "name") <- NULL
        attr(x[[i]]$formula, "name") <- NULL
        attr(x[[i]]$formula, ".Environment") <- NULL
        if("formula" %in% names(x[[i]])) {
          if(is.character(x[[i]]$formula))
            cat(x[[i]]$formula, "\n")
          else
            print(x[[i]]$formula, showEnv = FALSE)
        } else print(x[[i]])
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
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "tx4", "la", "n", "h", "lf", "af", "lf.vd", "re", "peer", "fpc", "lin", "rb", "tree"))
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
  for(j in specials) {
    i <- grep2(paste(j, "(", sep = ""), tl, fixed = TRUE)
    if(length(i)) {
      for(ii in i) {
        s1 <- strsplit(tl[ii], "")[[1]]
        s2 <- strsplit(paste(j, "(", sep = ""), "")[[1]]
        s1 <- paste(s1[1:length(s2)], collapse = "")
        s2 <- paste(s2, collapse = "")
        if(s1 == s2)
          sid <- c(sid, ii)
      }
    }
  }
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
    for(j in specials) {
      i <- grep2(paste(j, "(", sep = ""), tl, fixed = TRUE)
      if(length(i)) {
        for(ii in i) {
          s1 <- strsplit(tl[ii], "")[[1]]
          s2 <- strsplit(paste(j, "(", sep = ""), "")[[1]]
          s1 <- paste(s1[1:length(s2)], collapse = "")
          s2 <- paste(s2, collapse = "")
          if(s1 == s2)
            sid <- c(sid, ii)
        }
      }
    }
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
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "tx4", "la", "n", "h", "lf", "af", "lf.vd", "re", "peer", "fpc", "lin", "rb", "tree"))
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
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "tx4", "la", "n", "h", "lf", "af", "lf.vd", "re", "peer", "fpc", "lin", "rb", "tree"))
  if(inherits(x, "formula"))
    x <- terms(x, specials = specials)
  if(!inherits(x, "terms"))
    stop("x must be a 'terms' object!")
  return(length(unlist(attr(x, "specials"))) > 0)
}

has_pterms <- function(x, specials = NULL)
{
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "tx4", "la", "n", "h", "lf", "af", "lf.vd", "re", "peer", "fpc", "lin", "rb", "tree"))
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
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "tx4", "la", "n", "h", "lf", "af", "lf.vd", "re", "peer", "fpc", "lin", "rb", "tree"))
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
  specials <- unique(c(specials, "s", "te", "t2", "sx", "s2", "rs", "ti", "tx", "tx2", "tx3", "tx4", "la", "n", "h", "lf", "af", "lf.vd", "re", "peer", "fpc", "lin", "rb", "tree"))
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

  x$formula <- as.formula(x$formula)

  FUN <- list(...)$FUN

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
            if(is.null(obj$smooth.construct[[j]]$mono))
              obj$smooth.construct[[j]]$mono <- 0
            if(!is.null(obj$smooth.construct[[j]]$margin)) {
              for(mj in seq_along(obj$smooth.construct[[j]]$margin)) {
                if(is.null(obj$smooth.construct[[j]]$margin[[mj]]$mono))
                  obj$smooth.construct[[j]]$margin[[mj]]$mono <- 0
              }
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
            if(!is.null(obj$smooth.construct[[j]][["X.dim"]])) {
            b <- paste(id, "s", j, paste("b", 1:obj$smooth.construct[[j]][["X.dim"]], sep = ""), sep = ".")
            } else {
              state <- obj$smooth.construct[[j]][["state"]]
              b <- names(state$parameters)
              b <- b[!grepl("tau2", b)]
            }
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
            FUN = FUN, snames = snames, data = if(!is.null(obj$smooth.construct[[j]]$model.frame)) {
              obj$smooth.construct[[j]]$model.frame
            } else mf[, tn, drop = FALSE], grid = grid)
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
  burnin = NULL, thin = NULL, coef.only = FALSE, ...)
{
  if(!inherits(object, "bamlss") & !inherits(object, "bamlss.frame"))
    stop("object is not a 'bamlss' object!")
  if(is.null(object$samples))
    stop("no samples to extract!")
  tx <- terms(object, drop = FALSE)
  x <- object$samples
  x <- process.chains(x, combine, drop = FALSE, burnin, thin)

  if(coef.only) {
    cdrop <- c(".accepted", ".alpha", "logLik", "logPost", "AIC", 
      "BIC", "DIC", "pd", ".edf")
    for(d in cdrop) {
      for(k in seq_along(x)) {
        if(length(j <- grep(d, colnames(x[[k]]), fixed = TRUE))) {
          x[[k]] <- x[[k]][, -j, drop = FALSE]
        }
      }
    }
  }

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
      x[[k]] <- x[[k]][, j, drop = FALSE]
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
    drop <- c(".tau2", ".lambda", ".edf", ".accepted", if(!summary) ".alpha" else NULL,
      "logLik", "logPost", "AIC", "BIC", "DIC", "pd")
    if(is.null(FUN))
      FUN <- function(x) { mean(x, na.rm = TRUE) }
  }
  if(!pterms)
    drop <- c(drop, ".p.")
  if(!sterms)
    drop <- c(drop, ".s.") 
  par <- samps <- NULL
  rval <- list()
  mm <- if(is.null(list(...)$mm)) FALSE else TRUE
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
      rval$samples <- apply(rval$samples, 2, function(x, ...) { FUN(na.omit(x), ...) })
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
    pedf <- rval$parameters[grep(".p.edf", names(rval$parameters), fixed = TRUE)]
    if(length(di <- grep2(drop, names(rval$parameters), fixed = TRUE)))
      rval$parameters <- rval$parameters[-di]
    if(length(pedf) & !(".p." %in% drop))
      rval$parameters <- c(rval$parameters, pedf)
    if(summary) {
      if(!mm)
        rval$parameters <- rval$parameters[grep2(c(".tau2", ".edf"), names(rval$parameters), fixed = TRUE)]
    }
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
    foo <- function(x) {
      if(is.null(dim(x)))
        return(length(x))
      else
        return(ncol(x))
    }
    rd <- sapply(rval, foo)
    if(all(rd < 1)) {
      NULL
    } else {
      if(any(rd < 1))
        rval <- rval[!(rd < 1)]
      do.call("cbind", rval)
    }
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

    y <- NULL
    if(!is.null(object$y)) {
      y <- if(is.data.frame(object$y)) {
        if(ncol(object$y) < 2) {
          object$y[[1]]
        } else object$y
      } else {
        object$y
      }
    }

    if(!is.null(nd <- list(...)$newdata)) {
      rn <- response_name(object)
      y <- nd[[rn]]
      if(is.null(y))
        stop(paste("the response", rn , "is not available in newdata!"))
      rm(nd)
      n <- if(is.null(dim(y))) length(y) else nrow(y)
    }

    if(is.null(y)) {
      rn <- response_name(object)
      y <- model.frame(object)[[rn]]
    }

    if(is.null(y))
      stop("response variable is missing, cannot compute residuals!")

    par <- predict(object, nsamps = nsamps, drop = FALSE, ...)
    nas <- attr(par, "na.action")
    if(!is.null(nas)) {
      if(is.null(dim(y))) {
        y <- y[-nas]
      } else {
        y <- y[-nas, ]
      }
    }
    nod <- is.null(dim(par[[1L]]))
    for(j in family$names) {
      if(!nod)
        par[[j]] <- as.matrix(par[[j]])
      par[[j]] <- make.link2(family$links[j])$linkinv(par[[j]])
    }

    if(type == "quantile") {
      if(is.null(family$p)) {
        type <- "response"
        warning(paste("no $p() function in family '", family$family,
          "', cannot compute quantile residuals, computing response resdiuals instead!", sep = ""))
      } else {
        discrete <- FALSE
        if(!is.null(family$type)) {
          if(tolower(family$type) == "discrete")
            discrete <- TRUE
        }
        if(family$family == "binomial")
          discrete <- TRUE
        if(discrete) {
          ymin <- min(y, na.rm = TRUE)
          a <- family$p(ifelse(y == ymin, y, y - 1), par)
          a <- ifelse(y == ymin, 0, a)
          b <- family$p(y, par)
          u <- runif(length(y), a, b)
          u <- ifelse(u > 0.999999, u - 1e-16, u)
          u <- ifelse(u < 1e-06, u + 1e-16, u)
          res <- qnorm(u)
#	  a <- family$p(y - 1, par)
#	  b <- family$p(y, par)
#	  u <- runif(n = length(y), min = a, max = b)
#	  res <- qnorm(u)
        } else {
          prob <- family$p(y, par)
          res <- qnorm(prob)
          if(any(isnf <- !is.finite(res))) {
            warning("non finite quantiles from probabilities, set to NA!")
            res[isnf] <- NA
          }
        }
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

  if(any(j <- !is.finite(res)))
    res[j] <- NA

  return(res)
}


## Residuals plotting functions.
plot.bamlss.residuals <- function(x, which = c("hist-resid", "qq-resid", "wp"), spar = TRUE, ...)
{
  ## What should be plotted?
  which.match <- c("hist-resid", "qq-resid", "wp")
  if(!is.character(which)) {
    if(any(which > 3L))
      which <- which[which <= 3L]
    which <- which.match[which]
  } else which <- which.match[pmatch(tolower(which), which.match)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")

  if(is.null(dim(x)))
    x <- matrix(x, ncol = 1)
  nc <- ncol(x)
  cn <- colnames(x)

  if(nc > 10) {
    nc <- 1
    cn <- NULL
  }

  add <- list(...)$add
  if(is.null(add))
    add <- FALSE

  if(add)
    spar <- FALSE

  if(spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = n2mfrow(length(which) * nc))
  }

  type <- attr(x, "type")

  for(j in 1:nc) {
    for(w in which) {
      args <- list(...)
      if(w == "hist-resid") {
        if(ncol(x) > 1) {
          x2 <- rowMeans(x, na.rm = TRUE)
        } else {
          x2 <- x
        }
        rdens <- density(as.numeric(x2), na.rm = TRUE)
        rh <- hist(as.numeric(x2), plot = FALSE)
        args$ylim <- c(0, max(c(rh$density, rdens$y)))
#        if(is.null(args$xlim)) {
#          args$xlim <- range(x[is.finite(x)], na.rm = TRUE)
#          args$xlim <- c(-1, 1) * max(args$xlim)
#        }
        args$freq <- FALSE
        args$x <- as.numeric(x2)
        args <- delete.args("hist.default", args, package = "graphics", not = c("xlim", "ylim"))
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
        if(ncol(x) > 1) {
          x2 <- t(apply(x, 1, c95))

          args$x <- NULL
          args$plot.it <- FALSE
          args <- delete.args("qqnorm.default", args, package = "stats", not = c("col", "pch", "cex"))
          if(is.null(args$main))
            args$main <- paste("Normal Q-Q plot", if(!is.null(cn[j])) paste(":", cn[j]) else NULL)

          args$y <- x2[, "Mean"]
          mean <- do.call(qqnorm, args)
          args$y <- x2[, "2.5%"]
          lower <- do.call(qqnorm, args)
          args$y <- x2[, "97.5%"]
          upper <- do.call(qqnorm, args)

          ylim <- range(c(as.numeric(mean$y), as.numeric(lower$y), as.numeric(upper$y)),
            na.rm = TRUE)
          args$plot.it <- TRUE
          if(is.null(args$ylim))
            args$ylim <- ylim
          args$y <- x2[, "Mean"]
          mean <- do.call(qqnorm, args)

          if(is.null(args$ci.col))
            args$ci.col <- 1
          if(is.null(args$ci.lty))
            args$ci.lty <- 2

          lines(lower$x[order(lower$x)], lower$y[order(lower$x)], lty = args$ci.lty, col = args$ci.col)
          lines(upper$x[order(upper$x)], upper$y[order(upper$x)], lty = args$ci.lty, col = args$ci.col)

          args$y <- x2[, "Mean"]
          qqline(args$y)
        } else {
          args$y <- as.numeric(x)
#          if(is.null(args$ylim)) {
#            args$ylim <- range(x[is.finite(x)], na.rm = TRUE)
#            args$ylim <- c(-2.5, 2.5) * max(args$ylim)
#          }
#          if(is.null(args$xlim)) {
#            args$xlim <- range(x[is.finite(x)], na.rm = TRUE)
#            args$xlim <- c(-2.5, 2.5) * max(args$xlim)
#          }
          args$x <- NULL
          args <- delete.args("qqnorm.default", args, package = "stats", not = c("col", "pch", "xlim", "ylim", "cex"))
          if(is.null(args$main))
            args$main <- paste("Normal Q-Q plot", if(!is.null(cn[j])) paste(":", cn[j]) else NULL)
          args$plot.it <- !add
          ok <- try(do.call(qqnorm, args))
          if(add) {
            args <- delete.args("points.default", list(...), package = "graphics",
               not = c("col", "pch", "cex"))
            points(ok$x, ok$y, pch = args$pch, col = args$col, cex = args$cex)
          } else {
            if(!inherits(ok, "try-error"))
              qqline(args$y) ## abline(0,1)
          }
        }
      }
      if(w == "wp") {
        xlo <- xup <- NULL
        if(ncol(x) > 1) {
          x2 <- t(apply(x, 1, c95))
          xlo <- x2[, "2.5%"]
          xup <- x2[, "97.5%"]
          x2 <- x2[, "Mean"]
        } else {
          x2 <- x
        }
        d <- qqnorm(x2, plot = FALSE)

        probs <- c(0.25, 0.75)
        y3 <- quantile(x2, probs, type = 7, na.rm = TRUE)
        x3 <- qnorm(probs)
        slope <- diff(y3)/diff(x3)
        int <- y3[1L] - slope * x3[1L]
        d$y <- d$y - (int + slope * d$x)
        ##d$y <- d$y - d$x
        if(!is.null(xlo)) {
          d2 <- qqnorm(xlo, plot = FALSE)
          d$ylo <- d2$y - d2$x
          d$xlo <- d2$x
          d2 <- qqnorm(xup, plot = FALSE)
          d$yup <- d2$y - d2$x
          d$xup <- d2$x
        }
        level <- 0.95
        xlim <- max(abs(d$x), na.rm = TRUE)
        xlim <- c(-xlim, xlim)
        ylim <- max(abs(c(as.numeric(d$y), as.numeric(d$ylo), as.numeric(d$yup))), na.rm = TRUE)
        ylim <- c(-ylim, ylim)
        if(!is.null(args$ylim2))
          ylim <- args$ylim2
        if(!is.null(args$xlim2))
          xlim <- args$xlim2
        z <- seq(xlim[1] - 10, xlim[2] + 10, 0.25)
        p <- pnorm(z)
        se <- (1/dnorm(z)) * (sqrt(p * (1 - p)/length(d$y)))
        low <- qnorm((1 - level)/2) * se
        high <- -low
        args <- list(...) 
        if(is.null(args$col))
          args$col <- 1
        if(is.null(args$pch))
          args$pch <- 1
        if(is.null(args$cex))
          args$cex <- 1
        if(is.null(args$ylab))
          args$ylab <- "Deviation"
        if(is.null(args$xlab))
          args$xlab <- "Unit normal quantile"
        if(add) {
          points(d$x, d$y, col = args$col, pch = args$pch, cex = args$cex)
        } else {
          if(is.null(args$main))
            args$main <- paste("Worm plot", if(!is.null(cn[j])) paste(":", cn[j]) else NULL)
          plot(d$x, d$y, ylab = args$ylab, xlab = args$xlab, main = args$main,
            xlim = xlim, ylim = ylim, col = NA, type = "n")
          grid(lty = "solid")
          abline(0, 0, lty = 2, col = "lightgray")
          abline(0, 1e+05, lty = 2, col = "lightgray")
          lines(z, low, lty = 2)
          lines(z, high, lty = 2)
          points(d$x, d$y, col = args$col, pch = args$pch, cex = args$cex)
        }
        if(!is.null(xlo)) {
          if(is.null(args$ci.col))
            args$ci.col <- 4
          if(is.null(args$ci.lty))
            args$ci.lty <- 2
          i <- order(d$xlo)
          lines(d$ylo[i] ~ d$xlo[i], lty = args$ci.lty, col = args$ci.col)
          i <- order(d$xup)
          lines(d$yup[i] ~ d$xup[i], lty = args$ci.lty, col = args$ci.col)
        }
      }
    }
  }

  return(invisible(NULL))
}


c.bamlss.residuals <- function(...) {
  res <- list(...)
  Call <- match.call()
  mn <- as.character(Call)[-1L]
  names(res) <- if(is.null(names(res))) mn else names(res)
  class(res) <- "bamlss.residuals.list"
  return(res)
}

plot.bamlss.residuals.list <- function(x, ...) {
  class(x) <- "list"
  x <- as.data.frame(x)
  args <- list(...)
  ylim <- args$ylim
  if(is.null(ylim))
    ylim <- range(x, na.rm = TRUE)
  col <- args$col
  if(is.null(col))
    col <- 1:ncol(x)
  col <- rep(col, length.out = ncol(x))
  for(j in 1:ncol(x)) {
    plot(x[[j]], ylim = ylim, spar = FALSE, add = j > 1, col = col[j], ...)
  }
  legend <- args$legend
  if(is.null(legend))
    legend <- TRUE
  if(legend) {
    pos <- args$pos
    if(is.null(pos))
      pos <- "topleft"
    cex2 <- args$cex2
    if(is.null(cex2))
      cex2 <- 1
    pch <- args$pch
    if(is.null(pch))
      pch <- 1
    bty <- args$bty
    if(is.null(bty))
      bty <- "o"
    bg <- args$bg
    if(is.null(bg))
      bg <- "white"
    legend(pos, names(x), bty = bty, col = col, pch = pch, cex = cex2, bg = bg)
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
  p <- try(chol(x), silent = TRUE)
  p <- if(inherits(p, "try-error")) {
    try(solve(x), silent = TRUE)
  } else {
    try(chol2inv(p), silent = TRUE)
  }
  if(inherits(p, "try-error")) {
    x <- x + diag(0.0001, ncol(x))
    p <- try(solve(x), silent = TRUE)
  }
  if(inherits(p, "try-error") & force) {
    p <- diag(ncol(x))
  }
  if(is.null(dim(p))) {
    p <- matrix(p, 1, 1)
    rownames(p) <- rn[1]
    colnames(p) <- cn[1]
  } else {
    rownames(p) <- rn
    colnames(p) <- cn
  }
  return(p)
}


#if(FALSE) {
#  set.seed(1234)
#  a <- list()
#  sparse <- list()
#  for (i in 1:4){
#    a[[i]] <- crossprod(matrix(rnorm(i^2), nrow = i))
#  }

#  A <- as.matrix(do.call(bdiag, a))

#  sparse$block.index <- list(
#    1L,
#    2L:3L,
#    4L:6L,
#    7L:10L
#  )
#  sparse$is.diagonal <- FALSE

#  # calculate inverse
#  inv1 <- bamlss:::matrix_inv(A, sparse)
#  attr(inv1, "dimnames") <- NULL
#  inv2 <- bamlss:::matrix_inv(A)

#  all.equal(inv1, inv2)
#  inv1
#  inv2
#}


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


#if(FALSE) {
#  x <- runif(3000, -3, 3)
#  f <- bamlss:::simfun("pick")
#  y <- sin(x) + rnorm(3000, sd = scale2(f(scale2(x, 0, 1)), 0.01, 0.3))
#  plot(x, y)
#  b <- bamlss(list(y ~ n(x), ~ n(x)), sampler = FALSE)
#}

## Most likely transformations.
h <- function(...)
{
  ret <- te(...)
  class(ret$margin[[1]]) <- "ispline.smooth.spec"
  ret$label <- gsub("te(", "h(", ret$label, fixed = TRUE)
  ret$special <- TRUE
  class(ret) <- "mlt.smooth.spec"
  ret
}

smooth.construct.ispline.smooth.spec <- function(object, data, knots, ...)
{
  stopifnot(requireNamespace("splines2"))
  class(object) <- "ps.smooth.spec"
  object <- smooth.construct(object, data, knots)
  knots <- object$knots
  knots <- knots[(knots > min(data[[object$term]])) &
    (knots < max(data[[object$term]]))]
  object$knots <- knots
  object$X <- splines2::iSpline(data[[object$term]],
    knots = knots, degree = 3, intercept = FALSE, derivs = 0L)
  object$S <- list(crossprod(diff(diag(ncol(object$X)), differences = 2)))
  object$derivMat <- splines2::iSpline(data[[object$term]],
    knots = knots, degree = 3, intercept = FALSE, derivs = 1L)
  class(object) <- "ispline.smooth"
  return(object)
}

smooth.construct.mlt.smooth.spec <- function(object, data, knots, ...)
{
  stopifnot(requireNamespace("splines2"))
  object$margin[[1]] <- smooth.construct(object$margin[[1]], data, knots)
  if(length(object$margin) > 1) {
    for(j in 2:length(object$margin))
      object$margin[[j]] <- smoothCon(object$margin[[j]], data, knots, absorb.cons = TRUE,scale.penalty=TRUE)[[1]]
  }
  object$X <- tensor.prod.model.matrix(lapply(object$margin, function(x) { x$X } ))
  dX <- list(object$margin[[1]]$derivMat)
  if(length(object$margin) > 1)
    dX <- c(dX, lapply(object$margin[2:length(object$margin)], function(x) { x$X } ))
  object$dX <- tensor.prod.model.matrix(dX)
  object$S <- tensor.prod.penalties(sapply(object$margin, function(x) { x$S } ))
  object$force.positive <- TRUE
  class(object) <- "mlt.smooth"
  object
}

Predict.matrix.ispline.smooth <- function(object, data)
{
  stopifnot(requireNamespace("splines2"))
  X <- splines2::iSpline(data[[object$term]],
    knots = object$knots, degree = 3, intercept = FALSE, derivs = 0L)
  attr(X, "deriv") <- splines2::iSpline(data[[object$term]],
    knots = object$knots, degree = 3, intercept = FALSE, derivs = 1L)
  X
}

Predict.matrix.mlt.smooth <- function(object, data)
{
  X <- lapply(object$margin, function(x) { PredictMat(x, as.data.frame(data)) })
  X <- tensor.prod.model.matrix(X)
  X
}

#bamlss_bridge <- function(object, ...)
#{
#  stopifnot(requireNamespace("bridgesampling"))

#  if(!inherits(object, "bamlss"))
#    stop("object must be a 'bamlss' object!")
#  if(is.null(object$samples))
#    stop("object does not contain any samples!")

#  samps <- samples(object, coef.only = TRUE, drop = FALSE)

#  if(is.data.frame(object$y)) {
#    if(ncol(object$y) < 2)
#      object$y <- object$y[[1]]
#  }

#  if(is.null(attr(object$x, "bamlss.engine.setup")))
#    object$x <- bamlss.engine.setup(object$x)

#  logPost <- function(par, ...) {
#    log_posterior(par, object$x, object$y, object$family, verbose = FALSE)
#  }

#  lb <- rep(-Inf, ncol(samps[[1]]))
#  ub <- rep(Inf, ncol(samps[[1]]))
#  names(lb) <- names(ub) <- colnames(samps[[1]])
#  if(any(grepl("tau2", colnames(samps[[1]]))))
#    lb[grep("tau2", names(lb))] <- 0

#  bs <- bridge_sampler(samples = samps,
#    log_posterior = logPost,
#    data = list("n" = NA),
#    lb = lb, ub = ub, ...)

#  bs
#}

smooth.construct.re2.smooth.spec <- function(object, data, knots, ...)
{
  isn <- NULL
  for(j in object$term) {
    isn <- c(isn, is.numeric(data[[j]]))
  }
  if(any(isn) & FALSE) {
    center <- scale <- NULL
    for(j in which(isn)) {
      m <- mean(data[[j]])
      sd <- sd(data[[j]])
      data[[j]] <- (data[[j]] - m) / sd
      center <- c(center, m)
      scale <- c(scale, sd)
    }
    object$isn <- isn
    object$center <- center
    object$scale <- scale
  }
  object <- smooth.construct.re.smooth.spec(object, data, knots)
  class(object) <- if (inherits(object, "tensor.smooth.spec")) 
      c("random2.effect", "tensor.smooth", "random.effect")
    else c("random2.effect", "random.effect")
  object
}

Predict.matrix.random2.effect <- function(object, data) 
{
  if(!is.null(object$isn)) {
    i <- 1
    for(j in which(object$isn)) {
      data[[j]] <- (data[[j]] - object$center[i]) / object$scale[i]
      i <- i + 1
    }
  }
  X <- model.matrix(object$form, data)
  X
}

smooth.construct.sr.smooth.spec <- function(object, data, knots, ...)
{
  if(object$dim > 1)
    stop("more than one variable not supported!")

  class(object) <- "ps.smooth.spec"

  object <- smoothCon(object, as.data.frame(data), knots, absorb.cons = TRUE,scale.penalty=TRUE)[[1]]

  ev <- eigen(object$S[[1]], symmetric = TRUE)
  null.rank <- object$df - object$rank
  p.rank <- object$rank
  if(p.rank > ncol(object$X)) 
    p.rank <- ncol(object$X)
  object$xt$trans.U <- ev$vectors
  object$xt$trans.D <- c(ev$values[1:p.rank], rep(1, null.rank))
  object$xt$trans.D <- 1/sqrt(object$xt$trans.D)
  UD <- t(t(object$xt$trans.U) * object$xt$trans.D)
  object$X <- object$X %*% UD
  object$X <- object$X[, 1:p.rank, drop = FALSE]
  object$S <- list(diag(1, ncol(object$X)))

  object$xt$p.rank <- p.rank
  object$xt$df <- 1
  object$PredictMat <- Predict.matrix.srand.smooth
  object$X.dim <- ncol(object$X)

  class(object) <- c("srand.smooth", "no.mgcv")

  return(object)
}

Predict.matrix.srand.smooth <- function(object, data) 
{
  class(object) <- "pspline.smooth"
  X <- PredictMat(object, as.data.frame(data))
  UD <- t(t(object$xt$trans.U) * object$xt$trans.D)
  X <- X %*% UD
  X <- X[, 1:object$xt$p.rank, drop = FALSE]
  X
}


## Model fitting shortcuts.
boost2 <- function(...) {
  bamlss(..., sampler = FALSE, optimizer = opt_boost)
}

lasso2 <- function(...) {
  bamlss(..., sampler = FALSE, optimizer = opt_lasso)
}

bayesx2 <- function(...) {
  bamlss(..., sampler = sam_BayesX, optimizer = FALSE)
}

bboost <- function(..., data, type = 1, cores = 1, n = 2, prob = 0.623,
  fmstop = NULL, trace = TRUE, drop = FALSE, replace = FALSE)
{
  if(is.null(fmstop)) {
    fmstop <- function(model, data) {
      y <- response.name(model)
      p <- predict(model, newdata = data, model = "mu")
      mse <- NULL
      for(i in 1:nrow(model$parameters))
        mse <- c(mse, sqrt(mean((data[[y]] - p[, i])^2)))
      list("rMSE" = mse, "mstop" = which.min(mse))
    }
  }

  nobs <- nrow(data)
  ind <- 1:nobs

  if(type > 1) {
    size <- ceiling(nobs * prob)

    foo <- function(j) {
      if(trace)
        cat("... starting bootstrap sample", j, "\n")
      i <- sample(ind, size = size, replace = TRUE)
      d0 <- data[i, , drop = FALSE]
      d1 <- data[!(ind %in% i), , drop = FALSE]
      b <- bamlss(..., data = d0, plot = FALSE, sampler = FALSE,
        optimizer = opt_boost, boost.light = FALSE, light = TRUE)
      attr(b, "mstop") <- fmstop(b, d1)
      if(drop)
        b$parameters <- b$parameters[attr(b, "mstop")$mstop, ]
      if(trace)
        cat("... finished bootstrap sample", j, "\n")
      return(b)
    }
  } else {
    foo <- function(j) {
      if(trace)
        cat("... starting bootstrap sample", j, "\n")
      if(replace) {
        i <- sample(ind, size = nobs, replace = TRUE)
      } else {
        i <- sample(ind, size = ceiling(nobs * prob), replace = FALSE)
      }
      d0 <- data[i, , drop = FALSE]
      b <- bamlss(..., data = d0, plot = FALSE, boost.light = TRUE,
        light = TRUE, sampler = FALSE, optimizer = opt_boost)
      attr(b, "mstop") <- list(
        "logLik" = b$model.stats$optimizer$boost_summary$ic,
        "mstop" = nrow(b$parameters))
      b$parameters <- b$parameters[nrow(b$parameters), ]
      if(trace)
        cat("... finished bootstrap sample", j, "\n")
      return(b)
    }

    drop <- TRUE
  }

  m <- parallel::mclapply(1:n, foo, mc.cores = cores)
  class(m) <- "bboost"
  attr(m, "drop") <- drop
  m
}

bamlss.sl <- function(object, data, ...) {
  fam <- object[[1]]$family
  p <- predict.bboost(object, newdata = data)
  eta <- beta <- list()
  for(j in fam$names) {
    eta[[j]] <- do.call("cbind", lapply(p, function(x) { x[, grep(paste0(".", j), colnames(x))] }))
    beta[[j]] <- rep(0, ncol(eta[[j]]))
    colnames(eta[[j]]) <- paste0(j, 1:ncol(eta[[j]]))
  }
  k <- length(beta[[1]])
  beta <- unlist(beta)
  yname <- response.name(object[[1]])
  objfun <- function(beta) {
    par <- list()
    pen <- 0
    for(j in fam$names) {
      beta1 <- beta[paste0(j, 1:k)]
      if(any(beta1 > 0))
        beta1 <- beta1 / sum(beta1)
      par[[j]] <- drop(eta[[j]] %*% beta1)
      pen <- pen + t(beta1) %*% diag(1e-05, length(beta1)) %*% beta1
    }
    -1 * fam$loglik(data[[yname]], fam$map2par(par)) + pen
  }
  opt <- optim(beta, objfun, method = "L-BFGS-B", lower = 0)
  beta <- opt$par
  fit <- list()
  for(j in fam$names) {
    beta1 <- beta[paste0(j, 1:k)]
    if(any(beta1 > 0))
      beta1 <- beta1 / sum(beta1)
    else
      warning("All models have zero weight!", call. = FALSE)
    beta[paste0(j, 1:k)] <- beta1
    fit[[j]] <- drop(eta[[j]] %*% beta1)
  }
  rval <- list("coefficients" = beta, "fitted.values" = fit,
    "converged" = opt$convergence == 0L, "logLik" = -1 * opt$value)
  class(rval) <- "bamlss.sl"
  rval
}

predict.bamlss.sl <- function(object, models, newdata,
  type = c("link", "parameter"), ...)
{
  type <- match.arg(type)
  fam <- models[[1]]$family
  p <- predict.bboost(models, newdata = newdata, ...)
  coefficients <- object$coefficients
  fit <- list()
  for(j in fam$names) {
    eta <- do.call("cbind", lapply(p, function(x) { x[, grep(paste0(".", j), colnames(x))] }))
    fit[[j]] <- drop(eta %*% coefficients[paste0(j, 1:ncol(eta))])
    if(type != "link")
      fit[[j]] <- fam$linkfun(fit[[j]])
  }
  return(fit)
}

bboost_plot <- function(object, col = NULL)
{
  ncrit <- names(attr(object[[1]], "mstop"))
  ncrit <- ncrit[!(ncrit %in% "mstop")]
  crit <- mstops <- NULL
  for(j in 1:length(object)) {
    crit <- cbind(crit, attr(object[[j]], "mstop")[[ncrit]])
    mstops <- c(mstops, attr(object[[j]], "mstop")$mstop)
  }
  if(is.null(col))
    col <- rainbow_hcl(length(mstops))
  colnames(crit) <- paste0(ncrit, 1:ncol(crit))
  matplot(crit, type = "l", lty = 1, xlab = "Boosting iteration", ylab = ncrit, col = col)
  abline(v = mstops, col = col, lty = 2)
  return(invisible(crit))
}

plot.bboost <- function(...) {
  bboost_plot(...)
}

predict.bboost <- function(object, newdata, ..., cores = 1, pfun = NULL)
{
  if(is.null(pfun))
    pfun <- predict.bamlss
  n <- length(object)
  drop <- attr(object, "drop")
  foo <- function(j) {
    if(is.null(object[[j]]))
      return(NULL)
    if(drop) {
      p <- pfun(object[[j]], newdata = newdata, ...)
    } else {
      mstop <- attr(object[[j]], "mstop")
      if(is.list(mstop))
        mstop <- mstop$mstop
      p <- pfun(object[[j]], newdata = newdata, mstop = mstop, ...)
    }
    if(is.list(p)) {
      p <- do.call("cbind", p)
      colnames(p) <- paste0(paste0("m", j), ".", colnames(p))
    } else
      p <- as.numeric(p)
    return(p)
  }
  pred <- parallel::mclapply(1:n, foo, mc.cores = cores)
  if(is.null(dim(pred[[1]])))
    pred <- do.call("cbind", pred)
  return(pred)
}


vfun <- function(gamma, constr) {
  v <- diff(drop(gamma))
  if(constr < 2)
    v <- (v < 0) * 1
  else
    v <- (v > 0) * 1
  v
}


pathplot <- function(object, ...)
{
  if(is.null(object$model.stats$optimizer)) {
    warning("there is nothing to plot!")
  }
  if(!is.null(object$model.stats$optimizer$boost_summary))
    boost_plot(object, ...)
  if(!is.null(object$model.stats$optimizer$lasso.stats))
    lasso_plot(object, ...)
  if(!is.null(object$model.stats$optimizer$parpaths))
    bbfit_plot(object, ...)
  return(invisible(NULL))
}

smooth.construct.ms.smooth.spec <- function(object, data, knots, ...)
{
  class(object) <- "ps.smooth.spec"
  object <- smooth.construct.ps.smooth.spec(object, data, knots)
  if(is.null(object$xt$constr))
    object$xt$constr <- 1

  object$xt$force.center <- TRUE

  object$boost.fit <- function(x, y, nu, hatmatrix = FALSE, weights = NULL, nthreads = 1, ...) {
    ## process weights.
    if(is.null(weights))
      weights <- rep(1, length = length(y))

    ## Compute reduced residuals.
    xbin.fun(x$binning$sorted.index, weights, y, x$weights, x$rres, x$binning$order)
  
    ## Compute mean and precision.
    XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)

    S <- 0
    tau2 <- get.state(x, "tau2")
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S0[[j]]
    P <- matrix_inv(XWX + S, index = x$sparse.setup)
  
    ## New parameters.
    g <- nu * drop(P %*% crossprod(x$X, x$rres))

    D <- diff(diag(ncol(x$X)))

    d <- 1
    while(d > 0.0001) {
      v <- diag(vfun(g, constr = x$xt$constr))
      g <- drop(matrix_inv(XWX + S + 1e+10 * t(D) %*% v %*% D) %*% crossprod(x$X, x$rres))
      d <- sum((v - diag(vfun(g, constr = x$xt$constr)))^2)
    }
  
    ## Finalize.
    x$state$parameters <- set.par(x$state$parameters, g, "b")
    x$state$fitted.values <- x$fit.fun(x$X, get.state(x, "b"))

    x$state$rss <- sum((x$state$fitted.values - y)^2 * weights)

    if(hatmatrix)
      x$state$hat <- nu * x$X %*% P %*% t(x$X)
  
    return(x$state)
  }

  object$update <- function(x, family, y, eta, id, weights, criterion, ...) {
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
    ng0 <- names(g0)
  
    ## Compute mean and precision.
    XWX <- do.XWX(x$X, 1 / x$weights, x$sparse.setup$matrix)

    D <- diff(diag(ncol(x$X)))

    tau2 <- get.state(x, "tau2")

    S <- 0
    for(j in seq_along(x$S))
      S <- S + 1 / tau2[j] * x$S0[[j]]
    P <- matrix_inv(XWX + S, index = x$sparse.setup)
    g <- drop(P %*% crossprod(x$X, x$rres))

    d <- 1
    while(d > 0.0001) {
      v <- diag(vfun(g, constr = x$xt$constr))
      g <- drop((P <- matrix_inv(XWX + S + 1e+10 * t(D) %*% v %*% D)) %*% crossprod(x$X, x$rres))
      d <- sum((v - diag(vfun(g, constr = x$xt$constr)))^2)
    }
    names(g) <- ng0
  
    ## Compute fitted values.
    if(any(is.na(g)) | any(g %in% c(-Inf, Inf))) {
      x$state$parameters <- set.par(x$state$parameters, rep(0, length(get.state(x, "b"))), "b")
    } else {
      x$state$parameters <- set.par(x$state$parameters, g, "b")
    }
    x$state$fitted.values <- x$fit.fun(x$X, g)

    x$state$edf <- sum_diag(XWX %*% P)

    if(!is.null(x$prior)) {
      if(is.function(x$prior))
        x$state$log.prior <- x$prior(x$state$parameters)
    }
  
    return(x$state)
  }

  K <- object$S[[1]]
  object$S0 <- object$S
  object$S <- list()
  object$S[[1]] <- function(beta, ...) {
    b <- get.par(beta, "b")
    D <- diff(diag(length(b)))
    v <- diag(vfun(b, constr = object$xt$constr))
    return(K + 1e+10 * t(D) %*% v %*% D)
  }
  attr(object$S[[1]], "npar") <- ncol(object$X)

  return(object)
}


WAIC <- function(object, ..., newdata = NULL)
{
  object <- list(object, ...)
  rval <- NULL
  for(j in seq_along(object)) {
    rval <- rbind(rval, compute_WAIC(object[[j]], newdata))
  }
  Call <- match.call()
  row.names(rval) <- if(nrow(rval) > 1) as.character(Call[-1L]) else ""
  rval
}


compute_WAIC <- function(object, newdata = NULL) {
  if(is.null(object$samples)) {
    warning("cannot compute WAIC, return NULL!")
    return(NULL)
  }
  fam <- family(object)
  par <- predict.bamlss(object, newdata = newdata,
    type = "parameter", FUN = function(x) { return(x) }, drop = FALSE)
  y <- if(is.null(newdata)) {
    model.response2(object)
  } else {
    newdata[, response.name(object)]
  }
  if(!is.null(dim(y))) {
    if(ncol(y) < 2L)
      y <- y[, 1L]
  }
  nsamps <- ncol(par[[1L]])
  d <- matrix(NA, nrow(par[[1L]]), nsamps)
  for(i in 1L:nsamps) {
    tpar <- list()
    for(j in names(par))
      tpar[[j]] <- par[[j]][, i]
    d[, i] <- fam$d(y, tpar)
  }
  dm <- rowMeans(d, na.rm = TRUE)
  ldm <- sum(log(dm), na.rm = TRUE)

  a <- log(dm)
  b <- rowMeans(log(d), na.rm = TRUE)

  p1 <- sum(2 * (a - b))
  p2 <- sum(apply(log(d), 1, var))

  waic1 <- -2 * (ldm - p1)
  waic2 <- -2 * (ldm - p2)

  rval <- data.frame(
    "WAIC1" = waic1,
    "WAIC2" = waic2,
    "p1" = p1,
    "p2" = p2
  )

  return(rval)
}


smooth_check <- function(object, newdata = NULL, model = NULL, term = NULL, ...)
{
  if(!inherits(object, "bamlss")) {
    warning("nothing to do!")
    return(NULL)
  }
  samps <- !is.null(object$samples)
  nx <- names(object$x)
  if(!is.null(model))
    nx <- nx[grep2(model, nx, fixed = TRUE)]
  eff <- list()
  for(i in nx) {
    if(!is.null(object$x[[i]]$smooth.construct)) {
      eff[[i]] <- list()
      nt <- names(object$x[[i]]$smooth.construct)
      if(!is.null(term))
        nt <- nt[grep2(term, nt, fixed = TRUE)]
      for(j in nt) {
        p <- as.matrix(predict(object, newdata = newdata, model = i, term = j, intercept = FALSE, FUN = c95, ...))
        minp <- min(p[, "Mean"])
        maxp <- max(p[, "Mean"])
        pp <- (maxp - minp)
        if(samps) {
          se <- as.matrix(samples(object, model = i, term = j))
          secheck <- if(nrow(unique(se)) < 2L) FALSE else TRUE
        } else {
          secheck <- FALSE
        }
        eff[[i]][[j]] <- if(secheck) {
          mean(!((p[, "2.5%"] <= 0) & (p[, "97.5%"] >= 0)) & (pp > 1e-10))
        } else {
          mean(!((p[, "Mean"] <= 1e-20) & (p[, "Mean"] >= -1e-20)))
        }
      }
      eff[[i]] <- unlist(eff[[i]])
    }
  }
  if(length(eff) < 2L)
    eff <- eff[[1L]]
  return(eff)
}


ecdf_transform <- function(data, trans = NULL, notrans = NULL) {
  nt <- is.null(trans)
  if(nt)
    trans <- list()
  for(j in names(data)) {
    if(!(j %in% notrans)) {
      if(nt)
        trans[[j]] <- identity
      if(is.numeric(data[[j]]) | is.integer(data[[j]]) | !nt) {
        if((length(unique(data[[j]])) > 50L) | !nt) {
          if(nt)
            trans[[j]] <- ecdf(data[[j]])
          data[[j]] <- trans[[j]](data[[j]])
        }
      }
    }
  }
  attr(data, "ecdf_transform") <- trans
  return(data)
}

if(FALSE) {
  library("bamlss")
  library("MASS")

  n <- 1000
  x <- sort(runif(n, -pi, pi))
  y <- 2 + 0.5 * x + sin(x) + rnorm(n, sd = 0.3)

  ## P-spline design matrix.
  sm <- smooth.construct(n(~x,k=200,rint=0.1,sint=100,orthc=TRUE), list(x=x), NULL)
  Z <- sm$X
  S <- sm$S[[1]]

  ## Smooth X.
  X <- smooth.construct(s(x), list(x=x), NULL)$X

  ## Orthogonal complement of subspace.
  R <- cbind(X)
  A <- diag(n) - R %*% solve(t(R) %*% R) %*% t(R)
  C <- A %*% Z

  i <- fixDependence(R, C)
  if(!is.null(i)) {
    C <- Z[, -i]
    S <- S[-i, -i]
  }

  ## Centering.
#  Q <- qr.Q(qr(crossprod(C, rep(1, length = nrow(C)))), complete = TRUE)[, -1]
#  C <- C %*% Q
#  K <- crossprod(Q, S) %*% Q
  C <- scale(C)
  K <- diag(1, ncol(C))

  ## Plot basis functions.
  par(mfrow = c(2, 2))
  matplot(x, C, type = "l", lty = 1, col = 1)

  ## Final big design matrix used for estimation.
  G <- cbind(1, X, C)

  ## Estimate coefficients.
  K <- as.matrix(Matrix::bdiag(list(matrix(0, 2, 2), K)))

  GG <- crossprod(G)
  tG <- t(G)

  gcv <- function(lambda) {
    S <- G %*% solve(GG + lambda * K) %*% tG
    yhat <- S %*% y
    trS <- sum(diag(S))
    rss <- sum((y - yhat)^2)
    drop(rss * n / (n - trS)^2)
  }

  lambda <- optimize(gcv, lower = 1e-20, upper = 1e+10)$minimum

  beta <- solve(t(G) %*% G + lambda * K) %*% t(G) %*% y
  fit <- G %*% beta
  fitl <- G[, 1:2] %*% beta[1:2]
  fits <- G[, -c(1:2)] %*% beta[-c(1:2)]

  plot(x, y, main = "fit", ylim = range(c(y, fitl, fits)))
  lines(fit ~ x, col = 2, lwd = 2)
  lines(fitl ~ x, col = 4, lwd = 2)
  lines(fits ~ x, col = 4, lwd = 2)

  plot(x, fitl, type = "l", main = "linear", ylim = range(c(0.5*x, fitl)))
  lines(2 + 0.5*x ~ x, col = 2)

  plot(x, fits, type = "l", main = "smooth", ylim = range(c(sin(x), fits)))
  lines(sin(x) ~ x, col = 2)
}


.engines <- function(family)
{
  family <- bamlss.family(family)
  bayesx <- !is.null(family$bayesx)
  jags <- !is.null(family$bugs)
  optimizer <- !is.null(family$optimizer)
  sampler <- !is.null(family$sampler)
  tab <- c(
    "family" = family$family,
    "opt_bfit()" = !optimizer,
    "opt_boost()" = !optimizer,
    "opt_bbfit()" = !optimizer,
    "sam_GMCMC()" = !sampler,
    "sam_BayesX()" = bayesx,
    "sam_JAGS()" = jags,
    "special_opt()" = optimizer,
    "special_sam()" = sampler
  )
  return(tab)
}

engines <- function(family, ...)
{
  if(missing(family)) {
    family <- ls("package:bamlss")
    family <- grep("_bamlss", family, fixed = TRUE, value = TRUE)
    fam <- NULL
    for(f in family) {
      foo <- get(f)
      foo <- try(foo(), silent = TRUE)
      if(!inherits(foo, "try-error")) {
        if(inherits(foo, "family.bamlss"))
          fam <- c(fam, f)
      }
    }
    family <- fam
  }
  family <- list(family, ...)
  tab <- list()
  call <- match.call()
  fn <- as.character(call)[-1L]
  for(i in seq_along(family)) {
    tab[[i]] <- .engines(family[[i]])
  }
  tab <- do.call("rbind", tab)
  nt <- colnames(tab)
  rn <- tab[, 1]
  tab <- t(tab[, -1])
  mode(tab) <- "logical"
  rownames(tab) <- nt[-1]
  colnames(tab) <- rn
  return(as.data.frame(tab))
}


CRPS <- function(object, newdata = NULL, interval = c(-Inf, Inf), FUN = mean, term = NULL, ...) {
  yname <- response_name(object)
  fam <- family(object)
  if(!is.null(fam$type)) {
    if(tolower(fam$type) != "continuous")
      stop("CRPS only for continuous responses!")
  }
  if(is.null(fam$p))
    stop("no p() function in family object!")
  if(is.null(newdata))
    newdata <- model.frame(object)
  par <- as.data.frame(predict(object, newdata = newdata, type = "parameter", drop = FALSE))
  if(!is.null(fam$valid.response)) {
    vd <- rep(NA, 2)
    ty <- c(-0.0001, 0.0001)
    for(i in seq_along(ty)) {
      vd[i] <- fam$valid.response(ty[i])
    }
    if(!vd[1L])
      interval[1L] <- 1e-20
    if(!vd[2L])
      interval[2L] <- -1e-20
  }
  crps <- if(is.null(fam$crps)) {
    .CRPS(newdata[[yname]], par, fam, interval)
  } else {
    fam$crps(newdata[[yname]], par)
  }
  return(FUN(crps, ...))
}

.CRPS <- function(y, par, family, interval = c(-Inf, Inf)) {
  if(is.function(family))
    family <- family()
  if(inherits(family, "gamlss.family"))
    family <- tF(family)
  family <- complete.bamlss.family(family)
  n <- length(y)
  crps <- rep(0, n)
  for(i in 1:n) {
    foo1 <- function(x) {
      family$p(x, par[i, , drop = FALSE])^2
    }
    foo2 <- function(x) {
      (family$p(x, par[i, , drop = FALSE]) - 1)^2
    }
    int1 <- integrate(foo1, lower = interval[1L], upper = y[i])$value
    int2 <- integrate(foo2, lower = y[i], upper = interval[2L])$value
    crps[i] <- int1 + int2
  }
  return(crps)
}


GramSchmidt <- function(X) {
  itcpt <- if(is.null(colnames(X))) {
    FALSE
  } else {
    any(colnames(X) == "(Intercept)")
  }
  js <- if(itcpt) 2 else 1
  tX <- X
  coef <- list()
  for(j in js:ncol(X)) {
    b <- lm.fit(tX[, 1:(j - 1), drop = FALSE], X[, j])
    coef[[paste0("c", j)]] <- b$coefficients
    tX[, j] <- b$residuals
  }
  return(list("coefficients" = coef, "X" = tX))
}

