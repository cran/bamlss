######################
## BayesX interface ##
######################
BayesX.control <- function(n.iter = 1200, thin = 1, burnin = 200,
  seed = NULL, predict = "light", model.name = "bamlss", data.name = "d",
  prg.name = NULL, dir = NULL, verbose = FALSE, show.prg = TRUE, modeonly = FALSE, ...)
{
  if(is.null(seed))
    seed <- '##seed##'
  stopifnot(burnin < n.iter)
  if(is.null(model.name))
    model.name <- 'bamlss'
  if(is.null(data.name))
    data.name <- 'd'
  if(is.null(prg.name))
    prg.name <- paste(model.name, 'prg', sep = '.')
  if(!grepl(".prg", prg.name))
    prg.name <- paste(prg.name, "prg", sep = ".")
  if(is.null(dir)) {
    dir.create(dir <- tempfile())
    attr(dir, "unlink") <- TRUE
  } else dir <- path.expand(dir)
  if(!file.exists(dir)) dir.create(dir)
  cores <- 1

  cvals <- list(
    "prg" = list(
      "iterations" = n.iter, "burnin" = burnin, "step" = thin,
      "setseed" = seed, "predict" = predict, "modeonly" = modeonly
    ),
    "setup" = list(
      "main" = c(rep(FALSE, 3), rep(TRUE, 2)), "model.name" = model.name, "data.name" = data.name,
      "prg.name" = prg.name, "dir" = dir, "verbose" = verbose, "show.prg" = show.prg, "cores" = cores
    )
  )

  cvals
}


sam_BayesX <- BayesX <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
  data = NULL, control = BayesX.control(...), ...)
{
  cmd <- paste0('stopifnot(require', 'Name', 'space("Bayes', 'X', 'src"))')
  eval(parse(text = cmd))

  if(is.null(family$bayesx))
    stop("BayesX specifications missing in family object, cannot set up model!")

  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, update = update, parametric2smooth = FALSE, nodf = TRUE, ...)

  if(!is.null(start))
    x <- set.starting.values(x, start) ## FIXME: starting values for parametric parts?

  fbx <- family$bayesx

  model.name <- control$setup$model.name
  data.name <- rmf(control$setup$data.name)
  prg.name <- control$setup$prg.name
  dir <- control$setup$dir
  cores <- control$setup$cores
  if(is.null(cores)) cores <- 1

  if(!file.exists(dir)) {
    dir.create(dir)
    on.exit(unlink(dir))
  }

  if(is.null(data)) {
    data <- try(get('bf', envir = parent.frame())$model.frame, silent = TRUE)
    if(inherits(data, "try-error"))
      stop("cannot find the model.frame for creating BayesX model term objects!")
  }
  if(is.null(data))
    stop("no data available for creating BayesX model term objects!")

  for(j in names(data)) {
    if(is.factor(data[[j]])) {
      data <- cbind(data, as.data.frame(model.matrix(as.formula(paste("~ -1 +", j)), data = data)))
    }
  }

  if(is.data.frame(y)) {
    if(ncol(y) < 2) {
      y <- y[, 1, drop = FALSE]
    }
  }

  cny <- colnames(y)
  for(j in seq_along(cny))
    cny[j] <- all.vars(as.formula(paste("~", cny[j])))
  colnames(y) <- cny
  yname <- colnames(y)[1]

  for(j in names(y)) {
    if(is.factor(y[[j]])) {
      if(nlevels(y[[j]]) < 3) {
        y[[j]] <- as.integer(as.integer(y[[j]]) - 1L)
      }
    }
  }

  is.user <- function(x) {
    iu <- inherits(x, "userdefined.smooth.spec") | inherits(x, "tensorX.smooth") | inherits(x, "userdefined.smooth") | inherits(x, "tensor.smooth")
    if(!is.null(x$sx.construct))
      iu <- FALSE
    if(is.null(x$sx.construct) & !iu)
      iu <- TRUE
    return(iu)
  }

  is.tx <- function(x) {
    inherits(x, "tensorX.smooth") | inherits(x, "tensorX3.smooth")
  }

  single_eqn <- function(x, y, id) {
    rhs <- dfiles <- prgex <- sdata <- NULL

    if(!is.null(x$model.matrix)) {
      cn <- rmf(colnames(x$model.matrix))
      colnames(x$model.matrix) <- cn
      cn <- cn[cn != "Intercept"]
      if(length(cn)) {
        rhs <- c(rhs, cn)
        sdata <- as.data.frame(x$model.matrix[, cn, drop = FALSE])
      }
      if("Intercept" %in% colnames(x$model.matrix)) {
        sdata <- cbind("Intercept" = rep(1, if(is.null(dim(y))) length(y) else nrow(y)), sdata)
        rhs <- c("const", rhs)
      }
    }

    if(!is.null(sdata))
      sdata <- as.data.frame(sdata)

    if(!is.null(x$smooth.construct)) {
      for(j in names(x$smooth.construct)) {
        if(is.null(x$smooth.construct[[j]]$sx.construct))
          class(x$smooth.construct[[j]]) <- c("userdefined.smooth.spec", class(x$smooth.construct[[j]]))
        sxc <- sx.construct(x$smooth.construct[[j]], data, id = c(id, j), dir = dir, mcmcreg = TRUE)
        if(!is.null(attr(sxc, "write")))
          prgex <- c(prgex, attr(sxc, "write")(dir))
        rhs <- c(rhs, sxc)
        tl <- if(is.null(x$smooth.construct[[j]]$tx.term)) {
          x$smooth.construct[[j]]$term
        } else {
          x$smooth.construct[[j]]$tx.term
        }
        otl <- x$smooth.construct[[j]]$term
        if((x$smooth.construct[[j]]$by != "NA") & is.user(x$smooth.construct[[j]])) {
          tl <- c(tl, x$smooth.construct[[j]]$by)
          otl <- c(otl, x$smooth.construct[[j]]$by)
        }
        if((length(tl) > 1) & is.user(x$smooth.construct[[j]]) & !is.tx(x$smooth.construct[[j]]))
          tl <- paste(tl, collapse = "")
        if(is.null(sdata)) {
          if(is.user(x$smooth.construct[[j]])) {
            sdata <- match.index(data[, otl, drop = FALSE])$match.index
            colnames(sdata) <- tl
          } else {
            sdata <- data[, tl, drop = FALSE]
          }
        }
        if(!all(tl %in% colnames(sdata))) {
          for(tlj in tl) {
            if(tlj %in% names(data)) {
              sdata[[tlj]] <- data[[tlj]]
            } else {
              sdata[[tlj]] <- match.index(data[, otl, drop = FALSE])$match.index
            }
          }
        }
        if(x$smooth.construct[[j]]$by != "NA") {
          if(is.factor(data[[x$smooth.construct[[j]]$by]])) {
            mm <- model.matrix(as.formula(paste("~ -1 +", x$smooth.construct[[j]]$by)), data = data)
            for(tlj in colnames(mm))
              sdata[[tlj]] <- mm[, tlj]
          } else {
            sdata[[x$smooth.construct[[j]]$by]] <- data[[x$smooth.construct[[j]]$by]]
          }
        }
      }
    }

    rn <- response.name(as.formula(x$formula), hierarchical = FALSE)

    if(!family$family == "dirichlet") {
      if(is.null(family$cat)) {
        if(length(cny) > 1) {
          for(ii in 0:9) {
            if(grepl(ii, rn))
              rn <- grep(ii, cny, value = TRUE)
          }
          if(rn %in% family$names)
            rn <- yname
        } else {
          if(rn %in% family$names)
            rn <- NA
          if(is.na(rn))
            rn <- yname
        }
      }
    }

    eqn <- paste(rn, "=", paste(rhs, collapse = " + "))
    rval <- list("eqn" = eqn, "prgex" = prgex)

    if(!is.null(sdata)) {
      for(j in names(sdata)) {
        if(is.factor(sdata[[j]]))
          sdata[[j]] <- as.integer(as.character(sdata[[j]]))
      }
      if(nrow(sdata) == (if(is.null(dim(y))) length(y) else nrow(y)))
        sdata <- cbind(sdata, y)
      rval$dname <- paste(paste(id, collapse = "_"), data.name, sep = "_")
      write.table(sdata, file = file.path(dir, paste(rval$dname, ".raw", sep = "")),
        quote = FALSE, row.names = FALSE)
      rval$prgex <- c(
        paste("dataset", rval$dname),
        paste(rval$dname, ".infile using ",
          file.path(dir, paste(rval$dname, ".raw", sep = "")), sep = ""),
        rval$prgex
      )
    }

    rval
  }

  eqn <- list()
  prgex <- NULL
  n <- 1
  main <- if(is.null(family$bayesx$main)) c(TRUE, rep(FALSE, length(x) - 1)) else family$bayesx$main
  pcmd <- control$prg$predict
  control$prg$predict <- NULL
  control$prg$quantile <- family$bayesx$quantile
  modeonly <- control$prg$modeonly
  control$prg$modeonly <- NULL
  if(modeonly)
    control$prg <- control$prg[!(names(control$prg) %in% c("iterations", "burnin", "step"))]

  nrcat <- attr(family$bayesx, "nrcat")

  for(i in names(x)) {
    if(!all(c("fake.formula", "formula") %in% names(x[[i]]))) {
      stop("hierarchical models not supported yet!")
      eqn[[i]] <- list()
      k <- 1
      for(j in names(x[[i]])) {
        msp <- single_eqn(x[[i]][[j]], y, id = c(i, j))
        teqn <- paste(model.name, ".hregress ", msp$eqn,
          ", family=", if(k < 2) fbx[[i]][1] else "gaussian_re",
          " equationtype=", fbx[[i]][2],
          if(n == length(x) & k < 2) {
            paste(" ", paste(names(control$prg), "=", control$prg, sep = "", collapse = " "))
          } else NULL,
          if(main[n]) {
            paste(" predict=", pcmd, sep = "")
          } else NULL,
          if(!is.null(msp$dname)) paste(" using", msp$dname) else NULL, sep = "")
        eqn[[i]][[j]] <- teqn
        prgex <- c(prgex, msp$prgex)
        k <- k + 1
      }
      eqn[[i]] <- rev(eqn[[i]])
    } else {
      msp <- single_eqn(x[[i]], y, id = i)
      teqn <- paste(model.name, ".hregress ", msp$eqn, ", family=", fbx[[i]][1],
        if(!is.null(nrcat)) paste0(if(n == 1) " hlevel=1 " else NULL, " nrcat=", nrcat) else NULL,
        " equationtype=", if(!main[n]) fbx[[i]][length(fbx[[i]])] else fbx[[i]][2],
        if(n == length(x)) {
          paste(paste(" ", paste(names(control$prg), "=", control$prg, sep = "", collapse = " ")),
            if(modeonly) " modeonly" else "", sep = "")
        } else NULL,
        if(main[n]) {
          if(modeonly) " modeonly" else paste(" predict=", pcmd, sep = "")
        } else NULL,
        if(!is.null(msp$dname)) paste(" using", msp$dname) else NULL, sep = "")
      eqn[[i]] <- teqn
      prgex <- c(prgex, msp$prgex)
    }
    n <- n + 1
  }

  prg <- c(prgex, "", paste("mcmcreg", model.name), "")
  for(i in unlist(rev(eqn)))
    prg <- c(prg, i, "")

  if(!modeonly)
    prg <- c(prg, paste(model.name, "getsample", sep = "."))
  prg <- c(
    paste('%% BayesX program created by bamlss: ', as.character(Sys.time()), sep = ''),
    paste('%% usefile ', file.path(dir, prg.name), sep = ''), "",
    prg
  )
  prg <- unique(prg)

  if(any(grepl("##seed##", prg, fixed = TRUE)))
    prg <- gsub("##seed##", round(runif(1L) * .Machine$integer.max), prg, fixed = TRUE)

  if(control$setup$show.prg)
    writeLines(prg)

  prgf <- file.path(dir, prg.name)
  writeLines(prg, prgf)

  warn <- getOption("warn")
  options(warn = -1)
  cmd <- paste0("Bayes", "X", "src", rep(":", 2),
    "run.bayesx(prg = prgf, verbose = control$setup$verbose)")
  ok <- eval(parse(text = cmd))
  options("warn" = warn)
  if(length(i <- grep("error", ok$log, ignore.case = TRUE))) {
    errl <- gsub("^ +", "", ok$log[i])
    errl <- gsub(" +$", "", errl)
    errl <- encodeString(errl, width = NA, justify = "left")
    errl <- paste(" *", errl)
    errm <- paste("an error occurred running the BayesX binary! The following messages are returned:\n",
      paste(errl, collapse = "\n", sep = ""), sep = "")
    warning(errm, call. = FALSE)
  }
  if(length(i <- grep("-nan", ok$log, ignore.case = TRUE))) {
    warning("the BayesX engine returned NA samples, please check your model specification! In some cases it can be helpful to center continuous covariates!", call. = FALSE)
  }

  sfiles <- grep("_sample.raw", dir(file.path(dir, "output")), fixed = TRUE, value = TRUE)

  if(!length(sfiles)) {
    warning("BayesX did not return any samples!")
    return(NULL)
  }

  samples <- NULL
  for(i in names(x)) {
    if(!all(c("fake.formula", "formula") %in% names(x[[i]]))) {
      stop("hierarchical models not supported yet!")
    } else {
      if(!is.null(x[[i]]$model.matrix)) {
        sf <- grepl(paste("_", i, "_", sep = ""), sfiles, fixed = TRUE) & grepl(paste("LinearEffects", sep = ""), sfiles, fixed = TRUE)
        if(!any(sf) & (length(cny) > 1)) {
          for(ii in 0:9) {
            if(grepl(ii, i) & any(grepl(ii, cny))) {
              jj <- gsub(ii, "", i)
              ii <- grep(ii, cny, value = TRUE)
              break
            }
          }
          sf <- grepl(paste("_", jj, "_", sep = ""), sfiles, fixed = TRUE) & grepl(paste(ii, "_LinearEffects", sep = ""), sfiles, fixed = TRUE)
        }
        sf <- sfiles[sf]
        if(length(sf)) {
          samps <- as.matrix(read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE])
          colnames(samps) <- paste(i, ".p.", colnames(x[[i]]$model.matrix), sep = "")
          samples <- cbind(samples, samps)
        }
      }
      if(!is.null(x[[i]]$smooth.construct)) {
        for(j in seq_along(x[[i]]$smooth.construct)) {
          tl <- if(is.null(x[[i]]$smooth.construct[[j]]$tx.term)) {
            x[[i]]$smooth.construct[[j]]$term
          } else {
            x[[i]]$smooth.construct[[j]]$tx.term
          }
          if((x[[i]]$smooth.construct[[j]]$by != "NA") & is.user(x[[i]]$smooth.construct[[j]])) {
            tl <- if(is.tx(x[[i]]$smooth.construct[[j]])) {
              c(x[[i]]$smooth.construct[[j]]$by, tl)
            } else {
              c(tl, x[[i]]$smooth.construct[[j]]$by)
            }
          }
          if((length(tl) > 1) & is.user(x[[i]]$smooth.construct[[j]]) & !is.tx(x[[i]]$smooth.construct[[j]]))
            tl <- paste(tl, collapse = "")
          term <- paste(if(is.user(x[[i]]$smooth.construct[[j]])) "of" else NULL, "_", paste(tl, collapse = "_"), "_", "sample", sep = "")
          #term <- paste("of", term, "sample", sep = "")
          sf <- grepl(paste("_", i, "_", sep = ""), sfiles, fixed = TRUE) & grepl(term, sfiles, fixed = TRUE) & !grepl("_variance_", sfiles, fixed = TRUE)
          sf <- sfiles[sf]

          if((length(sf) < 1) & (length(cny) > 1)) {
            for(ii in 0:9) {
              if(grepl(ii, i) & any(grepl(ii, cny))) {
                jj <- gsub(ii, "", i)
                ii <- grep(ii, cny, value = TRUE)
                break
              }
            }
            sf <- grepl(paste("_", jj, "_", sep = ""), sfiles, fixed = TRUE) & grepl(term, sfiles, fixed = TRUE) & !grepl("_variance_", sfiles, fixed = TRUE) & grepl(paste0("_", ii, "_"), sfiles, fixed = TRUE)
            sf <- sfiles[sf]
          }

          if(inherits(x[[i]]$smooth.construct[[j]], "mrf.smooth")) {
            if(!inherits(x[[i]]$smooth.construct[[j]], "mgcv.smooth"))
              sf <- grep("_spatial_", sf, fixed = TRUE, value = TRUE)
          }
          if(inherits(x[[i]]$smooth.construct[[j]], "random.effect")) {
            if(!inherits(x[[i]]$smooth.construct[[j]], "mgcv.smooth"))
              sf <- grep("_random_", sf, fixed = TRUE, value = TRUE)
          }
          if(any(grepl("anisotropy", sf)))
            sf <- sf[-grep("anisotropy", sf)]
          tj <- grep("tensor", sf, fixed = TRUE)
          if(length(tj))
            sf <- if(is.tx(x[[i]]$smooth.construct[[j]]) & (length(x[[i]]$smooth.construct[[j]]$term) > 1)) sf[tj] else sf[-tj]
          if((x[[i]]$smooth.construct[[j]]$by != "NA") & !is.user(x[[i]]$smooth.construct[[j]]))
            sf <- grep(x[[i]]$smooth.construct[[j]]$by, sf, fixed = TRUE, value = TRUE)

          if(length(sf) < 1)
            next

          samps <- as.matrix(read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE])
          cn <- colnames(x[[i]]$smooth.construct[[j]]$X)
          if(is.null(cn))
            cn <- paste("b", 1:ncol(x[[i]]$smooth.construct[[j]]$X), sep = "")
          if(length(cn) != ncol(samps)) {
            cn <- paste("b", 1:ncol(samps), sep = "")
            warning(paste("number of returned parameters from BayesX is different from number of columns in the design matrix for term ",
              names(x[[i]]$smooth.construct)[j], "!", sep = ""))
          }
          colnames(samps) <- paste(i, ".s.", x[[i]]$smooth.construct[[j]]$label, ".", cn, sep = "")
          sf <- grepl(paste("_", i, "_", sep = ""), sfiles, fixed = TRUE) & grepl(term, sfiles, fixed = TRUE) & grepl("_variance_", sfiles, fixed = TRUE)
          sf <- sfiles[sf]
          term <- gsub("_sample", "_omega", term, fixed = TRUE)
          aniso <- grepl(paste("_", i, "_", sep = ""), sfiles, fixed = TRUE) & grepl(term, sfiles, fixed = TRUE) & grepl("_anisotropy_", sfiles, fixed = TRUE)
          aniso <- sfiles[aniso]
          if(length(aniso) & FALSE)
            sf <- aniso
          if(length(sf)) {
            tj <- grep("tensor", sf, fixed = TRUE)
            if(length(tj))
              sf <- if(is.tx(x[[i]]$smooth.construct[[j]]) & (length(x[[i]]$smooth.construct[[j]]$term) > 1)) sf[tj] else sf[-tj]
            if(inherits(x[[i]]$smooth.construct[[j]], "mrf.smooth")) {
              if(!inherits(x[[i]]$smooth.construct[[j]], "mgcv.smooth"))
                sf <- grep("_spatial_", sf, fixed = TRUE, value = TRUE)
            }
            if(inherits(x[[i]]$smooth.construct[[j]], "random.effect")) {
              if(!inherits(x[[i]]$smooth.construct[[j]], "mgcv.smooth"))
                sf <- grep("_random_", sf, fixed = TRUE, value = TRUE)
            }
            if((x[[i]]$smooth.construct[[j]]$by != "NA") & !is.user(x[[i]]$smooth.construct[[j]]))
              sf <- grep(x[[i]]$smooth.construct[[j]]$by, sf, fixed = TRUE, value = TRUE)
            vsamps <- as.matrix(read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE])
            colnames(vsamps) <- paste(i, ".s.", x[[i]]$smooth.construct[[j]]$label, ".", paste("tau2", 1:ncol(vsamps), sep = ""), sep = "")
            samps <- cbind(samps, vsamps)
          }
          samples <- cbind(samples, samps)
        }
      }
    }
  }

  if(FALSE) {
    sf <- grep("_DIC", dir(file.path(dir, "output")), fixed = TRUE, value = TRUE)
    dic <- read.table(file.path(dir, "output", sf), header = TRUE)[, -1, drop = FALSE]
    samples <- cbind(samples, "DIC" = dic$dic, "pd" = dic$pd)
  }

  as.mcmc(samples)
}


########################################
## (2) BayesX model term construction ##
########################################
sx <- function(x, z = NULL, bs = "ps", by = NA, ...)
{
  by <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  call <- match.call()

  available.terms <- c(
    "rw1", "rw2",
    "season",
    "ps", "psplinerw1", "psplinerw2", "pspline",
    "te", "pspline2dimrw2", "te1", "pspline2dimrw1",
    "kr", "kriging",
    "gk", "geokriging",
    "gs", "geospline",
    "mrf", "spatial",
    "bl", "baseline",
    "factor",    
    "ridge", "lasso", "nigmix",
    "re", "ra", "random",
    "cs", "catspecific",
    "offset",
    "generic",
    "rps", "hrandom_pspline"
  )
  if(!bs %in% available.terms) stop(paste("basis type", sQuote(bs), "not supported by BayesX"))

  if(bs %in% c("rsps", "hrandom_pspline")) {
    bs <- "rsps"
    x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
    rcall <- paste("r(x = ", by, ", bs = ", sQuote(bs), ", by = ", x, ", ...)", sep = "")
    rval <- eval(parse(text = rcall))
  } else {
    if(length(grep("~", term <- deparse(call$x))) && bs %in% c("re", "ra", "random")) {
      x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
      rcall <- paste("r(x = ", x, ", by = ", by, ", ...)", sep = "")
      rval <- eval(parse(text = rcall))
    } else {
      k <- -1
      m <- NA
      xt <- list(...)
      if("xt" %in% names(xt))
        xt <- xt[["xt"]]
      if(is.null(xt$lambda))
        xt$lambda <- 100
      if("m" %in% names(xt))
        stop("argument m is not allowed, please see function s() using this specification!")
      if("k" %in% names(xt))
        stop("argument k is not allowed, please see function s() using this specification!")
      if(!is.null(xt$xt))
        xt <- xt$xt
      warn <- getOption("warn")
      options("warn" = -1)
      if(by != "NA" && is.vector(by) && length(by) < 2L && !is.na(as.numeric(by))) {
        xt["b"] <- by
        by <- "NA"
      }
      options("warn" = warn)
      if(bs %in% c("pspline2dimrw1", "pspline2dimrw2", "te",
        "gs", "geospline", "kr", "gk", "kriging", "geokriging")) {
        if(by != "NA")
          stop(paste("by variables are not allowed for smooths of type bs = '", bs, "'!", sep = ""))
      }
      if(bs == "te1") bs <- "pspline2dimrw1"
      if(!is.null(xt$knots))
        xt$nrknots <- xt$knots
      if(bs %in% c("ps", "te", "psplinerw1", "psplinerw2", "pspline",
        "pspline2dimrw2", "pspline2dimrw1", "gs", "geospline")) {
        if(!is.null(xt$degree))
          m <- xt$degree
        if(!is.null(xt$order)) {
          if(is.na(m))
            m <- c(3L, xt$order)
          else
            m <- c(m[1L], xt$order)
        }
        if(length(m) < 2L && (bs %in% c("gs", "geospline")))
          m <- c(3L, 1L)
        if(length(m) < 2L && is.na(m))
          m <- c(3L, 2L)   
        if(is.null(xt$order) && length(m) < 2L)
          m <- c(m, 2L)
        if(is.null(xt$order) && length(m) < 2L)
          m <- c(m, 1L)
        m[1L] <- m[1L] - 1L
        if(!is.null(xt$nrknots))
          k <- xt$nrknots + m[1L]
        else {
          if(bs %in% c("ps", "psplinerw1", "psplinerw2", "pspline"))
            k <- 20L + m[1L]
          else
            k <- 10L + m[1L]
        }
      }
      if(bs %in% c("kr", "gk", "kriging", "geokriging")) {
        m <- c(1L, 1L)
        if(!is.null(xt$nrknots)) {
          k <- xt$nrknots
        } else k <- -1L
      }
      if(!is.null(xt$map))
        xt$map.name <- rmf(as.character(call$map))
      xt[c("degree", "order", "knots", "nrknots")] <- NULL
      if(!length(xt))
        xt <- NULL
      if(!is.null(call$z)) 
        term <- c(term, deparse(call$z))
      if(bs != "te") {
        rval <- s(x, z, k = k, bs = bs, m = m, xt = xt)
      } else {
        rval <- te(x, z, k = k, bs = "ps", m = m, xt = xt, mp = FALSE)
        rval$margin[[1]]$term <- term[1]
        rval$margin[[2]]$term <- term[1]
      }
      rval$by <- by
      rval$term <- term
      rval$dim <- length(term)
      rval$label <- paste("sx(", paste(term, collapse = ",", sep = ""),
        if(by != "NA") paste(",by=", by, sep = "") else NULL, ")", sep = "")
    }
  }

  if(is.null(rval$xt$prior))
    rval$xt$prior <- "ig"
  if(is.null(rval$xt$theta)) {
    rval$xt$theta <- switch(rval$xt$prior,
      "sd" = 0.00877812,
      "hc" = 0.01034553,
      "hn" = 0.1457644,
      "u" = 0.2723532
    )
  }
  if(is.null(rval$xt$scaletau2))
    rval$xt$scaletau2 <- rval$xt$theta
  if(is.null(rval$xt$hyperprior)) {
    rval$xt$hyperprior <- switch(rval$xt$prior,
      "ig" = "invgamma",
      "hn" = "hnormal",
      "sd" = "scaledep",
      "hc" = "hcauchy",
      "u" = "aunif"
    )
  }
  rval$xt$prior <- NULL
  rval$xt$theta <- NULL

  rval$special <- TRUE
  rval$sx.construct <- TRUE
  class(rval) <- c(class(rval), "no.mgcv")

  return(rval)
}

sx.construct <- function(object, data, ...)
{
  UseMethod("sx.construct")
}

sx.construct.default <- function(object, data, ...) 
{
  cl <- grep(".smooth.spec", class(object), value = TRUE, fixed = TRUE)
  bs <- gsub(".smooth.spec", "", cl, fixed = TRUE)
  if(length(cl)) {
    info <- if(cl == paste(bs, "smooth.spec", sep = ".")) {
      paste("with basis", sQuote(bs))
    } else {
      paste("of class", sQuote(cl))
    }
  } else info <- paste("of class", paste(sQuote(class(object)), collapse = ", "))
  stop(paste("BayesX does not support smooth terms ", info,
    ", it is recommended to use sx() for specifying smooth terms",
    sep = ""))
}

do.xt <- function(term, object, not = NULL, noco = FALSE)
{
  if(!is.null(object$xt)) {
    names.xt <- names(object$xt)
    if(is.null(not))
      not <- "not"
    count <- 1
    co <- ","
    if(noco)
      co <- NULL
    for(name in names.xt) {
      if(count > 1)
        co <- ","
      if(!name %in% not) {
        if(name %in% c("full", "catspecific", "center", "derivative", "nofixed") || 
          is.logical(object$xt[[name]])) {
          if(is.logical(object$xt[[name]])) {
            if(object$xt[[name]])
              term <- paste(term, co, name, sep = "")
          } else term <- paste(term, co, name, sep = "")
          count <- count + 1
        } else {
          term <- paste(term, co, name, "=", object$xt[name], sep = "")
          count <- count + 1
        }
      }
    }
  }

  return(term)
}

sx.construct.userdefined.smooth.spec <- sx.construct.tensorX.smooth <- function(object, data, id = NULL, dir = NULL, ...)
{
  if(is.null(object$xt$hyperprior)) {
    if(is.null(object$xt$prior))
      object$xt$prior <- "ig"
  } else {
    object$xt$prior <- switch(object$xt$hyperprior,
       "invgamma" = "ig",
       "hnormal" = "hn",
       "scaledep" = "sd",
       "hcauchy" = "hc",
       "aunif" = "u"
    )
  }
  if(is.null(object$xt$scaletau2)) {
    if(is.null(object$xt$theta)) {
      object$xt$theta <- switch(object$xt$prior,
        "sd" = 0.00877812,
        "hc" = 0.01034553,
        "hn" = 0.1457644,
        "u" = 0.2723532
      )
    }
  } else {
    object$xt$theta <- object$xt$scaletau2
  }
  object$xt$scaletau2 <- object$xt$theta
  object$xt$hyperprior <- switch(object$xt$prior,
    "ig" = "invgamma",
    "hn" = "hnormal",
    "sd" = "scaledep",
    "hc" = "hcauchy",
    "u" = "aunif"
  )
  if(!is.null(object$xt[["pSigma"]])) {
    cat("yes!\n")
    object$S <- object$sx.S <- object$xt[["pSigma"]]
  }
  object$state <- NULL
  if(!is.null(object$sx.S))
    object$S <- object$sx.S
  if(!is.null(object$sx.rank))
    object$rank <- object$sx.rank
  is.tx <- inherits(object, "tensorX.smooth") | inherits(object, "tensorX3.smooth")
  if(inherits(object, "tensor.smooth")) {
    S <- 0
    for(j in seq_along(object$S))
      S <- S + object$S[[j]]
    object$S <- list(S)
  }
  if(!is.null(object$xt$doC))
    object$C <- Cmat(object)
  if(!is.null(object$C)) {
    if(nrow(object$C) < 1)
      object$C <- NULL
  }
  if(is.null(object$C) & !is.tx)
    object$xt$nocenter <- TRUE
  if(is.null(id))
    id <- "t"
  id <- paste(rmf(id), collapse = "_")
  term <- if(length(object$term) > 1) {
    paste(if(is.null(object$tx.term)) object$term else object$tx.term, collapse = if(!is.tx) "" else "*")
  } else object$term
  by <- if(object$by != "NA") object$by else NULL
  if(!is.null(by)) {
    term <- if(is.tx) paste(by, term, sep = "*") else paste(term, by, sep = "")
    Sn <- paste(id, by, "S", sep = "_")
    Sn <- paste(Sn, "", 1:length(object$S), sep = "")
    Xn <- paste(id, by, "X", sep = "_")
    Cn <- paste(id, by, "C", sep = "_")
    Pn <- paste(id, by, "P", sep = "_")
    Pm <- paste(id, by, "Pm", sep = "_")
  } else {
    Sn <- paste(id, "S", sep = "_")
    Sn <- paste(Sn, "", 1:length(object$S), sep = "")
    Xn <- paste(id, "X", sep = "_")
    Cn <- paste(id, "C", sep = "_")
    Pn <- paste(id, "P", sep = "_")
    Pm <- paste(id, "Pm", sep = "_")
  }
  if(is.null(object$rank))
    object$rank <- sapply(object$S, function(x) { qr(x)$rank })
  if(!is.null(object$xt$nocenter))
    object$xt$centermethod <- NULL
  if(!is.null(object$C)) {
    object$xt$centermethod <- NULL
    object$xt$nocenter <- NULL
  }
  term <- paste(term, if(is.tx & (length(object$S) > 1)) "(tensor," else "(userdefined,", sep = "")
  for(j in seq_along(object$S))
    term <- paste(term, paste("penmatdata", if(j < 2) "" else j, "=", sep = ""), Sn[j], ",", sep = "")
  noc_and_cm <- is.null(object$xt$nocenter) & is.null(object$xt$centermethod)
  if(length(object$S) > 1) {
    Xn <- paste(Xn, "", 1:length(object$S), sep = "")
    for(j in seq_along(object$S)) {
      term <- paste(term, paste("designmatdata", if(j < 2) "" else j, "=", sep = ""),
        Xn[j], if(j == length(object$S)) { if(noc_and_cm) "," else "" } else ",", sep = "")
    }
  } else {
    term <- paste(term, "designmatdata=", Xn,
      if(noc_and_cm) "," else "",
      sep = "")
  }
  if(!is.null(object$C))
    term <- paste(term, "constrmatdata=", Cn, sep = "")
  if(!is.null(object$state$parameters))
    term <- paste(term, ",betastart=", Pn, sep = "")
  if(!is.null(object$xt[["pmean"]]))
    term <- paste(term, ",priormeandata=", Pm, sep = "")
  if(is.null(object$xt$nocenter) & is.null(object$xt$centermethod) & !is.null(object$rank))
    term <- paste(term, ",rankK=", sum(object$rank), sep = "")
  term <- paste(do.xt(term, object,
    c("center", "before", "penalty", "polys", "map", "map.name", "nb", "gra", "ft", "prior", "theta", "pmean", "pSigma", "doC", "constraint", "binning")), ")", sep = "")

  write <- function(dir) {
    exists <- NULL
    for(j in seq_along(object$S)) {
      if(!file.exists(file.path(dir, paste(Sn[j], ".raw", sep = "")))) {
        colnames(object$S[[j]]) <- rownames(object$S[[j]]) <- NULL
        write.table(round(object$S[[j]], 5), file = file.path(dir, paste(Sn[j], ".raw", sep = "")),
          quote = FALSE, row.names = FALSE)
      } else exists <- c(exists, file.path(dir, paste(Sn[j], ".raw", sep = "")))
    }
    if(is.tx) {
      for(j in seq_along(object$S)) {
        if(!file.exists(file.path(dir, paste(Xn[j], ".raw", sep = "")))) {
          write.table(round(object$margin[[j]]$X, 5), file = file.path(dir, paste(Xn[j], ".raw", sep = "")),
            quote = FALSE, row.names = FALSE)
        } else exists <- c(exists, file.path(dir, paste(Xn[j], ".raw", sep = "")))
      }
    } else {
      if(!file.exists(file.path(dir, paste(Xn, ".raw", sep = "")))) {
        write.table(round(object$X, 5), file = file.path(dir, paste(Xn, ".raw", sep = "")),
          quote = FALSE, row.names = FALSE)
      } else exists <- c(exists, file.path(dir, paste(Xn, ".raw", sep = "")))
    }
    if(!is.null(object$C)) {
      if(!file.exists(file.path(dir, paste(Cn, ".raw", sep = "")))) {
        write.table(object$C, file = file.path(dir, paste(Cn, ".raw", sep = "")),
          quote = FALSE, row.names = FALSE)
      } else exists <- c(exists, file.path(dir, paste(Cn, ".raw", sep = "")))
    }
    if(!is.null(object$state$parameters)) {
      spar <- as.matrix(matrix(get.par(object$state$parameters, "b"), ncol = 1))
      if(!file.exists(file.path(dir, paste(Pn, ".raw", sep = "")))) {
        write.table(spar, file = file.path(dir, paste(Pn, ".raw", sep = "")),
          quote = FALSE, row.names = FALSE)
      } else exists <- c(exists, file.path(dir, paste(Pn, ".raw", sep = "")))
    }
    if(!is.null(object$xt[["pmean"]])) {
      spar <- as.matrix(matrix(object$xt[["pmean"]], ncol = 1))
      if(!file.exists(file.path(dir, paste(Pm, ".raw", sep = "")))) {
        write.table(spar, file = file.path(dir, paste(Pm, ".raw", sep = "")),
          quote = FALSE, row.names = FALSE)
      } else exists <- c(exists, file.path(dir, paste(Pm, ".raw", sep = "")))
    }
    cmd <- NULL
    if(is.tx) {
      for(j in seq_along(object$S)) {
        if(!(file.path(dir, paste(Sn[j], ".raw", sep = "")) %in% exists)) {
          cmd <- c(cmd,
            paste("dataset", Sn[j]),
            paste(Sn[j], ".infile using ", file.path(dir, paste(Sn[j], ".raw", sep = "")), sep = "")
          )
        }
        if(!(file.path(dir, paste(Xn[j], ".raw", sep = "")) %in% exists)) {
          cmd <- c(cmd,
            paste("dataset", Xn[j]),
            paste(Xn[j], ".infile using ", file.path(dir, paste(Xn[j], ".raw", sep = "")), sep = "")
          )
        }
      }
    } else {
      if(!(file.path(dir, paste(Sn, ".raw", sep = "")) %in% exists)) {
        cmd <- c(cmd,
          paste("dataset", Sn),
          paste(Sn, ".infile using ", file.path(dir, paste(Sn, ".raw", sep = "")), sep = "")
        )
      }
      if(!(file.path(dir, paste(Xn, ".raw", sep = "")) %in% exists)) {
        cmd <- c(cmd,
          paste("dataset", Xn),
          paste(Xn, ".infile using ", file.path(dir, paste(Xn, ".raw", sep = "")), sep = "")
        )
      }
    }
    if(!is.null(object$C)) {
      if(!(file.path(dir, paste(Cn, ".raw", sep = "")) %in% exists)) {
        cmd <- c(cmd,
          paste("dataset", Cn),
          paste(Cn, ".infile using ", file.path(dir, paste(Cn, ".raw", sep = "")), sep = "")
        )
      }
    }
    if(!is.null(object$state$parameters)) {
      if(!(file.path(dir, paste(Pn, ".raw", sep = "")) %in% exists)) {
        cmd <- c(cmd,
          paste("dataset", Pn),
          paste(Pn, ".infile using ", file.path(dir, paste(Pn, ".raw", sep = "")), sep = "")
        )
      }
    }
    if(!is.null(object$xt[["pmean"]])) {
      if(!(file.path(dir, paste(Pm, ".raw", sep = "")) %in% exists)) {
        cmd <- c(cmd,
          paste("dataset", Pm),
          paste(Pm, ".infile using ", file.path(dir, paste(Pm, ".raw", sep = "")), sep = "")
        )
      }
    }
    return(cmd)
  }

  attr(term, "write") <- write

  term
}

sx.construct.pspline.smooth <- sx.construct.ps.smooth.spec <- sx.construct.psplinerw1.smooth.spec <-
sx.construct.psplinerw2.smooth.spec <- sx.construct.pspline.smooth.spec <-
function(object, data, mcmcreg = FALSE, ...)
{
  if(length(object$p.order) == 1L) 
    m <- rep(object$p.order, 2L)
  else 
    m <- object$p.order
  m[is.na(m)] <- 2L
  object$p.order <- m
  object$p.order[1L] <- object$p.order[1L] + 1L
  if(inherits(object, "psplinerw1.smooth.spec"))
    object$p.order[2L] <- 1L
  if(inherits(object, "psplinerw2.smooth.spec"))
    object$p.order[2L] <- 2L
  if(object$bs.dim < 0L)
    object$bs.dim <- 10L
  if(length(object$p.order) > 1L) {
    if(object$p.order[2L] > 2L) {
      warning("order of the difference penalty not supported by BayesX, set to 2!")
      object$p.order <- c(object$p.order[1L], 2L)
    }
  }
  nrknots <- object$bs.dim - object$p.order[1L] + 1L
  if(nrknots < 5L) {
    warning("number of inner knots smaller than 5 not supported by BayesX, set to 5!",
      call. = FALSE)
    nrknots <- 5L
  }
  termo <- object$term
  term <- if(mcmcreg) {
    paste(termo, "(pspline,nrknots=",
      nrknots, ",degree=", object$p.order[1L], ",difforder=", object$p.order[2L], sep = "")
  } else {
    paste(termo, "(psplinerw", object$p.order[2L], ",nrknots=",
      nrknots, ",degree=", object$p.order[1L], sep = "")
  }
  object$xt[c("knots", "nrknots", "degree", "difforder")] <- NULL
  term <- paste(do.xt(term, object, c("center", "before")), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)
  
  return(term)
}

sx.construct.tensor.smooth <- sx.construct.tensor.smooth.spec <- sx.construct.t2.smooth.spec <-
function(object, dir, prg, data, mcmcreg = FALSE, ...)
{
  by <- object$term[1L]
  term <- object$term[2L]
  object <- object$margin[[1L]]
  object$bs.dim <- as.integer(object$bs.dim^2)
  object$term <- term
  object$by <- by
  term <- sx.construct(object, dir, prg, data, mcmcreg = mcmcreg, ...)
  term <- gsub("(pspline", "(tensor", term, fixed = TRUE)
  if(!mcmcreg) {
    term <- gsub("psplinerw2", "pspline2dimrw2", term)
    term <- gsub("psplinerw1", "pspline2dimrw1", term)
  }
  return(term)
}

sx.construct.ra.smooth.spec <- sx.construct.re.smooth.spec <-
sx.construct.random.smooth.spec <- sx.construct.random.effect <- function(object, data, mcmcreg = FALSE, ...)
{
  term <- object$term
  if(is.null(object$ins))
    term <- paste(term, "(random", sep = "")
  else
    term <- paste(term, "(hrandom", sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

sx.construct.rps.smooth.spec <- function(object, data, ...)
{
  term <- paste(object$term, "(hrandom_pspline,centermethod=meansum2", sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  term <- paste(object$by, term , sep = "*")

  return(term)
}

sx.construct.kr.smooth.spec <- sx.construct.kriging.smooth.spec <- function(object, data, ...)
{
  termo <- object$term
  if(length(termo) < 2L)
    stop("kriging method needs two terms!")
  if(object$bs.dim < 0L)
    object$bs.dim <- 10L
	nrknots <- object$bs.dim
  xt <- object$xt
  if(is.null(xt$full))
    term <- paste(termo[1L], "*", termo[2L], "(kriging,nrknots=", nrknots, sep = "")
  else {
    term <- paste(termo[1L], "*", termo[2L], "(kriging,full", sep = "")    
    object$xt$full <- NULL
  }
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

construct.shrw <- function(object, data, what)
{
  term <- object$term
  term <- paste(term, "(", what, sep = "")
  term <- paste(do.xt(term, object, NULL), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  return(term)
}

sx.construct.offset.smooth.spec <- function(object, data, ...)
{
  return(construct.shrw(object, data, "offset"))
}

sx.construct.mrf.smooth <- sx.construct.mrf.smooth.spec <- sx.construct.spatial.smooth.spec <- function(object, data, ...)
{
  if(is.null(object$xt))
    stop("need to supply a map object in argument xt!")  
  map.name <- help.map.name(deparse(substitute(object, env = .GlobalEnv), 
    backtick = TRUE, width.cutoff = 500L))
  if(!is.null(object$xt$map.name))
    map.name <- object$xt$map.name
  if(!is.list(object$xt))
    object$xt <- list(object$xt)
  map.name <- rmf(gsub("\\s", "", paste(map.name, sep = "", collapse = "")))

  map <- object$xt$map
  if(is.null(map)) {
    if(!is.null(object$xt$polys))
      map <- object$xt$polys
    if(!is.null(object$xt$penalty))
      map <- object$xt$penalty
  }
  if(is.null(map))
    map <- object$xt$gra
  if(is.null(map)) {
    if(!is.list(object$xt[[1L]])) {
      if(inherits(object$xt[[1L]], "gra"))
        map <- object$xt[[1L]]
      else
        map <- object$xt
    } else map <- object$xt[[1L]]
    if(is.null(map)) {
      map <- object$xt
      if(is(map, "SpatialPolygonsDataFrame"))
        map <- BayesX::sp2bnd(map)
      if(is.null(map) || (!is.list(map) && !inherits(map, "bnd") || !inherits(map, "gra")))
        stop("need to supply a bnd or graph file object in argument xt!")
    }
  }
  if(is(map, "nb"))
    map <- BayesX::nb2gra(map)
  if(inherits(map, "SpatialPolygons"))
    map <- BayesX::sp2bnd(map)
  if(!inherits(map, "bnd") && !inherits(map, "gra")) {
    if(is.list(map))
      class(map) <- "bnd"
    else
      class(map) <- "gra"
  }

  write <- function(dir = NULL) {
    if(is.null(dir)) {
      dir.create(dir <- tempfile())
      on.exit(unlink(dir))
    } else dir <- path.expand(dir)
    if(!file.exists(dir)) dir.create(dir)

    counter <- NULL
    ok <- TRUE
    files <- list.files(dir)
    while(ok) {
      classm <- class(map)
      if(length(classm) > 1L)
        if("list" %in% classm)
          class(map) <- classm[classm != "list"]
      mapfile <- paste(map.name, counter, ".", class(map), sep = "")[1]
      if(any(grepl(mapfile, files))) {
        if(is.null(counter))
          counter <- 0L
        counter <- counter + 1L
      } else ok <- FALSE
    }
    mapfile <- file.path(dir, mapfile)

    prg <- paste("map", map.name)
    if(inherits(map, "bnd")) {
      if(!file.exists(mapfile))
        BayesX::write.bnd(map = map, file = mapfile, replace = TRUE)
      prg <- c(prg, paste(map.name, ".infile using ", mapfile, sep = ""))
    } else {
      if(!is.character(map)) {
        if(!file.exists(mapfile))
          BayesX::write.gra(map = map, file = mapfile, replace = TRUE)
        prg <- c(prg, paste(map.name, ".infile, graph using ", mapfile, sep = ""))
      } else {
        stopifnot(is.character(map))
        pos <- regexpr("\\.([[:alnum:]]+)$", map)
        fext <- ifelse(pos > -1L, substring(map, pos + 1L), "")
        if(fext == "gra")
          prg <- c(prg, paste(map.name, ".infile, graph using ", path.expand(map), sep = ""))
        else
          prg <- c(prg, paste(map.name, ".infile using ", path.expand(map), sep = ""))
      }
    }
    prg
  }

  term <- object$term
  term <- paste(term, "(spatial,map=", map.name, sep = "")
  term <- paste(do.xt(term, object, c("map", "polys", "penalty", "map.name", "nb")), ")", sep = "")
  if(object$by != "NA")
    term <- make_by(term, object, data)

  attr(term, "write") <- write
  attr(term, "map.name") <- map.name

  return(term)
}


make_by <- function(term, object, data)
{
  if(!missing(data) && !is.character(data)) {
    by <- data[[object$by]]
    if(is.factor(by) && nlevels(by) > 1) {
      nocenter <- paste(c("lasso", "nigmix", "ridge", "ra"), ".smooth.spec", sep = "")
      term <- paste(paste(rmf(object$by), rmf(levels(by)), sep = ""), "*", term, sep = "")
      if((k <- length(term)) > 1)
        for(j in 1:k) 
          if(!grepl("center", term[j]) && is.null(object$xt$center))
            if(!class(object) %in% nocenter)
              term[j] <- gsub(")", ",center)", term[j])
    } else term <- paste(rmf(object$by), "*", term, sep = "")
  } else term <- paste(rmf(object$by), "*", term, sep = "")

  return(term)
}

rmf <- function(x) 
{
  for(i in 1L:length(x)) {
    for(char in c("+", "-", "*", ":", "^", "/", " ", "(", ")", "]", "[",
      ",", ".", "<", ">", "?", "!", "'", "#", "~", "`", ";", "=", "&", "$", "@")) {
      x[i] <- gsub(char, "", x[i], fixed = TRUE)
    }
  }

  return(x)
}

help.map.name <- function(x)
{
  if(is.null(x))
    return("")
  x <- splitme(x)
  if(resplit(x[1L:2L]) == "s(") {
    x <- resplit(c("sfoofun", x[2L:length(x)]))
    x <- eval(parse(text = x), envir = parent.frame())
    if(is.null(x))
      x <- "map"
  } else  x <- "map"

  return(x)
}

sfoofun <- function(x, xt = NULL, ...)
{
  if(is.null(x) || is.null(xt))
    return(NULL)
  x <- rval <- deparse(substitute(xt), backtick = TRUE, width.cutoff = 500L)
  if(is.list(xt) && length(xt)>1)
    for(i in 1L:length(xt))
      if(inherits(xt[[i]], "bnd") || inherits(xt[[i]], "gra") || inherits(xt[[i]], "list")) {
        rval <- strsplit(x, ",", " ")[[1L]]
        if(length(rval) > 1L)
          rval <- rval[i]
      }
  rval <- splitme(rval)
  if(length(rval) > 5L)
    if(resplit(rval[1L:5L]) == "list(")
      rval <- rval[6L:length(rval)]
  if(rval[length(rval)] == ")")
    rval <- rval[1L:(length(rval) - 1)]
  if(any(grepl("=", rval)))
    rval <- rval[(grep("=", rval) + 2L):length(rval)]
  rval <- resplit(rval)
   
  return(rval)
}

splitme <- function(x) {
  return(strsplit(x, "")[[1L]])
}

resplit <- function(x) {
  if(!is.null(x))
    x <- paste(x, sep = "", collapse = "")
	
  return(x)
}


## Special tensor constructors.
te2 <- function (..., k = NA, bs = "cr", m = NA, d = NA, by = NA, fx = FALSE, 
  mp = TRUE, np = TRUE, xt = NULL, id = NULL, sp = NULL, pc = NULL)
{
  vars <- as.list(substitute(list(...)))[-1]
  dim <- length(vars)
  by.var <- deparse(substitute(by), backtick = TRUE)
  term <- deparse(vars[[1]], backtick = TRUE)
  if(dim > 1) 
    for(i in 2:dim) term[i] <- deparse(vars[[i]], backtick = TRUE)
  for(i in 1:dim) term[i] <- attr(terms(stats::reformulate(term[i])), 
    "term.labels")
  if(sum(is.na(d)) || is.null(d)) {
    n.bases <- dim
    d <- rep(1, dim)
  }
  else {
    d <- round(d)
    ok <- TRUE
    if(sum(d <= 0)) 
      ok <- FALSE
    if(sum(d) != dim) 
      ok <- FALSE
    if(ok) 
      n.bases <- length(d)
    else {
      warning("something wrong with argument d.")
      n.bases <- dim
      d <- rep(1, dim)
    }
  }
  if(sum(is.na(k)) || is.null(k)) 
    k <- 5^d
  else {
    k <- round(k)
    ok <- TRUE
    if(sum(k < 3)) {
      ok <- FALSE
      warning("one or more supplied k too small - reset to default")
    }
    if(length(k) == 1 && ok) 
      k <- rep(k, n.bases)
    else if(length(k) != n.bases) 
      ok <- FALSE
    if(!ok) 
      k <- 5^d
  }
  if(sum(is.na(fx)) || is.null(fx)) 
    fx <- rep(FALSE, n.bases)
  else if(length(fx) == 1) 
    fx <- rep(fx, n.bases)
  else if(length(fx) != n.bases) {
    warning("dimension of fx is wrong")
    fx <- rep(FALSE, n.bases)
  }
  xtra <- list()
  if(is.null(xt) || length(xt) == 1)
    for(i in 1:n.bases) xtra[[i]] <- xt
  else if(length(xt) == n.bases) 
    xtra <- xt
  else stop("xt argument is faulty.")
  if(length(bs) == 1) 
    bs <- rep(bs, n.bases)
  if(length(bs) != n.bases) {
    warning("bs wrong length and ignored.")
    bs <- rep("cr", n.bases)
  }
  bs[d > 1 & (bs == "cr" | bs == "cs" | bs == "cp")] <- "tp"
  if(!is.list(m) && length(m) == 1) 
    m <- rep(m, n.bases)
  if(length(m) != n.bases) {
    warning("m wrong length and ignored.")
    m <- rep(0, n.bases)
  }
  if(!is.list(m)) 
    m[m < 0] <- 0
  if(length(unique(term)) != dim) 
    stop("Repeated variables as arguments of a smooth are not permitted")
  j <- 1
  margin <- list()
  for(i in 1:n.bases) {
    j1 <- j + d[i] - 1
    if(is.null(xt)) 
      xt1 <- NULL
    else xt1 <- xtra[[i]]
    stxt <- "s("
    for(l in j:j1) stxt <- paste(stxt, term[l], ",", sep = "")
    stxt <- paste(stxt, "k=", deparse(k[i], backtick = TRUE), 
      ",bs=", deparse(bs[i], backtick = TRUE), ",m=", deparse(m[[i]], 
        backtick = TRUE), ",xt=xt1", ")")
    margin[[i]] <- eval(parse(text = stxt))
    j <- j1 + 1
  }
  if(mp) 
    mp <- TRUE
  else mp <- FALSE
  if(np) 
    np <- TRUE
  else np <- FALSE
  full.call <- paste("te(", term[1], sep = "")
  if(dim > 1) 
    for(i in 2:dim) full.call <- paste(full.call, ",", term[i], 
      sep = "")
  label <- paste(full.call, ")", sep = "")
  if(!is.null(id)) {
    if(length(id) > 1) {
      id <- id[1]
      warning("only first element of `id' used")
    }
    id <- as.character(id)
  }
  ret <- list(margin = margin, term = term, by = by.var, fx = fx, 
    label = label, dim = dim, mp = mp, np = np, id = id, 
    sp = sp, inter = FALSE)
  if(!is.null(pc)) {
    if(length(pc) < d) 
      stop("supply a value for each variable for a point constraint")
    if(!is.list(pc)) 
      pc <- as.list(pc)
    if(is.null(names(pc))) 
      names(pc) <- unlist(lapply(vars, all.vars))
    ret$point.con <- pc
  }
  class(ret) <- "tensor.smooth.spec"
  ret
}

tx <- function(..., bs = "ps", k = -1,
  ctr = c("center", "main", "both", "both1", "both2",
    "none", "meanf", "meanfd", "meansimple", "nullspace"),
  xt = NULL, special = TRUE)
{
  if(length(k) < 2) {
    if(k < 0)
      k <- 10
  }
  object <- te2(..., bs = bs, k = k)
  object$constraint <- match.arg(ctr)
  object$label <- gsub("te(", "tx(", object$label, fixed = TRUE)
  object$special <- special
  object$xt <- xt
  if(any(i <- sapply(object$margin, class) == "mrf.smooth.spec")) {
    xt <- c(object$xt, object$margin[[i]]$xt)
    object$xt <- object$margin[[i]]$xt <- xt
  }
  class(object) <- "tensorX.smooth.spec"
  object
}

tx2 <- function(...)
{
  object <- tx(..., special = FALSE)
  object$label <- gsub("tx(", "tx2(", object$label, fixed = TRUE)
  object
}

tx3 <- function(..., bs = "ps", k = c(10, 5), ctr = c("main", "center"), xt = NULL, special = TRUE)
{
  vars <- as.character(unlist(as.list(substitute(list(...)))[-1]))
  if(length(vars) != 3L)
    stop("3 variables are necessary for the space-time interaction term!")
  bs <- rep(bs, length.out = 2)
  k <- rep(k, length.out = 2)
  ctr <- rep(ctr, length.out = 2)
  m1 <- paste('tx(', vars[1], ',bs="', bs[1], '",k=', k[1], ',ctr="', ctr[1],'")', sep = '')
  m2 <- paste('tx(', vars[2], ',', vars[3], ',bs="', bs[2], '",k=', k[2], ',ctr="', ctr[2],'")', sep = '')
  object <- list()
  object$margin <- list(
    eval(parse(text = m1)),
    eval(parse(text = m2))
  )
  object$term <- vars
  object$by <- "NA"
  object$fx <- FALSE
  object$label <- paste("tx3(", paste(vars, collapse = ","), ")", sep = "")
  object$dim <- 3
  object$special <- special
  object$constraint <- ctr[1]
  object$xt <- xt
  class(object) <- "tensorX3.smooth.spec"
  object
}

smooth.construct.tensorX3.smooth.spec <- function(object, data, knots, ...)
{
  object$margin[[1]] <- smooth.construct.tensorX.smooth.spec(object$margin[[1]], data, knots, ...)
  object$margin[[2]] <- smooth.construct.tensorX.smooth.spec(object$margin[[2]], data, knots, ...)

  object$margin[[1]]$S <- object$margin[[1]]$margin[[1]]$S
  object$margin[[2]]$X <- tensor.prod.model.matrix(list(object$margin[[2]]$margin[[1]]$X, object$margin[[2]]$margin[[2]]$X))
  object$margin[[2]]$S <- list(
    kronecker(diag(1, ncol(object$margin[[2]]$margin[[1]]$S[[1]])), object$margin[[2]]$margin[[1]]$S[[1]]) +
    kronecker(object$margin[[2]]$margin[[2]]$S[[1]], diag(1, ncol(object$margin[[2]]$margin[[2]]$S[[1]])))
  )
  object$X <- tensor.prod.model.matrix(list(object$margin[[1]]$X, object$margin[[2]]$X))
  object$S <- tensor.prod.penalties(list(object$margin[[1]]$S[[1]], object$margin[[2]]$S[[1]]))

  p1 <- ncol(object$margin[[1]]$X)
  p2 <- ncol(object$margin[[2]]$X)

  if(object$constraint %in% c("meanf", "meanfd", "meansimple", "none", "nullspace")) {
    if(object$constraint == "none")
      object$xt$nocenter <- TRUE
    else
      object$xt$centermethod <- object$constraint
  } else {
    if(object$constraint == "main") {
      A1 <- matrix(rep(1, p1), ncol = 1)
      A2 <- matrix(rep(1, p2), ncol = 1)
      I1 <- diag(p1); I2 <- diag(p2)

      A <- cbind(kronecker(A1, I2), kronecker(I1,A2))

      i <- match.index(t(A))
      A <- A[, i$nodups, drop = FALSE]

      k <- 0
      while((qr(A)$rank < ncol(A)) & (k < 100)) {
        i <- sapply(1:ncol(A), function(d) { qr(A[, -d])$rank })
        j <- which(i == qr(A)$rank)
        if(length(j))
          A <- A[, -j[1], drop = FALSE]
        k <- k + 1
      }
      if(k == 100)
        stop("rank problems with constraint matrix!")

      object$C <- t(A)
    } else {
      object$C <- matrix(1, ncol = p1 * p2)
    }
    attr(object$C, "always.apply") <- TRUE
 }

  if(is.null(object$xt$hyperprior)) {
    if(is.null(object$xt$prior))
      object$xt$prior <- "ig"
  } else {
    object$xt$prior <- switch(object$xt$hyperprior,
       "invgamma" = "ig",
       "hnormal" = "hn",
       "scaledep" = "sd",
       "hcauchy" = "hc",
       "aunif" = "u"
    )
  }
  if(is.null(object$xt$scaletau2)) {
    if(is.null(object$xt$theta)) {
      object$xt$theta <- switch(object$xt$prior,
        "sd" = 0.00877812,
        "hc" = 0.01034553,
        "hn" = 0.1457644,
        "u" = 0.2723532
      )
    }
  } else {
    object$xt$theta <- object$xt$scaletau2
  }
  object$xt$scaletau2 <- object$xt$theta
  object$xt$hyperprior <- switch(object$xt$prior,
    "ig" = "invgamma",
    "hn" = "hnormal",
    "sd" = "scaledep",
    "hc" = "hcauchy",
    "u" = "aunif"
  )

  if(!is.null(object$xt[["ft"]])) {
    if(object$xt[["ft"]]) {
      stopifnot(requireNamespace("sdPrior"))
      if(length(object$margin) > 1) {
        nraniso <- if(is.null(object$xt$nraniso)) 11 else object$xt$nraniso
        minaniso <- if(is.null(object$xt$minaniso)) 0.05 else object$xt$minaniso
        omegaseq <- seq(from = minaniso, to = 1 - minaniso, length = nraniso)
        omegaseq[5] <- 0.4099
        omegaprob <- rep(1 / nraniso, nraniso)
        object$xt$theta <- try(hyperpar_mod2(object$X, object$margin[[1]]$S[[1]], object$margin[[2]]$S[[1]], A = object$C, c = 3,
          alpha = 0.1, omegaseq = omegaseq, omegaprob = omegaprob, R = 1000,
          type = toupper(object$xt$prior), lowrank=if(ncol(object$X) > 100) TRUE else FALSE,
          k = min(c(50, ceiling(0.5 * ncol(object$X))))), silent = TRUE)
      } else {
        object$xt$theta <- try(hyperpar_mod2(object$X, object$S[[1]], NULL, A = object$C, c = 3,
          alpha = 0.1, omegaseq = 1, omegaprob = 1, R = 1000, type = toupper(object$xt$prior),
          lowrank=if(ncol(object$X) > 100) TRUE else FALSE,
          k = min(c(50, ceiling(0.5 * ncol(object$X))))), silent = TRUE)
      }
      if(inherits(object$xt$theta, "try-error")) {
        print(object$xt$theta)
        print(object$label)
        stop("problems with sdPrior!")
      }
      object$xt$scaletau2 <- object$xt$theta
    }
  }

  object$tx.term <- paste(object$term, collapse = "")
  object$sx.S <- lapply(object$margin, function(x) { x$S[[1]] })
  object$sx.rank <- qr(do.call("+", object$S))$rank
  if(!(object$constraint %in% c("center", "meanf", "meanfd", "meansimple", "none", "nullspace")))
    object$sx.rank <- object$sx.rank - nrow(object$C)
  object$side.constrain <- if(object$special) FALSE else TRUE

  class(object) <- "tensorX3.smooth"

  object
}

Predict.matrix.tensorX3.smooth <- function(object, data) 
{
  tensor.prod.model.matrix(list(Predict.matrix(object$margin[[1]], data),
    Predict.matrix(object$margin[[2]]$margin[[1]], data),
    Predict.matrix(object$margin[[2]]$margin[[2]], data)))
}


tx4 <- function(..., ctr = c("center", "main", "both", "both1", "both2"))
{
  rval <- te(...)
  rval$special <- TRUE
  rval$mp <- TRUE
  rval$xt$doC <- TRUE
  rval$xt$constraint <- match.arg(ctr)
  rval$label <- gsub("te(", "tx4(", rval$label, fixed = TRUE)
  rval
}


smooth.construct.tensorX.smooth.spec <- function(object, data, knots, ...)
{
  if(length(object$margin) > 2)
    stop("more than two variables in tx() currently not supported!")

  side.constrain <- if(object$special) FALSE else TRUE

  object$np <- FALSE
  object$by.done <- TRUE
  if(is.null(object$inter))
    object$inter <- FALSE
  object <- smooth.construct.tensor.smooth.spec(object, data, knots)
  if(object$mp)
    object$sx.S <- lapply(object$margin, function(x) { x$S[[1]] })
  object$sx.rank <- qr(do.call("+", object$S))$rank

  if(object$by != "NA")
    object$label <- paste(object$label, object$by, sep = ":")

  object$side.constrain <- side.constrain

  ref <- sapply(object$margin, function(x) { inherits(x, "random.effect") })

  if(length(ref) < 2) {
    if(ref)
      object$constraint <- "center"
  }

  if(length(object$margin) > 1) {
    if(!object$mp)
      object$constraint <- "none"
  }

  if(object$constraint %in% c("meanf", "meanfd", "meansimple", "none", "nullspace")) {
    if(object$constraint == "none")
      object$xt$nocenter <- TRUE
    else
      object$xt$centermethod <- object$constraint
  } else {
    if(length(object$margin) < 2) {
      p <- ncol(object$margin[[1]]$X)
      object$C <- matrix(1, ncol = p)
      if(object$constraint == "main") {
        object$C <- t(cbind(1, 1:p))
      }
    } else {
      p1 <- ncol(object$margin[[2]]$X); p2 <- ncol(object$margin[[1]]$X)
      I1 <- diag(p1); I2 <- diag(p2)

      if(object$constraint == "center") {
        object$C <- matrix(1, ncol = p1 * p2)
      } else {
        if(object$constraint == "main") {
          ## Remove main effects only.
          A1 <- matrix(rep(1, p1), ncol = 1)
          A2 <- matrix(rep(1, p2), ncol = 1)
        }
        if(object$constraint == "both") {
          ## Remove main effects and varying coefficients.
          A1 <- if(ref[1]) rep(0, p1) else cbind(rep(1, p1), 1:p1)
          A2 <- if(ref[2]) rep(0, p2) else cbind(rep(1, p2), 1:p2)
        }
        if(object$constraint == "both1") {
          ## Remove main effects and varying coefficients.
          A1 <- matrix(rep(1, p1), ncol = 1)
          A2 <- cbind(rep(1, p2), 1:p2)
        }
        if(object$constraint == "both2") {
          ## Remove main effects and varying coefficients.
          A1 <- cbind(rep(1, p1), 1:p1)
          A2 <- matrix(rep(1, p2), ncol = 1)
        }

        if(ref[1])
          A1 <- matrix(rep(1, p1), ncol = 1)
        if(ref[2])
          A2 <- matrix(rep(1, p2), ncol = 1)

        A <- cbind(kronecker(A1, I2), kronecker(I1,A2))

        i <- match.index(t(A))
        A <- A[, i$nodups, drop = FALSE]

        k <- 0
        while((qr(A)$rank < ncol(A)) & (k < 100)) {
          i <- sapply(1:ncol(A), function(d) { qr(A[, -d])$rank })
          j <- which(i == qr(A)$rank)
          if(length(j))
            A <- A[, -j[1], drop = FALSE]
          k <- k + 1
        }
        if(k == 100)
          stop("rank problems with constraint matrix!")

        object$C <- t(A)
      }
    }
  }

  if(length(object$margin) > 1) {
    object$tx.term <- paste(unlist(lapply(object$margin, function(x) {
      paste(x$term, collapse = "")
    })), collapse = "")
  }

  if(object$by != "NA") {
    object$X <- data[[object$by]] * object$X
  }

  if(is.null(object$xt$hyperprior)) {
    if(is.null(object$xt$prior))
      object$xt$prior <- "ig"
  } else {
    object$xt$prior <- switch(object$xt$hyperprior,
       "invgamma" = "ig",
       "hnormal" = "hn",
       "scaledep" = "sd",
       "hcauchy" = "hc",
       "aunif" = "u"
    )
  }
  if(is.null(object$xt$scaletau2)) {
    if(is.null(object$xt$theta)) {
      object$xt$theta <- switch(object$xt$prior,
        "sd" = 0.00877812,
        "hc" = 0.01034553,
        "hn" = 0.1457644,
        "u" = 0.2723532
      )
    }
  } else {
    object$xt$theta <- object$xt$scaletau2
  }
  object$xt$scaletau2 <- object$xt$theta
  object$xt$hyperprior <- switch(object$xt$prior,
    "ig" = "invgamma",
    "hn" = "hnormal",
    "sd" = "scaledep",
    "hc" = "hcauchy",
    "u" = "aunif"
  )

  if(!is.null(object$xt[["ft"]])) {
    if(object$xt[["ft"]]) {
      stopifnot(requireNamespace("sdPrior"))
      if(length(object$margin) > 1) {
        nraniso <- if(is.null(object$xt$nraniso)) 11 else object$xt$nraniso
        minaniso <- if(is.null(object$xt$minaniso)) 0.05 else object$xt$minaniso
        omegaseq <- seq(from = minaniso, to = 1 - minaniso, length = nraniso)
        omegaseq[5] <- 0.4099
        omegaprob <- rep(1 / nraniso, nraniso)

        object$xt$theta <- try(hyperpar_mod2(unique(object$X), object$margin[[1]]$S[[1]], object$margin[[2]]$S[[1]], A = object$C, c = 3,
          alpha = 0.1, omegaseq = omegaseq, omegaprob = omegaprob, R = 1000, type = toupper(object$xt$prior),
          lowrank=if(ncol(object$X) > 100) TRUE else FALSE, k = min(c(50, ceiling(0.5 * ncol(object$X))))), silent = TRUE)
      } else {
        object$xt$theta <- try(hyperpar_mod2(unique(object$X), object$S[[1]], NULL, A = object$C, c = 3,
          alpha = 0.1, omegaseq = 1, omegaprob = 1, R = 1000, type = toupper(object$xt$prior),
          lowrank=if(ncol(object$X) > 100) TRUE else FALSE, k = min(c(50, ceiling(0.5 * ncol(object$X))))), silent = TRUE)
      }
      if(inherits(object$xt$theta, "try-error")) {
        print(object$xt$theta)
        print(object$label)
        stop("problems with sdPrior!")
      }
      object$xt$scaletau2 <- object$xt$theta
    }
  }

  if(!is.null(object$C))
    attr(object$C, "always.apply") <- TRUE
  
  class(object) <- "tensorX.smooth"

  if(length(object$margin) > 1) {
    if(!(object$constraint %in% c("center", "meanf", "meanfd", "meansimple", "none", "nullspace")))
      object$sx.rank <- object$sx.rank - nrow(object$C)
    if(!object$mp) {
      class(object) <- "tensor.smooth"
    }
  }

  return(object)
}

Predict.matrix.tensorX.smooth <- function(object, data) 
{
  Predict.matrix.tensor.smooth(object, data)
}


## Constraint matrices.
Cmat <- function(x)
{
  if(!is.null(x$margin)) {
    x$xt <- list()
    if(!is.null(x$margin[[1]]$xt$ctr))
      x$margin[[1]]$xt$constraint <- x$margin[[1]]$xt$ctr
    if(!is.null(x$margin[[1]]$xt$con))
      x$margin[[1]]$xt$constraint <- x$margin[[1]]$xt$con
    if(is.null(x$margin[[1]]$xt$constraint))
      x$margin[[1]]$xt$constraint <- "center"
    x$xt$constraint <- x$margin[[1]]$xt$constraint
  } else {
    if(!is.null(x$xt$ctr))
      x$xt$constraint <- x$xt$ctr
    if(!is.null(x$xt$con))
      x$xt$constraint <- x$xt$con
    if(is.null(x$xt$constraint))
      x$xt$constraint <- "center"
  }
  ref <- sapply(x$margin, function(z) { inherits(z, "random.effect") })
  if(length(ref)) {
    if(length(ref) < 2) {
      if(ref)
        x$xt$constraint <- "center"
    }
  }
  if(length(x$margin) < 2) {
    p <- if(is.null(x$margin)) ncol(x$X) else ncol(x$margin[[1]]$X)
    C <- matrix(1, ncol = p)
    if(x$xt$constraint == "main")
      C <- t(cbind(1, 1:p))
  } else {
    p1 <- ncol(x$margin[[2]]$X); p2 <- ncol(x$margin[[1]]$X)
    if(x$xt$constraint == "center") {
      C <- matrix(1, ncol = p1 * p2)
    } else {
      I1 <- diag(p1); I2 <- diag(p2)
      if(x$xt$constraint == "main") {
        ## Remove main effects only.
        A1 <- matrix(rep(1, p1), ncol = 1)
        A2 <- matrix(rep(1, p2), ncol = 1)
      }
      if(x$xt$constraint == "both") {
        ## Remove main effects and varying coefficients.
        A1 <- if(ref[1]) rep(0, p1) else cbind(rep(1, p1), 1:p1)
        A2 <- if(ref[2]) rep(0, p2) else cbind(rep(1, p2), 1:p2)
      }
      if(x$xt$constraint == "both1") {
        ## Remove main effects and varying coefficients.
        A1 <- matrix(rep(1, p1), ncol = 1)
        A2 <- cbind(rep(1, p2), 1:p2)
      }
      if(x$xt$constraint == "both2") {
        ## Remove main effects and varying coefficients.
        A1 <- cbind(rep(1, p1), 1:p1)
        A2 <- matrix(rep(1, p2), ncol = 1)
      }

      if(ref[1])
        A1 <- matrix(rep(1, p1), ncol = 1)
      if(ref[2])
        A2 <- matrix(rep(1, p2), ncol = 1)

      A <- cbind(kronecker(A1, I2), kronecker(I1,A2))

      i <- match.index(t(A))
      A <- A[, i$nodups, drop = FALSE]

      k <- 0
      while((qr(A)$rank < ncol(A)) & (k < 100)) {
        i <- sapply(1:ncol(A), function(d) { qr(A[, -d])$rank })
        j <- which(i == qr(A)$rank)
        if(length(j))
        A <- A[, -j[1], drop = FALSE]
        k <- k + 1
      }
      if(k == 100)
        stop("rank problems with constraint matrix!")

      C <- t(A)
    }
  }
  return(C)
}


## Download the newest version of BayesXsrc.
get_BayesXsrc <- function(dir = NULL, install = TRUE) {
  owd <- getwd()
  if(is.null(dir)) {
    dir.create(dir <- tempfile())
    on.exit(unlink(dir))
  }
  setwd(dir)
  system("svn checkout svn://scm.r-forge.r-project.org/svnroot/bayesr/pkg/BayesXsrc BayesXsrc")

  devel <- c(
    'REPOS=http://svn.gwdg.de/svn/bayesx/trunk',
    'USER=guest',
    'PASSWD=guest',
    'DIRS="adaptiv alex andrea bib dag graph leyre mcmc psplines samson structadd"',
    'FILES="export_type.h main.cpp values.h"',
    'mkdir -p src/bayesxsrc',
    'cd src/bayesxsrc',
    'for i in $DIRS ; do',
    '  svn checkout --username "${USER}" --password "${PASSWD}" $REPOS/$i $i',
    'done',
    'for i in $FILES ; do',
    '  svn export --username "${USER}" --password "${PASSWD}" $REPOS/$i $i',
    'done',
    'cd ..',
    'cp dev-Makefile Makefile',
    'cp dev-Makefile.win Makefile.win'
  )

  writeLines(devel, file.path("BayesXsrc", "bootstrap-devel.sh"))

  sh <- c(
    'cd BayesXsrc',
    'sh ./bootstrap-devel.sh',
    'cd ..',
    'R CMD build BayesXsrc',
    if(install) 'R CMD INSTALL BayesXsrc_*.tar.gz' else NULL
  )

  writeLines(sh, "get_BayesXsrc.sh")
  ok <- system('sh ./get_BayesXsrc.sh')

  setwd(owd)
  return(ok)
}


hyperpar_mod2 <- function(...)
{
  requireNamespace("sdPrior")
  if("hyperpar_mod" %in% ls(getNamespace("sdPrior"))) {
    hpf <- eval(parse(text = paste0("sdPrior", ":", ":", "hyperpar_mod")))
    return(hpf(...))
  } else stop("cannot find hyperpar_mod() in sdPrior!")
}

