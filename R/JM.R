#######################################
## Joint model model fitting engine. ##
#######################################
## (1) The family object.
jm_bamlss <- function(...)
{
  links = c(lambda = "log", gamma = "log", mu = "identity", sigma = "log",
            alpha = "identity", dalpha = "identity")
  
  rval <- list(
    "family" = "jm",
    "names" = c("lambda", "gamma", "mu", "sigma", "alpha", "dalpha"),
    "links" = links,
    "transform" = function(x, jm.start = NULL, timevar = NULL, idvar = NULL, init.models = FALSE, 
                           plot = FALSE, nonlinear = FALSE, edf_alt = FALSE, start_mu = NULL, k_mu = 6, ...) {
      rval <- jm_transform(x = x$x, y = x$y, terms = x$terms, knots = x$knots,
                formula = x$formula, family = x$family, data = x$model.frame,
                jm.start = jm.start, timevar = timevar, idvar = idvar, 
                nonlinear = nonlinear, edf_alt = edf_alt, start_mu = start_mu, k_mu = k_mu, ...)
      if(init.models) {
        x2 <- rval$x[c("mu", "sigma")]
        for(i in names(x2)) {
          for(j in seq_along(x2[[i]]$smooth.construct))
            x2[[i]]$smooth.construct[[j]]$update <- bfit_iwls
        }
        attr(x2, "bamlss.engine.setup") <- TRUE
        cat2("generating starting model for longitudinal part...\n")
        gpar <- opt_bfit(x = x2, y = rval$y[[1]][, "obs", drop = FALSE],
                     family = gF2(gaussian), start = NULL, maxit = 100)$parameters
        if(plot) {
          b <- list("x" = x2, "model.frame" = x$model.frame, "parameters" = gpar)
          class(b) <- "bamlss"
          plot(results.bamlss.default(b))
        }
        x2 <- rval$x[c("lambda", "gamma")]
        attr(x2, "bamlss.engine.setup") <- TRUE
        y <- rval$y[[1]][attr(rval$y[[1]], "take"), , drop = FALSE]
        attr(y, "width") <- attr(rval$y[[1]], "width")
        attr(y, "subdivisions") <- attr(rval$y[[1]], "subdivisions")
        cat2("generating starting model for survival part...\n")
        spar <- opt_Cox(x = x2, y = list(y), family = gF2("cox"),
                         start = NULL, maxit = 100)$parameters
        if(plot) {
          b <- list("x" = x2, "model.frame" = x$model.frame, "parameters" = spar)
          class(b) <- "bamlss"
          plot(results.bamlss.default(b))
        }
        rval$x <- set.starting.values(rval$x, start = c(unlist(gpar), unlist(spar)))
      }
      
      return(rval)
    },
    "optimizer" = opt_JM,
    "sampler" = sam_JM,
    "predict" = jm_predict
  )
  
  class(rval) <- "family.bamlss"
  rval
}


## (2) The transformer function.
jm_transform <- function(x, y, data, terms, knots, formula, family,
  subdivisions = 25, timedependent = c("lambda", "mu", "alpha", "dalpha"), 
  timevar = NULL, idvar = NULL, alpha = .Machine$double.eps, mu = NULL, sigma = NULL,
  sparse = TRUE, nonlinear = FALSE, edf_alt = FALSE, start_mu = NULL, k_mu = 6, ...)
{
  rn <- names(y)
  y0 <- y
  y <- y[[rn]]
  ntd <- timedependent
  if(!all(ntd %in% names(x)))
    stop("the time dependent predictors specified are different from family object names!")
  
  ## Create the time grid.  
  grid <- function(upper, length){
    seq(from = 0, to = upper, length = length)
  }
  take <- NULL
  if(is.character(idvar)) {
    if(!(idvar %in% names(data)))
      stop(paste("variable", idvar, "not in supplied data set!"))
    y2 <- cbind(y[, "time"], data[[idvar]])
  } else {
    did <- as.factor(idvar)
    y2 <- cbind(y[, "time"], did)
    idvar <- deparse(substitute(idvar), backtick = TRUE, width.cutoff = 500)
    idvar <- strsplit(idvar, "$", fixed = TRUE)[[1]]
    if(length(idvar) > 1)
      idvar <- idvar[2]
    data[[idvar]] <- did
  }
  colnames(y2) <- c("time", idvar)
  take <- !duplicated(y2)
  take_last <- !duplicated(y2, fromLast = TRUE)
  y2 <- y2[take, , drop = FALSE]
  nobs <- nrow(y2)
  grid <- lapply(y2[, "time"], grid, length = subdivisions)
  width <- rep(NA, nobs)
  for(i in 1:nobs)
    width[i] <- grid[[i]][2]
  y[y[, "obs"] == Inf, "obs"] <- .Machine$double.xmax * 0.9
  y[y[, "obs"] == -Inf, "obs"] <- .Machine$double.xmin * 0.9
  attr(y, "width") <- width
  attr(y, "subdivisions") <- subdivisions
  attr(y, "grid") <- grid
  attr(y, "nobs") <- nobs
  attr(y, "take") <- take
  attr(y, "take_last") <- take_last
  yname <- all.names(x[[ntd[1]]]$formula[2])[2]
  timevar_mu <- timevar
  if(is.null(timevar_mu))
    stop("the time variable is not specified, needed for mu!")
  timevar <- yname
  attr(y, "timevar") <- c("lambda" = timevar, "mu" = timevar_mu)
  attr(y, "idvar") <- idvar
  attr(y, "nonlinear") <- nonlinear
  attr(y, "edf_alt") <- edf_alt
  attr(y, "tp") <- tp <-  FALSE
  attr(y, "fac") <- fac <-  FALSE
  attr(y, "obstime") <- data[[timevar_mu]]
  attr(y, "ids") <- data[[idvar]]
  dalpha <- has_pterms(x$dalpha$terms) | (length(x$dalpha$smooth.construct) > 0)
  
  ## update alpha for nonlinear term
  if(nonlinear){
    if(is.null(start_mu)){
      yobs <- y[, "obs"]
      center <- median(yobs)
      # scale y towards center to account for outliers
      data$mu <- (0.9 * (yobs - center)) + center
    } else {
      data$mu <- start_mu
    }

    ## set up knots for spline in nonlinear
    nk <- k_mu - 2   # 2 : second order P-spline basis (cubic spline)
    bs_mu <- "ps"
    xu <- max(data$mu)
    xl <- min(data$mu)
    xr <- xu - xl
    xl <- xl - xr * 0.05
    xu <- xu + xr * 0.05
    dx <- (xu - xl)/(nk - 1)
    k <- seq(xl - dx * (2 + 1), xu + dx * (2 + 1), length = nk + 2 * 2 + 2)
    ygrid <- quantile(y[, "obs"], probs = seq(0.025, 0.975, 0.025))
    ## constraint matrix for sum-to-zero constraint on a fixed grid
    C <- matrix(apply(modSplineDesign(k, ygrid, derivs = 0), 2, mean), nrow = 1)

    nonlinear_vars <- all.vars(x$alpha$formula)[!all.vars(x$alpha$formula) %in% "alpha"]
    if(!is.null(x$alpha$smooth.construct)) {
      attr(y, "tp") <- tp <- TRUE
    } else {
      if(length(nonlinear_vars) > 0){
        if(length(nonlinear_vars) > 1) 
           stop("Only one covariate is allowed in alpha for associations nonlinear in mu.\n Please adjust your model specification.")
        if(!is.factor(data[[nonlinear_vars]])) 
          stop("Please provide the covariate in alpha as a factor or as a smooth term.")
        attr(y, "fac") <- fac <- TRUE
      }
    }

    ## updating the respective s() term
    if(!tp & !fac){
      ## smooth function of mu
      f_alpha_mu <- paste0("s(mu, bs = \"ps\", k = ", k_mu, ")")
      alpha_mu <- eval(parse(text = f_alpha_mu))
      alpha_mu$C <- C
      x$alpha$smooth.construct <- smoothCon(alpha_mu, data = data, knots = list(mu = k), 
                                            absorb.cons = TRUE, n = nrow(data), scale.penalty=FALSE) 
    }
    if(fac){
      ## smooth function of mu per factor
      f_alpha_mu <- paste0("s(mu, bs = \"ps\", by = ", nonlinear_vars, ", k = ", k_mu,  ")")
      alpha_mu <- eval(parse(text = f_alpha_mu))
      alpha_mu$C <- C
      x$alpha$smooth.construct <- smoothCon(alpha_mu, data = data, knots = list(mu = k), 
                                            absorb.cons = TRUE, n = nrow(data), scale.penalty=FALSE) 
      # insert factor again
      f_alpha_mu <- paste0(nonlinear_vars, " + s(mu, bs = \"ps\", by = ", nonlinear_vars, ", k = ", k_mu, ")")
    }
    if(tp){
      ## tensor spline over mu and covariate: under construction, not tested yet
      term <- x$alpha$smooth.construct[[1]]$term
      k_alpha <- x$alpha$smooth.construct[[1]]$bs.dim
      bs_alpha <- substr(unlist(strsplit(class(x$alpha$smooth.construct[[1]])[1], "[.]"))[1], 1, 2)
      f_alpha_mu <- paste0("te(mu, ", term, ", bs = c(\"", bs_mu, "\",\"", bs_alpha, 
                           "\"), k = c(", k_mu,",", k_alpha, "))")
      alpha_mu <- eval(parse(text = f_alpha_mu))
      alpha_mu$C <- C
      x$alpha$smooth.construct <- smoothCon(alpha_mu, data = data, knots = list(mu = k), 
                                            absorb.cons = FALSE, n = nrow(data), scale.penalty=FALSE) 
    }
    names(x$alpha$smooth.construct) <- alpha_mu$label
    formula$alpha$formula <- update(formula$alpha$formula, paste0("alpha ~ ", f_alpha_mu))
    terms$alpha <- terms(update(terms$alpha, paste0("alpha ~ ", f_alpha_mu)), specials = c("s", "te", "t2", "s2", "ti"))
  }

  ## Recompute design matrixes for lambda, gamma, alpha.
  for(j in c("lambda", "gamma", "alpha", if(dalpha) "dalpha" else NULL)) {
    if(nonlinear & (j =="alpha")){
      x[[j]]<- design.construct(terms, data = data[take_last, , drop = FALSE], knots = list(mu = k),
                 model.matrix = TRUE, smooth.construct = TRUE, model = j,
                 scale.x = FALSE, absorb.cons = TRUE, C = C)[[j]]
    } else {
      x[[j]]<- design.construct(terms, data = data[take, , drop = FALSE], knots = knots,
                 model.matrix = TRUE, smooth.construct = TRUE, model = j,
                 scale.x = FALSE)[[j]]
    }
  }
  
  ## Degrees of freedom to 1 for lambda, gamma and alpha smooths.
  for(j in c("lambda", "gamma", "alpha", if(dalpha) "dalpha" else NULL)) {
    for(sj in names(x[[j]]$smooth.construct)) {
      if(sj != "model.matrix")
        x[[j]]$smooth.construct[[sj]]$xt$df <- 3
    }
  }
  
  ## dalpha dmu initialize.
  if(dalpha) {
    terms <- c(terms, "dmu" = terms$mu)
    terms$dmu <- terms(update(terms$dmu, dmu ~ .), specials = c("s", "te", "t2", "s2", "ti"))
    formula$dmu <- formula$mu
    formula$dmu$formula <- update(formula$dmu$formula, dmu ~ .)
    i <- which(!grepl(timevar_mu, labels(terms$dmu), fixed = TRUE))
    if(length(i)) {
      terms$dmu <- drop.terms(terms$dmu, i, keep.response = TRUE)
      formula$dmu$formula <- formula(terms$dmu)
    }
    x$dmu <- design.construct(terms, data = data, knots = knots,
                              model.matrix = TRUE, smooth.construct = TRUE, model = "mu",
                              scale.x = FALSE)$mu
    if(has_sterms(terms$dmu)) {
      i <- grep(timevar_mu, names(x$dmu$smooth.construct))
      x$dmu$smooth.construct <- x$dmu$smooth.construct[i]
    }
    if(has_intercept(terms$dmu)) {
      x$dmu$model.matrix <- x$dmu$model.matrix[, -grep("Intercept", colnames(x$dmu$model.matrix)), drop = FALSE]
      x$dmu$terms <- terms(update(x$dmu$terms, . ~ -1 + .), specials = c("s", "te", "t2", "s2", "ti"))
      terms$dmu <- terms(update(terms$dmu, dmu ~ -1 + .), specials = c("s", "te", "t2", "s2", "ti"))
      formula$dmu$formula <- update(formula$dmu$formula, dmu ~ -1 + .)
      if(ncol(x$dmu$model.matrix) < 1)
        x$dmu$model.matrix <- NULL
    }
    i <- which(!grepl(timevar_mu, labels(x$dmu$terms), fixed = TRUE))
    if(length(i)) {
      x$dmu$terms <- drop.terms(x$dmu$terms, i, keep.response = TRUE)
      x$dmu$formula <- formula(x$dmu$terms)
    }
    ntd <- c(ntd, "dmu")
  }

  ## The basic setup.
  if(is.null(attr(x, "bamlss.engine.setup")))
    x <- bamlss.engine.setup(x, ...)

  ## Remove intercept from lambda and alpha.
  for (j in c("lambda", if(nonlinear) "alpha" else NULL)){
    if(!is.null(x[[j]]$smooth.construct$model.matrix)) {
      attr(terms[[j]], "intercept") <- 0
      cn <- colnames(x[[j]]$smooth.construct$model.matrix$X)
    if("(Intercept)" %in% cn)
        x[[j]]$smooth.construct$model.matrix$X <- x[[j]]$smooth.construct$model.matrix$X[, cn != "(Intercept)", drop = FALSE]
      if(ncol(x[[j]]$smooth.construct$model.matrix$X) < 1) {
        x[[j]]$smooth.construct$model.matrix <- NULL
        x[[j]]$terms <- drop.terms.bamlss(x[[j]]$terms, pterms = FALSE, keep.intercept = FALSE)
    } else {
        x[[j]]$smooth.construct$model.matrix$term <- gsub("(Intercept)+", "",
                                                          x[[j]]$smooth.construct$model.matrix$term, fixed = TRUE)
        x[[j]]$smooth.construct$model.matrix$state$parameters <- x[[j]]$smooth.construct$model.matrix$state$parameters[-1]
        x[[j]]$terms <- drop.terms.bamlss(x[[j]]$terms, keep.intercept = FALSE) # instead of attr(x[[j]]$terms, "intercept") <- 0 
        # additional changes to remove intercept
        x[[j]]$smooth.construct$model.matrix$label <- gsub("(Intercept)+", "",
                                                           x[[j]]$smooth.construct$model.matrix$label, fixed = TRUE)
        x[[j]]$smooth.construct$model.matrix$bs.dim <- as.integer(x[[j]]$smooth.construct$model.matrix$bs.dim - 1)

        pid <- !grepl("tau", names(x[[j]]$smooth.construct$model.matrix$state$parameters)) &
          !grepl("edf", names(x[[j]]$smooth.construct$model.matrix$state$parameters))
        x[[j]]$smooth.construct$model.matrix$pid <- list("b" = which(pid), "tau2" = which(!pid))
        if(!length(x[[j]]$smooth.construct$model.matrix$pid$tau2))
          x[[j]]$smooth.construct$model.matrix$pid$tau2 <- NULL
        x[[j]]$smooth.construct$model.matrix$sparse.setup$matrix <- NULL
      }
    }
  }
  
  ## Assign time grid predict functions.
  for(i in seq_along(ntd)) {
    if(has_pterms(x[[ntd[i]]]$terms)) {
      if(ntd[i] %in% c("lambda", if(nonlinear) "alpha")) {
        x[[ntd[i]]]$smooth.construct$model.matrix <- param_time_transform2(x[[ntd[i]]]$smooth.construct$model.matrix,
                                                                           drop.terms.bamlss(x[[ntd[i]]]$terms, sterms = FALSE, keep.response = FALSE), data, grid, yname,
                                                                           timevar_mu, take, timevar2 = timevar_mu, idvar = idvar, delete.intercept = TRUE)
      } else {
        x[[ntd[i]]]$smooth.construct$model.matrix <- param_time_transform(x[[ntd[i]]]$smooth.construct$model.matrix,
                                                                          drop.terms.bamlss(x[[ntd[i]]]$terms, sterms = FALSE, keep.response = FALSE), data, grid, yname, 
                                                                          if(ntd[i] != "mu" & ntd[i] != "dmu") timevar else timevar_mu, take, derivMat = (ntd[i] == "dmu"))
      }
    }
    if(length(x[[ntd[i]]]$smooth.construct)) {
      for(j in names(x[[ntd[i]]]$smooth.construct)) {
        if(j != "model.matrix") {
          xterm <- x[[ntd[i]]]$smooth.construct[[j]]$term
          by <- if(x[[ntd[i]]]$smooth.construct[[j]]$by != "NA") x[[ntd[i]]]$smooth.construct[[j]]$by else NULL
          x[[ntd[i]]]$smooth.construct[[j]] <- sm_time_transform2(x[[ntd[i]]]$smooth.construct[[j]],
                                                                  data[, unique(c(xterm, yname, by, timevar, timevar_mu, idvar)), drop = FALSE], grid, yname,
                                                                  if(ntd[i] != "mu" & ntd[i] != "dmu") timevar else timevar_mu, take_last, derivMat = (ntd[i] == "dmu"))
        }
      }
    }
  }
  ## Assign time grid predict function to margin of alpha
  if(tp){
    xterm <- x[["alpha"]]$smooth.construct[[1]]$margin[[1]]$term
    by <- NULL
    x[["alpha"]]$smooth.construct[[1]]$margin[[1]] <- sm_time_transform2(x[["alpha"]]$smooth.construct[[1]]$margin[[1]],
                                                                         data[, unique(c(xterm, yname, by, timevar, timevar_mu, idvar)), drop = FALSE], grid, yname,
                                                                         timevar, take_last, derivMat = FALSE)
    
    xterm <- x[["alpha"]]$smooth.construct[[1]]$margin[[2]]$term
    by <- NULL
    x[["alpha"]]$smooth.construct[[1]]$margin[[2]] <- sm_time_transform2(x[["alpha"]]$smooth.construct[[1]]$margin[[2]],
                                                                         data[, unique(c(xterm, yname, by, timevar, timevar_mu, idvar)), drop = FALSE], grid, yname,
                                                                         timevar, take_last, derivMat = FALSE)
  }


  ## Set alpha/mu/sigma intercept starting value.
  if(!is.null(x$alpha$smooth.construct$model.matrix)) {
    if(alpha == 0)
      alpha <- .Machine$double.eps
    x$alpha$smooth.construct$model.matrix$state$parameters[1] <- alpha
    x$alpha$smooth.construct$model.matrix$state$fitted.values <- x$alpha$smooth.construct$model.matrix$X %*% x$alpha$smooth.construct$model.matrix$state$parameters
  }
  if(!is.null(x$mu$smooth.construct$model.matrix)) {
    if(!is.null(mu)) {
      if(mu == 0)
        mu <- .Machine$double.eps
    }
    x$mu$smooth.construct$model.matrix$state$parameters[1] <- if(is.null(mu)) mean(y[, "obs"], na.rm = TRUE) else mu
    x$mu$smooth.construct$model.matrix$state$fitted.values <- x$mu$smooth.construct$model.matrix$X %*% x$mu$smooth.construct$model.matrix$state$parameters
  }
  if(!is.null(x$sigma$smooth.construct$model.matrix)) {
    x$sigma$smooth.construct$model.matrix$state$parameters[1] <- if(is.null(sigma)) log(sd(y[, "obs"], na.rm = TRUE)) else sigma
    x$sigma$smooth.construct$model.matrix$state$fitted.values <- x$sigma$smooth.construct$model.matrix$X %*% x$sigma$smooth.construct$model.matrix$state$parameters
  }
  
  ## Make fixed effects term the last term in mu.
  if(!is.null(x$mu$smooth.construct$model.matrix) & (length(x$mu$smooth.construct) > 1)) {
    nmu <- names(x$mu$smooth.construct)
    nmu <- c(nmu[nmu != "model.matrix"], "model.matrix")
    x$mu$smooth.construct <- x$mu$smooth.construct[nmu]
    if(dalpha) {
      if(!is.null(x$dmu$smooth.construct$model.matrix) & (length(x$dmu$smooth.construct) > 1)) {
        nmu <- names(x$dmu$smooth.construct)
        nmu <- c(nmu[nmu != "model.matrix"], "model.matrix")
        x$dmu$smooth.construct <- x$dmu$smooth.construct[nmu]
      }
    }
  }
  
  ## Sparse matrix setup for mu/dmu.
  for(j in c("mu", if(dalpha) "dmu" else NULL)) {
    for(sj in seq_along(x[[j]]$smooth.construct)) {
      x[[j]]$smooth.construct[[sj]] <- sparse_Matrix_setup(x[[j]]$smooth.construct[[sj]], sparse = sparse, take = take, nonlinear = nonlinear)
    }
  }
  
  ## Update prior/grad/hess functions
  for(j in names(x)) {
    for(sj in names(x[[j]]$smooth.construct)) {
      priors <- make.prior(x[[j]]$smooth.construct[[sj]])
      x[[j]]$smooth.construct[[sj]]$prior <- priors$prior
      x[[j]]$smooth.construct[[sj]]$grad <- priors$grad
      x[[j]]$smooth.construct[[sj]]$hess <- priors$hess
    }
  }
  
  ## Add small constant to penalty diagonal for stabilization
  for(j in names(x)) {
    for(sj in names(x[[j]]$smooth.construct)) {
      if(length(x[[j]]$smooth.construct[[sj]]$S)){
        for(k in seq_along(x[[j]]$smooth.construct[[sj]]$S)){
          if(!is.list(x[[j]]$smooth.construct[[sj]]$S[[k]]) & !is.function(x[[j]]$smooth.construct[[sj]]$S[[k]])) {
            nc <- ncol(x[[j]]$smooth.construct[[sj]]$S[[k]])
            x[[j]]$smooth.construct[[sj]]$S[[k]] <- x[[j]]$smooth.construct[[sj]]$S[[k]] + diag(1e-08, nc, nc)
          }
        }
      }
    }
  }
  

  y0[[rn]] <- y
  
  family$p2logLik <- function(par, logPost = FALSE, ...) {
    eta <- eta_timegrid <- list()
    lprior <- 0
    nx <- c("lambda", "mu", "gamma", "sigma", "alpha")
    if(dalpha)
      nx <- c(nx, "dalpha", "dmu")
    for(j in nx) {
      eta[[j]] <- 0
      if(!(j %in% c("gamma", "sigma")))
        eta_timegrid[[j]] <- 0
      for(sj in names(x[[j]]$smooth.construct)) {
        pn <- paste(j, if(sj != "model.matrix") "s" else "p", sep = ".")
        cn <- colnames(x[[j]]$smooth.construct[[sj]]$X)
        if(is.null(cn))
          cn <- paste("b", 1:ncol(x[[j]]$smooth.construct[[sj]]$X), sep = "")
        pn0 <- paste(pn, sj, sep = ".")
        pn <- paste(pn0, cn, sep = ".")
        if(all(is.na(par[pn]))) {
          if(sj == "model.matrix")
            pn <- gsub(".p.model.matrix.", ".p.", pn, fixed = TRUE)
        }
        eta[[j]] <- eta[[j]] + x[[j]]$smooth.construct[[sj]]$fit.fun(x[[j]]$smooth.construct[[sj]]$X, par[pn])
        if(!(j %in% c("gamma", "sigma")))
          eta_timegrid[[j]] <- eta_timegrid[[j]] + x[[j]]$smooth.construct[[sj]]$fit.fun_timegrid(par[pn])
        if(logPost) {
          pn2 <- paste(pn0, "tau2", sep = ".")
          tpar <- par[grep2(c(pn, pn2), names(par), fixed = TRUE)]
          lprior <- lprior + x[[j]]$smooth.construct[[sj]]$prior(tpar)
        }
      }
    }
    eta_timegrid <- if(dalpha) {
      eta_timegrid$lambda + eta_timegrid$alpha * eta_timegrid$mu + eta_timegrid$dalpha * eta_timegrid$dmu
    } else {
      if(nonlinear){
        eta_timegrid$lambda + eta_timegrid$alpha
      } else {
      eta_timegrid$lambda + eta_timegrid$alpha * eta_timegrid$mu
    }
    }
    eeta <- exp(eta_timegrid)
    int <- width * (0.5 * (eeta[, 1] + eeta[, subdivisions]) + apply(eeta[, 2:(subdivisions - 1)], 1, sum))
    logLik <- sum((eta_timegrid[, ncol(eta_timegrid)] + eta$gamma) * y[attr(y, "take"), "status"], na.rm = TRUE) -
      exp(eta$gamma) %*% int + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
    if(logPost)
      logLik <- logLik + lprior
    return(drop(logLik))
  }
  
  return(list("x" = x, "y" = y0, "family" = family, "terms" = terms, "formula" = formula))
}


sparse_Matrix_setup <- function(x, sparse = TRUE, force = FALSE, take, nonlinear)
{
  if(sparse & ((ncol(x$sparse.setup$crossprod) < (ncol(x$X) * 0.5)) | force)) {
      x$X <- Matrix(x$X, sparse = TRUE)
      x$XT <- Matrix(x$XT, sparse = TRUE)
      for(j in seq_along(x$S))
        x$S[[j]] <- Matrix(x$S[[j]], sparse = TRUE)
    if(nonlinear){
      x$update <- update_jm_mu_nonlin_Matrix
      x$propose <- propose_jm_mu_nonlin_Matrix
    } else{
      x$update <- update_jm_mu_Matrix
      x$propose <- propose_jm_mu_Matrix
    }
  } else {
    if(nonlinear){
      x$update <- update_jm_mu_nonlin
      x$propose <- propose_jm_mu_nonlin
    } else{
    x$update <- update_jm_mu
    x$propose <- propose_jm_mu_simple
  }
  }
  if(!is.null(x$sparse.setup$matrix)) {
    x$sparse.setup[["mu.matrix"]] <- x$sparse.setup$matrix[take, , drop = FALSE]
  }
  return(x)
}


opt_JM <- jm_mode <- function(x, y, start = NULL, weights = NULL, offset = NULL,
  criterion = c("AICc", "BIC", "AIC"), maxit = c(100, 1),
  nu = c("lambda" = 0.1, "gamma" = 0.1, "mu" = 0.1, "sigma" = 0.1, "alpha" = 0.1, "dalpha" = 0.1),
  update.nu = FALSE, eps = 0.0001, alpha.eps = 0.001, ic.eps = 1e-08, nback = 40,
  verbose = TRUE, digits = 4, ...)
{
  ## Hard coded.
  fix.lambda <- FALSE
  fix.gamma <- FALSE
  fix.alpha <- FALSE
  fix.mu <- FALSE
  fix.sigma <- FALSE
  fix.dalpha <- FALSE
  nu.eps <- -Inf
  edf.eps <- Inf
  
  dalpha <- has_pterms(x$dalpha$terms) | (length(x$dalpha$smooth.construct) > 0)
  
  if(!is.null(start)) {
    if(is.matrix(start)) {
      if(any(i <- grepl("Mean", colnames(start))))
        start <- start[, i]
      else stop("the starting values should be a vector not a matrix!")
    }
    if(dalpha) {
      if(!any(grepl("dmu.", names(start), fixed = TRUE))) {
        start.dmu <- start[grep("mu.", names(start), fixed = TRUE)]
        names(start.dmu) <- gsub("mu.", "dmu.", names(start.dmu), fixed = TRUE)
        start <- c(start, start.dmu)
      }
    }
    x <- set.starting.values(x, start)
  }
  
  criterion <- match.arg(criterion)
  ia <- interactive()
  
  ## Names of parameters/predictors.
  nx <- names(x)
  
  ## Extract y.
  y <- y[[1]]
  
  ## Number of observations.
  nobs <- attr(y, "nobs")
  
  ## Number of subdivions used for the time grid.
  sub <- attr(y, "subdivisions")
  
  ## The interval width from subdivisons.
  width <- attr(y, "width")
  
  ## Subject specific indicator
  take <- attr(y, "take")
  take_last <- attr(y, "take_last")
  nlong <- length(take)
  
  ## Extract the status for individual i.
  status <- y[take, "status"]
  
  ## nonlinear setup
  nonlinear <- attr(y, "nonlinear")
  tp <- attr(y, "tp") 
  fac <- attr(y, "fac")
  edf_alt <- attr(y, "edf_alt")
  
  ## Make id for individual i.
  id <- which(take)
  id <- append(id[-1], nlong + 1) - id 
  id <- rep(1:nobs, id)
  
  ## Compute additive predictors.
  eta <- get.eta(x, expand = FALSE)
  
  ## Correct gamma predictor if NULL.
  if(!length(x$gamma$smooth.construct)) {
    eta$gamma <- rep(0, length = nobs)
  }
  
  nu <- rep(nu, length.out = 6)
  if(is.null(names(nu)))
    names(nu) <- c("lambda", "gamma", "mu", "sigma", "alpha", "dalpha")
  
  ## For the time dependent part, compute
  ## predictors based on the time grid.
  eta_timegrid_mu <- 0
  if(length(x$mu$smooth.construct)) {
    for(j in names(x$mu$smooth.construct)) {
      b <- get.par(x$mu$smooth.construct[[j]]$state$parameters, "b")
      x$mu$smooth.construct[[j]]$state$fitted_timegrid <- x$mu$smooth.construct[[j]]$fit.fun_timegrid(b)
      eta_timegrid_mu <- eta_timegrid_mu + x$mu$smooth.construct[[j]]$fit.fun_timegrid(b)
      x$mu$smooth.construct[[j]]$state$nu <- nu["mu"]
    }
  }
  
  eta_timegrid_lambda <- 0
  if(length(x$lambda$smooth.construct)) {
    for(j in names(x$lambda$smooth.construct)) {
      b <- get.par(x$lambda$smooth.construct[[j]]$state$parameters, "b")
      eta_timegrid_lambda <- eta_timegrid_lambda + x$lambda$smooth.construct[[j]]$fit.fun_timegrid(b)
      x$lambda$smooth.construct[[j]]$state$nu <- nu["lambda"]
    }
  }
  
  eta_timegrid_dmu <- eta_timegrid_dalpha <- 0
  if(dalpha) {
    if(length(x$dalpha$smooth.construct)) {
      for(j in names(x$dalpha$smooth.construct)) {
        b <- get.par(x$dalpha$smooth.construct[[j]]$state$parameters, "b")
        eta_timegrid_dalpha <- eta_timegrid_dalpha + x$dalpha$smooth.construct[[j]]$fit.fun_timegrid(b)
        x$dalpha$smooth.construct[[j]]$state$nu <- nu["dalpha"]
      }
    }
    if(length(x$dmu$smooth.construct)) {
      for(j in names(x$dmu$smooth.construct)) {
        b <- get.par(x$mu$smooth.construct[[j]]$state$parameters, "b")
        eta_timegrid_dmu <- eta_timegrid_dmu + x$dmu$smooth.construct[[j]]$fit.fun_timegrid(b)
        eta$dmu <- x$dmu$smooth.construct[[j]]$fit.fun(x$dmu$smooth.construct[[j]]$X, b)
      }
    }
  }
  
  eta_timegrid_alpha <- 0
  if(length(x$alpha$smooth.construct)) {
    for(j in names(x$alpha$smooth.construct)) {
      b <- get.par(x$alpha$smooth.construct[[j]]$state$parameters, "b")
      # update alpha according to mu
      if(nonlinear & !("model.matrix" %in% attr(x$alpha$smooth.construct[[j]], "class"))){
        Xmu <- as.vector(t(eta_timegrid_mu))
        if(!tp) {
          X <- x$alpha$smooth.construct[[j]]$update(Xmu, "mu")    
        } else {
          Xalpha <- x$alpha$smooth.construct[[j]]$margin[[1]]$update(Xmu, "mu")
          Xalpha2 <- x$alpha$smooth.construct[[j]]$margin[[2]]$fit.fun_timegrid(NULL)
          X <- rowTensorProduct(Xalpha, Xalpha2) 
        }
        fit_timegrid <- matrix(drop(X %*% b), nrow = nrow(eta_timegrid_mu), 
                               ncol = ncol(eta_timegrid_mu), byrow = TRUE)
        x$alpha$smooth.construct[[j]]$state$fitted.values <- fit_timegrid[, ncol(fit_timegrid)]
        x$alpha$smooth.construct[[j]]$state$fitted_timegrid <- fit_timegrid
        eta_timegrid_alpha <- eta_timegrid_alpha + fit_timegrid
      } else {
        x$alpha$smooth.construct[[j]]$state$fitted_timegrid <- x$alpha$smooth.construct[[j]]$fit.fun_timegrid(b)
        eta_timegrid_alpha <- eta_timegrid_alpha + x$alpha$smooth.construct[[j]]$fit.fun_timegrid(b)   
      }
      x$alpha$smooth.construct[[j]]$state$nu <- nu["alpha"]
    }
    eta$alpha <- eta_timegrid_alpha[, ncol(eta_timegrid_alpha)]
  } 

  # for nonlinear effect eta_mu is within eta_alpha
  if(nonlinear){
    eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
  } else{
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
  }
  knots <- if(nonlinear){
              if(tp){x$alpha$smooth.construct[[1]]$margin[[1]]$knots} 
              else {x$alpha$smooth.construct[[1]]$knots}} 
           else {NULL}
  
  for(k in c("gamma", "sigma")) {
    if(length(x[[k]]$smooth.construct)) {
      for(j in names(x[[k]]$smooth.construct)) {
        x[[k]]$smooth.construct[[j]]$state$nu <- nu[k]
      }
    }
  }
  
  ## Set edf to zero 4all terms.
  for(k in names(x)) {
    if(length(x[[k]]$smooth.construct)) {
      for(j in names(x[[k]]$smooth.construct)) {
        x[[k]]$smooth.construct[[j]]$state$edf <- 0
      }
    }
  }
  
  ## Extract current value of the log-posterior.
  get_LogPost <- function(eta_timegrid, eta, log.priors) {
    eeta <- exp(eta_timegrid)
    int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
    logLik <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
      drop(exp(eta$gamma)) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
    logPost <- as.numeric(logLik + log.priors)
    return(logPost)
  }
  
  eps <- rep(eps, length.out = 2)
  if(length(maxit) < 2)
    maxit <- c(maxit, 1)
  
  ## Start the backfitting algorithm.
  eps0 <- eps0_surv <- eps0_long <- eps0_alpha <- eps0_dalpha <- eps[1] + 1
  iter <- 0; ic_contrib <- NULL
  edf <- 0; ok.alpha <- FALSE
  ptm <- proc.time()
  while((eps0 > eps[1]) & (iter < maxit[1])) {
    if(eps0 < nu.eps)
      update.nu <- FALSE
    
    eta0_surv <- do.call("cbind", eta[c("lambda", "gamma")])
    eta0_long <- rowsum(do.call("cbind", eta[c("mu", "sigma")]), id)
    eta0_alpha <- if(is.null(dalpha)) eta$alpha else do.call("cbind", eta[c("alpha", "dalpha")])
    
    if(max(c(eps0_surv, eps0_long)) < alpha.eps) 
      ok.alpha <- TRUE
    #checks# Setup for checks while running
    myplots <- FALSE
    if(myplots){
      times <- matrix(unlist(attr(y, "grid")), ncol = ncol(eta_timegrid), byrow = TRUE)
      obstimes <- attr(y, "obstime")
      ids <- attr(y, "ids")
      par(mfrow = c(2,2))
    }
    mydiag <- FALSE
    
    ################
    ## Alpha part ##
    ################
    if(!fix.alpha & ok.alpha) {
      eps00 <- eps[2] + 1; iter00 <- 0
      eta00 <- eta$alpha
      while(eps00 > eps[2] & iter00 < maxit[2]) {
        for(sj in seq_along(x$alpha$smooth.construct)) {
          # #checks#
          # cat("\n iteration: ", iter, ", predictor: alpha", ", term: ", sj)
          if(nonlinear){
            state <- update_jm_alpha_nonlin(x$alpha$smooth.construct[[sj]], eta, eta_timegrid,
                                           eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                                           status, update.nu, width, criterion, get_LogPost, nobs, eps0_alpha < edf.eps, edf = edf, edf_alt, 
                                           tp = tp, ...)
          } else {
          state <- update_jm_alpha(x$alpha$smooth.construct[[sj]], eta, eta_timegrid,
                                   eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                                   status, update.nu, width, criterion, get_LogPost, nobs, eps0_alpha < edf.eps, edf = edf, ...)
          }
          eta_timegrid_alpha <- eta_timegrid_alpha - x$alpha$smooth.construct[[sj]]$state$fitted_timegrid + state$fitted_timegrid
          if(nonlinear){
            eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
          } else {
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
          }       
          eta$alpha <- eta$alpha - fitted(x$alpha$smooth.construct[[sj]]$state) + fitted(state)
          #checks#
          if(myplots){
            plot(eta_timegrid_mu, state$fitted_timegrid, 
                 main=paste0("alpha: iteration ", iter, ", edf ", round(state$edf, 2), ", term: ", names(x$alpha$smooth.construct)[sj]))
            abline(h = 0, col = 2)
            abline(v = median(y[, 3]), col = 2)
            plot(eta_timegrid_mu, eta_timegrid_alpha, main=paste0("alpha: iteration ", iter, ", edf ", round(state$edf, 2)))
            plot(eta_timegrid_mu[, ncol(eta_timegrid_mu)], eta$alpha)
            matplot(t(times), t(eta_timegrid_alpha), main=paste0("alpha: iteration ", iter, ", edf ", round(state$edf, 2)), type = "l")
          }
          if(mydiag){
            if(min(eta_timegrid_alpha) < -10 | max(eta_timegrid_alpha) > 10) browser()
          }
          edf <- edf - x$alpha$smooth.construct[[sj]]$state$edf + state$edf
          x$alpha$smooth.construct[[sj]]$state <- state
        }
        
        if(maxit[2] > 1) {
          eps00 <- mean(abs((eta$alpha - eta00) / eta$alpha), na.rm = TRUE)
          if(is.na(eps00) | !is.finite(eps00)) eps00 <- eps[2] + 1
        }
        
        iter00 <- iter00 + 1
      }
    }
    
    if(!fix.dalpha & dalpha & ok.alpha) {
      eps00 <- eps[2] + 1; iter00 <- 0
      eta00 <- eta$dalpha
      while(eps00 > eps[2] & iter00 < maxit[2]) {
        for(sj in seq_along(x$dalpha$smooth.construct)) {
          state <- update_jm_dalpha(x$dalpha$smooth.construct[[sj]], eta, eta_timegrid,
                                    eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                                    status, update.nu, width, criterion, get_LogPost, nobs, eps0_alpha < edf.eps, edf = edf, ...)
          eta_timegrid_dalpha <- eta_timegrid_dalpha - x$dalpha$smooth.construct[[sj]]$state$fitted_timegrid + state$fitted_timegrid
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
          eta$dalpha <- eta$dalpha - fitted(x$dalpha$smooth.construct[[sj]]$state) + fitted(state)
          edf <- edf - x$dalpha$smooth.construct[[sj]]$state$edf + state$edf
          x$dalpha$smooth.construct[[sj]]$state <- state
        }
        
        if(maxit[2] > 1) {
          eps00 <- mean(abs((eta$dalpha - eta00) / eta$dalpha), na.rm = TRUE)
          if(is.na(eps00) | !is.finite(eps00)) eps00 <- eps[2] + 1
        }
        
        iter00 <- iter00 + 1
      }
    }
    
    #######################
    ## Longitudinal part ##
    #######################
    eps00 <- eps[2] + 1; iter00 <- 0
    eta00 <- do.call("cbind", eta[c("mu", "sigma")])
    while(eps00 > eps[2] & iter00 < maxit[2]) {
      if(!fix.mu) {
        for(sj in names(x$mu$smooth.construct)) {
          # #checks#
          # cat("\n iteration: ", iter, ", predictor: mu", ", term: ", sj)
          state <- x$mu$smooth.construct[[sj]]$update(x$mu$smooth.construct[[sj]], y, eta, eta_timegrid,
                                                      eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                                                      status, update.nu, width, criterion, get_LogPost, nobs, eps0_long < edf.eps, edf = edf,
                                                      dx = if(dalpha) x$dmu$smooth.construct[[sj]] else NULL, 
                                                      xsmalpha = if(nonlinear) x$alpha$smooth.construct else NULL, 
                                                      knots = knots, tp = if(nonlinear) tp else FALSE, fac = if(nonlinear) fac else FALSE, ...)
          eta_timegrid_mu <- eta_timegrid_mu - x$mu$smooth.construct[[sj]]$state$fitted_timegrid + state$fitted_timegrid

          if(nonlinear){
            for (i in names(x$alpha$smooth.construct)){
              if(i != "model.matrix"){   # only smooth.constructs need to be updated
                g_a <- get.par(x$alpha$smooth.construct[[i]]$state$parameters, "b")
                Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
                Xalpha <- constrain(x$alpha$smooth.construct[[i]], Xalpha)
                if(fac)
                  Xalpha <- Xalpha * x$alpha$smooth.construct[[i]]$by_timegrid
                if(tp)  
                  Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
                
                alpha_state <- matrix(Xalpha %*% g_a, nrow = nrow(eta_timegrid), ncol = ncol(eta_timegrid), byrow = TRUE)
                eta_timegrid_alpha <- eta_timegrid_alpha - x$alpha$smooth.construct[[i]]$state$fitted_timegrid + alpha_state
                x$alpha$smooth.construct[[i]]$state$fitted_timegrid <- alpha_state
                eta$alpha <- eta$alpha - fitted(x$alpha$smooth.construct[[i]]$state) + alpha_state[, ncol(eta_timegrid)]
                x$alpha$smooth.construct[[i]]$state$fitted.values <- alpha_state[, ncol(eta_timegrid)]
              } 
            }

            eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
          } else {
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
          }
          eta$mu <- eta$mu - fitted(x$mu$smooth.construct[[sj]]$state) + fitted(state)
          edf <- edf - x$mu$smooth.construct[[sj]]$state$edf + state$edf
          x$mu$smooth.construct[[sj]]$state <- state
          #checks# 
          if(myplots){
            matplot(t(times),t(state$fitted_timegrid), type = "l", main = paste0("mu: iteration ", iter, ", edf ", round(state$edf, 2), " ", sj))
            matplot(t(times), t(eta_timegrid_mu), type = "l", main = paste0("mu: iteration ", iter, ", edf ", round(state$edf, 2), " ", sj))
            plot(eta_timegrid_mu, eta_timegrid_alpha, main=paste0("mu: iteration ", iter, ", edf ", round(state$edf, 2), "alpha"))
            plot(eta_timegrid_mu[, ncol(eta_timegrid_mu)], eta$alpha, main=paste0("mu: iteration ", iter, ", edf ", round(state$edf, 2), "alpha"))
            # print(xyplot(eta$mu ~ obstimes, group = ids, type = "l", main = paste0("mu: iteration ", iter, ", edf ", round(state$edf, 2), " ", sj)))
          }
          
          if(mydiag){
            if(min(eta_timegrid_alpha) < -10 | max(eta_timegrid_alpha) > 10) browser()
            if(min(eta_timegrid_mu) < -10 | max(eta_timegrid_alpha) > 20) browser()
          }
          
          #!# Danger: dalpha and nonlinear not properly set up
          if(dalpha & (sj %in% names(x$dmu$smooth.construct))) {
            state <- update_jm_dmu(x$dmu$smooth.construct[[sj]], x$mu$smooth.construct[[sj]])
            eta_timegrid_dmu <- eta_timegrid_dmu - x$dmu$smooth.construct[[sj]]$state$fitted_timegrid + state$fitted_timegrid
            eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
            eta$dmu <- eta$dmu - fitted(x$dmu$smooth.construct[[sj]]$state) + fitted(state)
            x$dmu$smooth.construct[[sj]]$state <- state
          }
        }
      }
      
      if(!fix.sigma) {
        for(sj in seq_along(x$sigma$smooth.construct)) {
          state <- update_jm_sigma(x$sigma$smooth.construct[[sj]], y, eta,
                                   eta_timegrid, update.nu, criterion, get_LogPost, nobs, eps0_long < edf.eps, edf = edf, ...)
          eta$sigma <- eta$sigma - fitted(x$sigma$smooth.construct[[sj]]$state) + fitted(state)
          edf <- edf - x$sigma$smooth.construct[[sj]]$state$edf + state$edf
          x$sigma$smooth.construct[[sj]]$state <- state
        }
      }
      
      if(maxit[2] > 1) {
        eps00 <- do.call("cbind", eta[c("mu", "sigma")])
        eps00 <- mean(abs((eps00 - eta00) / eps00), na.rm = TRUE)
        if(is.na(eps00) | !is.finite(eps00)) eps00 <- eps[2] + 1
      }
      
      iter00 <- iter00 + 1
    }
    
    eps00 <- eps[2] + 1; iter00 <- 0
    eta00 <- do.call("cbind", eta[c("lambda", "gamma")])
    while(eps00 > eps[2] & iter00 < maxit[2]) {
      ###################
      ## Survival part ##
      ###################
      eeta <- exp(eta_timegrid)
      int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
      
      if(!fix.lambda) {
        for(sj in seq_along(x$lambda$smooth.construct)) {
          # cat("\n iteration: ", iter, ", predictor: lambda", ", term: ", sj)
          state <- update_jm_lambda(x$lambda$smooth.construct[[sj]], eta, eta_timegrid,
                                    eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                                    status, update.nu, width, criterion, get_LogPost, nobs, eps0_surv < edf.eps, edf = edf, nonlinear, ...)
          eta_timegrid_lambda <- eta_timegrid_lambda - x$lambda$smooth.construct[[sj]]$state$fitted_timegrid + state$fitted_timegrid
          if(nonlinear){
            eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
          } else {
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
          }
          eta$lambda <- eta$lambda - fitted(x$lambda$smooth.construct[[sj]]$state) + fitted(state)
          edf <- edf - x$lambda$smooth.construct[[sj]]$state$edf + state$edf
          x$lambda$smooth.construct[[sj]]$state <- state
        }
      }
      if(mydiag){
        if(min(eta_timegrid) < -100 | max(eta_timegrid) > 100) browser()
      }
      
      if(!fix.gamma) {
        if(length(x$gamma$smooth.construct)) {
          for(sj in seq_along(x$gamma$smooth.construct)) {
            # cat("\n iteration: ", iter, ", predictor: gamma", ", term: ", sj)
            state <- update_jm_gamma(x$gamma$smooth.construct[[sj]], eta, eta_timegrid, int0,
                                     status, update.nu, criterion, get_LogPost, nobs, eps0_surv < edf.eps, edf = edf, ...)
            eta$gamma <- eta$gamma - fitted(x$gamma$smooth.construct[[sj]]$state) + fitted(state)
            edf <- edf - x$gamma$smooth.construct[[sj]]$state$edf + state$edf
            x$gamma$smooth.construct[[sj]]$state <- state
          }
        }
      }
      
      if(maxit[2] > 1) {
        eps00 <- do.call("cbind", eta[c("lambda", "gamma")])
        eps00 <- mean(abs((eps00 - eta00) / eps00), na.rm = TRUE)
        if(is.na(eps00) | !is.finite(eps00)) eps00 <- eps[2] + 1
      }
      
      iter00 <- iter00 + 1
    }
    
    eps0_surv <- do.call("cbind", eta[c("lambda", "gamma")])
    eps0_long <- rowsum(do.call("cbind", eta[c("mu", "sigma")]), id)
    eps0_alpha <- if(is.null(dalpha)) eta$alpha else do.call("cbind", eta[c("alpha", "dalpha")])
    
    eps0_surv <- if(!(fix.lambda & fix.gamma)) {
      mean(abs((eps0_surv - eta0_surv) / eps0_surv), na.rm = TRUE)
    } else 0
    eps0_long <- if(!(fix.mu & fix.sigma)) {
      mean(abs((eps0_long - eta0_long) / eps0_long), na.rm = TRUE)
    } else 0
    eps0_alpha <- if(!fix.alpha) {
      if(ok.alpha) {
        mean(abs((eps0_alpha - eta0_alpha) / eps0_alpha), na.rm = TRUE)
      } else 1 + eps[1]
    } else 0
    
    eps0 <- mean(c(eps0_surv, eps0_long, eps0_alpha))
    if(is.na(eps0) | !is.finite(eps0)) eps0 <- eps[1] + 1
    
    eeta <- exp(eta_timegrid)
    int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
    logLik <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
      exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
    edf <- get.edf(x, type = 2)
    IC <- get.ic2(logLik, edf, length(eta$mu), criterion)
    
    ic_contrib <- c(ic_contrib, IC)
    if(iter > nback) {
      adic <- abs(diff(tail(ic_contrib, nback))) / tail(ic_contrib, nback - 1)
      if(all(adic < ic.eps))
        eps0 <- 0
    }
    
    iter <- iter + 1
    
    if(verbose) {
      cat(if(ia) "\r" else "\n")
      vtxt <- paste(
        criterion, " ", fmt(IC, width = 8, digits = digits),
        " logLik ", fmt(logLik, width = 8, digits = digits),
        " edf ", fmt(edf, width = 6, digits = digits),
        " eps ", paste(fmt(eps0, width = 6, digits = digits + 2), "|S",
                       fmt(eps0_surv, width = 6, digits = digits + 2), "|L",
                       fmt(eps0_long, width = 6, digits = digits + 2), "|A",
                       fmt(eps0_alpha, width = 6, digits = digits + 2), sep = ""),
        " iteration ", formatC(iter, width = nchar(maxit[1])), sep = ""
      )
      cat(vtxt)
      if(.Platform$OS.type != "unix" & ia) flush.console()
    }
  }
  
  elapsed <- c(proc.time() - ptm)[3]
  
  if(verbose) {
    et <- if(elapsed > 60) {
      paste(formatC(format(round(elapsed / 60, 2), nsmall = 2), width = 5), "min", sep = "")
    } else paste(formatC(format(round(elapsed, 2), nsmall = 2), width = 5), "sec", sep = "")
    cat("\nelapsed time: ", et, "\n", sep = "")
  }
  
  eeta <- exp(eta_timegrid)
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  logLik <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  logPost <- as.numeric(logLik + get.log.prior(x))
  
  if(iter == maxit[1])
    warning("the backfitting algorithm did not converge, please check argument eps and maxit!")
  
  return(list("fitted.values" = eta, "parameters" = get.all.par(x),
              "logLik" = logLik, "logPost" = logPost, "hessian" = get.hessian(x),
              "converged" = iter < maxit[1], "time" = elapsed))
}


### Backfitting updating functions.
update_jm_gamma <- function(x, eta, eta_timegrid, int,
                            status, update.nu, criterion, get_LogPost, nobs, do.optim2, edf, ...)
{
  ## Gradient and hessian.
  tint <- xhess0 <- 0
  for(i in 1:nobs) {
    tint <- tint + exp(eta$gamma[i]) * x$X[i, , drop = FALSE] * int[i]
    xhess0 <- xhess0 + exp(eta$gamma[i]) * x$X[i, ] %*% t(x$X[i, ]) * int[i]
  }
  xgrad <- drop(t(status) %*% x$X - tint)
  
  env <- new.env()
  
  if((!(!x$state$do.optim | x$fixed | x$fxsp)) & do.optim2) {
    par <- x$state$parameters
    
    edf0 <- if(is.null(edf)) {
      0
    } else {
      edf - x$state$edf
    }
    
    objfun1 <- function(tau2) {
      par[x$pid$tau2] <- tau2
      xgrad <- xgrad + x$grad(score = NULL, par, full = FALSE)
      xhess <- xhess0 + x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(xhess, index = x$sparse.setup)
      Hs <- Sigma %*% xgrad
      g <- par[x$pid$b]
      if(update.nu) {
        objfun.nu <- function(nu) {
          g2 <- drop(g + nu * Hs)
          par[x$pid$b] <- g2
          fit <- x$fit.fun(x$X, g2)
          eta$gamma <- eta$gamma - fitted(x$state) + fit
          lp <- get_LogPost(eta_timegrid, eta, x$prior(par))
          return(-1 * lp)
        }
        nu <- optimize(f = objfun.nu, interval = c(0, 1))$minimum
      } else {
        nu <- x$state$nu
      }
      g2 <- drop(g + nu * Hs)
      fit <- x$fit.fun(x$X, g2)
      eta$gamma <- eta$gamma - fitted(x$state) + fit
      edf1 <- sum_diag(xhess0 %*% Sigma)
      edf <- edf0 + edf1
      
      logLik <- get_LogPost(eta_timegrid, eta, 0)
      ic <- get.ic2(logLik, edf, length(eta$mu), criterion)
      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          par[x$pid$b] <- g2
          opt_state <- list("parameters" = par,
                            "fitted.values" = fit, "edf" = edf1, "hessian" = xhess,
                            "nu" = nu, "do.optim" = x$state$do.optim)
          assign("state", opt_state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }
    
    assign("ic00_val", objfun1(get.state(x, "tau2")), envir = env)
    tau2 <- tau2.optim(objfun1, start = get.state(x, "tau2"))
    
    if(!is.null(env$state))
      return(env$state)
    
    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
  }
  
  xgrad <- drop(xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE))
  xhess <- xhess0 + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(xhess, index = NULL)
  Hs <- Sigma %*% xgrad
  
  ## Update regression coefficients.
  g <- get.state(x, "b")
  
  if(update.nu) {
    objfun2 <- function(nu) {
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      x$state$parameters <- set.par(x$state$parameters, g2, "b")
      
      ## Update additive predictors.
      fit <- x$fit.fun(x$X, g2)
      eta$gamma <- eta$gamma - fitted(x$state) + fit
      lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
      return(-1 * lp)
    }
    
    x$state$nu <- optimize(f = objfun2, interval = c(0, 1))$minimum
  }
  
  g2 <- drop(g + x$state$nu * Hs)
  names(g2) <- names(g)
  x$state$parameters <- set.par(x$state$parameters, g2, "b")
  
  ## Update additive predictors.
  fit <- x$fit.fun(x$X, g2)
  x$state$fitted.values <- fit
  x$state$edf <- sum_diag(xhess0 %*% Sigma)
  x$state$hessian <- xhess
  
  return(x$state)
}


update_jm_lambda <- function(x, eta, eta_timegrid,
                             eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                             status, update.nu, width, criterion, get_LogPost, nobs, do.optim2, edf, nonlinear, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  
  ## Timegrid predictor.
  eeta <- exp(eta_timegrid)
  
  ## Compute gradient and hessian integrals.
  # int <- survint(X, eeta, width, exp(eta$gamma), index = x$sparse.setup$matrix) # throws error with time-varying covariates
  int <- survint(X, eeta, width, exp(eta$gamma))
  # returns negative Hessian
  
  xgrad <- drop(t(status) %*% x$XT - int$grad)
  
  env <- new.env()
  
  if((!(!x$state$do.optim | x$fixed | x$fxsp)) & do.optim2) {
    par <- x$state$parameters
    
    edf0 <- if(is.null(edf)) {
      0
    } else {
      edf - x$state$edf
    }
    
    objfun1 <- function(tau2) {
      par <- set.par(par, tau2, "tau2")
      xgrad <- xgrad + x$grad(score = NULL, par, full = FALSE)
      xhess <- int$hess + x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(xhess, index = NULL)
      Hs <- Sigma %*% xgrad
      g <- get.par(par, "b")
      if(update.nu) {
        objfun.nu <- function(nu) {
          g2 <- drop(g + nu * Hs)
          names(g2) <- names(g)
          x$state$parameters <- set.par(x$state$parameters, g2, "b")
          x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
          
          ## Update additive predictors.
          fit_timegrid <- x$fit.fun_timegrid(g2)
          eta_timegrid_lambda <- eta_timegrid_lambda - x$state$fitted_timegrid + fit_timegrid
          if(nonlinear){
            eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
          } else {
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
          }          
          lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
          return(-1 * lp)
        }
        nu <- optimize(f = objfun.nu, interval = c(0, 1))$minimum
      } else {
        nu <- x$state$nu
      }
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      fit <- x$fit.fun(x$X, g2, expand = FALSE)
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta$lambda <- eta$lambda - fitted(x$state) + fit
      eta_timegrid_lambda <- eta_timegrid_lambda - x$state$fitted_timegrid + fit_timegrid
      if(nonlinear){
        eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
      } else {
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
      } 
      edf1 <- sum_diag(int$hess %*% Sigma)
      edf <- edf0 + edf1
      logLik <- get_LogPost(eta_timegrid, eta, 0)
      ic <- get.ic2(logLik, edf, length(eta$mu), criterion)
      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          opt_state <- list("parameters" = set.par(par, g2, "b"),
                            "fitted.values" = fit, "fitted_timegrid" = fit_timegrid,
                            "edf" = edf1, "hessian" = xhess,
                            "nu" = nu, "do.optim" = x$state$do.optim)
          assign("state", opt_state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }
    
    assign("ic00_val", objfun1(get.state(x, "tau2")), envir = env)
    tau2 <- tau2.optim(objfun1, start = get.state(x, "tau2"))
    # print(round(tau2, 4))
    
    if(!is.null(env$state))
      return(env$state)
    
    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
  }
  
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(xhess, index = NULL)
  Hs <- Sigma %*% xgrad
  
  ## Update regression coefficients.
  g <- get.state(x, "b")
  
  if(update.nu) {
    objfun2 <- function(nu) {
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      x$state$parameters <- set.par(x$state$parameters, g2, "b")
      
      ## Update additive predictors.
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta_timegrid_lambda <- eta_timegrid_lambda - x$state$fitted_timegrid + fit_timegrid
      if(nonlinear){
        eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
      } else {
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
      }      
      lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
      return(-1 * lp)
    }
    x$state$nu <- optimize(f = objfun2, interval = c(0, 1))$minimum
  }
  
  g2 <- drop(g + x$state$nu * Hs)
  names(g2) <- names(g)
  
  x$state$parameters <- set.par(x$state$parameters, g2, "b")
  
  ## Update additive predictors.
  x$state$fitted_timegrid <- x$fit.fun_timegrid(g2)
  x$state$fitted.values <- x$fit.fun(x$X, g2, expand = FALSE)
  x$state$edf <- sum_diag(int$hess %*% Sigma)
  x$state$hessian <- xhess
  
  return(x$state)
}


update_jm_mu <- function(x, y, eta, eta_timegrid,
                         eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                         status, update.nu, width, criterion, get_LogPost, nobs, do.optim2, edf, 
                         dx = NULL, xsmalpha = NULL, knots = NULL, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  dX <- if(!is.null(dx)) dx$fit.fun_timegrid(NULL) else NULL
  
  ## Timegrid predictor.
  eeta <- exp(eta_timegrid)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_alpha, width, exp(eta$gamma),
                 eeta * eta_timegrid_alpha^2, index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  
  xgrad <- drop(t(x$X) %*% drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2) +
                  t(x$XT) %*% drop(eta$alpha * status) - int$grad)
  if(!is.null(dx))
    xgrad <- drop(xgrad + t(dx$XT) %*% drop(eta$dalpha * status))
  XWX <- if(is.null(x$sparse.setup$matrix)) {
    crossprod(x$X * (1 / exp(eta$sigma)^2), x$X)
  } else do.XWX(x$X, exp(eta$sigma)^2, index = x$sparse.setup$matrix)
  xhess0 <- -1 * XWX - int$hess
  
  env <- new.env()
  
  if((!(!x$state$do.optim | x$fixed | x$fxsp)) & do.optim2) {
    par <- x$state$parameters
    
    edf0 <- if(is.null(edf)) {
      0
    } else {
      edf - x$state$edf
    }
    
    objfun1 <- function(tau2) {
      par <- set.par(par, tau2, "tau2")
      xgrad <- xgrad + x$grad(score = NULL, par, full = FALSE)
      # negative Hessian in xhess0
      xhess <- xhess0 - x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(-1 * xhess, index = NULL)
      Hs <- Sigma %*% xgrad
      g <- get.par(par, "b")
      if(update.nu) {
        objfun.nu <- function(nu) {
          g2 <- drop(g + nu * Hs)
          names(g2) <- names(g)
          x$state$parameters <- set.par(x$state$parameters, g2, "b")
          x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
          fit_timegrid <- x$fit.fun_timegrid(g2)
          eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
          if(!is.null(dx))
            eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
          fit <- x$fit.fun(x$X, g2)
          eta$mu <- eta$mu - fitted(x$state) + fit
          lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
          return(-1 * lp)
        }
        nu <- optimize(f = objfun.nu, interval = c(0, 1))$minimum
      } else {
        nu <- x$state$nu
      }
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      fit <- x$fit.fun(x$X, g2)
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta$mu <- eta$mu - fitted(x$state) + fit
      eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
      if(!is.null(dx))
        eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
      edf1 <- sum_diag((-1 * xhess0) %*% Sigma)
      edf <- edf0 + edf1
      logLik <- get_LogPost(eta_timegrid, eta, 0)
      ic <- get.ic2(logLik, edf, length(eta$mu), criterion)
      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          opt_state <- list("parameters" = set.par(par, g2, "b"),
                            "fitted.values" = fit, "fitted_timegrid" = fit_timegrid,
                            "edf" = edf1, "hessian" = -1 * xhess,
                            "nu" = nu, "do.optim" = x$state$do.optim)
          assign("state", opt_state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }
    
    assign("ic00_val", objfun1(get.state(x, "tau2")), envir = env)
    tau2 <- tau2.optim(objfun1, start = get.state(x, "tau2"), maxit = 1)
    
    if(!is.null(env$state))
      return(env$state)
    
    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
  }
  
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- xhess0 - x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(-1 * xhess, index = NULL)
  Hs <- Sigma %*% xgrad
  
  ## Update regression coefficients.
  g <- get.state(x, "b")
  
  if(update.nu) {
    objfun2 <- function(nu) {
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      x$state$parameters <- set.par(x$state$parameters, g2, "b")
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
      if(!is.null(dx))
        eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
      fit <- x$fit.fun(x$X, g2)
      eta$mu <- eta$mu - fitted(x$state) + fit
      lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
      return(-1 * lp)
    }
    x$state$nu <- optimize(f = objfun2, interval = c(0, 1))$minimum
  }
  
  g2 <- drop(g + x$state$nu * Hs)
  names(g2) <- names(g)
  x$state$parameters <- set.par(x$state$parameters, g2, "b")
  
  ## Update fitted values..
  x$state$fitted_timegrid <- x$fit.fun_timegrid(g2)
  x$state$fitted.values <- x$fit.fun(x$X, g2)
  x$state$edf <- sum_diag((-1 * xhess0) %*% Sigma)
  x$state$hessian <- -1 * xhess
  
  return(x$state)
}


update_jm_mu_Matrix <- function(x, y, eta, eta_timegrid,
                                eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                                status, update.nu, width, criterion, get_LogPost, nobs, do.optim2, edf, 
                                dx = NULL, xsmalpha = NULL, knots = NULL, tp = FALSE, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  dX <- if(!is.null(dx)) dx$fit.fun_timegrid(NULL) else NULL
  
  ## Timegrid predictor.
  eeta <- exp(eta_timegrid)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_alpha, width, exp(eta$gamma),
                 eeta * eta_timegrid_alpha^2, index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  
  xgrad <- crossprod(x$X, drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2)) +
    crossprod(x$XT,  drop(eta$alpha * status)) - int$grad
  if(!is.null(dx))
    xgrad <- xgrad + crossprod(dx$XT, drop(eta$dalpha * status))
  XWX <- crossprod(Diagonal(x = 1 / exp(eta$sigma)^2) %*% x$X, x$X)
  xhess0 <- -1 * XWX - int$hess
  
  env <- new.env()
  
  if((!(!x$state$do.optim | x$fixed | x$fxsp)) & do.optim2) {
    par <- x$state$parameters
    
    edf0 <- if(is.null(edf)) {
      0
    } else {
      edf - x$state$edf
    }
    
    objfun1 <- function(tau2) {
      par[x$pid$tau2] <- tau2
      xgrad <- xgrad + x$grad(score = NULL, par, full = FALSE)
      xhess <- xhess0 - x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(-1 * xhess, index = NULL)
      Hs <- Sigma %*% xgrad
      g <- par[x$pid$b]
      if(update.nu) {
        objfun.nu <- function(nu) {
          g2 <- drop(g + nu * Hs)
          par[x$pid$b] <- g2
          eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + x$fit.fun_timegrid(g2)
          if(!is.null(dx))
            eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
          fit <- x$fit.fun(x$X, g2)
          eta$mu <- eta$mu - fitted(x$state) + fit
          lp <- get_LogPost(eta_timegrid, eta, x$prior(par))
          return(-1 * lp)
        }
        nu <- optimize(f = objfun.nu, interval = c(0, 1))$minimum
      } else {
        nu <- x$state$nu
      }
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      fit <- x$fit.fun(x$X, g2)
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta$mu <- eta$mu - fitted(x$state) + fit
      eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
      if(!is.null(dx))
        eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
      edf1 <- sum_diag((-1 * xhess0) %*% Sigma)
      edf <- edf0 + edf1
      logLik <- get_LogPost(eta_timegrid, eta, 0)
      ic <- get.ic2(logLik, edf, length(eta$mu), criterion)
      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          par[x$pid$b] <- g2
          opt_state <- list("parameters" = par,
                            "fitted.values" = fit, "fitted_timegrid" = fit_timegrid,
                            "edf" = edf1, "hessian" = -1 * xhess,
                            "nu" = nu, "do.optim" = x$state$do.optim)
          assign("state", opt_state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }
    
    assign("ic00_val", objfun1(get.state(x, "tau2")), envir = env)
    tau2 <- tau2.optim(objfun1, start = get.state(x, "tau2"), maxit = 1)
    
    if(!is.null(env$state))
      return(env$state)
    
    x$state$parameters[x$pid$tau2] <- tau2
  }
  
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- xhess0 - x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(-1 * xhess, index = NULL)
  Hs <- Sigma %*% xgrad
  
  ## Update regression coefficients.
  g <- get.state(x, "b")
  
  if(update.nu) {
    objfun2 <- function(nu) {
      g2 <- drop(g + nu * Hs)
      x$state$parameters[x$pid$b] <- g2
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
      if(!is.null(dx))
        eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
      fit <- x$fit.fun(x$X, g2)
      eta$mu <- eta$mu - fitted(x$state) + fit
      lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
      return(-1 * lp)
    }
    x$state$nu <- optimize(f = objfun2, interval = c(0, 1))$minimum
  }
  
  g2 <- drop(g + x$state$nu * Hs)
  x$state$parameters[x$pid$b] <- g2
  
  ## Update fitted values..
  x$state$fitted_timegrid <- x$fit.fun_timegrid(g2)
  x$state$fitted.values <- x$fit.fun(x$X, g2)
  x$state$edf <- sum_diag((-1 * xhess0) %*% Sigma)
  x$state$hessian <- -1 * xhess
  
  return(x$state)
}


update_jm_dmu <- function(dx, x)
{
  b <- get.state(x, "b")
  state <- x$state
  state$fitted.values <- dx$fit.fun(dx$X, b)
  state$fitted_timegrid <- dx$fit.fun_timegrid(b)
  state$edf <- 0
  return(state)
}


update_jm_sigma <- function(x, y, eta, eta_timegrid,
                            update.nu, criterion, get_LogPost, nobs, do.optim2, edf, ...)
{
  xgrad <- crossprod(x$X, -1 + (y[, "obs"] - eta$mu)^2 / exp(eta$sigma)^2)
  xhess0 <- -2 * crossprod(x$X*drop((y[, "obs"] - eta$mu) / exp(eta$sigma)^2),
                           x$X*drop(y[, "obs"] - eta$mu))
  
  if((!(!x$state$do.optim | x$fixed | x$fxsp)) & do.optim2) {
    par <- x$state$parameters
    
    edf0 <- if(is.null(edf)) {
      0
    } else {
      edf - x$state$edf
    }
    
    env <- new.env()
    
    objfun1 <- function(tau2) {
      par <- set.par(par, tau2, "tau2")
      xgrad <- xgrad + x$grad(score = NULL, par, full = FALSE)
      xhess <- xhess0 - x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(-1 * xhess, index = NULL)
      Hs <- Sigma %*% xgrad
      g <- get.par(par, "b")
      if(update.nu) {
        objfun.nu <- function(nu) {
          g2 <- drop(g + nu * Hs)
          names(g2) <- names(g)
          x$state$parameters <- set.par(x$state$parameters, g2, "b")
          x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
          fit <- x$fit.fun(x$X, g2)
          eta$sigma <- eta$sigma - fitted(x$state) + fit
          lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
          return(-1 * lp)
        }
        nu <- optimize(f = objfun.nu, interval = c(0, 1))$minimum
      } else {
        nu <- x$state$nu
      }
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      fit <- x$fit.fun(x$X, g2)
      eta$sigma <- eta$sigma - fitted(x$state) + fit
      edf1 <- sum_diag((-1 * xhess0) %*% Sigma)
      edf <- edf0 + edf1
      logLik <- get_LogPost(eta_timegrid, eta, 0)
      ic <- get.ic2(logLik, edf, length(eta$mu), criterion)
      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          par[x$pid$b] <- g2
          opt_state <- list("parameters" = par,
                            "fitted.values" = fit,
                            "edf" = edf1, "hessian" = -1 * xhess,
                            "nu" = nu, "do.optim" = x$state$do.optim)
          assign("state", opt_state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }
    
    assign("ic00_val", objfun1(get.state(x, "tau2")), envir = env)
    tau2 <- tau2.optim(objfun1, start = get.state(x, "tau2"), maxit = 1)
    
    if(!is.null(env$state))
      return(env$state)
    
    x$state$parameters[x$pid$tau2] <- tau2
  }
  
  xhess <- xhess0 - x$hess(score = NULL, x$state$parameters, full = FALSE)
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(-1 * xhess, index = NULL)
  Hs <- Sigma %*% xgrad
  
  ## Update regression coefficients.
  g <- get.state(x, "b")
  
  if(update.nu) {
    objfun2 <- function(nu) {
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      x$state$parameters <- set.par(x$state$parameters, g2, "b")
      fit <- x$fit.fun(x$X, g2)
      eta$sigma <- eta$sigma - fitted(x$state) + fit
      lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
      return(-1 * lp)
    }
    
    x$state$nu <- optimize(f = objfun2, interval = c(0, 1))$minimum
  }
  
  
  g2 <- drop(g + x$state$nu * Hs)
  names(g2) <- names(g)
  x$state$parameters <- set.par(x$state$parameters, g2, "b")
  
  ## Update fitted values.
  x$state$fitted.values <- x$fit.fun(x$X, g2)
  x$state$edf <- sum_diag((-1 * xhess0) %*% Sigma)
  x$state$hessian <- -1 * xhess
  
  return(x$state)
}


update_jm_alpha <- function(x, eta, eta_timegrid,
                            eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                            status, update.nu, width, criterion, get_LogPost, nobs, do.optim2, edf, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  
  ## Timegrid predictor.
  eeta <- exp(eta_timegrid)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_mu, width, exp(eta$gamma),
                 eeta * (eta_timegrid_mu^2), index = x$sparse.setup$matrix)
  xgrad <- t(x$XT) %*% drop(eta_timegrid_mu[, ncol(eeta)] * status) - int$grad
  
  env <- new.env()
  
  if((!(!x$state$do.optim | x$fixed | x$fxsp)) & do.optim2) {
    par <- x$state$parameters
    
    edf0 <- if(is.null(edf)) {
      0
    } else {
      edf - x$state$edf
    }
    
    objfun1 <- function(tau2) {
      par <- set.par(par, tau2, "tau2")
      xgrad <- xgrad + x$grad(score = NULL, par, full = FALSE)
      xhess <- int$hess + x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(xhess, index = NULL)
      Hs <- Sigma %*% xgrad
      g <- get.par(par, "b")
      if(update.nu) {
        objfun.nu <- function(nu) {
          g2 <- drop(g + nu * Hs)
          names(g2) <- names(g)
          x$state$parameters <- set.par(x$state$parameters, g2, "b")
          x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
          fit_timegrid <- x$fit.fun_timegrid(g2)
          eta_timegrid_alpha <- eta_timegrid_alpha - x$state$fitted_timegrid + fit_timegrid
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
          lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
          return(-1 * lp)
        }
        nu <- optimize(f = objfun.nu, interval = c(0, 1))$minimum
      } else {
        nu <- x$state$nu
      }
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      fit <- x$fit.fun(x$X, g2)
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta$alpha <- eta$alpha - fitted(x$state) + fit
      eta_timegrid_alpha <- eta_timegrid_alpha - x$state$fitted_timegrid + fit_timegrid
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
      edf1 <- sum_diag(int$hess %*% Sigma)
      edf <- edf0 + edf1
      logLik <- get_LogPost(eta_timegrid, eta, 0)
      ic <- get.ic2(logLik, edf, length(eta$mu), criterion)
      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          opt_state <- list("parameters" = set.par(par, g2, "b"),
                            "fitted.values" = fit, "fitted_timegrid" = fit_timegrid,
                            "edf" = edf1, "hessian" = xhess,
                            "nu" = nu, "do.optim" = x$state$do.optim)
          assign("state", opt_state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }
    
    assign("ic00_val", objfun1(get.state(x, "tau2")), envir = env)
    tau2 <- tau2.optim(objfun1, start = get.state(x, "tau2"))
    
    if(!is.null(env$state))
      return(env$state)
    
    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
  }
  
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(xhess, x$sparse.setup)
  Hs <- Sigma %*% xgrad
  
  ## Update regression coefficients.
  g <- get.state(x, "b")
  
  if(update.nu) {
    objfun2 <- function(nu) {
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      x$state$parameters <- set.par(x$state$parameters, g2, "b")
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta_timegrid_alpha <- eta_timegrid_alpha - x$state$fitted_timegrid + fit_timegrid
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
      lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
      return(-1 * lp)
    }
    x$state$nu <- optimize(f = objfun2, interval = c(0, 1))$minimum
  }
  
  g2 <- drop(g + x$state$nu * Hs)        
  names(g2) <- names(g)
  x$state$parameters <- set.par(x$state$parameters, g2, "b")
  
  ## Update fitted values.
  x$state$fitted_timegrid <- x$fit.fun_timegrid(g2)
  x$state$fitted.values <- x$fit.fun(x$X, g2, expand = FALSE)
  x$state$edf <- sum_diag(int$hess %*% Sigma)

  x$state$hessian <- xhess
  
  return(x$state)
}


update_jm_dalpha <- function(x, eta, eta_timegrid,
                             eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                             status, update.nu, width, criterion, get_LogPost, nobs, do.optim2, edf, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  
  ## Timegrid predictor.
  eeta <- exp(eta_timegrid)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_dmu, width, exp(eta$gamma),
                 eeta * (eta_timegrid_dmu^2), index = x$sparse.setup$matrix)
  xgrad <- t(x$XT) %*% drop(eta_timegrid_dmu[, ncol(eeta)] * status) - int$grad
  
  env <- new.env()
  
  if((!(!x$state$do.optim | x$fixed | x$fxsp)) & do.optim2) {
    par <- x$state$parameters
    
    edf0 <- if(is.null(edf)) {
      0
    } else {
      edf - x$state$edf
    }
    
    objfun1 <- function(tau2) {
      par <- set.par(par, tau2, "tau2")
      xgrad <- xgrad + x$grad(score = NULL, par, full = FALSE)
      xhess <- int$hess + x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(xhess, index = NULL)
      Hs <- Sigma %*% xgrad
      g <- get.par(par, "b")
      if(update.nu) {
        objfun.nu <- function(nu) {
          g2 <- drop(g + nu * Hs)
          names(g2) <- names(g)
          x$state$parameters <- set.par(x$state$parameters, g2, "b")
          x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
          fit_timegrid <- x$fit.fun_timegrid(g2)
          eta_timegrid_dalpha <- eta_timegrid_dalpha - x$state$fitted_timegrid + fit_timegrid
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
          lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
          return(-1 * lp)
        }
        nu <- optimize(f = objfun.nu, interval = c(0, 1))$minimum
      } else {
        nu <- x$state$nu
      }
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      fit <- x$fit.fun(x$X, g2)
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta$dalpha <- eta$dalpha - fitted(x$state) + fit
      eta_timegrid_dalpha <- eta_timegrid_dalpha - x$state$fitted_timegrid + fit_timegrid
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
      edf1 <- sum_diag(int$hess %*% Sigma)
      edf <- edf0 + edf1
      logLik <- get_LogPost(eta_timegrid, eta, 0)
      ic <- get.ic2(logLik, edf, length(eta$mu), criterion)
      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          opt_state <- list("parameters" = set.par(par, g2, "b"),
                            "fitted.values" = fit, "fitted_timegrid" = fit_timegrid,
                            "edf" = edf1, "hessian" = xhess,
                            "nu" = nu, "do.optim" = x$state$do.optim)
          assign("state", opt_state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }
    
    assign("ic00_val", objfun1(get.state(x, "tau2")), envir = env)
    tau2 <- tau2.optim(objfun1, start = get.state(x, "tau2"))
    
    if(!is.null(env$state))
      return(env$state)
    
    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
  }
  
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(xhess, x$sparse.setup)
  Hs <- Sigma %*% xgrad
  
  ## Update regression coefficients.
  g <- get.state(x, "b")
  
  if(update.nu) {
    objfun2 <- function(nu) {
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      x$state$parameters <- set.par(x$state$parameters, g2, "b")
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta_timegrid_dalpha <- eta_timegrid_dalpha - x$state$fitted_timegrid + fit_timegrid
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
      lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
      return(-1 * lp)
    }
    x$state$nu <- optimize(f = objfun2, interval = c(0, 1))$minimum
  }
  nu <- x$state$nu
  
  g2 <- drop(g + nu * Hs)
  
  names(g2) <- names(g)
  x$state$parameters <- set.par(x$state$parameters, g2, "b")
  
  ## Update fitted values.
  x$state$fitted_timegrid <- x$fit.fun_timegrid(g2)
  x$state$fitted.values <- x$fit.fun(x$X, g2, expand = FALSE)
  x$state$edf <- sum_diag(int$hess %*% Sigma)
  x$state$hessian <- xhess
  
  return(x$state)
}


## (5) Joint model MCMC.
sam_JM <- jm_mcmc <- function(x, y, family, start = NULL, weights = NULL, offset = NULL,
  n.iter = 1200, burnin = 200, thin = 1, verbose = TRUE, 
  digits = 4, step = 20, ...)
{
  ## Hard coded.
  fixed <- NULL
  slice <- NULL
  
  dalpha <- has_pterms(x$dalpha$terms) | (length(x$dalpha$smooth.construct) > 0)
  
  nu <- 1
  
  if(!is.null(start)) {
    if(is.matrix(start)) {
      if(any(i <- grepl("Mean", colnames(start))))
        start <- start[, i]
      else stop("the starting values should be a vector not a matrix!")
    }
    if(dalpha) {
      if(!any(grepl("dmu.", names(start), fixed = TRUE))) {
        start.dmu <- start[grep("mu.", names(start), fixed = TRUE)]
        names(start.dmu) <- gsub("mu.", "dmu.", names(start.dmu), fixed = TRUE)
        start <- c(start, start.dmu)
      }
    }
    x <- set.starting.values(x, start)
  }
  
  ## Names of parameters/predictors.
  nx <- names(x)
  
  ## Extract y.
  y <- y[[1]]
  
  ## Number of observations.
  nobs <- attr(y, "nobs")
  
  ## Number of subdivions used for the time grid.
  sub <- attr(y, "subdivisions")
  
  ## The interval width from subdivisons.
  width <- attr(y, "width")
  
  ## Subject specific indicator
  take <- attr(y, "take")
  nlong <- length(take)
  
  ## Extract the status for individual i.
  status <- y[take, "status"]
  
  ## nonlinear setup
  nonlinear <- attr(y, "nonlinear")
  tp <- attr(y, "tp") 
  fac <- attr(y, "fac")
  
  ## Make id for individual i.
  id <- which(take)
  id <- append(id[-1], nlong + 1) - id 
  id <- rep(1:nobs, id)
  
  ## Compute additive predictors.
  eta <- get.eta(x, expand = FALSE) # potentially include alpha_mu
  
  ## Correct gamma predictor if NULL.
  if(!length(x$gamma$smooth.construct)) {
    eta$gamma <- rep(0, length = nobs)
  }
  
  ## For the time dependent part, compute
  ## predictors based on the time grid.
  eta_timegrid_mu <- 0
  if(length(x$mu$smooth.construct)) {
    for(j in names(x$mu$smooth.construct)) {
      b <- get.par(x$mu$smooth.construct[[j]]$state$parameters, "b")
      x$mu$smooth.construct[[j]]$state$fitted_timegrid <- x$mu$smooth.construct[[j]]$fit.fun_timegrid(b)
      eta_timegrid_mu <- eta_timegrid_mu + x$mu$smooth.construct[[j]]$state$fitted_timegrid
    }
  }
  
  
  eta_timegrid_alpha <- 0
  if(length(x$alpha$smooth.construct)) {
    for(j in names(x$alpha$smooth.construct)) {
      b <- get.par(x$alpha$smooth.construct[[j]]$state$parameters, "b")
      # update alpha according to mu
      if(nonlinear & !("model.matrix" %in% attr(x$alpha$smooth.construct[[j]], "class"))){
        Xmu <- as.vector(t(eta_timegrid_mu))
        if(!tp) {
          X <- x$alpha$smooth.construct[[j]]$update(Xmu, "mu")    
        } else {
          Xalpha <- x$alpha$smooth.construct[[j]]$margin[[1]]$update(Xmu, "mu")
          Xalpha2 <- x$alpha$smooth.construct[[j]]$margin[[2]]$fit.fun_timegrid(NULL)
          X <- rowTensorProduct(Xalpha, Xalpha2) 
        }
        fit_timegrid <- matrix(drop(X %*% b), nrow = nrow(eta_timegrid_mu), 
                               ncol = ncol(eta_timegrid_mu), byrow = TRUE)
        x$alpha$smooth.construct[[j]]$state$fitted.values <- fit_timegrid[, ncol(fit_timegrid)]
        x$alpha$smooth.construct[[j]]$state$fitted_timegrid <- fit_timegrid
        eta_timegrid_alpha <- eta_timegrid_alpha + fit_timegrid
      } else {
        x$alpha$smooth.construct[[j]]$state$fitted_timegrid <- x$alpha$smooth.construct[[j]]$fit.fun_timegrid(b)
        eta_timegrid_alpha <- eta_timegrid_alpha + x$alpha$smooth.construct[[j]]$fit.fun_timegrid(b)   
      }
      x$alpha$smooth.construct[[j]]$state$nu <- nu["alpha"]
    }
    eta$alpha <- eta_timegrid_alpha[, ncol(eta_timegrid_alpha)]
  }     
  
  
  eta_timegrid_lambda <- 0
  if(length(x$lambda$smooth.construct)) {
    for(j in names(x$lambda$smooth.construct)) {
      b <- get.par(x$lambda$smooth.construct[[j]]$state$parameters, "b")
      x$lambda$smooth.construct[[j]]$state$fitted_timegrid <- x$lambda$smooth.construct[[j]]$fit.fun_timegrid(b)
      eta_timegrid_lambda <- eta_timegrid_lambda + x$lambda$smooth.construct[[j]]$state$fitted_timegrid
    }
  }
  
  eta_timegrid_dmu <- eta_timegrid_dalpha <- 0
  if(dalpha) {
    if(length(x$dalpha$smooth.construct)) {
      for(j in names(x$dalpha$smooth.construct)) {
        b <- get.par(x$dalpha$smooth.construct[[j]]$state$parameters, "b")
        x$dalpha$smooth.construct[[j]]$state$fitted_timegrid <- x$dalpha$smooth.construct[[j]]$fit.fun_timegrid(b)
        eta_timegrid_dalpha <- eta_timegrid_dalpha + x$dalpha$smooth.construct[[j]]$fit.fun_timegrid(b)
      }
    }
    if(length(x$dmu$smooth.construct)) {
      for(j in names(x$dmu$smooth.construct)) {
        b <- get.par(x$mu$smooth.construct[[j]]$state$parameters, "b")
        x$dmu$smooth.construct[[j]]$state$fitted_timegrid <- x$dmu$smooth.construct[[j]]$fit.fun_timegrid(b)
        eta_timegrid_dmu <- eta_timegrid_dmu + x$dmu$smooth.construct[[j]]$fit.fun_timegrid(b)
      }
    }
  }
  
  # for nonlinear effect eta_mu is within eta_alpha
  if(nonlinear){
    eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
  } else{
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
  }
  knots <- if(nonlinear){
    if(tp){x$alpha$smooth.construct[[1]]$margin[[1]]$knots} 
    else {x$alpha$smooth.construct[[1]]$knots}} 
  else {NULL}
  myplots <- FALSE
  if(myplots){
    times <- matrix(unlist(attr(y, "grid")), ncol = ncol(eta_timegrid), byrow = TRUE)
    obstimes <- attr(y, "obstime")
    ids <- attr(y, "ids")
    par(mfrow = c(2,2))
  }
  
  ## Porcess iterations.
  if(burnin < 1) burnin <- 1
  if(burnin > n.iter) burnin <- floor(n.iter * 0.1)
  if(thin < 1) thin <- 1
  iterthin <- as.integer(seq(burnin, n.iter, by = thin))
  
  ## Samples.
  samps <- list()
  for(i in nx) {
    samps[[i]] <- list()
    for(j in names(x[[i]]$smooth.construct)) {
      samps[[i]][[j]] <- list(
        "samples" = matrix(NA, nrow = length(iterthin), ncol = length(x[[i]]$smooth.construct[[j]]$state$parameters)),
        "edf" = rep(NA, length = length(iterthin)),
        "alpha" = rep(NA, length = length(iterthin)),
        "accepted" = rep(NA, length = length(iterthin))
      )
      colnames(samps[[i]][[j]]$samples) <- names(x[[i]]$smooth.construct[[j]]$state$parameters)
    }
  }
  logLik.samps <- logPost.samps <- rep(NA, length = length(iterthin))
  
  foo <- function(x) {
    if(is.null(x)) return(0)
    if(is.na(x)) return(0)
    x <- exp(x)
    if(x < 0)
      x <- 0
    if(x > 1)
      x <- 1
    x
  }
  
  ## Integrals.
  eeta <- exp(eta_timegrid)
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  
  ## Slice sampling?
  if(is.null(slice)) {
    slice <- c("lambda" = FALSE, "gamma" = FALSE, "mu" = FALSE,
               "sigma" = FALSE, "alpha" = FALSE, "dmu" = FALSE, "dalpha" = FALSE)
  } else {
    if(is.null(names(slice))) {
      slice <- rep(slice, length.out = 7)
      names(slice) <- c("lambda", "gamma", "mu", "sigma", "alpha", "dmu", "dalpha")
    } else {
      npar <- c("lambda", "gamma", "mu", "sigma", "alpha", "dmu", "dalpha")
      if(!all(j <- npar %in% names(slice))) {
        nnot <- npar[!j]
        for(j in npar) {
          if(j %in% nnot) {
            nslice <- c(names(slice), j)
            slice <- c(slice, FALSE)
            names(slice) <- nslice
          }
        }
      }
    }
  }
  
  nx2 <- if(dalpha) {
    c("lambda", "gamma", "mu", "sigma", "alpha", "dalpha")
  } else {
    c("lambda", "gamma", "mu", "sigma", "alpha")
  }
  
  ## Start sampling.
  if(verbose) {
    cat2("Starting the sampler...")
    if(!interactive())
      cat("\n")
  }
  
  nstep <- step
  step <- floor(n.iter / step)
  
  ptm <- proc.time()
  for(iter in 1:n.iter) {
    if(save <- iter %in% iterthin)
      js <- which(iterthin == iter)
    
    for(i in nx2) {
      if(i == "gamma") {
        eeta <- exp(eta_timegrid)
        int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
      }
      
      prop_fun <- get_jm_prop_fun(i, slice[i], nonlinear)
      # 3jump
      for(sj in names(x[[i]]$smooth.construct)) {
        
        p.state <- prop_fun(x[[i]]$smooth.construct[[sj]],
                            y, eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                            eta_timegrid_dmu, eta_timegrid_dalpha,
                            width, sub, nu, status, id = i, int0, nobs,
                            dx = if(dalpha & (i == "mu")) x[["dmu"]]$smooth.construct[[sj]] else NULL,
                            xsmalpha = if(nonlinear & (i == "mu")) x$alpha$smooth.construct else NULL, 
                            knots = if(nonlinear & (i == "mu")) knots else NULL, tp = tp, fac = fac, ...)
        ## If accepted, set current state to proposed state.
        accepted <- if(is.na(p.state$alpha)) FALSE else log(runif(1)) <= p.state$alpha
        
        #checks# if(i == "mu") cat("\n", sj, round(exp(p.state$alpha), 2))
        if(i %in% fixed)
          accepted <- FALSE
        
        if(accepted) {
          if(i %in% c("lambda", "mu", "alpha", "dalpha")) {
            if(i == "lambda")
              eta_timegrid_lambda <- eta_timegrid_lambda - x[[i]]$smooth.construct[[sj]]$state$fitted_timegrid + p.state$fitted_timegrid
            if(i == "mu") {
              # cat("\n iteration: ", iter, ", predictor: ", i, ", term: ", sj)
              eta_timegrid_mu <- eta_timegrid_mu - x[[i]]$smooth.construct[[sj]]$state$fitted_timegrid + p.state$fitted_timegrid
              if(nonlinear){
                for (j in names(x$alpha$smooth.construct)){
                  if(j != "model.matrix"){   # only smooth.constructs need to be updated
                    g_a <- get.par(x$alpha$smooth.construct[[j]]$state$parameters, "b")
                    Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
                    Xalpha <- constrain(x$alpha$smooth.construct[[j]], Xalpha)
                    if(fac)
                      Xalpha <- Xalpha * x$alpha$smooth.construct[[j]]$by_timegrid
                    if(tp)  
                      Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
                    alpha_state <- matrix(Xalpha %*% g_a, nrow = nrow(eta_timegrid), ncol = ncol(eta_timegrid), byrow = TRUE)
                    eta_timegrid_alpha <- eta_timegrid_alpha - x$alpha$smooth.construct[[j]]$state$fitted_timegrid + alpha_state
                    x$alpha$smooth.construct[[j]]$state$fitted_timegrid <- alpha_state
                    eta$alpha <- eta$alpha - fitted(x$alpha$smooth.construct[[j]]$state) + alpha_state[, ncol(eta_timegrid)]
                    x$alpha$smooth.construct[[j]]$state$fitted.values <- alpha_state[, ncol(eta_timegrid)]
                  } 
                }
              }
              if(dalpha & (sj %in% names(x[["dmu"]]$smooth.construct))) {
                p.state.dmu <- update_jm_dmu(x[["dmu"]]$smooth.construct[[sj]], x[["mu"]]$smooth.construct[[sj]])
                eta_timegrid_dmu <- eta_timegrid_dmu - x[["dmu"]]$smooth.construct[[sj]]$state$fitted_timegrid + p.state.dmu$fitted_timegrid
                eta[["dmu"]] <- eta[["dmu"]] - fitted(x[["dmu"]]$smooth.construct[[sj]]$state) + fitted(p.state.dmu)
                x[["dmu"]]$smooth.construct[[sj]]$state <- p.state.dmu
              }
            }
            if(i == "alpha"){
              eta_timegrid_alpha <- eta_timegrid_alpha - x[[i]]$smooth.construct[[sj]]$state$fitted_timegrid + p.state$fitted_timegrid
            }
            if(i == "dalpha")
              eta_timegrid_dalpha <- eta_timegrid_dalpha - x[[i]]$smooth.construct[[sj]]$state$fitted_timegrid + p.state$fitted_timegrid
            if(nonlinear){
              eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
            } else {
            eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
          }
          }
          eta[[i]] <- eta[[i]] - fitted(x[[i]]$smooth.construct[[sj]]$state) + fitted(p.state)
          if(i == "alpha"){
            if(myplots & iter == round(iter, -1)){
              par(mfrow=c(2,2))
              plot(eta_timegrid_mu, eta_timegrid_alpha, main=paste0("alpha: iteration ", iter, ", edf ", round(p.state$edf, 2)))
              plot(eta_timegrid_mu[, ncol(eta_timegrid_mu)], eta$alpha, main=paste0("alpha: iteration ", iter, ", edf ", round(p.state$edf, 2)))
              plot(eta_timegrid_mu, p.state$fitted_timegrid, 
                   main=paste0("alpha: iteration ", iter, ", edf ", round(p.state$edf, 2), ", term: ", sj))
              abline(h = 0, col = 2)
              abline(v = median(y[, 3]), col = 2)
              matplot(t(times), t(eta_timegrid_alpha), main=paste0("alpha: iteration ", iter, ", edf ", round(p.state$edf, 2)), type = "l")
             Sys.sleep(2)
              }
          }
            if(i == "mu"){
              if(myplots & iter == round(iter, -1)){
                par(mfrow=c(2,2))
                matplot(t(times),t(p.state$fitted_timegrid), type = "l", main = paste0("mu: iteration ", iter, ", edf ", round(p.state$edf, 2), " ", sj))
                matplot(t(times), t(eta_timegrid_mu), type = "l", main = paste0("mu: iteration ", iter, ", edf ", round(p.state$edf, 2), " ", sj))
                plot(eta_timegrid_mu, eta_timegrid_alpha, main=paste0("mu: iteration ", iter, ", edf ", round(p.state$edf, 2), "alpha"))
                plot(eta_timegrid_mu[, ncol(eta_timegrid_mu)], eta$alpha, main=paste0("mu: iteration ", iter, ", edf ", round(p.state$edf, 2), "alpha"))
                Sys.sleep(1)
              }
            }
          x[[i]]$smooth.construct[[sj]]$state <- p.state 
        }
        
        ## Save the samples and acceptance.
        if(save) {
          samps[[i]][[sj]]$samples[js, ] <- x[[i]]$smooth.construct[[sj]]$state$parameters
          samps[[i]][[sj]]$edf[js] <- x[[i]]$smooth.construct[[sj]]$state$edf
          samps[[i]][[sj]]$alpha[js] <- foo(p.state$alpha)
          samps[[i]][[sj]]$accepted[js] <- accepted
          if(dalpha & (i == "mu") & (sj %in% names(x[["dmu"]]$smooth.construct))) {
            samps[["dmu"]][[sj]]$samples[js, ] <- x[["dmu"]]$smooth.construct[[sj]]$state$parameters
            samps[["dmu"]][[sj]]$edf[js] <- x[["dmu"]]$smooth.construct[[sj]]$state$edf
            samps[["dmu"]][[sj]]$alpha[js] <- foo(p.state$alpha)
            samps[["dmu"]][[sj]]$accepted[js] <- accepted
          }
        }
      }
    }
    
    if(save) {
      logLik.samps[js] <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
        exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
      logPost.samps[js] <- as.numeric(logLik.samps[js] + get.log.prior(x))
    }
    
    if(verbose) barfun(ptm, n.iter, iter, step, nstep)
  }
  
  if(verbose) cat("\n")
  
  for(i in names(samps)) {
    for(j in names(samps[[i]])) {
      samps[[i]][[j]] <- do.call("cbind", samps[[i]][[j]])
      cn <- if(j == "model.matrix") {
        paste(i, "p", j, colnames(samps[[i]][[j]]), sep = ".")
      } else {
        paste(i, "s", j, colnames(samps[[i]][[j]]), sep = ".")
      }
      colnames(samps[[i]][[j]]) <- cn
    }
    samps[[i]] <- do.call("cbind", samps[[i]])
  }
  samps$logLik <- logLik.samps
  samps$logPost <- logPost.samps
  samps <- do.call("cbind", samps)
  
  ## Compute DIC. #
  dev <- -2 * logLik.samps
  mpar <- apply(samps, 2, mean, na.rm = TRUE)
  names(mpar) <- colnames(samps)
  ll <- family$p2logLik(mpar)
  mdev <- -2 * ll
  pd <- mean(dev) - mdev
  DIC <- mdev + 2 * pd
  samps <- cbind(samps,
                 "DIC" = rep(DIC, length.out = nrow(samps)),
                 "pd" = rep(pd, length.out = nrow(samps))
  )
  
  return(as.mcmc(samps))
}


## Get proposal function.
get_jm_prop_fun <- function(i, slice = FALSE, nonlinear = FALSE) {
  function(...) {
    prop_fun <- if(slice) {
      propose_jm_slice
    } else {
      if(nonlinear & (i == "alpha")){
        get("propose_jm_alpha_nonlin")
      } else {
      get(paste("propose_jm", i, sep = "_"))
    }
    }
    state <- prop_fun(...)
    if(inherits(state, "try-error")) {
      warning(paste("problems sampling parameter ", i, "!", sep = ""))
      return(list("alpha" = -Inf))
    } else {
      return(state)
    }
  }
}


## JM proposal functions.
uni.slice_beta_logPost <- function(g, x, family, y = NULL, eta = NULL, id,
                                   eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                                   eta_timegrid_dmu, eta_timegrid_dalpha,
                                   width, sub, status, dx = NULL, nonlinear = FALSE, ...)
{
  if(id %in% c("lambda", "mu", "alpha", "dalpha")) {
    if(id == "lambda")
      eta_timegrid_lambda <- eta_timegrid_lambda - x$state$fitted_timegrid + x$fit.fun_timegrid(g)
    if(id == "mu") {
      eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + x$fit.fun_timegrid(g)
      if(!is.null(dx)) {
        eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g)
      }
    }
    if(id == "alpha")
      eta_timegrid_alpha <- eta_timegrid_alpha - x$state$fitted_timegrid + x$fit.fun_timegrid(g)
    if(id == "dalpha")
      eta_timegrid_dalpha <- eta_timegrid_dalpha - x$state$fitted_timegrid + x$fit.fun_timegrid(g)
    if(nonlinear){
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
    } else {
    eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
  }
  }
  
  eta[[id]] <- eta[[id]] - fitted(x$state) + x$fit.fun(x$X, g, expand = FALSE)
  
  eeta <- exp(eta_timegrid)
  int <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  
  ll <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  lp <- x$prior(g)
  
  return(ll + lp)
}


uni.slice_tau2_logPost <- function(g, x, family, y = NULL, eta = NULL, id, ll = 0)
{
  lp <- x$prior(g)
  return(ll + lp)
}


propose_jm_slice <- function(x, y,
                             eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                             eta_timegrid_dmu, eta_timegrid_dalpha,
                             width, sub, nu, status, id, dx = NULL, nonlinear = FALSE, ...)
{
  for(j in seq_along(get.par(x$state$parameters, "b"))) {
    x$state$parameters <- uni.slice(x$state$parameters, x, family, y,
                                    eta, id = id, j, logPost = uni.slice_beta_logPost, eta_timegrid = eta_timegrid,
                                    eta_timegrid_lambda = eta_timegrid_lambda, eta_timegrid_mu = eta_timegrid_mu,
                                    eta_timegrid_alpha = eta_timegrid_alpha, eta_timegrid_dmu = eta_timegrid_dmu,
                                    eta_timegrid_dalpha = eta_timegrid_dalpha, width = width, sub = sub, status = status,
                                    dx = dx, nonlinear = nonlinear)
  }
  
  g <- get.par(x$state$parameters, "b")
  fit <- drop(x$X %*% g)
  eta[[id]] <- eta[[id]] - fitted(x$state) + fit
  x$state$fitted.values <- fit
  
  if(id %in% c("lambda", "mu", "alpha", "dalpha"))
    x$state$fitted_timegrid <- x$fit.fun_timegrid(g)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if((length(x$S) < 2) & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(g, x$S[[1]]) %*% g + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "mu", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- 0
  
  return(x$state)
}


propose_jm_lambda <- function(x, y,
                              eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                              eta_timegrid_dmu, eta_timegrid_dalpha,
                              width, sub, nu, status, id, ...)
{
  nonlinear <- attr(y, "nonlinear")
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  
  ## Timegrid predictor.
  eeta <- exp(eta_timegrid)
  
  ## Old logLik and prior.
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibeta <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  p1 <- x$prior(x$state$parameters)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta, width, exp(eta$gamma), index = x$sparse.setup$matrix)
  xgrad <- drop(t(status) %*% x$XT - int$grad)
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(xhess, index = NULL)
  
  ## Save old coefficients.
  g0 <- get.state(x, "b")
  
  ## Get new position.
  mu <- drop(g0 + nu * Sigma %*% xgrad)
  
  ## Sample new parameters.
  g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma, method="chol"))
  names(g) <- names(g0)
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  
  ## Compute log priors.
  p2 <- x$prior(x$state$parameters)
  qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)
  
  ## Update additive predictors.
  fit_timegrid <- x$fit.fun_timegrid(g)
  eta_timegrid_lambda <- eta_timegrid_lambda - x$state$fitted_timegrid + fit_timegrid
  if(nonlinear){
    eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
  } else {
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
  }
  x$state$fitted_timegrid <- fit_timegrid
  
  fit <- drop(x$X %*% g)
  eta$lambda <- eta$lambda - fitted(x$state) + fit
  x$state$fitted.values <- fit
  
  ## New logLik.
  eeta <- exp(eta_timegrid)
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibetaprop <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  
  ## Prior prob.
  int <- survint(X, eeta, width, exp(eta$gamma), index = x$sparse.setup$matrix)
  xgrad <- drop(t(status) %*% x$XT - int$grad)
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  Sigma2 <- matrix_inv(xhess, index = NULL)
  mu2 <- drop(g + nu * Sigma2 %*% xgrad)
  qbeta <- dmvnorm(g0, mean = mu2, sigma = Sigma2, log = TRUE)
  
  ## Save edf.
  x$state$edf <- sum_diag(int$hess %*% Sigma2)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if((length(x$S) < 2) & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(g, x$S[[1]]) %*% g + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "mu", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  # cat(paste("\n",pibeta, pibetaprop, qbeta, qbetaprop, p1, p2, x$state$alpha, sep=";"))
  return(x$state)
}


propose_jm_mu <- function(x, ...)
{
  return(x$propose(x, ...))
}

propose_jm_mu_simple <- function(x, y,
                                 eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                                 eta_timegrid_dmu, eta_timegrid_dalpha,
                                 width, sub, nu, status, id, dx = NULL, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  dX <- if(!is.null(dx)) dx$fit.fun_timegrid(NULL) else NULL
  
  ## Timegrid lambda.
  eeta <- exp(eta_timegrid)
  
  ## Old logLik and prior.
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibeta <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  p1 <- x$prior(x$state$parameters)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_alpha, width, exp(eta$gamma),
                 eeta * eta_timegrid_alpha^2, index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  
  xgrad <- drop(t(x$X) %*% drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2) +
                  t(x$XT) %*% drop(eta$alpha * status) - int$grad)
  if(!is.null(dx))
    xgrad <- drop(xgrad + t(dx$XT) %*% drop(eta$dalpha * status))
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  XWX <- if(is.null(x$sparse.setup$matrix)) {
    crossprod(x$X * (1 / exp(eta$sigma)^2), x$X)
  } else do.XWX(x$X, exp(eta$sigma)^2, index = x$sparse.setup$matrix)
  xhess <- -1 * XWX - int$hess
  xhess <- xhess - x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(-1 * xhess, index = NULL)
  
  ## Save old coefficients.
  g0 <- get.state(x, "b")
  
  ## Get new position.
  mu <- drop(g0 + nu * Sigma %*% xgrad)
  
  ## Sample new parameters.
  g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma, method="chol"))
  names(g) <- names(g0)
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  
  ## Compute log priors.
  p2 <- x$prior(x$state$parameters)
  qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)
  
  ## Update additive predictors.
  fit_timegrid <- x$fit.fun_timegrid(g)
  eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
  if(!is.null(dx))
    eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g)
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
  x$state$fitted_timegrid <- fit_timegrid
  
  fit <- drop(x$X %*% g)
  eta$mu <- eta$mu - fitted(x$state) + fit
  x$state$fitted.values <- fit
  
  ## New logLik.
  eeta <- exp(eta_timegrid)
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibetaprop <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  
  ## Prior prob.
  int <- survint(X, eeta * eta_timegrid_alpha, width, exp(eta$gamma),
                 eeta * eta_timegrid_alpha^2, index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  
  xgrad <- drop(t(x$X) %*% drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2) +
                  t(x$XT) %*% drop(eta$alpha * status) - int$grad)
  if(!is.null(dx))
    xgrad <- drop(xgrad + t(dx$XT) %*% drop(eta$dalpha * status))
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  XWX <- if(is.null(x$sparse.setup$matrix)) {
    crossprod(x$X * (1 / exp(eta$sigma)^2), x$X)
  } else do.XWX(x$X, exp(eta$sigma)^2, index = x$sparse.setup$matrix)
  xhess <- -1 * XWX - int$hess
  xhess <- xhess - x$hess(score = NULL, x$state$parameters, full = FALSE)
  Sigma2 <- matrix_inv(-1 * xhess, index = NULL)
  mu2 <- drop(g + nu * Sigma2 %*% xgrad)
  qbeta <- dmvnorm(g0, mean = mu2, sigma = Sigma2, log = TRUE)
  
  ## Save edf.
  x$state$edf <- sum_diag((-1 * int$hess) %*% Sigma2)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if((length(x$S) < 2) & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(g, x$S[[1]]) %*% g + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "mu", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  
  return(x$state)
}


propose_jm_mu_Matrix <- function(x, y,
                                 eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                                 eta_timegrid_dmu, eta_timegrid_dalpha,
                                 width, sub, nu, status, id, dx = NULL, ...)
{
  
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  dX <- if(!is.null(dx)) dx$fit.fun_timegrid(NULL) else NULL
  
  ## Timegrid lambda.
  eeta <- exp(eta_timegrid)
  
  ## Old logLik and prior.
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibeta <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  p1 <- x$prior(x$state$parameters)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_alpha, width, exp(eta$gamma),
                 eeta * eta_timegrid_alpha^2, index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  
  xgrad <- crossprod(x$X, drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2)) +
    crossprod(x$XT,  drop(eta$alpha * status)) - int$grad
  if(!is.null(dx))
    xgrad <- xgrad + crossprod(dx$XT, drop(eta$dalpha * status))
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  XWX <- crossprod(Diagonal(x = 1 / exp(eta$sigma)^2) %*% x$X, x$X)
  xhess0 <- -1 * XWX - int$hess
  xhess <- xhess0 - x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(-1 * xhess, index = NULL)
  
  ## Save old coefficients.
  g0 <- get.state(x, "b")
  
  ## Get new position.
  mu <- drop(g0 + Sigma %*% xgrad)
  Sigma <- as.matrix(Sigma)
  
  ## Sample new parameters using blockdiagonal structure.
  if(!is.null(x$sparse.setup$block.index)){
    lg <- lapply(1:length(x$sparse.setup$block.index), function(i){
      tmp <- x$sparse.setup$block.index[[i]]
      if(x$sparse.setup$is.diagonal){
        if(Sigma[tmp,tmp] < 0) Sigma[tmp,tmp] <- -1 * Sigma[tmp,tmp]
        drop(rnorm(n = 1, mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp])))
      } else{
        if(length(tmp) == 1){
          if(Sigma[tmp,tmp] < 0) Sigma[tmp,tmp] <- -1 * Sigma[tmp,tmp]
          drop(rnorm(n = 1, mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp])))
        } else {
          drop(rmvnorm(n = 1, mean = mu[tmp], sigma = Sigma[tmp,tmp], method="chol"))
        }
      } 
    })
    g <- unlist(lg)
  } else {
    g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma, method="chol"))
  }
  
  if(all(is.na(g))) {
    x$state$alpha <- 1
    return(x$state)
  }
  
  names(g) <- names(g0)
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  
  ## Compute log priors using blockdiagonal structure.
  p2 <- x$prior(x$state$parameters)
  if(!is.null(x$sparse.setup$block.index)){
    lqbetaprop <- lapply(1:length(x$sparse.setup$block.index), function(i){
      tmp <- x$sparse.setup$block.index[[i]]
      if(x$sparse.setup$is.diagonal){
        if(Sigma[tmp,tmp] < 0) Sigma[tmp,tmp] <- -1 * Sigma[tmp,tmp]
        drop(dnorm(g[tmp], mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp]), log = TRUE))
      } else{
        if(length(tmp) == 1){
          if(Sigma[tmp,tmp] < 0) Sigma[tmp,tmp] <- -1 * Sigma[tmp,tmp]
          drop(dnorm(g[tmp], mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp]), log = TRUE))
        } else {
          dmvnorm(g[tmp], mean = mu[tmp], sigma = Sigma[tmp,tmp], log = TRUE)
        }
      } 
    })
    qbetaprop <- sum(unlist(lqbetaprop))
  } else {
    qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)
  }
  
  ## Update additive predictors.
  fit_timegrid <- x$fit.fun_timegrid(g)
  eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
  if(!is.null(dx))
    eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g)
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
  x$state$fitted_timegrid <- fit_timegrid
  
  fit <- drop(x$X %*% g)
  eta$mu <- eta$mu - fitted(x$state) + fit
  x$state$fitted.values <- fit
  
  ## New logLik.
  eeta <- exp(eta_timegrid)
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibetaprop <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  
  ## Prior prob.
  int <- survint(X, eeta * eta_timegrid_alpha, width, exp(eta$gamma),
                 eeta * eta_timegrid_alpha^2, index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  
  xgrad <- crossprod(x$X, drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2)) +
    crossprod(x$XT,  drop(eta$alpha * status)) - int$grad
  if(!is.null(dx))
    xgrad <- xgrad + crossprod(dx$XT, drop(eta$dalpha * status))
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  XWX <- crossprod(Diagonal(x = 1 / exp(eta$sigma)^2) %*% x$X, x$X)
  xhess0 <- -1 * XWX - int$hess
  xhess <- xhess0 - x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  Sigma2 <- matrix_inv(-1 * xhess, index = NULL)
  mu2 <- drop(g + nu * Sigma2 %*% xgrad)
  Sigma2 <- as.matrix(Sigma2) 
  
  if(!is.null(x$sparse.setup$block.index)){
    lqbeta <- lapply(1:length(x$sparse.setup$block.index), function(i){
      tmp <- x$sparse.setup$block.index[[i]]
      if(x$sparse.setup$is.diagonal){
        if(Sigma2[tmp,tmp] < 0) Sigma2[tmp,tmp] <- -1 * Sigma2[tmp,tmp]
        drop(dnorm(g0[tmp], mean = mu2[tmp], sd = sqrt(Sigma2[tmp,tmp]), log = TRUE))
      } else{
        if(length(tmp) == 1){
          if(Sigma2[tmp,tmp] < 0) Sigma2[tmp,tmp] <- -1 * Sigma2[tmp,tmp]
          drop(dnorm(g0[tmp], mean = mu2[tmp], sd = sqrt(Sigma2[tmp,tmp]), log = TRUE))
        } else {
          dmvnorm(g0[tmp], mean = mu2[tmp], sigma = Sigma2[tmp,tmp], log = TRUE)
        }
      } 
    })
    qbeta <- sum(unlist(lqbeta))
  } else {
    qbeta <- dmvnorm(g0, mean = mu2, sigma = as.matrix(Sigma2), log = TRUE)
  }
  
  
  ## Save edf.
  x$state$edf <- sum_diag((-1 * int$hess) %*% Sigma2)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if((length(x$S) < 2) & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- drop(0.5 * crossprod(g, x$S[[1]]) %*% g + x$b)
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "mu", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  
  return(x$state)
}


propose_jm_alpha <- function(x, y,
                             eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                             eta_timegrid_dmu, eta_timegrid_dalpha,
                             width, sub, nu, status, id, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  
  ## Timegrid lambda.
  eeta <- exp(eta_timegrid)
  
  ## Old logLik and prior.
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibeta <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  p1 <- x$prior(x$state$parameters)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_mu, width, exp(eta$gamma),
                 eeta * eta_timegrid_mu^2, index = x$sparse.setup$matrix)
  xgrad <- t(x$XT) %*% drop(eta_timegrid_mu[, ncol(eeta)] * status) - int$grad
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(xhess, index = NULL)
  
  ## Save old coefficients.
  g0 <- get.state(x, "b")
  
  ## Get new position.
  mu <- drop(g0 + nu * Sigma %*% xgrad)
  
  ## Sample new parameters.
  g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma, method="chol"))
  names(g) <- names(g0)
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  
  ## Compute log priors.
  p2 <- x$prior(x$state$parameters)
  qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)
  
  ## Update additive predictors.
  fit_timegrid <- x$fit.fun_timegrid(g)
  eta_timegrid_alpha <- eta_timegrid_alpha - x$state$fitted_timegrid + fit_timegrid
  x$state$fitted_timegrid <- fit_timegrid
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
  
  fit <- x$fit.fun(x$X, g, expand = FALSE)
  eta$alpha <- eta$alpha - fitted(x$state) + fit
  x$state$fitted.values <- fit
  
  ## New logLik.
  eeta <- exp(eta_timegrid)
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibetaprop <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  
  ## Prior prob.
  int <- survint(X, eeta * eta_timegrid_mu, width, exp(eta$gamma),
                 eeta * eta_timegrid_mu^2, index = x$sparse.setup$matrix)
  xgrad <- t(x$XT) %*% drop(eta_timegrid_mu[, ncol(eeta)] * status) - int$grad
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  Sigma2 <- matrix_inv(xhess, index = NULL)
  mu2 <- drop(g + nu * Sigma2 %*% xgrad)
  qbeta <- dmvnorm(g0, mean = mu2, sigma = Sigma2, log = TRUE)
  
  ## Save edf.
  x$state$edf <- sum_diag(int$hess %*% Sigma2)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if((length(x$S) < 2) & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(g, x$S[[1]]) %*% g + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "alpha", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  return(x$state)
}


propose_jm_dalpha <- function(x, y,
                              eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                              eta_timegrid_dmu, eta_timegrid_dalpha,
                              width, sub, nu, status, id, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  
  ## Timegrid lambda.
  eeta <- exp(eta_timegrid)
  
  ## Old logLik and prior.
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibeta <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  p1 <- x$prior(x$state$parameters)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_dmu, width, exp(eta$gamma),
                 eeta * eta_timegrid_dmu^2, index = x$sparse.setup$matrix)
  xgrad <- t(x$XT) %*% drop(eta_timegrid_dmu[, ncol(eeta)] * status) - int$grad
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(xhess, index = NULL)
  
  ## Save old coefficients.
  g0 <- get.state(x, "b")
  
  ## Get new position.
  mu <- drop(g0 + nu * Sigma %*% xgrad)
  
  ## Sample new parameters.
  g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma, method="chol"))
  names(g) <- names(g0)
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  
  ## Compute log priors.
  p2 <- x$prior(x$state$parameters)
  qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)
  
  ## Update additive predictors.
  fit_timegrid <- x$fit.fun_timegrid(g)
  eta_timegrid_dalpha <- eta_timegrid_dalpha - x$state$fitted_timegrid + fit_timegrid
  x$state$fitted_timegrid <- fit_timegrid
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha * eta_timegrid_mu + eta_timegrid_dalpha * eta_timegrid_dmu
  
  fit <- x$fit.fun(x$X, g, expand = FALSE)
  eta$alpha <- eta$alpha - fitted(x$state) + fit
  x$state$fitted.values <- fit
  
  ## New logLik.
  eeta <- exp(eta_timegrid)
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibetaprop <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  
  ## Prior prob.
  int <- survint(X, eeta * eta_timegrid_dmu, width, exp(eta$gamma),
                 eeta * eta_timegrid_dmu^2, index = x$sparse.setup$matrix)
  xgrad <- t(x$XT) %*% drop(eta_timegrid_dmu[, ncol(eeta)] * status) - int$grad
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  Sigma2 <- matrix_inv(xhess, index = NULL)
  mu2 <- drop(g + nu * Sigma2 %*% xgrad)
  qbeta <- dmvnorm(g0, mean = mu2, sigma = Sigma2, log = TRUE)
  
  ## Save edf.
  x$state$edf <- sum_diag(int$hess %*% Sigma2)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if((length(x$S) < 2) & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(g, x$S[[1]]) %*% g + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "alpha", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  return(x$state)
}


propose_jm_gamma <- function(x, y,
                             eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                             eta_timegrid_dmu, eta_timegrid_dalpha,
                             width, sub, nu, status, id, int, nobs, ...)
{
  ## Timegrid lambda.
  eeta <- exp(eta_timegrid)
  
  ## Old logLik and prior.
  pibeta <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  p1 <- x$prior(x$state$parameters)
  
  ## Compute gradient and hessian integrals.
  tint <- xhess <- 0
  for(i in 1:nobs) {
    tint <- tint + exp(eta$gamma[i]) * x$X[i, , drop = FALSE] * int[i]
    xhess <- xhess + exp(eta$gamma[i]) * x$X[i, ] %*% t(x$X[i, ]) * int[i]
  }
  xgrad <- t(status) %*% x$X - tint
  xgrad <- drop(xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE))
  xhess <- xhess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(xhess, index = NULL)
  
  ## Save old coefficients.
  g0 <- get.state(x, "b")
  
  ## Get new position.
  mu <- drop(g0 + nu * Sigma %*% xgrad)
  
  ## Sample new parameters.
  g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma, method="chol"))
  names(g) <- names(g0)
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  
  ## Compute log priors.
  p2 <- x$prior(x$state$parameters)
  qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)
  
  ## Update additive predictors.
  fit <- x$fit.fun(x$X, g, expand = FALSE)
  eta$gamma <- eta$gamma - fitted(x$state) + fit
  x$state$fitted.values <- fit
  
  ## New logLik.
  pibetaprop <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  
  ## Prior prob.
  tint <- xhess <- 0
  for(i in 1:nobs) {
    tint <- tint + exp(eta$gamma[i]) * x$X[i, , drop = FALSE] * int[i]
    xhess <- xhess + exp(eta$gamma[i]) * x$X[i, ] %*% t(x$X[i, ]) * int[i]
  }
  xgrad <- t(status) %*% x$X - tint
  xgrad <- drop(xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE))
  xhess0 <- xhess
  xhess <- xhess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  
  Sigma2 <- matrix_inv(xhess, index = NULL)
  mu2 <- drop(g + nu * Sigma2 %*% xgrad)
  qbeta <- dmvnorm(g0, mean = mu2, sigma = Sigma2, log = TRUE)
  
  ## Save edf.
  x$state$edf <- sum_diag(xhess0 %*% Sigma2)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if((length(x$S) < 2) & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(g, x$S[[1]]) %*% g + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "gamma", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  
  return(x$state)
}


propose_jm_sigma <- function(x, y,
                             eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                             eta_timegrid_dmu, eta_timegrid_dalpha,
                             width, sub, nu, status, id, int, ...)
{
  ## Timegrid lambda.
  eeta <- exp(eta_timegrid)
  
  ## Old logLik and prior.
  pibeta <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  p1 <- x$prior(x$state$parameters)
  
  ## Compute gradient and hessian integrals.
  xgrad <- crossprod(x$X, -1 + (y[, "obs"] - eta$mu)^2 / exp(eta$sigma)^2)
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  
  xhess <- -2*crossprod(x$X*drop((y[, "obs"] - eta$mu)/exp(eta$sigma)^2),
                        x$X*drop(y[, "obs"] - eta$mu))  
  xhess <- xhess - x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(-1 * xhess, index = NULL)
  
  ## Save old coefficients.
  g0 <- get.state(x, "b")
  
  ## Get new position.
  mu <- drop(g0 + nu * Sigma %*% xgrad)
  
  ## Sample new parameters.
  g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma, method="chol"))
  names(g) <- names(g0)
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  
  ## Compute log priors.
  p2 <- x$prior(x$state$parameters)
  qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)
  
  ## Update additive predictors.
  fit <- x$fit.fun(x$X, g, expand = FALSE)
  eta$sigma <- eta$sigma - fitted(x$state) + fit
  x$state$fitted.values <- fit
  
  ## New logLik.
  pibetaprop <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  
  ## Prior prob.
  xgrad <- crossprod(x$X, -1 + (y[, "obs"] - eta$mu)^2 / exp(eta$sigma)^2)
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  
  xhess <- xhess0 <- -2*crossprod(x$X*drop((y[, "obs"] - eta$mu)/exp(eta$sigma)^2),
                                  x$X*drop(y[, "obs"] - eta$mu))  
  xhess <- xhess - x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  Sigma2 <- matrix_inv(-1 * xhess, index = NULL)
  mu2 <- drop(g + nu * Sigma2 %*% xgrad)
  qbeta <- dmvnorm(g0, mean = mu2, sigma = Sigma2, log = TRUE)
  
  ## Save edf.
  x$state$edf <- sum_diag((-1 * xhess0) %*% Sigma2)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if((length(x$S) < 2) & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(g, x$S[[1]]) %*% g + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "sigma", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  
  return(x$state)
}


## Time varying random slopes.
smooth.construct.Re.smooth.spec <- function(object, data, knots)
{
  isf <- sapply(data[object$term], is.factor)
  id <- data[[object$term[isf]]]
  if(object$bs.dim < 0)
    object$bs.dim <- 5
  xobj <- eval(as.call(c(as.symbol("s"),
                         as.symbol(object$term[!isf]),
                         k=object$bs.dim,xt=list(object$xt),
                         bs="ps")))
  xobj <- smoothCon(xobj, data, knots, absorb.cons = FALSE, n = length(id))[[1]]
  xl <- levels(id)
  object$X <- list()
  for(j in seq_along(xl)) {
    object$X[[j]] <- xobj$X[id == xl[j], , drop = FALSE]
  }
  object$X <- as.matrix(do.call("bdiag", object$X))
  object$zeros <- apply(object$X, 2, function(x) { all(x == 0) })
  object$X <- object$X[, !object$zeros]
  object$isf <- isf
  object$bs.dim <- ncol(object$X)
  object$S <- list(diag(object$bs.dim))
  object$rank <- object$bs.dim
  object$null.space.dim <- 0
  object$C <- matrix(0, 0, ncol(object$X))
  object$side.constrain <- FALSE
  object$plot.me <- TRUE
  object$te.ok <- 2
  object$random <- TRUE
  object$xt$force.center <- TRUE
  object$xobj <- xobj
  class(object) <- "Random.effect"
  object
}

Predict.matrix.Random.effect <- function(object, data) 
{
  id <- data[[object$term[object$isf]]]
  Xd <- PredictMat(object$xobj, data, n = length(id))
  X <- list()
  xl <- levels(id)
  for(j in seq_along(xl))
    X[[j]] <- Xd[id == xl[j], , drop = FALSE]
  X <- as.matrix(do.call("bdiag", X))
  X <- X[, !object$zeros]
  X
}

center.X <- function(X, S)
{
  QR <- qr(crossprod(X, rep(1, length = nrow(X))))
  Q2 <- qr.Q(QR, complete = TRUE)[, -1]
  X <- X %*% Q2
  S <- crossprod(Q2, S) %*% Q2
  return(list("X" = X, "S" = S, "Q2" = Q2))
}


smooth.construct.Re2.smooth.spec <- function(object, data, knots)
{
  isf <- sapply(data[object$term], is.factor)
  id <- data[[object$term[isf]]]
  if(object$bs.dim < 0)
    object$bs.dim <- 5
  xobj <- eval(as.call(c(as.symbol("s"),
                         as.symbol(object$term[!isf]),
                         k=object$bs.dim,xt=list(object$xt),
                         bs="ps")))
  xobj <- smoothCon(xobj, data, knots, absorb.cons = FALSE, n = length(id))[[1]]
  xl <- levels(id)
  object$X <- object$S <- object$Q2 <- list()
  for(j in seq_along(xl)) {
    Xt <- xobj$X[id == xl[j], , drop = FALSE]
    Xt <- center.X(Xt, diag(ncol(Xt)))
    object$X[[j]] <- Xt$X
    object$S[[j]] <- Xt$S
    object$Q2[[j]] <- Xt$Q2
  }
  object$X <- as.matrix(do.call("bdiag", object$X))
  object$S <- as.matrix(do.call("bdiag", object$S))
  object$zeros <- apply(object$X, 2, function(x) { all(x == 0) })
  object$X <- object$X[, !object$zeros]
  object$S <- list(object$S[!object$zeros, !object$zeros])
  object$isf <- isf
  object$bs.dim <- ncol(object$X)
  object$rank <- object$bs.dim
  object$null.space.dim <- 0
  object$C <- matrix(0, 0, ncol(object$X))
  object$side.constrain <- FALSE
  object$plot.me <- TRUE
  object$te.ok <- 2
  object$random <- TRUE
  object$xt$force.center <- TRUE
  object$xobj <- xobj
  class(object) <- "Random2.effect"
  object
}

Predict.matrix.Random2.effect <- function(object, data) 
{
  id <- data[[object$term[object$isf]]]
  Xd <- PredictMat(object$xobj, data, n = length(id))
  X <- list()
  xl <- levels(id)
  for(j in seq_along(xl)) {
    X[[j]] <- Xd[id == xl[j], , drop = FALSE]
    X[[j]] <- X[[j]] %*% object$Q2[[j]]
  }
  X <- as.matrix(do.call("bdiag", X))
  X <- X[, !object$zeros]
  X
}


# time transform for time-varying survival covariates
param_time_transform2 <- function (x, formula, data, grid, yname, timevar, take, 
                                   derivMat = FALSE, eps = 1e-07, 
                                   timevar2=NULL, idvar=NULL, delete.intercept = FALSE){
  if (derivMat) 
    ddata <- data
  id <- data[[idvar]]
  X <- Xn <- tvar <- NULL
  for (j in names(data)) {
    if ((!grepl("Surv(", j, fixed = TRUE) & !grepl("Surv2(", j, fixed = TRUE)) & 
        (j != yname) & (j != timevar)) {
      # split data per subject
      idata <- split(data[[j]], id)
      # check if timevarying variable
      temp <- lapply(1:length(idata), function(i){length(unique(idata[[i]])) > 1})
      if(any(unlist(temp))){
        tvar <- c(tvar, j)
        # extract unique time-varying values
        values <- lapply(1:length(idata), function(i){unique(idata[[i]])})
        # extract break points
        breaks <- lapply(1:length(idata), function(i){
          split(data[[timevar2]], id)[[i]][c(TRUE, diff(idata[[i]]) != 0)]})
        # transfer break points to evaluation grid
        igrid <- lapply(1:length(idata), function(i){
          if(length(breaks[[i]]) > 1){
            g <- cut(grid[[i]], breaks[[i]], labels=FALSE, include.lowest = TRUE)
            g[is.na(g)] <- max(g, na.rm=TRUE) + 1
            g
          } else {
            rep(1, length(grid[[i]]))
          }})
        # evaluate variable on that grid
        evalgrid <- lapply(1:length(idata), function(i){
          values[[i]][igrid[[i]]]})
        df <- data.frame(unlist(evalgrid))
        names(df) <- j
        X <- if (is.null(X)) 
          df
        else cbind(X, df)
        Xn <- c(Xn, j)
      }
    }
  }
  if (!is.null(take)) 
    data <- data[take, , drop = FALSE]
  for (j in names(data)) {
    if ((!grepl("Surv(", j, fixed = TRUE) & !grepl("Surv2(", 
                                                   j, fixed = TRUE)) & (j != yname) & (j != timevar) & !(j %in% tvar)) {
      df <- data.frame(rep(data[[j]], each = length(grid[[1]])))
      names(df) <- j
      X <- if (is.null(X)) 
        df
      else cbind(X, df)
      Xn <- c(Xn, j)
    }
  }
  if (!is.null(X)) 
    colnames(X) <- Xn
  X <- if (is.null(X)) 
    data.frame(unlist(grid))
  else cbind(X, unlist(grid))
  colnames(X)[ncol(X)] <- yname
  if (timevar != yname) {
    X <- cbind(X, unlist(grid))
    colnames(X)[ncol(X)] <- timevar
  }
  dX <- NULL
  if (derivMat) {
    dX <- X
    for (j in colnames(dX)) {
      if (!is.factor(dX[[j]]) & (grepl(timevar, j, fixed = TRUE)) & 
          (timevar %in% c(x$term, x$by))) {
        dX[[j]] <- dX[[j]] + eps
        ddata[[j]] <- ddata[[j]] + eps
      }
    }
  }

  if(delete.intercept){
    formula <- update(formula, ~ . + 1)
    X <- model.matrix(formula, data = X)
    if("(Intercept)" %in% colnames(X))
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  formula <- update(formula, ~ . - 1)
  } else {
  X <- model.matrix(formula, data = X)
  }
  

  gdim <- c(length(grid), length(grid[[1]]))
  x$XT <- extract_XT(X, gdim[1], gdim[2])
  if (derivMat) {
    dX <- model.matrix(formula, data = dX)
    X <- -1 * (X - dX)/eps
    x$XT <- extract_XT(X, gdim[1], gdim[2])
    x$X <- -1 * (x$X - model.matrix(formula, data = ddata))/eps
  }
  x$fit.fun_timegrid <- function(g) {
    if (is.null(g)) 
      return(X)
    g <- get.par(g, "b")
    f <- drop(X %*% g)
    f <- matrix(f, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    f
  }
  x$state$fitted_timegrid <- x$fit.fun_timegrid(get.state(x, "b"))
  class(x) <- c(class(x), "deriv.model.matrix")
  x
}

# simulate data
################################################################################
simJM <- function(nsub = 300, times = seq(0, 120, 1), probmiss = 0.75,
                  long_setting = "functional", alpha_setting = if(nonlinear) "linear" else "nonlinear", 
                  dalpha_setting = "zero", sigma = 0.3, long_df = 6, tmax=NULL,    
                  seed = NULL, full = FALSE, file = NULL, nonlinear = FALSE, fac = FALSE){
  
  if(is.null(tmax)){
    tmax <- max(times)
  }
  
  if(nonlinear & alpha_setting == "nonlinear2") fac <- TRUE
  if(nonlinear & fac & alpha_setting == "nonlinear") alpha_setting <- "nonlinear2"
  
  
  ## specify censoring function (aus CoxFlexBoost) changed into uniformly (as too much censoring)
  ## added censoring at tmax
  censoring <- function(time, tmax){
    ## censoring times are independent uniformly distributed
    censor_time <- runif(n = length(time), min=0, max=1.5*tmax)
    censor_time <- ifelse(censor_time > tmax, tmax, censor_time)
    event <- (time <= censor_time)
    survtime <- apply(cbind(time, censor_time), 1, min)
    ## return matrix of observed survival times and event indicator
    return(cbind(survtime, event))
  }
  
  
  ## introduce random missings 
  ## (excluding the first measurement to ensure all subjects remain in the data)
  miss_fct <- function(data, prop, obstime="obstime"){
    select <- which(data[[obstime]] > 0) 
    n <- length(select)
    n_miss <- round(prop*n, 0)
    miss <- sample(select, n_miss)
    data <- data[-miss,]
    return(data)
  }
  
  
  ## generate baseline covariates 
  ## (x1=survival covariate, x2=longitudinal covariate)
  gen_x <- function(nsub){
    x1 <- runif(nsub, -3, 3) 
    x2 <- runif(nsub, -3, 3) 
    x3 <- rbinom(nsub, 1, 0.5)
    cbind(x1, x2, x3)
  } 
  
  ## generate random effects 
  ## (r1=random intercept, r2=random slope)
  gen_r <- function(nsub){
    r1 <- r2 <- matrix(NA, nrow=nsub, ncol=2)
    r1 <- rnorm(nsub, 0, 0.25) 
    r2 <- rnorm(nsub, 0, 0.4) 
    cbind(r1, r2)
  }  
  
  ## function written by Fabian Scheipl for random functional effect
  ## (changed into random functional intercepts)
  gen_b <- function(times, nsub, long_df, pen = 2, l = c(1,1), seed = NULL) {
    if(!is.null(seed)) set.seed(seed)
    # Recursion for difference operator matrix
    makeDiffOp <- function(degree, dim){
      if(degree == 0){
        return(diag(dim))  
      } else {
        return(diff(makeDiffOp(degree - 1, dim)))
      }    
    }
    # Kronecker Product penalty from marginal
    Pi <- l[1] * kronecker(diag(nsub), diag(long_df)) 
    Pt <- l[2] * kronecker(diag(nsub), crossprod(makeDiffOp(pen[1], long_df)))
    P <- .1*diag(nsub*long_df) + Pt + Pi
    
    coef <- matrix(rmvnorm(nsub*long_df, sigma = solve(P), method="chol"), ncol = long_df, nrow = nsub)
    colnames(coef) <- paste0("b", 1:long_df)
    bt <- bs(times, df = long_df, intercept = FALSE)
    b_set <- list(knots = attr(bt, "knots"), Boundary.knots = attr(bt, "Boundary.knots"),
                  degree = attr(bt, "degree"), intercept = attr(bt, "intercept"))
    return(list(coef, b_set))
  }
  
  ## numerical derivative for B-Spline
  ## modified from JMbayes code
  dbs <-  function (x, df = NULL, knots = NULL, degree=3, intercept = FALSE, 
                    Boundary.knots = range(x), eps = 1e-07) {
    ex <- pmax(abs(x), 1)
    x1 <- x + eps * ex
    x2 <- x - eps * ex
    bs.xeps1 <- suppressWarnings(bs(x1, df, knots, degree, intercept, Boundary.knots))
    bs.xeps2 <- suppressWarnings(bs(x2, df, knots, degree, intercept, Boundary.knots))
    out <- (bs.xeps1 - bs.xeps2) / c(x1 - x2)
    out
  }
  
  
  ## compute predictors
  ## individual longitudinal trajectories
  mu <-  function(time, x, r, long_df, b_set, long_setting){
    
    # allow to always extract coefficients as column
    if(is.null(dim(r))){
      r <- matrix(r, nrow=1)
    } 
    
    beta <- r[, -c(1:2)]
    # duplicate vector beta for the multiple integration points
    if(is.null(dim(beta))){
      beta <- matrix(beta, nrow = length(time), ncol = long_df, byrow=TRUE)
    }
    # TODO: Suppress warnings
    switch(long_setting,
           "linear" = (1.25 + r[, 1] + 0.6*sin(x) + (-0.01)*time + r[, 2]*0.02*time),
           "nonlinear" = (0.5 + r[, 1] + 0.6*sin(x) + 0.1*(time+2)*exp(-0.075*time)),
           "functional" = (0.5 + r[, 1] + 0.6*sin(x) + 0.1*(time+2)*exp(-0.075*time) + 
                             apply(bs(time, long_df, b_set$knots, b_set$degree,
                                               b_set$intercept, b_set$Boundary.knots) * beta, 1, sum)))
  }
  
  ## derivative of individual longitudinal trajectories
  dmu <-  function(time, r, long_df, b_set, long_setting){
    
    # allow to always extract coefficients as column
    if(is.null(dim(r))){
      r <- matrix(r, nrow = 1)
    } 
    
    beta <- r[, -c(1:2)]
    # duplicate vector beta for the multiple integration points
    if(is.null(dim(beta))){
      beta <- matrix(beta, nrow = length(time), ncol = long_df, byrow = TRUE)
    }
    switch(long_setting,
           "simple" = (-0.02) + 0*time,
           "linear" = (-0.01 + r[, 2]*0.02) + 0*time,
           "nonlinear" = (0.085 - 0.0075*time)*exp(-0.075*time),    
           "functional" = (0.085 - 0.0075*time)*exp(-0.075*time) + 
             apply(dbs(time, long_df, b_set$knots, b_set$degree,
                       b_set$intercept, b_set$Boundary.knots) * beta,1,sum))
  }
  
  ## association between mu and log-hazard
  alpha <- function(time, alpha_setting){
    switch(alpha_setting,
           "zero" = 0*time,
           "constant" = 0*time + 1,
           "linear" = 1 - 0.015*time,
           "nonlinear" = cos((time-60)/20) + 1,
           "nonlinear2" = cos((time-33)/33))
  }
  
  alpha_nonlin <- function(time, alpha_setting, x2, x3, r, long_df, b_set, long_setting){
    switch(alpha_setting,
           "zero" = 0*time,
           "constant" = 0*time + 1,
           "linear" = 1 * mu(time, x2, r, long_df, b_set, long_setting),
           "nonlinear" = -0.1 * (mu(time, x2, r, long_df, b_set, long_setting) + 3)^2 +
             mu(time, x2, r, long_df, b_set, long_setting) + 1.8,
           "nonlinear2" = x3 * (-0.1 * (mu(time, x2, r, long_df, b_set, long_setting) + 3)^2 +
                                  mu(time, x2, r, long_df, b_set, long_setting) + 1.8) + 
                          (1 - x3) * (0.1 * (mu(time, x2, r, long_df, b_set, long_setting) - 3)^2 +
                                  0.75 * mu(time, x2, r, long_df, b_set, long_setting) - 0.8))
  }
  
  # used only for setting up final data
  alpha_nonlin_simple <- function(alpha_setting, x3, mu){
    switch(alpha_setting,
           "zero" = 0,
           "constant" = 1,
           "linear" = 1 * mu,
           "nonlinear" = -0.1 * (mu + 3)^2 + mu + 1.8,
           "nonlinear2" = x3 * (-0.1 * (mu + 3)^2 +  mu + 1.8) + 
                          (1 - x3) * (0.1 * (mu - 3)^2 + 0.75 * mu - 0.8))
  }
  
  ## association between dmu and log-hazard
  dalpha <- function(time, dalpha_setting){
    switch(dalpha_setting,
           "zero" = 0*time,
           "constant" = 0*time + 10,
           "linear" = 6 - 0.015*time,
           "nonlinear" = 50 + 55*sin((time)/20),
           "nonlinear2" = 50 + 55*sin((time)/20))
  }
  
  ## baseline hazard
  lambda <-  function(time){
    1.4*log((time + 10)/1000) - 1.5
  }
  
  ## baseline covariate
  gamma <-  function(x1, nonlinear = TRUE){
    if(nonlinear){
      # 0.3*x1 + 0.1*(1 - x3) 
      0.3*x1 
    } else {
      sin(x1)
    }
  }
  
  ## full hazard (including nonlinear effect)
  hazard <-  function(time, x, r, ...){
    if(nonlinear){
      exp(lambda(time) + gamma(x[1], nonlinear) + 
            alpha_nonlin(time, alpha_setting, x[2], x[3], r, long_df, b_set, long_setting) +
            dalpha(time, dalpha_setting)*dmu(time, r, long_df, b_set, long_setting))
    } else{
    exp(lambda(time) + gamma(x[1]) + 
          alpha(time, alpha_setting)*mu(time, x[2], r, long_df, b_set, long_setting) + 
          dalpha(time, dalpha_setting)*dmu(time, r, long_df, b_set, long_setting))
  }
  }
  
  
  # generate input
  id <- rep(1:nsub, each=length(times))
  if(!is.null(seed)){
    set.seed(seed)
  }
  r <- gen_r(nsub)
  x <- gen_x(nsub)
  if(!fac) x[, 3] <- rep(1, nsub)
  
  temp <- gen_b(times, nsub, long_df=long_df, l=c(1,5))
  r <- cbind(r, temp[[1]])
  b_set <- temp[[2]]
  
  data_short <- rJM(hazard, censoring, x, r, tmin = times[1], tmax = tmax) 
  
  data_long <- cbind(id, data_short[id,], obstime = rep(times, nsub))
  data_grid <- data.frame(survtime = times,
                          mu = seq(-0.5, 2.5, length.out = length(times)))
  
  i <- !duplicated(data_long$id)
  
  # Save true predictors
  data_long$dalpha <- dalpha(data_long$survtime, dalpha_setting) 
  data_grid$dalpha <- dalpha(data_grid$survtime, dalpha_setting) 
  
  # gamma and lambda have only joint intercept which is estimated in predictor gamma
  data_long$mu <- mu(data_long$obstime, data_long$x2, r[id,], long_df, b_set, long_setting)
  f_lambda <- lambda(data_long$survtime)[i]   
  f_gamma <- gamma(data_long$x1, nonlinear)[i]
  data_long$lambda <- lambda(data_long$survtime) - mean(f_lambda)
  data_grid$lambda <- lambda(data_grid$survtime) - mean(f_lambda)
  if(nonlinear){
    data_long$gamma <- gamma(data_long$x1, nonlinear) + mean(f_lambda)
    data_long$surv_mu <- mu(data_long$survtime, data_long$x2, r[id,], long_df, b_set, long_setting)
    data_long$alpha <- alpha_nonlin_simple(alpha_setting, data_long$x3, data_long$surv_mu)
    data_long$alpha_l <- alpha_nonlin_simple(alpha_setting, data_long$x3, data_long$mu)
    data_grid$mu <- seq(-0.5, 2.5, length.out = nrow(data_grid))
    if(fac){
      data_grid$alpha1 <- alpha_nonlin_simple(alpha_setting, rep(1, nrow(data_grid)), data_grid$mu)
      data_grid$alpha0 <- alpha_nonlin_simple(alpha_setting, rep(0, nrow(data_grid)), data_grid$mu)
    } else {
      data_grid$alpha <- alpha_nonlin_simple(alpha_setting, rep(1, nrow(data_grid)), data_grid$mu)
    }

  } else {
    data_long$alpha <- alpha(data_long$survtime, alpha_setting)
    data_grid$alpha <- alpha(data_grid$survtime, alpha_setting)
    data_long$gamma <- gamma(data_long$x1, nonlinear) + mean(f_lambda)
  }
  
  data_long$dmu <- dmu(data_long$obstime, r[id,], long_df, b_set, long_setting)
  data_long$id <- as.factor(data_long$id)
  data_long$sigma <- rep(log(sigma), nrow(data_long))
  data_long$x3 <- as.factor(data_long$x3)
  
  # censoring                   
  data_long <- data_long[data_long$obstime <= data_long$survtime,]
  
  # saving data without longitudinal missings
  data_full <- data_long
  
  # inducing longitudinal missings
  data_long <- miss_fct(data_long, probmiss)
  
  # Draw longitudinal observations
  data_long$y <- rnorm(nrow(data_long), data_long$mu, sigma)
  
  ygrid <- quantile(data_long$y, probs = seq(0.025, 0.975, 0.025))
  # adjust predictions with constraint median(y)
  if(nonlinear){
   # constraint acts on group 0. to be checked
    alpha_constraint <- mean(alpha_nonlin_simple(alpha_setting, 0, ygrid))
    data_long$alpha <- data_long$alpha - alpha_constraint
    data_long$alpha_l <- data_long$alpha_l - alpha_constraint
    if(fac){
      data_grid$alpha1 <- data_grid$alpha1 - alpha_constraint
      data_grid$alpha0 <- data_grid$alpha0 - alpha_constraint
    } else {
      data_grid$alpha <- data_grid$alpha - alpha_constraint
    }
    data_long$gamma <- data_long$gamma + alpha_constraint
    data_long$alpha_constraint <- alpha_constraint
  } 
  
  if(full){
    d <- list(data=data_long, data_grid = data_grid, data_full = data_full)
  } else {
    d <- data_long
  }
  if(!is.null(file)) { 
    save(d, file = file) 
    invisible(d) 
  } else { 
    return(d) 
  } 
}


rJM <- function(hazard, censoring, x, r, 
                subdivisions = 1000, tmin = 0, tmax, file = NULL, ...){
  ## compute hazard for every person i at time
  nsub <- nrow(x)
  
  time <- rep(NA, nsub) 
  
  Hazard <- function(hazard, time, x, r) { 
    integrate(hazard, 0, time, x = x, r = r,
              subdivisions = subdivisions)$value 
  } 
  
  InvHazard <- function(Hazard, hazard, x, r, tmin, tmax) { 
    negLogU <- -log(runif(1, 0, 1)) 
    # check if first Lambda value is smaller than sample
    rootfct <- function(time) { 
      negLogU - Hazard(hazard, time, x, r) 
    } 
    if(rootfct(tmin)<0){
      return(0)
    } else {
      root <- try(uniroot(rootfct, interval = c(0, tmax))$root, silent=TRUE)
      root <- if(inherits(root, "try-error")) {
        # if root not within [0, tmax] --> error --> set it to tmax + 0.01 (will be censored)
        tmax + 0.01
      } else {root}
    }
    return(root)
  }
  
  # Finding Survival Times
  cumhaz <- rep(NA, nsub)
  survprob <- rep(NA, nsub)
  for(i in 1:nsub) { 
    time[i] <- InvHazard(Hazard, hazard, x[i,], r[i,], tmin, tmax)
    cumhaz[i] <- Hazard(hazard, time[i], x[i,], r[i,])
    survprob[i] <- exp((-1)*cumhaz[i])
  } 
  
  time_event <- censoring(time, tmax) 
  
  # Make data (long format)
  data_short <- data.frame(survtime = time_event[, 1], event = time_event[, 2], 
                           x, r, cumhaz = cumhaz)
  names(data_short) <- gsub(".", "", names(data_short), fixed = TRUE) 
  
  return(data_short)
}


## Prediction.
jm_predict <- function(object, newdata, type = c("link", "parameter", "probabilities", "cumhaz", "loglik"),
                       dt, steps, id, FUN = function(x) { mean(x, na.rm = TRUE) }, subdivisions = 100, cores = NULL,
                       chunks = 1, verbose = FALSE,  ...)
{
  if(attr(object$y[[1]], "nonlinear") & any(type %in% c("probabilities", "cumhaz"))){
    stop("Computation of cumulative hazard and survival probabilities are currently 
          under construction for associations which are nonlinear in mu.")
  }
  if(missing(dt)) dt <- 0
  if(missing(steps)) steps <- 1
  if(missing(id)) i <- NULL else i <- id
  idvar <- attr(object$y[[1]], "idvar")
  
  if(length(type) > 1)
    type <- type[1]
  type <- match.arg(type)
  loglik <- type == "loglik"
  if(loglik) {
    type <- "cumhaz"
    dt <- 0
  }
  
  if(type == "probabilities"){
    if(!is.null(newdata)){
      warning("The provided newdata will be ignored for the prediction of conditional survival probabilities.")
      newdata <- NULL
    }
    if(dt == 0){
      stop("Please specify a time window for the prediction of conditional survival probabilities.")
    }
  }
  
  if(type == "cumhaz" & !is.null(newdata) & dt > 0){
    warning("The provided newdata will be ignored for the prediction of conditional survival t + dt.")
    newdata <- NULL
  }
  
  calculate_mu <- FALSE
  if(attr(object$y[[1]], "nonlinear") &                            # only if nonlinear effect
     (("alpha" %in% list(...)$model) | is.null(list(...)$model)) & # only if alpha amongst predicted effects
     (is.null(newdata)|is.null(newdata$mu))                        # only if mu is not explicitely specified in newdata
  ){
    calculate_mu <- TRUE
  }
  
  if(is.null(newdata)){
    newdata <- model.frame(object) 
  }
  
  if(calculate_mu){
    object$family$predict <- NULL
    long_timevar <- attr(object$y[[1]], "timevar")["mu"]
    surv_timevar <- attr(object$y[[1]], "timevar")["lambda"]
    tempdata <- newdata
    tempdata[, long_timevar] <- tempdata[, surv_timevar]
    mu <- predict.bamlss(object, model = "mu", newdata = tempdata, type = type,
                         cores = cores, chunks = chunks, verbose = verbose)
    newdata$mu <- mu
  }
  
  if(!(type %in% c("probabilities", "cumhaz"))) {
    object$family$predict <- NULL
    return(predict.bamlss(object, newdata = newdata, type = type,
                          FUN = FUN, cores = cores, chunks = chunks, verbose = verbose, ...))
  }
  
  if(object$family$family != "jm")
    stop("object must be a joint-model!")
  
    dalpha <- has_pterms(object$x$dalpha$terms) | (length(object$x$dalpha$smooth.construct) > 0)
  
  timevar <- attr(object$y[[1]], "timevar")
  tmax_model <- max(newdata[,timevar["lambda"]])
  
  if(!is.null(i)){
    if(!is.character(i))
      i <- levels(newdata[[idvar]])[i]
    newdata <- subset(newdata, newdata[[idvar]] %in% i)
  }
  
  tmax_pred <- max(newdata[,timevar["mu"]]) + dt
  
  if(tmax_pred > tmax_model){
    warning("Predictions should not be made beyond the modelled time range. 
            Please adjust the time window dt accordingly.")
  }
  
  
  
  ## Create the time grid.  
  grid <- function(lower, upper, length) {
    seq(from = lower, to = upper, length = length)
  }
  
  jm_probs <- function(data) {
    
    if(dt == 0){
      take <- !duplicated(data[, c(timevar["lambda"], idvar)])
      dsurv <- subset(data, take)
      timegrid <- lapply(dsurv[[timevar["lambda"]]], function(x){grid(0, x, subdivisions)})
    } else {
      take <- !duplicated(data[, c(timevar["lambda"], idvar)], fromLast = TRUE)
      dsurv <- subset(data, take)
      timegrid <- lapply(dsurv[[timevar["mu"]]], function(x){grid(x, x+dt, subdivisions)})
    }
    nobs <- nrow(dsurv)
    gdim <- c(length(timegrid), length(timegrid[[1]]))
    width <- rep(NA, nobs)
    for(i in 1:nobs)
      width[i] <- timegrid[[i]][2] - timegrid[[i]][1]
    
    pred.setup <- predict.bamlss(object, data, type = "link",
                                 get.bamlss.predict.setup = TRUE, ...)

    enames <- pred.setup$enames
    
    pred_gamma <- with(pred.setup, .predict.bamlss("gamma",
                                                   object$x$gamma, samps, enames$gamma, intercept,
                                                   nsamps, dsurv))
    
    pred_lambda <- with(pred.setup, .predict.bamlss.surv.td("lambda",
                                                            object$x$lambda$smooth.construct, samps, enames$lambda, intercept,
                                                            nsamps, dsurv, timevar["lambda"], timegrid,
                                                            drop.terms.bamlss(object$x$lambda$terms, sterms = FALSE, keep.response = FALSE),
                                                            type = 2))
    #!# Todo: Fix prediction for nonlinear
    pred_mu <- with(pred.setup, .predict.bamlss.surv.td("mu",
                                                        object$x$mu$smooth.construct, samps, enames$mu, intercept,
                                                        nsamps, dsurv, timevar["mu"], timegrid,
                                                        drop.terms.bamlss(object$x$mu$terms, sterms = FALSE, keep.response = FALSE)))
    
    pred_alpha <- with(pred.setup, .predict.bamlss.surv.td("alpha",
                                                           object$x$alpha$smooth.construct, samps, enames$alpha, intercept,
                                                           nsamps, dsurv, timevar["lambda"], timegrid,
                                                           drop.terms.bamlss(object$x$alpha$terms, sterms = FALSE, keep.response = FALSE)))
    
    if(dalpha) {
      pred_dalpha <- with(pred.setup, .predict.bamlss.surv.td("dalpha",
                                                              object$x$dalpha$smooth.construct, samps, enames$dalpha, intercept,
                                                              nsamps, dsurv, timevar["lambda"], timegrid,
                                                              drop.terms.bamlss(object$x$dalpha$terms, sterms = FALSE, keep.response = FALSE)))
      
      pred_dmu <- with(pred.setup, .predict.bamlss.surv.td("dmu",
                                                           object$x$dmu$smooth.construct, samps, enames$dmu, intercept,
                                                           nsamps, dsurv, timevar["mu"], timegrid,
                                                           drop.terms.bamlss(object$x$dmu$terms, sterms = FALSE, keep.response = FALSE),
                                                           derivMat = TRUE))
    }
    if(attr(object$y[[1]], "nonlinear")){
      eta_timegrid <- if(dalpha) {
        pred_lambda + pred_alpha + pred_dalpha * pred_dmu
      } else {
        pred_lambda + pred_alpha
      }
    } else {
    eta_timegrid <- if(dalpha) {
      pred_lambda + pred_alpha * pred_mu + pred_dalpha * pred_dmu
    } else {
      pred_lambda + pred_alpha * pred_mu
      }
    }

    if(loglik) {
      eta_gamma <- predict.bamlss(object, data[take, , drop = FALSE], model = "gamma")
      eta_mu <- predict.bamlss(object, data, model = "mu")
      eta_sigma <- predict.bamlss(object, data, model = "sigma")
      mf <- model.frame(object, data = data)
      y <- mf[, grep("Surv", names(mf), fixed = TRUE)]
      eta_timegrid <- matrix(eta_timegrid, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
      eeta <- exp(eta_timegrid)
      int <- width * (0.5 * (eeta[, 1] + eeta[, subdivisions]) + apply(eeta[, 2:(subdivisions - 1)], 1, sum))
      ll <- sum((eta_timegrid[, ncol(eta_timegrid)] + eta_gamma) * y[take, "status"], na.rm = TRUE) -
        exp(eta_gamma) %*% int + sum(dnorm(y[, "obs"], mean = eta_mu, sd = exp(eta_sigma), log = TRUE))
      return(drop(ll))
    }
    
    if(dt == 0){
      probs <- NULL
      for(i in 1:ncol(eta_timegrid)) {
        eta <- matrix(eta_timegrid[, i], nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
        eeta <- exp(eta)
        int <- width * (0.5 * (eeta[, 1] + eeta[, subdivisions]) + apply(eeta[, 2:(subdivisions - 1), drop = FALSE], 1, sum))
        probs <- if(type == "probabilities") {
          cbind(probs, exp(-1 * exp(pred_gamma[, i]) * int))
        } else {
          cbind(probs, exp(pred_gamma[, i]) * int)
        }
      }
      
      if(!is.null(FUN)) {
        if(is.matrix(probs)) {
          if(ncol(probs) > 1)
            probs <- apply(probs, 1, FUN)
          probs <- t(probs) 
        } 
      }
      
    } else {
      lprobs <- lapply(1:steps, function(x){
        probs <- NULL
        sub <- round(subdivisions/steps * x)
        for(i in 1:ncol(eta_timegrid)) {
          eta <- matrix(eta_timegrid[, i], nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
          eeta <- exp(eta)
          int <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1), drop = FALSE], 1, sum))
          probs <- if(type == "probabilities") {
            cbind(probs, exp(-1 * exp(pred_gamma[, i]) * int))
          } else {
            cbind(probs, exp(pred_gamma[, i]) * int)
          }
        }
        
        if(!is.null(FUN)) {
          if(is.matrix(probs)) {
            if(ncol(probs) > 1)
              probs <- apply(probs, 1, FUN)
            probs <- t(probs) 
          } 
        }
        probs
      }) 
      
      probs <- lprobs
      names(probs) <- paste("Time after last longitudinal observation:", (1:steps)*dt/steps) 
    }
    
    
    return(probs)
  }
  
  ia <- interactive()
  
  if(is.null(cores)) {
    if(chunks < 2) {
      probs <- jm_probs(newdata)
    } else {
      id <- sort(rep(1:chunks, length.out = nrow(newdata)))
      newdata <- split(newdata, id)
      chunks <- length(newdata)
      probs <- NULL
      for(i in 1:chunks) {
        if(verbose) {
          cat(if(ia) "\r" else "\n")
          cat("predicting chunk", i, "of", chunks, "...")
          if(.Platform$OS.type != "unix" & ia) flush.console()
        }
        if(i < 2) {
          probs <- jm_probs(newdata[[i]])
        } else {
          if(is.null(dim(probs))) {
            probs <- c(probs, jm_probs(newdata[[i]]))
          } else {
            probs <- rbind(probs, jm_probs(newdata[[i]]))
          }
        }
      }
      if(verbose) cat("\n")
    }
  } else {
    parallel_fun <- function(i) {
      if(chunks < 2) {
        pr <- jm_probs(newdata[[i]])
      } else {
        idc <- sort(rep(1:chunks, length.out = nrow(newdata[[i]])))
        nd <- split(newdata[[i]], idc)
        chunks <- length(nd)
        pr <- NULL
        for(j in 1:chunks) {
          if(j < 2) {
            pr <- jm_probs(nd[[j]])
          } else {
            if(is.null(dim(pr))) {
              pr <- c(pr, jm_probs(nd[[j]]))
            } else {
              pr <- rbind(pr, jm_probs(nd[[j]]))
            }
          }
        }
      }
      return(pr)
    }
    
    id <- sort(rep(1:cores, length.out = nrow(newdata)))
    newdata <- split(newdata, id)
    cores <- length(newdata)
    probs <- parallel::mclapply(1:cores, parallel_fun, mc.cores = cores)
    
    probs <- if(is.matrix(probs[[1]])) {
      do.call("rbind", probs)
    } else {
      do.call("c", probs)
    }
  }
  
  return(probs)
}



jm_survplot <- function(object, id = 1, dt = NULL, steps = 10, 
                        points = TRUE, rug = !points)
{
  on.exit(par(par(no.readonly = TRUE)))
  
  FUN <- function(x) { mean(x, na.rm = TRUE) }
  idvar <- attr(object$y[[1]], "idvar")
  timevar <- attr(object$y[[1]], "timevar")
  mf <- model.frame(object)
  i <- id
  
  if(!is.character(i))
    i <- levels(mf[[idvar]])[i]
  ii <- mf[[idvar]] == i
  
  tmax <- max(mf[[timevar["mu"]]][ii])
  tmax_model <- max(mf[[timevar["lambda"]]])
  if(is.null(dt))
    dt <- 0.4 * tmax
  maxtime <- tmax + dt
  if(maxtime > tmax_model){
    maxtime <- tmax_model
    dt <- maxtime - tmax
  }
    
  time <- seq(tmax, maxtime, length.out = steps)
  
  if(all(time == time[1])) {
    stop(paste("Not enough time points available for individual ",
               i, "!", sep = ""))
  }
  
  p_surv <- predict(object, type = "probabilities", dt = dt, steps = steps-1, id = id, FUN = c95)
  p_surv <- do.call(rbind, p_surv)
  # fix the last observed longitudinal timepoint to p_surv = 1
  p_surv <- rbind(c(1,1,1), p_surv)
  
  time2 <- sort(c(tmax, seq(0, maxtime, length = steps - 1)))
  nd <- data.frame(time2, time2)
  names(nd) <- timevar
  vars <- names(mf)[!(names(mf) %in% c(timevar, idvar))]
  for(j in vars) {
    if(!is.factor(mf[[j]])) {
      nd[[j]] <- FUN(mf[[j]][ii])
    } else {
      ftab <- table(mf[[j]][ii])
      ntab <- names(ftab)
      nd[[j]] <- factor(ntab[which.max(as.vector(ftab))], levels = levels(mf[[j]]))
    }
  }
  nd[[idvar]] <- factor(i, levels = levels(mf[[idvar]]))
  
  p_long <- predict(object, newdata = nd, model = "mu", FUN = c95)
  
  s2.col <- rev(rgb(0, 0, 1, alpha = seq(0.1, 0.01, length = 50)))
  
  par(mfrow = c(2, 1), mar = rep(0, 4),
      oma = c(4.1, 4.1, 1.1, 4.1))
  plot2d(p_surv ~ time, fill.select = c(0, 1, 0, 1),
         scheme = 2, axes = FALSE, ylim = c(0, 1),
         xlim = c(0, maxtime), s2.col = s2.col, xlab = "", ylab = "")
  abline(v = tmax, lty = 2)
  axis(2)
  box()
  mtext("Prob(T > t + dt |T > t)", side = 2, line = 2.5)
  plot2d(p_long[time2 <= tmax, ] ~ time2[time2 <= tmax], fill.select = c(0, 1, 0, 1),
         scheme = 2, axes = FALSE, xlim = c(0, maxtime),
         xlab = "", ylab = "",
         ylim = if(points) range(c(p_long, object$y[[1]][ii, "obs"])) else range(p_long))
  plot2d(p_long[time2 >= tmax, ] ~ time2[time2 >= tmax], fill.select = c(0, 1, 0, 1),
         scheme = 2, axes = FALSE, add = TRUE, s2.col = s2.col,
         xlab = "", ylab = "")
  abline(v = tmax, lty = 2)
  if(points)
    points(mf[[timevar["mu"]]][ii], object$y[[1]][ii, "obs"])
  if(rug)
    rug(mf[[timevar["mu"]]][ii])
  axis(1)
  axis(4)
  box()
  mtext("Effect of time", side = 4, line = 2.5)
  mtext("Time", side = 1, line = 2.5)
  
  invisible(NULL)
}



.predict.bamlss.jm.td <- function(id, x, samps, enames, intercept, nsamps, newdata, env,
                                  yname, grid, formula, type = 1, derivMat = FALSE)
{
  snames <- colnames(samps)
  enames <- gsub("p.Intercept", "p.(Intercept)", enames, fixed = TRUE)
  has_intercept <- any(grepl(paste(id, "p", "(Intercept)", sep = "."), snames, fixed = TRUE))
  if(!has_intercept)
    has_intercept <- any(grepl(paste(id, "p.model.matrix", "(Intercept)", sep = "."), snames, fixed = TRUE))
  if(intercept & has_intercept)
    enames <- c("p.(Intercept)", enames)
  enames <- unique(enames)
  ec <- sapply(enames, function(x) {
    paste(strsplit(x, "")[[1]][1:2], collapse = "")
  })
  enames2 <- sapply(enames, function(x) {
    paste(strsplit(x, "")[[1]][-c(1:2)], collapse = "")
  })
  
  eta <- 0
  if(length(i <- grep("p.", ec))) {
    for(j in enames2[i]) {
      if(j != "(Intercept)") {
        f <- as.formula(paste("~", if(has_intercept) "1" else "-1", "+", j), env = env)
        X <- param_Xtimegrid(f, newdata, grid, yname, type = type, derivMat = derivMat)
        if(has_intercept)
          X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
        sn <- snames[grep2(paste(id, "p", j, sep = "."), snames, fixed = TRUE)]
        if(!length(sn))
          sn <- snames[grep2(paste(id, "p.model.matrix", j, sep = "."), snames, fixed = TRUE)]
        eta <- if(is.null(eta)) {
          fitted_matrix(X, samps[, sn, drop = FALSE])
        } else {
          eta + fitted_matrix(X, samps[, sn, drop = FALSE])
        }
      } else {
        if(has_intercept) {
          sn <- snames[grep2(paste(id, "p", j, sep = "."), snames, fixed = TRUE)]
          if(!length(sn))
            sn <- snames[grep2(paste(id, "p.model.matrix", j, sep = "."), snames, fixed = TRUE)]
          eta <- eta + fitted_matrix(matrix(1, nrow = length(grid[[1]]) * nrow(newdata), ncol = 1),
                                     samps[, sn, drop = FALSE])
        }
      }
    }
  }
  if(length(i <- grep("s.", ec))) {
    for(j in enames2[i]) {
      if(!inherits(x[[j]], "no.mgcv") & !inherits(x[[j]], "special")) {
        X <- sm_Xtimegrid(x[[j]], newdata, grid, yname, derivMat = derivMat)
        sn <- snames[grep2(paste(id, "s", j, sep = "."), snames, fixed = TRUE)]
        eta <- if(is.null(eta)) {
          fitted_matrix(X, samps[, sn, drop = FALSE])
        } else {
          eta + fitted_matrix(X, samps[, sn, drop = FALSE])
        }
      } else {
        stop("no predictions for special terms available yet!")
      }
    }
  }
  
  eta
}


# Adjusted function for the fitting of nonlinear effects
sm_time_transform2 <- function(x, data, grid, yname, timevar, take, derivMat = FALSE, eps = 1e-7,
                               update_data = NULL)
{
  if(derivMat)
    ddata <- data
  if(!is.null(take))
    data <- data[take, , drop = FALSE]
  X <- NULL
  terms <- x$term[x$term %in% names(data)]
  for(j in terms) {
    if((j != yname) & (j != timevar)) {
      df <- data.frame(rep(data[[j]], each = length(grid[[1]])))
      names(df) <- j
      X <- if(is.null(X)) df else cbind(X, df)
    }
  }
  if(!is.null(X))
    colnames(X) <- x$term[!(x$term %in% c(yname, timevar))]
  
  X <- if(is.null(X)) data.frame(unlist(grid)) else cbind(X, unlist(grid))
  colnames(X)[ncol(X)] <- yname
  if(timevar != yname) {
    X <- cbind(X, unlist(grid))
    colnames(X)[ncol(X)] <- timevar
  }
  if(x$by != "NA" & x$by != yname){
    X[[x$by]] <- rep(data[[x$by]], each = length(grid[[1]]))
    if(!is.null(x$by.level)) 
      x$by_timegrid <- as.numeric(X[[x$by]] == x$by.level)
  }
  
  dX <- NULL
  if(derivMat) {
    dX <- X
    for(j in colnames(dX)) {
      if(!is.factor(dX[[j]]) & (grepl(timevar, j, fixed = TRUE)) & (timevar %in% c(x$term, x$by))) {
        dX[[j]] <- dX[[j]] + eps
        ddata[[j]] <- ddata[[j]] + eps
      }
    }
  }
  
  Xkeep <- X
  
  # make update function for later updates
  x$update <- function(update_data, name){
    Xkeep[[name]] <- update_data
    X <- PredictMat(x, Xkeep)
    X
  }
  x$propose <- NULL
  
  X <- PredictMat(x, X)
  gdim <- c(length(grid), length(grid[[1]]))
  
  x$XT <- extract_XT(X, gdim[1], gdim[2])
  
  if(derivMat) {
    dX <- PredictMat(x, dX)
    X <- -1 * (X - dX) / eps
    x$XT <- extract_XT(dX, gdim[1], gdim[2])
    x$X <- -1 * (x$X - PredictMat(x, ddata)) / eps
    gm <- sparse.matrix.index(X)
    if(!is.null(gm))
      x$sparse.setup$grid.matrix <- gm
  } else {
    x$sparse.setup$grid.matrix <- sparse.matrix.index(X)
  }
  ff <- make.fit.fun(x, type = 2)
  
  x$all_zero <- all(X == 0) & all(x$X == 0)
  x$timevar <- timevar
  
  x$fit.fun_timegrid <- function(g, Xupdate = NULL) {
    if(!is.null(Xupdate)) X <- Xupdate 
    if(is.null(g)) return(X)
    if(x$all_zero) return(rep(0, gdim[1]))
    f <- ff(X, g, expand = FALSE)
    f <- matrix(f, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    f
  }
  
  x$state$fitted_timegrid <- x$fit.fun_timegrid(get.state(x, "b"))
  x$state$optimize <- FALSE
  
  if(derivMat) {
    x$xt$orig.class <- class(x)
    x$xt$deriv.eps <- eps
    class(x) <- "deriv.smooth"
  }
  
  x
}


update_jm_alpha_nonlin <- function(x, eta, eta_timegrid,
                                  eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                                  status, update.nu, width, criterion, get_LogPost, nobs, do.optim2, edf, edf_alt = FALSE, tp = FALSE, ...)
{
  ## The time-dependent design matrix for the grid.
  if("model.matrix" %in% attr(x, "class")){
    X <- x$fit.fun_timegrid(NULL)
  } else {
    Xmu <- as.vector(t(eta_timegrid_mu))
    if(!tp) {
      X <- x$update(Xmu, "mu")    
    } else {
      Xalpha <- x$margin[[1]]$update(Xmu, "mu")
      Xalpha2 <- x$margin[[2]]$fit.fun_timegrid(NULL)
      X <- rowTensorProduct(Xalpha, Xalpha2) 
    }
  }

  gdim <- c(nrow(eta_timegrid), ncol(eta_timegrid))
  x$XT <- extract_XT(X, gdim[1], gdim[2])
  
  ## Timegrid predictor.
  eeta <- exp(eta_timegrid)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta, width, exp(eta$gamma), index = NULL)
  
  ## check for tensor-product, if sparse.setup helps
  xgrad <- drop(t(x$XT) %*% status - int$grad)
  
  env <- new.env()
  
  if((!(!x$state$do.optim | x$fixed | x$fxsp)) & do.optim2) {
    par <- x$state$parameters
    
    edf0 <- if(is.null(edf)) {
      0
    } else {
      edf - x$state$edf
    }
    
    objfun1 <- function(tau2) {
      par <- set.par(par, tau2, "tau2")
      xgrad <- xgrad + x$grad(score = NULL, par, full = FALSE)
      xhess <- int$hess + x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(xhess, index = NULL)
      Hs <- Sigma %*% xgrad
      g <- get.par(par, "b")
      if(update.nu) {
        objfun.nu <- function(nu) {
          g2 <- drop(g + nu * Hs)
          names(g2) <- names(g)
          x$state$parameters <- set.par(x$state$parameters, g2, "b")
          x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
          
          ## Update additive predictors
          # fit_timegrid <- x$fit.fun_timegrid(g2, X)
          fit_timegrid <- matrix(drop(X %*% g2), nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
          eta_timegrid_alpha <- eta_timegrid_alpha - x$state$fitted_timegrid + fit_timegrid
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
          lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
          if(lp == -Inf) lp <- -1 * .Machine$double.xmax
          return(-1 * lp)
        }
        nu <- optimize(f = objfun.nu, interval = c(0, 1))$minimum
      } else {
        nu <- x$state$nu
      }
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)

      fit <- x$fit.fun(x$XT, g2, no.sparse.setup = TRUE) 
      #!# check if still problematic: not equivalent to X%*% g2 due to sparse index in C-function
      # fit_timegrid <- x$fit.fun_timegrid(g2, X)
      fit_timegrid <- matrix(drop(X %*% g2), nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
      eta$alpha <- eta$alpha - fitted(x$state) + fit
      eta_timegrid_alpha <- eta_timegrid_alpha - x$state$fitted_timegrid + fit_timegrid
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
      S <- int$hess %*% Sigma
      if(edf_alt){
        edf1 <- sum_diag(2*S - S%*%S)
      } else {
        edf1 <- sum_diag(S)
      }
      #checks# if(edf1 < 0) browser()
      edf <- edf0 + edf1
      logLik <- get_LogPost(eta_timegrid, eta, 0)
      ic <- get.ic2(logLik, edf, length(eta$mu), criterion)
      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          opt_state <- list("parameters" = set.par(par, g2, "b"),
                            "fitted.values" = fit, "fitted_timegrid" = fit_timegrid,
                            "edf" = edf1, "hessian" = xhess,
                            "nu" = nu, "do.optim" = x$state$do.optim)
          assign("state", opt_state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }
    
    assign("ic00_val", objfun1(get.state(x, "tau2")), envir = env)
    tau2 <- tau2.optim(objfun1, start = get.state(x, "tau2"))
    # cat("\n objfun --- ", x$label, ", tau: ", tau2, ", edf: ", edf)
    if(!is.null(env$state))
      return(env$state)
    
    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
  }

  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)

  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(xhess, x$sparse.setup)
  Hs <- Sigma %*% xgrad

  ## Update regression coefficients.
  g <- get.state(x, "b")
  
  if(update.nu) {
    objfun2 <- function(nu) {
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      x$state$parameters <- set.par(x$state$parameters, g2, "b")
      
      ## Update additive predictors.
      # fit_timegrid <- x$fit.fun_timegrid(g2, X)
      fit_timegrid <- matrix(drop(X %*% g2), nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
      eta_timegrid_alpha <- eta_timegrid_alpha - x$state$fitted_timegrid + fit_timegrid
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
      lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
      if(lp == -Inf) lp <- -1 * .Machine$double.xmax
      return(-1 * lp)
    }
    x$state$nu <- optimize(f = objfun2, interval = c(0, 1))$minimum
  }

  g2 <- drop(g + x$state$nu * Hs)        
  names(g2) <- names(g)
  x$state$parameters <- set.par(x$state$parameters, g2, "b")
  ## Update fitted values.
  # x$state$fitted_timegrid <- x$fit.fun_timegrid(g2, X)
  x$state$fitted_timegrid <- matrix(drop(X %*% g2), nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
  x$state$fitted.values <- x$state$fitted_timegrid[, gdim[2]]
  # x$state$fitted.values <- x$fit.fun(x$XT, g2, expand = FALSE, no.sparse.setup = TRUE)
  x$state$edf <- sum_diag(int$hess %*% Sigma)
  x$state$hessian <- xhess
  
  return(x$state)
}

propose_jm_alpha_nonlin <- function(x, y,
                                   eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                                   eta_timegrid_dmu, eta_timegrid_dalpha,
                                   width, sub, nu, status, id, tp = FALSE, fac = FALSE, ...)
{
  ## The time-dependent design matrix for the grid.
  if("model.matrix" %in% attr(x, "class")){
    X <- x$fit.fun_timegrid(NULL)
  } else {
    Xmu <- as.vector(t(eta_timegrid_mu))
    if(!tp) {
      X <- x$update(Xmu, "mu")    
    } else {
      Xalpha <- x$margin[[1]]$update(Xmu, "mu")
      Xalpha2 <- x$margin[[2]]$fit.fun_timegrid(NULL)
      X <- rowTensorProduct(Xalpha, Xalpha2) 
    }
  }
  
  gdim <- c(nrow(eta_timegrid), ncol(eta_timegrid))
  x$XT <- extract_XT(X, gdim[1], gdim[2])
  
  ## Timegrid lambda.
  eeta <- exp(eta_timegrid)
  
  ## Old logLik and prior.
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibeta <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  p1 <- x$prior(x$state$parameters)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta, width, exp(eta$gamma), index = NULL)
  intold <- int
  xgrad <- drop(t(x$XT) %*% status - int$grad)
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(xhess, index = NULL)
  
  ## Save old coefficients.
  g0 <- get.state(x, "b")
  
  ## Get new position.
  mu <- drop(g0 + nu * Sigma %*% xgrad)
  
  ## Sample new parameters.
  g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma, method="chol"))
  names(g) <- names(g0)
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  
  ## Compute log priors.
  p2 <- x$prior(x$state$parameters)
  qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)
  
  ## Update additive predictors.
  fit_timegrid <- matrix(drop(X %*% g), nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
  # fit_timegrid <- x$fit.fun_timegrid(g)
  eta_timegrid_alpha <- eta_timegrid_alpha - x$state$fitted_timegrid + fit_timegrid
  x$state$fitted_timegrid <- fit_timegrid
  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
  
  fit <- fit_timegrid[, gdim[2]]
  eta$alpha <- eta$alpha - fitted(x$state) + fit
  x$state$fitted.values <- fit
  
  ## New logLik.
  eeta <- exp(eta_timegrid)
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibetaprop <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  
  ## Prior prob.
  int <- survint(X, eeta, width, exp(eta$gamma), index = NULL)
  xgrad <- drop(t(x$XT) %*% status - int$grad)
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- int$hess + x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  Sigma2 <- matrix_inv(xhess, index = NULL)
  if(inherits(Sigma2, "try-error")) {
    x$state$alpha <- 1
    return(x$state)
  }
  mu2 <- drop(g + nu * Sigma2 %*% xgrad)
  qbeta <- dmvnorm(g0, mean = mu2, sigma = Sigma2, log = TRUE)
  
  ## Save edf.
  x$state$edf <- sum_diag(int$hess %*% Sigma2)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if(length(x$S) < 2) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(g, x$S[[1]]) %*% g + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "alpha", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  return(x$state)
}


modSplineDesign <- function(knots, x, ord = 4, derivs = 0, sparse = FALSE, outer.ok = TRUE){
  # modified from mgcv:::Predict.matrix.pspline.smooth
  # ensures a more moderate behaviour if values are outside of the inner knots
  m <- ord - 1
  ll <- knots[m + 1]
  ul <- knots[length(knots) - m]
  m <- m + 1
  n <- length(x)
  ind <- x <= ul & x >= ll
  if (sum(ind) == n) {
    X <- splines::spline.des(knots, x, m, rep(derivs, n), sparse = sparse)$design
  }
  else {
    D <- splines::spline.des(knots, c(ll, ll, ul, ul), m, c(0, 1, 0, 1), sparse = sparse)$design
    X <- matrix(0, n, ncol(D))
    nin <- sum(ind)
    if (nin > 0) 
      X[ind, ] <- splines::spline.des(knots, x[ind], m, rep(derivs, nin), sparse = sparse)$design
    if (derivs < 2) {
      ind <- x < ll
      if (sum(ind) > 0) 
        X[ind, ] <- if (derivs == 0) 
          cbind(1, x[ind] - ll) %*% D[1:2, ]
      else matrix(D[2, ], sum(ind), ncol(D), byrow = TRUE)
      ind <- x > ul
      if (sum(ind) > 0) 
        X[ind, ] <- if (derivs == 0) 
          cbind(1, x[ind] - ul) %*% D[3:4, ]
      else matrix(D[4, ], sum(ind), ncol(D), byrow = TRUE)
    }
  }
  X
}

constrain <- function(sm, X) {
  # apply smooth.construct constrain to new design matrix
  # adapted from mgcv::smoothCon
  ncolC <- ncol(X)
  j <- attr(sm, "nCons")
  qrc <- attr(sm, "qrc")
  indi <- attr(sm, "indi")
  
  if (j > 0) {
    nx <- length(indi)
    if (nx < ncolC & !is.null(indi)) {
      nc <- j
      nz <- nx - nc
      X[, indi[1:nz]] <- t(qr.qty(qrc, t(X[, indi, drop = FALSE]))[(nc + 1):nx, ])
      X <- X[, -indi[(nz + 1):nx]]
    }
    else {
      {
        k <- ncol(X)
        X <- t(qr.qty(qrc, t(X))[(j + 1):k, ])
      }
    }
  }
  X
}

update_jm_mu_nonlin <- function(x, y, eta, eta_timegrid,
                               eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                               status, update.nu, width, criterion, get_LogPost, nobs, do.optim2, edf,
                               dx = NULL, xsmalpha = NULL, knots = NULL, tp = FALSE, fac = FALSE, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  dX <- if(!is.null(dx)) dx$fit.fun_timegrid(NULL) else NULL
  gdim <- c(nrow(eta_timegrid), ncol(eta_timegrid))

  ## Timegrid predictor.
  eeta <- exp(eta_timegrid)
  
  eta_timegrid_alpha1 <- eta_timegrid_alpha2 <- 0 
  ## evaluations of alpha
  for (i in names(xsmalpha)){
    if(i != "model.matrix"){      # derivative for model.matrix is  0
      g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
      dXalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 1) 
      dXalpha <- constrain(xsmalpha[[i]], dXalpha)
      d2Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 2)
      d2Xalpha <- constrain(xsmalpha[[i]], d2Xalpha)
      
      if(fac){
        dXalpha <- dXalpha * xsmalpha[[i]]$by_timegrid
        d2Xalpha <- d2Xalpha * xsmalpha[[i]]$by_timegrid
      }
      
      if(tp){
        Xalpha2 <- xsmalpha[[1]]$margin[[2]]$fit.fun_timegrid(NULL)
        dXalpha <- rowTensorProduct(dXalpha, Xalpha2) 
        d2Xalpha <- rowTensorProduct(d2Xalpha, Xalpha2) 
      }
      eta_timegrid_alpha1 <- eta_timegrid_alpha1 + matrix(dXalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE) 
      eta_timegrid_alpha2 <- eta_timegrid_alpha2 + matrix(d2Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    }
 }
  
  eta_alpha1 <- eta_timegrid_alpha1[, gdim[2]]
  eta_alpha2 <- eta_timegrid_alpha2[, gdim[2]]

  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_alpha1, width, exp(eta$gamma),
                 eeta * ((eta_timegrid_alpha1)^2 + eta_timegrid_alpha2),
                 index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  xgrad <- drop(t(x$X) %*% drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2) +
                  t(x$XT) %*% drop(eta_alpha1 * status) - int$grad)
  # if(!is.null(dx))
  #   xgrad <- drop(xgrad + t(dx$XT) %*% drop(eta$dalpha * status))
  XWX <- if(is.null(x$sparse.setup$matrix)) {
    crossprod(x$X * (1 / exp(eta$sigma)^2), x$X)
  } else do.XWX(x$X, exp(eta$sigma)^2, index = x$sparse.setup$matrix)
  XAX <- crossprod(x$XT * drop(eta_alpha2 * status), x$XT)
  xhess0 <- -1 * XWX + XAX - int$hess

    env <- new.env()
  
  if((!(!x$state$do.optim | x$fixed | x$fxsp)) & do.optim2) {
    par <- x$state$parameters
    
    edf0 <- if(is.null(edf)) {
      0
    } else {
      edf - x$state$edf
    }
    objfun1 <- function(tau2) {
      par <- set.par(par, tau2, "tau2")
      xgrad <- xgrad + x$grad(score = NULL, par, full = FALSE)
      xhess <- xhess0 - x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(-1 * xhess, index = NULL)
      Hs <- Sigma %*% xgrad
      g <- get.par(par, "b")
      if(update.nu) {
        objfun.nu <- function(nu) {
          g2 <- drop(g + nu * Hs)
          names(g2) <- names(g)
          x$state$parameters[grepl("b", names(x$state$parameters))] <- g2
          x$state$parameters[grepl("tau2", names(x$state$parameters))] <- tau2
          fit_timegrid <-x$fit.fun_timegrid(g2)
          
          eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
          eta_timegrid_alpha <- 0 
          for (i in names(xsmalpha)){
            g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
            if(i != "model.matrix"){
              Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
              Xalpha <- constrain(xsmalpha[[i]], Xalpha)
              if(fac)
                Xalpha <- Xalpha * xsmalpha[[i]]$by_timegrid
              if(tp)  
                Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
              
              eta_timegrid_alpha <- eta_timegrid_alpha + matrix(Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
            } else {
              eta_timegrid_alpha <- eta_timegrid_alpha + xsmalpha[[i]]$state$fitted_timegrid 
            }
          }

          # if(!is.null(dx))
          #   eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
          fit <- x$fit.fun(x$X, g2)
          eta$mu <- eta$mu - fitted(x$state) + fit
          lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
          if(lp == -Inf) lp <- -1 * .Machine$double.xmax
          return(-1 * lp)
        }
        nu <- optimize(f = objfun.nu, interval = c(0, 1))$minimum
      } else {
        nu <- x$state$nu
      }
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      fit <- x$fit.fun(x$X, g2)
      fit_timegrid <- x$fit.fun_timegrid(g2)
      
      eta$mu <- eta$mu - fitted(x$state) + fit 
      eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid

      eta_timegrid_alpha <- 0 
      for (i in names(xsmalpha)){
        g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
        if(i != "model.matrix"){
          Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
          Xalpha <- constrain(xsmalpha[[i]], Xalpha)
          if(fac)
            Xalpha <- Xalpha * xsmalpha[[i]]$by_timegrid
          if(tp)  
            Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
          
          eta_timegrid_alpha <- eta_timegrid_alpha + matrix(Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
        } else {
          eta_timegrid_alpha <- eta_timegrid_alpha + xsmalpha[[i]]$state$fitted_timegrid 
        }
      }

      # if(!is.null(dx))
      #   eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
      edf1 <- sum_diag((-1 * xhess0) %*% Sigma)
      # if(edf1 < 0) browser()
      edf <- edf0 + edf1
      logLik <- get_LogPost(eta_timegrid, eta, 0)
      ic <- get.ic2(logLik, edf, length(eta$mu), criterion)

      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          opt_state <- list("parameters" = set.par(par, g2, "b"),
                            "fitted.values" = fit, "fitted_timegrid" = fit_timegrid,
                            "edf" = edf1, "hessian" = -1 * xhess,
                            "nu" = nu, "do.optim" = x$state$do.optim)
          assign("state", opt_state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }
    assign("ic00_val", objfun1(get.state(x, "tau2")), envir = env)
    tau2 <- tau2.optim(objfun1, start = get.state(x, "tau2"), maxit = 1)

    if(!is.null(env$state))
      return(env$state)
    
    x$state$parameters[grepl("tau2", names(x$state$parameters))] <- tau2
  }
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- xhess0 - x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(-1 * xhess, index = NULL)
  Hs <- Sigma %*% xgrad
  
  ## Update regression coefficients.
  g <- get.state(x, "b")
  
  if(update.nu) {
    objfun2 <- function(nu) {
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      x$state$parameters <- set.par(x$state$parameters, g2, "b")
      fit_timegrid <- matrix(X %*% g2, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
      
      eta_timegrid_alpha <- 0 
      for (i in names(xsmalpha)){
        g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
        if(i != "model.matrix"){
          Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
          Xalpha <- constrain(xsmalpha[[i]], Xalpha)
          if(fac)
            Xalpha <- Xalpha * xsmalpha[[i]]$by_timegrid
          if(tp)  
            Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
          
          eta_timegrid_alpha <- eta_timegrid_alpha + matrix(Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
        } else {
          eta_timegrid_alpha <- eta_timegrid_alpha + xsmalpha[[i]]$state$fitted_timegrid 
        }
      }

      # if(!is.null(dx))
      #   eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
      fit <- x$fit.fun(x$X, g2)
      eta$mu <- eta$mu - fitted(x$state) + fit
      lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
      if(lp == -Inf) lp <- -1 * .Machine$double.xmax
      return(-1 * lp)
    }
    x$state$nu <- optimize(f = objfun2, interval = c(0, 1))$minimum
  }
  
  g2 <- drop(g + x$state$nu * Hs)
  names(g2) <- names(g)

  x$state$parameters <- set.par(x$state$parameters, g2, "b")
  
  ## Update fitted values
  x$state$fitted_timegrid <- x$fit.fun_timegrid(g2)
  x$state$fitted.values <- x$fit.fun(x$X, g2)
  x$state$edf <- sum_diag((-1 * xhess0) %*% Sigma)
  x$state$hessian <- -1 * xhess
  
  return(x$state)
}

update_jm_mu_nonlin_Matrix <- function(x, y, eta, eta_timegrid,
                                      eta_timegrid_lambda, eta_timegrid_alpha, eta_timegrid_mu, eta_timegrid_dalpha, eta_timegrid_dmu,
                                      status, update.nu, width, criterion, get_LogPost, nobs, do.optim2, edf, 
                                      dx = NULL, xsmalpha = NULL, knots = NULL, tp = FALSE, fac = FALSE, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  dX <- if(!is.null(dx)) dx$fit.fun_timegrid(NULL) else NULL
  gdim <- c(nrow(eta_timegrid), ncol(eta_timegrid))
  if(tp) Xalpha2 <- xsmalpha$margin[[2]]$fit.fun_timegrid(NULL) 
  
  ## Timegrid predictor.
  eeta <- exp(eta_timegrid)
  
  ## evaluations of alpha
  eta_timegrid_alpha1 <- eta_timegrid_alpha2 <- 0 
  ## evaluations of alpha
  for (i in names(xsmalpha)){
    if(i != "model.matrix"){      # derivative for model.matrix is  0
      g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
      dXalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 1) 
      dXalpha <- constrain(xsmalpha[[i]], dXalpha)
      d2Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 2)
      d2Xalpha <- constrain(xsmalpha[[i]], d2Xalpha)
      
      if(fac){
        dXalpha <- dXalpha * xsmalpha[[i]]$by_timegrid
        d2Xalpha <- d2Xalpha * xsmalpha[[i]]$by_timegrid
      }
      
      if(tp){
        Xalpha2 <- xsmalpha[[1]]$margin[[2]]$fit.fun_timegrid(NULL)
        dXalpha <- rowTensorProduct(dXalpha, Xalpha2) 
        d2Xalpha <- rowTensorProduct(d2Xalpha, Xalpha2) 
      }
      eta_timegrid_alpha1 <- eta_timegrid_alpha1 + matrix(dXalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE) 
      eta_timegrid_alpha2 <- eta_timegrid_alpha2 + matrix(d2Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    }
  }
  
  eta_alpha1 <- eta_timegrid_alpha1[, gdim[2]]
  eta_alpha2 <- eta_timegrid_alpha2[, gdim[2]]
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_alpha1, width, exp(eta$gamma),
                 eeta * ((eta_timegrid_alpha1)^2 + eta_timegrid_alpha2),
                 index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  
  xgrad <- crossprod(x$X, drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2)) +
    crossprod(x$XT,  drop(eta_alpha1 * status)) - int$grad
  # if(!is.null(dx))
  #   xgrad <- xgrad + crossprod(dx$XT, drop(eta$dalpha * status))
  XWX <- crossprod(Diagonal(x = 1 / exp(eta$sigma)^2) %*% x$X, x$X)
  XAX <- crossprod(x$XT * drop(eta_alpha2 * status), x$XT)
  xhess0 <- -1 * XWX + XAX - int$hess

  env <- new.env()
  
  if((!(!x$state$do.optim | x$fixed | x$fxsp)) & do.optim2) {
    par <- x$state$parameters
    
    edf0 <- if(is.null(edf)) {
      0
    } else {
      edf - x$state$edf
    }
    objfun1 <- function(tau2) {
      par <- set.par(par, tau2, "tau2")
      xgrad <- xgrad + x$grad(score = NULL, par, full = FALSE)
      xhess <- xhess0 - x$hess(score = NULL, par, full = FALSE)
      Sigma <- matrix_inv(-1 * xhess, index = NULL)
      Hs <- Sigma %*% xgrad
      g <- get.par(par, "b")
      if(update.nu) {
        objfun.nu <- function(nu) {
          g2 <- drop(g + nu * Hs)
          names(g2) <- names(g)
          x$state$parameters <- set.par(x$state$parameters, g2, "b")
          x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")

          fit_timegrid <- x$fit.fun_timegrid(g2)
          eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
          if(!is.null(dx))
            eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
          
          eta_timegrid_alpha <- 0 
          for (i in names(xsmalpha)){
            g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
            if(i != "model.matrix"){
              Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
              Xalpha <- constrain(xsmalpha[[i]], Xalpha)
              if(fac)
                Xalpha <- Xalpha * xsmalpha[[i]]$by_timegrid
              if(tp)  
                Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
              
              eta_timegrid_alpha <- eta_timegrid_alpha + matrix(Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
            } else {
              eta_timegrid_alpha <- eta_timegrid_alpha + xsmalpha[[i]]$state$fitted_timegrid 
            }
          }
          
          eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
          fit <- x$fit.fun(x$X, g2)
          eta$mu <- eta$mu - fitted(x$state) + fit
          lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
          if(lp == -Inf) lp <- -1 * .Machine$double.xmax
          return(-1 * lp)
        }
        nu <- optimize(f = objfun.nu, interval = c(0, 1))$minimum
      } else {
        nu <- x$state$nu
      }
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      fit <- x$fit.fun(x$X, g2)
      fit_timegrid <- x$fit.fun_timegrid(g2)
      
      eta$mu <- eta$mu - fitted(x$state) + fit
      eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
      if(!is.null(dx))
        eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
      
      eta_timegrid_alpha <- 0 
      for (i in names(xsmalpha)){
        g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
        if(i != "model.matrix"){
          Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
          Xalpha <- constrain(xsmalpha[[i]], Xalpha)
          if(fac)
            Xalpha <- Xalpha * xsmalpha[[i]]$by_timegrid
          if(tp)  
            Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
          
          eta_timegrid_alpha <- eta_timegrid_alpha + matrix(Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
        } else {
          eta_timegrid_alpha <- eta_timegrid_alpha + xsmalpha[[i]]$state$fitted_timegrid 
        }
      }

      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
      edf1 <- sum_diag((-1 * xhess0) %*% Sigma)

      edf <- edf0 + edf1
      logLik <- get_LogPost(eta_timegrid, eta, 0)
      ic <- get.ic2(logLik, edf, length(eta$mu), criterion)

      if(!is.null(env$ic_val)) {
        if((ic < env$ic_val) & (ic < env$ic00_val)) {
          par[x$pid$b] <- g2
          opt_state <- list("parameters" = set.par(par, g2, "b"),
                            "fitted.values" = fit, "fitted_timegrid" = fit_timegrid,
                            "edf" = edf1, "hessian" = -1 * xhess,
                            "nu" = nu, "do.optim" = x$state$do.optim)
          assign("state", opt_state, envir = env)
          assign("ic_val", ic, envir = env)
        }
      } else assign("ic_val", ic, envir = env)
      return(ic)
    }
    
    assign("ic00_val", objfun1(get.state(x, "tau2")), envir = env)
    tau2 <- tau2.optim(objfun1, start = get.state(x, "tau2"), maxit = 1)
    if(!is.null(env$state))
      return(env$state)
    
    x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
  }
  
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  xhess <- xhess0 - x$hess(score = NULL, x$state$parameters, full = FALSE)

  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(-1 * xhess, index = NULL)
  Hs <- Sigma %*% xgrad
  
  ## Update regression coefficients.
  g <- get.state(x, "b")
  
  if(update.nu) {
    objfun2 <- function(nu) {
      g2 <- drop(g + nu * Hs)
      names(g2) <- names(g)
      x$state$parameters <- set.par(x$state$parameters, g2, "b")
      fit_timegrid <- x$fit.fun_timegrid(g2)
      eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
      if(!is.null(dx))
        eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g2)
      
      eta_timegrid_alpha <- 0 
      for (i in names(xsmalpha)){
        g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
        if(i != "model.matrix"){
          Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
          Xalpha <- constrain(xsmalpha[[i]], Xalpha)
          if(fac)
            Xalpha <- Xalpha * xsmalpha[[i]]$by_timegrid
          if(tp)  
            Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
          
          eta_timegrid_alpha <- eta_timegrid_alpha + matrix(Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
        } else {
          eta_timegrid_alpha <- eta_timegrid_alpha + xsmalpha[[i]]$state$fitted_timegrid 
        }
      }
      
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha  + eta_timegrid_dalpha * eta_timegrid_dmu
      fit <- x$fit.fun(x$X, g2)
      eta$mu <- eta$mu - fitted(x$state) + fit
      lp <- get_LogPost(eta_timegrid, eta, x$prior(x$state$parameters))
      if(lp == -Inf) lp <- -1 * .Machine$double.xmax
      return(-1 * lp)
    }
    x$state$nu <- optimize(f = objfun2, interval = c(0, 1))$minimum
  }
  
  g2 <- drop(g + x$state$nu * Hs)
  names(g2) <- names(g)
  x$state$parameters <- set.par(x$state$parameters, g2, "b")
  
  ## Update fitted values..
  fit_timegrid <- x$fit.fun_timegrid(g2)
  x$state$fitted.values <- x$fit.fun(x$X, g2)
  x$state$edf <- sum_diag((-1 * xhess0) %*% Sigma)
  x$state$hessian <- -1 * xhess
  
  return(x$state)
}

rowTensorProduct <- function(X1,X2)
{
  # from fdboost
  do.call("rbind", lapply(1:nrow(X1), function(i)
    rep(X1[i,], each=ncol(X2)) * rep(X2[i,], ncol(X1))
  ))
}

propose_jm_mu_nonlin <- function(x, y,
                                eta, eta_timegrid, eta_timegrid_lambda, eta_timegrid_mu, eta_timegrid_alpha,
                                eta_timegrid_dmu, eta_timegrid_dalpha,
                                width, sub, nu, status, id, dx = NULL, 
                                xsmalpha = NULL, knots = NULL, tp = FALSE, fac = FALSE, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  dX <- if(!is.null(dx)) dx$fit.fun_timegrid(NULL) else NULL
  gdim <- c(nrow(eta_timegrid), ncol(eta_timegrid))

  ## Timegrid lambda.
  eeta <- exp(eta_timegrid)
  
  ## evaluations of alpha
  eta_timegrid_alpha1 <- eta_timegrid_alpha2 <- 0 
  for (i in names(xsmalpha)){
    if(i != "model.matrix"){      # derivative for model.matrix is  0
      g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
      dXalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 1) 
      dXalpha <- constrain(xsmalpha[[i]], dXalpha)
      d2Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 2)
      d2Xalpha <- constrain(xsmalpha[[i]], d2Xalpha)
      
      if(fac){
        dXalpha <- dXalpha * xsmalpha[[i]]$by_timegrid
        d2Xalpha <- d2Xalpha * xsmalpha[[i]]$by_timegrid
      }
      
      if(tp){
        Xalpha2 <- xsmalpha[[1]]$margin[[2]]$fit.fun_timegrid(NULL)
        dXalpha <- rowTensorProduct(dXalpha, Xalpha2) 
        d2Xalpha <- rowTensorProduct(d2Xalpha, Xalpha2) 
      }
      eta_timegrid_alpha1 <- eta_timegrid_alpha1 + matrix(dXalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE) 
      eta_timegrid_alpha2 <- eta_timegrid_alpha2 + matrix(d2Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    }
  }
  
  eta_alpha1 <- eta_timegrid_alpha1[, gdim[2]]
  eta_alpha2 <- eta_timegrid_alpha2[, gdim[2]]
  
  ## Old logLik and prior.
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibeta <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  p1 <- x$prior(x$state$parameters)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_alpha1, width, exp(eta$gamma),
                 eeta * ((eta_timegrid_alpha1)^2 + eta_timegrid_alpha2),
                 index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  
  xgrad <- drop(t(x$X) %*% drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2) +
                  t(x$XT) %*% drop(eta_alpha1 * status) - int$grad)
  # if(!is.null(dx))
  #   xgrad <- drop(xgrad + t(dx$XT) %*% drop(eta$dalpha * status))
  
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  
  XWX <- if(is.null(x$sparse.setup$matrix)) {
    crossprod(x$X * (1 / exp(eta$sigma)^2), x$X)
  } else do.XWX(x$X, exp(eta$sigma)^2, index = x$sparse.setup$matrix)
  XAX <- crossprod(x$XT * drop(eta_alpha2 * status), x$XT)
  xhess <- -1 * XWX + XAX - int$hess 
  xhess <- xhess - x$hess(score = NULL, x$state$parameters, full = FALSE) 
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(-1 * xhess, index = NULL)

  ## Save old coefficients.
  g0 <- get.state(x, "b")
  
  ## Get new position.
  mu <- drop(g0 + 0.5 * nu * Sigma %*% xgrad)
  
  ## Sample new parameters.
  g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma, method="chol"))

  names(g) <- names(g0)
  x$state$parameters <- set.par(x$state$parameters, g, "b")
  
  ## Compute log priors.
  p2 <- x$prior(x$state$parameters)
  qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)
  
  ## Update additive predictors.
  fit_timegrid <- x$fit.fun_timegrid(g)
  eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid
  
  # if(!is.null(dx))
  #   eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g)

  eta_timegrid_alpha <- 0 
  for (i in names(xsmalpha)){
    g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
    if(i != "model.matrix"){
      Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
      Xalpha <- constrain(xsmalpha[[i]], Xalpha)
      if(fac)
        Xalpha <- Xalpha * xsmalpha[[i]]$by_timegrid
      if(tp)  
        Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
      
      eta_timegrid_alpha <- eta_timegrid_alpha + matrix(Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    } else {
      eta_timegrid_alpha <- eta_timegrid_alpha + xsmalpha[[i]]$state$fitted_timegrid 
    }
  }


  eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
  x$state$fitted_timegrid <- fit_timegrid

  fit <- drop(x$X %*% g)
  eta$mu <- eta$mu - fitted(x$state) + fit
  x$state$fitted.values <- fit
  
  ## New logLik.
  eeta <- exp(eta_timegrid)
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibetaprop <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  ## Prior prob.

  ## evaluations of alpha
  eta_timegrid_alpha1 <- eta_timegrid_alpha2 <- 0 
  for (i in names(xsmalpha)){
    if(i != "model.matrix"){      # derivative for model.matrix is  0
      g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
      dXalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 1) 
      dXalpha <- constrain(xsmalpha[[i]], dXalpha)
      d2Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 2)
      d2Xalpha <- constrain(xsmalpha[[i]], d2Xalpha)
      
      if(fac){
        dXalpha <- dXalpha * xsmalpha[[i]]$by_timegrid
        d2Xalpha <- d2Xalpha * xsmalpha[[i]]$by_timegrid
      }
      
      if(tp){
        Xalpha2 <- xsmalpha[[1]]$margin[[2]]$fit.fun_timegrid(NULL)
        dXalpha <- rowTensorProduct(dXalpha, Xalpha2) 
        d2Xalpha <- rowTensorProduct(d2Xalpha, Xalpha2) 
      }
      eta_timegrid_alpha1 <- eta_timegrid_alpha1 + matrix(dXalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE) 
      eta_timegrid_alpha2 <- eta_timegrid_alpha2 + matrix(d2Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    }
  }
  
  eta_alpha1 <- eta_timegrid_alpha1[, gdim[2]]
  eta_alpha2 <- eta_timegrid_alpha2[, gdim[2]]
  
  int <- survint(X, eeta * eta_timegrid_alpha1, width, exp(eta$gamma),
                 eeta * ((eta_timegrid_alpha1)^2 + eta_timegrid_alpha2),
                 index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  
  xgrad <- drop(t(x$X) %*% drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2) +
                  t(x$XT) %*% drop(eta_alpha1 * status) - int$grad)
  # if(!is.null(dx))
  #   xgrad <- drop(xgrad + t(dx$XT) %*% drop(eta$dalpha * status))
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  
  XWX <- if(is.null(x$sparse.setup$matrix)) {
    crossprod(x$X * (1 / exp(eta$sigma)^2), x$X)
  } else do.XWX(x$X, exp(eta$sigma)^2, index = x$sparse.setup$matrix)
  XAX <- crossprod(x$XT *  drop(eta_alpha2 * status), x$XT)
  xhess0 <- -1 * XWX + XAX - int$hess
  xhess <- xhess0 - x$hess(score = NULL, x$state$parameters, full = FALSE) 
  
  Sigma2 <- matrix_inv(-1 * xhess, index = NULL)
  if(inherits(Sigma2, "try-error")) {
    x$state$alpha <- 1
    return(x$state)
  }

  mu2 <- drop(g + 0.5 * nu * Sigma2 %*% xgrad)
  qbeta <- dmvnorm(g0, mean = mu2, sigma = Sigma2, log = TRUE)
  
  ## Save edf.
  x$state$edf <- sum_diag((-1 * xhess0) %*% Sigma2)

  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if(length(x$S) < 2 & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- 0.5 * crossprod(g, x$S[[1]]) %*% g + x$b
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "mu", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  return(x$state)
}


propose_jm_mu_nonlin_Matrix <- function(x, y,
                                       eta, eta_timegrid, eta_timegrid_lambda, 
                                       eta_timegrid_mu, eta_timegrid_alpha,
                                       eta_timegrid_dmu, eta_timegrid_dalpha,
                                       width, sub, nu, status, id, dx = NULL, 
                                       xsmalpha = NULL, knots = NULL, tp = FALSE, 
                                       fac = FALSE, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  dX <- if(!is.null(dx)) dx$fit.fun_timegrid(NULL) else NULL
  gdim <- c(nrow(eta_timegrid), ncol(eta_timegrid))
  
  ## Timegrid lambda.
  eeta <- exp(eta_timegrid)
  
  ## evaluations of alpha
  eta_timegrid_alpha1 <- eta_timegrid_alpha2 <- 0 
  for (i in names(xsmalpha)){
    if(i != "model.matrix"){      # derivative for model.matrix is  0
      g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
      dXalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 1) 
      dXalpha <- constrain(xsmalpha[[i]], dXalpha)
      d2Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 2)
      d2Xalpha <- constrain(xsmalpha[[i]], d2Xalpha)
      
      if(fac){
        dXalpha <- dXalpha * xsmalpha[[i]]$by_timegrid
        d2Xalpha <- d2Xalpha * xsmalpha[[i]]$by_timegrid
      }
      
      if(tp){
        Xalpha2 <- xsmalpha[[1]]$margin[[2]]$fit.fun_timegrid(NULL)
        dXalpha <- rowTensorProduct(dXalpha, Xalpha2) 
        d2Xalpha <- rowTensorProduct(d2Xalpha, Xalpha2) 
      }
      eta_timegrid_alpha1 <- eta_timegrid_alpha1 + matrix(dXalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE) 
      eta_timegrid_alpha2 <- eta_timegrid_alpha2 + matrix(d2Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    }
  }
  
  eta_alpha1 <- eta_timegrid_alpha1[, gdim[2]]
  eta_alpha2 <- eta_timegrid_alpha2[, gdim[2]]
  
  ## Old logLik and prior.
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibeta <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  p1 <- x$prior(x$state$parameters)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_alpha1, width, exp(eta$gamma),
                 eeta * ((eta_timegrid_alpha1)^2 + eta_timegrid_alpha2),
                 index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  
  xgrad <- crossprod(x$X, drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2)) +
    crossprod(x$XT,  drop(eta_alpha1 * status)) - int$grad
  # if(!is.null(dx))
  #   xgrad <- xgrad + crossprod(dx$XT, drop(eta$dalpha * status))
  
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  
  XWX <- crossprod(Diagonal(x = 1 / exp(eta$sigma)^2) %*% x$X, x$X)
  XAX <- crossprod(x$XT * drop(eta_alpha2 * status), x$XT)
  xhess <- -1 * XWX + XAX - int$hess
  
  xhess <- xhess - x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(-1 * xhess, index = NULL)
  
  ## Save old coefficients.
  g0 <- get.state(x, "b")
  
  ## Get new position.
  mu <- drop(g0 + nu * Sigma %*% xgrad)
  Sigma <- as.matrix(Sigma)
  sample.problem <- TRUE
  sample.iter <- 0
  xgrad_keep <- xgrad

  while(sample.problem & sample.iter < 4){
    sample.iter <- sample.iter + 1
      ## Sample new parameters using blockdiagonal structure.
      if(!is.null(x$sparse.setup$block.index)){
        lg <- lapply(1:length(x$sparse.setup$block.index), function(i){
          tmp <- x$sparse.setup$block.index[[i]]
          if(x$sparse.setup$is.diagonal){
            if(Sigma[tmp,tmp] < 0){
              Sigma[tmp,tmp] <- max(diag(Sigma))
            }
            drop(rnorm(n = 1, mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp])))
          } else{
            if(length(tmp) == 1){
              if(Sigma[tmp,tmp] < 0){
                Sigma[tmp,tmp] <- max(diag(Sigma))
              }
              drop(rnorm(n = 1, mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp])))
            } else {
              if(any(eigen(Sigma[tmp, tmp])[[1]] < 0))
                Sigma[tmp, tmp] <- Sigma[tmp, tmp] - min(eigen(Sigma[tmp, tmp])[[1]]) * diag(length(tmp))
              drop(rmvnorm(n = 1, mean = mu[tmp], sigma = Sigma[tmp,tmp], method="chol"))
            }
          } 
        })
        g <- unlist(lg)
      } else {
        g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma, method="chol"))
      }
      
      if(all(is.na(g))) {
        x$state$alpha <- 1
        return(x$state)
      }
      names(g) <- names(g0)
      x$state$parameters <- set.par(x$state$parameters, g, "b")
      
      ## Compute log priors.
      p2 <- x$prior(x$state$parameters)
      if(!is.null(x$sparse.setup$block.index)){
        lqbetaprop <- lapply(1:length(x$sparse.setup$block.index), function(i){
          tmp <- x$sparse.setup$block.index[[i]]
          if(x$sparse.setup$is.diagonal){
            if(Sigma[tmp,tmp] < 0){
              Sigma[tmp,tmp] <- -1 * Sigma[tmp,tmp]
            }
            drop(dnorm(g[tmp], mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp]), log = TRUE))
          } else{
            if(length(tmp) == 1){
              if(Sigma[tmp,tmp] < 0){
                Sigma[tmp,tmp] <- -1 * Sigma[tmp,tmp]
              }
              drop(dnorm(g[tmp], mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp]), log = TRUE))
            } else {
              #!# to be checked
              if(any(eigen(Sigma[tmp, tmp])[[1]] < 0))
                Sigma[tmp, tmp] <- Sigma[tmp, tmp] - min(eigen(Sigma[tmp, tmp])[[1]]) * diag(length(tmp))
              dmvnorm(g[tmp], mean = mu[tmp], sigma = Sigma[tmp,tmp], log = TRUE)
            }
          } 
        })
        qbetaprop <- sum(unlist(lqbetaprop))
      } else {
        qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)
      }
      
      ## Update additive predictors.
      fit_timegrid <- matrix(X %*% g, nrow = gdim[1], ncol = gdim[2], byrow = TRUE) # instead of fit_timegrid <- x$fit.fun_timegrid(g)
      eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid

      # if(!is.null(dx))
      #   eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g)
      eta_timegrid_alpha <- 0 
      for (i in names(xsmalpha)){
        g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
        if(i != "model.matrix"){
          Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
          Xalpha <- constrain(xsmalpha[[i]], Xalpha)
          if(fac)
            Xalpha <- Xalpha * xsmalpha[[i]]$by_timegrid
          if(tp)  
            Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
          
          eta_timegrid_alpha <- eta_timegrid_alpha + matrix(Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
        } else {
          eta_timegrid_alpha <- eta_timegrid_alpha + xsmalpha[[i]]$state$fitted_timegrid 
        }
      }
    
      eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
      x$state$fitted_timegrid <- fit_timegrid
      
      fit <- drop(x$X %*% g)
      eta$mu <- eta$mu - fitted(x$state) + fit
      x$state$fitted.values <- fit

      ## New logLik.
      eeta <- exp(eta_timegrid)
      int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
      pibetaprop <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
        exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
      
      ## Prior prob.
      ## evaluations of alpha
      eta_timegrid_alpha1 <- eta_timegrid_alpha2 <- 0 
      for (i in names(xsmalpha)){
        if(i != "model.matrix"){      # derivative for model.matrix is  0
          g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
          dXalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 1) 
          dXalpha <- constrain(xsmalpha[[i]], dXalpha)
          d2Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 2)
          d2Xalpha <- constrain(xsmalpha[[i]], d2Xalpha)
          
          if(fac){
            dXalpha <- dXalpha * xsmalpha[[i]]$by_timegrid
            d2Xalpha <- d2Xalpha * xsmalpha[[i]]$by_timegrid
          }
          
          if(tp){
            Xalpha2 <- xsmalpha[[1]]$margin[[2]]$fit.fun_timegrid(NULL)
            dXalpha <- rowTensorProduct(dXalpha, Xalpha2) 
            d2Xalpha <- rowTensorProduct(d2Xalpha, Xalpha2) 
          }
          eta_timegrid_alpha1 <- eta_timegrid_alpha1 + matrix(dXalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE) 
          eta_timegrid_alpha2 <- eta_timegrid_alpha2 + matrix(d2Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
        }
      }
      
      eta_alpha1 <- eta_timegrid_alpha1[, gdim[2]]
      eta_alpha2 <- eta_timegrid_alpha2[, gdim[2]]
      
      int <- survint(X, eeta * eta_timegrid_alpha1, width, exp(eta$gamma),
                     eeta * ((eta_timegrid_alpha1)^2 + eta_timegrid_alpha2),
                     index = x$sparse.setup[["mu.matrix"]],
                     dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                     if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
      
      xgrad <- crossprod(x$X, drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2)) +
        crossprod(x$XT,  drop(eta_alpha1 * status)) - int$grad
      # if(!is.null(dx))
      #   xgrad <- xgrad + crossprod(dx$XT, drop(eta$dalpha * status))
      xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
      
      XWX <- crossprod(Diagonal(x = 1 / exp(eta$sigma)^2) %*% x$X, x$X)
      XAX <- crossprod(x$XT * eta_alpha2, x$XT)
      xhess0 <- -1 * XWX + XAX - int$hess
      xhess <- xhess0 - x$hess(score = NULL, x$state$parameters, full = FALSE)
      
      # repeat if Hessian is not negative definite in block-samples
      if(is.null(x$sparse.setup$block.index)){
        sample.problem <- FALSE}
      else{
        eigen_blockhess <- sapply(1:length(x$sparse.setup$block.index), function(i){
          max(eigen(xhess[x$sparse.setup$block.index[[i]], x$sparse.setup$block.index[[i]] ])[[1]])})
        if(!any(eigen_blockhess > 0)){
          sample.problem <- FALSE
        } else {
          cases <- which(eigen_blockhess > 0)
          cases <- unlist(lapply(cases, function(i) x$sparse.setup$block.index[[i]]))
          cat("\n", x$label, "sample iteration: ", sample.iter, " cases: ", cases)
          # adjust steplength of problematic cases
          mu[cases] <- drop(g0[cases] + 0.1 * nu * Sigma[cases, cases] %*% xgrad_keep[cases])
        }
      }
  }
  
  Sigma2 <- matrix_inv(-1 * xhess, index = NULL)
  if(inherits(Sigma2, "try-error")) {
    x$state$alpha <- 0
    return(x$state)
  }
  
  mu2 <- drop(g + nu * Sigma2 %*% xgrad)
  Sigma2 <- as.matrix(Sigma2) 

  
  if(!is.null(x$sparse.setup$block.index)){
    lqbeta <- lapply(1:length(x$sparse.setup$block.index), function(i){
      tmp <- x$sparse.setup$block.index[[i]]
      if(x$sparse.setup$is.diagonal){
        if(Sigma2[tmp,tmp] < 0){
          Sigma2[tmp,tmp] <- max(diag(Sigma))
        }
        drop(dnorm(g0[tmp], mean = mu2[tmp], sd = sqrt(Sigma2[tmp,tmp]), log = TRUE))
      } else{
        if(length(tmp) == 1){
          drop(dnorm(g0[tmp], mean = mu2[tmp], sd = sqrt(Sigma2[tmp,tmp]), log = TRUE))
        } else {
          dmvnorm(g0[tmp], mean = mu2[tmp], sigma = Sigma2[tmp,tmp], log = TRUE)
        }
      } 
    })
    qbeta <- sum(unlist(lqbeta))
  } else {
    qbeta <- dmvnorm(g0, mean = mu2, sigma = as.matrix(Sigma2), log = TRUE)
  }
  
  ## Save edf.
  x$state$edf <- sum_diag((-1 * xhess0) %*% Sigma2)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if(length(x$S) < 2  & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- drop(0.5 * crossprod(g, x$S[[1]]) %*% g + x$b)
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "mu", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  return(x$state)
}

propose_jm_mu_nonlin_Matrix_block <- function(x, y,
                                       eta, eta_timegrid, eta_timegrid_lambda, 
                                       eta_timegrid_mu, eta_timegrid_alpha,
                                       eta_timegrid_dmu, eta_timegrid_dalpha,
                                       width, sub, nu, status, id, dx = NULL, 
                                       xsmalpha = NULL, knots = NULL, tp = FALSE, 
                                       fac = FALSE, ...)
{
  ## The time-dependent design matrix for the grid.
  X <- x$fit.fun_timegrid(NULL)
  dX <- if(!is.null(dx)) dx$fit.fun_timegrid(NULL) else NULL
  gdim <- c(nrow(eta_timegrid), ncol(eta_timegrid))
  
  ## Timegrid lambda.
  eeta <- exp(eta_timegrid)
  
  # if(sum(is.na(eta_timegrid_mu)) > 0) browser()
  ## evaluations of alpha
  eta_timegrid_alpha1 <- eta_timegrid_alpha2 <- 0 
  for (i in names(xsmalpha)){
    if(i != "model.matrix"){      # derivative for model.matrix is  0
      g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
      dXalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 1) 
      dXalpha <- constrain(xsmalpha[[i]], dXalpha)
      d2Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 2)
      d2Xalpha <- constrain(xsmalpha[[i]], d2Xalpha)
      
      if(fac){
        dXalpha <- dXalpha * xsmalpha[[i]]$by_timegrid
        d2Xalpha <- d2Xalpha * xsmalpha[[i]]$by_timegrid
      }
      
      if(tp){
        Xalpha2 <- xsmalpha[[1]]$margin[[2]]$fit.fun_timegrid(NULL)
        dXalpha <- rowTensorProduct(dXalpha, Xalpha2) 
        d2Xalpha <- rowTensorProduct(d2Xalpha, Xalpha2) 
      }
      eta_timegrid_alpha1 <- eta_timegrid_alpha1 + matrix(dXalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE) 
      eta_timegrid_alpha2 <- eta_timegrid_alpha2 + matrix(d2Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    }
  }
  
  eta_alpha1 <- eta_timegrid_alpha1[, gdim[2]]
  eta_alpha2 <- eta_timegrid_alpha2[, gdim[2]]
  
  ## Old logLik and prior.
  int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
  pibeta <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
    exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
  p1 <- x$prior(x$state$parameters)
  
  ## Compute gradient and hessian integrals.
  int <- survint(X, eeta * eta_timegrid_alpha1, width, exp(eta$gamma),
                 eeta * ((eta_timegrid_alpha1)^2 + eta_timegrid_alpha2),
                 index = x$sparse.setup[["mu.matrix"]],
                 dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                 if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
  
  xgrad <- crossprod(x$X, drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2)) +
    crossprod(x$XT,  drop(eta_alpha1 * status)) - int$grad
  # if(!is.null(dx))
  #   xgrad <- xgrad + crossprod(dx$XT, drop(eta$dalpha * status))
  
  xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
  
  XWX <- crossprod(Diagonal(x = 1 / exp(eta$sigma)^2) %*% x$X, x$X)
  XAX <- crossprod(x$XT * drop(eta_alpha2 * status), x$XT)
  xhess <- -1 * XWX + XAX - int$hess
  
  xhess <- xhess - x$hess(score = NULL, x$state$parameters, full = FALSE)
  
  ## Compute the inverse of the hessian.
  Sigma <- matrix_inv(-1 * xhess, index = NULL)
  
  ## Save old coefficients.
  g0 <- get.state(x, "b")
  
  ## Get new position.
  mu <- drop(g0 + nu * Sigma %*% xgrad)
  Sigma <- as.matrix(Sigma)
  sample.problem <- TRUE
  sample.iter <- 0
  xgrad_keep <- xgrad
  
  while(sample.problem & sample.iter < 4){
    sample.iter <- sample.iter + 1
    ## Sample new parameters using blockdiagonal structure.
    if(!is.null(x$sparse.setup$block.index)){
      lg <- lapply(1:length(x$sparse.setup$block.index), function(i){
        tmp <- x$sparse.setup$block.index[[i]]
        if(x$sparse.setup$is.diagonal){
          if(Sigma[tmp,tmp] < 0){
            Sigma[tmp,tmp] <- max(diag(Sigma))
          }
          drop(rnorm(n = 1, mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp])))
        } else{
          if(length(tmp) == 1){
            if(Sigma[tmp,tmp] < 0){
              Sigma[tmp,tmp] <- max(diag(Sigma))
            }
            drop(rnorm(n = 1, mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp])))
          } else {
            if(any(eigen(Sigma[tmp, tmp])[[1]] < 0))
              Sigma[tmp, tmp] <- Sigma[tmp, tmp] - min(eigen(Sigma[tmp, tmp])[[1]]) * diag(length(tmp))
            drop(rmvnorm(n = 1, mean = mu[tmp], sigma = Sigma[tmp,tmp], method="chol"))
          }
        } 
      })
      g <- unlist(lg)
    } else {
      g <- drop(rmvnorm(n = 1, mean = mu, sigma = Sigma, method="chol"))
    }
    
    if(all(is.na(g))) {
      x$state$alpha <- 1
      return(x$state)
    }
    names(g) <- names(g0)
    x$state$parameters <- set.par(x$state$parameters, g, "b")
    
    ## Compute log priors.
    p2 <- x$prior(x$state$parameters)
    if(!is.null(x$sparse.setup$block.index)){
      lqbetaprop <- lapply(1:length(x$sparse.setup$block.index), function(i){
        tmp <- x$sparse.setup$block.index[[i]]
        if(x$sparse.setup$is.diagonal){
          if(Sigma[tmp,tmp] < 0){
            Sigma[tmp,tmp] <- -1 * Sigma[tmp,tmp]
          }
          drop(dnorm(g[tmp], mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp]), log = TRUE))
        } else{
          if(length(tmp) == 1){
            if(Sigma[tmp,tmp] < 0){
              Sigma[tmp,tmp] <- -1 * Sigma[tmp,tmp]
            }
            drop(dnorm(g[tmp], mean = mu[tmp], sd = sqrt(Sigma[tmp,tmp]), log = TRUE))
          } else {
            #!# to be checked
            if(any(eigen(Sigma[tmp, tmp])[[1]] < 0))
              Sigma[tmp, tmp] <- Sigma[tmp, tmp] - min(eigen(Sigma[tmp, tmp])[[1]]) * diag(length(tmp))
            dmvnorm(g[tmp], mean = mu[tmp], sigma = Sigma[tmp,tmp], log = TRUE)
          }
        } 
      })
      qbetaprop <- sum(unlist(lqbetaprop))
    } else {
      qbetaprop <- dmvnorm(g, mean = mu, sigma = Sigma, log = TRUE)
    }
    
    ## Update additive predictors.
    #?#
    fit_timegrid <- matrix(X %*% g, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
    # fit_timegrid <- x$fit.fun_timegrid(g)
    eta_timegrid_mu <- eta_timegrid_mu - x$state$fitted_timegrid + fit_timegrid

    # if(!is.null(dx))
    #   eta_timegrid_dmu <- eta_timegrid_dmu - dx$state$fitted_timegrid + dx$fit.fun_timegrid(g)
    eta_timegrid_alpha <- 0 
    for (i in names(xsmalpha)){
      g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
      if(i != "model.matrix"){
        Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 0)
        Xalpha <- constrain(xsmalpha[[i]], Xalpha)
        if(fac)
          Xalpha <- Xalpha * xsmalpha[[i]]$by_timegrid
        if(tp)  
          Xalpha <- rowTensorProduct(Xalpha, Xalpha2) 
        
        eta_timegrid_alpha <- eta_timegrid_alpha + matrix(Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
      } else {
        eta_timegrid_alpha <- eta_timegrid_alpha + xsmalpha[[i]]$state$fitted_timegrid 
      }
    }
    
    eta_timegrid <- eta_timegrid_lambda + eta_timegrid_alpha + eta_timegrid_dalpha * eta_timegrid_dmu
    x$state$fitted_timegrid <- fit_timegrid
    
    fit <- drop(x$X %*% g)
    eta$mu <- eta$mu - fitted(x$state) + fit
    x$state$fitted.values <- fit
    
    ## New logLik.
    eeta <- exp(eta_timegrid)
    int0 <- width * (0.5 * (eeta[, 1] + eeta[, sub]) + apply(eeta[, 2:(sub - 1)], 1, sum))
    pibetaprop <- sum((eta_timegrid[,ncol(eta_timegrid)] + eta$gamma) * status, na.rm = TRUE) -
      exp(eta$gamma) %*% int0 + sum(dnorm(y[, "obs"], mean = eta$mu, sd = exp(eta$sigma), log = TRUE))
    
    ## Prior prob.
    ## evaluations of alpha
    eta_timegrid_alpha1 <- eta_timegrid_alpha2 <- 0 
    for (i in names(xsmalpha)){
      if(i != "model.matrix"){      # derivative for model.matrix is  0
        g_a <- get.par(xsmalpha[[i]]$state$parameters, "b")
        dXalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 1) 
        dXalpha <- constrain(xsmalpha[[i]], dXalpha)
        d2Xalpha <- modSplineDesign(knots, as.vector(t(eta_timegrid_mu)), derivs = 2)
        d2Xalpha <- constrain(xsmalpha[[i]], d2Xalpha)
        
        if(fac){
          dXalpha <- dXalpha * xsmalpha[[i]]$by_timegrid
          d2Xalpha <- d2Xalpha * xsmalpha[[i]]$by_timegrid
        }
        
        if(tp){
          Xalpha2 <- xsmalpha[[1]]$margin[[2]]$fit.fun_timegrid(NULL)
          dXalpha <- rowTensorProduct(dXalpha, Xalpha2) 
          d2Xalpha <- rowTensorProduct(d2Xalpha, Xalpha2) 
        }
        eta_timegrid_alpha1 <- eta_timegrid_alpha1 + matrix(dXalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE) 
        eta_timegrid_alpha2 <- eta_timegrid_alpha2 + matrix(d2Xalpha %*% g_a, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)
      }
    }
    
    eta_alpha1 <- eta_timegrid_alpha1[, gdim[2]]
    eta_alpha2 <- eta_timegrid_alpha2[, gdim[2]]
    
    int <- survint(X, eeta * eta_timegrid_alpha1, width, exp(eta$gamma),
                   eeta * ((eta_timegrid_alpha1)^2 + eta_timegrid_alpha2),
                   index = x$sparse.setup[["mu.matrix"]],
                   dX, if(!is.null(dX)) eeta * eta_timegrid_dalpha else NULL,
                   if(!is.null(dX)) eeta * eta_timegrid_dalpha^2 else NULL)
    
    xgrad <- crossprod(x$X, drop((y[, "obs"]  - eta$mu) / exp(eta$sigma)^2)) +
      crossprod(x$XT,  drop(eta_alpha1 * status)) - int$grad
    # if(!is.null(dx))
    #   xgrad <- xgrad + crossprod(dx$XT, drop(eta$dalpha * status))
    xgrad <- xgrad + x$grad(score = NULL, x$state$parameters, full = FALSE)
    
    XWX <- crossprod(Diagonal(x = 1 / exp(eta$sigma)^2) %*% x$X, x$X)
    XAX <- crossprod(x$XT * eta_alpha2, x$XT)
    xhess0 <- -1 * XWX + XAX - int$hess
    xhess <- xhess0 - x$hess(score = NULL, x$state$parameters, full = FALSE)
    
    # repeat if Hessian is not negative definite in block-samples
    if(is.null(x$sparse.setup$block.index)){
      sample.problem <- FALSE}
    else{
      eigen_blockhess <- sapply(1:length(x$sparse.setup$block.index), function(i){
        max(eigen(xhess[x$sparse.setup$block.index[[i]], x$sparse.setup$block.index[[i]] ])[[1]])})
      if(!any(eigen_blockhess > 0)){
        sample.problem <- FALSE
      } else {
        cases <- which(eigen_blockhess > 0)
        cases <- unlist(lapply(cases, function(i) x$sparse.setup$block.index[[i]]))
        cat("\n", x$label, "sample iteration: ", sample.iter, " cases: ", cases)
        # adjust steplength of problematic cases
        mu[cases] <- drop(g0[cases] + 0.1 * nu * Sigma[cases, cases] %*% xgrad_keep[cases])
      }
    }
  }
  
  Sigma2 <- matrix_inv(-1 * xhess, index = NULL)
  if(inherits(Sigma2, "try-error")) {
    x$state$alpha <- 0
    return(x$state)
  }
  
  mu2 <- drop(g + nu * Sigma2 %*% xgrad)
  Sigma2 <- as.matrix(Sigma2) 
  
  
  if(!is.null(x$sparse.setup$block.index)){
    lqbeta <- lapply(1:length(x$sparse.setup$block.index), function(i){
      tmp <- x$sparse.setup$block.index[[i]]
      if(x$sparse.setup$is.diagonal){
        if(Sigma2[tmp,tmp] < 0){
          Sigma2[tmp,tmp] <- max(diag(Sigma))
        }
        drop(dnorm(g0[tmp], mean = mu2[tmp], sd = sqrt(Sigma2[tmp,tmp]), log = TRUE))
      } else{
        if(length(tmp) == 1){
          drop(dnorm(g0[tmp], mean = mu2[tmp], sd = sqrt(Sigma2[tmp,tmp]), log = TRUE))
        } else {
          dmvnorm(g0[tmp], mean = mu2[tmp], sigma = Sigma2[tmp,tmp], log = TRUE)
        }
      } 
    })
    qbeta <- sum(unlist(lqbeta))
  } else {
    qbeta <- dmvnorm(g0, mean = mu2, sigma = as.matrix(Sigma2), log = TRUE)
  }
  
  ## Save edf.
  x$state$edf <- sum_diag((-1 * xhess0) %*% Sigma2)
  
  ## Sample variance parameter.
  if(!x$fixed & is.null(x$sp) & length(x$S)) {
    if(length(x$S) < 2  & (attr(x$prior, "var_prior") == "ig")) {
      g <- get.par(x$state$parameters, "b")
      a <- x$rank / 2 + x$a
      b <- drop(0.5 * crossprod(g, x$S[[1]]) %*% g + x$b)
      tau2 <- 1 / rgamma(1, a, b)
      x$state$parameters <- set.par(x$state$parameters, tau2, "tau2")
    } else {
      i <- grep("tau2", names(x$state$parameters))
      for(j in i) {
        x$state$parameters <- uni.slice(x$state$parameters, x, NULL, NULL,
                                        NULL, id = "mu", j, logPost = uni.slice_tau2_logPost, lower = 0, ll = 0)
      }
    }
  }
  
  ## Compute acceptance probablity.
  x$state$alpha <- drop((pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1))
  return(x$state)
}
