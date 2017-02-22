## From mgcv.
expand.t2.smooths <- function(sm) 
{
  m <- length(sm)
  not.needed <- TRUE
  for(i in 1:m) if (inherits(sm[[i]], "t2.smooth") && length(sm[[i]]$S) > 1) {
    not.needed <- FALSE
    break
  }
  if(not.needed) 
    return(NULL)
  smr <- list()
  k <- 0
  for(i in 1:m) {
    if(inherits(sm[[i]], "t2.smooth")) {
      smi <- split.t2.smooth(sm[[i]])
      comp.ind <- (k + 1):(k + length(smi))
      for(j in 1:length(smi)) {
        k <- k + 1
        smr[[k]] <- smi[[j]]
        smr[[k]]$comp.ind <- comp.ind
      }
    } else {
      k <- k + 1
      smr[[k]] <- sm[[i]]
    }
  }
  smr
}

split.t2.smooth <- function (object) 
{
  if(!inherits(object, "t2.smooth")) 
    return(object)
  ind <- 1:ncol(object$S[[1]])
  ind.para <- object$first.para:object$last.para
  sm <- list()
  sm[[1]] <- object
  St <- object$S[[1]] * 0
  for(i in 1:length(object$S)) {
    indi <- ind[diag(object$S[[i]]) != 0]
    label <- paste(object$label, ".frag", i, sep = "")
    sm[[i]] <- list(S = list(object$S[[i]][indi, indi]), 
      first.para = min(ind.para[indi]), last.para = max(ind.para[indi]), 
      fx = object$fx[i], fixed = object$fx[i], sp = object$sp[i], 
      null.space.dim = 0, df = length(indi), rank = object$rank[i], 
      label = label, S.scale = object$S.scale[i])
    class(sm[[i]]) <- "t2.frag"
    St <- St + object$S[[i]]
  }
  i <- length(object$S) + 1
  indi <- ind[diag(St) == 0]
  if(length(indi)) {
    label <- paste(object$label, ".frag", i, sep = "")
    sm[[i]] <- list(S = NULL, first.para = min(ind.para[indi]), 
      last.para = max(ind.para[indi]), fx = TRUE, fixed = TRUE, 
      null.space.dim = 0, label = label, df = length(indi))
    class(sm[[i]]) <- "t2.frag"
  }
  sm
}


smooth2random <- function(object, vnames, type = 1) {
  UseMethod("smooth2random")
}

smooth2random.t2.smooth <- function (object, vnames, type = 1) 
{
    if (object$fixed) 
        return(list(fixed = TRUE, Xf = object$X))
    fixed <- rep(TRUE, ncol(object$X))
    random <- list()
    diagU <- rep(1, ncol(object$X))
    ind <- 1:ncol(object$X)
    pen.ind <- ind * 0
    n.para <- 0
    for (i in 1:length(object$S)) {
        indi <- ind[diag(object$S[[i]]) != 0]
        pen.ind[indi] <- i
        X <- object$X[, indi, drop = FALSE]
        D <- diag(object$S[[i]])[indi]
        diagU[indi] <- 1/sqrt(D)
        X <- X %*% diag(diagU[indi])
        fixed[indi] <- FALSE
        term.name <- new.name("Xr", vnames)
        group.name <- new.name("g", vnames)
        vnames <- c(vnames, term.name, group.name)
        if (type == 1) {
            form <- as.formula(paste("~", term.name, "-1", sep = ""), 
                env = .GlobalEnv)
            random[[i]] <- pdIdnot(form)
            names(random)[i] <- group.name
            attr(random[[i]], "group") <- factor(rep(1, nrow(X)))
            attr(random[[i]], "Xr.name") <- term.name
            attr(random[[i]], "Xr") <- X
        }
        else {
            random[[i]] <- X
            names(random)[i] <- term.name
            attr(random[[i]], "s.label") <- object$label
        }
        n.para <- n.para + ncol(X)
    }
    if (sum(fixed)) {
        Xf <- object$X[, fixed, drop = FALSE]
    }
    else Xf <- matrix(0, nrow(object$X), 0)
    list(rand = random, trans.D = diagU, Xf = Xf, fixed = FALSE, 
        rind = 1:n.para, rinc = rep(n.para, n.para), pen.ind = pen.ind)
}


smooth2random.mgcv.smooth <- function (object, vnames, type = 1) 
{
    if (object$fixed) 
        return(list(fixed = TRUE, Xf = object$X))
    if (length(object$S) > 1) 
        stop("Can not convert this smooth class to a random effect")
    ev <- eigen(object$S[[1]], symmetric = TRUE)
    null.rank <- object$df - object$rank
    p.rank <- object$rank
    if (p.rank > ncol(object$X)) 
        p.rank <- ncol(object$X)
    U <- ev$vectors
    D <- c(ev$values[1:p.rank], rep(1, null.rank))
    D <- 1/sqrt(D)
    UD <- t(t(U) * D)
    X <- object$X %*% UD
    if (p.rank < object$df) 
        Xf <- X[, (p.rank + 1):object$df, drop = FALSE]
    else Xf <- matrix(0, nrow(object$X), 0)
    term.name <- new.name("Xr", vnames)
    if (type == 1) {
        form <- as.formula(paste("~", term.name, "-1", sep = ""), 
            env = .GlobalEnv)
        random <- list(pdIdnot(form))
        group.name <- new.name("g", vnames)
        names(random) <- group.name
        attr(random[[1]], "group") <- factor(rep(1, nrow(X)))
        attr(random[[1]], "Xr.name") <- term.name
        attr(random[[1]], "Xr") <- X[, 1:p.rank, drop = FALSE]
    }
    else {
        random <- list(X[, 1:p.rank, drop = FALSE])
        names(random)[1] <- term.name
        attr(random[[1]], "s.label") <- object$label
    }
    rind <- 1:p.rank
    pen.ind <- rep(0, ncol(object$X))
    pen.ind[rind] <- 1
    rinc <- rep(p.rank, p.rank)
    list(rand = random, Xf = Xf, trans.U = U, trans.D = D, fixed = FALSE, 
        rind = rind, rinc = rinc, pen.ind = pen.ind)
}

smooth2random.tensor.smooth <- function (object, vnames, type = 1) 
{
    if (type == 2) 
        stop("te smooths not useable with gamm4: use t2 instead")
    if (sum(object$fx) == length(object$fx)) 
        return(list(fixed = TRUE, Xf = object$X))
    else if (sum(object$fx) != 0) 
        warning("gamm can not fix only some margins of tensor product.")
    sum.S <- object$S[[1]]/mean(abs(object$S[[1]]))
    if (length(object$S) > 1) 
        for (l in 2:length(object$S)) {
            sum.S <- sum.S + object$S[[l]]/mean(abs(object$S[[l]]))
        }
    null.rank <- object$null.space.dim
    ev <- eigen(sum.S, symmetric = TRUE)
    p.rank <- ncol(object$X) - null.rank
    if (p.rank > ncol(object$X)) 
        p.rank <- ncol(object$X)
    U <- ev$vectors
    D <- c(ev$values[1:p.rank], rep(1, null.rank))
    if (sum(D <= 0)) 
        stop("Tensor product penalty rank appears to be too low: please email Simon.Wood@R-project.org with details.")
    U <- U
    X <- object$X %*% U
    if (p.rank < ncol(X)) 
        Xf <- X[, (p.rank + 1):ncol(X), drop = FALSE]
    else Xf <- matrix(0, nrow(X), 0)
    for (l in 1:length(object$S)) {
        object$S[[l]] <- (t(U) %*% object$S[[l]] %*% U)[1:p.rank, 
            1:p.rank]
        object$S[[l]] <- (object$S[[l]] + t(object$S[[l]]))/2
    }
    term.name <- new.name("Xr", vnames)
    form <- as.formula(paste("~", term.name, "-1", sep = ""), 
        env = .GlobalEnv)
    attr(form, "S") <- object$S
    random <- list(pdTens(form))
    group.name <- new.name("g", vnames)
    names(random) <- group.name
    attr(random[[1]], "group") <- factor(rep(1, nrow(X)))
    attr(random[[1]], "Xr.name") <- term.name
    attr(random[[1]], "Xr") <- X[, 1:p.rank, drop = FALSE]
    rind <- 1:p.rank
    rinc <- rep(p.rank, p.rank)
    list(rand = random, Xf = Xf, trans.U = U, trans.D = rep(1, 
        ncol(U)), fixed = FALSE, rind = rind, rinc = rinc)
}


.ringDirxy2 <- function(xy) 
{
    a <- xy[, 1]
    b <- xy[, 2]
    nvx <- length(b)
    if ((a[1] == a[nvx]) && (b[1] == b[nvx])) {
        a <- a[-nvx]
        b <- b[-nvx]
        nvx <- nvx - 1
    }
    if (nvx < 3) 
        return(1)
    tX <- 0
    dfYMax <- max(b)
    ti <- 1
    for (i in 1:nvx) {
        if (b[i] == dfYMax && a[i] > tX) 
            ti <- i
    }
    if ((ti > 1) & (ti < nvx)) {
        dx0 = a[ti - 1] - a[ti]
        dx1 = a[ti + 1] - a[ti]
        dy0 = b[ti - 1] - b[ti]
        dy1 = b[ti + 1] - b[ti]
    }
    else if (ti == nvx) {
        dx0 = a[ti - 1] - a[ti]
        dx1 = a[1] - a[ti]
        dy0 = b[ti - 1] - b[ti]
        dy1 = b[1] - b[ti]
    }
    else {
        dx1 = a[2] - a[1]
        dx0 = a[nvx] - a[1]
        dy1 = b[2] - b[1]
        dy0 = b[nvx] - b[1]
    }
    v3 = ((dx0 * dy1) - (dx1 * dy0))
    if (v3 > 0) 
        return(as.integer(1))
    else return(as.integer(-1))
}


plot.density2 <- function (x, main = NULL, xlab = NULL, ylab = "Density", type = "l", 
    zero.line = TRUE, ...) 
{
    if (is.null(xlab)) 
        xlab <- paste("N =", x$n, "  Bandwidth =", formatC(x$bw))
    if (is.null(main)) 
        main <- deparse(x$call)
    plot.default(x, main = main, xlab = xlab, ylab = ylab, type = type, 
        ...)
    if (zero.line) 
        abline(h = 0, lwd = 0.1, col = "gray")
    invisible(NULL)
}

