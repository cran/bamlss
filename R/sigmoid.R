exp2 <- function(x) {
  i <- x > 20
  x[i] <- 20
  x <- exp(x)
  x
}

sigmoid.fit.w <- function(X, y, weights = NULL, j = 2, ...) {
  if(any(is.na(y)) | any(is.na(X)))
    stop("NA values in data!")
  nr <- nrow(X)
  nc <- ncol(X)
  if(is.null(weights))
    weights <- rep(1, nr)

  X <- X[, -1, drop = FALSE]
  j <- j - 1

  y2 <- scale2(y, 0.01, 0.99)

  objfun <- function(w) {
    ws <- w[1]
    w <- w[-1]
    X[, j] <- X[, j] - ws
    wR <- 1 / (1 + exp(-(X %*% w)))
    wL <- 1 - wR
    X <- cbind(wL, wR)
    XX <- crossprod(X * weights)
    beta <- chol2inv(chol(XX + diag(1e-3, 2))) %*% t(X * weights) %*% y
    fit <- drop(X %*% beta)
    sum((y2 - fit)^2)
  }

  opt <- optim(c(sum(X[,j]*weights)/sum(weights), rep(0, nc - 1)), fn = objfun, gr = NULL, method = "L-BFGS-B")

  w <- opt$par
  ws <- w[1]
  w <- w[-1]
  X[, j] <- X[, j] - ws
  wR <- 1 / (1 + exp(-(X %*% w)))
  wL <- 1 - wR
  X <- cbind(wL, wR)
  XX <- crossprod(X * weights)
  beta <- drop(chol2inv(chol(XX + diag(1e-3, 2))) %*% t(X * weights) %*% y)
  fit <- drop(X %*% beta)

  rval <- list(
    "fitted.values" = fit,
    "residuals" = y - fit,
    "coefficients" = w,
    "w" = list("L" = wL, "R" = wR)
  )

  rval$residuals <- y - rval$fitted.values
  class(rval) <- "sigmoid"
  rval
}

sigmoid.fit <- function(X, y, weights = NULL) {
  if(any(is.na(y)) | any(is.na(X)))
    stop("NA values in data!")
  nr <- nrow(X)
  nc <- ncol(X)
  if(is.null(weights))
    weights <- rep(1, nr)
  objfun <- function(w) {
    sum((weights * (y - (w[1] + w[2] / (1 + exp2(-(X %*% w[-c(1:2)]))))))^2)
  }
  gradfun <- function(w) {
    ez <- exp2(-(X %*% w[-c(1:2)]))
    ez1 <- 1 + ez
    ez2 <- 1 / ez1
    eta <- w[1] + w[2] / ez1
    s1 <- -(2 * (y - eta))
    g <- matrix(0, nrow = nr, ncol = nc + 2)
    g[, 1] <- s1
    g[, 2] <- s1 * ez2
    for(j in 1:nc) {
      g[, j + 2] <- s1 * w[2] * ez2 * (1 - ez2) * X[, j]
    }
    colSums(g * weights)
  }
  opt <- optim(rep(0.1, ncol(X) + 2), fn = objfun, gr = gradfun, method = "L-BFGS-B")
  w <- opt$par
  rval <- list(
    "fitted.values" = drop(w[1] + w[2] / (1 + exp(-(X %*% w[-c(1:2)])))),
    "fit.fun" = function(X) drop(w[1] + w[2] / (1 + exp(-(X %*% w[-c(1:2)])))),
    "coefficients" = w[-c(1:2)],
    "weights" = weights
  )
  rval$rank <- 2 + ncol(X)
  rval$residuals <- y - rval$fitted.values
  class(rval) <- "sigmoid"
  rval
}

logLik.sigmoid <- function(object, ...) {
  res <- object$residuals
  w <- object$weights
  N <- length(res)
  val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w * res^2))))
  attr(val, "nall") <- N
  attr(val, "nobs") <- N
  attr(val, "df") <- object$rank + 1
  class(val) <- "logLik"
  val
}

sigmoid <- function(x, ...) {
  UseMethod("sigmoid")
}

sigmoid.formula <- function(formula, data, weights, ..., subset, na.action, contrasts = NULL) 
{
  class.ind <- function(cl) {
    n <- length(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1L:n) + n * (as.vector(unclass(cl)) - 1L)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    x
  }
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval.parent(m$data))) 
    m$data <- as.data.frame(data)
  m$... <- m$contrasts <- NULL
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  x <- model.matrix(Terms, m, contrasts)
  cons <- attr(x, "contrast")
  w <- model.weights(m)
  if(length(w) == 0L) 
    w <- rep(1, nrow(x))
  y <- model.response(m)
  res <- sigmoid.default(x, y, w, ...)
  res$terms <- Terms
  res$coefnames <- colnames(x)
  res$call <- match.call()
  res$na.action <- attr(m, "na.action")
  res$contrasts <- cons
  res$xlevels <- .getXlevels(Terms, m)
  class(res) <- c("sigmoid.formula", "sigmoid")
  res
}

sigmoid.default <- function(x, y, weights, ...) {
  sigmoid.fit(x, y, weights)
}

estfun.sigmoid <- function(x, ...) {
  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  if(any(alias <- is.na(coef(x)))) 
    xmat <- xmat[, !alias, drop = FALSE]
  wts <- weights(x)
  if(is.null(wts)) 
    wts <- 1
  res <- residuals(x)
  rval <- as.vector(res) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  return(rval)
}

stree <- function(x, y, k = 20, verbose = TRUE, plot = FALSE, min = 1, ...) {
  if(is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
    X <- cbind(1, x)
  } else {
    x <- as.matrix(x)
    has_itcpt <- apply(x, 2, function(x) all(x == 1) )
    if(!any(has_itcpt)) {
      X <- cbind(1, x)
    } else {
      x <- x[, -which(has_itcpt), drop = FALSE]
      X <- cbind(1, x)
    }
  }
  lw <- data.frame(rep(1, length(y)))
  K <- 0
  y0 <- y
  fit <- 0
  y <- y - fit
  while(K < k) {
    crit <- rep(NA, ncol(lw))
    for(j in 1:ncol(lw)) {
      if(min < sum(lw[, j])) {
        ws <- lw[, j]
        sfit <- sigmoid.fit.w(X, y, weights = ws)
        lw[, j] <- do.call("cbind", sfit$w) * ws
        b <- lm.fit(as.matrix(lw), y)
        crit[j] <- sum(residuals(b)^2)
        lw[, j] <- ws
      }
    }
    if(all(is.na(crit))) {
      warning("all Nas!")
      break
    }
    j <- which.min(crit)
    sfit <- sigmoid.fit.w(X, y, weights = lw[, j])
    lw[, j] <- do.call("cbind", sfit$w) * lw[, j]
    b <- lm.fit(as.matrix(lw), y)
    fit <- fit + fitted(b)
    y <- y - fitted(b)
    lw <- as.data.frame(as.matrix(lw))
    K <- ncol(lw)
  }
  lw <- as.matrix(lw)
  b <- lm.fit(lw, y0)
#  par(mfrow = c(3, 1))
#  plot(d$x, y)
#  plot(d$x, d$y)
#  plot2d(fitted(b) ~ d$x, add = TRUE, col.lines = 4, lwd = 2)
#  plot2d(lw ~ d$x, scheme = 1, col.lines = rainbow_hcl(ncol(lw)), lwd = 3)
  b$weights <- lw
  return(b)
}

build_net_w <- function(X, y, k = 10, weights = NULL, wts = NULL, maxit = 10, ...) {
  stopifnot(requireNamespace("nnet"))
  j <- grep("bw", names(wts))
  wts <- as.numeric(c(wts[j], runif(1, -0.5, 0.5), wts[-j]))
  m <- nnet::nnet(X[, -1, drop = FALSE], y, weights = weights,
    linout = TRUE, maxit = maxit, size = k, Wts = wts,
    MaxNWts = 20000, trace = FALSE, decay = 0.0001)
  w <- extract_nnet_weights(coef(m))
  return(w)
}


build_net_w1 <- function(X, y, k = 10, n = 10, weights = NULL, ...) {
  ind <- 1:nrow(X)
  tX <- t(X)
  w <- NULL
#plot(X[,2], y, ylim = c(0, 1))
  i <- -1
  n[which.min(n)] <- min(c(max(c(min(n), ceiling(ncol(X) * 6))), length(y)))
  H <- NULL
  rss <- 1e+20
  eps <- 1
  while(i < k) {
    j <- sample(ind, size = 1)
    tx <- as.numeric(X[j, , drop = FALSE])
    cs <- colSums((tX - tx)^2)
    if(length(n) < 2)
      n2 <- n
    else
      n2 <- sample(n, size = 1)
    take <- order(cs)[1:n2]
    yn <- scale2(y[take], 0.01, 0.99)
    xn <- X[take, , drop = FALSE]
    m <- glm.fit(xn, yn, family = binomial())
    wm <- coef(m)
    if(!any(is.na(wm))) {
      w <- cbind(w, wm)
      H2 <- H
      H2 <- cbind(H2, 1 / (1 + exp(-(X %*% wm))))
      b <- lm.wfit(H2, y, weights)
##hist(residuals(b), breaks = "Scott", xlim = c(-0.5, 0.5))
      rss1 <- sum(residuals(b)^2)
      eps1 <- (rss - rss1)/rss
      if(eps1 <= eps) {
        H <- cbind(H, 1 / (1 + exp(-(X %*% wm))))
        w <- cbind(w, wm)
        rss <- rss1
        eps <- eps * 0.99
##print(eps)
      }
    }
    if(!is.null(w))
      i <- ncol(w)
  }
  return(w)
}



build_net_w0 <- function(X, y, k = 10, n = 10, linear = FALSE, ...) {
  ind <- 1:nrow(X)
  tX <- t(X)
  w <- NULL
#plot(X[,2], y, ylim = c(0, 1))
  i <- -1
  n[which.min(n)] <- max(c(min(n), ncol(X) * 4))
  while(i < k) {
    j <- sample(ind, size = 1)
    tx <- as.numeric(X[j, , drop = FALSE])
    cs <- colSums((tX - tx)^2)
    if(length(n) < 2)
      n2 <- n
    else
      n2 <- sample(n, size = 1)
    take <- order(cs)[1:n2]
    yn <- y[take]
    xn <- X[take, , drop = FALSE]
    if(linear) {
      m <- lm.fit(xn, scale2(yn, 0.001, 0.999))
      cm <- coef(m)
      if(!any(is.na(cm))) {
        wm <- cm[-1] * 4
        wm <- c(-1 * sum(wm * tx[-1]), wm)
        w <- cbind(w, wm)
      }
    } else {
      m <- glm.fit(xn, scale2(yn, 0.01, 0.99), family = binomial())
      wm <- coef(m)
      if(!any(is.na(wm)))
        w <- cbind(w, wm)
    }
#X2 <- cbind(1, seq(-3, 3, length = 500))
#fit <- drop(1/(1 + exp(-c(X2 %*% wm))))
#plot2d(fit ~ X2[, 2], add = TRUE, col.lines = 4)
#points(xn[,2], yn, col = 2, pch = 16)
    if(!is.null(w))
      i <- ncol(w)
  }
  return(w)
}


extract_nnet_weights <- function(x)
{
  if(inherits(x, "nnet"))
    x <- coef(x)
  nw <- nw0 <- names(x)
  nw0 <- sapply(strsplit(nw0, "->", fixed = TRUE), function(x) { x[2] })
  nw0 <- paste0("->", nw0)
  nw <- grep("->h", nw, fixed = TRUE, value = TRUE)
  nw <- sapply(strsplit(nw, "->", fixed = TRUE), function(x) { x[2] })
  nw <- unique(nw)
  wts <- NULL
  for(j in nw) {
    xtmp <- x[nw0 == paste0("->", j)]
    wts <- cbind(wts, xtmp)
  }
  rownames(wts) <- NULL
  colnames(wts) <- nw
  return(wts)
}

