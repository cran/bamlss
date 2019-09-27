## Functions used for simulation.
simfun <- function(type = "sinus")
{
  ## Known function types.
  known_types <- c("linear", "quadratic", "unimodal", "double", "sinus",
    "cosinus", "pick", "complicated", "const", "spatial", "2d",
    "yue1", "yue2", "yue3")
  if(is.character(type))
    type <- match.arg(type, known_types)

  f <- switch(type,
    "linear" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      f <- -4.5 * x
      return(f - mean(f))
    },
    "quadratic" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      f <- 3.5 * (x - 0.5)^2
      return(f - mean(f))
    },
    "unimodal" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      f <- 120 * x * exp(-10 * x)
      return(f - mean(f))
    },
    "double" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      f <- 1.3 * (120 * x * exp(-10 * x) + 2.75 * x^2)
      return(f - mean(f))
    },
    "sinus" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      f <- sin(2 * pi * x)
      return(f - mean(f))
    },
    "cosinus" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      f <- cos(2 * pi * x)
      return(f - mean(f))
    },
    "pick" = function(x, min = 0, max = 1) { 
      x <- (x - min) / (max - min)
      f <- sin(2 * (4 * x - 2)) + 2 * exp(-16^2 * (x - 0.5)^2)
      return(f - mean(f))
    },
    "complicated" = function(x, min = 0, max = 1) {
      x <- (x - min) / (max - min)
      f <- exp(-400 * (x - 0.6)^2) + 5 * exp(-500 * (x - 0.75)^2) / 3 + 2 * exp(-500 * (x - 0.9)^2)
      return(f - mean(f))
    },
    "const" = function(..., const = 1.2) {
      const
    },
    "spatial" = function(id, f1 = sin, f2 = function(x) sin(x * pi * 2)) {
      n <- ceiling(sqrt(length(id)))
      co <- expand.grid("long" = seq(0, 1, length = n), "lat" = seq(0, 1, length = n))
      f <- f1(co[, 1]) * f2(co[, 2])
      f <- f - mean(f)
      f <- data.frame("long" = co[, 1], "lat" = co[, 2], "f" = f)
      f <- f[seq_along(id), ]
      return(f)
    },
    "2d" = function(x, y) {
      x <- scale2(x, -3, 3)
      y <- scale2(y, -3, 3)
      f <- sin(x) * cos(y)
      f <- f - mean(f)
      return(f)
    },
    "yue1" = function(x) {
      x <- scale2(x, 0, 1)
      knots <- c(0.2, 0.6, 0.7)
      B <- splineDesign(knots, x, ord = 3, outer.ok = TRUE)
      f <- B %*% c(20, 4, 6, 11, 6)
      f <- f - mean(f)
      return(f)
    },
    "yue2" = function(x) {
      x <- scale2(x, -2, 2)
      f <- sin(x) + 2 * exp(-30 * x^2)
      f <- f - mean(f)
      return(f)
    },
    "yue3" = function(x) {
      x <- scale2(x, 0, 1)
      e <- 0.125
      f <- sqrt(x * (1 - x)) * sin(2 * pi * (1 + e) / (x + e))
      f <- f - mean(f)
      return(f)
    }
  )
  if(!is.character(type))
    type <- known_types[type]
  attr(f, "type") <- type

  f
}

## Function for scaling.
scale2 <- function(x, lower = -1.5, upper = 1.5)
{
  x <- if(length(unique(x)) > 1) {
    (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)) * (upper - lower) + lower
  } else x
  x
}


## Artificial data set generating function.
GAMart <- function(n = 500, sd = 0.1, seed = FALSE, ti = c("none", "vcm", "main", "both"))
{
  if(seed) set.seed(111)

  ti <- match.arg(ti)

  n2 <- ceiling(sqrt(n))
  
  ## Other covariates.
  d <- data.frame(
    "x1" = runif(n, 0, 1),
    "x2" = runif(n, 0, 1),
    "x3" = runif(n, 0, 1),
    "fac" = factor(sample(1:3, size = n, replace = TRUE), labels = c("low", "medium", "high")),
    "lon" = sample(seq(0, 1, length = n2), size = n, replace = TRUE),
    "lat" = sample(seq(0, 1, length = n2), size = n, replace = TRUE)
  )

  i <- match.index(d[, c("lon", "lat")])
  d$id <- as.factor(i$match.index)

  ## Functions.
  ## Linear.
  f1 <- function(x) {
    scale2(-1.2 * x, -1, 1)
  }

  ## Doublemode.
  f2 <- function(x) {
    scale2(1.3 * (120 * x * exp(-10 * x) + 2.75 * x^2), -1, 1)
  }

  ## Quadratic. 
  f3 <- function(x) {
    scale2(3.5 * (x - mean(x))^2, -1, 1)
  }

  ## Spatial.
  f4 <- function(lon, lat) {
    scale2(sin(scale2(lon, -3, 3)) * cos(scale2(lat, -3, 3)), -1, 1)
  }

  ## Random.
  f5 <- function(id) {
    scale2(rnorm(length(unique(id)), sd = 0.2)[id], -0.2, 0.2)
  }

  ## Factor.
  f6 <- function(fac) {
    scale2(sort(rnorm(length(unique(fac)), sd = 1))[fac], -0.5, 0.5)
  }

  if(sd < 0) {
    sd <- exp(-2 + sin(d$x1 * 2 *pi - pi) - d$x2)
  }
  
  ## Response.
  if(ti == "none") {
    d$eta <- with(d, scale2(f1(x1) + f2(x2) + f3(x3) + f4(lon, lat) + f5(id) + f6(fac), -1, 1))
  } else {
    if(ti == "main")
      d$eta <- with(d, scale2(f1(lon) + f3(lat) + f4(lon, lat), -1, 1))
    if(ti == "both")
      d$eta <- with(d, scale2(0.7 * lon - 0.9 * lat + f1(lon) + f2(lat) + f4(lon, lat), -1, 1))
    if(ti == "vcm")
      d$eta <- with(d, 1.2 + 1.4 * x1 + x1 * f2(x2))
  }
  d$err <- rnorm(n, sd = sd)
  d$num <- with(d, eta + err)
  d$pnum <- d$num + abs(min(d$num)) + 1
  d$bnum <- scale2(d$num, 0.01, 0.99)
  d$cnum <- round(scale2(d$num, 0, 100))
  d$bin <- cut(d$num, quantile(d$num, probs = c(0, 0.5, 1)), labels = c("no", "yes"), include.lowest = TRUE)
  d$cat <- cut(d$num, quantile(d$num), labels = c("none", "low", "medium", "high"), include.lowest = TRUE)
  ystar <- d$eta + rnorm(n, sd = 1)
  d$cens <- ifelse(ystar > 0.0, ystar, 0.0)
  d <- d[, c("num", "pnum", "bnum", "cnum", "bin", "cat", "cens", "eta", "x1", "x2", "x3", "fac", "id", "lon", "lat", "err")]

  d
}


Volcano <- function(sd = 0.3) {
  eta <- scale2(t(get("volcano")), 0, 1)
  lon <- rev(scale2(1:nrow(eta), -36.877366, -36.877640))
  lat <- scale2(1:ncol(eta), 174.764712, 174.764879)
  d <- expand.grid("lon" = lon, "lat" = lat)
  d$y <- as.numeric(eta) + rnorm(length(eta), sd = sd)
  d
}


Crazy <- function(n = 1000) {
  d <- data.frame("x" = runif(n, -3, 3))
  d$eta <- sin(20 * exp(bamlss::scale2(d$x, 0, 1))) * bamlss::scale2(d$x, 0, 1)^2
  d$eta[d$x >= -1]  <- d$eta[d$x >= -1] - 1.5
  d$eta[d$x <= -2] <- - 1.5
  d$eta[d$x >= -2 & d$x <= -1] <- 0
  d$y <- d$eta + rnorm(n, sd = 0.1)
  return(d)
}

