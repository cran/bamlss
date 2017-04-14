## Plot map using longitude and latitude coordinates
## using the maps database.
xymap <- function(x, y, z, color = sequential_hcl(99, h = 100), raw.color = FALSE, symmetric = FALSE,
  swap = TRUE, p.cex = 0.01, pch = 22, legend = TRUE, add = FALSE, domar = TRUE,
  layout = TRUE, mmar = c(0, 0, 0, 0), lmar = c(1.3, 3, 1.3, 2), interp = FALSE,
  grid = 8, linear = FALSE, extrap = FALSE, duplicate = "mean", xlim = NULL,
  ylim = NULL, map = TRUE, boundary = TRUE, interior = TRUE, rivers = FALSE, mcol = NULL,
  contour.data = NULL, cgrid = 100, data = NULL, subset = NULL, box = FALSE,
  ireturn = FALSE, sort = TRUE, proj4string = CRS(as.character(NA)), eps = 0.0001, ...)
{
  ## projection = "+proj=longlat +ellps=WGS84 +datum=WGS84"

  if(!ireturn) {
    opar <- par(no.readonly = TRUE)
    on.exit(opar)
    if(domar & layout) {
      omar <- opar$mar
      par(mar = mmar)
    }
  }

  if(!is.null(data)) {
    subset <- deparse(substitute(subset), backtick = TRUE, width.cutoff = 500L)
    if(subset != "NULL") {
      subset <- eval(parse(text = subset), envir = data)
      data <- subset(data, subset)
    }
    x <- eval(parse(text = deparse(substitute(x), backtick = TRUE, width.cutoff = 500L)), envir = data)
    y <- eval(parse(text = deparse(substitute(y), backtick = TRUE, width.cutoff = 500L)), envir = data)
    z <- eval(parse(text = deparse(substitute(z), backtick = TRUE, width.cutoff = 500L)), envir = data)
  }

  data <- unique(data.frame("x" = as.numeric(x), "y" = as.numeric(y), "z" = as.numeric(z)))
  if(sort)
    data <- data[order(data$x), ]
  data <- na.omit(data)

  if(interp) {
    data <- na.omit(data)

    pp <- cbind(data$x, data$y)
    dx <- abs(diff(pp[, 1])); dy <- abs(diff(pp[, 2]))
    dx <- dx[abs(dx) > eps]
    dy <- dy[abs(dy) > eps]
    dx <- dx[dx != 0]; dy <- dy[dy != 0]
    res <- c(min(dx), min(dy)) / 2
    
    grid <- grid + 1

    px <- apply(pp, 1, function(x) {
      xs <- seq(x[1] - res[1], x[1] + res[1], length = grid)
      ys <- seq(x[2] - res[2], x[2] + res[2], length = grid)
      xs <- xs[-grid] + (xs[2] - xs[1]) / 2
      ys <- ys[-grid] + (ys[2] - ys[1]) / 2
      expand.grid("x" = xs, "y" = ys)
    })
    px <- do.call("rbind", px)

    data <- as.data.frame(mba.points(data, px, extend = TRUE, verbose = FALSE)$xyz.est)

    where <- maps::map.where("world", data$x, data$y)
    data <- data[!is.na(where), ]
    if(ireturn)
      return(data)
  }
  colors <- colorlegend(x = data$z, plot = FALSE, color = color,
    swap = swap, symmetric = symmetric, ...)
  col <- colors$map(data$z)
  p.cex <- if(is.null(p.cex)) {
    if(interp) 0.1 else 1
  } else p.cex

  coordinates(data) <- c("x", "y")
  proj4string(data) <- proj4string

  if(!add && legend && layout) {
    mar <- par()$mar
    par(mar = mar)
    w <- (4.5 + mar[4L]) * par("csi") * 3.7
  }
  if(!add && legend && layout)
    layout(matrix(c(1, 2), nrow = 1), widths = c(1, lcm(w)))

  pp <- coordinates(SpatialPoints(coordinates(data),
    proj4string = proj4string))
  dx <- abs(diff(pp[, 1])); dy <- abs(diff(pp[, 2]))
  dx <- dx[dx != 0]; dy <- dy[dy != 0]
  dx <- dx[abs(dx) > eps]
  dy <- dy[abs(dy) > eps]
  res <- c(min(dx), min(dy))

  if(!add) {
    plot(data, col = NA, bg = NA, xlim = xlim, ylim = ylim)
    if(box)
      box()
  }

  addmap <- NULL
  if(!is.logical(map)) {
    addmap <- map
    map <- FALSE
  }

  if(map) {
    m <- maps::map("world", add = TRUE, xlim = if(!add) NULL else xlim, ylim = if(!add) NULL else ylim,
      fill = if(is.null(mcol)) FALSE else TRUE, col = if(is.null(mcol)) gray(0.6) else mcol,
      boundary = boundary, interior = interior)
  }

  if(rivers) {
    stopifnot(requireNamespace("mapdata"))
    maps::map("rivers", add = TRUE, col = "lightblue")
  }

  rect(pp[, 1] - res[1] / 2, pp[, 2] - res[2] / 2, pp[, 1] + res[1] / 2, pp[, 2] + res[2] / 2,
    col = col, border = col, lwd = 0)

  if(rivers) {
    maps::map("rivers", add = TRUE, col = "lightblue")
  }

  if(!is.null(mcol))
    points(data, col = col, bg = col, pch = pch, cex = p.cex)
  if(!is.null(contour.data)) {
    if(is.logical(contour.data) & contour.data)
      contour.data <- data.frame("x" = as.numeric(x), "y" = as.numeric(y), "z" = as.numeric(z))
    contour.data <- unique(contour.data)
    x <- contour.data[, 1]; y <- contour.data[, 2]; z <- contour.data[, 3]
    xo <- seq(min(x), max(x), length = cgrid)
    yo <- seq(min(y), max(y), length = cgrid)
    fit <- interp2(x, y, z, xo = xo, yo = yo, grid = cgrid)

    if(!is.null(addmap)) {
      maptools::gpclibPermit()
      eg <- expand.grid("x" = xo, "y" = yo)
      nob <- length(slot(slot(addmap, "polygons")[[1]], "Polygons"))
      pip <- NULL
      for(j in 1:nob) {
        oco <- slot(slot(slot(addmap, "polygons")[[1]], "Polygons")[[j]], "coords")
        pip <- cbind(pip, point.in.polygon(eg$x, eg$y, oco[, 1L], oco[, 2L], mode.checked = FALSE))
      }
      pip <- apply(pip, 1, function(x) { any(x > 0) })
      fit[!pip] <- NA
      fit <- matrix(fit, cgrid, cgrid)
    }

    contour(xo, yo, fit, add = TRUE)
  }
  if(legend) {
    if(layout) {
      par(mar = lmar)
      colorlegend(x = data$z, full = TRUE,
        side.legend = 2, side.ticks = 2, color = color,
        swap = swap, symmetric = symmetric, ...)
    } else {
      colorlegend(x = data$z, plot = FALSE, add = TRUE,
        color = color, swap = swap, symmetric = symmetric, ...)
    }
  }
  
  invisible(colors)
}


## Function to compute a polygon map from
## a grid overlayed on x- and y-coordinates.
pixelmap <- function(x, y, size = 0.1, width = NULL, data = NULL,
  all = TRUE, at.xy = FALSE, yscale = TRUE, n = 20, ...)
{
  if(missing(x) & missing(y) & is.null(data)) {
    x <- expand.grid("x" = seq(0, 1, length = n), "y" = seq(0, 1, length = n))
    at.xy <- TRUE; size <- NA
  }
  if(missing(y) & (is.matrix(x) | is.data.frame(x) | is.list(x))) {
    x <- as.data.frame(x)
  } else {
    if(inherits(x, "formula")) {
      if(is.null(data))
        data <- environment(x)
      x <- model.frame(x, data = data)
      x <- x[, 2:1]
    } else {
      nx <- c(deparse(substitute(x), backtick = TRUE, width.cutoff = 500),
        deparse(substitute(y), backtick = TRUE, width.cutoff = 500))
      x <- data.frame(x, y)
      names(x) <- nx
    }
  }

  xr <- range(x[, 1], na.rm = TRUE)
  yr <- range(x[, 2], na.rm = TRUE)

  if(at.xy) {
    id <- rep(NA, nrow(x))
    xu <- unique(x)
    if(is.null(width)) {
      xd <- abs(diff(xu[, 1]))
      xd <- xd[xd > 0]
      xstep <- if(is.numeric(size)) {
        size * min(xd, na.rm = TRUE)
      } else {
        min(xd, na.rm = TRUE) / 2
      }
    } else xstep <- width / 2
    if(yscale) {
      p <- xstep / abs(diff(xr))
      ystep <- abs(diff(yr)) * p
    } else ystep <- xstep
    map <- list()
    n <- nrow(xu)
    for(j in 1:n) {
      map[[j]] <- cbind(
        "x" = c(xu[j, 1] - xstep, xu[j, 1] + xstep, xu[j, 1] + xstep, xu[j, 1] - xstep, x[j, 1] - xstep),
        "y" = c(xu[j, 2] - ystep, xu[j, 2] - ystep, xu[j, 2] + ystep, xu[j, 2] + ystep, x[j, 2] - ystep)
      )
      id[x[, 1] >= xu[j, 1] - xstep & x[, 1] <= xu[j, 1] + xstep & x[, 2] >= xu[j, 2] - ystep & x[, 2] <= xu[j, 2] + ystep] <- j
    }
    names(map) <- as.character(1:n)
  } else {
    xstep <- if(is.null(width)) abs(diff(xr)) * size else width / 2
    if(yscale) {
      p <- xstep / abs(diff(xr))
      ystep <- abs(diff(yr)) * p
    } else ystep <- xstep
    xstart <- xstart0 <- xr[1] - 0.5 * xstep
    ystart <- yr[1] - 0.5 * ystep
    map <- list(); k <- 1; id <- rep(NA, nrow(x))
    while(ystart < yr[2]) {
      xstart <- xstart0
      while(xstart < xr[2]) {
        map[[k]] <- cbind(
          "x" = c(xstart, xstart + xstep, xstart + xstep, xstart, xstart),
          "y" = c(ystart, ystart, ystart + ystep, ystart + ystep, ystart)
        )
        id[x[, 1] >= xstart & x[, 1] < xstart + xstep & x[, 2] >= ystart & x[, 2] < ystart + ystep] <- k
        xstart <- xstart + xstep
        k <- k + 1
      }
      ystart <- ystart + ystep
    }
    names(map) <- as.character(1:length(map))

    if(!all)
      map <- map[names(map) %in% as.character(id)]
  }

  class(map) <- c("bnd", "list")
  nmat <- neighbormatrix(map, ...)
  nn <- rowSums(nmat)
  nmat[nmat > 0] <- -1
  diag(nmat) <- nn
  id <- factor(as.character(id))
  lid <- levels(id)
  nmat <- nmat[lid, lid]

  return(list("map" = map, "nmat" = nmat, "data" = cbind(x, "id" = id)))
}


## Function to create a neighbormatrix from
## a list() of polygons or objects of class
## "SpatialPolygons".
neighbormatrix <- function(x, type = c("boundary", "dist", "delaunay", "knear"),
  k = 1, id = NULL, nb = FALSE, names = NULL, ...)
{
  if(inherits(x, "bnd") | inherits(x, "list")) {
    x <- list2sp(x)
  }

  type <- match.arg(type)

  adjmat <- if(!inherits(x, "nb")) {
    switch(type,
      "boundary" = spdep::poly2nb(x, ...),
      "dist" = spdep::dnearneigh(coordinates(x), ...),
      "delaunay" = spdep::tri2nb(coordinates(x), ...),
      "knear" = spdep::knn2nb(spdep::knearneigh(coordinates(x), k = k, ...), sym = TRUE))
  } else x

  if(!(spdep::is.symmetric.nb(adjmat, verbose = FALSE, force = TRUE))) {
    warning("neighbormatrix is not symmetric, will envorce symmetry!")
    adjmat <- spdep::make.sym.nb(adjmat)
  }

  if(!nb) {
    adjmat <- spdep::nb2mat(adjmat, style = "B", zero.policy = TRUE)
    if(inherits(x, "SpatialPolygonsDataFrame")) {
      names <- if(is.null(names)) {
        as.character(slot(x, "data")$OBJECTID)
      } else {
        try(as.character(slot(x, "data")[[names]]), silent = TRUE)
      }
    }
    if(!is.null(names) & class(names) != "try-error") {
      if(length(names) == nrow(adjmat)) {
        rownames(adjmat) <- names
        colnames(adjmat) <- names
      } else names <- rownames(adjmat)
    }

    if(!is.null(id)) {
      id <- if(is.factor(id)) {
        levels(id)
      } else {
        as.character(unique(id))
      }
      i <- rownames(adjmat) %in% id
      adjmat <- adjmat[i, i]
    }

    nn <- rowSums(adjmat)
    adjmat[adjmat > 0] <- -1
    diag(adjmat) <- nn

    colnames(adjmat) <- rownames(adjmat)
  }

  adjmat
}


## Function to plot the neighborhood relationship
## given some polygon list() map x.
plotneighbors <- function(x, add = FALSE, ...)
{
  nb <- neighbormatrix(x, nb = TRUE, ...)
  if(!add)
    plot(x, col = "lightgray")
  if(inherits(x, "bnd") | inherits(x, "list"))
    x <- list2sp(x)
  plot(nb, coordinates(x), add = TRUE, pch = 16)
  invisible(NULL)
}


### Function to create a spatial weight matrix
### of a polygon map.
#spatial.weights <- function(x, ...)
#{
#  requireNamespace("spdep")
#  nb <- neighbormatrix(x, nb = TRUE, ...)
#  weights <- listw2mat(nb2listw(nb, ...))
#  weights
#}


### Function to create spatial weight matrix
### from xy-coordinates
#spatial.weights2 <- function(x, y = NULL, d1 = 0, d2 = 0.5, W = FALSE, ...)
#{
#  if(!is.null(y))
#    x <- cbind(x, y)
#  if(!is.matrix(x))
#    x <- as.matrix(x)
#  weights <- nb2listw(dnearneigh(x, d1, d2), ...)
#  if(W) {
#    weights <- listw2mat(weights)
#  }
#  weights
#}


### Spatial weighted smooth constructor.
#smooth.construct.sws.smooth.spec <- function(object, data, knots) 
#{
#  requireNamespace("spdep")
#  xt <- object$xt
#  object$xt <- NULL
#  if(!is.null(xt$coords))
#    W <- spatial.weights2(xt$coords, W = FALSE)
#  if(!is.null(xt$weights))
#    W <- xt$weights
#  bs <- if(is.null(xt$bs)) "tp" else xt$bs
#  class(object) <- paste(bs, "smooth.spec", sep = ".")
#  object <- smooth.construct(object, data, knots)
#  object$X <- lag.listw(W, object$X,
#    zero.policy = if(is.null(xt$zero.policy)) TRUE else xt$zero.policy,
#    NAOK = TRUE)
#  object
#}


#smooth.construct.sws2.smooth.spec <- function(object, data, knots) 
#{
#  requireNamespace("spdep")
#  call <- if(length(object$term) > 1) {
#    paste("s(", paste(object$term[-1], collapse = ", "), ", bs = 'kr', k = ", object$bs.dim, ")", sep = "")
#  } else {
#    paste("s(", paste(object$term, collapse = ", "), ", bs = 'kr', k = ", object$bs.dim, ")", sep = "")
#  }
#  object$X <- smooth.construct(eval(parse(text = call)), data, NULL)$X
#  object$X <- object$X * data[[object$term[1]]]
#  object$S <- list()
#  object$rank <- 0
#  object$null.space.dim <- 0
#  object$fixed <- TRUE
#  object$plot.me <- FALSE
#  object
#}


#Predict.matrix.sws2.smooth <- function(object, data)
#{
#  call <- if(length(object$term) > 1) {
#    paste("s(", paste(object$term[-1], collapse = ", "), ", bs = 'kr', k = ", object$bs.dim, ")", sep = "")
#  } else {
#    paste("s(", paste(object$term, collapse = ", "), ", bs = 'kr', k = ", object$bs.dim, ")", sep = "")
#  }
#  object$X <- smooth.construct(eval(parse(text = call)), data, NULL)$X
#  object$X <- object$X * data[[object$term[1]]]
#  object
#}


## Compute centroids of polygons.
centroids <- function(x, id = NULL, verbose = FALSE, check.dups = TRUE)
{
  if(inherits(x, "SpatialPolygons")) {
    x <- BayesX::sp2bnd(x)
  }
  if(!is.list(x))
    stop("argument map must be a list() of matrix polygons!")

  n <- length(x)
  cp <- matrix(0, n, 2L)
  for(i in 1L:n) {
    cp[i,] <- centroidpos(na.omit(x[[i]]))
  }
  cp  <- as.data.frame(cp)
  nx <- nx0 <- names(x)
  if(any(i <- duplicated(nx))) {
    warning(paste("the following polygons are duplicated:", paste(nx[i], collapse = ", ")))
    for(j in nx[i]) {
      dups <- nx[nx == j]
      nx[nx == j] <- paste(dups, 1:length(dups), sep = ":")
    }
  }
  rownames(cp) <- nx
  colnames(cp) <- c("x", "y")

  if(!is.null(id)) {
    id <- as.character(unlist(id))
    cp2 <- matrix(NA, nrow = length(id), ncol = 2)
    for(j in unique(id)) {
      if(verbose) {
        cat("processing polygon:", j, "\n")
      }
      i <- which(nx0 == j)
      k <- which(id == j)
      pall <- list()
      take <- cp[i, ]
      for(l in 1:nrow(take))
        pall[[l]] <- as.numeric(take[l, ])
      pall <- rep(pall, length.out = length(k))
      pall <- do.call("rbind", pall)
      for(l in 1:length(k)) {
        cp2[k[l], ] <- pall[l, ]
      }
    }
    if(verbose) cat("creating data.frame\n")
    cp <- as.data.frame(cp2)
    if(check.dups) {
      if(any(i <- duplicated(id))) {
        for(j in id[i]) {
          if(verbose) cat("managing duplicates for region:", j, "\n")
          if(length(dups <- id[id == j]))
            id[id == j] <- paste(dups, 1:length(dups), sep = ":")
        }
      }
    }
    rownames(cp) <- id
    colnames(cp) <- c("x", "y")
  }

  return(cp)
}


centroidpos <- function(polygon) 
{
  polygon <- na.omit(polygon)
  p <- polygon
  np <- (nrow(p) - 1L)
  if(is.na(p[1L, 1L])) {
    p <- p[2L:(np + 1L),]
    np <- np - 1L
  }
  if((p[1L, 1L] != p[(np + 1L), 1L]) || (p[1L, 2L] != p[(np + 1L), 2L]))
    p[(np + 1L),] <- p[1L,]
	out <- cpos(p, np)

  return(out)
}

cpos <- function(p, np) 
{
  rval <- .Call("cpos",
    as.numeric(p),
    as.integer(np), PACKAGE = "bamlss")
  return(rval)
}

centroidtext <- function(polygon, poly.name = NULL, counter = "NA", cex = 1, ...) 
{
  pos <- centroidpos(polygon)		
  if(is.null(poly.name))
    txt <- paste(counter)
  else
    txt <- poly.name
  text(pos[1L], pos[2L], txt, cex = cex, ...)
  return(invisible(NULL))
}


## Function to drop data outside the polygon area.
drop2poly <- function(x, y, map, union = FALSE)
{
  if(inherits(map, "bnd") | inherits(map, "list")) {
    map <- list2sp(map)
  }
  if(union)
    map <- maptools::unionSpatialPolygons(map, rep(1L, length = length(map)), avoidGEOS  = TRUE)

  np <- length(map@polygons)
  pip <- NULL
  for(j in 1:np) {
    for(i in 1:length(map@polygons[[j]]@Polygons)) {
      oco <- map@polygons[[j]]@Polygons[[i]]@coords
      pip <- cbind(pip, point.in.polygon(x, y, oco[, 1L], oco[, 2L], mode.checked = FALSE))
    }
  }
  pip <- apply(pip, 1, any)

  return(pip)
}

xy2poly <- function(x, y, map, verbose = TRUE)
{
  id <- names(map)
  rval <- rep(NA, length(x))
  for(j in 1:length(map)) {
    if(verbose) {
      if(j > 1) cat("\r")
      cat("polygon", j)
    }
    tm <- map[j]
    class(tm) <- "bnd"
    i <- drop2poly(x, y, tm)
    rval[i] <- id[j]
  }
  if(verbose) cat("\n")
  rval
}
