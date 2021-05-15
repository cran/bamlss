plot2d <- function(x, residuals = FALSE, rug = FALSE, jitter = TRUE, 
  col.residuals = NULL, col.lines = NULL, col.polygons = NULL, 
  col.rug = NULL, c.select = NULL, fill.select = NULL, data = NULL,
  sep = "", month = NULL, year = NULL, step = 12,
  shift = NULL, trans = NULL, scheme = 2, s2.col = NULL, grid = 50, ...)
{
  rugp <- attr(x, "rug")
  if(is.null(x))
    return(invisible(NULL))
  if(is.character(x)) {
    stopifnot(file.exists(x <- path.expand(x)))
    x <- read.table(x, header = TRUE, sep = sep)
  }
  if(is.character(data)) {
    stopifnot(file.exists(data <- path.expand(data)))
    data <- read.table(data, header = TRUE, sep = sep)
  }
  if(inherits(x, "formula")) {
    if(is.null(data))
      data <- environment(x)
    else
      if(is.matrix(data))
        data <- as.data.frame(data)
    if(any(grep("+", as.character(x)[2L]))) {
      xch <- as.character(x)
      if(xch[2L] %in% names(data)) {
        if(inherits(data[[xch[2L]]], "data.frame")) {
          data[[xch[2L]]] <- as.matrix(data[[xch[2L]]])
        }
      }
      x <- try(model.frame(as.formula(paste("~", xch[2L])), data = data), silent = TRUE)
      if(inherits(x, "try-error")) {
        x <- model.frame(as.formula(paste0("~ as.matrix(", xch[2L], ")")), data = data)
      }
      x <- as.matrix(x)
      colnames(x) <- gsub(paste0(xch[2L], "."), "", colnames(x), fixed = TRUE)
      x <- cbind(model.frame(as.formula(paste("~", xch[3L])), data = data), x)
    } else x <- model.frame(x, data = data)
    if(ncol(x) < 2L)
      stop("formula is specified wrong!")
  }
  is.bayesx <- grepl(".bayesx", class(x))[1L]
  if(is.data.frame(x)) {
    if(!is.na(match("intnr", names(x))) & !is.null(c.select) & !is.character(c.select))
      c.select <- c.select - 1
    x <- df2m(x)
  }
  if(!is.list(x) && !is.matrix(x))
    stop("x must be a matrix!")
  if(!is.list(x) && ncol(x) < 2L)
    stop("x must have at least 2 columns!")
  if(ncol(x) < 3L)
    scheme <- 1
  if(ncol(x) > 10L)
    scheme <- 1
  if(scheme != 1) {
    if(is.null(col.lines)) {
      col.lines <- c(NA, "black", NA)
    } else {
      col.lines <- c(NA, col.lines[!is.na(col.lines)][1L], NA)
    }
  }
  args <- list(...)
  nc <- ncol(x)
  if(is.null(c.select)) {
    if(is.bayesx)
      c.select <- c(1L, 2L, 3L, 4L, 6L, 7L) 
    else 
      c.select <- 1L:nc
  }
  if(is.null(c.select))
    c.select <- 1L:nc
  if(length(c.select) > nc)
    c.select <- c.select[1L:nc]
  if(is.null(fill.select)) {
    if(is.bayesx)
      fill.select <- c(0L, 0L, 1L, 2L, 2L, 1L)
    if(all(c("Mean", "2.5%", "97.5%") %in% colnames(x)) & (ncol(x) == 4L))
      fill.select <- c(0, 1, 0, 1)
  }
  if(!is.bayesx && length(fill.select) < nc) {
    fill.select <- NULL
  }
  if(is.null(col.polygons))
    args$col.polygons <- rep(c("grey80", "grey70"), round(nc/2))
  else
    args$col.polygons <- col.polygons
  if(residuals && !is.null(pres <- attr(x, "residuals")))
    residuals <- TRUE
  else
    residuals <- FALSE
  by <- attr(x, "specs")$by
  if(is.null(by))
    by <- "NA"
  xnam <- attr(x, "specs")$term
  if(is.null(xnam))
    xnam <- colnames(x)[1L]
  if(is.null(xnam))
    xnam <- "x"
  if(by[1L] != "NA") {
    if(any(by == 0))
      x <- x[by != 0,]
    if(length(xnam) > 1L)	
      byname <- xnam[length(xnam)]
    else
      byname <- by
    xnam <- xnam[1L]
  }
  if(length(xnam) > 1L)
    xnam <- xnam[1L]
  if(is.null(args$xlab))
    args$xlab <- xnam
  if(is.null(args$ylab)) {
    if(is.null(attr(x, "specs")$label))
      args$ylab <- paste("Effect of", args$xlab)
    else
      args$ylab <- attr(x, "specs")$label
  }	
  if(is.character(c.select)) 
    c.select <- pmatch(c.select, colnames(x))
  x <- x[, c.select, drop = FALSE]
  if(!is.null(shift)) {
    shift <- as.numeric(shift[1])
    x[, 2:ncol(x)] <- x[, 2:ncol(x)] + shift
  }
  if(!is.null(trans)) {
    if(!is.function(trans)) stop("argument trans must be a function!")
    for(j in 2:ncol(x))
      x[, j] <- trans(x[, j])
  }
  if(residuals) {
    if(!is.null(shift)) pres[, 2L] <- pres[, 2L] + shift
    if(!is.null(trans)) pres[, 2L] <- trans(pres[, 2L])
    attr(x, "residuals") <- pres
  }
  if(is.null(args$ylim)) {
    ylim <- range(x[, -1], na.rm = TRUE)
    if(residuals)
      args$ylim <- range(c(ylim, pres[,2L]), na.rm = TRUE)
    else
      args$ylim <- range(ylim, na.rm = TRUE)
  }
  if(is.null(args$xlim))
    args$xlim <- base::range(x[,1L], na.rm = TRUE)
  if(!(!is.null(args$add) && args$add)) {
    graphics::plot(args$xlim, args$ylim, type = "n", axes = FALSE, 
      xlab = args$xlab, ylab = args$ylab, main = args$main)
  }
  args <- set.plot2d.specs(ncol(x) - 1L, args, col.lines, is.bayesx)
  if(!is.null(args$add)) {
    if(args$add) {
      if(is.null(args$axes))
        args$axes <- FALSE
    }
  }
  args$rugp <- rugp
  args$specs <- args
  args$residuals <- residuals
  args$col.residuals <- col.residuals
  args$col.rug <- col.rug
  args$fill.select <- fill.select
  args$pb <- FALSE
  args$rug <- rug
  args$jitter <- jitter
  args$x <- x
  args$specs$scheme <- scheme
  args$specs$s2.col <- s2.col
  args$specs$grid <- grid
  do.call(plot2d.default, delete.args(plot2d.default, args))
  if(is.null(args$type))
    box()
  else
    if(args$type != "n")
      box()
  if(is.null(args$axes)) {
    axis(2L, cex.axis = args$cex.axis)
    if(!is.null(month) & !is.null(year)) {
      start <- min(x[, 1], na.rm = TRUE) - month + 1
      stop <- max(x[, 1] + 1, na.rm=TRUE)
      pos <- seq(start, stop, step)
      label <- (pos - pos[1]) / step + year
      if(nrow(x) <= 24) {
        label2 <- month.abb[ifelse(step == 12, 1:12,
          ifelse(step == 4, c(1, 4, 7, 10),
          ifelse(step == 2, c(1, 7), FALSE)))]
        label2 <- rep(label2, length.out = nrow(x) + month - 1)
        label2 <- label2[month:(nrow(x) + month - 1)]
        start2 <- x[1, 1]
        stop2 <- max(x[, 1], na.rm = TRUE)
        pos2 <- seq(start2, stop2, 1)
        axis(side = 1, at = pos2, labels = label2, cex.axis = args$cex.axis)
      } else axis(side = 1, at = pos, labels = label, cex.axis = args$cex.axis)
    } else axis(1L, cex.axis = args$cex.axis)
  } else {
    if(args$axes) {
      axis(2L, cex.axis = args$cex.axis)
      axis(1L, cex.axis = args$cex.axis)
    }
  }

  return(invisible(NULL))
}


plot2d.default <- function(x, residuals, range, col.residuals = "black",
  fill.select = NULL, col.polygons = NULL, col.rug = NULL, pb = FALSE, 
  x.co = NULL, rug = FALSE, jitter = FALSE, specs)
{
  if(residuals && !is.null(pres <- attr(x, "residuals")))
    residuals <- TRUE
  else
    residuals <- FALSE
  if(nrow(x) > 1)
    x <- na.omit(x)
  if(!is.matrix(x))
    x <- matrix(x, nrow = 1L)
  if(residuals)
    e <- attr(x, "residuals")
  x <- unique(x)
  if(pb) {
    nc <- ncol(x)
    if(length(ux <- unique(x[,2L:nc])) < 3L) {
      fill.select <- NULL
      if(!is.matrix(ux))
        ux <- matrix(ux, nrow = 1L)
    } else ux <- matrix(unique(x[,2L:nc]), ncol = (nc - 1L))
    nux <- nrow(ux)
    if(nux < 2L) {
      nux <- 2L
      ux <- rbind(ux, ux)
    }
    x.co <- seq(x.co + range[1L], x.co - range[2L], length = nux)
    x <- cbind(x.co, ux)
    x <- rbind(x, x, x)
  }
  x <- x[order(x[,1L]), , drop = FALSE]
  if(!is.null(fill.select)) {      
    ufs <- unique(fill.select)
    ufs <- ufs[ufs != 0]
    nu <- length(ufs)
    if(!is.null(specs$poly.lty))
      specs$poly.lty <- rep(specs$poly.lty, length.out = nu)
    else
      specs$poly.lty <- rep(0, nu)
    if(is.null(specs$angle))
      specs$angle <- rep(45, nu)
    else
      specs$angle <- rep(specs$angle, length.out = nu)
    if(!is.null(specs$density))
      specs$density <- rep(specs$density, length.out = nu)
    else
      specs$density <- NULL
    if(!is.null(specs$border))
      specs$border <- rep(specs$border, length.out = nu)
    if(!is.null(specs$poly.lwd))
      specs$poly.lwd <- rep(specs$poly.lwd, length.out = nu)
    else
      specs$poly.lwd <- rep(1, nu)
    if(is.null(specs$scheme))
      specs$scheme <- 2
    for(k in 1L:nu) {
      check <- fill.select == ufs[k]
      if(length(check) == ncol(x)) {
        poly <- x[, check, drop = FALSE]
        if(specs$scheme == 1) {
          p1 <- poly[, 1L, drop = FALSE]
          p2 <- poly[, 2L, drop = FALSE]
          y.co <- c(p1, p2[length(p2):1L])
          x.co <- x[,1L]
          x.co <- c(x.co, x.co[length(x.co):1L])
          graphics::polygon(x = x.co, y = y.co, col = col.polygons[k], 
            lty = specs$poly.lty[k], border = specs$border[k], 
            density = specs$density[k], angle = specs$angle[k], 
            lwd = specs$poly.lwd[k])
        } else {
          grid <- if(is.null(specs$grid)) 30 else specs$grid
          mx <- grep("50%", colnames(x), fixed = TRUE)
          if(!length(mx))
            mx <- grep("mean", colnames(x), ignore.case = TRUE)
          if(!length(mx))
            mx <- grep("median", colnames(x), ignore.case = TRUE)
          if(!length(mx))
            mx <- grep("estimate", colnames(x), ignore.case = TRUE)
          if(length(mx)) {
            poly2 <- cbind(poly[, 1], x[, mx], poly[, 2])
            poly <- apply(poly2, 1, function(x) {
              x <- as.numeric(x)
              c(seq(x[1], x[2], length = ceiling(grid / 2)), seq(x[2], x[3], length = ceiling(grid / 2)))
            })
          } else {
            poly <- apply(poly, 1, function(x) { x <- as.numeric(x); seq(x[1], x[2], length = grid) })
          }
          if(is.null(specs$s2.col))
            specs$s2.col <- rev(gray.colors(grid, start = 0, end = 1, gamma = 1))
          if(is.function(specs$s2.col))
            specs$s2.col <- specs$s2.col(grid)
          specs$s2.col <- rep(specs$s2.col, length.out = grid / 2)
          for(pj in 1:(grid / 2)) {
            p1 <- poly[pj,]
            p2 <- poly[grid - pj + 1, ]
            y.co <- c(p1, p2[length(p2):1L])
            x.co <- x[,1L]
            x.co <- c(x.co, x.co[length(x.co):1L])
            graphics::polygon(x = x.co, y = y.co, col = specs$s2.col[pj], 
              lty = specs$poly.lty[k], border = specs$border[k], 
              density = specs$density[k], angle = specs$angle[k], 
              lwd = specs$poly.lwd[k])
          }
        }
      }
    }
  }    
  if(residuals) {
    pargs <- list()
    pargs$x <- pres[,1L]
    pargs$y <- pres[,2L]
    pargs$cex <- specs$cex
    pargs$type <- specs$type
    pargs$pch <- specs$pch
    pargs$col <- col.residuals
    do.call(graphics::points, pargs)
  }
  for(k in 2L:ncol(x)) {
    lines(x[,k] ~ x[,1L], lty = specs$lty[k - 1L], lwd = specs$lwd[k - 1L], 
      col = specs$col.lines[k - 1L])
  }
  if(rug) {
    specs$col <- col.rug
    rugp <- if(!is.null(specs$rugp)) specs$rugp else x[,1L]
    if(jitter)      
      specs$x <- jitter(rugp)
    else
      specs$x <- rugp
    do.call(graphics::rug, delete.args(graphics::rug, specs))
  }

  return(invisible(NULL))
}


plot3d <- function(x, residuals = FALSE, col.surface = NULL, 
  ncol = 99L, swap = FALSE, col.residuals = NULL, col.contour = NULL, 
  c.select = NULL, grid = 30L, image = FALSE, contour = FALSE, 
  legend = TRUE, cex.legend = 1, breaks = NULL, range = NULL, 
  digits = 2L, d.persp = 1L, r.persp = sqrt(3), 
  outscale = 0, data = NULL, sep = "",
  shift = NULL, trans = NULL,
  type = "mba", linear = FALSE, extrap = FALSE, k = 40, ...)
{
  if(is.null(x))
    return(invisible(NULL))
  if(is.character(x)) {
    stopifnot(file.exists(x <- path.expand(x)))
    x <- read.table(x, header = TRUE, sep = sep)
  }
  if(is.character(data)) {
    stopifnot(file.exists(data <- path.expand(data)))
    data <- read.table(data, header = TRUE, sep = sep)
  }
  if(inherits(x,"formula")) {
    if(is.null(data))
      data <- environment(x)
    else
      if(is.matrix(data))
        data <- as.data.frame(data)
    x <- model.frame(x, data = data)
    if(ncol(x) < 3L)
      stop("formula is specified wrong!")
    if(ncol(x) > 3L)
      x <- x[, c(2L, 3L, 1L, 4L:ncol(x))]
    else
      x <- x[, c(2L, 3L, 1L)]
  }
  if(is.data.frame(x)) {
    if(!is.na(match("intnr", names(x))) & !is.null(c.select) & !is.character(c.select))
      c.select <- c.select - 1
    x <- df2m(x)
  }
  if(!is.matrix(x))
    stop("x must be a matrix!")
  if(ncol(x) < 3)
    stop("x must have at least 3 columns!")
  args <- list(...)
  if(!is.null(shift))
    shift <- as.numeric(shift[1])
  e <- NULL
  if(!is.null(attr(x, "residuals"))) {
    e <- attr(x, "residuals")
    if(!is.null(shift))
      e[, 3L] <- e[, 3L] + shift
  }
  if(!is.null(e) && all(is.na(e)))
    residuals <- FALSE
  specs <- attr(x, "specs")
  by <- specs$by
  if(is.null(by))
    by <- "NA"
  else {
    if(!is.null(specs)  && length(specs$term) > 2L)
      by <- specs$term[length(specs$term)]
  }
  nx <- colnames(x)
  x <- x[order(x[, 1L]), ]
  X <- x[, 1L]
  z <- x[, 2L]
  xrd <- diff(range(X))
  zrd <- diff(range(z))
  xn <- seq(min(X) - outscale * xrd , max(X) + outscale * xrd, length = grid)
  zn <- seq(min(z) - outscale * zrd, max(z) + outscale * zrd, length = grid)
  fitted <- list(NA)
  if(!is.null(c.select)) {
    take <- NULL
    id <- 1L:length(nx)
    if(length(c.select) < 2L && c.select == 95) 
      c.select <- as.character(c.select)
    if(length(c.select) < 2L && c.select == 80)
      c.select <- as.character(c.select)
    is.se <- FALSE
    if(!is.na(pmatch("95", c.select))) {
      take <- id[nx %in% c("2.5%", "97.5%")]
      is.se <- TRUE
    }
    if(!is.na(pmatch("80", c.select))) {
      take <- id[nx %in% c("10%", "90%")]
      is.se <- TRUE
    }
    if(is.se) {
      take2 <- c("mean", "Mean", "MEAN", "estimate", 
        "Estimate", "ESTIMATE", "mean", "pmode")
      for(k in take2)
        if(any(nx %in% k))
          take <- c(take[1], id[nx %in% k][1], take[2])
    }
    if(!length(take)) take <- NULL
    if(is.null(take) && !is.character(c.select)) {
      if(min(c.select) < 3L)
        c.select <- c.select + 2L
      if(max(c.select) > ncol(x) || min(c.select) < 3L)
        stop("argument c.select is specified wrong!")
      take <- unique(c.select)
    }
    if(is.null(take) && is.character(c.select))
      for(k in c.select)
        for(i in 1L:length(nx))
          if(!is.na(pmatch(k, nx[i])))
            take <- c(take, i)
    if(is.null(take))
      stop("argument c.select is specified wrong!")
    for(k in 1:length(take)) {
      fitted[[k]] <- interp2(X, z, x[, take[k]], xo = xn, yo = zn,
        type = type, linear = linear, extrap = extrap, k = k)
    }
  }
  if(length(fitted[[1L]]) == 1L && is.na(fitted[[1L]][1L])) {
    fitted[[1L]] <- interp2(X, z, x[, 3L], xo = xn, yo = zn,
      type = type, linear = linear, extrap = extrap, k = k)
  }
  if(!is.null(range)) {
    for(k in 1L:length(fitted)) {
      if(min(range, na.rm = TRUE) > min(fitted[[k]], na.rm = TRUE))
        fitted[[k]][fitted[[k]] < min(range, na.rm = TRUE)] <- min(range, na.rm = TRUE)  
      if(max(range, na.rm = TRUE) < max(fitted[[k]], na.rm = TRUE))
        fitted[[k]][fitted[[k]] > max(range, na.rm = TRUE)] <- max(range, na.rm = TRUE)  
    }
  }
  if(!is.null(shift)) {
    for(k in 1L:length(fitted)) {
        fitted[[k]] <- fitted[[k]] + shift
    }
  }
  if(!is.null(trans)) {
    if(!is.function(trans)) stop("argument trans must be a function!")
    for(k in 1L:length(fitted)) {
      fitted[[k]] <- trans(fitted[[k]])
    }
  }
  names <- colnames(x)[1L:2L]
  if(residuals)
    zlimit <- range(c(unlist(fitted), e[, 3L]), na.rm = TRUE)
  else
    zlimit <- range(unlist(fitted), na.rm = TRUE)
  if(is.null(args$xlab))
    args$xlab <- names[1L]
  if(is.null(args$ylab))
    args$ylab <- names[2L]
  if(is.null(args$zlab)) {
    if(!is.null(specs) && is.null(specs$label))
      args$zlab <- "fitted"
    else
      args$zlab <- specs$label
  }
  if(is.null(args$zlab))
    args$zlab <- try(paste("f(", nx[1L], ",", nx[2L], ")", sep = ""))
  args$y <- substitute(zn)
  args$x <- substitute(xn)
  if(is.null(col.surface))
    col.surface <- colorspace::diverge_hcl
  if(!is.null(args$add) && args$add)
    par(new = TRUE)
  pmat0 <- NULL
  if(!image && !contour) {
    myfit <- matrix(fitted[[1L]], grid, grid)
    if(length(fitted) < 2L) {
      args$col <- make_pal(col = col.surface, ncol = ncol, data = myfit, 
        range = range, breaks = breaks, swap = swap, 
        symmetric = args$symmetric)$map(myfit)
    } else args$col <- col.surface
    args$z <- substitute(myfit)
    args$d <- d.persp
    args$r <- r.persp
    if(is.null(args$zlim))
      args$zlim <- zlimit
    args$zlim <- range(args$zlim)
    if(identical(args$zlim[1], args$zlim[2])) {
      args$zlim[1] <- args$zlim[1] - 0.1
      args$zlim[2] <- args$zlim[2] + 0.1
    }
    if(is.null(args$theta))
      args$theta <- 40
    if(is.null(args$phi))
      args$phi <- 40
    if(!is.null(c.select) && length(fitted) > 1L) {
      nf <- length(fitted)
      if(is.null(args$border))
        args$border <- c(3, 1, 2)
      if(is.function(args$col) || is.null(args$col))
        args$col <- NA
      color <- rep(args$col, length.out = nf)
      bcol <- rep(args$border, length.out = nf)
      args$col <- color[1L]
      args$border <- bcol[1L]
      pmat <- pmat0 <- do.call(graphics::persp, 
        delete.args("persp.default", args, c("lwd", "lty"), package = "graphics"))
      for(k in 2L:length(fitted)) {
        par(new = TRUE)
        args$col <- color[k]
        args$border <- bcol[k]
        if(k > 1) {
          args$xlab <- args$ylab <- args$zlab <- args$main <- NA
          args$box <- FALSE
          args$axes <- FALSE
        }
        myfit <- matrix(fitted[[k]], grid, grid)
        args$z <- substitute(myfit)
        pmat <- do.call(graphics::persp, 
          delete.args("persp.default", args, c("lwd", "lty"), package = "graphics"))
      }
    } else {
      pmat <- pmat0 <- do.call(graphics::persp, delete.args("persp.default",
        args, c("lwd", "lty"), package = "graphics"))
    }
    if(residuals && !is.null(e)) {
      t3d <- trans3d(e[,1L], e[,2L], e[,3L], pmat)
      if(is.null(col.residuals))
        col.residuals <- "black"
      points(x = t3d$x, y = t3d$y, cex = args$cex, pch = args$pch, col = col.residuals)
    }
  }
  if(image || contour) {
    myfit <- matrix(fitted[[1L]], grid[1L], grid[1L])
    pal <- make_pal(col = col.surface, ncol = ncol, data = myfit, 
      range = range, breaks = breaks, swap = swap, 
      symmetric = args$symmetric)
    args$col <- pal$colors
    args$breaks <- pal$breaks
    if(is.null(args$xlim))
      args$xlim <- range(xn)
    if(is.null(args$ylim))
      args$ylim <- range(zn)
    add <- FALSE
    args$z <- substitute(myfit)
    args$x <- xn
    args$y <- zn
    args$zlab <- NULL
    if(image) {
      if(legend & is.null(args$pos)) {
        mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
        mar <- mar.orig
        on.exit(par(par.orig))
        mar[4L] <- 0
        par(mar = mar)
        w <- (3 + mar[2L]) * par("csi") * 2.54
        layout(matrix(c(1, 2), nrow = 1), widths = c(1, lcm(w)))
      }
      do.call(graphics::image, 
        delete.args(graphics::image.default, args, 
        c("xlab", "ylab", "main", "axes")))
      if(!is.null(args$image.map)) {
          args2 <- args
          args2$map <- args$image.map
          args2$add <- TRUE
          args2$legend <- FALSE
          args2$x <- NULL
          args2$id <- NULL
          args2$col <- NULL
          do.call(plotmap, delete.args(plotmap, args2))
      }
      if(contour) {
        if(is.null(col.contour)) 
          args$col <- "black"
        else
          args$col <- col.contour
        args$add <- TRUE
        do.call(graphics::contour.default, 
          delete.args(graphics::contour.default, args, 
          c("xlab", "ylab", "main", "axes")))
        contour <- FALSE
      }
      if(legend) {
        if(is.null(args$pos)) {
          mar <- mar.orig
          mar[2L] <- 1
          mar[4L] <- 3.1
          par(mar = mar, xaxs = "i", yaxs = "i")
        }
        args2 <- args
        if(is.null(args$side.legend))
          args2$side.legend <- 2L
        if(is.null(args2$pos)) {
          if(is.null(args$distance.labels))
            args2$distance.labels <- 3L
          if(is.null(args$side.ticks))
            args2$side.ticks <- 2L
        }
        args2[["col"]] <- NULL
        args2$color <- col.surface
        args2$ncol <- ncol
        args2$x <- args$z
        args2$xlim <- range(xn)
        args2$ylim <- range(zn)
        args2$breaks <- breaks
        args2$swap <- swap
        args2$plot <- TRUE
        args2$digits <- digits
        args2$cex <- cex.legend
        args2$range <- range
        if(is.null(args$pos)) {
          args2$add <- FALSE
          args2$full <- TRUE
        } else {
          args2$add <- TRUE
          args2$full <- FALSE
          args2$plot <- FALSE
        }
        do.call(colorlegend, delete.args(colorlegend, args2, "font"))$map
      }
    }
    if(contour) {
      if(is.null(col.contour)) 
        args$col <- "black"
      else
        args$col <- col.contour
      do.call(graphics::contour.default, 
        delete.args(graphics::contour.default, args, 
        c("xlab", "ylab", "main", "axes")))
    }
    args$pal <- pal
  }
  args$pmat <- pmat0

  return(invisible(args))
}


plotblock <- function(x, residuals = FALSE, range = c(0.3, 0.3), 
  col.residuals = "black", col.lines = "black", c.select = NULL, 
  fill.select = NULL , col.polygons = NULL, data = NULL,
  shift = NULL, trans = NULL, labels = NULL, ...)
{
  if(is.null(x))
    return(invisible(NULL))
  if(inherits(x, "formula")) {
    if(is.null(data))
      data <- environment(x)
    else
      if(is.matrix(data))
        data <- as.data.frame(data)
    if(any(grep("+", as.character(x)[2]))) {
      xch <- as.character(x)
      x <- model.frame(as.formula(paste("~", xch[2L])), data = data)
      x <- cbind(model.frame(as.formula(paste("~", xch[3L])), data = data), x)
    } else x <- model.frame(x, data = data)
    if(ncol(x) < 2L)
      stop("formula is specified wrong!")
  }
  is.bayesx <- grepl(".bayesx", class(x))[1L]
  if(is.data.frame(x)) {
    labels <- as.character(x[, 1])
    x <- cbind(as.integer(x[, 1]), as.matrix(x[, -1]))
  }
  if(!is.list(x) && !is.matrix(x))
    stop("x must be a matrix!")
  if(!is.list(x) && ncol(x) < 2L)
    stop("x must have at least 2 columns!")
  args <- list(...)
  if(is.null(args$xlab))
    args$xlab <- attr(x, "specs")$term
  if(is.null(args$ylab)) {
    if(is.null(attr(x, "specs")$label))
      args$ylab <- paste("f(", args$xlab, ")", sep = "")
    else
      args$ylab <- attr(x, "specs")$label
  }
  if(!is.null(shift))
    shift <- as.numeric(shift[1])
  if(!is.list(x))
    nc <- ncol(x)
  else
    nc <- ncol(x[[1L]])
  if(is.null(c.select)) {
    if(is.bayesx)
      c.select <- c(1L, 2L, 3L, 4L, 6L, 7L) 
    else
      c.select <- 1L:nc
  }
  if(is.null(c.select))
    c.select <- 1L:nc
  if(length(c.select) > nc)
    c.select <- c.select[1L:nc]
  if(is.null(fill.select))
    if(is.bayesx)
      fill.select <- c(0L, 0L, 1L, 2L, 2L, 1L) 
  if(!is.bayesx && length(fill.select) < nc)
    fill.select <- NULL
  xnam <- attr(x, "specs")$term
  if(is.null(xnam))
    xnam <- "x"  
  partial.resids <- NULL
  if(is.null(range)) {
    dow <- 0.3
    up <- 0.3
  } else {
    dow <- range[1L]
    up <- range[2L]	
    if(is.na(dow))
      dow <- 0.3
    if(is.na(up))
      up <- 0.3
  }
  ylim <- NULL
  if(!is.list(x)) {
    if(is.null(e <- attr(x, "residuals")))
      residuals <- FALSE
    xu <- unique(x[,1L])
    n <- length(xu)      
    effects <- vector("list", n)
    for(i in 1L:n) {
      effects[[i]] <- x[x[,1L] == xu[i], c.select, drop = FALSE]
      if(!is.matrix(effects[[i]]))
        effects[[i]] <- matrix(effects[[i]], nrow = 1L)
      if(!is.matrix(effects[[i]]))
        effects[[i]] <- matrix(effects[[i]], nrow = 1L)
      if(!is.null(shift)) effects[[i]][, 2L:ncol(x[[i]])] <- effects[[i]][, 2L:ncol(x[[i]])] + shift
      if(!is.null(trans)) {
        if(!is.function(trans)) stop("argument trans must be a function!")
        for(j in 2:ncol(effects[[i]]))
          effects[[i]][, j] <- trans(effects[[i]][, j])
      }
      ylim <- c(ylim, effects[[i]][, 2L:ncol(effects[[i]])])
      colnames(effects[[i]]) <- rep(paste(xnam, xu[i], sep = ""), ncol(effects[[i]]))
      if(residuals) {
        if(length(pres <- e[e[,1L] == xu[i],])) {
          if(is.null(dim(pres)))
            pres <- matrix(pres, nrow = 1)
          if(!is.null(shift))
            pres[, 2L:ncol(pres)] <- pres[, 2L:ncol(pres)] + shift
          attr(effects[[i]], "residuals") <- pres
          ylim <- c(ylim, pres[, 2L:ncol(pres)])
        }
      }
    }
    x <- effects
  } else {
    n <- length(x)	
    for(i in 1L:n) {
      if(residuals && !is.null(pres <- attr(x[[i]], "residuals"))) {
        pres <- pres[pres[,1L] != 0 & pres[,1L] != -1,]
        if(!is.matrix(pres))
          pres <- matrix(pres, nrow = 1L)
        pres[,1L] <- i
        ylim <- c(ylim, pres[,2L:ncol(pres)])
      }
      if(is.data.frame(x[[i]]))
        x[[i]] <- df2m(x[[i]])
      cn <- colnames(x[[i]])
      cn <- cn[c.select]
      x[[i]] <- x[[i]][,c.select]
      if(!is.matrix(x[[i]]))
        x[[i]] <- matrix(x[[i]], nrow = 1L)
      x[[i]] <- x[[i]][x[[i]][,1L] != 0 & x[[i]][,1L] != -1,]
      if(!is.matrix(x[[i]]))
        x[[i]] <- matrix(x[[i]], nrow = 1L)
      if(nrow(x[[i]]) == 2L && x[[i]][1L, 1L] == -1)
        x[[i]] <- matrix(x[[i]][2L,], nrow = 1L)
      colnames(x[[i]]) <- cn
      if(residuals) {
        if(!is.null(shift)) {
          if(is.matrix(pres))
            pres[, 2L:ncol(pres)] <- pres[, 2L:ncol(pres)] + shift
          else
            pres <- pres + shift
        }
        attr(x[[i]], "residuals") <- pres
      }
      if(!is.null(shift)) x[[i]][, 2L:ncol(x[[i]])] <- x[[i]][, 2L:ncol(x[[i]])] + shift
      if(!is.null(trans)) {
        if(!is.function(trans)) stop("argument trans must be a function!")
        for(j in 2:ncol(x[[i]]))
          x[[i]][, j] <- trans(x[[i]][, j])
      }
      ylim <- c(ylim, x[[i]][,2L:ncol(x[[i]])])
    }
  }
  if(!is.null(labels))
    labels <- rep(labels, length.out = length(x))
  if(is.null(args$xlim))
    args$xlim <- c(0.5, n + 0.5)
  if(is.null(args$ylim))
    args$ylim <- base::range(ylim, na.rm = TRUE)
  if(is.null(args$xlab))
    args$xlab <- xnam	
  if(!is.null(args$add) && args$add)
    par(new = TRUE)
  graphics::plot(args$xlim, args$ylim, type = "n", axes = FALSE, 
    xlab = args$xlab, ylab = args$ylab, main = args$main)
  args <- set.plot2d.specs(ncol(x[[1L]]) - 1L, args, col.lines, is.bayesx)
  xnames <- NULL
  axn <- rep(NA, n)
  args$specs <- args
  args$residuals <- residuals
  args$range <- range
  args$col.residuals <- col.residuals
  args$fill.select <- fill.select
  if(is.null(col.polygons))
    args$col.polygons <- rep(c("grey70", "grey50"), round(nc / 2))
  else
    args$col.polygons <- col.polygons
  args$pb <- TRUE
  for(i in 1L:n) {
    args$x.co <- i
    if(length(unique(x[[i]][, -1])) < 2) {
      xvals <- sort(jitter(x[[i]][, -1], factor = .Machine$double.eps^0.5))
      x[[i]][, 2:ncol(x[[i]])] <- xvals
    }
    args$x <- x[[i]]
    if(!is.null(attr(args$x, "residuals")))
      attr(args$x, "residuals")[,1L] <- i
    do.call(plot2d.default, delete.args(plot2d.default, args))
    axn[i] <- if(is.null(labels)) colnames(x[[i]])[1L] else labels[i]
  }
  if(is.null(args$type))
    box()
  else
    if(args$type != "n")
      box()
  if(is.null(args$axes)) {
    axis(2L)
    axis(1L, at = 1L:n, labels = axn)
  } else
    if(args$axes) {
      axis(2L)
      axis(1L, at = 1L:n, labels = axn)
    }

  return(invisible(NULL))
}


colorlegend <- function(color = NULL, ncol = NULL, x = NULL, breaks = NULL, 
  pos = "center", shift = 0.02, side.legend = 1L, side.ticks = 1L, range = NULL, lrange = NULL, 
  width = 0.25, height = 0.05, scale = TRUE, xlim = NULL, ylim = NULL, plot = NULL, full = FALSE,
  add = FALSE, col.border = "black", lty.border = 1L, lwd.border = 1L, ticks = TRUE, 
  at = NULL, col.ticks = "black", lwd.ticks = 1L, lty.ticks = 1L, length.ticks = 0.3, 
  labels = NULL, distance.labels = 0, col.labels = "black", cex.labels = 1L, 
  digits = 2L, swap = FALSE, symmetric = TRUE, xpd = NULL,
  title = NULL, side.title = 2, shift.title = c(0, 0), cex.title = 1, ...)
{
  args <- list(...)
  if(is.null(xlim)) {
    if(add)
      xlim <- par("usr")[1L:2L]
    else
      xlim <- c(0L, 1L)
  }
  if(is.null(ylim)) {
    if(add)
      ylim <- par("usr")[3L:4L]
    else
      ylim <- c(0L, 1L)
  }
  if(!side.legend %in% c(1L, 2L)) {
    warning("argument side.legend is specified wrong, set to default!")
    side.legend <- 1L
  }
  if(full) {
    scale <- FALSE
    pos = c(0, 0)
    par(xaxs = "i")
    par(yaxs = "i")
    if(side.legend < 2L) {
      width <- diff(xlim)
      height <- diff(ylim)
    } else {
      height <- diff(xlim)
      width <- diff(ylim)
    }
  }
  if(is.null(plot) || plot == TRUE) {
    plot <- TRUE
    graphics::plot.default(xlim, ylim, type = "n", xlab = "", ylab = "", 
      axes = FALSE, xlim = xlim, ylim = ylim, asp = NA)
  } else plot <- FALSE
  if(is.null(xpd))
    xpd <- FALSE
  if(xpd)
    par(xpd = xpd)
  pos2 <- NULL
  postxt <- c("bottomleft", "topleft", "topright", "bottomright",
    "left", "right", "top", "bottom", "center")
  poscheck <- pmatch(pos, postxt)
  if(all(!is.na(poscheck)) && length(poscheck) > 0) {
    pos2 <- postxt[pmatch(pos, postxt)]
    pos <- c(0, 0)
  }
  if(is.null(pos)) {
    pos <- c(0.35, 0.15)
    if(side.legend < 2L)
      pos <- rev(pos)   
  }
  limits <- list(xlim, ylim)
  pos <- opos <- c(min(limits[[1L]], na.rm = TRUE) + pos[1L] * diff(limits[[1L]]), 
    min(limits[[2L]], na.rm = TRUE) + pos[2L] * diff(limits[[2L]])) 
  if(side.legend > 1L)
    limits <- rev(limits)
  if(scale) {
    width <- width * diff(limits[[1L]])
    height <- height * diff(limits[[2L]])
  }
  if(side.legend > 1L) {
    wi <- width
    width <- height
    height <- wi
  }
  if(full)
    shift <- 0
  shift <- rep(shift, length.out = 2)
  if(is.null(pos2)) {
    xlim <- range(c(pos[1L], pos[1L] + width, pos[1L] + width, pos[1L])) + shift[1] * width
    ylim <- range(c(pos[2L], pos[2L], pos[2L] + height, pos[2L] + height)) + shift[2] * height
  } else {
    pos2 <- dopos(pos2, limits, width, height, side.legend, shift)
    xlim <- pos2$xlim
    ylim <- pos2$ylim
  }
  if(!is.null(x)) {
    if(is.null(lrange)) {
      if(is.null(range)) {      
        lrange <- range(x, na.rm = TRUE)
        if(symmetric) {
          mar <- max(abs(lrange))
          lrange <- c(0 - mar, mar)
        }
      } else lrange <- range
    }
    x <- unique(na.omit(sort(x)))
  } else { 
    if(is.null(range))
      range <- xlim
    if(is.null(lrange))
      lrange <- range
  }
  if(is.null(color))
    color <- grDevices::gray.colors
  args$col <- color
  args$ncol <- ncol
  args$data <- x
  args$range <- range
  args$breaks <- breaks
  args$swap <- swap
  args$symmetric <- symmetric
  pal <- do.call(make_pal, delete.args(make_pal, args))
  if(plot || add) {
    if(!is.null(lrange)) {
      if(min(lrange) > min(pal$breaks)) 
        pal$breaks[pal$breaks <= min(lrange)] <- min(lrange)
      if(max(lrange) < max(pal$breaks))
        pal$breaks[pal$breaks >= max(lrange)] <- max(lrange)
    }
    br <- c(min(pal$breaks, lrange), pal$breaks, max(pal$breaks, lrange))
    cl <- c(head(pal$colors, 1L), pal$colors, tail(pal$colors, 1L))
    obs2legend <- function(x, xr) ((x - lrange[1L]) / diff(lrange)) * diff(xr) + xr[1L]
    if(side.legend < 2L) {
      graphics::rect(obs2legend(head(br, -1L), xlim), ylim[1L], obs2legend(tail(br, -1L), xlim),
        ylim[2L], col = cl, border = cl, xpd = xpd, lwd = 0.01)
    } else {
      graphics::rect(xlim[1L], obs2legend(head(br, -1L), ylim), xlim[2L], 
        obs2legend(tail(br, -1L), ylim), col = cl, border = cl, xpd = xpd, lwd = 0.01)
    }
    graphics::rect(xlim[1L], ylim[1L], xlim[2L], ylim[2L], 
      border = col.border, lwd = lwd.border, lty = lty.border, xpd = xpd)
    dl <- TRUE
    if(!is.null(labels)) {
      if(is.logical(labels)) {
        if(!labels)
          dl <- FALSE
      }
    }
    if(ticks || dl) {
      if(is.null(at)) {
        at <- pal$breaks
        if(abs(diff(lrange / max(lrange))) / length(at) < 0.2)
          at <- seq(min(lrange), max(lrange), length.out = 3L)
      }
      if(is.null(labels)) {
        justify <- "centre"
        if((side.legend == 2) & (side.ticks == 1))
          justify <- "right"
        if((side.legend == 2) & (side.ticks == 2))
          justify <- "left"
        labels <- format(gsub(" ", "", format(at, digits = digits, nsmall = digits)), justify = justify)
      }
      if(side.legend < 2L) {
        at <- obs2legend(at, xlim)
        length.ticks <- length.ticks * height
        if(any(at > max(xlim))) 
          at[at > max(xlim)] <- max(xlim)
        if(any(at < min(xlim)))
          at[at < min(xlim)] <- min(xlim)
      } else {
        at <- obs2legend(at, ylim)
        length.ticks <- length.ticks * width
        if(any(at > max(ylim))) 
          at[at > max(ylim)] <- max(ylim)
        if(any(at < min(ylim)))
          at[at < min(ylim)] <- min(ylim)
      }
      at <- unique(at)
      if(side.ticks > 1L)
        length.ticks <- (-1) * length.ticks
      nat <- length(at)
      lwd.ticks <- rep(lwd.ticks, length.out = nat)
      lty.ticks <- rep(lty.ticks, length.out = nat)
      col.ticks <- rep(col.ticks, length.out = nat)
      col.labels <- rep(col.labels, length.out = nat)
      cex.labels <- rep(cex.labels, length.out = nat)
      labels <- rep(labels, length.out = nat)
      if(!full) {
        for(i in 1L:nat) {
          if(side.legend < 2L) {
            if(ticks) {
              graphics::lines(c(at[i], at[i]), c(ylim[side.ticks], ylim[side.ticks] - length.ticks),
                lwd = lwd.ticks[i], lty = lty.ticks[i], col = col.ticks[i])
            }
            if(dl) {
              graphics::text(at[i], ylim[side.ticks] - length.ticks - (distance.labels * length.ticks * 2),
                labels = labels[i], col = col.labels[i], cex = cex.labels[i], pos = 1, ...)
            }
          } else {
            if(ticks) {
              graphics::lines(c(xlim[side.ticks], xlim[side.ticks] - length.ticks), c(at[i], at[i]),
                lwd = lwd.ticks[i], lty = lty.ticks[i], col = col.ticks[i]) 
            }
            if(dl) {
              graphics::text(xlim[side.ticks] - length.ticks - (distance.labels * length.ticks * 2), 
                at[i], labels = labels[i], col = col.labels[i], cex = cex.labels[i],
                pos = if(side.ticks < 2L) 2 else 4, ...)
            }
          }
        }
      } else {
        if(side.legend < 2L) {
          if(side.ticks < 2L) where <- 1L else where <- 3L
        } else {
          if(side.ticks < 2L) where <- 2L else where <- 4L
        }
      axis(where, at = at, labels = labels, col = col.labels, 
        tick = ticks, lty = lty.ticks, col.ticks = col.ticks, 
        lwd.ticks = lwd.ticks, cex.axis = cex.labels[1])
      }
    }
    if(!is.null(title)) {
      if(length(shift.title) < 2)
        shift.title <- c(shift.title, 0)
      if(!full) {
        xp <- xlim[1L] + shift.title[1] * diff(range(xlim)) + diff(range(xlim)) / 2
        yp <- ylim[2L] + shift.title[2] * diff(range(ylim))
        text(if(side.legend < 2) xp else yp,
          if(side.legend < 2) yp else xp, title, pos = 3,
          srt = if(side.legend == 2) 270 else 0, cex = cex.title, xpd = xpd)
      } else {
        mtext(title, side = side.title, cex = cex.title)
      }
    }
  }

  return(invisible(pal))
}


make_pal <- function(col, ncol = NULL, data = NULL, range = NULL, 
  breaks = NULL, swap = TRUE, symmetric = TRUE) 
{
  if(is.null(symmetric))
    symmetric <- TRUE
  if(is.null(col))
    col <- colorspace::diverge_hcl
  if(is.null(ncol) && is.null(breaks))
    ncol <- 99L
  if(is.null(ncol) && !is.null(breaks))
    ncol <- length(breaks) - 1L
  if(is.function(col))
    col <- col(ncol)    
  else 
    ncol <- length(col)
  if(swap) 
    col <- rev(col)
  if(all(is.null(data), is.null(range), is.null(breaks))) 
    stop("at least one needs to be specified")
  if(is.null(breaks)) {
    if(is.null(range)) {
      range <- range(data, na.rm = TRUE)
      if(symmetric) { 
        mar <- max(abs(range))
        range <- c(0 - mar, mar)
      }
    }
    if(diff(range) == 0)
      breaks <- seq(min(range) - 1, min(range) + 1, length.out = ncol + 1L)
    else
      breaks <- seq(range[1L], range[2L], length.out = ncol + 1L)
  } else stopifnot(length(breaks) == ncol + 1L)
  if(is.matrix(data)) {
    obs2col <- function(x) {
      hgt <- (x[-1L, -1L] + x[-1L, -(ncol(x) - 1L)] + 
        x[-(nrow(x) -1L), -1L] + x[-(nrow(x) -1L), -(ncol(x) - 1L)])/4
      c(col[1L], col, col[ncol])[cut(hgt, c(-Inf, breaks, Inf))]
      }
  } else {
    obs2col <- function(x) c(col[1L], col, col[ncol])[cut(x, c(-Inf, breaks, Inf))]
  }

  return(list(colors = col, breaks = breaks, map = obs2col))
}


dopos <- function(pos, limits, width, height, side.legend, shift)
{
  if(side.legend > 1L)
    limits <- rev(limits)
  shift <- rep(shift, length.out = 2)
  shift <- c(diff(limits[[1]]) * shift[1], diff(limits[[2]]) * shift[2])
  if(pos == "bottomleft") {
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + shift[1]
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + shift[2]
  }
  if(pos == "topleft") {
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + shift[1]
    ylim <- c(max(limits[[2L]], na.rm = TRUE) - height, max(limits[[2L]], na.rm = TRUE)) - shift[2]
  }
  if(pos == "topright") {
    xlim <- c(max(limits[[1L]], na.rm = TRUE) - width, max(limits[[1L]], na.rm = TRUE)) - shift[1]
    ylim <- c(max(limits[[2L]], na.rm = TRUE) - height, max(limits[[2L]], na.rm = TRUE)) - shift[2]
  }
  if(pos == "bottomright") {
    xlim <- c(max(limits[[1L]], na.rm = TRUE) - width, max(limits[[1L]], na.rm = TRUE)) - shift[1]
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + shift[2]
  }
  if(pos == "bottom") {
    m <- mean(limits[[1]] - min(limits[[1]], na.rm = TRUE), na.rm = TRUE) - 0.5 * width
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + m
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + shift[2]
  }
  if(pos == "top") {
    m <- mean(limits[[1]] - min(limits[[1]], na.rm = TRUE), na.rm = TRUE) - 0.5 * width
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + m
    ylim <- c(max(limits[[2L]], na.rm = TRUE) - height, max(limits[[2L]], na.rm = TRUE)) - shift[2]
  }
  if(pos == "left") {
    m <- mean(limits[[2]] - min(limits[[2]], na.rm = TRUE), na.rm = TRUE) - 0.5 * height
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + shift[1]
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + m
  }
  if(pos == "right") {
    m <- mean(limits[[2]] - min(limits[[2]], na.rm = TRUE), na.rm = TRUE) - 0.5 * height
    xlim <- c(max(limits[[1L]], na.rm = TRUE) - width, max(limits[[1L]], na.rm = TRUE)) - shift[1]
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + m
  }
  if(pos == "center") {
    mx <- mean(limits[[1]] - min(limits[[1]], na.rm = TRUE), na.rm = TRUE) - 0.5 * width
    my <- mean(limits[[2]] - min(limits[[2]], na.rm = TRUE), na.rm = TRUE) - 0.5 * height
    xlim <- c(min(limits[[1L]], na.rm = TRUE), min(limits[[1L]], na.rm = TRUE) + width) + mx
    ylim <- c(min(limits[[2L]], na.rm = TRUE), min(limits[[2L]], na.rm = TRUE) + height) + my
  }

  return(list(xlim = xlim, ylim = ylim))
}


list2sp <- function(x) 
{
  bndNames <- names(x)
  regions <- unique(bndNames)
  bndAttributes <- attributes(x)
  x <- lapply(x, FUN = function(polygon) {
    polygon <- if(!isTRUE(identical(polygon[1, ], polygon[nrow(polygon), ]))) {
      rbind(polygon, polygon[1, ])
    } else {
      polygon
    }
    as.matrix(na.omit(polygon))
  })
  ret <- list()
  for(id in regions) {
    idMatches <- which(id == bndNames)
    idPolygons <- lapply(na.omit(x[idMatches]), FUN = sp::Polygon, hole = FALSE)
    ret[[id]] <- sp::Polygons(srl = idPolygons, ID = id)
  }
  surrounding <- bndAttributes$surrounding
  whichAreInner <- which(sapply(surrounding, length) > 0L)
  for(innerInd in whichAreInner) {
    hole <- sp::Polygon(coords = x[[innerInd]], hole = TRUE)
    outerId <- surrounding[[innerInd]]
    outerPolys <- ret[[outerId]]@Polygons
    outerPolys <- c(outerPolys, hole)
    ret[[outerId]] <- sp::Polygons(srl = outerPolys, ID = outerId)
  }
  ret <- sp::SpatialPolygons(Srl = ret)
  return(ret)
}

table2df <- function(x)
{
  if(!inherits(x, "table"))
    return(x)
  data.frame("id" = as.factor(names(x)), "x" = as.numeric(x))
}


plotmap <- function(map, x = NA, id = NULL, select = NULL,
  legend = TRUE, names = FALSE, values = FALSE, ...)
{
  if(inherits(map, "bnd") | inherits(map, "list"))
    map <- list2sp(map)

  if(!is.null(id) & inherits(x, "table")) {
    if(length(id) == 1L) {
      x <- data.frame(as.numeric(x), names(x), stringsAsFactors = FALSE)
      names(x) <- c("plotmap_values", as.character(id))
      select <- "plotmap_values"
    }
  }

  if(!inherits(map, "SpatialPolygons"))
    stop("please supply a 'SpatialPolygons' object to argument map!")
  if(!inherits(map, "SpatialPolygonsDataFrame"))
    map <- as(map, "SpatialPolygonsDataFrame")
  if(!is.null(x)) {
    if(all(is.na(x)))
      x <- NULL
  }
  if(!is.null(x)) {
    if(is.data.frame(x) & is.character(id) & (length(id) == 1)) {
      if(!(id %in% names(x)))
        stop("id variable must be part of data x!")
      if(!(id %in% names(map@data)))
        stop("id variable must be part of data in map!")
      map@data[[".sortID"]] <- 1:nrow(map@data)
      map@data <- merge(map@data, x, by = id, all = TRUE)
      map@data <- map@data[order(map@data[[".sortID"]]), , drop = FALSE]
      map@data[[".sortID"]] <- NULL
      if(is.null(select)) {
        select <- names(x)
        select <- select[!(select %in% id)]
        select <- select[1]
      }
    } else {
      if(is.data.frame(x)) {
        if(any(is.f <- apply(x, 2, is.factor))) {
          id <- as.factor(x[, which(is.f)[1]])
          x <- x[, which(!is.f)[1]]
        }
        if(any(is.c <- apply(x, 2, is.character))) {
          id <- as.factor(x[, which(is.c)[1]])
          x <- x[, which(!is.c)[1]]
        }
      }
      if(is.matrix(x)) {
        if(ncol(x > 1)) {
          id <- as.factor(as.integer(x[, 2]))
          x <- x[, 1]
        }
      }
      if(is.null(id))
        stop("no region identifier found!")
      x <- as.numeric(x)
      if(length(x) != length(id))
        stop("length of x != length of id!")
      pol_id <- as.factor(as.character(sapply(slot(map, "polygons"), function(x) {
        slot(x, "ID")
      })))
      if(!any(id %in% pol_id))
        stop("region identifier does not match with polygon ID!")
      if(all(pol_id %in% id)) {
        map@data <- data.frame("ID" = pol_id, "x" = x)
      } else {
        map@data <- data.frame("ID" = pol_id, "x" = NA)
        for(j in seq_along(pol_id)) {
          if(pol_id[j] %in% id) {
            map@data$x[j] <- x[id %in% pol_id[j]]
          }
        }
      }

      select <- "x"
    }
  }

  nm <- names(map)
  args <- list(...)

  if((((length(nm)) < 2 & (nm[1] == "dummy")) | is.null(x)) & is.null(select)) {
    id <- as.character(1:nrow(map@data))
    plot(map, ...)
  } else {
    if(is.null(select))
      select <- nm[1]
    if(!is.character(select)) {
      select <- min(c(select, length(nm)))
      select <- nm[select]
    }
    map@data[[select[1]]] <- as.numeric(map@data[[select[1]]])

    if(is.null(args$color))
      args$color <- heat_hcl
    if(is.null(args$symmetric))
      args$symmetric <- FALSE
    if(is.null(args$swap))
      args$swap <- TRUE
    if(is.null(args$pos))
      args$pos <- "topleft"
    if(is.null(args$digits))
      args$digits <- 2
    if(is.null(args$angle))
      args$angle <- 45
    if(is.null(args$axes))
      args$axes <- FALSE
    if(is.null(args$lty))
      args$lty <- par("lty")
    if(is.null(args$add))
      args$add <- FALSE
    if(is.null(args$border))
      args$border <- 1
    if(is.null(args$lwd))
      args$lwd <- 1

    args$x <- map@data[[select[1]]]
    args$plot <- FALSE
    pal <- do.call("colorlegend", delete.args("colorlegend", args))
    missing <- any(is.na(map@data[[select[1]]]))
    if(missing) {
      if(is.null(args$mdensity))
        args$mdensity <- 20
      plot(map, density = args$mdensity, border = args$border, xlim = args$xlim,
        ylim = args$ylim, xpd = args$xpd, angle = args$angle, pbg = args$pbg,
        axes = args$axes, lty = args$lty, lwd = args$lwd, xlab = args$xlab,
        ylab = args$ylab, main = args$main)
    }
    plot(map, col = pal$map(map@data[[select[1]]]), border = args$border, add = args$add | missing,
      xlim = args$xlim, ylim = args$ylim, xpd = args$xpd, density = args$density,
      angle = args$angle, pbg = args$pbg, axes = args$axes, lty = args$lty,
      lwd = args$lwd, xlab = args$xlab, ylab = args$ylab, main = args$main)

    if(legend) {
      args$add <- TRUE
      pal <- do.call("colorlegend", delete.args("colorlegend", args))
    }

    if(values) {
      co <- coordinates(map)
      text(co[, 1], co[, 2],
        as.character(round(map@data[[select[1]]], digits = args$digits)),
        cex = args$cex)
    }
  }

  if(names) {
    co <- coordinates(map)
    if(is.null(args$names_id)) {
      if(length(id) == 1)
        id <- as.character(map@data[[id]])
    } else {
      id <- as.character(map@data[[args$names_id]])
    }
    text(co[, 1], co[, 2], id, cex = args$cex)
  }
}


.plotmap <- function(map, x = NULL, id = NULL, c.select = NULL, legend = TRUE,
  missing = TRUE, swap = FALSE, range = NULL, names = FALSE, values = FALSE, col = NULL,
  ncol = 100, breaks = NULL, cex.legend = 1, cex.names = 1, cex.values = cex.names,
  digits = 2L, mar.min = 2, add = FALSE, interp = FALSE, grid = 200, land.only = FALSE,
  extrap = FALSE, outside = FALSE, type = "akima", linear = FALSE, k = 40,
  p.pch = 15, p.cex = 1, shift = NULL, trans = NULL, scheme = 1, ...)
{
  if(missing(map))
    stop("map object is missing!")
  if(inherits(map, "SpatialPolygons")) {
    map <- BayesX::sp2bnd(map)
  }
  if(!is.list(map))
    stop("argument map must be a list() of matrix polygons!")
  args <- list(...)
  map.limits <- find.limits(map, mar.min, ...)
  if(is.null(args$asp)) {
    args$asp <- attr(map, "asp")
    if(is.null(args$asp))
      args$asp <- map.limits$asp
  }
  n <- length(map)
  if(is.null(x))
    legend <- FALSE
  poly.names.orig <- names(map)
  if(!any(is.na(poly.names <- x2int(names(map))))) {
    op <- order(poly.names)
    poly.names <- poly.names[op]
    poly.names <- as.character(poly.names)
  } else {
    poly.names <- names(map)
    op <- order(poly.names)
    poly.names <- poly.names[op]
  }
  poly.names.orig <- poly.names.orig[op]
  map <- map[op]
  if(length(upn <- unique(poly.names)) < length(poly.names)) {
    nn <- NULL
    for(i in upn) {
      j <- poly.names == i
      poly.names[j] <- paste(poly.names[j],
        if(sum(j) > 1) paste(".", 1:sum(j), sep = "") else NULL,
        sep = "")
    }
    names(map) <- poly.names
  }
  map <- map[poly.names]
  poly.names <- names(map)
  surrounding <- attr(map, "surrounding")
  inner.id <- which(sapply(surrounding, length) > 0L)
  if(length(inner.id)) {
    poly.names <- c(poly.names[- inner.id], poly.names[inner.id])
    map <- c(map[- inner.id], map[inner.id])
  }
  if(!is.null(args$ylim))
    map.limits$ylim <- args$ylim
  if(!is.null(args$xlim))
    map.limits$xlim <- args$xlim
  if(is.null(args$symmetric))
    symmetric <- TRUE
  else
    symmetric <- args$symmetric
  if(!is.null(x)) {
    if(is.null(col)) {
      col <- colorspace::diverge_hcl
      # col <- colorspace::diverge_hcl(ncol, h = c(130, 10), c = 250,
      #  l = c(30, 90), power = 1.5, gamma = 2.4, fixup = TRUE)
    }
    x <- compute.x.id(x, id, c.select, range, symmetric)
    if(!is.null(shift)) {
      shift <- as.numeric(shift[1])
      x$x <- x$x + shift
    }
    if(!is.null(trans)) {
      if(!is.function(trans)) stop("argument trans must be a function!")
      x$x <- trans(x$x)
    }
    map_fun <- make_pal(col = col, ncol = ncol, data = as.numeric(x$x), 
      range = range, breaks = breaks, swap = swap, 
      symmetric = symmetric)$map
    colors <- map_fun(as.numeric(x$x))
  } else {
    if(is.null(col))
      colors <- rep(NA, length.out = n)
    else {
      if(is.function(col))
        colors <- col(ncol)
      else colors <- col
      colors <- rep(colors, length.out = n)
    }
  }
  if(is.null(args$pos))
    args$pos <- "right"
  if(legend && !is.null(args$pos) && args$pos[1L] == "right") {
    par.orig <- par(c("mar", "las", "mfrow"))
    mar.orig <- mar <- par.orig$mar
    mar[4L] <- 0
    mar[c(1, 3)] <- 1
    on.exit(par(par.orig))
    par(mar = mar)
    w <- (3 + mar[2L]) * par("csi") * 2
    w <- max(c(2.84, w))
    layout(matrix(c(1, 2), nrow = 1), widths = c(1, lcm(w)))
  }
  if(!is.null(map.limits$mar) && is.null(args$asp) && !add)
    par(mar = map.limits$mar)
  args$x <- map.limits$x
  args$y <- map.limits$y
  if(is.null(args$type))
    args$type <- "n"
  if(is.null(args$axes))
    args$axes <- FALSE
  if(is.null(args$xlab))
    args$xlab <- ""
  if(is.null(args$ylab))
    args$ylab <- ""
  if(!add)
    do.call(graphics::plot.default, delete.args(graphics::plot.default, args))
  if(interp & !is.null(x)) {
    cdata <- data.frame(centroids(map), "id" = names(map))
    cdata <- merge(cdata, data.frame("z" = x$x, "id" = x$id), by = "id")
    cdata <- unique(cdata)

    xo <- seq(map.limits$x[1], map.limits$x[2], length = grid)
    yo <- seq(map.limits$y[1], map.limits$y[2], length = grid)
    ico <- interp2(x = cdata[["x"]], y = cdata[["y"]], z = cdata[["z"]],
      xo = xo,
      yo = yo,
      type = type, linear = linear, extrap = extrap,
      k = if(is.null(k)) ceiling(length(map) * 0.8) else as.integer(k))
    
    yco <- rep(yo, each = length(xo))
    xco <- rep(xo, length(yo))

    d4x <- abs(diff(xco))
    d4x <- min(d4x[d4x != 0], na.rm = TRUE)
    d4y <- abs(diff(yco))
    d4y <- min(d4y[d4y != 0], na.rm = TRUE)
    res <- c(d4x, d4y)
    pp <- NULL
    if(length(res))
      pp <- cbind(xco, yco)

    cvals <- as.numeric(ico)
    cvals[cvals < min(cdata$z)] <- min(cdata$z)
    cvals[cvals > max(cdata$z)] <- max(cdata$z)
    icolors <- map_fun(cvals)

    if(!outside) {
      maptools::gpclibPermit()
      class(map) <- "bnd"
      mapsp <- list2sp(map)
      ob <- maptools::unionSpatialPolygons(mapsp, rep(1L, length = length(mapsp)), avoidGEOS  = TRUE)

      nob <- length(slot(slot(ob, "polygons")[[1]], "Polygons"))
      pip <- NULL
      for(j in 1:nob) {
        oco <- slot(slot(slot(ob, "polygons")[[1]], "Polygons")[[j]], "coords")
        pip <- cbind(pip, point.in.polygon(xco, yco, oco[, 1L], oco[, 2L], mode.checked = FALSE) < 1L)
      }
      pip <- apply(pip, 1, function(x) { x < 1 })
    
      icolors[pip] <- NA
    }

    if(land.only) {
      icolors[is.na(maps::map.where("world", xco, yco))] <- NA
    }

    if(length(res)) {
     rect(pp[, 1] - res[1] / 2, pp[, 2] - res[2] / 2, pp[, 1] + res[1] / 2, pp[, 2] + res[2] / 2,
       col = icolors, border = NA, lwd = 0)
    } else {
      points(SpatialPoints(cbind(xco, yco)), col = icolors, pch = p.pch, cex = p.cex)
    }
    colors <- rep(NA, length = length(colors))
  }
  args$ylab <- args$xlab <- args$main <- ""
  args$type <- NULL
  args$axes <- NULL
  lwd.p <- if(!is.null(args$lwd)) rep(args$lwd, length.out = n) else NULL
  if(is.null(lwd.p))
    lwd.p <- rep(1, length.out = n)
  lty.p <- if(!is.null(args$lty)) rep(args$lty, length.out = n) else NULL
  if(is.null(lty.p))
    lty.p <- rep(1, length.out = n)
  border.p <- if(!is.null(args$border)) rep(args$border, length.out = n) else NULL
  if(is.null(border.p))
    border.p <- rep(if(scheme != 1) "white" else "black", length.out = n)
  density.p <- if(!is.null(args$density)) rep(args$density, length.out = n) else NULL
  angle.p <- if(!is.null(args$angle)) rep(args$angle, length.out = n) else NULL
  if(is.null(angle.p))
    angle.p <- rep(90, length.out = n)

  for(poly in unique(poly.names.orig)) {
    for(i in which(poly.names.orig == poly)) {
      args$x <- map[[i]][, 1L]
      args$y <- map[[i]][, 2L]
      args$border <- border.p[i]
      args$angle <- angle.p[i]
      args$lwd <- lwd.p[i]
      args$lty <- lty.p[i]
      if(!is.null(density.p))
        args$density <- density.p[i]
      if(!is.null(x)){ 
        if(!is.na(k <- pmatch(poly, x$id))) {
          args$col <- colors[k]
          args$density <- NULL
        } else {
          if(!missing) next
          args$col <- NULL
          if(is.null(args$density))
            args$density <- 20L
        }
      } else args$col <- colors[i]
      do.call(graphics::polygon, 
        delete.args(graphics::polygon, args, 
        c("lwd", "cex", "lty")))
      if(names && !values) {
        args$polygon <- map[[i]]
        args$poly.name <- poly.names.orig[i]
        args$counter <- i
        args$cex <- cex.names
        do.call(centroidtext, delete.args(centroidtext, args, "font"))
      }
      if(values && !names) {
        args$polygon <- map[[i]]
        args$poly.name <- as.character(round(x$x[k], digits = digits))
        args$counter <- k
        args$cex <- cex.values
        do.call(centroidtext, delete.args(centroidtext, args, "font"))
      }
    }
  }

  if(legend) {
    if(is.null(args$pos))
      args$pos <- "topleft"
    if(args$pos[1L] == "right") {
      args$full <- TRUE
      args$side.legend <- 2L
      args$side.ticks <- 2L
      mar <- mar.orig
      mar[2L] <- 0.5
      mar[4L] <- 3.1
      par(mar = mar, xaxs = "i", yaxs = "i")
      args$plot <- TRUE
      args$add <- FALSE
    } else {
      args$plot <- FALSE
      if(is.null(args$xpd))
        args$xpd <- TRUE
      args$add <- TRUE
    }
    args$shift <- args$legend.shift
    args$xlim <- map.limits$xlim
    args$ylim <- map.limits$ylim
    args$color <- col
    args$ncol <- ncol
    args$x <- x$x
    args$breaks <- breaks
    args$swap <- swap
    args$digits <- digits
    args$cex.labels <- cex.legend
    args$symmetric <- symmetric
    if(is.null(range)) {
      range <- range(args$x)
      if(diff(range) == 0)
        range <- unique(range) + c(-1, 1)
    }
    args$range <- range
    if(is.null(args$lrange))
      args$lrange <- args$range
    do.call(colorlegend, delete.args(colorlegend, args, c("font")))
  }
  if(!is.null(args$xlab))
    mtext(args$xlab, side = 1L)
  if(!is.null(args$ylab))
    mtext(args$ylab, side = 2L)

  return(invisible(NULL))
}


plot.bnd <- function(x, ...) {
  plotmap(x, ...)
}


compute.x.id <- function(x, id = NULL, c.select = NULL, range = NULL, symmetric = TRUE)
{ 
  if(is.null(id) && (is.vector(x) || is.array(x))) {
    if(!is.null(names(x))) {
      id <- names(x)
      x <- as.vector(x)
    }
  }
  if(is.factor(id))
    id <- as.character(id)
  if(is.array(x) && length(dim(x)) < 2L)
    x <- as.vector(x)
  if(is.null(dim(x)) && is.null(dim(id))) {
    if(length(x) != length(id))
      stop("arguments x and id are differing!")
  } else {
    x <- unclass(x)
    if(is.list(x)) 
      nx <- names(x)
    if(is.matrix(x)) {
      if(ncol(x) < 2 & !is.null(id)) {
        x <- data.frame("id" = id, "x" = as.numeric(x))
        nx <- names(x)
        c.select <- "x"
        id <- NULL
      } else {
        x <- as.list(as.data.frame(x))
        nx <- names(x)
        if(all(nx %in% paste("V", 1L:length(nx), sep = ""))) {
          nx[1L:2L] <- c("id", "x")
          c.select <- "x"
        }
      }
    }
    if(is.data.frame(x)) {
      x <- as.list(x)
      nx <- names(x)
    }
    if(is.null(id))
      id <- x[[1L]]
    else {
      if(is.character(id)) {
        if(is.na(id <- pmatch(id, nx)))
          stop("argument id is specified wrong!")
      } else {
        if(id > length(nx))
          stop("argument id is specified wrong!")
      }
      id <- x[[id]]
    }
    if(is.null(c.select)) {
      take <- c("mean", "Mean", "MEAN", "estimate", 
        "Estimate", "ESTIMATE", "mean", "pmode", "pmean_tot")
      did.take <- FALSE
      for(k in take) {
        if(!is.na(pmatch(k, nx)) & !did.take) {
          x <- x[[k]]
          did.take <- TRUE
        }
      }
     if(!did.take && length(x) > 1L)
       x <- x[[2L]]
    } else {
      if(is.character(c.select)) {
        k <- pmatch(c.select, nx)
      if(is.na(k))
        stop("argument c.select is specified wrong!")
      x <- x[[k]]
      } else {
        if(c.select > length(nx))
          stop("argument c.select is specified wrong!")
        x <- x[[c.select]]
      }
    }
  }
  xrange <- range(x, na.rm = TRUE)
  if(symmetric) {
    xrange <- c(-1 * max(abs(xrange)), max(abs(xrange))) 
    if(is.null(range)) {
      if(min(x) < 0)
        m <- (-1)
      else
        m <- 1
      if(abs(min(x)) > abs(max(x)))
        x <- c(x, abs(min(x)))
      if(abs(max(x)) > abs(min(x)))
        x <- c(x, m * abs(max(x)))
      id <- c(as.character(id), "added")
    } else {
      if(max(range) > max(x)) {
        x <- c(x, max(range))
        id <- c(as.character(id), "added")
      } else x[x > max(range)] <- max(range)
      if(min(range) < min(x)) {
        x <- c(x, min(range))
        id <- c(as.character(id), "added")
      } else x[x < min(range)] <- min(range)
    }
  }

  return(list(id = as.character(id), x = x, range = xrange))
}


df2m <- function(x)
{
  if(!is.null(x)) {
    xattr <- attributes(x)
    nxa <- names(xattr)
    x$intnr <- NULL
    cn <- colnames(x)
    x <- as.matrix(x)
    rownames(x) <- 1L:nrow(x)
    colnames(x) <- rep(cn, length.out = ncol(x))
    for(k in 1L:length(nxa)) 
      if(all(nxa[k] != c("dim", "dimnames", "class", "names", "row.names"))) {
        attr(x, nxa[k]) <- xattr[[k]]
      }
    }

  return(x)
}


find.limits <- function(map, mar.min = 2, ...)
{
  if(!is.list(map))
    stop("argument map must be a list() of matrix polygons!")
  n <- length(map)
  myrange <- function(x, c.select = 1L, ...) {
    if(!is.null(dim(x))) {
      return(na.omit(x[, c.select], ...))
    }
    else return(na.omit(x))
  }
  xlim <- range(unlist(lapply(map, myrange, c.select = 1L, ...)))
  ylim <- range(unlist(lapply(map, myrange, c.select = 2L, ...)))
  mar <- NULL
  asp <- attr(map, "asp")
  if(is.null(asp))
    asp <- (diff(ylim) / diff(xlim)) / cos((mean(ylim) * pi) / 180)
  if(!is.null(height2width <- attr(map, "height2width"))) {
    height2width <- height2width * 0.8
    if(!is.null(mar.min)) {
      if(height2width > 1) {
        side <- 17.5 * (1 - 1/height2width) + mar.min / height2width
        mar <- c(mar.min, side, mar.min, side)
      }
      else {
        top <- 17.5  * (1 - height2width) + mar.min * height2width
        mar <- c(top, mar.min, top, mar.min)
      }
    }
  }

  return(list(ylim = ylim, xlim = xlim, mar = mar, asp = asp))
}


set.plot2d.specs <- function(nc, args, col.lines, is.bayesx)
{
  lwd <- args$lwd
  if(is.null(lwd) || any(is.na(lwd)))
    lwd <- rep(1, nc)
  else
    lwd <- rep(lwd, length.out = nc)
  lty <- args$lty
  if((is.null(lty) || any(is.na(lty))) && !is.bayesx[1L])
    lty <- rep(1, nc)
  if((is.null(lty) || any(is.na(lty))) && is.bayesx[1L])
    lty <- c(1, 0, 0, 0, 0)
  if(length(lty) == 1L || length(lty) < nc)
    lty <- rep(lty, length.out = nc)
  if(is.null(col.lines))
    col.lines <- rep("black", nc)
  else
    col.lines <- rep(col.lines, length.out = nc)
  args$lty <- lty
  args$lwd <- lwd
  args$col.lines <- col.lines

  return(args)
}


interp2 <- function(x, y, z, xo = NULL, yo = NULL, grid = 30,
  type = c("mba", "akima", "mgcv", "gam", "raw"), linear = FALSE, extrap = FALSE, k = 40)
{
  type <- tolower(type)
  type <- match.arg(type)

  if(is.null(xo))
    xo <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = grid)
  if(is.null(yo))
    yo <- seq(min(y, na.rm = TRUE), max(y, na.rm = TRUE), length = grid)

  xgrid <- length(xo)
  ygrid <- length(yo)
  x <- as.numeric(x); y <- as.numeric(y); z <- as.numeric(z)

  if(type %in% c("mgcv", "gam")) {
    xo <- as.numeric(xo); yo <- as.numeric(yo)
    xr <- range(x, na.rm = TRUE)
    yr <- range(y, na.rm = TRUE)
    x <- (x - xr[1]) / diff(xr)
    y <- (y - yr[1]) / diff(yr)

    if(k > length(z))
      k <- ceiling(0.8 * length(z))

    b <- mgcv::gam(z ~ s(x, y, k = k))

    x2 <- (xo - xr[1]) / diff(xr)
    y2 <- (yo - yr[1]) / diff(yr)
    nd <- data.frame("x" = rep(x2, ygrid), "y" = rep(y2, rep(xgrid, xgrid)))
    fit <- as.vector(predict(b, newdata = nd))

    if(!extrap) {
      pid <- chull(X <- cbind(x, y))
      pol <- X[c(pid, pid[1]), ]
      pip <- point.in.polygon(nd$x, nd$y, pol[, 1], pol[, 2])
      fit[pip < 1] <- NA
    }
  }
  if(type == "mba") {
    fit <- MBA::mba.surf(data.frame("x" = x, "y" = y, "z" = z), xgrid, ygrid)$xyz.est$z
  }
  if(type == "akima") {
    if(isTRUE(getOption("use.akima"))) {
      stopifnot(requireNamespace("akima"))
    } else {
      if(requireNamespace("akima")) {
        cat("NOTE: Package 'akima' has an ACM license that restricts applications to non-commercial usage.\n")
      } else {
        stop(paste("plot3d() can only be used if the 'akima' package is installed. ",
          "Note that 'akima' has an ACM license that restricts applications to ",
          "non-commercial usage.", sep = ""))
      }
    }

    fit <- try(akima::interp(x, y, z, xo = xo, yo = yo, 
      duplicate = "mean", linear = linear, extrap = extrap)$z, silent = TRUE)
    if(inherits(fit, "try-error") | all(is.na(fit))) {
      cat("NOTE: akima::interp() is designed for irregular data points, the coordinates will be slightly jittered to obtain irregular spaced points.\n")
      fit <- try(akima::interp(jitter(x, amount = .Machine$double.eps),
        jitter(y, amount = .Machine$double.eps), z, xo = xo, yo = yo, 
        duplicate = "mean", linear = linear, extrap = extrap)$z, silent = TRUE)
    }
  }

  return(matrix(fit, xgrid, ygrid))
}


x2int <- function(x) 
{
  warn <- getOption("warn")
  options(warn = -1)
  rval <- as.integer(as.numeric(as.character(x)))
  options("warn" = warn)
  rval
}


sliceplot <- function(x, y = NULL, z = NULL, view = 1, c.select = NULL,
  values = NULL, probs = c(0.1, 0.5, 0.9), grid = 100,
  legend = TRUE, pos = "topright", digits = 2, data = NULL,
  rawdata = FALSE, type = "mba", linear = FALSE, extrap = FALSE,
  k = 40, rug = TRUE, rug.col = NULL, jitter = TRUE, ...)
{
  if(is.vector(x) & is.vector(y) & is.vector(z)) {
    nx <- c(
      deparse(substitute(x), backtick = TRUE, width.cutoff = 500),
      deparse(substitute(y), backtick = TRUE, width.cutoff = 500),
      deparse(substitute(z), backtick = TRUE, width.cutoff = 500)
    )
    x <- cbind(x, y, z)
    colnames(x) <- nx
  } else {
    if(inherits(x,"formula")) {
      if(is.null(data))
        data <- environment(x)
      else
        if(is.matrix(data))
          data <- as.data.frame(data)
      x <- model.frame(x, data = data)
      if(ncol(x) < 3L)
        stop("formula is specified wrong!")
      if(ncol(x) > 3L)
        x <- x[, c(2L, 3L, 1L, 4L:ncol(x))]
      else
        x <- x[, c(2L, 3L, 1L)]
    }
  }
  stopifnot(is.matrix(x) || is.data.frame(x))
  nx <- colnames(x)
  if(is.null(c.select))
    c.select <- 3
  if(c.select < 3)
    c.select <- if(c.select < 2) 3 else 4 
  if(c.select > ncol(x))
    stop("column number selected in c.select is larger than the number of existing columns in x!")
  view <- view[1]
  if(is.character(view))
    view <- grep(view, nx, ignore.case = TRUE)
  x <- x[order(x[, view]), ]
  noview <- if(view < 2) 2 else 1
  values <- if(is.null(values)) {
    quantile(x[, noview], probs = probs, type = 1)
  } else values
  if(!rawdata) {
    xo <- seq(min(x[, view]), max(x[, view]), length = grid)
    yo <- seq(min(x[, noview]), max(x[, noview]), length = grid)
    zi <- interp2(x[, view], x[, noview], x[, c.select],
      xo = xo,
      yo = yo,
      type = type, linear = linear, extrap = extrap, k = k)
    yg <- rep(yo, each = grid)
    zg <- as.vector(zi)
    slices <- xo
  } else {
    yg <- x[, noview]
    zg <- x[, c.select]
    slices <- unique(x[, view])
  }
  for(j in values) {
    val <- unique(yg[which.min(abs(yg - j))])
    slices <- cbind(slices, zg[yg == val])
  }
  k <- ncol(slices)
  args <- l.args <- list(...)
  args$lty <- if(is.null(args$lty)) 1:k else rep(args$lty, length.out = k)
  args$col <- if(is.null(args$col)) "black" else rep(args$col, length.out = k)
  args$lwd <- if(is.null(args$lwd)) 1 else rep(args$lwd, length.out = k)
  if(is.null(args$xlab))
    args$xlab <- nx[view]
  if(is.null(args$ylab))
    args$ylab <- paste("Effect of", nx[view])
  args$x <- slices[, 1]
  args$y <- slices[, 2:ncol(slices)]
  args$type = "l"
  do.call("matplot", delete.args("matplot", args, c("axes", "main", "xlab", "ylab")))
  if(legend) {
    l.args$x <- pos
    l.args$legend <- paste(nx[noview], "=", round(values, digits))
    l.args$lty <- args$lty
    l.args$col <- args$col
    l.args$lwd <- args$lwd
    if(is.null(l.args$bg))
      l.args$bg <- NA
    if(is.null(l.args$box.col))
      l.args$box.col <- NA
    do.call("legend", delete.args("legend", l.args))
  }
  if(rug) {
    args$x <- if(jitter) jitter(x[, view]) else x[, view]
    args$col <- rug.col
    do.call(graphics::rug, delete.args(graphics::rug, args))
  }
  invisible(args)
}

