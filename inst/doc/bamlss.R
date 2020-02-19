## ----preliminaries, echo=FALSE, message=FALSE, results="hide"-----------------
library("bamlss")
prefix <- "http://www.bamlss.org/articles/" ## ""
prefix2 <- "http://www.bamlss.org/reference/" ## ""

## ----installation-cran, eval=FALSE--------------------------------------------
#  install.packages("bamlss")

## ----installation-rforge, eval=FALSE------------------------------------------
#  install.packages("bamlss", repos = "http://R-Forge.R-project.org")

## -----------------------------------------------------------------------------
data("Golf", package = "bamlss")
head(Golf)

## -----------------------------------------------------------------------------
f <- price ~ age + kilometer + tia + abs + sunroof

## ---- message=FALSE, results="hide"-------------------------------------------
library("bamlss")

set.seed(111)

b1 <- bamlss(f, family = "gaussian", data = Golf)

## -----------------------------------------------------------------------------
summary(b1)

## ---- eval=FALSE--------------------------------------------------------------
#  plot(b1, which = "samples")

## ---- fig.width = 9, fig.height = 5, fig.align = "center", echo = FALSE, dev = "png", results = 'hide', message=FALSE----
bsp <- b1
bsp$samples <- bsp$samples[, c("mu.p.(Intercept)", "sigma.p.(Intercept)")]
plot(bsp, which = "samples")

## -----------------------------------------------------------------------------
confint(b1, prob = c(0.025, 0.975))

## ---- message=FALSE, results="hide"-------------------------------------------
set.seed(111)

f <- log(price) ~ age + kilometer + tia + abs + sunroof

b2 <- bamlss(f, family = "gaussian", data = Golf)

## -----------------------------------------------------------------------------
DIC(b1, b2)

## ---- message=FALSE, results="hide"-------------------------------------------
set.seed(222)

f <- log(price) ~ poly(age, 3) + poly(kilometer, 3) + poly(tia, 3) + abs + sunroof

b3 <- bamlss(f, family = "gaussian", data = Golf)

## -----------------------------------------------------------------------------
DIC(b1, b2, b3)

## -----------------------------------------------------------------------------
nd <- data.frame("age" = seq(min(Golf$age), max(Golf$age), length = 100))

nd$page <- predict(b3, newdata = nd, model = "mu", term = "age",
  FUN = c95, intercept = FALSE)

head(nd)

## -----------------------------------------------------------------------------
nd$kilometer <- seq(min(Golf$kilometer), max(Golf$kilometer), length = 100)
nd$tia <- seq(min(Golf$tia), max(Golf$tia), length = 100)

nd$pkilometer <- predict(b3, newdata = nd, model = "mu", term = "kilometer",
  FUN = c95, intercept = FALSE)
nd$ptia <- predict(b3, newdata = nd, model = "mu", term = "tia",
  FUN = c95, intercept = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  par(mfrow = c(1, 3))
#  ylim <- range(c(nd$page, nd$pkilometer, nd$ptia))
#  plot2d(page ~ age, data = nd, ylim = ylim)
#  plot2d(pkilometer ~ kilometer, data = nd, ylim = ylim)
#  plot2d(ptia ~ tia, data = nd, ylim = ylim)

## ---- fig.width = 8, fig.height = 2.4, fig.align = "center", dev = "png", results='hide', message=FALSE, echo=FALSE, out.width="100%"----
par(mfrow = c(1, 3), mar = c(4.1, 4.1, 0.1, 1.1))
ylim <- range(c(nd$page, nd$pkilometer, nd$ptia))
plot2d(page ~ age, data = nd, ylim = ylim)
plot2d(pkilometer ~ kilometer, data = nd, ylim = ylim)
plot2d(ptia ~ tia, data = nd, ylim = ylim)

## ----preliminaries2, echo=FALSE, message=FALSE, results="hide"----------------
data("mcycle", package = "MASS")
f <- list(accel ~ s(times, k = 20), sigma ~ s(times, k = 20))
set.seed(123)
b <- bamlss(f, data = mcycle, family = "gaussian")

## -----------------------------------------------------------------------------
data("mcycle", package = "MASS")
head(mcycle)

## -----------------------------------------------------------------------------
f <- list(accel ~ s(times, k = 20), sigma ~ s(times, k = 20))

## ---- eval=FALSE--------------------------------------------------------------
#  set.seed(123)
#  
#  b <- bamlss(f, data = mcycle, family = "gaussian")

## ---- eval=FALSE--------------------------------------------------------------
#  plot(b, model = c("mu", "sigma"))

## ---- fig.width = 7, fig.height = 3, fig.align = "center", echo = FALSE, dev = "png", results = 'hide', message=FALSE, out.width="80%"----
par(mar = c(4.1, 4.1, 1.1, 1.1), mfrow = c(1, 2))
plot(b, pages = 1, spar = FALSE, scheme = 2, grid = 100)

## -----------------------------------------------------------------------------
summary(b)

## ---- eval=FALSE--------------------------------------------------------------
#  plot(b, which = "samples")

## ---- fig.width = 9, fig.height = 5, fig.align = "center", echo = FALSE, dev = "png", results = 'hide', message=FALSE----
bsp <- b
bsp$samples <- bsp$samples[, c("mu.p.(Intercept)", "sigma.p.(Intercept)")]
plot(bsp, which = "samples")

## ---- eval = FALSE------------------------------------------------------------
#  plot(b, which = c("hist-resid", "qq-resid"))

## ---- fig.width = 7.5, fig.height = 4, fig.align = "center", echo = FALSE, dev = "png", results = 'hide', message=FALSE, out.width="80%"----
par(mfrow = c(1, 2))
plot(b, which = c("hist-resid", "qq-resid"), spar = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  nd <- data.frame("times" = seq(2.4, 57.6, length = 100))
#  nd$ptimes <- predict(b, newdata = nd, model = "mu", FUN = c95)
#  plot2d(ptimes ~ times, data = nd)

## ---- fig.width = 5, fig.height = 4, fig.align = "center", dev = "png", results='hide', message=FALSE, echo=FALSE, out.width="50%"----
par(mar = c(4.1, 4.1, 1.1, 1.1))
nd <- data.frame("times" = seq(2.4, 57.6, length = 100))
nd$ptimes <- predict(b, newdata = nd, model = "mu", FUN = c95)
plot2d(ptimes ~ times, data = nd)

## ---- fig.width = 5, fig.height = 4, fig.align = "center", dev = "png", out.width="50%"----
## Predict for the two scenarios.
nd <- data.frame("times" = c(10, 40))
ptimes <- predict(b, newdata = nd, FUN = identity, type = "parameter")

## Extract the family object.
fam <- family(b)

## Compute densities.
dens <- list("t10" = NULL, "t40" = NULL)
for(i in 1:ncol(ptimes$mu)) {
  ## Densities for times = 10.
  par <- list(
    "mu" = ptimes$mu[1, i, drop = TRUE],
    "sigma" = ptimes$sigma[1, i, drop = TRUE]
  )
  dens$t10 <- cbind(dens$t10, fam$d(mcycle$accel, par))

  ## Densities for times = 40.
  par <- list(
    "mu" = ptimes$mu[2, i, drop = TRUE],
    "sigma" = ptimes$sigma[2, i, drop = TRUE]
  )
  dens$t40 <- cbind(dens$t40, fam$d(mcycle$accel, par))
}

## Visualize.
par(mar = c(4.1, 4.1, 0.1, 0.1))
col <- rainbow_hcl(2, alpha = 0.01)
plot2d(dens$t10 ~ accel, data = mcycle,
  col.lines = col[1], ylab = "Density")
plot2d(dens$t40 ~ accel, data = mcycle,
  col.lines = col[2], add = TRUE)

