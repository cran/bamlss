\name{trans_AR1}
\alias{trans_AR1}
\alias{AR1}

\title{
  AR1 Transformer Function
}

\description{
  The transformer function takes a \code{\link{bamlss.frame}} object and transforms
  the response and the design matrices to account for lag 1 autocorrelation. The method
  is also known as Prais-Winsten estimation.
}

\usage{
trans_AR1(rho = 0.1)
AR1(rho = 0.1)
}

\arguments{
  \item{rho}{Specifies the correlation parameter at lag 1.}
}

\value{
  A transformer function which can be used in the \code{\link{bamlss}} call.
}

\references{
  Johnston, John (1972). Econometric Methods (2nd ed.). New York: McGraw-Hill. pp. 259--265.
}

\seealso{
  \code{\link{bamlss.frame}}, \code{\link{bamlss}}, \code{\link[mgcv]{smooth2random}}.
} 

\examples{\dontrun{## Simulate AR1 data.
set.seed(111)

n <- 240
d <- data.frame("t" = 1:n)

## Nonlinear function.
f <- function(x) {
  2 + sin(x / n * 2 * pi - pi)
}

## Correlated errors.
rho <- 0.8
e <- rnorm(n, sd = 0.1)
u <- c(e[1], rep(NA, n - 1))
for(i in 2:n){
  u[i] <- rho * u[i - 1] + e[i]
}

## Response.
d$y <- f(d$t) + u

## Plot time-series data.
plot(d, type = "l")

## Estimate models without and with AR1 transformation.
b0 <- bamlss(y ~ s(t,k=20), data = d, criterion = "BIC")
b1 <- bamlss(y ~ s(t,k=20), data = d, criterion = "BIC",
  transform = AR1(rho = 0.8))

## Estimate full AR1 model.
b2 <- bamlss(y ~ s(t,k=20), data = d, criterion = "BIC",
  family = "AR1")
rho <- predict(b2, model = "rho", type = "parameter")
print(range(rho))

## Estimated standard deviations.
sd0 <- predict(b0, model = "sigma", type = "parameter")
sd1 <- predict(b1, model = "sigma", type = "parameter")
sd2 <- predict(b2, model = "sigma", type = "parameter")
print(round(c(sd0[1], sd1[1], sd2[1]), 2))

## Plot fitted trends.
p0 <- predict(b0, model = "mu")
p1 <- predict(b1, model = "mu")
p2 <- predict(b2, model = "mu")

plot(d, type = "l")
lines(f(d$t) ~ d$t, col = 2, lwd = 2)
lines(p0 ~ d$t, col = 4, lwd = 2)
lines(p1 ~ d$t, col = 3, lwd = 3)
lines(p2 ~ d$t, col = 5, lwd = 3)
legend("topleft",
  c("no trans", "with trans", "AR1 model", "truth"),
  lwd = 2, col = c(4, 3, 5, 2), bty = "n")
}}

\keyword{regression}
\keyword{smooth}

