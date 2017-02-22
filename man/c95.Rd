\name{c95}
\alias{c95}

\title{
  Compute 95\% Credible Interval and Mean
}

\description{
  Small helper function that computes the 2.5\% and 97.5\% quantiles and the
  mean of a vector. Useful for example when using function \code{\link{predict.bamlss}}.
}

\usage{
c95(x)
}

\arguments{
  \item{x}{A numeric vector.}
}

\seealso{
  \code{\link{predict.bamlss}}, \code{\link{coef.bamlss}}
}

\examples{
x <- rnorm(100)
c95(x)

\dontrun{## Example computing predictions.
set.seed(123)
d <- data.frame("x" = seq(-3, 3, length = 30))
d$y <- sin(d$x) + rnorm(30, sd = 0.3)

## Estimate model and compute predictions.
## with c95().
b <- bamlss(y ~ s(x), data = d)
p <- predict(b, model = "mu", FUN = c95)
plot(d)
matplot(d$x, p, type = "l", lty = c(2, 1, 2),
  col = "black", add = TRUE)

## Example extracting coefficients.
coef(b, FUN = c95)
}
}

\keyword{regression}

