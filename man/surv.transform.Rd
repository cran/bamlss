\name{surv_transform}
\alias{surv_transform}

\title{
  Survival Model Transformer Function
}

\description{
  This function takes a \code{\link{bamlss.frame}} and computes design matrices of model terms
  based on a time grid for time-dependent structured additive predictors in a survival context.
  Note that this transformer function is usually used internally by function \code{\link{bamlss}}
  and is the default transformer function using the \code{\link{cox_bamlss}} family object.

  The time grid design matrices can be used to construct the full structured additive predictor
  for each time point. This way it is possible to solve the integrals that are part of, e.g.,
  a Newton-Raphson updating scheme, numerically.

  See the example section on how to extract the time grid design matrices.
}

\usage{
surv_transform(x, y, data, family,
  subdivisions = 100, timedependent = "lambda",
  timevar = NULL, idvar = NULL, is.cox = FALSE,
  alpha = 0.1, ...)
}

\arguments{
 \item{x}{The \code{x} list, as returned from function
    \code{\link{bamlss.frame}} and transformed by function \code{\link{surv_transform}},
    holding all model matrices and other information that is used for
    fitting the model.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{data}{The data.frame that should be used for setting up all matrices.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.
    In this case this is the \code{\link{cox_bamlss}} family object.}
  \item{subdivisions}{How many time points should be created for each individual.}
  \item{timedependent}{A character vector specifying the names of parameters in \code{x}
    that are time-dependent. Time grid design matrices are only computed for these parameters.}
  \item{timevar}{A character specifying the name of the survival time variable in the
    data set.}
  \item{idvar}{Depending on the type of data set, this is the name of the variable specifying
     identifier of individuals.}
  \item{is.cox}{Should the \code{\link{bamlss.frame}} be set up for a Cox type survival model.}
  \item{alpha}{A value for the intercept of a parameter names alpha. Typically the association
    parameter of a longitudinal and survival process in a joint model.}
  \item{\dots}{Arguments passed to function \code{\link{bamlss.engine.setup}}.}
}

\value{
  A \code{\link{bamlss.frame}} including the time grid design matrices.
}

\seealso{
  \code{\link{cox_bamlss}}, \code{\link{cox_mode}}, \code{\link{cox_mcmc}},
  \code{\link{simSurv}}, \code{\link{bamlss}}
}

\examples{
library("survival")
set.seed(111)

## Simulate survival data.
d <- simSurv(n = 20)

## Formula of the survival model, note
## that the baseline is given in the first formula by s(time).
f <- list(
  Surv(time, event) ~ s(time) + s(time, by = x3),
  gamma ~ s(x1) + s(x2)
)

## Create the bamlss.frame.
bf <- bamlss.frame(f, family = "cox", data = d)

## Lambda is the time-dependent parameter.
print(bf)

## Apply the transformer.
bf <- with(bf, surv_transform(x, y, data = model.frame,
  family = family, is.cox = TRUE, subdivisions = 25))

## Extract the time grid design matrix for term s(time).
X <- bf$x$lambda$smooth.construct[["s(time)"]]$fit.fun_timegrid(NULL)
dim(X)

## Compute fitted values for each time point.
grid <- attr(bf$y[[1]], "grid")
gdim <- c(length(grid), length(grid[[1]]))
b <- runif(ncol(X))
fit <- X %*% b
fit <- matrix(fit, nrow = gdim[1], ncol = gdim[2], byrow = TRUE)

plot(as.vector(fit) ~ unlist(grid), type = "n",
  xlab = "Survival time", ylab = "Effect")
for(j in seq_along(grid)) {
  lines(fit[j, ] ~ grid[[j]], lwd = 2, col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
  points(grid[[j]][gdim[2]], fit[j, gdim[2]], col = "red")
}
}

\keyword{regression}
\keyword{survival}

