\name{cox_predict}
\alias{cox_predict}

\title{
  Cox Model Prediction
}

\description{
  This function takes a fitted Cox model, i.e., a model estimated by \code{\link{cox_mode}}
  or \code{\link{cox_mcmc}} and computes predictions given a new data set or the
  original values. Survival probabilities are computed using numerical integration, therefore,
  computation may take some time. To avoid problems with computer memory, the prediction
  of survival probabilities can be split into chunks and computed parallel on different cores.
}

\usage{
cox_predict(object, newdata,
  type = c("link", "parameter", "probabilities"),
  FUN = function(x) { mean(x, na.rm = TRUE) },
  time = NULL, subdivisions = 100, cores = NULL,
  chunks = 1, verbose = FALSE, ...)
}

\arguments{
  \item{object}{A \code{"bamlss"} object as returned from function \code{\link{bamlss}}
    using the optimizer \code{\link{cox_mode}} or sampler function \code{\link{cox_mcmc}}.}
  \item{newdata}{A data frame or list containing the values of the model
    covariates at which predictions are required. If missing \code{newdata} is the
    \code{model.frame} of the provided model.}
  \item{type}{Specifies the type of predictions that should be computed.}
  \item{FUN}{A function that should be applied on each row of the samples
    of the additive predictor, parameter or probabilities. Per default
    the function computes means of samples, however, other functions like \code{\link{quantile}}
    can be supplied.}
  \item{time}{\code{numeric}, specifies the time for which survival probabilities should
    be computed if \code{type = "probabilities"}. Note that this overwrites survival times
    that are supplied in argument \code{newdata}.}
  \item{subdivisions}{How many time points should be created for each individual.}
  \item{cores}{Specifies the number of cores that should be used for prediction. The problem is
    split into \code{core} chunks, each chunk is then processed by one core.}
  \item{chunks}{The number of chunks that should be processed sequentially on one core.
    This way memory problems can be avoided when computing survival times for large problems.}
  \item{verbose}{Print progress information.}
  \item{\dots}{Arguments passed to \code{\link{predict.bamlss}}.}
}

\value{
  Depending on the type of function provided in argument \code{FUN}, a \code{numeric} vector
  or \code{matrix}.
}

\seealso{
  \code{\link{cox_mcmc}}, \code{\link{cox_bamlss}}, \code{\link{surv_transform}},
  \code{\link{simSurv}}, \code{\link{bamlss}}, \code{\link{predict.bamlss}}
}

\examples{
\dontrun{library("survival")
set.seed(123)

## Simulate survival data.
d <- simSurv(n = 500)

## Formula of the survival model, note
## that the baseline is given in the first formula by s(time).
f <- list(
  Surv(time, event) ~ s(time) + s(time, by = x3),
  gamma ~ s(x1) + s(x2)
)

## Cox model with continuous time.
## Note the the family object cox_bamlss() sets
## the default optimizer and sampler function!
## First, posterior mode estimates are computed
## using function cox_mode(), afterwards the
## sampler cox_mcmc() is started.
b <- bamlss(f, family = "cox", data = d)

## Predict survival probabilities P(T > t).
p <- predict(b, type = "probabilities",
  time = 3, subdivisions = 100, FUN = c95)
}
}

\keyword{regression}
\keyword{survival}

