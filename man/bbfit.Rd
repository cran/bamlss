\name{opt_bbfit}
\alias{bbfit}
\alias{opt_bbfit}
\alias{bbfitp}
\alias{opt_bbfitp}
\alias{contribplot}

\title{Batchwise Backfitting}

\description{
  Batchwise backfitting estimation engine for GAMLSS using very large data sets.
}

\usage{
## Batchwise backfitting engine.
opt_bbfit(x, y, family, shuffle = TRUE, start = NULL, offset = NULL,
  epochs = 1, nbatch = 10, verbose = TRUE, ...)

bbfit(x, y, family, shuffle = TRUE, start = NULL, offset = NULL,
  epochs = 1, nbatch = 10, verbose = TRUE, ...)

## Parallel version.
opt_bbfitp(x, y, family, mc.cores = 1, ...)

## Loglik contribution plot.
contribplot(x, ...)
}

\arguments{
   \item{x}{For function \code{bfit()} the \code{x} list, as returned from function
    \code{\link{bamlss.frame}}, holding all model matrices and other information that is used for
    fitting the model. For the updating functions an object as returned from function
    \code{\link[mgcv]{smooth.construct}} or \code{\link[mgcv]{smoothCon}}. For function
    \code{contribplot()}, a \code{"bamlss"} object using \code{bbfit()} with argument
    \code{select = TRUE}.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see \code{\link{family.bamlss}}.}
  \item{shuffle}{Should observations be shuffled?}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{epochs}{For how many epochs should the algorithm run?}
  \item{nbatch}{Number of batches. Can also be a number between 0 and 1, i.e., determining
    the fraction of observations that should be used for fitting.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{mc.cores}{On how many cores should estimation be started?}
  \item{\dots}{For \code{bbfitp()} all arguments to be passed to \code{bbfit()}.}
}

\details{
  The algorithm uses batch-wise estimation of regression coefficients and smoothing variances.
  The smoothing variances are estimated on an hold-out batch. This way, models for very large
  data sets can be estimated. Note, the algorithm can work in combination with the \pkg{ff}
  and \pkg{ffbase} package, i.e., the entire data is never in the computer RAM. Therefore, the
  data can either to be stored as comma separated file on disc or provided as \code{"ffdf"}
  data frame, see also the examples.

  The optimizer functions use additional arguments:
  \itemize{
    \item \code{batch_ids}. This argument can either be a list of indices specifying the
      batches that should be used for estimation, or a vector of length 2, where the
      first element specifies the number of observations that should be sampled for each
      batch and the second argument specifies the number of batches, see the example.
    \item \code{nu}, the step length control parameter. Defaults to \code{nu = 0.05}.
      If argument \code{slice = TRUE} then \code{nu = 1}.
    \item \code{loglik}, defaults to \code{loglik = FALSE}. If set to \code{loglik = TRUE}
      the "out-of-sample" log-likelihood is used for smoothing variance estimation.
    \item \code{aic}, defaults to \code{aic = FALSE}, If set to \code{aic = TRUE}
      the "out-of-sample" AIC is used for smoothing variance estimation.
    \item \code{eps_loglik}, defaults to \code{eps_loglik = 0.01}. This argument specifies
      the relative change in the "out-of-sample" log-likelihood that is needed such that
      a model term gets updated.
    \item \code{select}, defaults to \code{select = FALSE}. If set to \code{select = TRUE},
      the algorithm only selects the model term with the largest contribution in the
      "out-of-sample" log-likelihood for updating in each iteration/batch.
    \item \code{always}, defaults to \code{always = FALSE}. If set to \code{always = TRUE}
      no log-likelihood contribution checks will be used and model terms are always updated.
    \item \code{K}, defaults to \code{K = 2}. This argument controls the penalty on the
      degrees of freedom in the computation of the AIC.
    \item \code{slice}, defaults to \code{slice = FALSE}. If set to \code{slice = TRUE},
      slice sampling using the "out-of-sample" log-likelihood or AIC is used for smoothing
      variance estimation. Moreover, \code{always = TRUE}, \code{eps_loglik = -Inf} and
      \code{nu = 1}. If \code{slice} is an integer \code{n}, slice sampling is started after
      \code{n} iterations, before smoothing variances are optimized.
  }

  When using function \code{opt_bbfitp}, the parameter updates are stored as \code{"mcmc"}
  objects. In this case the traceplots can be visualized using \code{\link{plot.bamlss}}.
}

\value{
  For function \code{opt_bbfit()} a list containing the following objects:
  \item{fitted.values}{A named list of the fitted values of the modeled parameters
    of the selected distribution.}
  \item{parameters}{The estimated set regression coefficients and smoothing variances.}
  \item{shuffle}{Logical}
  \item{runtime}{The runtime of the algorithm.}
}

\seealso{
\code{\link{bamlss}}, \code{\link{bfit}}
}

\examples{
\dontrun{## Simulate data.
set.seed(123)
d <- GAMart(n = 27000, sd = -1)

## Write data to disc.
tf <- tempdir()
write.table(d, file.path(tf, "d.raw"), quote = FALSE, row.names = FALSE, sep = ",")

## Model formula.
f <- list(
  num ~ s(x1,k=40) + s(x2,k=40) + s(x3,k=40) + te(lon,lat,k=10),
  sigma ~ s(x1,k=40) + s(x2,k=40) + s(x3,k=40) + te(lon,lat,k=10)
)

## Specify 50 batches with 1000 observations.
batch_ids <- c("nobs" = 1000, "nbatch" = 50)

## Note, can also be a list of indices, e.g.
## batch_ids <- lapply(1:50, function(i) { sample(1:nrow(d), size = 1000) })

## Different flavors:
## (1) Using "out-of-sample" aic for smoothing
##     variance estimation. Update is only accepted
##     if the "out-of-sample" log-likelihood is
##     increased. If data is a filepath, the data set is
##     read into R using package ff and model and
##     design matrices are processed with ff. This may
##     take some time depending on the size of the data.
set.seed(1)
b1 <- bamlss(f, data = file.path(tf, "d.raw"),
  sampler = FALSE, optimizer = opt_bbfit,
  batch_ids = batch_ids, nu = 0.1, aic = TRUE, eps_loglik = -Inf,
  always = FALSE)

## Plot estimated effects.
plot(b1)

## Plot coefficient paths for x3 in mu.
pathplot(b1, name = "mu.s.s(x3).b")

## (2) Same but always update, this mimics the classic SGD.
##     Note, for prediction only the last iteration is
##     used in this case. To use more iterations use opt_bbfitp(),
##     Then iterations are stored as "mcmc" object and we can
##     predict using the burnin argment, e.g.,
##     p <- predict(b2, model = "mu", burnin = 35)
set.seed(2)
b2 <- bamlss(f, data = file.path(tf, "d.raw"),
  sampler = FALSE, optimizer = opt_bbfit,
  batch_ids = batch_ids, nu = 0.1, aic = TRUE, eps_loglik = -Inf,
  always = TRUE)

## Plot coefficient paths for x3 in mu.
pathplot(b2, name = "mu.s.s(x3).b")

## (3) Boosting type flavor, only update model term with
##     the largest contribution in the "out-of-sample"
##     log-likelihood. In this case, if edf = 0 during
##     runtime of the algorithm, no model has an additional
##     contribution and the algorithm converges. This
##     behavior is controlled by argument eps_loglik, the
##     higher eps_loglik, the more restrictive is the
##     updating step.

## Initialize intercepts.
set.seed(0)
b0 <- bamlss(num ~ 1, data = d, sampler = FALSE, optimizer = opt_bbfitp,
  batch_ids = batch_ids, slice = TRUE)

## Compute starting values, remove the first
## 10 iterates and compute the mean of the
## remaining iterates.
start <- coef(b0, FUN = mean, burnin = 10)

## Start boosting, only update if change in
## "out-of-sample" log-likelihood is 0.1%, i.e.,
## eps_loglik = 0.001.
set.seed(3)
b3 <- bamlss(f, data = d, sampler = FALSE, optimizer = opt_bbfit,
    batch_ids = batch_ids, nu = 0.1, aic = TRUE, eps_loglik = 0.001,
    select = TRUE, always = FALSE, start = start)

## Plot log-likelihood contributions.
contribplot(b3)
## In this case, the algorithm did not converge,
## we need more iterations/batches.

## Note, prediction uses last iterate.
p3 <- predict(b3, model = "mu")

## (4) Use slice sampling under the "out-of-sample"
##     log likelihood for estimation of smoothing
##     variances. In this case model terms are always
##     updated ad parameter paths behave like a MCMC
##     chain. Therefore, use opt_bbfitp(), which stores
##     parameter paths as "mcmc" objects and we can
##     inspect using traceplots. Note nu = 1 if
##     slice = TRUE.
set.seed(4)
b4 <- bamlss(f, data = d, sampler = FALSE, optimizer = opt_bbfitp,
  batch_ids = batch_ids, aic = TRUE, slice = TRUE)

plot(b4)

## Plot parameter updates/samples.
plot(b4, which = "samples")

## Predict with burnin and compute mean
## prediction of the last 20 iterates.
p4 <- predict(b4, model = "mu", burnin = 30, FUN = mean)
}
}

\keyword{regression}

