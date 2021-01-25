\name{simJM}
\alias{simJM}
\alias{rJM}

\title{
  Simulate longitudinal and survival data for joint models
}

\description{
  Simulates longitudinal data with normal error and (Cox-type) survival times
  using the inversion method. The function \code{simJM()} is a wrapper specifying
  all predictors and the resulting data sets. The wrapper calls \code{rJM()} to sample
  the survival times, a modified version of \code{rSurvtime()} from the \R package
  \pkg{CoxFlexBoost}.
}

\usage{
simJM(nsub = 300, times = seq(0, 120, 1), probmiss = 0.75,
  long_setting = "functional",
  alpha_setting = if(nonlinear) "linear" else "nonlinear",
  dalpha_setting = "zero", sigma = 0.3, long_df = 6, tmax = NULL,
  seed = NULL, full = FALSE, file = NULL, nonlinear = FALSE,
  fac = FALSE)
  
rJM(hazard, censoring, x, r, 
  subdivisions = 1000, tmin = 0, tmax,
  file = NULL, ...)
}

\arguments{
  \item{nsub}{number of individuals for which longitudinal data and survival times should be simulated.}
  \item{times}{vector of time points at which longitudinal measurements are "sampled".}
  \item{probmiss}{proportion of longitudinal measurements to be set to missing. Used to induce
    sparsity in the longitudinal measurements.}
  \item{long_setting}{Specification of the longitudinal trajectories of the sampled subjects.
    Preset specifications are \code{"linear"}, \code{"nonlinear"} and \code{"functional"}. See Details.}
  \item{alpha_setting}{specification of the association between survival and longitudinal. Preset
    specifications are \code{"simple"}, \code{"linear"}, \code{"nonlinear"} and
    \code{"nonlinear2"}. See Details.}
  \item{dalpha_setting}{specification of the association between survival and the derivative of the
    longitudinal. Work in progress.}
  \item{sigma}{standard deviation of the normal error around the true longitudinal measurements.}
  \item{long_df}{number of basis functions from which functional random intercepts are sampled.}
  \item{tmax}{For function \code{simJM()}, longest possible survival time, observations are censored
    after that timepoint. Defaults to \code{max(times)} and should not be specified longer than
    \code{max(times)} for longitudinal setting "functional". For function \code{rJM()}, latest
    time point to sample a survival time.}
  \item{seed}{numeric scalar setting the random seed.}
  \item{full}{logical indicating if only the longitudinal data set should be returned (\code{FALSE})
    or additionally also the data for the survival part evaluated on a regular time grid and the
    longitudinal data set without longitudinal missings (\code{TRUE}). }
  \item{file}{name of the data file the generated data set should be stored into
    (e.g., "simdata.RData") or NULL if the dataset should directly be returned in \R.}
  \item{nonlinear}{If set to \code{TRUE}, a nonlinear interaction between \code{alpha} and
    \code{mu} is simulated.}
  \item{fac}{If set to \code{TRUE}, a smooth interaction that varies by a factor is simulated.}
  \item{hazard}{complete hazard function to specify the joint model. Time must be the
    first argument.}
  \item{censoring}{function to compute (random) censoring. }
  \item{x}{matrix of sampled covariate values.}
  \item{r}{matrix of sampled random coefficients.}
  \item{subdivisions}{the maximum number of subintervals for the integration.}
  \item{tmin}{earliest time point to sample a survival time.}
  \item{\dots}{further arguments to be passed to \code{hazard} or \code{censoring}.}
}

\details{
  The function simulates longitudinal data basing on the given specification at given \code{times}.
  The full hazard is built from all joint model predictors \eqn{\eta_{\mu}}, \eqn{\eta_{\sigma}},
  \eqn{\eta_{\lambda}}, \eqn{\eta_{\gamma}}, \eqn{\eta_{\alpha}} as presented in
  Koehler, Umlauf, and Greven (2016), see also {\code{\link{jm_bamlss}}}. Survival times are sampled using the inversion
  method (cf. Bender, Augustin, & Blettner, 2005). Additional censoring and missingness is
  introduced. The longitudinal information is censored according to the survival information. The
  user can also specify own predictors and use only \code{rJM} to simulate survival times
  accordingly.  
  
  Pre-specified functions for \eqn{\eta_{\mu}} in \code{long_setting} are for \code{linear}
  \deqn{\eta_{\mu i}(t) = 1.25 + r_{1i} + 0.6 \sin(x_{2i}) + (-0.01) t + 0.02 r_{2i} t},
  for \code{nonlinear} \deqn{\eta_{\mu i}(t) = 0.5 + r_{1i} + 0.6 \sin(x_{2i}) + 0.1 (t+1) \exp(-0.075 t)}
  and for \code{functional}
  \deqn{\eta_{\mu i}(t) = 0.5 + r_{1i} + 0.6 \sin(x_{2i}) + 0.1 (t+1) \exp(-0.075 t) + \sum_k \beta_{ki} B(t)},
  where \eqn{B(.)} denotes a B-spline basis function and \eqn{\beta_{ki}} are the sampled penalized
  coefficients from \code{gen_b} per person. 

  Prespecified functions for \eqn{\eta_{\alpha}} in \code{alpha_setting} are for \code{constant}
  \deqn{\eta_{\alpha}(t) = 1}, for \code{linear} \deqn{\eta_{\alpha}(t) = 1 - 0.015 t}, for
  \code{nonlinear} \deqn{\eta_{\alpha}(t) = \cos((time-20)/20)}, and for \code{nonlinear}
  \deqn{\eta_{\alpha}(t) = \cos((time-33)/33)}. 

  Additionally the fixed functions for \eqn{\eta_{\lambda} = 0.1(t+2)\exp(-0.075t)} and
  \eqn{\eta_{\lambda} = 0.1(t+2)\exp(-0.075t)} are employed.
}


\value{
  For \code{full = TRUE} a list of the three \code{data.frame}s is returned:
    \item{data}{Simulated dataset in long format including all longitudinal and survival covariates.}
    \item{data_grid}{Dataset of the time-varying survival predictors which are not subject specific,
      evaluated at a grid of fixed time points.}
    \item{data_full}{Simulated data set prior to generating longitudinal missings. Useful to
      assess the longitudinal fit.}
  For \code{full = FALSE} only the first dataset is returned. 

  Covariates within these datasets include a subject identifier \code{id}, the sampled survival
  times \code{survtime}, the event indicator \code{event}, the time points of longitudinally
  "observed" measurements  \code{obstime},   the longitudinal response \code{y}, the cumulative
  hazard at the survival time \code{cumhaz}, as well as covariates  \code{x1, x2}, random effects
  \code{r1, r2, b1, ...},  and the true predictors \code{alpha, lambda, gamma, mu, sigma}. 
}

\references{
  Hofner, B (2016). \pkg{CoxFlexBoost}: Boosting Flexible Cox Models
  (with Time-Varying Effects). \R package version 0.7-0.

  Bender, R., Augustin, T., and Blettner, M. (2005).
  Generating Survival Times to Simulate Cox Proportional Hazards Models.
  \emph{Statistics in Medicine}, \bold{24}, 1713-1723.
  
   Koehler N, Umlauf N, Beyerlein, A., Winkler, C., Ziegler, A., and Greven S (2016). Flexible Bayesian Additive Joint Models with an
    Application to Type 1 Diabetes Research. \emph{(submitted)}
}

\seealso{
  \code{\link{jm_bamlss}}, \code{\link{opt_JM}}, \code{\link{sam_JM}}, \code{\link{bamlss}}.
}

\examples{
\dontrun{## Simulate survival data
## with functional random intercepts and a nonlinear effect 
## of time, time-varying association alpha.
d <- simJM(nsub = 300)
head(d)

## Simulate survival data
## with random intercepts/slopes and a linear effect of time,
## constant association alpha.
d <- simJM(nsub = 200, long_setting = "linear", 
  alpha_setting = "constant")
head(d)
}
}

\keyword{survival}
\keyword{distribution}

