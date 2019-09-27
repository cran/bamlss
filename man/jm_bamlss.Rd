\name{jm_bamlss}
\alias{jm_bamlss}
\alias{jm_mode}
\alias{jm_mcmc}
\alias{jm_transform}
\alias{jm_survplot}
\alias{jm_predict}

\title{
  Fit Flexible Additive Joint Models
}

\description{
  Family object to fit a flexible additive joint model for longitudinal and survival
  data under a Bayesian approach as presented in Koehler et al. (2017a, b). 
  All parts of the joint model can be specified as structured additive predictors. 
  See the details and examples.
}

\usage{
## JM family object.
jm_bamlss(...)

## "bamlss.frame" transformer function
## to set up joint models.
jm_transform(x, y, data, terms, knots, formula, family, subdivisions = 25,
  timedependent = c("lambda", "mu", "alpha", "dalpha"), timevar = NULL, 
  idvar = NULL, alpha = .Machine$double.eps, mu = NULL, sigma = NULL, 
  sparse = TRUE, nonlinear = FALSE, edf_alt = FALSE, start_mu = NULL, 
  k_mu = 6, ...)
                    
## Posterior mode optimizing engine.     
jm_mode(x, y, start = NULL, weights = NULL, offset = NULL, 
  criterion = c("AICc", "BIC", "AIC"), maxit = c(100, 1), 
  nu = c("lambda" = 0.1, "gamma" = 0.1, "mu" = 0.1, "sigma" = 0.1, 
  "alpha" = 0.1, "dalpha" = 0.1),
  update.nu = FALSE, eps = 0.0001, alpha.eps = 0.001, ic.eps = 1e-08, 
  nback = 40, verbose = TRUE, digits = 4, ...)

## Sampler function.
jm_mcmc(x, y, family, start = NULL, weights = NULL, offset = NULL, 
  n.iter = 1200, burnin = 200, thin = 1, verbose = TRUE, digits = 4, 
  step = 20, ...)
  
## Predict function, set to default in jm_bamlss().
jm_predict(object, newdata,
  type = c("link", "parameter", "probabilities", "cumhaz", "loglik"),
  dt, steps, id, FUN = function(x) { mean(x, na.rm = TRUE) }, 
  subdivisions = 100, cores = NULL, chunks = 1,
  verbose = FALSE,  ...)
  
## Survival plot.
jm_survplot(object, id = 1, dt = NULL, steps = 10, 
  points = TRUE, rug = !points)
}

\arguments{
  \item{x}{The \code{x} list, as returned from function
    \code{\link{bamlss.frame}} (and transformed by function \code{jm_transform()}),
    holding all model matrices and other information that is used for
    fitting the model.}
  \item{y}{The model response, as returned from function \code{\link{bamlss.frame}}.}
  \item{data}{A \code{\link{data.frame}} or \code{\link{list}} containing the model
    response variable(s) and covariates specified in the \code{formula} in long format.
    By default the variables are taken from \code{environment(formula)}:
    typically the environment from which \code{bamlss} is called.}
  \item{terms}{The corresponding \code{\link{terms.bamlss}} object needed for processing.}
  \item{knots}{An optional list containing user specified knots, see the documentation of
    function \code{\link[mgcv]{gam}}.}
  \item{formula}{The corresponding \code{\link{bamlss.formula}}.}
  \item{family}{The \code{\link{bamlss.family}} object.}
  \item{subdivisions}{How many time points should be created for each individual.}
  \item{timedependent}{A character vector specifying the names of parameters in \code{x}
    that are time-dependent. Time grid design matrices are only computed for these parameters.}
  \item{timevar}{A character specifying the name of the survival time variable in the
    data set.}
  \item{idvar}{Depending on the type of data set, this is the name of the variable specifying
     identifier of individuals.}
  \item{alpha}{Numeric, a starting value for the intercept of the association parameter alpha.}
  \item{mu}{Numeric, a starting value for the intercept of the mu parameter.}
  \item{sigma}{Numeric, a starting value for the intercept of the sigma parameter.}
  \item{sparse}{Logical, indicating if sparse matrix structures are used for updating and sampling
    of mu parameter model terms.}
  \item{nonlinear}{Logical, indicating if association is nonlinear in mu. See Details on the different
    model specifications.}
  \item{edf_alt}{Logical, indicating if an alternative computation of estimated degrees of freedom 
    for penalized model terms should be used.}
  \item{start_mu}{Starting values for the computation of mu. For estimating associations which are 
    nonlinear in mu, knot placement is based on these starting values which can improve stability. }
  \item{k_mu}{Number of knots for spline basis of association nonlinear in mu. Reducing this number
    improves stability of the estimation.}
  \item{start}{A named numeric vector containing possible starting values, the names are based on
    function \code{\link{parameters}}.}
  \item{weights}{Currently not supported.}
  \item{offset}{Currently not supported.}
  \item{criterion}{Information criterion to be used, e.g., for smoothing
    variance selection. Options are the corrected AIC \code{"AICc"} (see Details), the \code{"BIC"} and
    \code{"AIC"}. Defaults to \code{"AICc"}?}
  \item{maxit}{Vector containing the maximum number of iterations for the backfitting
    algorithm with \code{maxit[1]} defining the iterations for the full model and \code{maxit[2]}
    the iterations within each predictor. \code{maxit[2]} defaults to 1 if only one value is
    specified.}
  \item{nu}{Vector of step lengths for parameter updates of one Newton-Raphson update for 
    each predictor of the joint model.}
  \item{update.nu}{Should the updating step length be optimized in each iteration
    of the backfitting algorithm? Uses \code{nu} as starting value if set to \code{TRUE}.}
  \item{eps}{The relative convergence tolerance of the backfitting algorithm.}
  \item{alpha.eps}{The relative convergence tolerance of the backfitting algorithm
    for predictor alpha.}
  \item{ic.eps}{The relative convergence tolerance of the information criterion used, e.g.,
    for smoothing variance selection.}
  \item{nback}{For computing \code{ic.eps}, how many iterations back should be included
    when computing relative convergence tolerance of the information criterion.}
  \item{verbose}{Print information during runtime of the algorithm.}
  \item{digits}{Set the digits for printing when \code{verbose = TRUE}.}
  \item{n.iter}{the number of MCMC iterations.}
  \item{burnin}{the burn-in phase of the sampler, i.e., the number of starting samples that
    should be removed.}
  \item{thin}{the thinning parameter for MCMC simulation. E.g., \code{thin = 10} means,
    that only every 10th sampled parameter will be stored.}
  \item{step}{How many times should algorithm runtime information be printed, divides \code{n.iter}.}
  \item{object}{A \code{"bamlss"} object processed with the JM optimizer function
    \code{jm_mode()} ans/or sampler function \code{jm_mcmc()} for which the survival plot
    should be created.}
  \item{newdata}{Dataset for which to create predictions. Not needed for conditional survival probabilities.}
  \item{type}{Character string indicating which type of predictions to compute. \code{link} returns estimates
    for all predictors with the respective link functions applied, \code{"parameter"} returns the estimates
    for all pedictors, \code{"probabilities"} returns the survival probabilities conditional on the
    survival up to the last longitudinal measurement, and \code{"cumhaz"} return the cumulative hazard
    up to the survival time or for a time window after the last longitudinal measurement. If \code{type}
    is set to \code{"loglik"}, the log-likelihood of the joint model is returned.}
  \item{id}{Integer or character, that specifies the individual for which the plot should be created.}
  \item{dt}{The time window after the last observed measurement for which predictions should be computed. 
    The default is \code{0.4 * max(obstime)} and \code{obstime} are the individual's longitudinal measurement times.}
  \item{steps}{Integer, the number of steps for which to evaluate the conditional survival probability 
    up to \code{dt}.}
  \item{FUN}{A function that should be applied on the samples of predictors or
    parameters, depending on argument \code{type}.}
  \item{cores}{Specifies the number of cores that should be used for prediction. Note that
    this functionality is based on the \code{\link[parallel]{parallel}} package.}
  \item{chunks}{Should computations be split into \code{chunks}? Prediction is then processed
    sequentially.}
  \item{points}{Should longitudinal observations be added to the plot.}
  \item{rug}{Should longitudinal observed time points be added on the x-axis to the plot.}
  \item{\dots}{Currently not used.}
}

\details{
We refer to the papers of Koehler et al. (2017a, b) for details on the flexible
additive joint model. In short, we model the hazard of subject \eqn{i} an event at time
\eqn{t} as
\deqn{h_{i}(t)= \exp [\eta_{\lambda i}(t)+
\eta_{\gamma i}+\eta_{\alpha i}(\eta_{\mu i}(t), t) ]}
with predictor \eqn{\eta_{\lambda}} for all survival covariates that are time-varying or have a
time-varying coefficient (including the log baseline hazard), predictor \eqn{\eta_{\gamma}} for
baseline survival covariates, predictor \eqn{\eta_{\alpha}} representing the potentially
time-varying or nonlinear association between the longitudinal marker \eqn{\eta_{\mu}} and the hazard. The
longitudinal response \eqn{y_{ij}} at time points \eqn{t_{ij}} is modeled as
\deqn{y_{ij}=\eta_{\mu i}(t_{ij})+e_{ij}} with independent normal errors \eqn{N(0, \exp[\eta_{\sigma i}(t_{ij})]^2)}. 

Each predictor \eqn{\eta_{ki}} is a structured additive predictor, i.e. a sum of functions of
covariates \eqn{\eta_{ki} =  \sum_{m=1}^{M_k} f_{km}(\bm{x}_{ki})}. Each of these functions can be
modeled parametrically or using basis function evaluations from the smooth constructors in
\pkg{mgcv} such as \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}} and \code{\link[mgcv]{ti}} and can
include smooth time-varying, random or spatial effects. For the Bayesian estimation of these effects
we specify corresponding priors: For linear or parametric terms we use vague normal priors, smooth
and random effect terms are regularized by placing generic multivariate normal priors on the
coefficients and for anisotropic smooths, when multiple smoothing variance parameters are involved,
more complex prior are in place (cf. Koehler et al., 2017a). We use inverse Gamma
hyper-priors, i.e. IG(0.001, 0.001) to obtain an inverse Gamma full conditional for the variance
parameters. We estimate the posterior mode by maximizing the log-posterior of the model using a
Newton-Raphson procedure, the posterior mean is obtained via derivative-based Metropolis-Hastings
sampling. We recommend to use posterior mode estimates for a quick model assessment. In order to
draw correct inferences from the model, posterior mean estimates should be computed. 
We approximate integration in the survival part of the likelihood using trapezoidal rule. For
posterior mode estimation.

A variety specifications of the association \eqn{\eta_{\alpha i}(\eta_{\mu i}(t), t)} are possible with an 
important distinction between associations which are nonlinear in \eqn{\eta_{\mu}}  for \code{nonlinear = TRUE} (Koehler et al. 2017b) or linear where \eqn{\eta_{\alpha i}(\eta_{\mu i}(t), t) = \eta_{\alpha i}(t)\eta_{\mu i}(t)} for \code{nonlinear = FALSE} (Koehler et al. 2017a).
}

\note{The indicator \code{nonlinear} for associations with are linear or nonlinear in \eqn{\eta_{\mu}} was named \code{interaction]} in earlier versions stages of the development.}

\references{
  Koehler M, Umlauf N, Beyerlein, A., Winkler, C. Ziegler, A.-G., Greven S (2017). Flexible
  Bayesian Additive Joint Models with an Application to Type 1 Diabetes Research.
  \emph{Biometrical Journal}, \url{http://dx.doi.org/10.1002/bimj.201600224}.
  
  Meike Koehler, Nikolaus Umlauf, and Sonja Greven (2018). Nonlinear association structures in
  flexible Bayesian additive joint models. \emph{Statistics in Medicine},
  \url{http://dx.doi.org/10.1002/sim.7967}.
}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}}.
}

\examples{
\dontrun{

set.seed(123)
## Simulate survival data
## with random intercepts/slopes and a linear effect of time,
## constant association alpha and no effect of the derivative
d <- simJM(nsub = 200, long_setting = "linear", 
  alpha_setting = "constant",
  dalpha_setting = "zero", full = FALSE)

## Formula of the according joint model
f <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime, bs = "ps"),
  gamma ~ s(x1, bs = "ps"),
  mu ~ obstime + s(id, bs = "re") + 
    s(id, obstime, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

## Joint model estimation
## jm_bamlss() sets the default optimizer and sampler function.
## First, posterior mode estimates are computed using function
## jm_mode(), afterwards the sampler jm_mcmc() is started.
b <- bamlss(f, data = d, family = "jm",
           timevar = "obstime", idvar = "id")

## Plot estimated effects.
plot(b)

## Predict event probabilities for two individuals
## at 12 time units after their last longitudinal measurement.
## The event probability is conditional on their survival
## up to their last observed measurement.
p <- predict(b, type = "probabilities", id = c(1, 2), dt = 12, FUN = c95)
print(p)

## Plot of survival probabilities and
## corresponding longitudinal effects
## for individual id.
jm_survplot(b, id = 3)
jm_survplot(b, id = 30)


## Simulate survival data
## with functional random intercepts and a nonlinear effect 
## of time, time-varying association alpha and no effect 
## of the derivative.
## Note: This specification is the simJM default.
d <- simJM(nsub = 200, full = FALSE)

## Formula of the according joint model
## number of knots for the smooth nonlinear effect of time
klong <- 8
f <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime, bs = "ps"),
  gamma ~ s(x1, bs = "ps"),
  mu ~ ti(id, bs = "re") + 
    ti(obstime, bs = "ps", k = klong) + 
    ti(id, obstime, bs = c("re", "ps"), 
      k = c(nlevels(d$id), klong)) + 
    s(x2, bs = "ps"),
  sigma ~ 1,
  alpha ~ s(survtime, bs = "ps"),
  dalpha ~ -1
)

## Estimating posterior mode only using jm_mode()
b_mode <- bamlss(f, data = d, family = "jm",
                 timevar = "obstime", idvar = "id", 
                 sampler = FALSE)

## Estimating posterior means using jm_mcmc() 
##  with starting values generated from posterior mode 
b_mean <- bamlss(f, data = d, family = "jm",
                 timevar = "obstime", idvar = "id", optimizer = FALSE,
                 start = parameters(b_mode), results = FALSE)

## Plot effects.
plot(b_mean, model = "alpha")



## Simulate survival data
## with functional random intercepts and an association nonlinear in mu
set.seed(234)
d <- simJM(nsub = 300, long_setting = "functional", alpha_setting = "nonlinear", 
           nonlinear = TRUE, full = FALSE, probmiss = 0.9)

## Calculate longitudinal model to obtain starting values for mu
long_df <- 7
f_start <- y ~ ti(id, bs = "re") + ti(obstime, bs = "ps", k = long_df) + 
               ti(id, obstime, bs = c("re", "ps"), k = c(nlevels(d$id), long_df)) + 
               s(x2, bs = "ps")
b_start <- bamlss(f_start, data = d, sampler = FALSE)
mu <- predict(b_start)$mu

## Fit joint model with nonlinear association (nonlinear = TRUE)
f <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime, bs = "ps"),
  gamma ~ x1,
  mu ~ ti(id, bs = "re") + ti(obstime, bs = "ps", k = long_df) + 
       ti(id, obstime, bs = c("re", "ps"), k = c(nlevels(d$id), long_df)) + 
       s(x2, bs = "ps"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
b <- bamlss(f, data = d, family = "jm", timevar = "obstime", idvar = "id", 
            nonlinear = TRUE, start_mu = mu, 
            n.iter = 6000, burnin = 2000, thin = 2)
plot(b)
samplestats(b$samples)
}
}

\keyword{regression}
\keyword{models}

