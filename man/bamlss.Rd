\name{bamlss}
\alias{bamlss}

\title{
  Fit Bayesian Additive Models for Location Scale and Shape (and Beyond)
}

\description{
  This is the main model fitting function of the package. Function \code{bamlss()}
  is a wrapper function that parses the \code{data} and the model \code{formula}, or
  extended \code{\link{bamlss.formula}}, as well as the \code{\link{bamlss.family}}
  into a \code{\link{bamlss.frame}}. The \code{\link{bamlss.frame}} then holds all model
  matrices and information that is needed for setting up estimation engines.
  The model matrices are based on \code{\link{mgcv}} infrastructures, i.e.,
  smooth terms are constructed using \code{\link[mgcv]{smooth.construct}} and
  \code{\link[mgcv]{smoothCon}}. Therefore, all \code{\link{mgcv}} model term constructors like
  \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}}, \code{\link[mgcv]{t2}} and \code{\link[mgcv]{ti}}
  can be used. Identifiability conditions are imposed using function \code{\link[mgcv]{gam.side}}.

  After the \code{\link{bamlss.frame}} is set up function \code{bamlss()} applies optimizer
  and/or sampling functions. These functions can also be provided by the user. See the details
  below on how to create new engines to be used with function \code{bamlss()}.

  Finally, the estimated parameters and/or samples are used to create model output results like
  summary statistics or effect plots. The computation of results may also be controlled by the user.
}

\usage{
bamlss(formula, family = "gaussian", data = NULL,
  start = NULL, knots = NULL, weights = NULL,
  subset = NULL, offset = NULL, na.action = na.omit,
  contrasts = NULL, reference = NULL, transform = NULL,
  optimizer = NULL, sampler = NULL, samplestats = NULL,
  results = NULL, cores = NULL, sleep = NULL,
  combine = TRUE, model = TRUE, x = TRUE,
  light = FALSE, ...)
}

\arguments{
  \item{formula}{A formula or extended formula, i.e., the \code{formula} can be a
    \code{\link{list}} of formulas where each list entry specifies the details of one parameter
    of the modeled response distribution, see \code{\link{bamlss.formula}}. For incorporating
    smooth terms, all model term constructors implemented in \code{\link{mgcv}} such as
    \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}} and \code{\link[mgcv]{ti}} can be used, amongst others.}
  \item{family}{A \code{\link{bamlss.family}} object, specifying the details of the modeled
    distribution such as the parameter names, the density function, link functions, etc.
    Can be a character without the \code{"_bamlss"} extension of the
    \code{\link{bamlss.family}} name.}
  \item{data}{A \code{\link{data.frame}} or \code{\link{list}} containing the model
    response variable(s) and covariates specified in the \code{formula}.
    By default the variables are taken from \code{environment(formula)}:
    typically the environment from which \code{bamlss} is called.}
  \item{start}{A named numeric vector containing starting values to be send to the \code{optimizer}
    and/or \code{sampler} function. For a possible naming convention for the parameters see
    function \code{\link{parameters}}, but this is not restrictive and engine specific.}
  \item{knots}{An optional list containing user specified knots, see the documentation of
    function \code{\link[mgcv]{gam}}.}
  \item{weights}{Prior weights on the data.}
  \item{subset}{An optional vector specifying a subset of observations to be
    used in the fitting process.}
  \item{offset}{Can be used to supply model offsets for use in fitting.}
  \item{na.action}{A function which indicates what should happen when the data
    contain \code{NA}'s. The default is set by the \code{na.action}
    setting of \code{\link{options}}, and is \code{\link{na.omit}} if set to \code{NULL}.}
  \item{contrasts}{An optional list. See the \code{contrasts.arg} of 
    \code{\link[stats]{model.matrix.default}}.}
  \item{reference}{A \code{character} specifying a reference category, e.g., when
    fitting a multinomial model.}
  \item{transform}{A transformer function that is applied on the \code{\link{bamlss.frame}}.
    See, e.g., function \code{\link{randomize}} and \code{\link{bamlss.engine.setup}}.}
  \item{optimizer}{An optimizer function that returns, e.g., posterior mode estimates
    of the parameters as a named numeric vector. The default optimizer function is
    \code{\link{bfit}}. If set to \code{FALSE}, no optimizer function will be used.}
  \item{sampler}{A sampler function that returns a matrix of samples, the columns represent the
    parameters, the rows the iterations. The returned matrix must be coerced to an object of
    class \code{"mcmc"}, see \code{\link[coda]{as.mcmc}}. The default sampler function is
    \code{\link{GMCMC}}. If set to \code{FALSE}, no sampler function will be used.}
  \item{samplestats}{A function computing statistics from samples, per default function
    \code{\link{samplestats}} is used. If set to \code{FALSE}, no \code{samplestats} function
    will be used. Note that this option is crucial for very large datasets, as computing
    statistics from samples this way may be very time consuming!}
  \item{results}{A function computing results from the parameters and/or samples, e.g., for
    creating effect plots, see function \code{link{results.bamlss.default}}. If set \code{FALSE}
    no \code{results} function will be used.}
  \item{cores}{An integer specifying the number of cores that should be used for the sampler
    function. This is based on function \code{\link[parallel]{mclapply}} of the \pkg{parallel}
    package.}
  \item{sleep}{Time the system should sleep before the next core is started.}
  \item{combine}{If samples are computed on multiple cores, should the samples be combined into
    one \code{\link[coda]{mcmc}} matrix?}
  \item{model}{If set to \code{FALSE} the model frame used for modeling is not part of the
    return value.}
  \item{x}{If set to \code{FALSE} the model matrices are not part of the return value.}
  \item{light}{Should the returned object be lighter, i.e., if \code{light = TRUE} the returned
    object will not contain the model.frame and design and penalty matrices are deleted.}
  \item{\dots}{Arguments passed to the \code{transformer}, \code{optimizer}, \code{sampler},
    \code{results} and \code{samplestats} function.}
}

\details{
  The main idea of this function is to provide infrastructures that make it relatively easy
  to create estimation engines for new problems, or write interfaces to existing software
  packages.

  The steps that are performed within the function are:
  \itemize{
    \item First, the function parses the \code{data}, the \code{formula} or the extended
      \code{\link{bamlss.formula}} as well as the \code{\link{bamlss.family}} into a model frame
      like object, the \code{\link{bamlss.frame}}. This object holds all necessary model matrices
      and information that is needed for subsequent model fitting engines. Per default,
      all package \code{\link{mgcv}} smooth term constructor functions like
      \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}}, \code{\link[mgcv]{t2}} and
      \code{\link[mgcv]{ti}} can be used (see also function \code{\link[mgcv]{smooth.construct}}),
      however, even special user defined constructors can be included, see the manual of
      \code{\link{bamlss.frame}}.
    \item In a second step, the \code{\link{bamlss.frame}} can be transformed, e.g., if a mixed
      model representation of smooth terms is needed, see function \code{\link{randomize}}.
    \item Then an optimizer function is started, e.g., a function that finds posterior mode estimates
      of the parameters. A convention for model fitting engines is that such functions should have
      the following arguments:

      \code{optimizer(x, y, family, start, weights, offset, ...)}

      Internally, function \code{bamlss()} will send the \code{x} object that holds all
      model matrices, the response \code{y} object, the \code{family} object, \code{start}ing
      values for the parameters, possible \code{weights} and \code{offset}s of the created
      \code{\link{bamlss.frame}} to the
      optimizer function (see the manual of \code{\link{bamlss.frame}} for more details on the
      \code{x}, \code{y} and other objects). The job of the optimizer is to return a named numeric
      vector of optimum parameters. The names of the parameters should be such that they can be
      uniquely mapped to the corresponding model matrices in \code{x}. See function
      \code{\link{parameters}} for more details on parameter names. The default optimizer function
      is \code{\link{bfit}}. The optimizer can return more information than only the optimum
      parameters. It is possible to return a list, the convention here is that an element named
      \code{"parameters"} then holds the named vector of estimated parameters. Possible other return
      values could be fitted values, the Hessian matrix, information criteria or information
      about convergence of the algorithm, etc. Note that the parameters are simply added to the
      \code{\link{bamlss.frame}} in an (list) entry named \code{parameters}.
    \item After the optimization step, a \code{sampler} function is started. The arguments of such
      sampler functions are the same as for the \code{optimizer} functions

      \code{sampler(x, y, family, start, weights, offset, ...)}

      Sampler functions must return a matrix of samples, each row represents one iteration and the matrix
      can be coerced to \code{\link[coda]{mcmc}} objects. The function may return a list of samples,
      e.g., if multiple chains are returned each list entry then holds one sample matrix of
      one chain. The column names of the sample matrix should be the same as the names of estimated
      parameters. For a possible naming convention see function \code{\link{parameters}}, which
      ensures unique mapping of samples with the model matrices in the \code{x} object of the
      \code{\link{bamlss.frame}}. The samples are added to the \code{\link{bamlss.frame}}
      in an (list) entry named \code{samples}.
    \item Next, the \code{samplestats} function is applied. This function can compute any quantity
      from the samples and the \code{x} object, the arguments of such functions are

      \code{samplestats(samples, x, y, family, ...)}

      where argument \code{samples} are the samples returned from the \code{sampler} function,
      and \code{x}, \code{y} and \code{family} are the same objects as passed to the optimizer
      and or sampler functions. For example, the default function in \code{bamlss()} for this task
      is also called \code{\link{samplestats}} and returns the mean of the log-likelihood and the
      log-posterior computed of all samples, as well as the DIC.
    \item The last step is to compute more complex information about the model using the
      \code{results} function. The arguments of such \code{results} functions are

      \code{results(bamlss.frame, ...)}

      here, the full \code{\link{bamlss.frame}} including possible \code{parameters} and
      \code{samples} is passed to the function within \code{bamlss()}. The default function
      for this task is \code{\link{results.bamlss.default}} which returns an object of class
      \code{"bamlss.results"} for which generic plotting functions are and a \code{\link{summary}}
      function is provided. Hence, the user can control the output of the model, the plotting
      and summary statistics, too.
  }

  Note that function \code{transform()}, \code{optimizer()}, \code{sampler()}, \code{samplestats()}
  and \code{results()} can be provided from the \code{\link{bamlss.family}} object, e.g.,
  if a \code{\link{bamlss.family}} object has an element named \code{"optimizer"}, which
  represents a valid optimizer function such as \code{\link{bfit}}, exactly this optimizer
  function will be used as a default when using the family.
}

\value{
  An object of class \code{"bamlss"}. The object is in principle only a slight extension
  of a \code{\link{bamlss.frame}}, i.e., if an \code{optimizer} is applied it will hold the
  estimated parameters in an additional element named \code{"parameters"}. If a sampler function
  is applied it will additionally hold the samples in an element named \code{"samples"}.
  The same mechanism is used for \code{results} function.

  If the \code{optimizer} function computes additional output next to the parameters, this will
  be saved in an element named \code{"model.stats"}. If a \code{samplestats} function is applied,
  the output will also be saved in the \code{"model.stats"} element.

  Additionally, all functions that are called are saved as attribute \code{"functions"} in the
  returned object.
}

\references{
  Umlauf N, Klein N, Zeileis A (2016). Bayesian Additive Models for Location
    Scale and Shape (and Beyond). \emph{(to appear)}
}

\author{
  Nikolaus Umlauf, Nadja Klein, Achim Zeileis.
}

\seealso{
  \code{\link{bamlss.frame}}, \code{\link{family.bamlss}}, \code{\link{bamlss.formula}},
  \code{\link{randomize}}, \code{\link{bamlss.engine.setup}},
  \code{\link{bfit}}, \code{\link{GMCMC}}, \code{\link{continue}},
  \code{\link{coef.bamlss}}, \code{\link{parameters}}, \code{\link{predict.bamlss}},
  \code{\link{plot.bamlss}}
}

\examples{
\dontrun{## Simulated data example.
d <- GAMart()
f <- num ~ s(x1) + s(x2) + s(x3) + te(lon, lat)
b <- bamlss(f, data = d)
summary(b)
plot(b)
plot(b, which = 3:4)
plot(b, which = "samples")

## Use of optimizer and sampler functions:
## * first run optimizer,
b1 <- bamlss(f, data = d, optimizer = bfit, sampler = FALSE)
print(b1)
summary(b1)

## * afterwards, start sampler with staring values,
b2 <- bamlss(f, data = d, start = coef(b1), optimizer = FALSE, sampler = GMCMC)
print(b2)
summary(b2)

## Continue sampling.
b3 <- continue(b2, n.iter = 12000, burnin = 0, thin = 10)
plot(b3, which = "samples")
plot(b3, which = "max-acf")
plot(b3, which = "max-acf", burnin = 500, thin = 4)
}
}

\keyword{models}
\keyword{regression}

