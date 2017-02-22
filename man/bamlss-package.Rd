\name{bamlss-package}
\alias{bamlss-package}
\docType{package}

\title{
  Bayesian Additive Models for Location Scale and Shape (and Beyond)
}

\description{
  The \pkg{bamlss} package is a general tool for complex Bayesian regression modeling
  with structured additive predictors based on Markov chain Monte Carlo simulation.
  The design of this package substantially focuses on maximum flexibility and easy
  integration of new code and/or standalone systems. The package makes heavy use of
  \code{\link{mgcv}} infrastructures to build up all necessary model matrices and
  information from which it is relatively easy for the user to construct estimation
  algorithms or interfaces to existing software packages. The package can also be
  seen as an harmonized framework for regression modeling since it does not restrict
  to any type of problem. The main function in this package is \code{\link{bamlss}},
  which is a wrapper function that calls optimizer and/or sampling functions for fitting
  Bayesian additive models for location scale and shape (and beyond). These model fitting
  functions can be exchanged by the user. Moreover, the package contains numerous
  functions for creating post-estimation results like summary statistics and
  effect plots etc.
}

\author{
  Nikolaus Umlauf, Nadja Klein, Achim Zeileis.
}

\keyword{package}

\seealso{
  \code{\link{bamlss}}, \code{\link{bamlss.frame}}
}

\keyword{regression}
\keyword{MCMC}
