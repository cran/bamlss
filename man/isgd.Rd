\name{isgd}
\alias{isgd}
\title{
Implicit Stochastic Gradient Descent Optimizer
}
\description{
This optimizer performs an implicit stochastic gradient descent
algorithm. It is mainly used within a \code{\link{bamlss}} call.
}
\usage{
isgd(x, y, family, weights = NULL, offset = NULL,
     gammaFun = function(i) 1/(1 + i), shuffle = TRUE,
     CFun = function(beta) diag(length(beta)),
     start = NULL, i.state = 0)
}
\arguments{
  \item{x}{For function \code{boost()} the \code{x} list, as returned from
           function \code{\link{bamlss.frame}}, holding all model matrices
           and other information that is used for fitting the model.}
  \item{y}{The model response, as returned from function
           \code{\link{bamlss.frame}}.}
  \item{family}{A \pkg{bamlss} family object, see
                \code{\link{family.bamlss}}.}
  \item{weights}{Prior weights on the data, as returned from function \code{\link{bamlss.frame}}.}
  \item{offset}{Can be used to supply model offsets for use in fitting,
    returned from function \code{\link{bamlss.frame}}.}
  \item{gammaFun}{Function specifying the step length.}
  \item{shuffle}{Should the data be shuffled?}
  \item{CFun}{Hessian approximating function.}
  \item{start}{Vector of starting values.}
  \item{i.state}{Added to \code{gammaFUN()}.}
}
\details{
tpf
}
\value{
  For function \code{isgd()} a list containing the following objects:
  \item{fitted.values}{A named list of the fitted values based on the last iteration
    of the modeled parameters of the selected distribution.}
  \item{parameters}{A matrix, each row corresponds to the parameter values of one iteration.}
  \item{sgd.summary}{The summary of the stochastic gradient descent algorithm which can be printed and plotted.}
}
\references{
Toulis, P and Airoldi, EM (2015): Scalable estimation strategies based on
stochastic approximations: Classical results and new insights.
Statistics and Computing, 25, no. 4, 781--795.
doi: 10.1007/s11222-015-9560-y
}
\author{
Thorsten Simon
}
\note{
Motivated by the lecture 'Regression modelling with large data sets'
given by Ioannis Kosmidis in Innsbruck, January 2017.
}

\section{Warning }{
CAUTION: Arguments \code{weights} and \code{offset} are not implemented yet!
}

\seealso{
\code{\link{bamlss.frame}}, \code{\link{bamlss}}
}
\examples{
\dontrun{
set.seed(111)
d <- GAMart(n = 10000)
f <- num ~ s(x1) + s(x2) + s(x3) + te(lon, lat)
b <- bamlss(f, data = d, optimizer = bamlss:::isgd, sampler = FALSE)
plot(b, ask = F)

## loop over observations a 2nd time
b <- bamlss(f, data = d, optimizer = bamlss:::isgd, sampler = FALSE, start = parameters(b),
            i.state = b$model.stats$optimizer$sgd.summary$i.state)
plot(b, ask = F)

## try differeent gammaFuns, e.g.,
# gammaFun <- function(i) .3/sqrt((1+i)) + 0.001

## testing some families
f2 <- bin ~ s(x1) + s(x2) + s(x3) + te(lon, lat)
b2 <- bamlss(f2, data = d, optimizer = isgd, sampler = FALSE, family = "binomial")

f3 <- cens ~ s(x1) + s(x2) + s(x3) + te(lon, lat)
b3 <- bamlss(f3, data = d, optimizer = isgd, sampler = FALSE, family = "cnorm")
}
}
\keyword{regression}
