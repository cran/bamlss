\docType{data}
\name{simdata}
\alias{simdata}
\title{Reference data.}
\format{
An object of class \code{list} of length 3.
}
\usage{
data("simdata")
}
\description{
Simulated data to test the implementation of the bamlss families.
}

\seealso{
  \code{\link{mvnchol_bamlss}}
}

\examples{
\dontrun{## Reproducing code.
set.seed(111)
n <- 2000

## build orthogonal rotation matrix
thetax <- pi/4
thetay <- pi/4
thetaz <- pi/4
Rx <- matrix( c(1,0,0, 0,cos(thetax),sin(thetax), 0,-sin(thetax),cos(thetax) ), 3, 3 )
Ry <- matrix( c(cos(thetay),0,-sin(thetay), 0,1,0,  sin(thetay),0,cos(thetay) ), 3, 3 )
Rz <- matrix( c(cos(thetaz),sin(thetaz),0, -sin(thetaz),cos(thetaz),0, 0,0,1 ), 3, 3 )
R <- Rx \%*\% Ry \%*\% Rz

## non-linear functions
f1 <- function(x) (sin(pi * x))^2
f2 <- function(x) (cos(pi * x))^2

## random derivitates  
x <- runif(n)

## eigenvalues
val1 <- f1(x) 
val2 <- f2(x)
val3 <- rep(0, n)

## initialize vectors for parameter lists
p12 <- NULL
p13 <- NULL
p23 <- NULL
sig <- matrix(0, n, 3)

lamdiag <- matrix(0, n, 3)
lambda <- matrix(0, n, 3)

y <- matrix(0, n, 3)
log_dens_ref <- rep(0, n)

tau <- .1  ## offset on diagonal
l <- 0     ## count occasions with invertible cv 
dens1 <- NULL
for ( ii in seq(n) ) {
    mu <- rep(0, 3)
  
    val <- diag( c(val1[ii], val2[ii], val3[ii]) ) + diag(tau, 3)
    ## compute covariance matrix from rotation matrix and eigenvalues
    cv <- R \%*\% val \%*\% t(R)
  
    ## compute parameters for parameter list
    sig[ii,] <- sqrt(diag(cv))
    p12[ii] <- cv[1,2]
    p13[ii] <- cv[1,3]
    p23[ii] <- cv[2,3]
  
    ## compute paramters for Cholesky family
    chol_cv <- solve(chol(cv))  # lambdas come from L^-1 not L
    lamdiag[ii,] <- diag(chol_cv)
    lambda[ii,] <- chol_cv[upper.tri(chol_cv)]
  
    ## Check if cv is invertible 
    if ( !is.matrix(try(chol(cv))) ) l <- l + 1
  
    y[ii,] <- mvtnorm::rmvnorm(1, mu, cv)
  
    log_dens_ref[ii] <- mvtnorm::dmvnorm(y[ii,], mu, cv, log = TRUE)
}
print(l)

## Data
d <- as.data.frame(y)
names(d) <- paste0("y", 1:3)
d$x <- x

## make parameter list for mvn chol family
par <- list()
par[["mu1"]] <- rep(0,n)
par[["mu2"]] <- rep(0,n)
par[["mu3"]] <- rep(0,n)
par[["lamdiag1"]] <- lamdiag[,1]
par[["lamdiag2"]] <- lamdiag[,2]
par[["lamdiag3"]] <- lamdiag[,3]
par[["lambda12"]] <- lambda[,1]
par[["lambda13"]] <- lambda[,2]
par[["lambda23"]] <- lambda[,3]

simdata <- list(
    d   = d,
    par = par,
    y   = y
)

## save(simdata, file = "simdata.rda")
## End of simulation
}
}
\keyword{datasets}
