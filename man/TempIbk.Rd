\docType{data}
\name{TempIbk}
\alias{TempIbk}
\title{Temperature data.}
\format{
An object of class \code{data.frame} with 1798 rows and 17 columns.
}
\usage{
data("TempIbk")
}
\description{
Temperature Data for Innsbruck Airport
}
\details{
Numerical weather predictions (NWP) and observations of
2 meter temperature at Innsbruck Airport.
The observations from the SYNOP station 11120 cover 5 years from
2015-01-01 to 2019-31-12.
The NWP data are derived from GEFS reforecasts (Hamill et al. 2013).
The data contain following variables:

\itemize{
  \item \code{init}: Time of initialization of the NWP model.
  \item \code{obs_*}: Observations for lead time \code{*}.
  \item \code{mean_ens_*}: NWP ensemble mean for lead time \code{*}.
  \item \code{logsd_ens_*}: NWP logarithm of ensemble standard deviation for lead time \code{*}.
  \item \code{yday}: Yearday.
}
}

\seealso{
  \code{\link{mvnchol_bamlss}}
}

\examples{
\dontrun{## Innsbruck temperature data.
data("TempIbk", package = "bamlss")

## Five lead times.
lead <- seq(192, 216, by = 6)

## Set up formulas.
f <- c(
  ## mu equations
  sprintf('obs_\%s ~ s(yday, bs = "cc") + s(yday, bs = "cc", by = mean_ens_\%s)', lead, lead),

  ## lambda diag equations
  sprintf('lamdiag\%s ~ s(yday, bs = "cc") + s(yday, bs = "cc", by = logsd_ens_\%s)', 1:5, lead),

  ## lambda off-diag equations
  sprintf('lambda\%s ~ s(yday, bs = "cc")', apply(combn(1:5, 2), 2, paste, collapse = ""))
)
f <- lapply(f, as.formula)

## Multivariate normal family with basic Cholesky parameterization.
fam <- mvnchol_bamlss(k = 5, type = "basic")

## Fit model.
set.seed(123)
b <- bamlss(f, family = fam, data = TempIbk, optimizer = opt_boost, maxit = 1000)

## Show estimated effects.
par(mfrow = c(2, 2))
plot(b, model = "mu1", scale = 0, spar = FALSE)
plot(b, model = "lamdiag2", term = "s(yday)", spar = FALSE)
plot(b, model = "lambda12")

## Predict sample case.
nd <- subset(TempIbk, format(init, "\%Y-\%m-\%d") \%in\% c("2015-01-03", "2015-10-10"))
fit <- predict(b, newdata = nd, type = "parameter")

## Plot correlation matrix for GEFS initialization 2015-10-10.
plot_cor <- function(i) {
    image(lead, lead, fam$correlation(fit)[[i]][5:1, ], zlim = c(0, 1),
    	 col = hcl.colors(10, "Blues 3", rev = TRUE), axes = FALSE,
    	 xlab = "lead time in hours", ylab = "lead time in hours",
    	 main = sprintf("Correlation matrix fitted for \%s", nd[i, "init"]))
    axis(1, lead)
    axis(2, lead, rev(lead))
    box()
}
par(mfrow = c(1, 2))
plot_cor(1)
plot_cor(2)

## Plot means and standard deviations.
plot_ms <- function(i) {
	stdev <- fam$stdev(fit)[[i]]
	means <- fam$means(fit)[[i]]
	lower <- means - stdev
	upper <- means + stdev
	
	plot(lead, means, type = 'b', cex = 2, lwd = 1, lty = 2, axes = FALSE,
		 ylim = c(-6, 16), # c(min(lower), max(upper)),
		 ylab = expression("Temperature in " * degree * "C"),
		 xlab = "lead time in hours",
		 main = sprintf("Means +/- one st. dev. for \%s", nd[i, "init"]))
	segments(lead, y0 = lower, y1 = upper)
	axis(1, lead)
	axis(2)
	box()
}
par(mfrow = c(1, 2))
plot_ms(1)
plot_ms(2)
}

}
\references{
Hamill TM, Bates GT, Whitaker JS, Murray DR, Fiorino M, Galarneau Jr TJ,
    Zhu Y, Lapenta W (2013). NOAA's Second-Generation Global Medium-Range
    Ensemble Reforecast Data Set. \emph{Bulletin of the American Meteorological
    Society}, 94(10), 1553-1565.
}
\keyword{datasets}
