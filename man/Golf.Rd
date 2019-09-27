\name{Golf}
\alias{Golf}

\title{Prices of Used Cars Data}

\description{
  This dataset is taken from the Regression Book and is about prices of used VW Golf cars.
}

\usage{data("Golf")}

\format{
The \code{Golf} data contains the following variables:
\describe{
  \item{price:}{Numeric, sale price in 1000 Euro.}
  \item{age:}{Numeric, age of the car in month.}
  \item{kilometer:}{Numeric, kilometer reading in 1000 kilometers.}
  \item{tia:}{Numeric, month until the next TIA appointment (German TUEV).}
  \item{abs:}{Factor, does the car have abs?}
  \item{sunroof:}{Factor, does the car have a sunroof?}
}
}

\references{
  Fahrmeir, L., Kneib, T., Lang, S. and Marx, B. (2013).
  Regression - Models, Methods and Applications, Springer.
  \url{https://www.uni-goettingen.de/de/books/550514.html}.
}

\examples{
data("Golf")
plot(price ~ age, data = Golf)
plot(price ~ kilometer, data = Golf)
}

\keyword{datasets}

