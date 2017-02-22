\name{neighbormatrix}
\alias{neighbormatrix}
\alias{plotneighbors}

\title{
  Compute a Neighborhood Matrix from Spatial Polygons
}

\description{
  The function takes a \code{\link[sp]{SpatialPolygonsDataFrame}} and computes
  the neighbor penalty matrix that can be used to fit a Markov
  random field, e.g., using the smooth constructor
  \code{\link[mgcv]{smooth.construct.mrf.smooth.spec}}.
}

\usage{
## Compute the neighborhood matrix.
neighbormatrix(x, type = c("boundary", "dist", "delaunay", "knear"),
  k = 1, id = NULL, nb = FALSE, names = NULL, ...)

## Plot neighborhood structure.
plotneighbors(x, add = FALSE, ...)
}

\arguments{
  \item{x}{An object of class \code{"SpatialPolygons"} or \code{"SpatialPolygonsDataFrame"}.}
  \item{type}{Which type of neighborhood structure should be used,
    \code{"boundary"} uses \code{\link[spdep]{poly2nb}}, \code{"dist"}
    uses function \code{\link[spdep]{dnearneigh}}, \code{"delaunay"}
    uses function \code{\link[spdep]{tri2nb}} and option
    \code{"knear"} applies function \code{\link[spdep]{knn2nb}}.}
  \item{k}{For \code{type = "knear"}, specifies number of nearest neighbors.}
  \item{id}{An identifier variable for which the penalty matrix should be computed.}
  \item{nb}{Should only the neighborhood object be returned.}
  \item{names}{Specifies the column where the regions names are provided in the data
    slot in the \code{"SpatialPolygonsDataFrame"}}
  \item{add}{Should the neighborhood structure be added to an existing plot?}
  \item{\dots}{Arguments to be passed to function \code{\link[spdep]{poly2nb}},
    \code{\link[spdep]{dnearneigh}}, \code{\link[spdep]{tri2nb}} or
    \code{\link[spdep]{knn2nb}}.}
}

\seealso{
  \code{\link[mgcv]{smooth.construct.mrf.smooth.spec}}, \code{\link[spdep]{dnearneigh}},
  \code{\link[spdep]{tri2nb}}, \code{\link[spdep]{knn2nb}}.
}

\examples{
data("LondonFire")

## Compute polygon boundary based
## neighborhood matrix.
nm <- neighbormatrix(LondonBoroughs)
print(nm)

## Plot neighborhood structures.
plotneighbors(LondonBoroughs)
plotneighbors(LondonBoroughs, type = "delaunay")
plotneighbors(LondonBoroughs, type = "dist", d1 = 0, d2 = 0.15)
}

\keyword{models}
\keyword{regression}
\keyword{aplot}

