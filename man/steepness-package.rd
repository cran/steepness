\name{steepness-package}
\alias{steepness}
\docType{package}
\title{
Testing Steepness of Dominance Hierarchies
}
\description{
Steepness is a package that computes steepness as a property of dominance hierarchies. Steepness is defined as the absolute slope of the straight line fitted to the normalized David's scores. The normalized David's scores can be obtained on the basis of dyadic dominance indices corrected for chance or from the matrix of win proportions. Given an observed sociomatrix, it computes hierarchy's steepness and estimates statistical significance by means of a randomization test (see de Vries, Stevens and Vervaecke, 2006).
}
\details{
\tabular{ll}{
Package: \tab steepness\cr
Version: \tab 0.2-2\cr
Date: \tab 2014-29-09\cr
Depends: \tab >= 3.1.0\cr
License: \tab GPL version 2 or newer\cr
}

Index:
\preformatted{
getDij            Dyadic dominance index corrected for chance -Dij-
getDS             David's scores -DS-
getNormDS         Normalized David's scores -NormDS-
getOrderedMatrix  Ordered matrix according to NormDS values
getPij            Matrix of proportions of wins -Pij-
getStp            Hierarchy's steepness measure -Stp-
getwl             Several win and loss measures at individual level
steeptest         Statistical significance for steepness statistic
}
}
\author{
David Leiva <dleivaur@ub.edu> & Han de Vries <J.deVries1@uu.nl>.

Maintainer: David Leiva <dleivaur@ub.edu>
}
\references{
de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}
\keyword{package}
\seealso{
For more information see: \code{\link{getDij}}, \code{\link{getDS}}, \code{\link{getNormDS}}, \code{\link{getOrderedMatrix}}, \code{\link{getPij}}, \code{\link{getStp}}, \code{\link{getwl}}, \code{\link{steeptest}}.
}
