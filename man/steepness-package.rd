\name{steepness-package}
\alias{steepness}
\docType{package}
\title{
Testing Steepness of Dominance Hierarchies
}
\description{
Steepness is a package that computes steepness as a property of dominance hierarchies. Steepness is defined as the absolute slope of the straight line fitted to the normalized David's scores. The normalized David's scores are obtained on the basis of dyadic dominance indices corrected for chance. Given an observed sociomatrix, it computes hierarchy's steepness and estimates statistical significance by means of a randomization test (see de Vries, Stevens and Vervaecke, 2006).
}
\details{
\tabular{ll}{
Package: \tab steepness\cr
Version: \tab 0.1\cr
Date: \tab 2008-20-05\cr
Depends: \tab >= 2.6.1\cr
License: \tab GPL version 2 or newer\cr
}

Index:
\preformatted{
getDij                  Dyadic dominance index corrected for chance -Dij-
getDs                   David's scores based on Dij -DS-
getinterc               Intercept of the fitted line based on Dij -interc-
getl1                   Individuals' sum of Dji -l1-
getl2                   Individuals' weighted sum of Dji -l2-
getNormDS               Normalized David's scores based on Dij -NormDS-
getStp                  Hierarchy's steepness based on Dij -Stp-
getw1                   Individuals' sum of Dij -w1-
getw2                   Individuals' weighted sum of Dij -w2-
getplot                 Steepness plot
steeptest               Statistical significance for steepness statistic
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
For more information see \code{RcmdrPlugin.steepness}.
}
