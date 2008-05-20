\name{getinterc}

\alias{getinterc}

\title{Intercept of the fitted line based on Dij -interc-}

\description{
Function to obtain intercept of the fitted line based on dyadic dominance indices by means of the observed sociomatrix.}

\usage{
   getinterc(X)
}

\arguments{
  \item{X}{Original sociomatrix.}
}

\details{
\code{getinterc} is the intercept of the best-fitted line between the normalized David's scores and the rank dominance in a decreasing order. The regression is obtained by Ordinary Least Squares method.
}

\value{

\item{getinterc}{Intercept based on dyadic dominance indices corrected for chance.}
}

\references{

de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}

\seealso{ \code{\link{getDij}}, \code{\link{getNormDS}}, \code{\link{getStp}}.
}

\examples{

  X <- matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,byrow=TRUE)
}

\keyword{misc}

\keyword{htest}
