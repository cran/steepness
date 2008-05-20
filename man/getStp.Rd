\name{getStp}

\alias{getStp}

\title{Hierarchy's steepness based on Dij -Stp-}

\description{
Function to obtain hierarchy's steepness measure based on dyadic dominance indices by means of the observed sociomatrix.}

\usage{
   getStp(X)
}

\arguments{
  \item{X}{Original sociomatrix.}
}

\details{
\code{getStp} is the absolute value of the slope of the best-fitted line between the normalized David's scores and the rank dominance in a decreasing order. The regression is obtained by Ordinary Least Squares method.
}

\value{

\item{getStp}{Steepness measure based on dyadic dominance indices corrected for chance.}
}

\references{

de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}

\seealso{ \code{\link{getDij}}, \code{\link{getNormDS}}.
}

\examples{

 X <- matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,byrow=TRUE)
 
}

\keyword{misc}

\keyword{htest}

