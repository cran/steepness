\name{getNormDS}

\alias{getNormDS}

\title{Normalized David's scores based on Dij -NormDS-}

\description{
Function to obtain normalized David's scores based on dyadic dominance indices by means of the observed sociomatrix.}

\usage{
   getNormDS(X)
}

\arguments{
  \item{X}{Original sociomatrix.}
}

\details{
\code{getNormDS} is obtained by means of the following expression: \eqn{NormDS = (DS + N(N-1)/2)/N}
}

\value{

\item{NormDS}{Normalized David's scores based on dyadic dominance indices corrected for chance.}
}

\references{

de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}

\seealso{ \code{\link{getDij}}, \code{\link{getDS}}.
}

\examples{

  X <- matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,byrow=TRUE)
}

\keyword{misc}

\keyword{htest}

