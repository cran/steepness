\name{getDij}

\alias{getDij}

\title{Dyadic dominance index corrected for chance -Dij-}

\description{
Function to obtain matrix of dyadic dominance indices corrected for chance by means of the observed sociomatrix.
}

\usage{
   getDij(X)
}

\arguments{
  \item{X}{Original sociomatrix.}
}

\details{
\code{getDij} is only applied for square matrices in which the set of \emph{n} actors is also the set of \emph{n} partners.
}

\value{

\item{Dij}{Matrix of observed dyadic dominance indices corrected for chance.}
}

\references{

de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}

\seealso{ \code{\link{steeptest}}.
}

\examples{

  X <- matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,byrow=TRUE)
}

\keyword{misc}

\keyword{htest}


\keyword{misc}

\keyword{htest}
