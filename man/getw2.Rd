\name{getw2}

\alias{getw2}

\title{Individuals' weighted sum of Dij -w2-}

\description{
Function to obtain individuals' weighted sums of dyadic dominance indices by means of the observed sociomatrix.}

\usage{
   getw2(X)
}

\arguments{
  \item{X}{Original sociomatrix.}
}

\details{
\code{getw2} is necessary for the computation of David's scores.
}

\value{

\item{w2}{Individuals' weighted sum of dyadic dominance indices corrected for chance.}
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

