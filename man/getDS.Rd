\name{getDS}

\alias{getDS}

\title{David's scores based on Dij -DS-}

\description{
Function to obtain David's scores based on dyadic dominance indices by means of the observed sociomatrix.
}

\usage{
   getDS(X)
}

\arguments{
  \item{X}{Original sociomatrix.}
}

\details{
\code{getDS} is obtained by means of the following expression: \eqn{DS = w1 + w2 - l1 - l2} where \emph{w1} is the sum of i's dyadic dominance indices corrected for chance; \emph{w2} is the weighted sum of i`s dyadic dominance indices corrected for chance; \emph{l1} is the sum of j's dyadic dominance indices and \emph{l2} is the sum of j's dyadic dominance indices corrected for chance.
}

\value{

\item{DS}{David's scores based on dyadic dominance indices.}
}

\references{

de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}

\seealso{ \code{\link{getDij}}, \code{\link{getw1}}, \code{\link{getw2}}, \code{\link{getl1}}, \code{\link{getl2}}.
}

\examples{

  X <- matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,byrow=TRUE)
}

\keyword{misc}

\keyword{htest}

