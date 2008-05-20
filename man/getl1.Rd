\name{getl1}

\alias{getl1}

\title{Individuals' sum of Dji -l1-}

\description{
Function to obtain individuals' sums of dyadic dominance indices by rows by means of the observed sociomatrix.}

\usage{
   getl1(X)
}

\arguments{
  \item{X}{Original sociomatrix.}
}

\details{
\code{getl1} is necessary for the computation of David's scores.
}

\value{

\item{l1}{Individuals' sum of dyadic dominance indices by rows.}
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

