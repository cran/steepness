\name{getl2}

\alias{getl2}

\title{Individuals' weighted sum of Dji -l2-}

\description{
Function to obtain individuals' weighted sums of dyadic dominance indices by rows by means of the observed sociomatrix.}

\usage{
   getl2(X)
}

\arguments{
  \item{X}{Original sociomatrix.}
}

\details{
\code{getl2} is necessary for the computation of David's scores.
}

\value{

\item{l2}{Individuals' weighted sum of dyadic dominance indices by rows.}
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

