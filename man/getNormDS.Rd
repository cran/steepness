\name{getNormDS}

\alias{getNormDS}

\title{Normalized David's scores -NormDS-}

\description{
Function to obtain normalized David's scores from the observed sociomatrix.}

\usage{
   getNormDS(X, names=NULL, method=c("Dij","Pij"))
}

\arguments{
  \item{X}{Empirical sociomatrix containing wins-losses frequencies in dyadic encounters. The matrix must be square and numeric.}
  \item{names}{Character vector with the names of individuals. This vector is NULL by default}
  \item{method}{A character string indicating which dyadic dominance measure is to be used for the computation of David's scores. One of "Dij" or "Pij", can be abbreviated.}
}

\details{
\code{getNormDS} is obtained by means of the following expression: \eqn{NormDS = (DS + N(N-1)/2)/N}
}

\value{

\item{NormDS}{Normalized David's scores based on dyadic dominance indices corrected for chance or based on the win proportions, depending on the method specified.}
}

\references{
David, H. A. (1988). \emph{The Method of Paired Comparisons}. London: C. Griffin.

de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}

\seealso{ \code{\link{getDij}}, \code{\link{getPij}}, \code{\link{getDS}}.
}

\examples{

##############################################################################
###               Example taken from Vervaecke et al. (2007):              ###
##############################################################################

X <- matrix(c(0,58,50,61,32,37,29,39,25,8,0,22,22,9,27,20,10,48,
              3,3,0,19,29,12,13,19,8,5,8,9,0,33,38,35,32,57,
              4,7,9,1,0,28,26,16,23,4,3,0,0,6,0,7,6,12,
              2,0,4,1,4,4,0,5,3,0,2,1,1,5,8,3,0,10,3,1,3,0,0,4,1,2,0),
              nrow=9,byrow=TRUE)

individuals <- c("V","VS","B","FJ","PR","VB","TOR","MU","ZV")

res <- getNormDS(X,names=individuals,method="Dij")

print(res,digits=3)

}

\keyword{misc}

\keyword{htest}

