\name{steeptest}

\alias{steeptest}

\title{Statistical significance for steepness of dominance hierarchies statistic}

\description{
Estimates statistical significance for steepness measure on the basis of dyadic dominance indices corrected for chance \emph{Dij} or based on proportions of wins \emph{Pij}.
}


\usage{
   steeptest(X, rep, names=NULL, method=c("Dij","Pij"), order=TRUE)
}

\arguments{
  \item{X}{Empirical sociomatrix containing wins-losses frequencies in dyadic encounters. The matrix must be square and numeric.}
  \item{rep}{Number of simulations for carrying out the randomization test.}
  \item{names}{Character vector with individuals' names.}
  \item{method}{A character string indicating which dyadic dominance measure is to be used for the computation of David's scores. One of "Dij" or "Pij", can be abbreviated.}
  \item{order}{Logical, if TRUE, results for Dij, DS and NormDS are ordered according to the individuals' NormDS values. TRUE by default.}
}

\details{
\code{steeptest} estimates statistical significance for steepness measures based on dyadic dominance index corrected for chance \emph{Dij} or based on the matrix of win proportions \emph{Pij}, depending on the \code{method} specified. This procedure simulates a number of sociomatrices under a uniform distribution by means of callings to C routine \emph{steep}, then computes steepness based on \emph{Dij} or \emph{Pij}. Specifically, it computes normalized David's scores, see \code{\link{getNormDS}} for more details. Then it computes the steepness measure based on these indices, see \code{\link{getStp}}. After \code{rep} simulations the sampling distribution for the statistic (\emph{Stp}) is estimated. Then statistical significance is computed as follows when results are shown by means of \code{summary} method:
\eqn{p=NS+1/NOS+1}
Where \emph{NS} is computed as:
\enumerate{
\item The number of times that simulated values are greater than or equal to the empirical value, if right-tailed p value is calculated.
\item And the number of times that simulated values are lower than or equal to the empirical value, if left-tailed p value is calculated.
}
And \emph{NOS} represents the number of simulated values.
}

\value{
\code{steeptest} returns an object of class steeptest containing the following components:
  \item{call }{Function call.}
\item{names}{Character vector with individuals' names.}
\item{method}{A character string indicating which dyadic dominance measure is used for the computation of David's scores.}
\item{rep}{Number of simulations for carrying out the randomization test.}
\item{matdom}{If \code{method} is set to be \emph{Dij} the function returns the matrix of observed dyadic dominance indices corrected for chance. If \code{method} is \emph{Pij} the matrix of proportions of wins is returned as a part of the output.}
\item{DS}{David's scores based on \emph{Dij} or \emph{Pij}, depending on the specification of the \code{method}.}
\item{NormDS}{Normalized David's scores based on dyadic dominance indices corrected for chance or on proportions of wins in dyadic encounters.}
\item{Stp}{Steepness value based on Normalized David's scores.}
\item{interc}{Intercept of the fitted line based on Normalized David's scores.}
\item{Stpsim}{The function provides results of the randomization procedure for the steepness measure based on \emph{NormDS}.}
}

\references{
David, H. A. (1988). \emph{The Method of Paired Comparisons}. London: C. Griffin.

de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}

\seealso{ \code{\link{getDij}}, \code{\link{getPij}}, \code{\link{getNormDS}}
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

 STP <- steeptest(X, rep=9999, names=individuals, method="Dij", order=TRUE)
 summary(STP)
 plot(STP)
}

\keyword{misc}

