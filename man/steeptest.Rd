\name{steeptest}

\alias{steeptest}

\title{Statistical significance for steepness of dominance hierarchies statistic}

\description{
Estimates statistical significance for steepness measure based on dyadic dominance indices corrected for chance \emph{Dij}.
}


\usage{
   steeptest(X, rep, names=NULL, option.console=FALSE)
}

\arguments{
  \item{X}{Original sociomatrix.}
  \item{rep}{Number of simulations for carrying out the randomization test.}
  \item{names}{Character vector with individuals' names.}
  \item{option.console}{Logical, if TRUE, results of the procedure will be printed in R console. FALSE by default.}
}

\details{
\code{steeptest} estimates statistical significance for steepness measures based on dyadic dominance index corrected for chance \emph{Dij}. This procedure simulates a number of sociomatrices under a uniform distribution by means of callings to C routine \emph{steep}, then computes steepness based on \emph{Dij}. Specifically, it computes normalized David's scores based on \emph{Dij}, see \code{\link{getNormDS}} for more details. Then it computes the steepness measure based on these indices, see \code{\link{getStp}}. After \emph{rep} simulations the sampling distribution for the statistic (\emph{Stp}) is estimated. Then statistical significance is computed as follows:
\eqn{p=NS+1/NOS+1}
Where \emph{NS} is computed as:
\enumerate{
\item The number of times that simulated values are greater than or equal to the empirical value, if right-tailed p value is calculated.
\item And the number of times that simulated values are lower than or equal to the empirical value, if left-tailed p value is calculated.
}
And \emph{NOS} represents the number of simulated values.
}

\value{
\code{steeptest} returns a list containing the following components:

\item{dyadic.dominance}{Matrix of observed dyadic dominance indices corrected for chance.}
\item{david.scores}{David's scores based on dyadic dominance indices corrected for chance.}
\item{norm.david.scores}{Normalized David's scores based on proportions of wins in dyadic encounters.}
\item{steepness}{Steepness value based on Normalized David's scores.}
\item{steep.right.pvalue}{Statistical significance for the empirical steepness value under the null hypothesis of random wins.}
\item{steep.left.pvalue}{Statistical significance for the empirical steepness value under the null hypothesis of random wins.}
\item{intercept}{Intercept of the fitted line based on Normalized David's scores.}
\item{results}{The function provides results of the randomization procedure for the steepness measure based on \emph{NormDS}.}
}

\references{

de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}

\seealso{ \code{\link{getDij}}, \code{\link{getNormDS}}
}

\examples{

  X <- matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,byrow=TRUE)
  rep <- 10000
  names <- c("Ag1","Ag2","Ag3")
}

\keyword{misc}

