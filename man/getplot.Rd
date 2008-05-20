\name{getplot}

\alias{getplot}

\title{Steepness plot}

\description{
Function to plot NormDS against Rank order and fit linear regression line
}

\usage{
   getplot(X,names)
}

\arguments{
  \item{X}{Original sociomatrix.}
  \item{names}{Character vector with the names of individuals. This vector is NULL by default}
}

\details{
\code{getplot} plots NormDS against rank scores.
}

\references{

de Vries, H., Stevens, J. M. G., & Vervaecke, H. (2006). Measuring and testing the steepness of dominance hierarchies. \emph{Animal Behaviour}, \emph{71}, 585-592.
}

\author{
David Leiva \email{dleivaur@ub.edu} & Han de Vries \email{J.deVries1@uu.nl}.
}

\seealso{ \code{\link{getNormDS}}, \code{\link{steeptest}}.
}

\examples{

  X <- matrix(c(0, 1, 2, 1, 0, 1, 3, 1, 0),nrow=3,byrow=TRUE)
  names <- c("Ag1","Ag2","Ag3")
}

\keyword{misc}

\keyword{htest}

