# Function to obtain the steepness measure based on dyadic dominance indices #

getStp <- function(X){
NormDS <- getNormDS(X);
descend <- array(nrow(NormDS):1);
sortNDS <- sort(NormDS);
NDS <- array(dim=c(nrow(NormDS),1),0.);
for (i in 1:nrow(NormDS))
{point <- descend[i];
NDS[i] <- sortNDS[point]}
NormDS <- NDS;
rnk <- array(dim=c(nrow(NormDS),1),0.);
for (i in 1:nrow(NormDS))
{rnk[i] <- i}
crossp <- sum(NormDS*rnk);
crosss <- (sum(NormDS))*(sum(rnk));
sumsq <- sum(rnk*rnk);
sqsum <- sum(rnk)**2;
Slp <- (nrow(X)*crossp- crosss)/(nrow(X)*sumsq-sqsum);
Stp <- abs(Slp);
return(Stp)
}
