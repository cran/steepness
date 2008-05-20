# Function to obtain the intercept based on dyadic dominance indices */

getinterc <- function(X){
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
interc <- (sum(NormDS) -(Slp*sum(rnk)))/nrow(X);
return(interc)
}
