# Function to obtain the normalized David's scores (NormDS) based on dyadic dominance indices #

getNormDS <- function(X){
Dij <- getDij(X);
DS <- getDS(X);	
NormDS <- array(dim=c(nrow(Dij),1),0.);
maxDS <- nrow(X)*(nrow(X)-1)/2;
for (i in 1:nrow(X))
NormDS[i] <- (DS[i] + maxDS)/nrow(X);
return(NormDS)
}
