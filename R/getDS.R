# Function to obtain the David's scores (DS) based on dyadic dominance indices #

getDS <- function(X){
Dij <- getDij(X);
w1 <- getw1(X);
w2 <- getw2(X);
l1 <- getl1(X);
l2 <- getl2(X);
DS <- array(dim=c(nrow(Dij),1),0.);
for (i in 1:nrow(X))
DS[i] <- w1[i] + w2[i] - l1[i] - l2[i];
return(DS)
}
