# Function to obtain the weighted sum of individuals' dyadic dominance indices by columns #

getl2 <- function(X){
Dij <- getDij(X);
l1 <- getl1(X);
l2mat <- array(dim=c(nrow(Dij),ncol(Dij)),0.);
for (i in 1:nrow(X))
for (j in 1:ncol(X))
l2mat[i,j] <- Dij[i,j]*l1[i];
l2 <- colSums(l2mat);
return(l2)
}
