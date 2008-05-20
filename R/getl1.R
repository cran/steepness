# Function to obtain the sum of individuals' dyadic dominance indices by columns #

getl1 <- function(X){
Dij <- getDij(X);
l1 <-colSums(Dij);
return(l1)
}
