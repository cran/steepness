# Function to obtain the sum of individuals' dyadic dominance indices by rows #

getw1 <- function(X){
Dij <- getDij(X);
w1 <- rowSums(Dij);
return(w1)
}
