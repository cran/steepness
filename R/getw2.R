# Function to obtain the weighted sum of individuals' dyadic dominance indices by rows #

getw2 <- function(X){
Dij <- getDij(X);
w2mat <- array(dim=c(nrow(Dij),ncol(Dij)),0.);
w1 <- getw1(X);
for (i in 1:nrow(X))
for (j in 1:ncol(X))
w2mat[i,j] <- Dij[i,j]*w1[j];
w2 <- rowSums(w2mat);
return(w2)
}
