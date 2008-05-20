# Function to obtain the matrix of dyadic dominance indices -Dij- #

getDij <- function(X){
dyadc <- X + t(X);
Dij <- array(dim=c(nrow(X),ncol(X)),0.);
for (i in 1:nrow(X))
for (j in 1:ncol(X))
if (i!=j){
if (dyadc[i,j]!=0)
{Dij[i,j] <- (X[i,j]/dyadc[i,j])-(((X[i,j]/dyadc[i,j])-0.5)/(dyadc[i,j]+1))}
else
Dij[i,j] <- 0.}
return(Dij)
}
