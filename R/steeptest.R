# Function to estimate p values for the steepness index based on the Dij measures #

steeptest <- function(X,rep,names=NULL,option.console=FALSE){

# Is matrix X square? #

  if (nrow(X) != ncol(X))
    return("Error: Matrix X is not square and can not be analyzed")
   
# Limit the number of replications #

  if ((rep < 1) | (rep > 1000000))
    return("Error: Number of replications must be between 1 and 1000000")

# Compute steepness measures for the original matrix #

Dij <- getDij(X);
DS <- getDS(X);
NormDS <- getNormDS(X);
interc <- getinterc(X);
Stp <- getStp(X);

# Carrying out the statistical test by means a C program #

vecX <- array(dim=c(nrow(X)*(ncol(X)-1)/2,1));
m <- 0.;
for (i in 1:nrow(X))
for (j in 1:ncol(X))
{
m <- m+1;
vecX[m] <- X[i,j];}

out <- .C("steep",as.double(vecX),
          as.integer(nrow(X)),
          as.integer(rep),
          res1=double(rep),
          PACKAGE="steepness")

Stpsim <- out$res1;

# Computation of some summary statistics and print results #

if (is.null(names)) names <- paste('Ind.',1:nrow(X))
  else names <- names

dimnames(Dij) <- list(c(names),c(names))
dimnames(DS) <- list(c(names),"DS Values")
dimnames(NormDS) <- list(c(names),"NormDS Values")

data <- array(dim=c(rep,1))
data[,1] <- Stpsim
colnames(data) <- c("Stpsim")
Stp_prightvalue <- (sum(Stp <= data[,"Stpsim"])+1)/(rep+1)
Stp_pleftvalue <- (sum(Stp >= data[,"Stpsim"])+1)/(rep+1)
results <- array((c(Stp, Stp_prightvalue,Stp_pleftvalue,rep,mean(data[,"Stpsim"]),var(data[,"Stpsim"]),
min(data[,"Stpsim"]),quantile(data[,"Stpsim"],.25,names=F),quantile(data[,"Stpsim"],.50,names=F),
quantile(data[,"Stpsim"],.75,names=F),max(data[,"Stpsim"]))),dim=c(11,1))
dimnames(results) <- list(c("Original value", "p-right value", "p-left value", "N simulations", "Mean",
"Variance","Minimum", "25th Pctl","50th Pctl", "75th Pctl","Maximum"),"Stp")

if (option.console == TRUE){
  cat("    ","\n")
  cat("    ","\n")
  cat("    ","\n")
  cat("RESULTS OF STEEPNESS ANALYSIS OF THE MATRIX OF DYADIC DOMINANCES CORRECTED FOR CHANCE","\n")
  cat("=====================================================================================","\n")
  cat("    ","\n")
  cat("    ","\n")
  cat("Dij","\n")
  cat("===","\n")
  cat("    ","\n")
  print(Dij)
  cat("    ","\n")
  cat("    ","\n")
  cat("DAVID'S SCORES","\n")
  cat("==============","\n")
  cat("    ","\n")
  print(DS)
  cat("    ","\n")
  cat("    ","\n")
  cat("NORMALIZED DAVID'S SCORES","\n")
  cat("=========================","\n")
  cat("    ","\n")
  print(NormDS)
  cat("    ","\n")
  cat("    ","\n")
  cat("STEEPNESS","\n")
  cat("=========","\n")
  cat("    ","\n")
  cat("Slope (absolute) = ",Stp,"\n")
  cat("    ","\n")
  cat("    ","\n")
  cat("INTERCEPT","\n")
  cat("=========","\n")
  cat("    ","\n")
  cat("Intercept = ",interc,"\n")
  cat("    ","\n")
  cat("    ","\n")
  cat("    ","\n")
  cat("SUMMARY STATISTICS OF THE RANDOMIZATION PROCEDURE FOR TESTING STEEPNESS BASED ON THE Dij MEASURES","\n")
  cat("=================================================================================================","\n")
  print(results)
}
  else
    list(dyadic.dominance=Dij,david.scores=DS,norm.david.scores=NormDS,steepness=Stp,steep.right.pvalue=Stp_prightvalue,
         steep.left.pvalue=Stp_pleftvalue,intercept=interc,results=results)
}
