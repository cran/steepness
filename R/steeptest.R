# Function to estimate p values for the steepness index based on the Dij or Pij measures #

steeptest <- function(X,rep,names=NULL,method=c("Dij","Pij"),order=TRUE,option.console=FALSE){

# Is matrix X square? #

  if (nrow(X) != ncol(X))
    return("Error: Matrix X is not square and can not be analyzed")
   
# Limit the number of replications #

  if ((rep < 1) | (rep > 1000000))
    return("Error: Number of replications must be between 1 and 1000000")

# Compute steepness measures for the original matrix #

  method <- match.arg(method)
  dyadc <- X + t(X);
  if (method == "Dij"){
    Dij <- X/dyadc-(((X/dyadc)-0.5)/(dyadc+1))
    Dij[is.nan(Dij)] <- 0.
    w1 <- rowSums(Dij);
    w2 <- Dij%*%w1;
    l1 <-colSums(Dij);
    l2 <- t(l1)%*%Dij;
  }
  if (method == "Pij"){
    Pij <- array(dim=c(nrow(X),ncol(X)),0.);
    Pij <- X/dyadc;
    Pij[is.nan(Pij)] <- 0.
    w1 <- rowSums(Pij);
    w2 <- Pij%*%w1;
    l1 <-colSums(Pij);
    l2 <- t(l1)%*%Pij;
  }
  DS <- w1 + w2 - l1 - t(l2);
  maxDS <- nrow(X)*(nrow(X)-1)/2;
  NormDS <- (DS + maxDS)/nrow(X);
  SortNormDS <- sort(NormDS,decreasing=TRUE,index.return=TRUE)
  rnk <- 1:nrow(X)
  Stp <- abs(lm(SortNormDS$x ~ rnk)$coefficients[2])
  names(Stp)<-NULL
  interc <- lm(SortNormDS$x ~ rnk)$coefficients[1]
  names(interc)<-NULL

# Carrying out the statistical test by means a C program #

  vecX <- c(t(X))

  if (method == "Dij")
  out <- .C("steep",as.double(vecX),
          as.integer(nrow(X)),
          as.integer(rep),
          res1=double(rep),
          PACKAGE="steepness")
          
  if (method == "Pij")
  out <- .C("steep2",as.double(vecX),
          as.integer(nrow(X)),
          as.integer(rep),
          res1=double(rep),
          PACKAGE="steepness")

  Stpsim <- out$res1;

# Computation of some summary statistics and print results #

  if (is.null(names)) names <- paste('Ind.',1:nrow(X))
  else names <- names
  
  if (order == TRUE){
    names <- names[SortNormDS$ix]
    DS <- array(DS[SortNormDS$ix],dim=c(nrow(X),1))
    NormDS <- array(SortNormDS$x,dim=c(nrow(X),1))
  }

  if (method == "Dij")
    dimnames(Dij) <- list(c(names),c(names))
  if (method == "Pij")
    dimnames(Pij) <- list(c(names),c(names))
  dimnames(DS) <- list(c(names),"DS Values")
  dimnames(NormDS) <- list(c(names),"NormDS Values")

  data <- array(dim=c(rep,1))
  data[,1] <- Stpsim
  colnames(data) <- c("Stpsim")
  Stp_rightpvalue <- (sum(Stp <= data[,"Stpsim"])+1)/(rep+1)
  Stp_leftpvalue <- (sum(Stp >= data[,"Stpsim"])+1)/(rep+1)
  results <- array((c(Stp, Stp_rightpvalue,Stp_leftpvalue,rep,mean(data[,"Stpsim"]),var(data[,"Stpsim"]),
  min(data[,"Stpsim"]),quantile(data[,"Stpsim"],.25,names=F),quantile(data[,"Stpsim"],.50,names=F),
  quantile(data[,"Stpsim"],.75,names=F),max(data[,"Stpsim"]))),dim=c(11,1))
  dimnames(results) <- list(c("Empirical value", "Right p-value", "Left p-value", "N simulations", "Mean",
  "Variance","Minimum", "25th Pctl","50th Pctl", "75th Pctl","Maximum"),"Stp")
  options(digits=6,scipen=999)
  results <-round(as.data.frame(results),round(log(results[4,],10)))
  if (option.console == TRUE){
    cat("    ","\n")
    cat("    ","\n")
    cat("    ","\n")
    if (method == "Dij"){
    cat("RESULTS OF STEEPNESS ANALYSIS OF THE MATRIX OF DYADIC DOMINANCES CORRECTED FOR CHANCE","\n")
    cat("=====================================================================================","\n")}
    if (method == "Pij"){
    cat("RESULTS OF STEEPNESS ANALYSIS OF THE MATRIX OF WIN PROPORTIONS","\n")    
    cat("==============================================================","\n")}
    cat("    ","\n")
    cat("    ","\n")
    if (method == "Dij")
    cat("Dij","\n")
    if (method == "Pij")
    cat("Pij","\n")    
    cat("===","\n")
    cat("    ","\n")
    if (method == "Dij")
    print(Dij)
    if (method == "Pij")
    print(Pij)
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
    if (method == "Dij")
    cat("SUMMARY STATISTICS OF THE RANDOMIZATION PROCEDURE FOR TESTING STEEPNESS BASED ON THE Dij MEASURES","\n")
    if (method == "Pij")
    cat("SUMMARY STATISTICS OF THE RANDOMIZATION PROCEDURE FOR TESTING STEEPNESS BASED ON THE Pij MEASURES","\n")
    cat("=================================================================================================","\n")
    print(results)
  }
  else
  {
    if (method == "Dij") matdom <- Dij
    else matdom <- Pij
    list(dyadic.dominance=matdom,david.scores=DS,norm.david.scores=NormDS,steepness=Stp,steep.right.pvalue=Stp_rightpvalue,
         steep.left.pvalue=Stp_leftpvalue,intercept=interc,results=results)
  }
}
