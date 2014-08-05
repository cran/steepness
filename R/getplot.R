# Function to plot NormDS against Rank order and fit linear regression line #

getplot <- function (X, names=NULL, method=c("Dij","Pij")) {

if (nrow(X) != ncol(X)) 
  return("Error: Sociomatrix must be square");
if ( is.na(X) || !is.numeric(X))
  return("Error: Sociomatrix must be numeric");

if ((is.null(names)) && (length(names)<=15)) names <- paste('Ind.',1:nrow(X))
if ((is.null(names)) && (length(names)>15)) names <- 1:nrow(X)
if ((nchar(names)<10) && (length(names)<=15)) names <- names
if ((nchar(names)>10) || (length(names)>15)) names <- 1:nrow(X)
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
names <- names[SortNormDS$ix]
NormDS <- array(SortNormDS$x,dim=nrow(X),dimnames=list(names))
rnk <- array(1:nrow(X))

plot(rnk,NormDS, type = "b", lty = 2, ylim = c(0,(nrow(X)-1)), xlab="Individuals in Rank Order", 
     ylab="Normalized David's Scores", main=paste("NormDS (based on",  if (method == 'Dij')
     "Dij)" else "Pij)", "plotted against rank order"), pch=15, col="blue",axes=FALSE,xaxs = "r",yaxs = "i")
axis(1, at=1:length(names),names,tcl = 0.3, adj=c(0.5))
axis(2, at=0:(nrow(X)-1),tcl = 0.3, srt = 1)
box()

myline.fit <- lm(NormDS ~ rnk)

abline(myline.fit,col="green",lwd=2)
lgnd <- paste("Fitted line:  Y = ",format(myline.fit$coefficients[2],digits=4),
"X", "+" ,format(myline.fit$coefficients[1],digits=4))
legend("topright",as.graphicsAnnot(lgnd),lty=1,col="green",lwd=2,box.lty=0)
}
