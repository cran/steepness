# Function to plot NormDS against Rank order and fit linear regression line #

getplot <- function (X, names=NULL) {

if ((is.null(names)) & (length(names)<=15)) names <- paste('Ind.',1:nrow(X))
if ((is.null(names)) & (length(names)>15)) names <- 1:nrow(X)
if ((nchar(names)<10) & (length(names)<=15)) names <- names
if ((nchar(names)>10) | (length(names)>15)) names <- 1:nrow(X)

matord <- array(dim=c(nrow(X),ncol(X)),0)
NormDS <- getNormDS(X);
namestemp <- names
descend <- array(nrow(X):1);
sortNDS <- sort(NormDS);
for (i in 1:length(names))
{namestemp[rank(NormDS)[i]] <- names[i]}
names <- namestemp
for (i in 1:(ncol(X)-1))
for (j in (i+1):ncol(X))
{matord[rank(NormDS)[i],rank(NormDS)[j]]=X[i,j]
matord[rank(NormDS)[j],rank(NormDS)[i]]=X[j,i]}
X <- matord 
NDS <- array(dim=c(nrow(NormDS),1),0.);
for (i in 1:nrow(NormDS))
{point <- descend[i]
NDS[i] <- sortNDS[point]
namestemp[i] <- names[point]}
NormDS <- NDS;
names <- namestemp
rownames(NormDS) <- names
for (i in 1:(ncol(X)-1))
for (j in (i+1):ncol(X))
{matord[i,j] <- X[descend[i],descend[j]]
matord[j,i] <- X[descend[j],descend[i]]}
X <- matord
rownames(X) <- names
colnames(X) <- names
rnk <- array(1:nrow(X));

plot(rnk,NormDS, type = "b", lty = 2, ylim = c(0,(nrow(X)-1)), xlab="Individuals in Rank Order", 
     ylab="Normalized David's Scores", main="NormDS plotted against rank order", 
     pch=15, col="blue",axes=FALSE,xaxs = "r",yaxs = "i")
axis(1, at=1:length(names),names,tcl = 0.3, adj=c(0.5))
axis(2, at=0:(nrow(X)-1),tcl = 0.3, srt = 1)
box()

myline.fit <- lm(NormDS ~ rnk)

abline(myline.fit,col="green",lwd=2) 
 
lgnd <- paste("Fitted line:  Y = ",format(myline.fit$coefficients[2],digits=4),"X"," + ",format(myline.fit$coefficients[1],digits=4))
legend("topright",as.graphicsAnnot(lgnd),lty=1,col="green",lwd=2,box.lty=0)
}
