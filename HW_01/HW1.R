# Author: Conner Engle
# Johns Hopkins: AS.410.671.82.SU22

# Analysis of GEO: GDS2880

#######################
#          1          #
####################### 

# Download data and verify 22 arrays of 22,283 probesets

setwd("\\Users\\Conner\\Desktop\\JHUDataSets\\")
Dataset <- read.table("/Users/Conner/Desktop/renal_cell_carcinoma.txt", header=TRUE, row.names=1, strip.white=TRUE)

dim(Dataset)


#######################
#          2          #
#######################

# Label header columns to add tumor vs. normal identity

identity <- rep(c('Normal','Tumor'),each=11)
colnames(Dataset) <- paste(colnames(Dataset),identity,sep=”_”)

colnames(Dataset)


#######################
#          3          #
#######################

# Correlation Plot

data.cor <- cor(Dataset, use='pairwise.complete.obs')

layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))
cx <- rev(colorpanel(25,"yellow","black","blue"))
leg <- seq(min(data.cor,na.rm=T),max(data.cor,na.rm=T),length=10)
image(data.cor,main="Correlation plot of normal vs renal cell carcinoma tissue",axes=F,col=cx)
axis(1,at=seq(0,1,length=ncol(data.cor)),label=dimnames(data.cor)[[2]],cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(data.cor)),label=dimnames(data.cor)[[2]],cex.axis=0.9,las=2)

image(as.matrix(leg),col=cx,axes=F)
tmp <- round(leg,2)
axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=1)


# Hierarchical Clustering Dendrogram

dat <- t(Dataset)
dat.dist <- dist(dat,method='euclidian')
dat.clust <- hclust(dat.dist,method='single')
plot(dat.clust,labels=names(dat),cex=0.75)


# CV vs. mean plot

dat.mean <- apply(log2(Dataset),2,mean)
dat.sd <- sqrt(apply(log2(Dataset),2,var))
dat.cv <- dat.sd/dat.mean

plot(dat.mean,dat.cv,main="GDS2880 Normal vs. Renal Cell Carcinoma Dataset\nSample CV vs. Mean",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(dat.mean,dat.cv,bg="lightblue",col=1,pch=21)
text(dat.mean,dat.cv,label=dimnames(Dataset)[[2]],pos=1,cex=0.5)


# Average Correlation Plot

dat.avg <- apply(data.cor,1,mean)
par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of Tumor/Normal samples",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(Dataset)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")


#######################
#          4          #
#######################

# Install and load impute library

source("http://bioconductor.org/biocLite.R") 
biocLite("impute") 

Library(impute)


#######################
#          5          #
#######################

# Remove outlier samples identified in question 3

dat <- subset(Dataset,select=-c(GSM146798_Normal,GSM146799_Tumor))
length(colnames(dat))


#######################
#          6          #
#######################

# Extract KNG1 and AQP2 probesets then plot expression intensity vs. samples

KNG1.probe1 <- dat["206054_at",]
KNG1.probe2 <- dat["217512_at",]
AQP2 <- dat["206672_at",]

plot(c(1, ncol(KNG1.probe1)), range(KNG1.probe1), type = 'n', main="Profile plot of KNG1 206054_at probe", xlab = "", ylab = "Expression Intensity", axes = F)
axis(side=1, at=c(1:20), labels = dimnames(KNG1.probe1)[[2]], cex.axis = 0.4, las = 2)
axis(side=2)
for(i in 1:length(KNG1.probe1)) { 
    dat.kng.1 <- as.numeric(KNG1.probe1[i, ]) 
		lines(c(1:ncol(KNG1.probe1)), dat.kng.1, col = i, lwd = 2) 
}
	
plot(c(1, ncol(KNG1.probe2)), range(KNG1.probe2), type = 'n', main="Profile plot of KNG1 217512_at probe", xlab = "", ylab = "Expression Intensity", axes = F)
axis(side=1, at=c(1:20), labels = dimnames(KNG1.probe2)[[2]], cex.axis = 0.4, las = 2)
axis(side=2)
for(i in 1:length(KNG1.probe2)) { 
    dat.kng.2 <- as.numeric(KNG1.probe2[i, ]) 
		lines(c(1:ncol(KNG1.probe2)), dat.kng.2, col = i, lwd = 2) 
}
	
plot(c(1, ncol(AQP2)), range(AQP2), type = 'n', main="Profile plot of AQP2", xlab = "", ylab = "Expression Intensity", axes = F)
axis(side=1, at=c(1:20), labels = dimnames(AQP2)[[2]], cex.axis = 0.4, las = 2)
axis(side=2)
for(i in 1:length(AQP2)) { 
    dat.aqp2 <- as.numeric(AQP2[i, ]) 
		lines(c(1:ncol(AQP2)), dat.aqp2, col = i, lwd = 2) 
}


#######################
#          7          #
#######################

# Assess accuracy of missing value imputation

dat.kng1 <- as.matrix(dat)
original <- dat.kng1['206054_at','GSM146784_Normal']


#######################
#          8          #
#######################

# Estimate missing values using 6 enarest neighbors and Euclidean distance

dat.knn <- impute.knn(dat.kng1,6)
imputed <- dat.knn$data['206054_at','GSM146784_Normal']
imputed


#######################
#          9          #
#######################

# Calculate relative error for imputed value using its actual value

rel.error <- abs(original-imputed)/original
rel.error


#######################
#          10         #
#######################

# Impute values using SVD imputation method

dat.svd <- as.matrix(dat)
dat.svd['206054_at','GSM146784_Normal'] <- NA
dat.pca <- pca(dat.svd, method=’svdImpute’, nPcs=9)
dat.output <- completeObs(dat.pca)


#######################
#          11         #
#######################

# Plot gene profile plot where two imputed values are represented as different colored points and actual value is the third point

knn <- dat.knn$data['206054_at',]
svd <- dat.output['206054_at',]
actual <- dat['206054_at',]
combined <- rbind(actual, knn, svd)

plot(c(1, ncol(combined)), range(combined), type = "n", main = "Profile plot of actual vs. imputed values of\nGSM146784 normal tissue in probeset 206054_at", xlab = "", ylab = "Expression Intensity", axes = F)
axis(side = 1, at = 1:length(combined), labels = dimnames(combined)[[2]], cex.axis = 0.4, las = 2)
axis(side = 2)
for(i in 1:length(combined)) { 
    dat.y <- as.numeric(combined[i, ]) 
    lines(c(1:ncol(combined)), dat.y, col = i, lwd = 2) 
} 
legend("topright", legend=c("Actual", "knn","SVD"),col=c(1,2,3),lty=1,cex=0.8)




dat.kng1['206054_at','GSM146784_Normal'] <- NA


