# Author: Conner Engle
# Johns Hopkins: AS.410.671.82.SU22

# Power and Sample Size

# Alizadeh, A. A., Eisen, M. B., Davis, R. E., Ma, C., Lossos, I. S., Rosenwald, A., Boldrick, J. C., Sabet, H., Tran, T., Yu, X., Powell, J. I., Yang, L., 
# Marti, G. E., Moore, T., Hudson, J., Jr, Lu, L., Lewis, D. B., Tibshirani, R., Sherlock, G., Chan, W. C., … Staudt, L. M. 
# (2000). 
# Distinct types of diffuse large B-cell lymphoma identified by gene expression profiling. 
# Nature, 403(6769), 503–511. 
# https://doi.org/10.1038/35000501


#######################
#          1          #
####################### 

# Load Eisen dataset

data <- read.table("/Users/Conner/Desktop/eisen.txt", header=TRUE, row.names=1, blank.lines.skip=FALSE, na.strings=”NA”)


#######################
#          2          #
####################### 

# Load the class label file

classes <- read.table("/Users/Conner/Desktop/eisenClasses.txt", header=TRUE)


#######################
#          3          #
####################### 

# Subset DF with class labels so you know where one class ends and the other begins

cl <- as.character(classes[,2]) 
dat <- data[,cl]


#######################
#          4          #
####################### 

# Pick a gene, remove cells that have N/A and plot the values for both classes with a
# 1) boxplot
# 2) histogram

gc <- cl[1:19]
act <- cl[20:39]
x <- as.numeric(dat[12,gc]) y <- as.numeric(dat[12,act]) x <- x[!is.na(x)]
y <- y[!is.na(y)]

# boxplot
boxplot(x,y,col=c('red','blue'),main='Example gene DLCL.0037 from DLBCL cDNA 2- channel dataset',axes=F,ylab="log2(ratio intensity)")
axis(2)
axis(1,at=c(1,2),c("GC","ACT"))

# histogram
par(mfrow=c(2,1)) 
hist(x,col='red') 
hist(y,col='blue')


#######################
#          5          #
####################### 

# Calculate pooled variance and calculate the minimum sample size necessary to detect a 1.5 fold difference (at 80% power and 99% confidence)

nx <- length(x)
ny <- length(y)
pool.var <- (((nx-1)*var(x)) + ((ny-1)*var(y)))/(nx+ny-2) 
pool.var

dif.1.5fold <- log2(1.5)/sqrt(pool.var)
pl.ss1.5 <- pwr.t.test(d=dif.1.5fold,sig.level=.01,power=0.8,type="two.sample")
pl.ss1.5


#######################
#          6          #
####################### 

# Now calculate the sample size required for the same gene selected in #5 using the empirically determined delta between the two groups, assuming 99% conf and 80% power

dif <- abs(mean(x)-mean(y))/sqrt(pool.var)
delt.ss1.5 <- pwr.t.test(d=dif,sig.level=.01,power=0.8,type="two.sample")
delt.ss1.5


#######################
#          7          #
####################### 

# Load ssize and gdata libraries, calculate the standard deviation for each gene in the matrix and plot a histogram of the sd's

library(ssize)
library(gdata)

dat.sd <-apply(data,1,sd,na.rm=T)
hist(dat.sd,n=20, col="cyan", border="blue", main="", xlab="Standard Deviation (for data on the log2 scale)")
dens <- density(dat.sd)
lines(dens$x, dens$y*par("usr")[4]/max(dens$y),col="red",lwd=2)
title("Histogram of Standard Deviations for genes in DLBCL Dataset")


#######################
#          8          #
####################### 

# Calculate and plot a proportion of genes vs. sample size graph to get an idea of the number of genes that have an adequate sample size for 95% confidence,
# effect size 3 (log2 scale), and 80% power.

fold.change <- 3
power <- 0.8
sig.level <- 0.05
all.size <- ssize(sd=dat.sd, delta=log2(fold.change), sig.level=sig.level, power=power) ssize.plot(all.size, lwd=2, col="magenta", xlim=c(1,20))
xmax <- par("usr")[2]-1;
ymin <- par("usr")[3] + 0.05
legend(x=xmax, y=ymin, legend= strsplit( paste("fold change=",fold.change,",", "alpha=", sig.level,",", "power=",power,",", "# genes=", length(exp.sd), sep=''), "," )[[1]], xjust=1, yjust=0, cex=1.0)
title("Sample Size to Detect 3-Fold Change")
