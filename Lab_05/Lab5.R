# Author: Conner Engle
# Johns Hopkins: AS.410.671.82.SU22

# Differential Expression Analysis - Rat Ketogenic Diet Samples on Affymetrix Arrays (RAE230A)

#######################
#          1          #
####################### 

# Load dataset into R

dat <- read.table(“/Users/Conner/Desktop/rat_KD.txt”, header=TRUE, row.names=1)


#######################
#          2          #
####################### 

# log2 transform data, then use Student's t-test to calculate changing genes between the control diet vs. ketogenic diet classes

log.dat <- log2(dat)
con <- names(log.dat[1:6])
exp <- names(log.dat[7:11])

t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

p.val <- apply(log.dat, 1, t.test.all.genes, s1 = con, s2 = exp)


#######################
#          3          #
####################### 

# Plot a histogram of the p-vals and determine how many probesets have p<0.05 and p<0.01. Then divide an alpha of 0.05 by the total number of probesets and
# determine how many probesets have a p-value less than this value (Bonferroni correction).

hist(p.val,col="lightblue",xlab="p-values",main="P-value dist’n between\nc=Control and Ketogenic Diet groups",cex.main=0.9)
abline(v=.05,col=2,lwd=2)

length(p.val[p.val < 0.05])
# 5160 probesets

length(p.val[p.val <0.01])
# 2414 probesets

bonf.limit <- 0.05/length(p.val)
length(p.val[p.val < bonf.limit])
# 12 probesets


#######################
#          4          #
#######################

# Calculate mean for each gene, and calculate fold change between the groups (control vs. ketogenic diet)

con.m <- apply(log.dat[,con],1,mean,na.rm=T)
exp.m <- apply(log.dat[,exp],1,mean,na.rm=T)
fold <- con.m-exp.m


#######################
#          5          #
#######################

# What is the maximum and minimum fold change value on linear scale. Now determine probesets with a p-val < Bonferroni threshold and |fold change| > 2

2**min(fold)
# 0.082

2**max(fold)
# 55.16

names(p.val[p.val < bonf.limit & 2**abs(fold) >=2])
# 1367553_x_at, 1370239_at, 1370240_x_at, 1371102_x_at, 1371245_a_at, 1388608_x_at


#######################
#          6          #
#######################

# Transform p-val (-1*log10(p-val)) and create a volcano plot with the p-value and fold change vectors. Make sure to use a log10 tranformation for the p-val
# and a log2 transformation for the fold change. Plot lines for each fold change threshold

p.trans <- -1 * log10(p.val
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\ncontrol and Ketogenic Diet group differences')

points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(.05)&fold>log2(2))],fold[(p.trans> - log10(.05)&fold>log2(2))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(.05)&fold< -log2(2))],fold[(p.trans> -log10(.05)&fold< - log2(2))],col=1,bg=3,pch=21)

abline(v= -log10(.05))
abline(h= -log2(2))
abline(h=log2(2))
