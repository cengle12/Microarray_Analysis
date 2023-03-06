# Author: Conner Engle
# Johns Hopkins: AS.410.671.82.SU22

# Golub Analysis


#######################
#          1          #
####################### 

# Load golub data in  multtest library. Also load Biobase and annotate libraries.

library(multtest)
data(golub)
library(Biobase)
library(annotate)


#######################
#          2          #
####################### 

# Cast matrtix to DF and label the gene names as numbers

dat <- as.data.frame(golub)
rows <- paste('g',dimnames(dat)[[1]],sep="")
dimnames(dat)[[1]] <- rows


#######################
#          3          #
####################### 

# Get sample labels and set the sample labels to the DF

ann.dat <- golub.cl
dimnames(dat)[[2]] <- ann.dat


#######################
#          4          #
####################### 

# Perform WMW test on all genes in dataset

wilcox.test.all.genes <- function(x,s1,s2) { x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  w.out <- wilcox.test(x1, x2, exact = FALSE, alternative = "two.sided", correct = TRUE,
var.equal = TRUE)
out <- as.numeric(w.out$statistic) return(out)
}

original.wmw.run <- apply(dat, 1, wilcox.test.all.genes, s1 = colnames(dat) == 0, s2 = colnames(dat) == 1)


#######################
#          5          #
####################### 

# Write for loop to iterate 500 times, wehere in each iteration the columns of DF are shuffled, the WMW test is calcuated on all genes, and max test statistic is saved.

max.vals <- list(1:500)
for (i in 1:500){
  colnames(dat) <- sample(colnames(dat))
  new.vals <- apply(dat, 1, wilcox.test.all.genes, s1 = colnames(dat) == 0, s2 = colnames(dat) == 1)
  max.vals[i] <- max(new.vals)
  }


#######################
#          6          #
####################### 

# Once list of test statistics has been created get the 95% value test statistic. Subset original list of values with only those that have higher test statistic
# than the 95% value. Print the gene names and test statistics out.

max.vals <- t(max.vals)
threshold <- length(max.vals)*0.95
threshold.val <- sort(as.character(max.vals))[threshold]
dat.subset <- original.wmw.run[original.wmw.run > threshold.val]

dat.subset


#######################
#          7          #
####################### 

# Compare results to those using empirical Bayes method in limma package. Load library and calculate p-values for the same dataset using the eBayes() function. 

library(limma)

group.1 <- dat[colnames(dat) == 0]
group.2 <- dat[colnames(dat) ==1]
frame <- cbind(a = 1, b = c(rep(1, length(group.1)), rep(0, length(group.2))))

fit.dat <- lmFit(dat,frame)
dat.bayes <- eBayes(fit.dat)

p.bayes <- dat.bayes$p.value[,-1]


#######################
#          8          #
####################### 

# Sort the empirical Bayes p-values and acquire the lowest n p-values, where n is defined as the number of significant test statistics that you found in problem 6. 

p.bayes.sort <- sort(p.bayes)
n <- length(dat.subset)
p.bayes.sort <- sort(p.bayes)
bayes.cutoff <- p.bayes.sort[1:n]
overlaps <- intersect(names(bayes.cutoff),names(dat.subset))

n

length(overlaps)

# Of top 660 genes for each set, 340 of them overlap between the two methods


#######################
#          9          #
####################### 

# Compare the results from a Student's t-test with the empirical Bayes method. First calculate a two sample (two-tailed) Student's t-test on all genes.
# Then extract only those genes with p-value less than 0.01 from this test. Plot the gene p-values < 0.01 for the Student's t-test vs. the same gens in the
# empirical Bayes method.

t.test.all.genes <- function(x, s1, s2){ x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1, x2, alternative="two.sided", var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

t.res <- apply(dat, 1, t.test.all.genes, s1 = colnames(dat) == 0, s2 = colnames(dat) == 1)
t.sig <- t.res[t.res < 0.01]
t.sig <- as.matrix(t.sig)

bayes.data <- as.matrix(p.bayes)
common.names <- intersect(dimnames(t.sig)[[1]],dimnames(bayes.data)[[1]])
bayes.subset <- as.matrix(bayes.data[common.names,])

merged <- cbind(t.sig,bayes.subset)
colnames(merged) <- header

plot(c(1, nrow(merged)), range(merged), type = "n", xaxt ="n", ylab = "P-Value", xlab="Genes")
points(1:nrow(merged), col = "Red", merged[ , 1])
points(1:nrow(merged), col = "Blue", merged[ , 2])
title(main="P-Value Plot\nStudents T-Test Vs. Empirical Bayes")
legend('topright' , colnames(merged), col = c("Red", "Blue"), pch = 15, cex = 0.9)
