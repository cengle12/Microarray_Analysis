# Author: Conner Engle
# Johns Hopkins: AS.410.671.82.SU22

# Analysis of GEO: GSE12050


#######################
#          1          #
####################### 

# Load marray library and 4 GenePix files, extracting forground and background median values from Cy5 and Cy3 channels

library(limma)
library(marray)
dir.path <- "/Users/Conner/Desktop/GSE12050_amend.gpix/"
dat <- read.GenePix(path= dir.path, name.Gf = "F532 Median", name.Gb = "B532 Median", name.Rf = "F635 Median", name.Rb = "B635 Median", name.W ="Flags")


#######################
#          2          #
####################### 

# Normalize each using median global, loess, and print-tip-group loess methods. Plot MvA plots of all 4 arrays comparing no normalization to 3 normalized approaches

sub.1 <- dat[,1]
sub.2 <- dat[,2]
sub.3 <- dat[,3]
sub.4 <- dat[,4]

# Subject 1
sub.1.med <- maNorm(sub.1, norm = c("median"))
sub.1.loess <- maNorm(sub.1, norm = c("loess"))
sub.1.PT <- maNorm(sub.1, norm = c("printTipLoess"))

par(mfrow = c(4, 1))
maPlot(sub.1, main = "Subject 1 - Raw Data", lines.func = NULL, legend.func = NULL)
maPlot(sub.1.med, main = "Subject 1 - Global Median Normalization", lines.func = NULL, legend.func = NULL)
maPlot(sub.1.loess, main = "Subject 1 - Loess Normalization", lines.func = NULL, legend.func = NULL)
maPlot(sub.1.PT, main = "Subject 1 - Print-Tip-Group Loess Normalization", lines.func = NULL, legend.func = NULL)

# Subject 2
sub.2.med <- maNorm(sub.2, norm = c("median"))
sub.2.loess <- maNorm(sub.2, norm = c("loess"))
sub.2.PT <- maNorm(sub.2, norm = c("printTipLoess"))

par(mfrow = c(4, 1))
maPlot(sub.2, main = "Subject 2 - Raw Data", lines.func = NULL, legend.func = NULL)
maPlot(sub.2.med, main = "Subject 2 - Global Median Normalization", lines.func = NULL, legend.func = NULL)
maPlot(sub.2.loess, main = "Subject 2 - Loess Normalization", lines.func = NULL, legend.func = NULL)
maPlot(sub.2.PT, main = "Subject 2 - Print-Tip-Group Loess Normalization", lines.func = NULL, legend.func = NULL)

# Subject 3
sub.3.med <- maNorm(sub.3, norm = c("median"))
sub.3.loess <- maNorm(sub.3, norm = c("loess"))
sub.3.PT <- maNorm(sub.3, norm = c("printTipLoess"))

par(mfrow = c(4, 1))
maPlot(sub.3, main = "Subject 3 - Raw Data", lines.func = NULL, legend.func = NULL)
maPlot(sub.3.med, main = "Subject 3 - Global Median Normalization", lines.func = NULL, legend.func = NULL)
maPlot(sub.3.loess, main = "Subject 3 - Loess Normalization", lines.func = NULL, legend.func = NULL)
maPlot(sub.3.PT, main = "Subject 3 - Print-Tip-Group Loess Normalization", lines.func = NULL, legend.func = NULL)

# Subject 4
sub.4.med <- maNorm(sub.4, norm = c("median"))
sub.4.loess <- maNorm(sub.4, norm = c("loess"))
sub.4.PT <- maNorm(sub.4, norm = c("printTipLoess"))

par(mfrow = c(4, 1))
maPlot(sub.4, main = "Subject 4 - Raw Data", lines.func = NULL, legend.func = NULL)
maPlot(sub.4.med, main = "Subject 4 - Global Median Normalization", lines.func = NULL, legend.func = NULL)
maPlot(sub.4.loess, main = "Subject 4 - Loess Normalization", lines.func = NULL, legend.func = NULL)
maPlot(sub.4.PT, main = "Subject 4 - Print-Tip-Group Loess Normalization", lines.func = NULL, legend.func = NULL)


#######################
#          3          #
####################### 

# Plot density plots of log ratio values for each normalization (and pre normalization) for Subject 4

plot(
  density(na.omit(maM(sub.4))), 
 	main = "Density Plot for Array 4 (pre- and post-normalized)", 
 	ylim=c(0, 0.9),xlim=c(-6.6, 11.7),
 	col="black"
	)
lines(density(na.omit(maM(sub.4.med))), col = "red") 
lines(density(na.omit(maM(sub.4.loess))), col = "blue") 
lines(density(na.omit(maM(sub.4.PT))), col = "green")

legend.key <- c(
 	"Pre Normalization", 
 	"Global Median Normalization", 
 	"Loess Normalization", 
 	"Print-Tip-Group Loess Normalization"
 	)
legend(1, 0.85, legend = legend.key, lty = c(1, 1), lwd = c(2.5, 2.5),col = c("black", "red", "blue", "green"))


#######################
#          4          #
####################### 

# Extract Cy5 foreground values and log2 transform them. Calculate global median normalization on these 4 arrays using background subtracted Cy5 values.
# Median values should be 1 (or 0 when log2 transformed)

for(i in 1:4){
   	id <- paste("subject", i, sep = ".")
   	background <- maRb(dat[ , i])
  	foreground <- maRf(dat[ , i])
   	diff <- foreground - background
   	assign(id, log2(diff))
   }

raw <- cbind(subject.1, subject.2, subject.3, subject.4)
median  <- apply(raw, 2, median, na.rm = T)
norm <- sweep(raw, 2, median)

header <- c('Array.1','Array.2','Array.3','Array.4')
colnames(norm) <- header


#######################
#          5          #
####################### 

# Calculate Spearman's rank correlation between all 4 arrays from #4 and do the same for M values from loess normalized data from #2. Then plot scatter plot matrix
# for each of the two normalizations. Print correlation coefficients to the screen.

M1 <- maM(sub.1.loess)
M2 <- maM(sub.2.loess)
M3 <- maM(sub.3.loess)
M4 <- maM(sub.4.loess)

dat.loess <- cbind(M1,M2,M3,M4)
header <- c('Array.1','Array.2','Array.3','Array.4')
colnames(dat.loess) <- header

dat.loess.spear <- cor(dat.loess, use = "complete.obs", method = "spearman")
dar.spear <- cor(norm, use = "complete.obs", method = "spearman")

dat.loess.spear
dat.spear

pairs(norm, main = "Scatterplot Matrix\nGlobal Median Normalization")
pairs(dat.loess, main = "Scatterplot Matrix\nGlobal Median Normalization")


#######################
#          6          #
####################### 

# Compare normalizations to quantile normalized data. 
# 1) Subtract the forground - background for each of the 4 chips for only the Cy5 channel
# 2) Sort each column independently in this new matrix
# 3) Calculate row means for the sorted matrix
# 4) Create a new matrix with each row having the same values as the sorted row mean vectors
# 5) Rank the columns independently on the original bacground subtracted matrix
# 6) Reorder the columns in the new matrix from step 4 using the ranks from step 5

for(i in 1:4){
 	id <- paste("sub", i, sep = "")
 	background <- maRb(dat[, i])
 	foreground <- maRf(dat[, i])
 	assign(id, foreground - background)
 }

subjects <- cbind(sub1,sub2,sub3,sub4)
colnames(subjects) <- c('Subject1','Subject2','Subject3','Subject4')
sorted <- apply(subjects,2,sort)
r.means <- rowMeans(sorted,na.rm=FALSE)

ord1 <- order(sub1)
ord2 <- order(sub2)
ord3 <- order(sub3)
ord4 <- order(sub4)

for(i in 1:4){
  name <- paste("norm", i, sep = "")
 	assign(name, rep(NA, nrow(subjects)))
}

norm1[ord1] <- r.means
norm2[ord2] <- r.means
norm3[ord3] <- r.means
norm4[ord4] <- r.means
quantile.norm <- cbind(norm1, norm2, norm3, norm4)

# Verify all columns for quantile normalized data have same median value
median(quantile.norm[,1])
median(quantile.norm[,2])
median(quantile.norm[,3])
median(quantile.norm[,4])

# Plot
par(mfrow = c(4, 1))
hist(log2(quantile.norm[, 1]))
hist(log2(quantile.norm[, 2]))
hist(log2(quantile.norm[, 3]))
hist(log2(quantile.norm[, 4]))


#######################
#          7          #
####################### 

# Log2 transform new matrix from step 6 and calculate spearman's ranmk correlation between the 4 arrays and plot a scatter plot matrix. Print correlation coefficients.

log.quant.norm <- log2(quantile.norm)
quant.spearman <- cor(log.quant.norm, use = "complete.obs", method = "spearman")
quant.spearman

pairs(log.quant.norm, main = "Scatterplot Matrix\nQuantile Normalization")


#######################
#          8          #
####################### 

# Using normalized qRT-PCR data, condunct Spearman's rank correlation to find which two patients are most correlated. Plot two patients as a scatter plot

dat <- read.table("/Users/Conner/Desktop/fold_chg_matrix.txt", header=T,row.names=1,fill=TRUE, na.strings = c())
dat.comp <- dat[1:32]
t.dat.comp <- t(dat.comp)
cor(t.dat.comp,use="complete.obs",method=c("spearman"))

# Plot of patient 6 vs. 7
pat.7 <- as.numeric(dat[5,1:32])
pat.6 <- as.numeric(dat[8,1:32])
plot(pat.6, pat.7, main="Scatterplot of Patient 6 vs Patient 7",
     xlab="Patient 6 ", ylab="Patient 7 ", pch=21)

