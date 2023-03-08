# Author: Conner Engle
# Johns Hopkins: AS.410.671.82.SU22

# Normalization and Bioconductor

#######################
#          1          #
####################### 

# Load marray library and swirl data set

library(marray)
dat <- swirl


#######################
#          2          #
####################### 

# Plot an MvA plot of array 3 without any stratified lines.

array.3 <- dat[,3]
maPlot(array.3, main='Pre-normalization MvA plot for Array 3', lines.func=NULL, legend.func=NULL)


#######################
#          3          #
####################### 

# Normalize array 3 by global median location normalization

med.norm <- maNorm(array.3, norm=c('median'))


#######################
#          4          #
#######################

# Plot an MvA plot of the normalized array without the stratified lines or legend

maPlot(med.norm, main='Post-normalization MvA plot for Array 3 - Median Normalized', lines.func=NULL, legend.func=NULL)


#######################
#          5          #
#######################

# Repeat #3 and #4 applying loess global intensity normalization

loess.norm <- maNorm(array.3, norm=c('loess'))
maPlot(loess.norm, main='Post-normalization MvA plot for Array 3 - Loess Normalized', lines.func=NULL, legend.func=NULL)


#######################
#          6          #
#######################

# Using a.cdna object, normalize both arrays and provide MvA plots for each array normalized by the following 3 methods:
# 1) no normalization
# 2) print-tip loess normalization
# 3) scale print-tip normalization
# Put all 3 plots for a single patient array on the same page

pat.1 <- a.cdna[,1]
pat.2 <- a.cdna[,2]

pat.1.loess <- maNorm(pat.1, norm=c('printTipLoess'))
pat.1.MAD <- maNorm(pat.1, norm=c('scalePrintTipMAD'))

pat.2.loess <- maNorm(pat.2, norm=c('printTipLoess'))
pat.2.MAD <- maNorm(pat.2, norm=c('scalePrintTipMAD'))

par(mfrow=c(3,1))
maPlot(pat.1, main='Pre-normalization MvA plot for Patient 1', lines.func=NULL, legend.func=NULL)
maPlot(pat.1.loess, main='Post-normalization MvA plot for Patient 1 – Print Tip Loess Normalized', lines.func=NULL, legend.func=NULL)
maPlot(pat.1.MAD, main='Post-normalization MvA plot for Patient 1 – Scale Print Tip Normalized', lines.func=NULL, legend.func=NULL)

par(mfrow=c(3,1))
maPlot(pat.2, main='Pre-normalization MvA plot for Patient 2', lines.func=NULL, legend.func=NULL)
maPlot(pat.2.loess, main='Post-normalization MvA plot for Patient 2 – Print Tip Loess Normalized', lines.func=NULL, legend.func=NULL)
maPlot(pat.2.MAD, main='Post-normalization MvA plot for Patient 2 – Scale Print Tip Normalized', lines.func=NULL, legend.func=NULL)


#######################
#          7          #
#######################

# Create a data matrix that can be written out to a file with 19,200 rows and 2 columns. Using the functions maM(), maGnames(), and maLabels(),
# figure out how to create the data matrix, get the probe IDs, and assign the probe IDs to the row names
# Do this for the 2 normalized metadata objects that you created in #6 above

ids <- a.cdna@maGnames@maLabels
header <- c(‘Patient_1’, ‘Patient_2’)

dat <- cbind(pat.1.loess, pat.2.loess)
dat.loess <- maM(dat)
dimnames(dat.loess)[[2]] <- header
dimnames(dat.loess)[[1]] <- ids
 
dat.2 <- cbind(pat.1.MAD, pat.2.MAD)
dat.MAD <- maM(dat.2)
dimnames(dat.MAD)[[2]] <- header
dimnames(dat.MAD)[[1]] <- ids


#######################
#          8          #
#######################

# Load the following libraries: affy, limma, affydata, affyPLM, and fpc

library(affy)
library(affydata)
library(affyPLM)
library(limma)
library(fpc)


#######################
#          9          #
#######################

# Read in 3 raw Affymetrix .CEL files and normalize them with 2 different algorithms.

dir.path <- "/Users/Conner/Desktop/HGU133plus_files"
fns <- sort(list.celfiles(path=dir.path,full.names=TRUE))
data.affy <- ReadAffy(filenames=fns,phenoData=NULL)


#######################
#          10         #
#######################

# Using the function: expresso in addition to exprs(), create the normalized data matrices with 54,675 rows and 3 columns for the 2 different norm. algorithms

dat.RMA <- expresso(Dilution,bgcorrect.method="rma",normalize.method="quantiles",pmcorrect.met hod="pmonly",summary.method="medianpolish")

dat.MAS <- expresso(data.affy,bgcorrect.method="mas",normalize.method="quantiles",pmcorrect.me thod="pmonly",summary.method="medianpolish")


#######################
#          11         #
#######################

# Now use the cor() function to calculate the correlation between the 3 arrays for both normalized data matrices

cor.rma <- cor(exprs(dat.RMA))
cor.rma

cor.mas <- cor(exprs(dat.MAS))
cor.mas

# As expected, high level of correlation between the 3 healthy subjects across all genes
