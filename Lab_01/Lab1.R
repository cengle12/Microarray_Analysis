# Author: Conner Engle
# Johns Hopkins: AS.410.671.82.SU22

# Basic R syntax/plots with data solutions


#######################
#          1          #
####################### 

# Unzip .txt file and read into R

dat <- read.table(file="/Users/Conner/Documents/JH/Gene_Expression/Datasets/sle_b_cell.txt", header=T, row.names=1)


#######################
#          2          #
####################### 

# Print table dimensions to the screen to verify 26 samples present

ncol(dat)


#######################
#          3          #
####################### 

# Print sample names to screen

colnames(dat)


#######################
#          4          #
####################### 

# Plot second SLE patient sample versus first normal control samples in an xy scatterplot. Label x and y axes 'Normal' and 'SLE' respectively and title plot.

x <- dat[,18]
y <- dat[,2]
dat.plot <- plot(x,y,xlab=’Normal’,ylab=’SLE’,main=’SLE B cell sample vs. Normal B cell sample – all probesets’)
grid(dat.plot, col='gray')


#######################
#          5          #
####################### 

# Create same plot but pick only the first 20 probesets. Change the shape and color of the points.

x.20 <- dat[1:20,18]
y.20 <- dat[1:20,2]
dat.plot.20 <- plot(x.20,y.20,xlab='Normal',ylab='SLE', main='SLE B cell sample vs. Normal B cell sample – first 20 subsets',pch=15,col='blue')
grid(dat.plot.20, col='gray')


#######################
#          6          #
####################### 

# Plot a gene profile plot of IGLJ3 (probeset ID 211881_x_at)

gene.vector <- as.numeric(dat["211881_x_at",])
gene.profile <- plot(range(1:26),range(gene.vector),type='n', xlab='Samples', ylab='Intensity',main='Gene Profile Plot for IGLJ3')
lines(gene.vector)
grid(gene.profile,col='gray')


#######################
#          7          #
####################### 

# Create boxplot with single distribution per condition. First create factor vector to indicate which sample is disease vs. control

gene.vector <- as.numeric(dat["211881_x_at",])
f <- c(rep("SLE",17),rep("Control",9))
boxplot(gene.vector~f, xlab='Condition', ylab='Intensity',main='Gene Profile by Condition for IGLJ3')
stripchart(gene.vector~f,add = TRUE, vertical=TRUE, method=’jitter’)


