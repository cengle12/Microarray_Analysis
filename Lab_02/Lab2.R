# Author: Conner Engle
# Johns Hopkins: AS.410.671.82.SU22

# Data visualization using Spellman data

# Spellman, P. T., Sherlock, G., Zhang, M. Q., Iyer, V. R., Anders, K., Eisen, M. B., Brown, P. O., Botstein, D., & Futcher, B. (1998). 
# Comprehensive identification of cell cycle-regulated genes of the yeast Saccharomyces cerevisiae by microarray hybridization. 
# Molecular biology of the cell, 9(12), 3273–3297. 
# https://doi.org/10.1091/mbc.9.12.3273


#######################
#          1          #
####################### 

# Unzip .txt file and read into R

dat <- read.table("/Users/Conner/Desktop/spellman.txt",header=TRUE, row.names=1, strip.white=TRUE)


#######################
#          2          #
####################### 

# Verify that there are 6,178 genes and 77 arrays/samples

dim(dat)


#######################
#          3          #
####################### 

# Isolate only the cdc15 experiment (samples 23-46)

cdc15 <- dat[23:46]


#######################
#          4          #
####################### 

# Calculate a correlation matrix between the time points (use Pearson's). Title the plot, label axes, and provide a legend. 
# Use value=pairwise.complete.obs since all of these arrays have at least one missing value

cdc15.cor <- cor(cdc15,use="pairwise.complete.obs")

layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow = TRUE))
par(oma=c(5,7,1,1))
cx <- rev(colorpanel(25,"yellow","black","blue"))
leg <- seq(min(cdc15.cor,na.rm=T),max(cdc15.cor,na.rm=T),length=10) 
image(cdc15.cor,main="Correlation plot cdc15 temp senstive mutant data by timepoint",axes=F,col=cx) 
axis(1,at=seq(0,1,length=ncol(cdc15.cor)),label=dimnames(cdc15.cor)[[2]],cex.axis=0.9,l as=2) 
axis(2,at=seq(0,1,length=ncol(cdc15.cor)),label=dimnames(cdc15.cor)[[2]],cex.axis=0.9,l as=2)

image(as.matrix(leg),col=cx,axes=F)
tmp <- round(leg,2) axis(1,at=seq(0,1,length=length(leg)),labels=tmp,cex.axis=1)


#######################
#          5          #
####################### 

# Select the gene YAL002W (VPS8). Impute the missing values with the row mean. Make sure to cast the gene to numeric

cdc15 <- dat[23:46]
gene1 <- as.numeric(cdc15['YAL002W',]) 
gene1[is.na(gene1)] <- mean(gene1, na.rm=TRUE)


#######################
#          6          #
####################### 

# Generate a profile plot of the same gene. Title the plot, label the axes, and provide time points on x-axis for each array. 

gene1.x <- c(1:length(gene1))
plot(gene1.x,gene1,type='n',main="Profile plot of YAL002W gene across cdc15 samples",xlab="Timepoint",ylab="Expression",axes=F) 
points(gene1.x,gene1,col='red',pch='*',cex=3)
lines(gene1.x,gene1,col='red',lwd=2,lty=2) 
axis(side=1,at=c(1:24),labels=c(10,30,50,70,80,90,100,110,120,130,140,150,160,170,180 ,190,200,210,220,230,240,250,270,290),cex.axis=0.4,las=2)
axis(side=2)


#######################
#          7          #
####################### 

# Create a simple shiny app which allows the user to select and correlate any time point vs another time point across all genes. 


# The following code establishes the dataset of all cdc15 info across all genes. Missing values have been imputed row-wise and rows containing
# no cdc information have been removed

cdc15 <- Dataset[23:46] 
for (i in 1:nrow(cdc15)) {
  num <- as.numeric(unlist(cdc15[i,])) 
  num[is.na(num)] <- mean(num, na.rm=TRUE) 
  cdc15[i,] <- num
}
cdc15 <- na.omit(cdc15)


# The following code establishes the ui and server functions allowing users to interogate the data

ui <- fluidPage(
  titlePanel("Relative Gene Expression by cdc15 Timepoint"),
  sidebarLayout( position = "left", 
  sidebarPanel("Choose your timepoint samples and color",
    selectInput('xcol', 'X Variable', names(cdc15),selected=names(cdc15)[[2]]), 
    selectInput('ycol', 'Y Variable', names(cdc15),selected=names(cdc15)[[2]]), 
    selectInput('color', 'Color', choices=c('Red','Blue','Black','Green')),
  ),
  mainPanel(plotOutput("plot1") 
 )
))

server <- function(input,output) { 
  selectedData <- reactive({
    cdc15[,c(input$xcol, input$ycol)] 
  })
  output$plot1 <- renderPlot({
  plot(selectedData(),xlab=input$xcol,ylab=input$ycol,col=input$color,pch=19,cex =0.5)
  points(selectedData(), col=input$color,pch=1, cex=0.25, lwd=4)
  }) 
}

shinyApp(ui= ui, server=server)
