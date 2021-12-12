library(devtools)
library(roxygen2)
setwd("C:\\Users\\marco.bee\\Dropbox\\LognPareto")

setwd("./LNPar")
document()

setwd('..')
install('LNPar')
library('LNPar')

y <- read.delim('C:\\Users\\marco.bee\\Dropbox\\LognPareto\\WebCodes\\impreseTN2016.csv',header=F,sep=',')
y <- y[,3]
y <- y[y>0]
n <- length(y)
minRank <- 33
nbootMLE <- 2

resFit <- LPfit(y,n-minRank,nbootMLE)
# resTest <- LPtest(100,0,1,0.67,100,90)
