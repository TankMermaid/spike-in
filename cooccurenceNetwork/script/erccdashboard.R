rm(list = ls())

# source("https://bioconductor.org/biocLite.R")
# biocLite("erccdashboard")
library(erccdashboard)

data(SEQC.Example)
exDat <- initDat(datType="array", isNorm=FALSE, 
                 exTable=UHRR.HBRR.arrayDat,
                 filenameRoot="testRun", sample1Name="UHRR",
                 sample2Name="HBRR", erccmix="RatioPair", 
                 erccdilution = 1, spikeVol = 50, 
                 totalRNAmass = 2.5*10^(3), choseFDR=0.01)



exDat <- est_r_m(exDat)
exDat <- dynRangePlot(exDat)

exDat <- geneExprTest(exDat)

exDat <- erccROC(exDat)

exDat$Figures$rocPlot
