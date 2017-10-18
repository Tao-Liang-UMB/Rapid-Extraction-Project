library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
# test spectra of new method v2, those spectra are generated from 30min, pH 4.0 condition

source("R scripts/Will_AnalyzeSpectra.R")   # change to your path!
set.seed(1234)

# Imports spectra, normalizes and performs peak picking. Check 
# <AnalyzeSpectra.R> for details

features.raw <- AnalyzeSpectra("data/New method v2 Correlation 30min 003 Caroff", snr = 5)  # change to your path
massrange<-as.numeric(colnames(features.raw))         # change to column names to numeric
massrange1<- massrange[massrange>=1400]           # choose the mass range that you want
massrange2<-length(massrange)-length(massrange1)  
features<-features.raw[,massrange2:length(massrange)]  # select the mass range from the processed raw data

# split into New method matrix and Caroff matrix
features1 <- features[1:19,]
features2 <- features[20:38,]   # change according to how many spectra you imported!

# rename the rowname of features, to simply the bacterial species name

rownames(features1) <- c("AI A. baumannii mcr*",
                         "A. baumannii WT",
                         "C. sakazakii**",
                         "E. cloacae",
                        "E.coli mcr*",
                        "E.coli WT",
                        "F.novicida 25C",
                        "F.novicida 37C",
                        "K. pneumoniae mcr*",
                        "K. pneumoniae WT",
                        "K. pneumoniae Colistin R.**",
                        "M. morganii**",
                        "P. aeruginosa mcr*",
                        "P. aeruginosa WT",
                        "P. mirabilis**",
                        "S. typhi",
                        "S. marcescens**",
                        "Y. pestis 25C",
                        "Y. pestis 37C")

rownames(features2) <- c("SA A. baumannii mcr*",
                         "A. baumannii WT",
                         "C. sakazakii**",
                         "E. cloacae",
                         "E.coli mcr*",
                         "E.coli WT",
                         "F.novicida 25C",
                         "F.novicida 37C",
                         "K. pneumoniae mcr*",
                         "K. pneumoniae WT",
                         "K. pneumoniae Colistin R.**",
                         "M. morganii**",
                         "P. aeruginosa mcr*",
                         "P. aeruginosa WT",
                         "P. mirabilis**",
                         "S. typhi",
                         "S. marcescens**",
                         "Y. pestis 25C",
                         "Y. pestis 37C")


# Pearson's correlation/ similarity calculation -----------------------------------------------------
features1t<-t(features1)
features2t<-t(features2)
sim <- cor(features1t, features2t, method = "pearson")

# Use corrplot package for similarity --------------------------------------------------
library(corrplot)
library(Hmisc)
col<- colorRampPalette(c( "white","deepskyblue3","indianred3"))(256)
#col <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
                           # "cyan", "#007FFF", "blue","#00007F"))(256)  # another good color

# Computing the p-value of Pearson's correlations
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# combine New method and Caroff original data together
features12<- rbind2(features1, features2)
features12t<-t(features12)
p.mat.all <- cor.mtest(features12t)  # calculate the p value from above in-house function, has to use original dataset
p.mat<- p.mat.all[1:19, 20:38]  # change according to how many spectra you imported!

# Generate and export the correlation figure
tiff("R figures/correlation figures SNR5 sig001/New Method 30min 003 vs Caroff Pearson.tiff",
width = 12, height = 12, units = 'in', res = 300)  # change to your path

corrplot(sim, method="circle", main=NULL,
        bg="ghostwhite", tl.col="black", tl.srt = 50, col=col, tl.cex = 1.1,
         number.digits = 1,cl.cex = 1.0, mar=c(0,0,1,0), 
         p.mat = p.mat, sig.level = 0.01, insig="blank", addCoef.col = "black", 
         number.cex = 1.1, addgrid.col = "gray23")  # make the plot

dev.off()

# use Use chart.Correlation(): Draw scatter plots
#library("PerformanceAnalytics")
#chart.Correlation(t(features12), histogram=TRUE, pch=19)

#Output files section ---------------------------------------------------------------------------
# for reprudicibility ANALYSIS!!!
# Extract the diagonal correlation efficient values from the sim matrix

diagdata<-as.data.frame(round(diag(sim),2))
row.names(diagdata)<-rownames(features1)
diagdata
write.csv(diagdata, "R output files/New method v2 30min 003 Pearson vs Caroff diagdata SNR5 sig001.csv") #change to your path

# needs to be flipped around to match heatmap
simwrite <- as.data.frame(sim[nrow(sim):1, ])  # dfwrite is a variable name
write.csv(simwrite, "R output files/New method v2 30min 003 Pearson vs Caroff SNR5 sig001.csv")  #change to your path

