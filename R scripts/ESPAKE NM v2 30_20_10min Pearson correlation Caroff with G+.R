library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
# compare ESKAPE spectra of new method v2 with Caroff spectra, those spectra are generated from 10, 20, 30min, pH 4.0 condition

#-------------------------------------------------------------------------------------------------------
#import and process New method spectra and Caroff spectra
source("R scripts/AnalyzeSpectra new method v2 ESKAPE spectra.R")   # change to your path!
set.seed(1234)
features <- AnalyzeSpectra("Data/NM ESKAPE 10min Pearson Corr Caroff", snr = 8)  # change to your path!!!
massrange<-as.numeric(colnames(features))         # change to column names to numeric
massrange1<- massrange[massrange>=1250]           # choose the mass range that you want
massrange2<-length(massrange)-length(massrange1)  
features<-features[,massrange2:length(massrange)]  # select the mass range from the processed raw data

# split into New method matrix and Caroff matrix
features.Car <- features[1:10,]
features.NM <- features[11:20,]   # change according to how many spectra you imported!
#--------------------------------------------------------------------------------------------------------


# rename the rowname of features, to simply the bacterial species name

rownames(features.Car) <- c("Caroff E. cloacae",
                            "S. aureus",
                            "K. pneumoniae mcr*",
                            "K. pneumoniae WT",
                            "K. pneumoniae Colistin R.**",
                            "A. baumannii mcr*",
                            "A. baumannii WT",
                            "P. aeruginosa mcr*",
                            "P. aeruginosa WT",
                            "E. faecium")

rownames(features.NM) <- c("NM E. cloacae",
                           "S. aureus",
                           "K. pneumoniae mcr*",
                           "K. pneumoniae WT",
                           "K. pneumoniae Colistin R.**",
                           "A. baumannii mcr*",
                           "A. baumannii WT",
                           "P. aeruginosa mcr*",
                           "P. aeruginosa WT",
                           "E. faecium")
#---------------------------------------------------------------------------------------------------

# Pearson's correlation/ similarity calculation
features.Cart<-t(features.Car)
features.NMt<-t(features.NM)
sim <- cor(features.Cart, features.NMt, method = "pearson")
sim<-sim[nrow(sim):1, ] #flipped to change down-to-up orientation

# Use corrplot package for visulaization --------------------------------------------------
library(corrplot)
library(Hmisc)
col<- colorRampPalette(c( "blue","white","firebrick2"))(256)
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
features12<- rbind2(features.Car, features.NM)
features12t<-t(features12)
p.mat.all <- cor.mtest(features12t)  # calculate the p value from above in-house function, has to use original dataset
p.mat<- p.mat.all[1:10, 11:20]  # change according to how many spectra you imported!
p.mat<-p.mat[nrow(p.mat):1, ] #flipped to change down-to-up orientation

# Generate and export the correlation figure
tiff("Figures/ESKAPE New Method vs Caroff Pearson 10min noLabel with G+.tiff",
width = 8, height = 8, units = 'in', res = 300)  # change to your path

corrplot(sim, method="circle", main=NULL,
        bg="ghostwhite", tl.col="black", tl.srt = 50, col=col, tl.cex = 1.1, tl.pos = "n",
         number.digits = 1,cl.cex = 1.0, mar=c(0,0,1,0), cl.lim = c(-0.3,1), addCoef.col = "black",
         p.mat = p.mat, sig.level = 0.01, insig="blank", 
         number.cex = 1.1, addgrid.col = "gray23")

dev.off()

diagdata<-as.data.frame(round(diag(sim[nrow(sim):1, ]),2))
row.names(diagdata)<-rownames(sim[nrow(sim):1, ])
diagdata
write.csv(diagdata, "Results/ESKAPE New Method vs Caroff 10min Pearson diag.csv")
#change to your path
