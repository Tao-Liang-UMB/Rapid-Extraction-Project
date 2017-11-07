library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
# compare ESKAPE spectra of new method v2 with Caroff spectra, those spectra are generated from 30min, pH 4.0 condition

#-------------------------------------------------------------------------------------------------------
#import and process New method spectra and Caroff spectra
source("R scripts/AnalyzeSpectra.R")   # change to your path!
set.seed(1234)

features1 <- AnalyzeSpectra("Data/NM ESKAPE 10min 001 Pearson Caroff prism", snr = 8)
features2 <- AnalyzeSpectra("Data/NM ESKAPE 10min 002 Pearson Caroff prism", snr = 8)
features3 <- AnalyzeSpectra("Data/NM ESKAPE 10min 003 Pearson Caroff prism", snr = 8)
features4 <- AnalyzeSpectra("Data/NM ESKAPE 20min 001 Pearson Caroff prism", snr = 8)
features5 <- AnalyzeSpectra("Data/NM ESKAPE 20min 002 Pearson Caroff prism", snr = 8)
features6 <- AnalyzeSpectra("Data/NM ESKAPE 20min 003 Pearson Caroff prism", snr = 8)
features7 <- AnalyzeSpectra("Data/NM ESKAPE 30min 001 Pearson Caroff prism", snr = 8)
features8 <- AnalyzeSpectra("Data/NM ESKAPE 30min 002 Pearson Caroff prism", snr = 8)
features9 <- AnalyzeSpectra("Data/NM ESKAPE 30min 003 Pearson Caroff prism", snr = 8)

#----below is the function for pearson diag data grabbing-----------------------------------------------
diaggrab <- function(feature){
  
  massrange<-as.numeric(colnames(feature))         
  massrange1<- massrange[massrange>=1250]           # choose the mass range that you want
  massrange2<-length(massrange)-length(massrange1)  
  feature<-feature[,massrange2:length(massrange)]  # select the mass range from the processed raw data
  
  feature.Car <- feature[1:10,]
  feature.NM <- feature[11:20,]   # change according to how many spectra you imported!
  #--------------------------------------------------------------------------------------------------------
  
  feature.Cart<-t(feature.Car)
  feature.NMt<-t(feature.NM)
  sim <- cor(feature.Cart, feature.NMt, method = "pearson")
  sim<-sim[nrow(sim):1, ] #flipped to change down-to-up orientation
  
  diagdata<-as.data.frame(round(diag(sim[nrow(sim):1, ]),2))
  row.names(diagdata)<-rownames(sim[nrow(sim):1, ])
  
  return(diagdata)

}
#------------------------------------------------------------------------------------------------------

diag1<-diaggrab(features1)
diag2<-diaggrab(features2)
diag3<-diaggrab(features3)
diag4<-diaggrab(features4)
diag5<-diaggrab(features5)
diag6<-diaggrab(features6)
diag7<-diaggrab(features7)
diag8<-diaggrab(features8)
diag9<-diaggrab(features9)

all.diag<-cbind(diag1,diag2,diag3,diag4,diag5,diag6,diag7,diag8,diag9)
colnames(all.diag) <- c("10min 001", "10min 002","10min 003",
                        "20min 001","20min 002","20min 003",
                        "30min 001","30min 002","30min 003")

write.csv(all.diag, "Results/NM ESKAPE Pearson All diag.csv")
#change to your path
