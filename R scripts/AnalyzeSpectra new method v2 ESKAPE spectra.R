## Description
# Takes input mzML spectra and runs them through a MALDIQuant pipeline for processing
# MALDI spectra. Outputs a peak list for each spectra.

## Input
# Path - a folder containing MzXML spectra
# snr - acceptable SNR ratio for peaks

## Output
# Features - a peak list of detected from the paraeters provided

AnalyzeSpectra <- function(path,snr){
  
  library(MALDIquant)
  library(MALDIquantForeign)

  # Read in MzML files --------------------------------------------------------
  file.list <- list.files(path, pattern=".mzML$",full.names=T)
  
  spectra <- importMzMl(file.list)
  
  # Normalize MS intensities ---------------------------------------------------
  spectra2 <- transformIntensity(spectra, method = "sqrt")
  spectra3 <- smoothIntensity(spectra2, method="SavitzkyGolay", halfWindowSize = 20)
  
  # Remove Baseline ------------------------------------------------------------
  spectra4 <- removeBaseline(spectra3, method = "SNIP", iterations = 100)
  
  spectra5 <- calibrateIntensity(spectra4, method = "TIC")
  
  # Align MS -------------------------------------------------------------------
  spectra6_0 <- alignSpectra(spectra5,reference = spectra5[[1]])
  
  # average technical replicates -----------------------------------------------
  # change based on how many spectra you imported
  average.sequence<- factor(c("1","2","3","4","5","6","7","8","9","10","11",
                               "12","12","12","13","13","13","14","14","14","15","15","15","16",
                              "16","16","17","17","17","18","18","18","19","19","19","20","20",
                              "20","21","21","21","22","22","22"),
                            levels = c("1","2","3","4","5","6","7","8","9","10",
                                       "11","12","13","14","15","16","17","18","19","20","21","22"))
  spectra6<-averageMassSpectra(spectra6_0, labels = average.sequence, method="mean") # Average technical replicates
  
  # Pick Peaks -----------------------------------------------------------------
  peaks <- detectPeaks(spectra6, SNR = snr, halfWindowSize = 10, method = "MAD")
  peaks <- binPeaks(peaks, tolerance = 0.5)
  
  # Preparing Spectra for dot product ------------------------------------------
  # This small section extracts the sample file name and adjusts it to be the sample name.
  sample <- sapply(spectra6, function(x) metaData(x)$file)
  sample <- gsub(".*\\\\", "", sample)
  sample <- gsub(".mzML$","",sample)
  sample <- factor(sample) 
  
  # Retrieves identified peaks and intensities as a matrix (row = sample, col = m/z, value = intensity)
  features <- intensityMatrix(peaks, spectra6)
  rownames(features) <- sample
  
  return(features)
}