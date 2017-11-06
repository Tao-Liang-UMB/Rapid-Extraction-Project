## Description
# Takes input mzML spectra and runs them through a MALDIQuant pipeline for processing
# MALDI spectra. Outputs a long tibble containing peak list for each spectra.

## Input
# fileDf - a data frame containing a minimum of these 2 named columns: 
#           -> path - The paths and filenames of the mzML files.
#           -> samp - The samp the file belongs to for import. Files to be 
#                   averaged together should have the same samp
# snr - acceptable SNR ratio for peaks

## Output
# Features - a peak list of detected from the parameters provided

AnalyzeSpectra <- function(fileDf ,snr){
  library(tidyverse)
  library(MALDIquant)
  library(MALDIquantForeign)

  # Read in MzML files --------------------------------------------------------
  spectra <- importMzMl(as.character(fileDf$path))
  
  # Normalize MS intensities ---------------------------------------------------
  spectra2 <- transformIntensity(spectra, method = "sqrt")
  spectra3 <- smoothIntensity(spectra2, method="SavitzkyGolay", halfWindowSize = 20)
  
  # Remove Baseline ------------------------------------------------------------
  spectra4 <- removeBaseline(spectra3, method = "SNIP", iterations = 100)
  
  spectra5 <- calibrateIntensity(spectra4, method = "TIC")
  
  # Align MS -------------------------------------------------------------------
  spectra6 <- alignSpectra(spectra5,reference = spectra5[[1]])
  
  # average technical replicates -----------------------------------------------
  spectra7 <- averageMassSpectra(spectra6, fileDf$samp, method="mean")
  
  # Pick Peaks -----------------------------------------------------------------
  peaks <- detectPeaks(spectra7, SNR = snr, halfWindowSize = 10, method = "MAD")
  peaks <- binPeaks(peaks, tolerance = 0.5)
  
  # Preparing Spectra for dot product ------------------------------------------
  # This small section extracts the sample file name and adjusts it to be the sample name.
  samp <- sapply(spectra7, function(x) metaData(x)$file)
  samp <- gsub(".*\\\\", "", samp)
  samp <- gsub(".mzML$","",samp)
  samp <- samp[1, ]
  samp <- gsub("001","",samp)  #since triplicate, so only select "001" sample for renaming
  samp <- factor(samp) 
  
  # Retrieves identified peaks and intensities as a matrix (row = sample, col = m/z, value = intensity)
  features <- intensityMatrix(peaks, spectra7)
  rownames(features) <- samp
  
  # Tranform wide matrix to long dataframe -------------------------------------
  featTbl <- as.tibble(rownames_to_column(as.data.frame(features))) %>%
      mutate(colistin = ifelse(str_detect(rowname, "WT"), "susceptible", "resistant"),
             samp = str_match(rowname, "/([^/]+) $")[ , 2],
             species = str_match(samp, "^(.+?) ")[ , 2]) %>%
      select(-rowname) %>%
      gather(mz, relInt, -colistin, -samp, -species) %>%
      group_by(samp) %>%
      mutate(mz = as.numeric(mz),
             relInt = relInt / max(relInt),
             relInt = ifelse(colistin == "susceptible", -relInt, relInt))
  
  return(featTbl)
}