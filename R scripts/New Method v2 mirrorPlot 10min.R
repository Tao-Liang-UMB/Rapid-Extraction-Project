library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
# this script is for mirror plot generation.

source("R scripts/AnalyzeSpectra mirrorPlot.R")   # change to your path!
set.seed(1234)

# Imports spectra, normalizes and performs peak picking. Check 
# <AnalyzeSpectra.R> for details

features <- AnalyzeSpectra("Data/mirror Plot mzML/10min mirror Plot", snr = 5)  # change to your path


