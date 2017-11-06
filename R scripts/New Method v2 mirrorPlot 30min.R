# this script is for mirror plot generation.
library(tidyverse)
library(stringr)

source("R scripts/AnalyzeSpectra mirrorPlot.R")
source("R scripts/makeMirrorPlot.R")
set.seed(1234)

# Imports spectra, normalizes and performs peak picking. Check 
# <AnalyzeSpectra.R> for details
dir <- "Data/mirror Plot mzML/30min mirror Plot"
suffix <- "_30min"

fileDf <- data.frame(path = list.files(dir, pattern=".mzML$",full.names=T)) %>%
    mutate(samp = str_match(as.character(path), "/([^/]+) [^/]+\\.mzML$")[ , 2],
           species = str_match(samp, "^(.+?) ")[ , 2]) 

speciesKey <- list("Kp" = "K. pneumoniae (mcr1)",
                   "Ab" = "A. baumannii",
                   "Ecoli" = "E. coli",
                   "Pa" = "P. aeuruginosa")

features <- fileDf %>% 
    group_by(species) %>%
    do(AnalyzeSpectra(., snr = 5)) %>%
    mutate(speciesLab = map_chr(species, ~ speciesKey[[.]]),
           speciesLab = ifelse(str_detect(samp, "ColR"), "K. pneumoniae (Ara4N)", speciesLab))

# make copy of Kp WT
features <- features %>% 
    filter(str_detect(samp, "Kp.+WT")) %>%
    mutate(speciesLab = "K. pneumoniae (Ara4N)") %>%
    full_join(features)

# Create plots
plots <- features %>%
    group_by(speciesLab) %>%
    do(plot = makeMirrorPlot(.)) %>%
    mutate(resultName = str_replace_all(speciesLab, "\\.* ", "_"),
           resultName = str_replace_all(resultName, "[\\(\\)]", ""))

# Save Plots
plots %>%
    do(ignore = ggsave(paste0("results/", .$resultName, suffix, ".tiff"),
              plot = .$plot, width = 140, height = 70, unit = "mm"))

# How to access plots individually:
plots$plot[[1]] # Ab

# or plot them all:
map(plots$plot, ~ .)
