library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
# spectra of new method v2, those spectra are generated from 10min, pH 4 condition
#Use pearson's method for spectrum similarity calculation

source("R scripts/AnalyzeSpectra heatmap.R")   # change to your path
# this R scripts calculate the average spectra of techniqical replicates
set.seed(1234)

# Imports spectra, normalizes and performs peak picking. Check 
# <AnalyzeSpectra.R> for details

features <- AnalyzeSpectra("Data/New method v2 heatmap 30min", snr = 8)  # change to your path

# rename the rowname of features, to simply the bacterial species name
rownames(features) <- c("A. baumannii *",
                        "A. baumannii",
                        "C. sakazakii**",
                        "E.coli *",
                        "E.coli",
                        "F.novicida 25°C",
                        "F.novicida 37°C",
                        "K. pneumoniae *",
                        "K. pneumoniae",
                        "K. pneumoniae **",
                        "M. morganii**",
                        "P. aeruginosa *",
                        "P. aeruginosa",
                        "P. mirabilis**",
                        "S. typhi",
                        "S. marcescens**",
                        "Y. pestis 25°C",
                        "Y. pestis 37°C",
                        "S. aureus",
                        "C. albicans")

# Pearson correlation ------------------------------------------------------------------
features.t<-t(features)
sim <- cor(features.t, features.t, method = "pearson")  # correlation calculation

#-----------------------------------------------------------------------------------
# convert to data.frame
df <- as.data.frame(sim)
row_names <- rownames(features)
df$names <- factor(row.names(df), levels = row_names)

# Making a heatmap -------------------------------------------------------------
# Reformats wide dataframe to long form for plotting with ggplot
dfMelt <- melt(df)
dfMelt$names <- factor(dfMelt$names)
dfMelt$variable <- factor(dfMelt$variable)
dfMelt$names <- factor(dfMelt$names, levels = row_names)
dfMelt$variable <- factor(dfMelt$variable, levels = row_names)

# The Plotting
p <- ggplot(dfMelt, aes(variable,names)) + 
  geom_tile(aes(fill=value), color = "black", size = 1) +
  scale_fill_gradientn(name = "", limits=c(-0.15,1),
                       colours = colorRampPalette(c("white","white", "coral", "black"))(256),
                       values =rescale(c(-0.15, 0, 0.5, 1))) +
  theme(axis.text.x = element_text(angle=55, hjust=1),
        axis.text = element_text(color="black", , face = "bold.italic"), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 18, color = "black", face = "bold")) +
  guides(fill = guide_colorbar(ticks = F )) +
  #scale_x_discrete(labels = labels) +
  #scale_y_discrete(labels = labels) +
  coord_fixed(ratio = 1)  
p

ggsave("Figures/New method v2 Pearson's heatmap 30min G+-fungi.tiff",width=10,height=10,units="in")


# file output seection----------------------------------------------------------------------
# needs to be flipped around to match heatmap
dfwrite <- as.data.frame(sim[nrow(sim):1, ])  # dfwrite is a variable name
write.csv(dfwrite, "Results/new method v2 Pearson's heatmap 30min G+-fungi.csv")  #change to your path
