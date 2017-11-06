# Creates a mirror plot with the colistin-resistant spectrum on top and the 
# susceptible spectrum on bottom.
#
# Input:
#  df - A long formatted data frame that must contain the following columns:
#       -> mz - The m/z values to plot
#       -> relInt - The relative intensity at the corresposding mz. relInt from
#                   a susceptible spectrum must be negative values.
#       -> colistin - A factor indicating colistin resistance
#
#  mzRange - A vector of length 2 indicating the m/z range to plot over
#
# Output: A ggplot2 object containing the mirrored spectra.

makeMirrorPlot <- function(df, 
                           mzRange = c(1000, 2400)) {
    
    resLabs <- data.frame(mz = rep(mzRange[2], 2), 
                          relInt = c(1, -1),
                          txt = c("Resistant", "Susceptible"))
    
    # make the plot
    ggplot(df, aes(x = mz, ymax = relInt, ymin = 0)) +
        geom_linerange(size = 0.25, aes(color = colistin)) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text = element_text(color = "black"),
              text = element_text(size = 8),
              legend.position = "none") +
        geom_text(data = resLabs, aes(x = mz, y = relInt,label = txt), 
                  hjust = "inward", vjust = "inward") +
        xlab(expression(italic("m/z"))) +
        ggtitle(df$speciesLab) +
        xlim(mzRange)
}