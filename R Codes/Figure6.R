#library("tidyverse") # includes ggplot2
library("reshape2") # needed for melt()
library("ggplot2")

# Generate sample data
load("Figure.RData")
#load("Figure6.RData")

df <- data.frame(index = rep(c(1,2,3), 12),
                  SI = as.vector(t(SPPV_Result[seq(1,34,length=12),])),
                  NO = as.vector(t(SPPV_Result[seq(2,35,length=12),])),
                  SNR = as.vector(t(SPPV_Result[seq(3,36,length=12),])),
                  ruralgender = rep(rep(c("Rural Female", 
                                          "Rural Male", 
                                          "Urban Female", 
                                          "Urban Male"), each=3), 3),
                  measure = rep(c("Confidence", "Anxiety", "Desperation"), each=12))

# Melt converts a dataframe from wide to long format.
# Need to specify a data frame, the id variables (which will be left at their settings)
# and the measured variables (columns of data) to be stacked.
# This is necessary to ensure that only a single column is used for the "y" variable in ggplot.
melted_df <- melt(df, id.vars=c("index", "measure", "ruralgender")) # wide to long

#Transform into Factor to Ensure Ordering
melted_df$measure.f<-factor(melted_df$measure,levels=c("Confidence", "Anxiety", "Desperation"))


# Initialize a ggplot object. Recall that variables depending on the data
# must be wrapped within "aes" whereas static variables (such as alpha, for transparency)
# go outside of the "aes" call.
p6<-ggplot(melted_df, aes(x=index, y=value)) + 
  # Add points, whose color and shape varies with the "variable" column
  geom_point(aes(color=variable, shape=variable), size=2, alpha=0.5) +
  # Provide breakpoints and respective labelings for the x-axis
  scale_x_continuous(breaks=c(1,2,3), labels=c("NT", "CP", "AT")) +
  scale_y_continuous(breaks=c(0,1), labels = c(0, 1), limits = c(0, 1)) +
  # Lay out plots in a grid fomat, with "measure" used as the vertical
  # facet group and "ruralgender" used as the horizontal facet group
  facet_grid(measure.f~ruralgender) +
  # Add horizontal lines at 0.05 and 0.95
  geom_hline(yintercept=c(0.05,0.95), linetype="dashed") +
  # Add labels. Note that to rename the legend, you have to rename both the
  # "shape" and the "color" variables.
  labs(title = "", x = "", y="",
       shape="", color="") +
  # Choose a grayscale palette
  scale_color_grey() +
  # Remove the default grey background
  theme_bw() +
  # Customizations. Center plot title, and remove background lines
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")
#pdf("SPPV.pdf", width = 8.5, height = 5)
p6
#dev.off()

