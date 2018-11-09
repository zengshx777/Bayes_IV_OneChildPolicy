library("reshape2") # needed for melt()
library("ggplot2")
#Codes Below Produce Figure 4 in the Paper
load("Figure.RData")
#load("Figure4.RData")

#Transform From List to Long Vector
MTE=MTE_up=MTE_low=NULL
for (res in 1:3)
{
  for(k in 1:4)
  {
    MTE<-c(MTE,mte[[res]][[k]][1,])
    MTE_up<-c(MTE_up,mte[[res]][[k]][2,])
    MTE_low<-c(MTE_low,mte[[res]][[k]][3,])
  }
}

df <- data.frame(Grid=rep(seq(0.01,0.99,length=100),12),
                 MTE = MTE,
                 MTE_low = MTE_low,
                 MTE_up = MTE_up,
                 ruralgender = rep(rep(c("Rural Female", 
                                         "Rural Male", 
                                         "Urban Female", 
                                         "Urban Male"), each=100), 3),
                 measure = rep(c("Confidence", "Anxiety", "Desperation"), each=400))
df$measure.f=factor(df$measure,levels=c("Confidence","Anxiety","Desperation"))


p4<-ggplot(df, aes(Grid)) + 
  # Draw MTE lines, as well as the Confidence Interval
  geom_line(aes(y=MTE),size=0.9,linetype="solid") +
  geom_line(aes(y=MTE_up),size=0.9,linetype="dashed") +
  geom_line(aes(y=MTE_low),size=0.9,linetype="dashed") +
  # Lay out plots in a grid fomat, with "measure" used as the vertical
  # facet group and "ruralgender" used as the horizontal facet group
  facet_grid(measure.f~ruralgender) +
  # Add horizontal lines at 0
  geom_hline(yintercept=0, linetype="dotdash",alpha=0.5) +
  labs(title = "", x = "", y="",
       shape="", color="") +
  # Choose a grayscale palette
  scale_color_grey() +
  # Remove the default grey background
  theme_bw() +
  # Customizations. Center plot title, and remove background lines
  theme(plot.title = element_text(hjust = 0.5),
        panel.spacing.x = unit(1, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")
#pdf(file = "Expected_MTE.pdf",height=6,width=8.5)
p4
#dev.off()



