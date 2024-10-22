########################################################################################################################
##############################################      Compare mutation rate       ########################################
########################################################################################################################
library(ggplot2)
library(scales) # to access break formatting functions
library(RColorBrewer)

setwd("C:/Users/xli/Documents/Emily/MA_analysis/MS/0_MA_GBE_final_all_data/mutation.rate.summary")

mut <- read.csv("Mutation_rate_comparition_adjusted.csv",header=T)
myColors <- c("blue", "black", "green", "orange", "Magenta", "red", "darkgreen")
names(myColors) <- levels(mut$Species)
colScale <- scale_colour_manual(name = "Species",values = myColors)
shapeScale <- scale_shape_manual(name = "Species",values = c(16, 17, 15, 3, 8, 22, 8))


# x and y axis are transformed and formatted
pdf('Mut3.pdf', width=7, height=7)

# Show log tick marks
p2 <- ggplot(mut, aes(x = Rpts, y = Mut, colour = Species, shape = Species)) + geom_point(size=2.5) + labs(x="Repeat number",  y = "Mutation rate/locus/cell division")+ 
     scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)), limits = c(-8, -2)) +
     theme_bw(base_size = 16)+theme(legend.text = element_text(face = "italic"))+ annotation_logticks() + colScale + shapeScale 

 p2 + stat_smooth(method="lm", se=TRUE, fill=NA, formula = y ~ poly(x), size=0.8, linetype = "dashed")
		

dev.off()

