library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(WGCNA)
library(ggplot2)
enableWGCNAThreads()
allowWGCNAThreads()


modTrait <- readxl::read_xlsx("./Module_Trait_P35_CX.xlsx", sheet = 1)


mod_plot <- ggplot(modTrait, aes(x=Module,y=Group, label=text, width = 1)) +
        geom_tile(aes(fill=signedLog10fdr),color="grey60") +
        scale_fill_gradient2(low = "blue", high = "red","[beta]\nsigned\n-log10FDR\n") +
        geom_text(size=1.5, color="black")  +
        ggtitle(paste0("soft_power_", 18,"_module_trait_module_size_", 100, "_DeepSplit_",0)) +
        theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=10),
        plot.margin=unit(c(1,1,-0.5,1),"mm"),
        panel.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10)
        ) 

colorbar <- ggplot(data.frame(x = modTrait$Module,y=0.1),aes(x,y)) +
  geom_tile(fill = modTrait$Module, colour = "white") + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        plot.margin=unit(c(1,1,-0.5,1),"mm"),
        text = element_text(size=10),
        axis.text.x = element_text(angle=60,hjust = 1),panel.border = ele
        )

gA <- ggplotGrob(mod_plot)
gB <- ggplotGrob(colorbar)
#gA$widths <- gB$widths
gB$widths <- gA$widths


final_plot <- grid.arrange(gA,gB, layout_matrix = cbind(c(1,1,1,1,1,1,1,1,1,1,1,2,2)))
ggsave(final_plot, filename = "./Module_Trait_P35_CX_sft_18_mm_20_ds_0.pdf",limitsize = F)

g <- arrangeGrob(gA,gB,nrow = 2, heights = c(7/8,1/8))
grid.draw(g)


