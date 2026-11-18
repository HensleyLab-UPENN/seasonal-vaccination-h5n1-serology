# Load packages 
library(tidyverse)
library(ggplot2)
library(viridis)
library(corrplot)
library(grid)
library(cowplot)
library(ggpubr)

df <- read.csv("/Users/jort/Documents/R/seasonal-vaccination-h5n1-serology/Data/neut_depletion_efficiencies.csv")

df$Condition <- factor(df$Condition, levels=c("Mock","Stalk","H5 HA","N1 NA","HA+NA"))
df$Depletion <- pmax(df$Depletion, 0)

df <- df %>%
  mutate(Color = Titer == 100)
df$Color = factor(df$Color, levels=c("TRUE", "FALSE"))

for (p in unique(df$Participant)) {
  for (a in unique(df$Antigen)) {
    plot <- ggplot(df %>% filter(Participant == p & Antigen == a), aes(x=Condition, y=Depletion, fill=Color, shape=factor(Replicate))) +
      geom_point(size=3) +
      scale_shape_manual(values = c("1" = 21, "2" = 22), name="Replicate") +
      scale_fill_manual(values = c("TRUE" = "#BBBBBB", "FALSE" = "white"), name="Below LOD", labels = c("TRUE" = "True", "FALSE" = "False")) +
      # geom_text(data = df %>% filter(Participant == p & Antigen == a & Titer == 100),
      #           aes(x=Condition, y=Depletion, label="*"),
      #           inherit.aes = FALSE,
      #           hjust = -1,
      #           vjust = 0.75,
      #           size = 6,
      #           color = "black") +
      # scale_y_continuous(trans="log2", breaks=c(0.78125,1.5625,3.125,6.25,12.5,25,50,100), limits=c(0.390625,200)) +
      scale_y_continuous(limits=c(0,100)) +
      theme_classic(base_size=12) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            rect = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5),
            axis.text.y = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5, size = 7),
            # axis.text = element_text(size = 9),
            # axis.title = element_text(size = 9),
            panel.border = element_rect(color = "black", linewidth = 0.6, fill = NA),
            axis.line = element_line(linewidth = 0),
            strip.text = element_text(size = 12),
            plot.margin = margin(10, 30, 10, 10)
      ) +
      labs(x="Depletion condition", y="% Depletion") +
      scale_x_discrete(labels = c("Mock"="Mock",
                                  "Stalk"="HL H5",
                                  "H5 HA"="FL H5",
                                  "N1 NA"="N1 NA",
                                  "HA+NA"="HA+NA"))
    
    # show(plot)
    path <- paste0("seasonal-vaccination-h5n1-serology/Figures/",p,"_",a,"_neutDepletion.pdf")
    ggsave(path, width = 2, height = 1.5, dpi = 300, scale = 2, bg = "transparent")
    
    }
}


legend_plot <- plot +
  theme(legend.position = "right") +
  guides(
    shape = guide_legend(order=1, override.aes = list(fill = "#FFFFFF", color = "black", size = 3)),
    fill  = guide_legend(order=2, override.aes = list(shape = 21, color = "black", size = 3))
  )

legend <- cowplot::ggdraw(cowplot::get_legend(legend_plot))
path <- "seasonal-vaccination-h5n1-serology/Figures/neutDepletion-separateLegend.pdf"
ggsave(path, legend, width = 2.5, height = 2)
