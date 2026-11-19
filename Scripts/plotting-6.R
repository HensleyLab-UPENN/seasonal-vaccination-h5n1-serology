# Load packages 
library(tidyverse)
library(ggplot2)
library(viridis)
library(corrplot)
library(grid)
library(cowplot)
library(ggpubr)

df <- read.csv("../Data/mouse-Ab-titers-long.csv", na.strings = character(0))
df$Group <- factor(df$Group, levels=c("Plasma", "HA-depleted", "NA-depleted", "PBS"))
df$Replicate <- factor(df$Replicate)

ggplot(df, aes(x = Group, y = Titer, shape = Exclude, color = Replicate)) +
  geom_jitter(size = 3, width=0.2, height=0) +
  facet_wrap(~Antigen, scales = "free_x") +
  scale_shape_manual(
    values = c("Y" = 4, "N" = 1),  # 4 = X, 1 = unfilled circle
    labels = c("Y" = "False", "N" = "True")
  ) +
  scale_color_manual(values = c("1" = "#000000", "2" = "grey70")) +
  scale_x_discrete(labels = c("Plasma" = "Pooled plasma", "HA-depleted" = "HA-depleted plasma", "NA-depleted" = "NA-depleted plasma", "PBS" = "PBS")) +
  scale_y_continuous(trans = "log2",
                     breaks = c(12.5, 25, 50, 100, 200, 400, 800, 1600),
                     labels = function(x) {
                       sub("\\.0+$", "", as.character(x))
                     }
  ) +
  geom_hline(
    yintercept = 25,
    linetype = 2,
    color = "grey70"
  ) +
  labs(
    x = "Group",
    y = "Endpoint titer",
    shape = "Included in analysis",
    color = "Experiment",
    title = ""
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold",
                              margin = margin(b = -10)),
    legend.position = "right",
    panel.background = element_blank(),
    plot.background  = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    strip.background = element_blank(),
    axis.line = element_line(linewidth = 0.35),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(shape = 1)))

path <- paste0("../Figures/mouse-passive-transfer-titers.pdf")
ggsave(path, width = 2.5, height = 2, dpi = 300, scale = 2, bg = "transparent")

