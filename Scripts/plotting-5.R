# Load packages 
library(tidyverse)
library(ggplot2)
library(viridis)
library(corrplot)
library(grid)
library(cowplot)
library(ggpubr)

df <- read.csv("/Users/jort/Documents/R/seasonal-vaccination-h5n1-serology/Data/human_sera_escape_dms.csv")
site_mapping_df <- read.csv("/Users/jort/Documents/R/seasonal-vaccination-h5n1-serology/Data/site_numbering_map_dms.csv")

combined_df <- left_join(df, site_mapping_df, by=c("site" = "reference_site"))

sum_df <- combined_df %>%
  mutate(escape = pmax(escape, 0), ha_domain = ifelse(region == "HA2", "HA2", "HA1")) %>%
  # group_by(antibody, site, sequential_site, mature_H5_site, region) %>%
  summarize(sum = sum(escape), .by = c(antibody, site, sequential_site, mature_H5_site, region, ha_domain)) %>%  #, sequential_site = sequential_site, region = region) %>%
  ungroup()
sum_df

for (p in unique(sum_df$antibody)) {
  p_lab <- str_sub(p, -4)
  
  plot_df <- sum_df %>% filter(antibody == p) %>% arrange(sequential_site)
  plot_df$ln_grp <- cumsum(c(1, diff(plot_df$sequential_site) != 1))
  
  height <- 0.6
  
  plot <- ggplot(plot_df, aes(x=sequential_site, y=sum, group=ln_grp)) +
    geom_point(size=0.4) +
    geom_line(linewidth=0.15) +
    theme_classic(base_size = 12) +
    scale_y_continuous(limits = c(0,10),
                       breaks = seq(0,10,2)) +
    scale_x_continuous(limits = c(0,NA),
                       breaks = NULL) +
    coord_cartesian(xlim=c(20,547),
                    ylim = c(0.41,9.5)) +
    labs(#title = p_lab,
         x = NULL,
         y = "Escape") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold",
                                margin = margin(b = -10)),
      legend.position = "none",
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      strip.background = element_blank(),
      axis.line = element_line(linewidth = 0.35),
    )
  
  if (p_lab == '0002') {
    plot <- plot +
      geom_rect(
        data = plot_df %>% 
          group_by(ha_domain) %>% 
          summarise(xmin = min(sequential_site, na.rm = TRUE),
                    xmax = max(sequential_site, na.rm = TRUE),
                    .groups = "drop"),
        aes(xmin = xmin, xmax = xmax, ymin = 8.25, ymax = 9.75, fill = ha_domain),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.3
      ) +
      scale_fill_manual(values = c("#F59E0B", "#D97706")) +
      geom_text(
        data = plot_df %>% 
          group_by(ha_domain) %>% 
          summarise(x = (min(sequential_site, na.rm = TRUE) + max(sequential_site, na.rm = TRUE))/2,
                    .groups = "drop"),
        aes(x = x, y = 9, label = ha_domain),
        inherit.aes = FALSE,
        size = 3
      )
  }
  
  if (p_lab == '0210') {
    plot <- plot +
      scale_x_continuous(
        limits = c(0, NA),
        breaks = seq(0,500,100)
      ) +
      labs(x = "HA resiude")
    height <- 0.8
  }
  
  # show(plot)
  path <- paste0("seasonal-vaccination-h5n1-serology/Figures/", p, "-dms_escape.pdf")
  ggsave(path, width = 2.5, height = height, dpi = 300, scale = 2, bg = "transparent")
}



path <- "seasonal-vaccination-h5n1-serology/Figures/dms_escape-protein_legend.pdf"
pdf(file = path, width=2.5, height=4)
plot.new()

colorlegend(
  colbar = colorRampPalette(c("#FFFFFF", "#DB5F57"))(200),
  labels = c("No escape", "High escape"),
  vertical = TRUE,
  ratio.colbar = 0.25,
  xlim = c(0, 1),
  ylim = c(0, 1),
  align = "l",
  cex = 0.8
)

dev.off()

