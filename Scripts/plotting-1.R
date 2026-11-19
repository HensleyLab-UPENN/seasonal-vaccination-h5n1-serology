# Load packages 
library(tidyverse)
library(ggplot2)
library(viridis)
library(corrplot)
library(grid)
library(cowplot)
library(ggpubr)


# Load results
wide_df <- read_csv("../Data/compiled_titers.csv")
long_df <- pivot_longer(wide_df, cols = -c(Participant, Sex, YOB, Visit),
                             names_to = "Test",
                             values_to = "Titer")

highlight_participants <- c(2, 59, 65, 164, 210)

wide_df$Highlight <- wide_df$Participant %in% highlight_participants
long_df$Highlight <- long_df$Participant %in% highlight_participants

fc_df <- long_df %>%
  filter(Visit %in% c("v1", "v3")) %>%
  pivot_wider(
    id_cols = c(Participant, Sex, YOB, Test, Highlight),
    names_from = Visit,
    values_from = Titer
  ) %>%
  mutate(FoldChange = v3 / v1) %>%
  dplyr::select(-v1, -v3)



yobs <- wide_df[,'YOB']

# plot YOBs as a histogram
plot <- ggplot(yobs, aes(x=YOB)) +
  geom_histogram(
    binwidth=2,
    boundary=min(yobs$YOB),
    fill='#506E82',
    alpha=0.4,
    color='#3c5261',
    linewidth=0.25
  ) +
  scale_x_continuous(
    limits = c(1944, 2004),
    breaks = seq(1930, 2020, 20),
    minor_breaks = c(1940, 1960, 1980, 2000, 2020),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 12.5),
    breaks = seq(0,12,2)) +
  coord_cartesian(ylim = c(0,12.5)) +
  guides(x = guide_axis(minor.ticks = TRUE)) +
  labs(x = 'Year of birth', y = 'Number of individuals') +
  geom_vline(xintercept = c(1957, 1968, 1977), linetype ="dashed", 
             color = "black", linewidth = 0.2) +
  theme_classic(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        rect = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 7),
        # axis.text = element_text(size = 10),
        # axis.title = element_text(size = 10),
        panel.border = element_rect(color = "black", linewidth = 0.6, fill = NA),
        axis.line = element_line(linewidth = 0),
        strip.text = element_text(size = 12),
        plot.margin = margin(10, 30, 10, 10),
        # axis.ticks = element_line(linewidth = 0.25),
        # axis.minor.ticks.length.x = unit(0.6, "mm"),
        # axis.minor.ticks.length.y = unit(0.6, "mm")
  )

ggsave("../Figures/YOB_histogram.pdf", width = 1.5, height = 1, dpi = 300, scale = 2, bg = "transparent")


tests <- unique(long_df$Test)
visits <- unique(long_df$Visit)

# for each test, plot V1 and V3 titer data
# then, plot titer fold-change data
for (test in tests) {
  if (grepl("NAI50", test)) {
    y_lim <- c(8, 2048)
    labels <- c(10,20,40,80,160,320,640,"1,280")
    breaks <- c(10,20,40,80,160,320,640,1280)
    y_label <- "NAI50 titer"
  } else if (grepl("NT90", test)) {
    y_lim <- c(8, 512)
    labels <- c(10,20,40,80,160,320)
    breaks <- c(10,20,40,80,160,320)
    y_label <- "H5N1 NT90"
  } else {
    y_lim <- c(512,81920)
    labels <- c(640,"1,280","2,560","5,120","10,240","20,480","40,960","81,920")
    breaks <- c(640,1280,2560,5120,10240,20480,40960,81920)
    y_label <- paste(
      case_match(
        test,
        "H5_FL" ~ "Full-length H5",
        "H5_stalk" ~ "Headless H5",
        "H1_stalk" ~ "Headless H1"
      ),
      "titer")
  }
  sig_y <- (log2(max(y_lim)) - log2(min(y_lim))) * .9 + log2(min(y_lim))
  fc_y <- 2**((log2(max(y_lim)) - log2(min(y_lim))) * .875 + log2(min(y_lim)))
  
  label_height <- max(y_lim)
  
  for (visit in visits) {
    
    data <- long_df[(long_df[, "Test"] == test) & (long_df[, "Visit"] == visit),]
    
    plot = ggplot(data, aes(x = YOB, y = Titer)) +
      geom_vline(xintercept = c(1957, 1968, 1977), linetype ="dashed", 
                 color = "black", linewidth = 0.2) +
      geom_smooth(col = "#506E82", fill = "#506E82",
                  alpha = 0.4, span = 0.5) +
      geom_jitter(data = data[!data$Highlight,],
                  shape = 16, size = 1.5,
                  color = "#506E82",
                  width = 0.1, height = 0.1,
                  alpha = 0.6, stroke = 0.1) +
      geom_jitter(data = data[data$Highlight,],
                  shape = 0, size = 1,
                  color = "#506E82",
                  width = 0.1, height = 0.1,
                  alpha = 0.8, stroke = .5) +
      scale_x_continuous(breaks = seq(1930, 2010, 20),
                         minor_breaks = c(1940, 1960, 1980, 2000)) +
      coord_cartesian(ylim = y_lim) +
      scale_y_continuous(trans = "log2",
                         labels = labels,
                         breaks = breaks) + 
      guides(x = guide_axis(minor.ticks = TRUE)) +
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
            panel.border = element_rect(color = "black", linewidth = 0.6, fill = NA),
            axis.line = element_line(linewidth = 0),
            strip.text = element_text(size = 12),
            plot.margin = margin(10, 30, 10, 10),
      ) +
      labs(x = "Year of birth",
           y = ifelse(visit=="v1", y_label, "")
      )
    
    path <- paste0("../Figures/",test,"_",visit,".pdf")
    ggsave(path, width = 1.5, height = 1, dpi = 300, scale = 2, bg = "transparent")
  
  }
  
  fc_data <- fc_df[fc_df[, "Test"] == test,]
  
  plot = ggplot(fc_data, aes(x = YOB, y = FoldChange)) +
    geom_vline(xintercept = c(1957, 1968, 1977), linetype ="dashed", 
               color = "black", linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype ="dashed", 
               color = "black", linewidth = 0.2) +
    geom_smooth(col = "#506E82", fill = "#506E82",
                alpha = 0.4, span = 0.5) +
    geom_jitter(data = fc_data[!fc_data$Highlight,],
                shape = 16, size = 1.5,
                color = "#506E82",
                width = 0.1, height = 0.1,
                alpha = 0.6, stroke = 0.1) +
    geom_jitter(data = fc_data[fc_data$Highlight,],
                shape = 0, size = 1,
                color = "#506E82",
                width = 0.1, height = 0.1,
                alpha = 0.8, stroke = .5) +
    scale_x_continuous(breaks = seq(1930, 2010, 20),
                       minor_breaks = c(1940, 1960, 1980, 2000)) +
    coord_cartesian(ylim = c(.0625,16)) +
    scale_y_continuous(trans = "log2",
                       labels = c(expression(1/16),expression(1/8),expression(1/4),expression(1/2),1,2,4,8,16),
                       breaks = c(.0625,.125,.25,.5,1,2,4,8,16)) + 
    guides(x = guide_axis(minor.ticks = TRUE)) +
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
          axis.line = element_line(linewidth = 0.1),
          strip.text = element_text(size = 12),
          plot.margin = margin(10, 30, 10, 10),
          # axis.ticks = element_line(linewidth = 0.25),
          # axis.minor.ticks.length.x = unit(0.6, "mm"),
          # axis.minor.ticks.length.y = unit(0.6, "mm")
    ) +
    labs(x = "Year of birth",
         y = paste(y_label, "\nfold-change")
    )
  
  path <- paste0("../Figures/",test,"_",visit,"-fold_change",".pdf")
  ggsave(path, width = 1.6, height = 1, dpi = 300, scale = 2, bg = "transparent")
  
  # --- Pre/post plot ---
  
  pp_data <- long_df[long_df[, "Test"] == test,]
  pp_data$Visit <- factor(pp_data$Visit, levels = c("v1", "v3"))
  pp_data_log <- pp_data
  pp_data_log$Titer <- log2(pp_data_log$Titer)
  gmt_fc <- 2 ** (mean(pp_data_log[pp_data_log$Visit == "v3", ]$Titer) - mean(pp_data_log[pp_data_log$Visit == "v1", ]$Titer))
  gmt_fc_label <- paste0(sprintf("%.2f", round(gmt_fc,3)),"×")
  
  plot = ggplot(pp_data, aes(x = Visit, y = Titer)) +
    geom_jitter(data = pp_data[!pp_data$Highlight,],
                shape = 16, size = 1.5,
                color = "#506E82",
                width = 0.08, height = 0,
                alpha = 0.6, stroke = 0.1) +
    geom_jitter(data = pp_data[pp_data$Highlight,],
                shape = 0, size = 1,
                color = "#506E82",
                width = 0.08, height = 0,
                alpha = 0.8, stroke = .5) +
    stat_summary(fun = mean, geom = "crossbar",
                 width = 0.25, linewidth = 0.3, color = "black") +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
                 width = 0.12, linewidth = 0.3, color = "black") +
    ggsignif::geom_signif(comparisons = list(c("v1", "v3")),
                          test = "wilcox.test",
                          test.args = list(paired = TRUE),
                          textsize = 3, size = 0.3,
                          tip_length = 0,
                          color = "black",
                          y_position = sig_y,
                          vjust = -0.2,
                          map_signif_level = function(p) {
                            if (p < 0.001) {
                              sprintf("%.2e", p)
                            } else {
                              sprintf("%.3f", p)
                            }
                          }
                          ) +
    annotate("text", label=gmt_fc_label, x=1.5, y=fc_y, size=3) +
    coord_cartesian(ylim = y_lim) +
    scale_y_continuous(trans = "log2",
                       labels = labels,
                       breaks = breaks) +
    theme_classic(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          rect = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5),
          axis.text.y = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5, size = 7),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 0.3),#0.1),
          strip.text = element_text(size = 12),
          plot.margin = margin(10, 30, 10, 10),
    ) +
    scale_x_discrete(labels = c("v1" = "Pre", "v3" = "Post")) +
    labs(x = "Vaccination status",
         y = y_label)
  
  path <- paste0("../Figures/",test,"_pre-vs-post",".pdf")
  ggsave(path, width = 1.2, height = 1, dpi = 300, scale = 2, bg = "transparent")
  
  plot <- plot + labs(y="")
  path <- paste0("../Figures/",test,"_pre-vs-post_no-Y-label",".pdf")
  ggsave(path, width = 1.2, height = 1, dpi = 300, scale = 2, bg = "transparent")
  
}





sorted_wide_df <- wide_df %>% arrange(NT90) # sort so that higher NT90 values are drawn on top

x_labels <- c(20,80,320,"1,280")
x_breaks <- c(20,80,320,1280)

y_labels <- c(640,"2,560","10,240","40,960")
y_breaks <- c(640,2560,10240,40960)

for (visit in visits) {
  data <- sorted_wide_df[sorted_wide_df[, "Visit"] == visit,]
  
  plot = ggplot(data, aes(x=NAI50, y=H5_stalk, color=NT90)) +
    geom_jitter(data=data[!data$Highlight,],
                shape=16, width = 0.1, height = 0.1, alpha = 0.9, size=2) +
    geom_jitter(data=data[data$Highlight,],
                shape=15, width = 0.1, height = 0.1, alpha = 0.9, size=2) +
    scale_x_continuous(trans="log2",
                       breaks=x_breaks,
                       labels=x_labels,
                       limits=c(10*.9,1280/.9)) +
    scale_y_continuous(trans="log2",
                       breaks=y_breaks,
                       labels=y_labels,
                       limits=c(640*.9,40960/.9)) +
    scale_color_viridis(trans="log2", option="mako", end=0.8,
                        limits=c(10*.95,320/.95),
                        breaks=c(10,20,40,80,160,320)) +
    guides(fill = guide_colorbar(ticks = FALSE)) +
    theme_classic(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          rect = element_blank(),
          axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5),
          axis.text.y = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5, size = 7),
          panel.border = element_rect(color = "black", linewidth = 0.6, fill = NA),
          axis.line = element_line(linewidth = 0),
          strip.text = element_text(size = 12),
          plot.margin = margin(10, 30, 10, 10),
    ) +
    guides(x = guide_axis(minor.ticks = TRUE),
           y = guide_axis(minor.ticks = TRUE)) +
    labs(
      x = "NAI50 titer",
      y = ifelse(visit == "v1", "Headless H5 titer", "")
    )
  
  path <- paste0("../Figures/H5_stalk-and-NAI50-vs-NT90_",visit,".pdf")
  ggsave(path, width = 1.63, height = 1, dpi = 300, scale = 2, bg = "transparent")
  
  plot_no_legend <- plot + theme(legend.position = "none")
  path <- paste0("../Figures/H5_stalk-and-NAI50-vs-NT90_",visit,"-noLegend.pdf")
  ggsave(path, width = 1.27*1.1, height = 1.1, dpi = 300, scale = 2, bg = "transparent")
  
}

legend <- get_legend(plot)
path <- "../Figures/H5_stalk-and-NAI50-vs-NT90-separateLegend.pdf"
pdf(file = path, width = 1, height = 2)
grid.newpage()
grid.draw(legend)
dev.off()



# y=NT90 vs x=titer for all nonNT_tests
nonNT_tests <- unique(long_df[long_df$Test!="NT90",]$Test)
y_breaks <- c(10,20,40,80,160,320)
y_labels <- c(10,20,40,80,160,320)
y_lim_outer <- c(4, 1024)
y_lim <- c(8, 512)
for (visit in visits) {
  for (test in nonNT_tests) {
    if (grepl("NAI50", test)) {
      x_lim_outer <- c(4, 4096)
      x_lim <- c(8, 2048)
      x_labels <- c(10,40,160,640)
      x_breaks <- c(10,40,160,640)
      x_label <- "NAI50 titer"
    } else {
      x_lim_outer <- c(256,163840)
      x_lim <- c(512,81920)
      x_labels <- c(640,"2,560","10,240","40,960")
      x_breaks <- c(640,2560,10240,40960)
      x_label <- paste(
        case_match(
          test,
          "H5_FL" ~ "Full-length H5",
          "H5_stalk" ~ "Headless H5",
          "H1_stalk" ~ "Headless H1"
        ),
        "titer")
    }
    
    test <- sym(test)
    data <- wide_df[wide_df[, "Visit"] == visit,]
    
    plot <- ggplot(data, aes(x=!!test, y=NT90)) +
      geom_jitter(data = data[!data$Highlight,],
                  shape = 16, size = 1.5,
                  color = "#506E82",
                  width = 0.1, height = 0.1,
                  alpha = 0.6, stroke = 0.1) +
      geom_jitter(data = data[data$Highlight,],
                  shape = 0, size = 1,
                  color = "#506E82",
                  width = 0.1, height = 0.1,
                  alpha = 0.8, stroke = .5) +
      geom_smooth(
        method = "lm",
        col = "#506E82",
        fill = "#506E82",
        alpha = 0.4,
        span = .6
      ) +
      stat_cor(
        size = 2.7,
        label.y = log2(max(y_breaks)),
        label.x = log2(min(x_breaks)),
        method = "spearman"
      ) +
      scale_x_continuous(trans="log2",
                         breaks=x_breaks,
                         labels=x_labels,
                         limits=x_lim_outer) +
      scale_y_continuous(trans="log2",
                         breaks=y_breaks,
                         labels=y_labels,
                         limits=y_lim_outer) +
      coord_cartesian(ylim = y_lim, xlim = x_lim) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            rect = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5),
            axis.text.y = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5, size = 7),
            axis.text = element_text(size = 9),
            axis.title = element_text(size = 9),
            panel.border = element_rect(color = "black", fill = NA),
            axis.line = element_line(linewidth = 0.1),
            strip.text = element_text(size = 12),
            plot.margin = margin(10, 30, 10, 10),
            axis.ticks = element_line(linewidth = 0.3),
            axis.minor.ticks.length.x = unit(0.6, "mm"),
            axis.minor.ticks.length.y = unit(0.6, "mm"),
      ) +
      guides(x = guide_axis(minor.ticks = TRUE)) +
      labs(
        x = x_label,
        y = "NT90"
      )
    
    path <- paste0("../Figures/NT90_vs_",test,"_",visit,".pdf")
    ggsave(path, width = 1.4, height = 1, dpi = 300, scale = 2, bg = "transparent")
    
  }
}


# corrplots
wide_df_v1_data <- wide_df[wide_df$Visit == 'v1',] %>% dplyr::select("NT90","H1_stalk","H5_stalk","H5_FL","NAI50")
wide_df_v3_data <- wide_df[wide_df$Visit == 'v3',] %>% dplyr::select("NT90","H1_stalk","H5_stalk","H5_FL","NAI50")

new_names <- c("H5N1 NT90", "Headless H1", "Headless H5", "Full-length H5", "NAI50")
colnames(wide_df_v1_data) <- new_names
colnames(wide_df_v3_data) <- new_names

cor_matrix_v1 <- cor(wide_df_v1_data, method = "spearman")
cor_matrix_v3 <- cor(wide_df_v3_data, method = "spearman")

testRes_v1 = cor.mtest(wide_df_v1_data, conf.level = 0.95)
testRes_v3 = cor.mtest(wide_df_v3_data, conf.level = 0.95)

path <- "../Figures/Titer_correlation-plot_v1.pdf"
pdf(file = path)
corrplot(
  cor_matrix_v1,
  method = "color",
  type = "lower",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 1.6,          # row/column labels (was 0.8)
  number.cex = 2.0,      # correlation coefficients
  cl.cex = 2.0,          # color legend text
  addCoef.col = "black",
  col.lim = c(0, 1),
  cl.length = 6,
  cl.pos = "n",
  diag = FALSE,
  col = colorRampPalette(c("#FFFFFF", "#FFFFFF", "#506E82"))(200),
)
dev.off()

path <- "../Figures/Titer_correlation-plot_v3.pdf"
pdf(file = path)
corrplot(
  cor_matrix_v3,
  method = "color",
  type = "lower",
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 1.6,          # row/column labels (was 0.8)
  number.cex = 2.0,      # correlation coefficients
  cl.cex = 2.0,          # color legend text
  addCoef.col = "black",
  col.lim = c(0, 1),
  cl.length = 6,
  cl.pos = "n",
  diag = FALSE,
  col = colorRampPalette(c("#FFFFFF", "#FFFFFF", "#506E82"))(200)
)
dev.off()

path <- "../Figures/Titer_correlation-plot_legend.pdf"
pdf(file = path, width=2, height=4)
plot.new()

colorlegend(
  colbar = colorRampPalette(c("#FFFFFF", "#506E82"))(200),
  labels = seq(0, 1, 0.2),
  vertical = TRUE,
  ratio.colbar = 0.25,
  xlim = c(0, 1),
  ylim = c(0, 1),
  align = "l",
  cex = 0.8
)

dev.off()









path <- "../Figures/H5_stalk-and-NAI50-vs-NT90_customLegend.pdf"
pdf(file = path, width=2, height=4)
plot.new()

colorlegend(
  colbar = viridis_pal(option = "mako", end=0.8)(200),
  labels = c(10,20,40,80,160,320),
  vertical = TRUE,
  ratio.colbar = 0.25,
  xlim = c(0, 1),
  ylim = c(0, 1),
  align = "l",
  cex = 0.8
)

dev.off()

