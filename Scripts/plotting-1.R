# Load packages 
library(tidyverse)
library(ggplot2)
library(viridis)
library(corrplot)
library(grid)
library(cowplot)
library(ggpubr)


# Load results
wide_df <- read_csv("seasonal-vaccination-h5n1-serology/Data/compiled_titers.csv")
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
  # annotate("text", x = 1960.5, y = 11.75, label = "1957",
  #          size = 2.5, color = "black") +
  # annotate("text", x = 1971.5, y = 11.75, label = "1968",
  #          size = 2.5, color = "black") +
  # annotate("text", x = 1980.5, y = 11.75, label = "1977",
  #          size = 2.5, color = "black") +
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

ggsave("seasonal-vaccination-h5n1-serology/Figures/YOB_histogram.pdf", width = 1.5, height = 1, dpi = 300, scale = 2, bg = "transparent")


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
                  color = "#506E82", # color = "#3c5261"
                  # fill = "#506E82",
                  width = 0.1, height = 0.1,
                  alpha = 0.6, stroke = 0.1) +
      geom_jitter(data = data[data$Highlight,],
                  # shape = 15, size = 1.5,
                  shape = 0, size = 1,
                  color = "#506E82",
                  # fill = "#E69F00", # color = "#ac7700",
                  # fill = "#E69F00",
                  width = 0.1, height = 0.1,
                  alpha = 0.8, stroke = .5) +
      # annotate("text", x = 1961.2, y = label_height, label = "1957",
      #          size = 2.5, color = "black") +
      # annotate("text", x = 1972.2, y = label_height, label = "1968",
      #          size = 2.5, color = "black") +
      # annotate("text", x = 1981.2, y = label_height, label = "1977",
      #          size = 2.5, color = "black") +
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
            # axis.text = element_text(size = 9),
            # axis.title = element_text(size = 9),
            panel.border = element_rect(color = "black", linewidth = 0.6, fill = NA),
            axis.line = element_line(linewidth = 0),
            strip.text = element_text(size = 12),
            plot.margin = margin(10, 30, 10, 10),
            # axis.ticks = element_line(linewidth = 0.25),
            # axis.minor.ticks.length.x = unit(0.6, "mm"),
            # axis.minor.ticks.length.y = unit(0.6, "mm")
      ) +
      labs(x = "Year of birth",
           y = ifelse(visit=="v1", y_label, "")
      )
    
    path <- paste0("seasonal-vaccination-h5n1-serology/Figures/",test,"_",visit,".pdf")
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
                color = "#506E82", # color = "#3c5261"
                # fill = "#506E82",
                width = 0.1, height = 0.1,
                alpha = 0.6, stroke = 0.1) +
    geom_jitter(data = fc_data[fc_data$Highlight,],
                # shape = 15, size = 1.5,
                shape = 0, size = 1,
                color = "#506E82",
                # fill = "#E69F00", # color = "#ac7700",
                # fill = "#E69F00",
                width = 0.1, height = 0.1,
                alpha = 0.8, stroke = .5) +
    # annotate("text", x = 1961.2, y = label_height, label = "1957",
    #          size = 2.5, color = "black") +
    # annotate("text", x = 1972.2, y = label_height, label = "1968",
    #          size = 2.5, color = "black") +
    # annotate("text", x = 1981.2, y = label_height, label = "1977",
    #          size = 2.5, color = "black") +
    scale_x_continuous(breaks = seq(1930, 2010, 20),
                       minor_breaks = c(1940, 1960, 1980, 2000)) +
    coord_cartesian(ylim = c(.0625,16)) +
    scale_y_continuous(trans = "log2",
                       # labels = c(.0625,.125,.25,.5,1,2,4,8,16),
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
  
  path <- paste0("seasonal-vaccination-h5n1-serology/Figures/",test,"_",visit,"-fold_change",".pdf")
  ggsave(path, width = 1.6, height = 1, dpi = 300, scale = 2, bg = "transparent")
  
  # --- Pre/post plot ---
  
  pp_data <- long_df[long_df[, "Test"] == test,]
  pp_data$Visit <- factor(pp_data$Visit, levels = c("v1", "v3"))
  pp_data_log <- pp_data
  pp_data_log$Titer <- log2(pp_data_log$Titer)
  gmt_fc <- 2 ** (mean(pp_data_log[pp_data_log$Visit == "v3", ]$Titer) - mean(pp_data_log[pp_data_log$Visit == "v1", ]$Titer))
  gmt_fc_label <- paste0(sprintf("%.2f", round(gmt_fc,3)),"×")
  
  plot = ggplot(pp_data, aes(x = Visit, y = Titer)) +
    # geom_line(aes(group = Participant),
    #           color = "#506E82", linewidth = 0.2, alpha = 0.4) +
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
    # stat_summary(fun.data = mean_se, geom = "errorbar",
    #              width = 0.12, linewidth = 0.3, color = "black") +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
                 width = 0.12, linewidth = 0.3, color = "black") +
    ggsignif::geom_signif(comparisons = list(c("v1", "v3")),
                          test = "wilcox.test",
                          test.args = list(paired = TRUE),
                          # map_signif_level = TRUE,
                          textsize = 3, size = 0.3,
                          tip_length = 0,
                          color = "black",
                          y_position = sig_y,
                          vjust = -0.2,
                          map_signif_level = function(p) {
                            if (p < 0.001) {
                              # sprintf("p = %.2e", p)
                              sprintf("%.2e", p)
                            } else {
                              # sprintf("p = %.3f", p)
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
          # axis.text = element_text(size = 9),
          # axis.title = element_text(size = 9),
          # panel.border = element_rect(color = "black", fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(linewidth = 0.3),#0.1),
          strip.text = element_text(size = 12),
          plot.margin = margin(10, 30, 10, 10),
          # axis.ticks = element_line(linewidth = 0.3),
    ) +
    scale_x_discrete(labels = c("v1" = "Pre", "v3" = "Post")) +
    labs(x = "Vaccination status",
         y = y_label)
  
  path <- paste0("seasonal-vaccination-h5n1-serology/Figures/",test,"_pre-vs-post",".pdf")
  ggsave(path, width = 1.2, height = 1, dpi = 300, scale = 2, bg = "transparent")
  
  plot <- plot + labs(y="")
  path <- paste0("seasonal-vaccination-h5n1-serology/Figures/",test,"_pre-vs-post_no-Y-label",".pdf")
  ggsave(path, width = 1.2, height = 1, dpi = 300, scale = 2, bg = "transparent")
  
}




# y=H5_stalk vs x=NAI50 with color=NT90
# crest <- c(
#   '#a5cd90', '#a4cc90', '#a3cc91', '#a2cb91', '#a0cb91', '#9fca91', '#9eca91', '#9dc991', '#9cc891', '#9bc891', '#9ac791', '#99c791', '#98c691', '#96c691', '#95c591', '#94c591', '#93c491', '#92c491', '#91c391', '#90c391', '#8fc291', '#8ec291', '#8dc191', '#8bc191', '#8ac091', '#89bf91', '#88bf91', '#87be91', '#86be91', '#85bd91', '#84bd91', '#82bc91', '#81bc91', '#80bb91', '#7fbb91', '#7eba91', '#7dba91', '#7cb991', '#7bb991', '#79b891', '#78b891', '#77b791', '#76b791', '#75b690', '#74b690', '#73b590', '#72b490', '#71b490', '#70b390', '#6fb390', '#6eb290', '#6db290', '#6cb190', '#6bb190', '#6ab090', '#69b090', '#68af90', '#67ae90', '#66ae90', '#65ad90', '#64ad90', '#63ac90', '#62ac90', '#62ab90', '#61aa90', '#60aa90', '#5fa990', '#5ea990', '#5da890', '#5ca890', '#5ba790', '#5ba690', '#5aa690', '#59a590', '#58a590', '#57a490', '#57a490', '#56a390', '#55a290', '#54a290', '#53a190', '#53a190', '#52a090', '#519f90', '#509f90', '#509e90', '#4f9e90', '#4e9d90', '#4e9d90', '#4d9c90', '#4c9b90', '#4b9b90', '#4b9a8f', '#4a9a8f', '#49998f', '#49988f', '#48988f', '#47978f', '#47978f', '#46968f', '#45958f', '#45958f', '#44948f', '#43948f', '#43938f', '#42928f', '#41928f', '#41918f', '#40918f', '#40908e', '#3f8f8e', '#3e8f8e', '#3e8e8e', '#3d8e8e', '#3c8d8e', '#3c8c8e', '#3b8c8e', '#3a8b8e', '#3a8b8e', '#398a8e', '#388a8e', '#38898e', '#37888e', '#37888d', '#36878d', '#35878d', '#35868d', '#34858d', '#33858d', '#33848d', '#32848d', '#31838d', '#31828d', '#30828d', '#2f818d', '#2f818d', '#2e808d', '#2d808c', '#2d7f8c', '#2c7e8c', '#2c7e8c', '#2b7d8c', '#2a7d8c', '#2a7c8c', '#297b8c', '#287b8c', '#287a8c', '#277a8c', '#27798c', '#26788c', '#25788c', '#25778c', '#24778b', '#24768b', '#23758b', '#23758b', '#22748b', '#22748b', '#21738b', '#21728b', '#20728b', '#20718b', '#20718b', '#1f708b', '#1f6f8a', '#1e6f8a', '#1e6e8a', '#1e6d8a', '#1e6d8a', '#1d6c8a', '#1d6c8a', '#1d6b8a', '#1d6a8a', '#1d6a8a', '#1c6989', '#1c6889', '#1c6889', '#1c6789', '#1c6689', '#1c6689', '#1c6589', '#1c6488', '#1c6488', '#1c6388', '#1d6388', '#1d6288', '#1d6188', '#1d6187', '#1d6087', '#1d5f87', '#1d5f87', '#1e5e87', '#1e5d86', '#1e5d86', '#1e5c86', '#1e5b86', '#1f5b86', '#1f5a85', '#1f5985', '#1f5985', '#205885', '#205784', '#205784', '#205684', '#215584', '#215583', '#215483', '#225383', '#225283', '#225282', '#225182', '#235082', '#235081', '#234f81', '#244e81', '#244e80', '#244d80', '#254c80', '#254c7f', '#254b7f', '#254a7f', '#26497e', '#26497e', '#26487e', '#27477d', '#27477d', '#27467c', '#27457c', '#28457c', '#28447b', '#28437b', '#28427a', '#29427a', '#29417a', '#294079', '#294079', '#2a3f78', '#2a3e78', '#2a3d78', '#2a3d77', '#2a3c77', '#2a3b76', '#2b3b76', '#2b3a76', '#2b3975', '#2b3875', '#2b3875', '#2b3774', '#2b3674', '#2c3574', '#2c3573', '#2c3473', '#2c3373', '#2c3272', '#2c3172', '#2c3172'
# )

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
                       labels=y_labels,#if(visit == "v1") y_labels else c("","","",""),
                       limits=c(640*.9,40960/.9)) +
    # coord_cartesian(ylim = y_lim) +
    # scale_color_continuous(trans="log2", low="#506E82", high="#E69F00") +
    scale_color_viridis(trans="log2", option="mako", end=0.8,
                        limits=c(10*.95,320/.95),
                        breaks=c(10,20,40,80,160,320)) +
    guides(fill = guide_colorbar(ticks = FALSE)) +
    # scale_color_gradientn(trans="log2", colors=rev(crest)) +
    theme_classic(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          rect = element_blank(),
          # legend.position = "none",
          axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5),
          axis.text.y = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5, size = 7),
          # axis.text = element_text(size = 9),
          # axis.title = element_text(size = 9),
          panel.border = element_rect(color = "black", linewidth = 0.6, fill = NA),
          axis.line = element_line(linewidth = 0),
          strip.text = element_text(size = 12),
          plot.margin = margin(10, 30, 10, 10),
          # axis.ticks = element_line(linewidth = 0.25),
          # axis.minor.ticks.length.x = unit(0.6, "mm"),
          # axis.minor.ticks.length.y = unit(0.6, "mm"),
    ) +
    guides(x = guide_axis(minor.ticks = TRUE),
           y = guide_axis(minor.ticks = TRUE)) +
    labs(
      x = "NAI50 titer",
      y = ifelse(visit == "v1", "Headless H5 titer", "")
    )
  
  path <- paste0("seasonal-vaccination-h5n1-serology/Figures/H5_stalk-and-NAI50-vs-NT90_",visit,".pdf")
  ggsave(path, width = 1.63, height = 1, dpi = 300, scale = 2, bg = "transparent")
  
  plot_no_legend <- plot + theme(legend.position = "none")
  path <- paste0("seasonal-vaccination-h5n1-serology/Figures/H5_stalk-and-NAI50-vs-NT90_",visit,"-noLegend.pdf")
  ggsave(path, width = 1.27*1.1, height = 1.1, dpi = 300, scale = 2, bg = "transparent")
  
}

legend <- get_legend(plot)
path <- "seasonal-vaccination-h5n1-serology/Figures/H5_stalk-and-NAI50-vs-NT90-separateLegend.pdf"
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
      # x_labels <- c(10,20,40,80,160,320,640,"1,280")
      # x_breaks <- c(10,20,40,80,160,320,640,1280)
      x_labels <- c(10,40,160,640)
      x_breaks <- c(10,40,160,640)
      x_label <- "NAI50 titer"
    } else {
      x_lim_outer <- c(256,163840)
      x_lim <- c(512,81920)
      # x_labels <- c(640,"1,280","2,560","5,120","10,240","20,480","40,960","81,920")
      # x_breaks <- c(640,1280,2560,5120,10240,20480,40960,81920)
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
      # geom_jitter(col="#506E82", shape=16, width = 0.1, height = 0.1, alpha = 0.9) +
      geom_jitter(data = data[!data$Highlight,],
                  shape = 16, size = 1.5,
                  color = "#506E82", # color = "#3c5261"
                  # fill = "#506E82",
                  width = 0.1, height = 0.1,
                  alpha = 0.6, stroke = 0.1) +
      geom_jitter(data = data[data$Highlight,],
                  # shape = 15, size = 1.5,
                  shape = 0, size = 1,
                  color = "#506E82",
                  # fill = "#E69F00", # color = "#ac7700",
                  # fill = "#E69F00",
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
        # color = "#506E82"
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
    
    path <- paste0("seasonal-vaccination-h5n1-serology/Figures/NT90_vs_",test,"_",visit,".pdf")
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

path <- "seasonal-vaccination-h5n1-serology/Figures/Titer_correlation-plot_v1.pdf"
pdf(file = path)
corrplot(
  cor_matrix_v1,
  method = "color",
  type = "lower",
  # tl.cex = 0.8,
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

path <- "seasonal-vaccination-h5n1-serology/Figures/Titer_correlation-plot_v3.pdf"
pdf(file = path)
corrplot(
  cor_matrix_v3,
  method = "color",
  type = "lower",
  # tl.cex = 0.8,
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

path <- "seasonal-vaccination-h5n1-serology/Figures/Titer_correlation-plot_legend.pdf"
pdf(file = path, width=2, height=4)
plot.new()

# colorlegend(
#   colbar = colorRampPalette(c("#FFFFFF", "#506E82"))(200),
#   labels = seq(0,1,.2)
#   # labels = seq(-1, 1, length.out = 5), 
#   # at = seq(-1, 1, 0.5), 
#   # xlim = c(2, 4), 
#   # ylim = c(-1, 1), 
#   # vertical = TRUE, 
#   # addlabels = TRUE
# )

# plot.new()
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
# mtext("Spearman\ncorrelation", side = 3, line = 0, cex = 1, font = 1)  # label above the bar

dev.off()









path <- "seasonal-vaccination-h5n1-serology/Figures/H5_stalk-and-NAI50-vs-NT90_customLegend.pdf"
pdf(file = path, width=2, height=4)
plot.new()

# colorlegend(
#   colbar = colorRampPalette(c("#FFFFFF", "#506E82"))(200),
#   labels = seq(0,1,.2)
#   # labels = seq(-1, 1, length.out = 5), 
#   # at = seq(-1, 1, 0.5), 
#   # xlim = c(2, 4), 
#   # ylim = c(-1, 1), 
#   # vertical = TRUE, 
#   # addlabels = TRUE
# )

# plot.new()
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
# mtext("NT90", side = 3, line = 0, cex = 1, font = 1)  # label above the bar

dev.off()












### model beta coefficients ###
model_table <- read.csv("model2_table.csv")
plot_df <- model_table %>%
  filter(term != "(Intercept)") %>%
  mutate(
    p.value = as.numeric(p.value),
    neg_log10_p = -log10(p.value),
    significant = p.value < 0.05,
    term = term
  )

plot_df$term = c(
    "Headless H5 titer", "NAI50 titer", "Visit (post-vaccination)", "Sex (male)", "Year of birth"
    )
plot_df$term = factor(plot_df$term, levels=plot_df %>% arrange(neg_log10_p) %>% pull(term))

# plot_df$term = factor(plot_df$term, levels=rev(plot_df$term))

ggplot(plot_df, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.15, color = "black") +
  # non-significant points: plain open circles, no color mapping
  # significant points: filled, color-scaled by -log10(p)
  geom_point(data = plot_df,
             aes(fill = neg_log10_p),
             shape = 21, size = 3.5, stroke = 1, color = "black") +
  scale_fill_gradient(
    low = "white", high = "#506E82",
    limits = c(0, NA),   # floor at the p=0.05 threshold
    name = expression(italic(P)),
    breaks = c(0,1,2,3,4,5),
    labels = 10**(-c(0,1,2,3,4,5)),
    na.value = "black"       # in case anything exceeds a set upper limit
  ) +
  labs(
    x = expression(beta ~ "(95% CI)"),
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

path <- "testModelCoefficients.pdf"
ggsave(path, width=2.5, height=1, dpi=300, scale=2, bg="transparent")


path <- "testModelCoefficients-customLegend.pdf"
pdf(file = path, width=2, height=4)
plot.new()

colorlegend(
  colbar = colorRampPalette(c("#FFFFFF", "#506E82"))(200),
  labels = 10**(-c(0,1,2,3,4,5)),
  vertical = TRUE,
  ratio.colbar = 0.25,
  xlim = c(0, 1),
  ylim = c(0, 1),
  align = "l",
  cex = 0.8
)

dev.off()




