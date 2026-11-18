# Load packages 
library(tidyverse)
library(ggplot2)
library(viridis)
library(corrplot)
library(grid)
library(cowplot)
library(ggpubr)
library(survival)
library(survminer)
library(glue)
library(ggtext)

# df <- read.csv("/Users/jort/Documents/R/seasonal-vaccination-h5n1-serology/Data/20260714_mouse_experiment.csv")
df <- read.csv("/Users/jort/Documents/R/seasonal-vaccination-h5n1-serology/Data/20260819_mouse_experiments_combined.csv")
df$group <- factor(
  df$group,
  levels = c("Unabsorbed", "HA absorbed", "NA absorbed", "PBS")
)

fit <- survfit(Surv(day, status) ~ group, data = df)

groups_to_test <- c("Unabsorbed", "HA absorbed", "NA absorbed")  # whatever your non-control groups are called

results <- lapply(groups_to_test, function(g) {
  sub <- subset(df, group %in% c("PBS", g))
  sd <- survdiff(Surv(day, status) ~ group, data = sub)
  p <- 1 - pchisq(sd$chisq, df = 1)
  data.frame(comparison = paste(g, "vs PBS"), chisq = sd$chisq, p_raw = p)
})

results_df <- do.call(rbind, results)
results_df$p_adj <- p.adjust(results_df$p_raw, method = "BH")
results_df

p <- ggsurvplot(
  fit,
  data = df,
  risk.table = FALSE,            # numbers-at-risk table below plot
  censor = FALSE,
  palette = c("#506E82", "#F59E0B", "#009E73", "#454545"),        # matches your usual color scheme
  xlab = "Days post-infection",
  ylab = "% Survival",
  legend.title = "Group",
  ggtheme = theme_classic(base_size = 12),
  # legend.labs = c(glue("Pooled plasma (*P* = {round(results_df$p_adj[1],3)})"),
  #                 glue("HA depleted plasma (*P* = {round(results_df$p_adj[2],3)})"),
  #                 glue("NA depleted plasma (*P* = {round(results_df$p_adj[3],3)})"),
  #                 "PBS")
  legend.labs = c("Pooled plasma",
                  "HA-depleted plasma",
                  "NA-depleted plasma",
                  "PBS")
)

plot <- p$plot +
  scale_x_continuous(
    limits = c(0, 14),
    breaks = c(0, 7, 14),
    minor_breaks = seq(0,14)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = c(0, 25, 50, 75, 100)
  ) +
  guides(x = guide_axis(minor.ticks = TRUE)) +
  theme(
    legend.position = "right",
    panel.background = element_blank(),
    plot.background  = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    strip.background = element_blank(),
    axis.line = element_line(linewidth = 0.35),
    legend.text = element_markdown()
  )

path <- "seasonal-vaccination-h5n1-serology/Figures/mouse_survival.pdf"
ggsave(path, width = 2.15, height = 1, dpi = 300, scale = 2, bg = "transparent")

plot_no_legend <- plot + theme(legend.position = "none")
path <- "seasonal-vaccination-h5n1-serology/Figures/mouse_survival-no_legend.pdf"
ggsave(path, width = 1.2, height = 1, dpi = 300, scale = 2, bg = "transparent")

## spacing issues with this method, use cowplot instead
# legend <- get_legend(plot)
# path <- "seasonal-vaccination-h5n1-serology/Figures/mouse_survival-separateLegend.pdf"
# pdf(file = path, width = 2.5, height = 2)
# grid.newpage()
# grid.draw(legend)
# dev.off()

legend <- cowplot::ggdraw(cowplot::get_legend(plot))
path <- "seasonal-vaccination-h5n1-serology/Figures/mouse_survival-separateLegend.pdf"
ggsave(path, legend, width = 2.5, height = 2)







# bw_df <- read.csv("/Users/jort/Documents/R/seasonal-vaccination-h5n1-serology/Data/20260714_mouse_experiment-body_weight.tsv", sep="\t")
bw_df <- read.csv("/Users/jort/Documents/R/seasonal-vaccination-h5n1-serology/Data/20260819_mouse_experiment-body_weight_combined.csv")
bw_df_long <- bw_df %>%
  pivot_longer(
    !c(experiment, mouse, mouse_unique, group), names_to = "day", names_prefix = "X",
               names_transform = list(day = as.integer), values_to = "weight_pct"
  )
bw_df_long

summary_bw_df <- bw_df_long %>%
  group_by(group, day) %>%
  summarise(
    mean_val = mean(weight_pct, na.rm = TRUE),
    sem = sd(weight_pct, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
summary_bw_df

summary_bw_df$group <- factor(summary_bw_df$group, levels=levels(df$group))

plot <- ggplot(summary_bw_df, aes(x = day, y = mean_val, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean_val - sem, ymax = mean_val + sem), 
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +  # baseline
  scale_color_manual(values = c("#506E82", "#F59E0B", "#009E73", "#454545")) +
  scale_fill_manual(values = c("#506E82", "#F59E0B", "#009E73", "#454545")) +
  labs(x = "Days post-infection", y = "% Body weight", color = "Group", fill = "Group") +
  scale_x_continuous(
    limits = c(0, 14),
    breaks = c(0, 7, 14),
    minor_breaks = seq(0,14)
  ) +
  scale_y_continuous(
    limits = c(80, 105),
    breaks = c(80, 85, 90, 95, 100, 105)
  ) +
  guides(x = guide_axis(minor.ticks = TRUE)) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    strip.background = element_blank(),
    axis.line = element_line(linewidth = 0.35),
    legend.text = element_markdown()
  )

path <- "seasonal-vaccination-h5n1-serology/Figures/mouse_body-weight.pdf"
ggsave(path, width = 1.2, height = 1, dpi = 300, scale = 2, bg = "transparent")









# cs_df <- read.csv("/Users/jort/Documents/R/seasonal-vaccination-h5n1-serology/Data/20260714_mouse_experiment-clinical_scores.tsv", sep="\t")
cs_df <- read.csv("/Users/jort/Documents/R/seasonal-vaccination-h5n1-serology/Data/20260819_mouse_experiment-clinical_scores_combined.csv")
cs_df_long <- cs_df %>%
  pivot_longer(
    !c(experiment, mouse, mouse_unique, group), names_to = "day", names_prefix = "X",
    names_transform = list(day = as.integer), values_to = "clinical_score"
  )
cs_df_long

summary_cs_df <- cs_df_long %>%
  group_by(group, day) %>%
  summarise(
    mean_val = mean(clinical_score, na.rm = TRUE),
    sem = sd(clinical_score, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
summary_cs_df

summary_cs_df$group <- factor(summary_cs_df$group, levels=levels(df$group))

plot <- ggplot(summary_cs_df, aes(x = day, y = mean_val, color = group, fill = group)) +
  geom_ribbon(aes(ymin = mean_val - sem, ymax = mean_val + sem), 
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(values = c("#506E82", "#F59E0B", "#009E73", "#454545")) +
  scale_fill_manual(values = c("#506E82", "#F59E0B", "#009E73", "#454545")) +
  labs(x = "Days post-infection", y = "Clinical score", color = "Group", fill = "Group") +
  # scale_y_continuous(limits=c(0,10)) +
  scale_x_continuous(
    limits = c(0, 14),
    breaks = c(0, 7, 14),
    minor_breaks = seq(0,14)
  ) +
  guides(x = guide_axis(minor.ticks = TRUE)) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    strip.background = element_blank(),
    axis.line = element_line(linewidth = 0.35),
    legend.text = element_markdown()
  )

path <- "seasonal-vaccination-h5n1-serology/Figures/mouse_clinical-scores.pdf"
ggsave(path, width = 1.2, height = 1, dpi = 300, scale = 2, bg = "transparent")

