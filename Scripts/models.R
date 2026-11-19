# Load packages 
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn)
library(dplyr)
library(tibble)
library(corrplot)

extract_mixed_model_data <- function(mod, model_name) {
  # Exact profile likelihood CIs
  ci <- as.data.frame(confint(mod, method = "profile")) %>%
    rownames_to_column("term") %>%
    rename(
      conf.low = `2.5 %`,
      conf.high = `97.5 %`
    )
  
  # Fixed effects estimates
  est <- data.frame(
    term = names(fixef(mod)),
    estimate = fixef(mod)
  )
  
  # p-values from lmerTest
  pvals <- data.frame(
    term = rownames(coef(summary(mod))),
    p.value = coef(summary(mod))[, "Pr(>|t|)"]
  )
  
  # Combine
  model_table <- est %>%
    left_join(ci, by = "term") %>%
    left_join(pvals, by = "term") %>%
    mutate(
      p.value = case_when(
        p.value < 0.001 ~ formatC(p.value, format = "e", digits = 2),
        TRUE ~ sprintf("%.3f", p.value)
      ),
      model = model_name
    )
  
  return(model_table)
}

# Load results
wide_df <- read_csv("../Data/compiled_titers.csv")

mod <- lmer(
  log2(NT90) ~ log2(H5_stalk) +
    log2(NAI50) +
    log2(H1_stalk) +
    log2(H5_FL) +
    Visit +
    (1|Participant),
  data=wide_df,
  REML=TRUE
)

mod_dem <- lmer(
  log2(NT90) ~ log2(H5_stalk) +
    log2(NAI50) +
    log2(H1_stalk) +
    log2(H5_FL) +
    Visit +
    Sex + YOB + # demography
    (1|Participant),
  data=wide_df,
  REML=TRUE
)

r.squaredGLMM(mod)
r.squaredGLMM(mod_dem)
summary(mod)
summary(mod_dem)
anova(mod, mod_dem)


# Extract data from both models
mod_data <- extract_mixed_model_data(mod, "Unadjusted")
mod_dem_data <- extract_mixed_model_data(mod_dem, "Adjusted")

# Combine dataframes
combined_df <- bind_rows(mod_data, mod_dem_data)
write_csv(combined_df, "../Tables/linear_mixed_models.csv")

# Mapping for term names
mapping <- c(
  "log2(H5_stalk)" = "Headless H5 titer",
  "log2(NAI50)" = "NAI50 titer",
  "log2(H1_stalk)" = "Headless H1 titer",
  "log2(H5_FL)" = "Full-length H5 titer",
  "Visitv3" = "Visit (post-vaccination)",
  "SexM" = "Sex (male)",
  "YOB" = "Year of birth"
)

# Process combined data
plot_df <- combined_df %>%
  filter(term != "(Intercept)") %>%
  mutate(
    p.value = as.numeric(p.value),
    neg_log10_p = -log10(p.value),
    significant = p.value < 0.05,
    term_pretty = mapping[term]
  )

# Create ordered factor for y-axis (order by unadjusted model's significance)
term_order <- plot_df %>%
  filter(model == "Adjusted") %>%
  arrange(neg_log10_p) %>%
  pull(term_pretty) %>%
  unique()

plot_df <- plot_df %>%
  mutate(term_pretty = factor(term_pretty, levels = term_order)) %>%
  mutate(model = factor(model, levels = c("Adjusted", "Unadjusted")))

# Create plot with shapes differentiated by model
plot <- ggplot(plot_df, aes(x = estimate, y = term_pretty, shape = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, color = "black", position = position_dodge(width = 0.8)) +
  geom_point(aes(fill = neg_log10_p),
             size = 3.5, stroke = 1, color = "black",
             position = position_dodge(width = 0.8)) +
  scale_shape_manual(
    values = c("Unadjusted" = 21, "Adjusted" = 22),
    name = "Model",
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_fill_gradient(
    low = "white", high = "#506E82",
    limits = c(0, 5),
    name = expression(italic(P)),
    breaks = c(0, 1, 2, 3, 4, 5),
    labels = 10**(-c(0, 1, 2, 3, 4, 5))
  ) +
  labs(
    x = expression(beta ~ "(95% CI)"),
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.justification = "top",
    legend.title=element_blank()
  ) +
  guides(fill = "none") +
  scale_x_continuous(limits = c(-0.6, 0.8))

show(plot)
path <- "../Figures/ModelCoefficients-without_with_adjustment.pdf"
ggsave(path, width=3.25, height=1.75, dpi=300, scale=2, bg="transparent")

path <- "../Figures/ModelCoefficients-without_with_adjustment-customLegend.pdf"
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






# single HA titer models
mod_h1stalk <- lmer(
  log2(NT90) ~ log2(NAI50) +
    log2(H1_stalk) +
    Visit +
    Sex + YOB + # demography
    (1|Participant),
  data=wide_df,
  REML=TRUE
)

mod_h5stalk <- lmer(
  log2(NT90) ~ log2(NAI50) +
    log2(H5_stalk) +
    Visit +
    Sex + YOB + # demography
    (1|Participant),
  data=wide_df,
  REML=TRUE
)

mod_h5fl <- lmer(
  log2(NT90) ~ log2(NAI50) +
    log2(H5_FL) +
    Visit +
    Sex + YOB + # demography
    (1|Participant),
  data=wide_df,
  REML=TRUE
)


# Extract data from all models
mod_h1stalk_data <- extract_mixed_model_data(mod_h1stalk, "H1 headless only")
mod_h5stalk_data <- extract_mixed_model_data(mod_h5stalk, "H5 headless only")
mod_h5fl_dem_data <- extract_mixed_model_data(mod_h5fl, "H5 full-length only")

combined_singleHA_df <- bind_rows(mod_h1stalk_data, mod_h5stalk_data, mod_h5fl_dem_data)

write_csv(combined_singleHA_df, "../Tables/linear_mixed_models-single_HA_titer.csv")






# fold-change models
extract_model_data <- function(mod, model_name) {
  # Exact profile likelihood CIs
  ci <- as.data.frame(confint(mod, method = "profile")) %>%
    rownames_to_column("term") %>%
    rename(
      conf.low = `2.5 %`,
      conf.high = `97.5 %`
    )
  
  # Fixed effects estimates
  est <- data.frame(
    term = names(coef(mod)),
    estimate = coef(mod)
  )
  
  # p-values from lmerTest
  pvals <- data.frame(
    term = rownames(coef(summary(mod))),
    p.value = coef(summary(mod))[, "Pr(>|t|)"]
  )
  
  # Combine
  model_table <- est %>%
    left_join(ci, by = "term") %>%
    left_join(pvals, by = "term") %>%
    mutate(
      p.value = case_when(
        p.value < 0.001 ~ formatC(p.value, format = "e", digits = 2),
        TRUE ~ sprintf("%.3f", p.value)
      ),
      model = model_name
    )
  
  return(model_table)
}

df_fc <- wide_df %>%
  pivot_wider(
    id_cols = c(Participant, Sex, YOB),
    names_from = Visit,
    values_from = c(H5_FL, H5_stalk, H1_stalk, NAI50, NT90),
    names_sep = "_"
  ) %>%
  mutate(
    FC_H5_FL = H5_FL_v3 / H5_FL_v1,
    FC_H5_stalk = H5_stalk_v3 / H5_stalk_v1,
    FC_H1_stalk = H1_stalk_v3 / H1_stalk_v1,
    FC_NAI50 = NAI50_v3 / NAI50_v1,
    FC_NT90 = NT90_v3 / NT90_v1
  )

fc_mod <- lm(
  log2(FC_NT90) ~ log2(FC_H5_stalk) +
    log2(FC_NAI50) +
    log2(FC_H1_stalk) +
    log2(FC_H5_FL),
  data=df_fc
)

fc_mod_dem <- lm(
  log2(FC_NT90) ~ log2(FC_H5_stalk) +
    log2(FC_NAI50) +
    log2(FC_H1_stalk) +
    log2(FC_H5_FL) +
    YOB + Sex, # demography
  data=df_fc
)

r.squaredGLMM(fc_mod)
r.squaredGLMM(fc_mod_dem)
summary(fc_mod)
summary(fc_mod_dem)
anova(fc_mod, fc_mod_dem)


# Extract data from both models
fc_mod_data <- extract_model_data(fc_mod, "Unadjusted")
fc_mod_dem_data <- extract_model_data(fc_mod_dem, "Adjusted")

# Combine dataframes
fc_combined_df <- bind_rows(fc_mod_data, fc_mod_dem_data)
write_csv(fc_combined_df, "../Tables/fc_linear_models.csv")

# Mapping for term names
fc_mapping <- c(
  "log2(FC_H5_stalk)" = "Headless H5 titer fold-change",
  "log2(FC_NAI50)" = "NAI50 titer fold-change",
  "log2(FC_H1_stalk)" = "Headless H1 titer fold-change",
  "log2(FC_H5_FL)" = "Full-length H5 titer fold-change",
  "Visitv3" = "Visit (post-vaccination)",
  "SexM" = "Sex (male)",
  "YOB" = "Year of birth"
)

# Process combined data
fc_plot_df <- fc_combined_df %>%
  filter(term != "(Intercept)") %>%
  mutate(
    p.value = as.numeric(p.value),
    neg_log10_p = -log10(p.value),
    significant = p.value < 0.05,
    term_pretty = fc_mapping[term]
  )

# Create ordered factor for y-axis (order by unadjusted model's significance)
fc_term_order <- fc_plot_df %>%
  filter(model == "Adjusted") %>%
  arrange(neg_log10_p) %>%
  pull(term_pretty) %>%
  unique()

fc_plot_df <- fc_plot_df %>%
  mutate(term_pretty = factor(term_pretty, levels = fc_term_order)) %>%
  mutate(model = factor(model, levels = c("Adjusted", "Unadjusted")))

# Create plot with shapes differentiated by model
plot <- ggplot(fc_plot_df, aes(x = estimate, y = term_pretty, shape = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, color = "black", position = position_dodge(width = 0.8)) +
  geom_point(aes(fill = neg_log10_p),
             size = 3.5, stroke = 1, color = "black",
             position = position_dodge(width = 0.8)) +
  scale_shape_manual(
    values = c("Unadjusted" = 21, "Adjusted" = 22),
    name = "Model",
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_fill_gradient(
    low = "white", high = "#506E82",
    limits = c(0, 3),
    name = expression(italic(P)),
    breaks = c(0, 1, 2, 3),
    labels = 10**(-c(0, 1, 2, 3))
  ) +
  labs(
    x = expression(beta ~ "(95% CI)"),
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.justification = "top",
    legend.title=element_blank()
  ) +
  guides(fill = "none") +
  scale_x_continuous(limits = c(-0.6, 0.8))

show(plot)
path <- "../Figures/ModelCoefficients-FoldChange-without_with_adjustment.pdf"
ggsave(path, width=3.25, height=1.75, dpi=300, scale=2, bg="transparent")

path <- "../Figures/ModelCoefficients-FoldChange-without_with_adjustment-customLegend.pdf"
pdf(file = path, width=2, height=4)
plot.new()

colorlegend(
  colbar = colorRampPalette(c("#FFFFFF", "#506E82"))(200),
  labels = 10**(-c(0,1,2,3)),
  vertical = TRUE,
  ratio.colbar = 0.25,
  xlim = c(0, 1),
  ylim = c(0, 1),
  align = "l",
  cex = 0.8
)

dev.off()
