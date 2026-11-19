# Load packages 
library(tidyverse)
library(ggplot2)
library(viridis)
library(corrplot)
library(grid)
library(cowplot)
library(ggpubr)
library(ggnewscale)
library(drc)

fit_4pl <- function(x, y) {
  # remove NAs if necessary
  data_clean <- data.frame(x = x, y = y) %>% filter(!is.na(x) & !is.na(y))
  
  if(nrow(data_clean) < 4) return(NULL)  # Need at least 4 points
  
  tryCatch({
    # fit using drc package (LL.4 = 4-parameter logistic)
    # use floor of 0 and ceiling of 100
    model <- drm(y ~ x, data = data_clean, fct = LL.4(fixed = c(NA, 0, 100, NA), names = c("b", "c", "d", "e")))
    return(model)
  }, error = function(e) {
    return(NULL)
  })
}

data_df <- read.csv("../Data/20260714_depletion_neut_data.csv")
data_df$participant <- str_pad(data_df$participant, width = 4, side = "left", pad = "0")
participants <- unique(data_df$participant)

# convert to long format
plot_df <- data_df %>%
  pivot_longer(
    cols = c(mock, stalk, ha, na, `ha.na`),
    names_to = "condition",
    values_to = "value"
  ) %>%
  mutate(
    log2_dilution = log2(dilution)
  )

condition_mapping = c(
  "mock" = "Mock",
  "stalk" = "Headless H5",
  "ha" = "Full-length H5",
  "na" = "N1 NA",
  "ha.na" = "Full-length H5 + N1 NA"
)


fits <- plot_df %>%
  nest(data = -c(participant, condition)) %>%
  mutate(
    model = map(data, ~fit_4pl(.x$log2_dilution, .x$value))
  ) %>%
  filter(!map_lgl(model, is.null))


predictions <- fits %>%
  mutate(
    pred_data = map2(model, data, ~{
      x_range <- seq(min(.y$log2_dilution), max(.y$log2_dilution), length.out = 1000)
      pred_y <- predict(.x, data.frame(x = x_range))
      data.frame(log2_dilution = x_range, dilution = 2**x_range, value = pred_y)
    })
  ) %>%
  dplyr::select(participant, condition, pred_data) %>%
  unnest(pred_data)

ic50_df <- predictions %>%
  group_by(participant, condition) %>%
  group_modify(~{
    
    pred <- .x %>% arrange(value)
    
    # Skip curves that never cross 50%
    if (max(pred$value) < 50)
      return(tibble())
    
    if (min(pred$value) > 50)
      return(tibble(
        dilution = 10,
        value = 50
      ))
    
    tibble(
      dilution = approx(
        x = pred$value,
        y = pred$dilution,
        xout = 50
      )$y,
      value = 50
    )
  }) %>%
  ungroup() %>%
  complete(
    participant,
    condition,
    fill = list(dilution = 10, value = 50)
  )



plot_df$condition <- factor(
  plot_df$condition,
  levels = c("mock", "stalk", "ha", "na", "ha.na")
)

predictions$condition <- factor(
  predictions$condition,
  levels = levels(plot_df$condition)
)

ic50_df$condition <- factor(
  ic50_df$condition,
  levels = levels(plot_df$condition)
)

plot_df <- plot_df %>%
  arrange(condition)
predictions <- predictions %>%
  arrange(condition)
ic50_df <- ic50_df %>%
  arrange(condition)

for (p in unique(plot_df$participant)) {
  pred_data <- predictions %>% filter(participant == p)
  # pred_data$dilution = 2**pred_data$log2_dilution
  # pred_data$value = pred_data$pred_value
  plot <- ggplot(
    filter(plot_df, participant == p),
    aes(
      dilution,
      value,
      color = condition,
      fill = condition,
      shape = condition,
      group = condition
    )
  ) +
    geom_line(
      data = filter(pred_data, participant == p),
      linewidth = 1
    ) +
    geom_point(
      size = 2,
      alpha = 0.4,
      stroke = 0
    ) +
    geom_hline(
      yintercept = 50,
      linetype = 2,
      color = "grey70"
    ) +
    geom_point(
      data = filter(ic50_df, participant == p),
      aes(dilution, value),
      size = 2,
      stroke = 1,
      alpha = 1,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        mock  = "#506E82",
        stalk = "#D97706",
        ha    = "#F59E0B",
        na    = "#009E73",
        ha.na = "#454545"
      ),
      name = "Depletion condition",
      labels = condition_mapping[levels(plot_df$condition)]
    ) +
    scale_fill_manual(
      values = c(
        mock  = "#506E82",
        stalk = "#D97706",
        ha    = "#F59E0B",
        na    = "#009E73",
        ha.na = "#454545"
      ),
      guide = "none"
    ) +
    scale_shape_manual(
      values = c(
        mock  = 21,
        stalk = 22,
        ha    = 23,
        na    = 24,
        ha.na = 25
      ),
      name = "Condition",
      labels = condition_mapping
    ) +
    scale_x_continuous(
      trans = "log2",
      limits = c(20, 640*2**0.5),
      breaks = c(20, 40, 80, 160, 320, 640)
      # labels = c(20, 40, 80, 160, 320, 640)
    ) +
    scale_y_continuous(
      limits = c(0, 125),
      breaks = c(0, 25, 50, 75, 100, 125)
    ) +
    coord_cartesian(ylim = c(0, 125)) +
    labs(
      title = p,
      x = "Plasma dilution",
      y = "% Infection"
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = c(21, 22, 23, 24, 25),
          color = c(
            "#506E82",
            "#D97706",
            "#F59E0B",
            "#009E73",
            "#454545"
          ),
          fill = c(
            "#506E82",
            "#D97706",
            "#F59E0B",
            "#009E73",
            "#454545"
          ),
          alpha = 1,
          size = 3,
          linewidth = 1.3
        )
      ),
      shape = "none"
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
      
    )
  
  plot_no_legend <- plot + theme(legend.position = "none")
  path <- paste0("../Figures/",p,"-depletion_neut_curves-noLegend.pdf")
  ggsave(path, width = 1.2, height = 1, dpi = 300, scale = 2, bg = "transparent")
  
}

legend <- get_legend(plot)
path <- "../Figures/depletion_neut_curves-separateLegend.pdf"
pdf(file = path, width = 2.5, height = 2)
grid.newpage()
grid.draw(legend)
dev.off()




ic50_df$condition <- factor(
  ic50_df$condition,
  levels = c("mock", "stalk", "ha", "na", "ha.na")
)

condition_mapping_axis = c(
  "mock" = "Mock",
  "stalk" = "Headless\nH5 HA",
  "ha" = "Full-length\nH5 HA",
  "na" = "N1 NA",
  "ha.na" = "Full-length\nH5 HA +\nN1 NA"
)


plot <- ggplot(
  ic50_df,
  aes(
    x = condition,
    y = dilution,
    group = participant,
    shape = condition
  )
) +
  geom_line(color = "grey75", alpha = 0.5) +
  geom_point(aes(fill = condition, color = condition), size = 2, stroke=1) +
  geom_hline(
    yintercept = 20,
    linetype = 2,
    color = "grey70"
  ) +
  scale_color_manual(
    values = c(
      mock  = "#506E82",
      stalk = "#D97706",
      ha    = "#F59E0B",
      na    = "#009E73",
      ha.na = "#454545"
    ),
    name = "Depletion",
    labels = condition_mapping
  ) +
  scale_x_discrete(labels = condition_mapping_axis) +
  scale_fill_manual(
    values = c(
      mock  = "#506E82",
      stalk = "#D97706",
      ha    = "#F59E0B",
      na    = "#009E73",
      ha.na = "#454545"
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      mock  = 21,
      stalk = 22,
      ha    = 23,
      na    = 24,
      ha.na = 25
    ),
    name = "Condition",
    labels = condition_mapping
  ) +
  scale_y_continuous(
    trans="log2",
    limits = c(10/2**0.5, 320*2**0.5),
    breaks = c(10, 20, 40, 80, 160, 320)
  ) +
  labs(
    x = "Depletion condition",
    y = "ID50",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold",
                              margin = margin(b = -10)),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background  = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    strip.background = element_blank(),
    axis.line = element_line(linewidth = 0.35)
  )

path <- "../Figures/depletion_neut_IC50s.pdf"
ggsave(path, width = 2.5, height = 1, dpi = 300, scale = 2, bg = "transparent")


