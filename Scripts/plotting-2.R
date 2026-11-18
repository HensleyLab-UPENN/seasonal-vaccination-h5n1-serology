# don't use

library(tidyverse)
library(minpack.lm)

df <- read_csv("/Users/jort/coding/GFP-neut-area-analysis/20260714_depletion_neut_data.csv")

# Convert to long format
plot_df <- df %>%
  pivot_longer(
    cols = c(mock, stalk, ha, na, `ha+na`),
    names_to = "condition",
    values_to = "value"
  ) %>%
  mutate(
    log2_dilution = log2(dilution)
  )

# Optional: remove the highest dilution from fitting only
fit_df <- plot_df %>%
  group_by(participant, condition) %>%
  filter(dilution != max(dilution)) %>%
  ungroup()

# Function to fit one 4PL
fit_4pl <- function(dat) {
  
  fit <- try(
    nlsLM(
      value ~ bottom + (top - bottom) /
        (1 + exp((ec50 - log2_dilution) * hill)),
      data = dat,
      start = list(
        bottom = min(dat$value),
        top = max(dat$value),
        ec50 = median(dat$log2_dilution),
        hill = 1
      ),
      control = nls.lm.control(maxiter = 500)
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(tibble())
  }
  
  xpred <- seq(
    min(dat$log2_dilution),
    max(dat$log2_dilution),
    length.out = 200
  )
  
  tibble(
    log2_dilution = xpred,
    dilution = 2^xpred,
    value = predict(
      fit,
      newdata = data.frame(log2_dilution = xpred)
    )
  )
}

# Fit curves
curve_df <- bind_rows(
  fit_df %>%
    group_by(participant, condition) %>%
    group_map(~{
      pred <- fit_4pl(.x)
      
      pred$participant <- .y$participant
      pred$condition <- .y$condition
      
      pred
    })
)

# Plot
ggplot(plot_df,
       aes(dilution, value, color = condition)) +
  geom_point(size = 2, alpha=0.2) +
  geom_line(
    data = curve_df,
    linewidth = 1
  ) +
  scale_x_continuous(
    trans = "log2",
    breaks = scales::log_breaks(base = 2)
  ) +
  facet_wrap(~participant) +
  labs(
    x = "Dilution",
    y = "% Infection",
    color = NULL
  ) +
  theme_bw()





library(purrr)

# Fit curves
curve_df <- bind_rows(
  fit_df %>%
    group_by(participant, condition) %>%
    group_map(~{
      pred <- fit_4pl(.x)
      
      pred$participant <- .y$participant
      pred$condition <- .y$condition
      
      pred
    })
)

# Plot
ggplot(plot_df,
       aes(dilution, value, color = condition)) +
  geom_point(size = 2, alpha=0.2) +
  geom_line(
    data = curve_df,
    linewidth = 1
  ) +
  scale_x_continuous(
    trans = "log2",
    breaks = scales::log_breaks(base = 2)
  ) +
  facet_wrap(~participant) +
  labs(
    x = "Dilution",
    y = "% Infection",
    color = NULL
  ) +
  theme_bw()

participant_plots <- map(
  unique(plot_df$participant),
  function(id) {
    
    ggplot(
      filter(plot_df, participant == id),
      aes(dilution, value, color = condition)
    ) +
      geom_point(size = 2, alpha = 0.4) +
      geom_line(
        data = filter(curve_df, participant == id),
        linewidth = 1.3
      ) +
      scale_x_continuous(
        trans = "log2",
        breaks = c(32, 64, 128, 256, 512, 1024, 2048, 4096)
      ) +
      coord_cartesian(ylim = c(0, 110)) +
      labs(
        title = id,
        x = "Dilution",
        y = "% Infection",
        color = NULL
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
  }
)

names(participant_plots) <- unique(plot_df$participant)









# with 50% inhibitory concentrations
library(dplyr)
library(purrr)
library(tidyverse)
library(minpack.lm)

df <- read_csv("/Users/jort/coding/GFP-neut-area-analysis/20260714_depletion_neut_data.csv")

# Convert to long format
plot_df <- df %>%
  pivot_longer(
    cols = c(mock, stalk, ha, na, `ha+na`),
    names_to = "condition",
    values_to = "value"
  ) %>%
  mutate(
    log2_dilution = log2(dilution)
  )

# # Function to fit one 4PL
# fit_4pl <- function(dat) {
#   
#   fit <- try(
#     nlsLM(
#       value ~ bottom + (top - bottom) /
#         (1 + exp((ec50 - log2_dilution) * hill)),
#       data = dat,
#       start = list(
#         bottom = min(dat$value),
#         top = max(dat$value),
#         ec50 = median(dat$log2_dilution),
#         hill = 1
#       ),
#       control = nls.lm.control(maxiter = 500)
#     ),
#     silent = TRUE
#   )
#   
#   if (inherits(fit, "try-error")) {
#     message(sprintf("Fit failed for %s / %s: %s",
#                     unique(dat$participant), unique(dat$condition),
#                     attr(fit, "condition")$message))
#     return(tibble())
#   }
#   
#   xpred <- seq(
#     min(dat$log2_dilution),
#     max(dat$log2_dilution),
#     length.out = 200
#   )
#   
#   tibble(
#     log2_dilution = xpred,
#     dilution = 2^xpred,
#     value = predict(
#       fit,
#       newdata = data.frame(log2_dilution = xpred)
#     )
#   )
# }

# Function to fit one 4PL, with hill slope bounded at 10
fit_4pl <- function(dat) {
  
  fit <- try(
    nlsLM(
      value ~ bottom + (top - bottom) / (1 + exp((ec50 - log2_dilution) * hill)),
      data = dat,
      start = list(
        bottom = min(dat$value, na.rm = TRUE),
        top = max(dat$value, na.rm = TRUE),
        ec50 = median(dat$log2_dilution, na.rm = TRUE),
        hill = 1
      ),
      lower = c(bottom = -Inf, top = -Inf, ec50 = -Inf, hill = 0.002),
      upper = c(bottom = Inf, top = Inf, ec50 = Inf, hill = ),
      control = nls.lm.control(maxiter = 500)
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    message(sprintf(
      "Fit failed for %s / %s: %s",
      unique(dat$participant), unique(dat$condition),
      attr(fit, "condition")$message
    ))
    return(tibble())
  }

  xpred <- seq(
    min(dat$log2_dilution),
    max(dat$log2_dilution),
    length.out = 200
  )

  tibble(
    log2_dilution = xpred,
    dilution = 2^xpred,
    value = predict(
      fit,
      newdata = data.frame(log2_dilution = xpred)
    )
  )
}

# Fit curves
curve_df <- bind_rows(
  plot_df %>%
    group_by(participant, condition) %>%
    group_map(~{
      pred <- fit_4pl(.x)
      
      pred$participant <- .y$participant
      pred$condition <- .y$condition
      
      pred
    })
)

ic50_df <- curve_df %>%
  group_by(participant, condition) %>%
  group_modify(~{
    
    pred <- .x %>% arrange(value)
    
    # Skip curves that never cross 50%
    if (min(pred$value) > 50 || max(pred$value) < 50)
      return(tibble())
    
    tibble(
      dilution = approx(
        x = pred$value,
        y = pred$dilution,
        xout = 50
      )$y,
      value = 50
    )
  }) %>%
  ungroup()

participant_plots <- map(
  unique(plot_df$participant),
  function(id) {
    ggplot(
      filter(plot_df, participant == id),
      aes(dilution, value, color = condition)
    ) +
      geom_point(size = 2, alpha = 0.4) +
      geom_line(
        data = filter(curve_df, participant == id),
        linewidth = 1.3
      ) +
      
      geom_hline(
        yintercept = 50,
        linetype = 2,
        color = "grey70"
      ) +
      geom_point(
        data = filter(ic50_df, participant == id),
        aes(dilution, value, color = condition),
        size = 3
      ) +
      geom_text(
        data = filter(ic50_df, participant == id),
        aes(
          dilution,
          value,
          label = sprintf("%.0f", dilution),
          color = condition
        ),
        nudge_y = 6,
        size = 3,
        show.legend = FALSE
      ) +
      geom_vline(
        data = filter(ic50_df, participant == id),
        aes(xintercept = dilution, color = condition),
        linetype = 3,
        linewidth = 0.6,
        show.legend = FALSE
      ) +
      
      scale_x_continuous(
        trans = "log2",
        breaks = c(32, 64, 128, 256, 512, 1024, 2048, 4096)
      ) +
      coord_cartesian(ylim = c(0, 110)) +
      labs(
        title = id,
        x = "Dilution",
        y = "% Infection",
        color = NULL
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
  }
)

names(participant_plots) <- unique(plot_df$participant)

walk(participant_plots, print)

