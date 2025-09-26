# Load required libraries
library(tidyverse)
library(readr)
library(ggplot2)
library(viridis)
library(gridExtra)
library(broom)
library(scales)
library(RColorBrewer)

# Define survey years and file paths
ebs_surveys <- c("200707", "200809", "200909", "201006", "201207", "201407", 
                 "201608", "201807", "202207", "202408")
ebs_result_files <- file.path(getwd(), "results", paste0("results_", ebs_surveys, ".csv"))
zdist_files <- file.path(getwd(), "results", paste0("zdists_", ebs_surveys, ".csv"))

# Read and process zdist files
zdists <- map_dfr(zdist_files, ~ read_csv(.x))
zdists <- zdists %>%
  mutate(zdist = str_replace(zdist, "Distributions\\.", ""))

zdists %>%
  group_by(zdist) %>%
  summarise(n = n(), .groups = 'drop')

zdists %>%
  group_by(class, zdist) %>%
  summarise(n = n(), .groups = 'drop')

# Read and process results files
results <- map_dfr(ebs_result_files, function(f) {
  survey <- str_sub(basename(f), 9, 14)
  cat(survey, "\n")
  
  df <- read_csv(f) %>%
    mutate(survey = survey,
           year = as.numeric(str_sub(survey, 1, 4))) %>%
    filter(age > 0) %>%
    pivot_longer(cols = c(n, biomass), names_to = "variable", values_to = "value")
  
  return(df)
})

# Calculate year-age CV
year_age <- results %>%
  filter(variable == "n") %>%
  group_by(year, age) %>%
  summarise(cv = sd(value) / mean(value), .groups = 'drop')

hist(year_age$cv)
quantile(year_age$cv, c(0.1, 0.9))

year_age_wide <- year_age %>%
  pivot_wider(names_from = year, values_from = cv, names_prefix = "year_")

# EVA-estimated CVs from cruise reports
eva <- data.frame(
  year = c(2007, 2008, 2009, 2010, 2012, 2014, 2016, 2018, 2022, 2024),
  cv_1d = c(3.8, 5.6, 6.9, 5.4, 3.4, 3.4, 1.9, 3.9, 6.8, 5.6),
  n = c(10.04, 5.24, 8.63, 12.97, 7.49, 19.51, 12.22, 5.57, 9.67, 11.55),
  biomass = c(2.28, 1.404, 1.331, 2.636, 2.279, 4.743, 4.838, 2.497, 3.834, 2.871)
)

# Calculate totals
totals <- results %>%
  group_by(survey, year, variable, i) %>%
  summarise(value = sum(value) / 1e9, .groups = 'drop')

# Create density plots
p1 <- totals %>%
  filter(variable == "n") %>%
  ggplot(aes(x = value, fill = factor(year))) +
  geom_density(alpha = 0.5) +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  labs(x = "Abundance (billions)", fill = "Year") +
  theme_minimal()

p2 <- totals %>%
  filter(variable == "biomass") %>%
  ggplot(aes(x = value, fill = factor(year))) +
  geom_density(alpha = 0.5) +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  labs(x = "Biomass (MT)", fill = "Year") +
  theme_minimal()

combined_density <- grid.arrange(p1, p2, ncol = 2)
ggsave(file.path(getwd(), "plots", "total_cv_distributions.png"), 
       combined_density, width = 10, height = 5, dpi = 300)

# QQ plots for abundance
create_qq_plots <- function(data, var_name, color_val) {
  years <- unique(data$year)
  plots <- map(years, function(yr) {
    data %>%
      filter(year == yr, variable == var_name) %>%
      ggplot(aes(sample = value)) +
      stat_qq(color = color_val) +
      stat_qq_line(color = "black") +
      ggtitle(yr) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  })
  return(plots)
}

qq_abundance <- create_qq_plots(totals, "n", "blue")
qq_abundance_combined <- do.call(grid.arrange, c(qq_abundance, ncol = 3))
ggsave(file.path(getwd(), "plots", "qq_normal_plots_abundance.png"), 
       qq_abundance_combined, width = 8, height = 8, dpi = 300)

qq_biomass <- create_qq_plots(totals, "biomass", "red")
qq_biomass_combined <- do.call(grid.arrange, c(qq_biomass, ncol = 3))
ggsave(file.path(getwd(), "plots", "qq_normal_plots_biomass.png"), 
       qq_biomass_combined, width = 8, height = 8, dpi = 300)

# Calculate annual statistics
annual <- totals %>%
  group_by(survey, year, variable) %>%
  summarise(
    upper = quantile(value, 0.975),
    lower = quantile(value, 0.025),
    std = sd(value),
    value = mean(value),
    .groups = 'drop'
  ) %>%
  mutate(cv = round(std / value * 100, 1)) %>%
  left_join(eva, by = "year") %>%
  mutate(cvstring = paste0(" ", cv, " (", cv_1d, ")")) %>%
  arrange(variable, year)

# Save annual uncertainty results
annual %>%
  select(year, variable, mean = value, std, cv, upper, lower) %>%
  write_csv(file.path(getwd(), "results", "annual_uncertainty.csv"))

# Summary statistics
annual %>%
  group_by(variable) %>%
  summarise(
    mean = mean(cv),
    min = min(cv),
    max = max(cv),
    .groups = 'drop'
  )

mean(eva$cv_1d)

# Create time series plots
p_n <- annual %>%
  filter(variable == "n") %>%
  ggplot(aes(x = year, y = value)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "blue") +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  geom_text(aes(label = str_extract(cvstring, "^[^(]+")), 
            hjust = 0, vjust = 0, size = 3) +
  scale_x_continuous(breaks = seq(2008, 2024, 2), limits = c(2006.5, 2026)) +
  scale_y_continuous(limits = c(0, 30)) +
  labs(x = "", y = "Abundance (billions)") +
  theme_minimal()

p_b <- annual %>%
  filter(variable == "biomass") %>%
  ggplot(aes(x = year, y = value)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "red") +
  geom_line(color = "red") +
  geom_point(color = "red") +
  geom_text(aes(label = cvstring), hjust = 0, vjust = 0, size = 3) +
  scale_x_continuous(breaks = seq(2008, 2024, 2), limits = c(2006.5, 2026)) +
  scale_y_continuous(limits = c(0, 6)) +
  labs(x = "Year", y = "Biomass (Mt)") +
  theme_minimal()

timeseries_plot <- grid.arrange(p_n, p_b, ncol = 1)
ggsave(file.path(getwd(), "plots", "timeseries.png"), 
       timeseries_plot, width = 8, height = 6, dpi = 300)

# Bootstrap vs 1D comparison
p_cv <- annual %>%
  ggplot(aes(x = year)) +
  geom_line(aes(y = cv, color = variable, linetype = "Bootstrap")) +
  geom_point(aes(y = cv, color = variable)) +
  geom_line(data = eva, aes(y = cv_1d, linetype = "1D geostatistical")) +
  geom_point(data = eva, aes(y = cv_1d)) +
  scale_color_manual(values = c("biomass" = "red", "n" = "blue"),
                     labels = c("Bootstrap (biomass)", "Bootstrap (abundance)")) +
  scale_linetype_manual(values = c("Bootstrap" = "solid", "1D geostatistical" = "dashed")) +
  scale_x_continuous(breaks = seq(2008, 2024, 2)) +
  scale_y_continuous(limits = c(0, 35)) +
  labs(x = "Year", y = "CV (%)", title = "(a)", color = "", linetype = "") +
  theme_minimal() +
  theme(plot.title.position = "plot")

# Linear regression for biomass
biomass_data <- annual %>% filter(variable == "biomass")
m_biomass <- lm(cv ~ cv_1d, data = biomass_data)
summary(m_biomass)

# Prediction for regression line
df_pred <- data.frame(cv_1d = seq(0, 8, 0.1))
pred_results <- predict(m_biomass, df_pred, interval = "confidence")
df_pred <- cbind(df_pred, pred_results)

coefs <- coef(m_biomass)

p_reg <- ggplot(df_pred, aes(x = cv_1d, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "red") +
  geom_line(color = "red") +
  geom_point(data = biomass_data, aes(x = cv_1d, y = cv), color = "red") +
  annotate("text", x = 5, y = 0, 
           label = paste0("y = ", round(coefs[1], 2), " + ", round(coefs[2], 2), "x"),
           color = "red", size = 4) +
  labs(x = "1D CV (%)", y = "Bootstrap CV (%)", title = "(b)") +
  theme_minimal() +
  theme(plot.title.position = "plot")

bootstrap_comparison <- grid.arrange(p_cv, p_reg, ncol = 2)
ggsave(file.path(getwd(), "plots", "bootstrap_vs_1d.png"), 
       bootstrap_comparison, width = 9, height = 3.5, dpi = 300)

# Calculate ratios
annual %>%
  group_by(variable) %>%
  summarise(ratio = mean(cv / cv_1d), .groups = 'drop')

# Annual age class analysis
annual_age <- results %>%
  filter(species_code == 21740) %>%
  group_by(year, variable, age) %>%
  summarise(
    mean = mean(value),
    std = sd(value),
    .groups = 'drop'
  ) %>%
  mutate(
    cv = std / mean,
    upper = mean + 2 * std,
    lower = mean - 2 * std
  )

write_csv(annual_age, file.path(getwd(), "results", "age_classes_uncertainty.csv"))

annual_age %>%
  group_by(variable) %>%
  summarise(
    mean = mean(cv),
    lower = quantile(cv, 0.05),
    upper = quantile(cv, 0.95),
    .groups = 'drop'
  )

# Age class bubble plots
pn_bubble <- annual_age %>%
  filter(variable == "n") %>%
  ggplot(aes(x = age, y = year)) +
  geom_point(aes(size = upper/5e8), stroke = 0.2, alpha = 0.7) +
  geom_point(aes(size = lower/5e8), color = "white", stroke = 0.2) +
  scale_x_continuous(breaks = 1:10, labels = c(as.character(1:9), "10+"),
                     limits = c(0, 11)) +
  scale_y_continuous(breaks = seq(2008, 2024, 2), limits = c(2006, 2025)) +
  labs(x = "Age class", y = "Year", title = "(a)") +
  theme_minimal() +
  theme(legend.position = "none", plot.title.position = "plot")

pb_bubble <- annual_age %>%
  filter(variable == "biomass") %>%
  ggplot(aes(x = age, y = year)) +
  geom_point(aes(size = upper/1e8), color = "red", stroke = 0.2, alpha = 0.7) +
  geom_point(aes(size = lower/1e8), color = "white", stroke = 0.2) +
  scale_x_continuous(breaks = 1:10, labels = c(as.character(1:9), "10+"),
                     limits = c(0, 11)) +
  scale_y_continuous(breaks = seq(2008, 2024, 2), limits = c(2006, 2025)) +
  labs(x = "Age class", y = "Year", title = "(b)") +
  theme_minimal() +
  theme(legend.position = "none", plot.title.position = "plot")

bubble_plots <- grid.arrange(pn_bubble, pb_bubble, ncol = 2)
ggsave(file.path(getwd(), "plots", "age_classes_bubbles.png"), 
       bubble_plots, width = 10, height = 5, dpi = 300)

# Error analysis
ebs_error_files <- file.path(getwd(), "results", paste0("stepwise_error_", ebs_surveys, ".csv"))

errors <- map_dfr(ebs_error_files, function(f) {
  survey <- str_sub(basename(f), 16, 21)
  cat(survey, "\n")
  
  df <- read_csv(f) %>%
    mutate(survey = survey,
           year = as.numeric(str_sub(survey, 1, 4))) %>%
    pivot_longer(cols = c(n, biomass), names_to = "variable", values_to = "value")
  
  return(df)
})

errors <- errors %>%
  mutate(variable = recode(variable, 
                          "n" = "Abundance", 
                          "biomass" = "Biomass"))

# Bootstrap function for error analysis
resample_df <- function(df, stochastic = TRUE) {
  n <- nrow(df)
  if (stochastic) {
    indices <- sample(1:n, n, replace = TRUE)
  } else {
    indices <- 1:n
  }
  return(df[indices, ])
}

# Bootstrap resampling (simplified - using 100 iterations for speed)
set.seed(123)
stds_boot <- map_dfr(1:100, function(i) {
  df <- resample_df(errors)
  df %>%
    group_by(year, error_label, variable) %>%
    summarise(cv = sd(value) / mean(value), .groups = 'drop') %>%
    mutate(iteration = i)
})

summary_boot <- stds_boot %>%
  group_by(year, error_label, variable) %>%
  summarise(
    cv = mean(cv),
    cv_se = sd(cv),
    .groups = 'drop'
  ) %>%
  mutate(error_label = factor(error_label, 
    levels = c("Calibration", "Spatial sampling", "Selectivity", "Resample catches",
               "Nearbottom coefs", "Trawl assignment", "TS models", "Age-length", 
               "Length-weight", "All"))) %>%
  arrange(year)

# Error series analysis
error_series <- errors %>%
  group_by(year, added_error, error_label, variable) %>%
  summarise(
    std = sd(value),
    mean = mean(value),
    .groups = 'drop'
  ) %>%
  mutate(cv = std / mean) %>%
  mutate(error_label = factor(error_label,
    levels = c("Calibration", "Spatial sampling", "Resample catches", "Selectivity",
               "Nearbottom coefs", "TS models", "Age-length", "Trawl assignment",
               "Length-weight", "All")))

# Error time series plots
pe1 <- error_series %>%
  filter(variable == "Abundance") %>%
  ggplot(aes(x = year, y = cv * 100, color = error_label)) +
  geom_line(size = 1) +
  geom_point() +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(2008, 2024, 2)) +
  labs(x = "", y = "Abundance CV (%)", color = "Error Source") +
  theme_minimal()

pe2 <- error_series %>%
  filter(variable == "Biomass") %>%
  ggplot(aes(x = year, y = cv * 100, color = error_label)) +
  geom_line(size = 1) +
  geom_point() +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = seq(2008, 2024, 2)) +
  labs(x = "Year", y = "Biomass CV (%)", color = "Error Source") +
  theme_minimal()

error_timeseries <- grid.arrange(pe1, pe2, ncol = 1)
ggsave(file.path(getwd(), "plots", "error_timeseries.png"), 
       error_timeseries, width = 8, height = 5, dpi = 300)

# Error summary statistics
error_summary <- error_series %>%
  group_by(error_label, variable) %>%
  summarise(
    cv_std = sd(cv),
    cv_max = quantile(cv, 0.75),
    cv_min = quantile(cv, 0.25),
    cv_se = sd(cv) / sqrt(n()),
    cv_mean = mean(cv),
    cv_med = median(cv),
    .groups = 'drop'
  ) %>%
  arrange(error_label)

# Boxplot of error sources
error_boxplot <- error_series %>%
  ggplot(aes(x = reorder(error_label, cv), y = cv * 100, fill = variable)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  scale_fill_manual(values = c("Abundance" = "blue", "Biomass" = "red")) +
  labs(x = "", y = "CV (%)", fill = "Variable") +
  theme_minimal() +
  theme(panel.grid.minor.y = TRUE)

ggsave(file.path(getwd(), "plots", "error_sources.png"), 
       error_boxplot, width = 5, height = 4, dpi = 300)

print("Analysis complete! Check the 'plots' directory for output files.")
