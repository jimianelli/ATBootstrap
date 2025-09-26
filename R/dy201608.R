library(readr)
library(dplyr)
library(ggplot2)
library(here)

# Source the ATBootstrap functions (assuming they've been ported to R)
source(file.path(here(), "R", "ATBootstrap.R"))

# Set survey and seed
survey <- "201608"
set.seed(as.numeric(survey))

# Define paths
surveydir <- file.path(here(), "surveydata", survey)

# Constants
km2nmi <- 1 / 1.852
resolution <- 10.0  # km
dA <- (resolution * km2nmi)^2

# Log ranges
log_ranges <- list(c(279.5, 8427.5), c(8430.5, 14259))

# Preprocess survey data
preprocess_survey_data(surveydir, dx = resolution, ebs = TRUE, log_ranges = log_ranges)

# Read survey files
survey_files <- read_survey_files(surveydir)
acoustics <- survey_files$acoustics
scaling <- survey_files$scaling
age_length <- survey_files$age_length
length_weight <- survey_files$length_weight
trawl_locations <- survey_files$trawl_locations
surveygrid <- survey_files$surveygrid

# Create scatter plot with acoustics data
p1 <- ggplot(acoustics, aes(x = x, y = y, color = class, size = nasc/500)) +
  geom_point(alpha = 0.5, stroke = 0) +
  coord_equal() +
  theme_minimal()

# Add trawl locations
p1 <- p1 + geom_point(data = trawl_locations, aes(x = x, y = y), 
                      color = "black", size = 2, inherit.aes = FALSE)

print(p1)

# Create survey data object
surveydata <- ATSurveyData(acoustics, scaling, age_length, length_weight, 
                          trawl_locations, surveygrid, dA)

# Create bootstrap problem
atbp <- ATBootstrapProblem(surveydata)

# Calculate simulation distances
sim_dists <- zdists(atbp)
sim_dists$survey <- survey

# Write results
write_csv(sim_dists %>% select(survey, class, zdist),
          file.path(here(), "results", paste0("zdists_", survey, ".csv")))

# Plot geosimulation statistics
plot_geosim_stats(atbp, surveydata, 500)
ggsave(file.path(here(), "plots", paste0("conditional_nasc_stats_", survey, ".png")))

# Plot variograms
plot_class_variograms(atbp)

# Plot simulated NASC
plot_simulated_nasc(atbp, surveydata, width = 1000, height = 600)

# Get first class problem and simulate
cp1 <- atbp$class_problems[[1]]
nasc <- simulate_nasc(cp1)
dom <- solution_domain(cp1)

# Create scatter plot of simulation
p2 <- ggplot(data.frame(x = dom$x, y = dom$y, nasc = nasc), 
             aes(x = x, y = y, color = nasc)) +
  geom_point(stroke = 0) +
  scale_color_viridis_c() +
  labs(title = "2016") +
  theme_minimal() +
  theme(legend.position = "none")

# Add acoustics points
p2 <- p2 + geom_point(data = acoustics, aes(x = x, y = y), 
                      color = "white", size = 2, inherit.aes = FALSE)

print(p2)

# Run bootstrap simulation
results <- simulate(atbp, surveydata, nreplicates = 500)

# Plot bootstrap results
plot_boot_results(results)

# Write results
write_csv(results, file.path(here(), "results", paste0("results_", survey, ".csv")))

# Summarize bootstrap results
n_summary <- summarize_bootstrap(results, "n")
biomass_summary <- summarize_bootstrap(results, "biomass")

# Stepwise error analysis
results_step <- stepwise_error(atbp, surveydata, nreplicates = 500)

# Plot error sources by age
plot_error_source_by_age(results_step, results, "n")

# Merge results
results_totals <- merge_results(results, results_step)

# Write stepwise error results
write_csv(results_totals, 
          file.path(here(), "results", paste0("stepwise_error_", survey, ".csv")))

# Plot error sources
plot_error_sources(results_totals, 
                   plot_title = survey,
                   xticks = seq(0, 0.15, 0.01),
                   width = 800, 
                   height = 500)
