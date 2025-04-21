# SDL "Regression Rebels"
# Research Question: Is air quality associated with prevalence of asthma?

# 02_statistical_analysis.R ----------------------------------------------
library(tidyverse)
library(gt)
library(dplyr)

# Load cleaned data
filtered_data_2022 <- read_csv("C:/Users/Owner/Downloads/cleaned_2022_data.csv")

# Calculate IQR bounds
Q1 <- quantile(filtered_data$rate_per_1000, 0.25, na.rm = TRUE)
Q3 <- quantile(filtered_data$rate_per_1000, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# Model 1: Poisson with population offset
poisson_model <- glm(
  total_prescriptions_2022 ~ avg_aqi_2022,
  data = filtered_data_2022,
  family = poisson(link = "log"),
  offset = log(POPESTIMATE2022)
)

summary(poisson_model)

# Add predictions
filtered_data_2022 <- filtered_data_2022 |>
  mutate(
    fitted = predict(poisson_model, type = "response"),
    observed_rate_per_1000 = rate_per_1000_2022,
    fitted_rate_per_1000 = (fitted / POPESTIMATE2022) * 1000
  )

# Plot Model 1
ggplot(filtered_data_2022, aes(x = avg_aqi_2022)) +
  geom_point(aes(y = observed_rate_per_1000), color = "blue") +
  geom_line(aes(y = fitted_rate_per_1000), color = "red") +
  labs(
    title = "Poisson Regression",
    subtitle = "Asthma Prescriptions per 1,000 vs. AQI",
    x = "Average AQI",
    y = "Prescriptions per 1,000 Residents"
  ) +
  theme_minimal()

# Model 2: Poisson Regression offset by all prescriptions

poisson_model_offset_all_prescriptions <- glm(
  total_prescriptions_2022 ~ avg_aqi_2022,
  data = filtered_data_2022,
  family = poisson(link = "log"),
  offset = log(all_prescriptions_2022)
)

summary(poisson_model_offset_all_prescriptions)

filtered_data_2022 <- filtered_data_2022 |>
  mutate(
    fitted = predict(poisson_model_offset_all_prescriptions, type = "response"),
    predicted_fraction = fitted / all_prescriptions_2022
  )

ggplot(filtered_data_2022, aes(x = avg_aqi_2022)) +
  geom_point(aes(y = asthma_fraction_2022), color = "blue") +
  geom_line(aes(y = predicted_fraction), color = "red") +
  labs(
    title = "Asthma Prescription Fraction vs. AQI",
    subtitle = "Poisson with Offset for All Prescriptions",
    x = "Average AQI",
    y = "Asthma Fraction"
  ) +
  theme_minimal()

# Log-Linear Model for Asthma Cost Fraction vs AQI -----------------------

# Log-transform the outcome
filtered_data_2022 <- filtered_data_2022 |>
  mutate(log_cost_fraction = log(asthma_cost_fraction_2022))

# Fit linear regression model
log_cost_model <- lm(log_cost_fraction ~ avg_aqi_2022, data = filtered_data_2022)

# Model summary
summary(log_cost_model)

# Predict and back-transform
filtered_data_2022 <- filtered_data_2022 |>
  mutate(
    predicted_log = predict(log_cost_model),
    predicted_fraction = exp(predicted_log)
  )

# Plotting
ggplot(filtered_data_2022, aes(x = avg_aqi_2022)) +
  geom_point(aes(y = asthma_cost_fraction_2022), color = "blue", alpha = 0.6) +
  geom_line(aes(y = predicted_fraction), color = "red", size = 1.2) +
  labs(
    title = "Asthma Drug Cost Fraction vs. AQI",
    subtitle = "Log-Transformed Linear Regression",
    x = "Average AQI (2022)",
    y = "Fraction of Total Drug Cost for Asthma"
  ) +
  theme_minimal()


# Temporal Analysis

# Load matched city-year data
matched_data <- read_csv("C:/Users/Owner/Downloads/matched_data_2019_2022.csv")

# Paired t-test
t_test_result <- t.test(
  matched_data$asthma_cost_fraction_2022,
  matched_data$asthma_cost_fraction_2019,
  paired = TRUE
)
print(t_test_result)

# Linear regression on changes
change_model <- lm(
  change_in_cost_fraction ~ change_in_aqi,
  data = matched_data
)

summary(change_model)

# Plot results
ggplot(matched_data, aes(x = change_in_aqi, y = change_in_cost_fraction)) +
  geom_point(color = "darkblue", alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Change in Asthma Cost Fraction vs. Change in AQI",
    subtitle = "2019 to 2022 by Matched City",
    x = "Change in Average AQI",
    y = "Change in Asthma Drug Cost Fraction"
  ) +
  theme_minimal()

# Summary Table of Statistical Tests ---------------------------------------

# Extract model summaries
poisson_summary <- summary(poisson_model)
poisson_offset_summary <- summary(poisson_model_offset_all_prescriptions)
log_cost_summary <- summary(log_cost_model)
change_model_summary <- summary(change_model)

# Create a summary table with 95% confidence intervals
stat_summary_table <- tibble::tibble(
  model = c(
    "Poisson: Total Asthma Prescriptions ~ Average AQI",
    "Poisson: Asthma Prescription Fraction ~ Average AQI",
    "Log-Linear: Cost Fraction ~ Average AQI",
    "Linear: ΔCost Fraction ~ ΔAQI (2019 ~ 2022)"
  ),
  estimate = c(
    poisson_summary$coefficients[2, 1],
    poisson_offset_summary$coefficients[2, 1],
    log_cost_summary$coefficients[2, 1],
    change_model_summary$coefficients[2, 1]
  ),
  std_error = c(
    poisson_summary$coefficients[2, 2],
    poisson_offset_summary$coefficients[2, 2],
    log_cost_summary$coefficients[2, 2],
    change_model_summary$coefficients[2, 2]
  ),
  p_value = c(
    poisson_summary$coefficients[2, 4],
    poisson_offset_summary$coefficients[2, 4],
    log_cost_summary$coefficients[2, 4],
    change_model_summary$coefficients[2, 4]
  )
) |>
  mutate(
    ci_lower = estimate - 1.96 * std_error,
    ci_upper = estimate + 1.96 * std_error,
    conf_int = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]")
  )

# View table
print(stat_summary_table)

# Format table for research letter
stat_summary_table |>
  mutate(
    estimate  = round(estimate, 3),
    p_value   = signif(p_value, 3)
  ) |>
  select(model, estimate, conf_int, p_value) |>
  gt() |>
  tab_header(
    title = "Summary of Statistical Test Results",
    subtitle = "Asthma Outcomes and Air Quality (2022 vs 2019)"
  ) |>
  cols_label(
    model     = "Model/Test",
    estimate  = "Estimate",
    conf_int  = "95% CI",
    p_value   = "p-value"
  ) |>
  tab_options(
    table.font.size = px(13),
    heading.title.font.size = px(16),
    heading.subtitle.font.size = px(14),
    column_labels.font.weight = "bold",
    table.border.top.width = px(1),
    table.border.bottom.width = px(1),
    table.border.bottom.color = "gray"
  )

stat_summary_table |>
  mutate(
    estimate  = round(estimate, 3),
    p_value   = ifelse(p_value < 2.2e-16, "< 2.2e-16", signif(p_value, 3))
  ) |>
  select(model, estimate, conf_int, p_value) |>
  gt() |>
  tab_header(
    title = "Summary of Statistical Test Results",
    subtitle = "Asthma Outcomes and Air Quality (2022 vs 2019)"
  ) |>
  cols_label(
    model     = "Model/Test",
    estimate  = "Estimate",
    conf_int  = "95% CI",
    p_value   = "p-value"
  )
