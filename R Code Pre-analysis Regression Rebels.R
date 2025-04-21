# SDL "Regression Rebels"
# Research Question: Is air quality associated with prevalence of asthma?

# ----------------------------------------------------------------------
# Load libraries -------------------------------------------------------
# ----------------------------------------------------------------------

library(tidyverse)
library(haven)
library(stringdist)
library(fuzzyjoin)

# ----------------------------------------------------------------------
# Shared objects & helper functions ------------------------------------
# ----------------------------------------------------------------------

# Asthma-related drugs
target_drugs <- c(
  "Albuterol", "Salmeterol", "Formoterol", "Ipratropium", "Tiotropium",
  "Fluticasone", "Budesonide", "Mometasone", "Beclomethasone", "Montelukast",
  "Zafirlukast", "Theophylline", "Omalizumab", "Mepolizumab", "Reslizumab",
  "Dupilumab", "Roflumilast", "Acetylcysteine", "Dornase alfa", "Guaifenesin"
)

# Clean city names
clean_city_name <- function(city) {
  city |>
    str_to_lower() |>
    str_replace_all("[^a-z ]", "") |>
    str_trim()
}

# Count prescriptions by city
count_prescriptions_by_city <- function(medicare_df, drug_list, suffix = "") {
  medicare_df |>
    mutate(gnrc_name = str_to_lower(Gnrc_Name)) |>
    filter(str_detect(gnrc_name, str_c(str_to_lower(drug_list), collapse = "|"))) |>
    group_by(Prscrbr_City, Prscrbr_State_Abrvtn) |>
    summarise(
      !!paste0("total_prescriptions", suffix) := sum(Tot_Clms, na.rm = TRUE),
      .groups = "drop"
    )
}

# Calculate average AQI
calculate_avg_aqi <- function(aqi_df, suffix = "") {
  aqi_df |>
    group_by(CBSA) |>
    summarise(
      !!paste0("avg_aqi", suffix) := mean(AQI, na.rm = TRUE),
      .groups = "drop"
    )
}

# Match drug cities to AQI cities
match_cities <- function(drug_data, aqi_data, suffix = "") {
  drug_data |>
    rowwise() |>
    mutate(
      matched = list(
        aqi_data |>
          filter(str_detect(get(paste0("aqi_city_clean", suffix)), fixed(get(paste0("drug_city_clean", suffix))))) |>
          select(matches("aqi_city.*"), matches("avg_aqi.*"))
      )
    ) |>
    unnest(matched) |>
    filter(str_detect(get(paste0("aqi_city", suffix)), fixed(Prscrbr_State_Abrvtn, ignore_case = FALSE)))
}

# Merge matched data with population info
merge_with_population <- function(drug_data, population_data, suffix = "") {
  population_data |>
    select(NAME, POPESTIMATE2022) |>
    crossing(drug_data) |>
    filter(str_detect(NAME, fixed(get(paste0("drug_city", suffix)), ignore_case = TRUE))) |>
    filter(str_detect(NAME, fixed(Prscrbr_State_Abrvtn, ignore_case = FALSE))) |>
    arrange(desc(get(paste0("avg_aqi", suffix)))) |>
    distinct(get(paste0("total_prescriptions", suffix)), .keep_all = TRUE)
}

# ----------------------------------------------------------------------
# 2022 Analysis --------------------------------------------------------
# ----------------------------------------------------------------------

# Load data
medicare_2022 <- read_csv("Medicare_Part_D_2022.csv")
aqi_2022 <- read_csv("daily_aqi_by_cbsa_2022.csv")
population_data <- read_csv("cbsa-est2024-alldata.csv")

# Prep AQI data
aqi_data_2022 <- calculate_avg_aqi(aqi_2022, "_2022") |>
  rename(aqi_city_2022 = CBSA) |>
  mutate(aqi_city_clean_2022 = clean_city_name(aqi_city_2022))

# Prep prescription data
drug_data_2022 <- count_prescriptions_by_city(medicare_2022, target_drugs, "_2022") |>
  rename(drug_city_2022 = Prscrbr_City) |>
  mutate(drug_city_clean_2022 = clean_city_name(drug_city_2022))

# Match AQI and prescription data by city and filter
combined_2022 <- drug_data_2022 |>
  rowwise() |>
  mutate(
    matched = list(
      aqi_data_2022 |>
        filter(str_detect(aqi_city_clean_2022, fixed(drug_city_clean_2022))) |>
        select(aqi_city_2022, avg_aqi_2022)
    )
  ) |>
  unnest(matched) |>
  filter(str_detect(aqi_city_2022, fixed(Prscrbr_State_Abrvtn, ignore_case = FALSE))) |>
  filter(total_prescriptions_2022 >= 2000) |>
  distinct(avg_aqi_2022, .keep_all = TRUE)

# Merge with population data
merged_2022 <- merge_with_population(combined_2022, population_data, "_2022") |>
  mutate(
    POPESTIMATE2022 = as.numeric(POPESTIMATE2022),
    total_prescriptions_2022 = as.numeric(total_prescriptions_2022)
  ) |>
  group_by(drug_city_2022, Prscrbr_State_Abrvtn, POPESTIMATE2022) |>
  summarise(
    total_prescriptions_2022 = sum(total_prescriptions_2022, na.rm = TRUE),
    avg_aqi_2022 = mean(avg_aqi_2022, na.rm = TRUE),
    .groups = "drop"
  )

# Add total prescription volume
total_2022 <- medicare_2022 |>
  group_by(Prscrbr_City, Prscrbr_State_Abrvtn) |>
  summarise(
    all_prescriptions_2022 = sum(Tot_Clms, na.rm = TRUE),
    .groups = "drop"
  )

merged_2022 <- merged_2022 |>
  left_join(total_2022, by = c("drug_city_2022" = "Prscrbr_City", "Prscrbr_State_Abrvtn")) |>
  mutate(
    asthma_fraction_2022 = total_prescriptions_2022 / all_prescriptions_2022,
    rate_per_1000_2022 = (total_prescriptions_2022 / POPESTIMATE2022) * 1000
  )

# Calculate total drug costs
cost_all_2022 <- medicare_2022 |>
  group_by(Prscrbr_City, Prscrbr_State_Abrvtn) |>
  summarise(
    total_cost_all_2022 = sum(GE65_Tot_Drug_Cst, na.rm = TRUE),
    .groups = "drop"
  )

cost_asthma_2022 <- medicare_2022 |>
  mutate(gnrc_name = str_to_lower(Gnrc_Name)) |>
  filter(str_detect(gnrc_name, str_c(str_to_lower(target_drugs), collapse = "|"))) |>
  group_by(Prscrbr_City, Prscrbr_State_Abrvtn) |>
  summarise(
    total_cost_asthma_2022 = sum(GE65_Tot_Drug_Cst, na.rm = TRUE),
    .groups = "drop"
  )

# Merge cost data into main dataset
merged_2022 <- merged_2022 |>
  left_join(cost_all_2022, by = c("drug_city_2022" = "Prscrbr_City", "Prscrbr_State_Abrvtn")) |>
  left_join(cost_asthma_2022, by = c("drug_city_2022" = "Prscrbr_City", "Prscrbr_State_Abrvtn")) |>
  mutate(
    asthma_cost_fraction_2022 = total_cost_asthma_2022 / total_cost_all_2022
  )

# Export final dataset

readr::write_csv(merged_2022, "C:/Users/Owner/Downloads/cleaned_2022_data.csv")

# ----------------------------------------------------------------------
# 2019 Analysis --------------------------------------------------------
# ----------------------------------------------------------------------

# Load 2019 data
medicare_2019 <- read_csv("C:/Users/Owner/Downloads/Medicare_Part_D_Prescribers_by_Provider_and_Drug_2019/Medicare_Part_D_Prescribers_by_Provider_and_Drug_2019.csv")
aqi_2019 <- read_csv("C:/Users/Owner/Downloads/daily_aqi_by_cbsa_2019/daily_aqi_by_cbsa_2019.csv")

# Prep data
aqi_data_2019 <- calculate_avg_aqi(aqi_2019, "_2019") |>
  rename(aqi_city_2019 = CBSA) |>
  mutate(aqi_city_clean_2019 = clean_city_name(aqi_city_2019))

drug_data_2019 <- count_prescriptions_by_city(medicare_2019, target_drugs, "_2019") |>
  rename(drug_city_2019 = Prscrbr_City) |>
  mutate(drug_city_clean_2019 = clean_city_name(drug_city_2019))

combined_2019 <- drug_data_2019 |>
  rowwise() |>
  mutate(
    matched = list(
      aqi_data_2019 |>
        filter(str_detect(aqi_city_clean_2019, fixed(drug_city_clean_2019))) |>
        select(aqi_city_2019, avg_aqi_2019)
    )
  ) |>
  unnest(matched) |>
  filter(str_detect(aqi_city_2019, fixed(Prscrbr_State_Abrvtn, ignore_case = FALSE))) |>
  filter(total_prescriptions_2019 >= 2000) |>
  distinct(avg_aqi_2019, .keep_all = TRUE)

merged_2019 <- merge_with_population(combined_2019, population_data, "_2019") |>
  mutate(
    POPESTIMATE2022 = as.numeric(POPESTIMATE2022),
    total_prescriptions_2019 = as.numeric(total_prescriptions_2019)
  ) |>
  group_by(drug_city_2019, Prscrbr_State_Abrvtn, POPESTIMATE2022) |>
  summarise(
    total_prescriptions_2019 = sum(total_prescriptions_2019, na.rm = TRUE),
    avg_aqi_2019 = mean(avg_aqi_2019, na.rm = TRUE),
    .groups = "drop"
  )


# Add total prescriptions
total_2019 <- medicare_2019 |>
  group_by(Prscrbr_City, Prscrbr_State_Abrvtn) |>
  summarise(all_prescriptions_2019 = sum(Tot_Clms, na.rm = TRUE), .groups = "drop")

merged_2019 <- merged_2019 |>
  left_join(total_2019, by = c("drug_city_2019" = "Prscrbr_City", "Prscrbr_State_Abrvtn")) |>
  mutate(asthma_fraction_2019 = total_prescriptions_2019 / all_prescriptions_2019)

# Add costs
cost_all_2019 <- medicare_2019 |>
  group_by(Prscrbr_City, Prscrbr_State_Abrvtn) |>
  summarise(total_cost_all_2019 = sum(GE65_Tot_Drug_Cst, na.rm = TRUE), .groups = "drop")

cost_asthma_2019 <- medicare_2019 |>
  mutate(gnrc_name = str_to_lower(Gnrc_Name)) |>
  filter(str_detect(gnrc_name, str_c(str_to_lower(target_drugs), collapse = "|"))) |>
  group_by(Prscrbr_City, Prscrbr_State_Abrvtn) |>
  summarise(total_cost_asthma_2019 = sum(GE65_Tot_Drug_Cst, na.rm = TRUE), .groups = "drop")

merged_2019 <- merged_2019 |>
  left_join(cost_all_2019, by = c("drug_city_2019" = "Prscrbr_City", "Prscrbr_State_Abrvtn")) |>
  left_join(cost_asthma_2019, by = c("drug_city_2019" = "Prscrbr_City", "Prscrbr_State_Abrvtn")) |>
  mutate(asthma_cost_fraction_2019 = total_cost_asthma_2019 / total_cost_all_2019)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Compare 2022 vs. 2019 ------------------------------------------------
# ----------------------------------------------------------------------

# Load cleaned 2022 data
filtered_2022 <- read_csv("C:/Users/Owner/Downloads/cleaned_2022_data.csv") |>
  mutate(
    city_id = str_c(drug_city_2022, Prscrbr_State_Abrvtn, sep = ", ")
  )

# Ensure 2019 dataset has city_id column
merged_2019 <- merged_2019 |>
  mutate(
    city_id = str_c(drug_city_2019, Prscrbr_State_Abrvtn, sep = ", ")
  )

# Merge datasets by city_id
matched_data_19_22 <- inner_join(
  filtered_2022 |>
    select(city_id, avg_aqi_2022 = avg_aqi_2022, asthma_cost_fraction_2022),
  merged_2019 |>
    select(city_id, avg_aqi_2019, asthma_cost_fraction_2019),
  by = "city_id"
) |>
  mutate(
    change_in_aqi = avg_aqi_2022 - avg_aqi_2019,
    change_in_cost_fraction = asthma_cost_fraction_2022 - asthma_cost_fraction_2019
  )
readr::write_csv(matched_data_19_22, "C:/Users/Owner/Downloads/matched_data_2019_2022.csv")


