# example_waste_adj.R
# -------------------
# Worked example: waste-adjusted LCA impacts for wheat bread and beef burger.
#
# Demonstrates the full pipeline from inventory data to consumed-unit impacts
# with Monte Carlo uncertainty propagation across both production uncertainty
# (from lcaStats-R) and waste rate uncertainty (from FAO FLW Database).
#
# Requirements:
#   - lcaStats-R scripts (pedigree_R.R, monte_carlo.R)
#   - foodWasteAdj.R
#   - Data_FAO.csv (FAO Food Loss and Waste Database)
#     Download: https://www.fao.org/platform-food-loss-waste/flw-data/en/
#   - R packages: ggplot2, openxlsx


# ---------------------------------------------------------------------------
# 0. Setup — adjust lca_path to your lcaStats-R folder
# ---------------------------------------------------------------------------

lca_path <- "/home/eric/projects/lcaStats-R/"  # <-- change this

source(paste0(lca_path, "pedigree_R.R"))
source(paste0(lca_path, "monte_carlo.R"))
source("foodWasteAdj.R")


# ---------------------------------------------------------------------------
# 1. Build inventory
# ---------------------------------------------------------------------------
# Each row is one exchange: a single input/output at a specific supply chain
# stage for a specific food item.
#
# Columns:
#   food_item                 : product name (used as label throughout)
#   stage                     : supply chain stage (Agriculture, Processing, etc.)
#   amount                    : mean impact value (e.g. kg CO2-eq per kg food)
#   reliability               : DQI score 1-5 (1 = best, 5 = worst)
#   completeness              : DQI score 1-5
#   temporal_correlation      : DQI score 1-5
#   geographical_correlation  : DQI score 1-5
#   technology_correlation    : DQI score 1-5
#   basic_var                 : basic variance of the measurement (default 0.0006)

inventory <- data.frame(
  food_item                = c("Wheat bread", "Wheat bread",
                               "Beef burger", "Beef burger"),
  stage                    = c("Agriculture", "Processing",
                               "Agriculture", "Processing"),
  amount                   = c(0.8, 0.3, 22.0, 1.5),
  reliability              = c(2L, 2L, 3L, 2L),
  completeness             = c(2L, 2L, 2L, 2L),
  temporal_correlation     = c(3L, 2L, 3L, 2L),
  geographical_correlation = c(1L, 1L, 2L, 1L),
  technology_correlation   = c(2L, 2L, 3L, 2L),
  basic_var                = rep(0.0006, 4)
)


# ---------------------------------------------------------------------------
# 2. Run Monte Carlo simulation (lcaStats-R)
# ---------------------------------------------------------------------------
# batch_from_dataframe() converts the inventory into lcaStats-R format.
# run_monte_carlo() draws n = 10,000 samples from lognormal distributions
# parameterised by the pedigree DQI scores.

set.seed(42)   # for reproducibility

results <- inventory |>
  batch_from_dataframe() |>
  run_monte_carlo(n = 10000)

# Inspect the production impact distributions
results$total          # mean, SD, 95% CI per food item
plot_mc_total(results, x_label = "GWP100 (kg CO2-eq per kg food produced)")


# ---------------------------------------------------------------------------
# 3. Load waste rate data (FAO FLW Database)
# ---------------------------------------------------------------------------
# load_waste_data() reads the FAO CSV, maps FAO stage names to our six
# standardised stages, and calculates mean and SD from observed study variation.
#
# min_studies = 3 (default): only include stages with at least 3 data points.
#
# Note: Latin American meat data is only available at "Whole supply chain"
# level in the FAO database. Global meat data is used here as a fallback.

waste_wheat <- load_waste_data(
  source    = "fao_flw",
  file      = "Data_FAO.csv",
  commodity = "Wheat"
)

waste_meat <- load_waste_data(
  source    = "fao_flw",
  file      = "Data_FAO.csv",
  commodity = "Meat of"
)

waste_combined <- rbind(waste_wheat, waste_meat)

# Inspect waste rates loaded
print(waste_combined)


# ---------------------------------------------------------------------------
# 4. Adjust for food waste
# ---------------------------------------------------------------------------
# adjust_for_waste() draws Beta-distributed waste rates for each stage and
# each of the 10,000 MC iterations, computes cumulative survival fractions,
# and divides produced impacts by survival fraction to get consumed impacts.
#
# food_category_map links food item names in the inventory to food category
# names in the waste data.

adjusted <- adjust_for_waste(
  mc_results        = results,
  waste_df          = waste_combined,
  food_category_map = c("Wheat bread" = "Wheat",
                        "Beef burger" = "Meat of"),
  seed              = 42
)


# ---------------------------------------------------------------------------
# 5. Inspect results
# ---------------------------------------------------------------------------

# Summary tables
adjusted$produced   # per kg produced: mean, SD, median, 95% CI, CV
adjusted$consumed   # per kg consumed: same structure

# Console summary with waste rates applied per stage
waste_summary_table(adjusted)


# ---------------------------------------------------------------------------
# 6. Visualise
# ---------------------------------------------------------------------------

plot_waste_comparison(
  adjusted,
  x_label = "GWP100 (kg CO2-eq per kg food)",
  title   = "Produced vs Consumed Impact — Wheat Bread and Beef Burger"
)

# Save plot
# library(ggplot2)
# ggsave("produced_vs_consumed.png", width = 8, height = 5, dpi = 150)


# ---------------------------------------------------------------------------
# 7. Export to Excel
# ---------------------------------------------------------------------------
# Creates a workbook with four sheets:
#   Produced            — summary stats per kg produced
#   Consumed            — summary stats per kg consumed
#   Waste rates applied — mean, SD, 95% CI of drawn waste rates per stage
#   Waste source data   — raw waste rate inputs used

export_waste_results(adjusted, file = "waste_results.xlsx")


# ---------------------------------------------------------------------------
# 8. Optional: compare with Gustavsson (2011) waste data
# ---------------------------------------------------------------------------
# Gustavsson et al. (2011) provides broad regional averages for 7 food
# categories. Useful as a sensitivity check or when FAO data is sparse.

waste_gustavsson <- rbind(
  load_waste_data("gustavsson_2011")[
    load_waste_data("gustavsson_2011")$food_category == "Cereals", ],
  load_waste_data("gustavsson_2011")[
    load_waste_data("gustavsson_2011")$food_category == "Meat", ]
)

# Remap categories to match food items
adjusted_g2011 <- adjust_for_waste(
  mc_results        = results,
  waste_df          = waste_gustavsson,
  food_category_map = c("Wheat bread" = "Cereals",
                        "Beef burger" = "Meat"),
  seed              = 42
)

# Compare consumed means: FAO vs Gustavsson
comparison <- data.frame(
  food_item      = adjusted$consumed$food_item,
  consumed_FAO   = round(adjusted$consumed$mean, 3),
  consumed_G2011 = round(adjusted_g2011$consumed$mean, 3)
)
print(comparison)
