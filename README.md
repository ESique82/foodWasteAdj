# foodWasteAdj

**Waste-adjusted LCA impacts for food systems — R module**

Adjusts per-production-unit environmental impacts to per-consumed-unit impacts
using food loss and waste rates, with full Monte Carlo uncertainty propagation.

Designed as a direct downstream module for
[lcaStats-R](https://github.com/ESique82/lcaStats-R).

---

## The problem it solves

Life cycle assessments of food typically report impacts **per kg produced**.
But food is lost and wasted at every stage of the supply chain — from farm to
consumer — so the true environmental cost **per kg eaten** is always higher.

This module bridges that gap. It takes Monte Carlo LCA output from `lcaStats-R`
and applies stage-level waste rates to produce consumed-unit impacts with
uncertainty propagated from both sources simultaneously:

```
Uncertainty in production impacts   (pedigree matrix / lcaStats-R)
          +
Uncertainty in waste rates          (observed variation across FAO studies)
          ↓
Combined uncertainty in consumed-unit impact
```

---

## Installation

No package installation required. Source the file directly:

```r
source("foodWasteAdj.R")
```

### Dependencies

| Package | Used for | Install |
|---|---|---|
| `ggplot2` | `plot_waste_comparison()` | `install.packages("ggplot2")` |
| `openxlsx` | `export_waste_results()` | `install.packages("openxlsx")` |

Both are optional — core adjustment functions work without them.

### Waste data

The FAO Food Loss and Waste Database CSV is required for `source = "fao_flw"`.
Download free from:
<https://www.fao.org/platform-food-loss-waste/flw-data/en/>

---

## Quick start

```r
# 1. Source dependencies
lca_path <- "path/to/lcaStats-R/"
source(paste0(lca_path, "pedigree_R.R"))
source(paste0(lca_path, "monte_carlo.R"))
source("foodWasteAdj.R")

# 2. Build inventory and run Monte Carlo (lcaStats-R)
inventory <- data.frame(
  food_item                = c("Wheat bread", "Wheat bread"),
  stage                    = c("Agriculture", "Processing"),
  amount                   = c(0.8, 0.3),
  reliability              = c(2L, 2L),
  completeness             = c(2L, 2L),
  temporal_correlation     = c(3L, 2L),
  geographical_correlation = c(1L, 1L),
  technology_correlation   = c(2L, 2L),
  basic_var                = rep(0.0006, 2)
)

results <- inventory |>
  batch_from_dataframe() |>
  run_monte_carlo(n = 10000)

# 3. Load waste data
waste_df <- load_waste_data("fao_flw",
                             file      = "Data_FAO.csv",
                             commodity = "Wheat")

# 4. Adjust for waste
adjusted <- adjust_for_waste(
  mc_results        = results,
  waste_df          = waste_df,
  food_category_map = c("Wheat bread" = "Wheat")
)

# 5. Explore results
adjusted$produced          # summary: per kg produced
adjusted$consumed          # summary: per kg consumed
waste_summary_table(adjusted)
plot_waste_comparison(adjusted, x_label = "GWP100 (kg CO2-eq per kg food)")
```

See `example_waste_adj.R` for a full two-product worked example.

---
---

## Standalone use — no lcaStats-R required

`adjust_for_waste_simple()` is an alternative entry point that does **not** require Monte Carlo output from lcaStats-R. Use it when you have impact figures from an external source (Agribalyse, Ecoinvent, a published LCA study, or a waste database) and want to adjust them from per-kg-produced to per-kg-consumed.

The production impact is treated as a fixed point estimate. Only waste rate uncertainty is propagated through the simulation. This is the correct approach when:
- You are exploring or sense-checking figures from a waste or LCA database
- You do not have pedigree scores or GSD² values for the impacts
- Your interest is specifically the sensitivity of consumed-unit impacts to waste rate uncertainty

### Minimal example

```r
source("foodWasteAdj.R")

# Impact figures from any external source
impacts <- data.frame(
  food_item   = c("Wheat bread", "Beef burger", "Whole milk"),
  impact_mean = c(1.124, 17.4, 1.87)   # kg CO2-eq per kg produced
)

# Map food items to waste categories in your chosen waste data source
food_category_map <- c(
  "Wheat bread" = "Cereals & bakery",
  "Beef burger" = "Meat",
  "Whole milk"  = "Dairy"
)

# Load waste rates
waste_df <- load_waste_data("gustavsson_2011")

# Run adjustment — only waste uncertainty propagated
results <- adjust_for_waste_simple(
  impacts           = impacts,
  waste_df          = waste_df,
  food_category_map = food_category_map,
  n                 = 10000,
  seed              = 42
)

# All downstream functions work identically to the full pipeline
waste_summary_table(results)
plot_waste_comparison(results, x_label = "GWP100 (kg CO2-eq per kg food)")
export_waste_results(results, file = "results.xlsx")
```

### Optional: propagate production uncertainty too

If you also have a coefficient of variation (CV = SD / mean) for the production impacts, you can propagate both sources of uncertainty jointly:

```r
impacts_with_cv <- data.frame(
  food_item   = c("Wheat bread", "Beef burger"),
  impact_mean = c(1.124, 17.4),
  impact_cv   = c(0.12, 0.25)   # 12% and 25% CV
)

results <- adjust_for_waste_simple(
  impacts                          = impacts_with_cv,
  waste_df                         = waste_df,
  food_category_map                = food_category_map,
  propagate_production_uncertainty = TRUE
)
```

### When to use which entry point

| Situation | Function to use |
|---|---|
| You have inventory data with pedigree scores | `adjust_for_waste()` via full lcaStats-R pipeline |
| You have point estimates from an external DB | `adjust_for_waste_simple()` |
| You have point estimates + CV from literature | `adjust_for_waste_simple(..., propagate_production_uncertainty = TRUE)` |
## Functions

| Function | Description |
|---|---|
| `load_waste_data()` | Load waste rates from FAO FLW, Gustavsson (2011), or custom data |
| `adjust_for_waste()` | Apply waste adjustment with Monte Carlo uncertainty propagation |
| `waste_summary_table()` | Print formatted console summary of waste rates and impact change |
| `plot_waste_comparison()` | Side-by-side produced vs consumed point-range plot |
| `export_waste_results()` | Export results to multi-sheet Excel workbook |

---

## Waste data sources

### FAO Food Loss and Waste Database (`"fao_flw"`)
Over 30,000 data points from 700+ published studies. Uncertainty is derived
from observed variation across studies — the most defensible approach for
LCA uncertainty propagation. Requires downloaded CSV.

### Gustavsson et al. (2011) (`"gustavsson_2011"`)
Built-in regional aggregates for 7 food categories × 6 supply chain stages.
No file required. Note: currently returns European values only.

### Custom (`"custom"`)
Supply your own data frame with columns:
`food_category`, `stage`, `waste_mean`, `waste_sd`.

---

## Supply chain stages

| Stage | Covers |
|---|---|
| Farm | Growing, harvesting, pre-harvest losses |
| Post-harvest | Drying, sorting, grading, packing, storage |
| Processing | Manufacturing and transformation losses |
| Distribution | Transport, wholesale, export |
| Retail | Supermarkets and markets |
| Consumer | Households and food service |

---

## Key methodological notes

- **Beta distribution** is used for waste rates — it is bounded [0, 1] and
  therefore correct for proportions. Lognormal (used for LCA exchanges) is
  not appropriate here.
- **Uncertainty is derived from study variation** in the FAO database, not
  from assumed coefficients of variation.
- **Stages with fewer than 3 studies** are excluded by default (`min_studies = 3`).
  Lower this threshold with caution.
- **Latin American meat data** in the FAO database is reported only at
  "Whole supply chain" level — no stage breakdown is available. Use global
  FAO meat data or Gustavsson (2011) as fallback for regional studies.

---

## References

Gustavsson, J., Cederberg, C., Sonesson, U., van Otterdijk, R., & Meybeck, A.
(2011). *Global food losses and food waste: extent, causes and prevention.*
FAO, Rome.

FAO (2021). *Food Loss and Waste Database.* Technical Platform on the
Measurement and Reduction of Food Loss and Waste.
<https://www.fao.org/platform-food-loss-waste/flw-data/en/>

Parfitt, J., Barthel, M., & Macnaughton, S. (2010). Food waste within food
supply chains: quantification and potential for change to 2050.
*Philosophical Transactions of the Royal Society B*, 365, 3065–3081.

---

## Related

- [lcaStats-R](https://github.com/ESique82/lcaStats-R) — pedigree matrix
  uncertainty and Monte Carlo simulation for LCA
- [food-innovation-HA-team](https://github.com/food-innovation-HA-team) —
  Harper Food Innovation: Digital toolset

---

*Harper Food Innovation: Digital · Harper Adams University*
