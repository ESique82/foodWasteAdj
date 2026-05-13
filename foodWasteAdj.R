# foodWasteAdj.R
# ---------------
# Waste-adjusted LCA impact calculation for food systems.
#
# Adjusts per-production-unit environmental impacts to per-consumed-unit
# impacts using food loss and waste rates, with full Monte Carlo uncertainty
# propagation across both production uncertainty (from lcaStats-R) and
# waste rate uncertainty (derived from observed study variation).
#
# Designed to connect directly to run_monte_carlo() output from lcaStats-R:
#
#   inventory |>
#     batch_from_dataframe() |>          # lcaStats-R / pedigree_R.R
#     run_monte_carlo(n = 10000) |>      # lcaStats-R / monte_carlo.R
#     adjust_for_waste()                 # this module
#
# WASTE DATA SOURCES
# ------------------
# Three built-in sources are supported:
#
#   "fao_flw"         FAO Food Loss and Waste Database (2021+)
#                     Most comprehensive. Uncertainty derived from observed
#                     variation across studies. Requires path to downloaded CSV.
#                     Download: https://www.fao.org/platform-food-loss-waste/flw-data/en/
#
#   "gustavsson_2011" Gustavsson et al. (2011) FAO global estimates.
#                     Regional aggregates, 7 food categories.
#                     Built-in — no file path required.
#
#   "custom"          User-supplied data frame in standard format.
#
# REFERENCES
# ----------
# Gustavsson, J., Cederberg, C., Sonesson, U., van Otterdijk, R., &
#   Meybeck, A. (2011). Global food losses and food waste: extent, causes
#   and prevention. FAO, Rome.
#
# FAO (2021). Food Loss and Waste Database. Technical Platform on the
#   Measurement and Reduction of Food Loss and Waste.
#   https://www.fao.org/platform-food-loss-waste/flw-data/en/
#
# Parfitt, J., Barthel, M., & Macnaughton, S. (2010). Food waste within
#   food supply chains: quantification and potential for change to 2050.
#   Philosophical Transactions of the Royal Society B, 365, 3065-3081.


# ---------------------------------------------------------------------------
# Stage mapping — FAO food_supply_stage to our six standardised stages
# ---------------------------------------------------------------------------

#' @keywords internal
FAO_STAGE_MAP <- c(
  "Farm"               = "Farm",
  "Pre-harvest"        = "Farm",
  "Harvest"            = "Farm",
  "Post-harvest"       = "Post-harvest",
  "Storage"            = "Post-harvest",
  "Stacking"           = "Post-harvest",
  "Grading"            = "Post-harvest",
  "Collector"          = "Post-harvest",
  "Packing"            = "Post-harvest",
  "Processing"         = "Processing",
  "Transport"          = "Distribution",
  "Distribution"       = "Distribution",
  "Wholesale"          = "Distribution",
  "Trader"             = "Distribution",
  "Export"             = "Distribution",
  "Retail"             = "Retail",
  "Market"             = "Retail",
  "Households"         = "Consumer",
  "Food Services"      = "Consumer"
)

#' @keywords internal
STAGES <- c("Farm", "Post-harvest", "Processing",
            "Distribution", "Retail", "Consumer")


# ---------------------------------------------------------------------------
# Beta distribution parameter helpers
# ---------------------------------------------------------------------------

#' Convert mean and SD to Beta distribution shape parameters
#'
#' Waste rates are proportions bounded [0,1]. Beta is the correct
#' distribution — unlike lognormal which is unbounded above.
#'
#' @param mean Numeric. Mean waste rate (0-1).
#' @param sd   Numeric. Standard deviation of waste rate (0-1).
#' @return Named numeric vector with elements alpha and beta.
#' @keywords internal
beta_params <- function(mean, sd) {
  if (mean <= 0 | mean >= 1) {
    stop(sprintf("mean must be in (0, 1), got %.4f", mean))
  }
  if (sd <= 0) {
    stop(sprintf("sd must be > 0, got %.4f", sd))
  }
  var   <- sd^2
  alpha <- mean * (mean * (1 - mean) / var - 1)
  beta  <- (1 - mean) * (mean * (1 - mean) / var - 1)
  if (alpha <= 0 | beta <= 0) {
    stop(sprintf(
      "Invalid Beta parameters (alpha=%.3f, beta=%.3f). SD may be too large relative to mean.",
      alpha, beta
    ))
  }
  c(alpha = alpha, beta = beta)
}


# ---------------------------------------------------------------------------
# load_waste_data()
# ---------------------------------------------------------------------------

#' Load and process food waste rate data
#'
#' Returns a standardised waste rate table with columns:
#' food_category, stage, waste_mean, waste_sd, n_studies, source.
#'
#' @param source Character. One of "fao_flw", "gustavsson_2011", or "custom".
#' @param file Character. Path to FAO FLW CSV. Required for "fao_flw".
#' @param commodity Character or NULL. Filter to commodity e.g. "Wheat".
#' @param country Character or NULL. Filter to country e.g. "Mexico".
#' @param region Character or NULL. Filter to region for gustavsson_2011.
#' @param min_studies Integer. Minimum studies per stage. Default 3.
#' @param custom_df Data frame. Required for source = "custom".
#' @return Standardised waste rate data frame.
#' @export
load_waste_data <- function(source,
                             file        = NULL,
                             commodity   = NULL,
                             country     = NULL,
                             region      = NULL,
                             min_studies = 3,
                             custom_df   = NULL) {

  source <- match.arg(source, c("fao_flw", "gustavsson_2011", "custom"))

  if (source == "fao_flw") {
    .load_fao_flw(file, commodity, country, min_studies)
  } else if (source == "gustavsson_2011") {
    .load_gustavsson(region)
  } else {
    .load_custom(custom_df)
  }
}


# ---------------------------------------------------------------------------
# Internal loaders
# ---------------------------------------------------------------------------

#' @keywords internal
.load_fao_flw <- function(file, commodity, country, min_studies) {

  if (is.null(file) || !file.exists(file)) {
    stop(
      "file path required for source = 'fao_flw'.\n",
      "Download from: https://www.fao.org/platform-food-loss-waste/flw-data/en/"
    )
  }

  df <- read.csv(file, stringsAsFactors = FALSE)
  
  df$country   <- trimws(df$country)
  df$commodity <- trimws(df$commodity)

  required <- c("commodity", "food_supply_stage", "loss_percentage", "country")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("FAO FLW file missing columns: %s",
                 paste(missing, collapse = ", ")))
  }

  if (!is.null(commodity)) {
    df <- df[grepl(commodity, df$commodity, ignore.case = TRUE), ]
    if (nrow(df) == 0) stop(sprintf("No data found for commodity '%s'.", commodity))
  }

  if (!is.null(country)) {
    df_country <- df[grepl(country, df$country, ignore.case = TRUE), ]
    if (nrow(df_country) == 0) {
      warning(sprintf("No data for country '%s'. Using all available data.", country))
    } else {
      df <- df_country
    }
  }

  df$stage <- FAO_STAGE_MAP[df$food_supply_stage]
  df <- df[!is.na(df$stage), ]

  df$loss_percentage <- suppressWarnings(as.numeric(df$loss_percentage))
  df <- df[!is.na(df$loss_percentage) & df$loss_percentage > 0, ]
  df$waste_rate <- df$loss_percentage / 100

  food_label <- if (!is.null(commodity)) commodity else "All commodities"

  result <- do.call(rbind, lapply(STAGES, function(s) {
    rows <- df[df$stage == s, ]
    n    <- nrow(rows)
    if (n < min_studies) return(NULL)
    data.frame(
      food_category = food_label,
      stage         = s,
      waste_mean    = mean(rows$waste_rate, na.rm = TRUE),
      waste_sd      = if (n > 1) sd(rows$waste_rate, na.rm = TRUE)
                      else rows$waste_rate * 0.2,
      n_studies     = n,
      source        = "FAO FLW Database",
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(result) || nrow(result) == 0) {
    stop(sprintf(
      "No stages had >= %d studies. Try lowering min_studies or broadening filters.",
      min_studies
    ))
  }

  result
}


#' @keywords internal
.load_gustavsson <- function(region) {

  gustavsson <- data.frame(
    food_category = rep(c(
      "Cereals", "Roots and tubers", "Oilseeds and pulses",
      "Fruits and vegetables", "Meat", "Dairy", "Fish and seafood"
    ), each = 6),
    stage = rep(STAGES, 7),
    waste_mean = c(
      0.020, 0.010, 0.010, 0.010, 0.020, 0.250,
      0.200, 0.050, 0.050, 0.020, 0.050, 0.170,
      0.030, 0.010, 0.020, 0.010, 0.020, 0.040,
      0.200, 0.050, 0.020, 0.050, 0.100, 0.190,
      0.030, 0.010, 0.050, 0.020, 0.050, 0.110,
      0.010, 0.010, 0.020, 0.010, 0.050, 0.070,
      0.090, 0.020, 0.120, 0.020, 0.060, 0.110
    ),
    waste_sd = c(
      0.005, 0.003, 0.003, 0.003, 0.005, 0.050,
      0.050, 0.010, 0.010, 0.005, 0.010, 0.040,
      0.007, 0.003, 0.005, 0.003, 0.005, 0.010,
      0.050, 0.010, 0.005, 0.010, 0.020, 0.040,
      0.007, 0.003, 0.010, 0.005, 0.010, 0.030,
      0.003, 0.003, 0.005, 0.003, 0.010, 0.020,
      0.020, 0.005, 0.030, 0.005, 0.015, 0.030
    ),
    n_studies = 1L,
    source    = "Gustavsson et al. (2011)",
    stringsAsFactors = FALSE
  )

  if (!is.null(region) && !grepl("europe", region, ignore.case = TRUE)) {
    warning(sprintf(
      "Region '%s' not yet built-in for Gustavsson. Returning Europe values.",
      region
    ))
  }

  gustavsson
}


#' @keywords internal
.load_custom <- function(custom_df) {

  if (is.null(custom_df)) stop("custom_df required when source = 'custom'.")

  required <- c("food_category", "stage", "waste_mean", "waste_sd")
  missing  <- setdiff(required, names(custom_df))
  if (length(missing) > 0) {
    stop(sprintf("custom_df missing columns: %s", paste(missing, collapse = ", ")))
  }

  invalid_stages <- setdiff(unique(custom_df$stage), STAGES)
  if (length(invalid_stages) > 0) {
    stop(sprintf("Invalid stage(s): %s. Must be one of: %s",
                 paste(invalid_stages, collapse = ", "),
                 paste(STAGES, collapse = ", ")))
  }

  if (!("n_studies" %in% names(custom_df))) custom_df$n_studies <- NA_integer_
  if (!("source"    %in% names(custom_df))) custom_df$source    <- "Custom"

  custom_df[, c("food_category", "stage", "waste_mean", "waste_sd",
                "n_studies", "source")]
}

# adjust_for_waste()
# Main function: applies waste adjustment to MC output
# ------------------------------------------------------------

adjust_for_waste <- function(mc_results,
                             waste_df,
                             food_category_map,
                             seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  food_items <- mc_results$meta$food_items
  n          <- mc_results$meta$n
  stages_mc  <- mc_results$meta$stages
  
  # Validate food_category_map covers all food items
  unmapped <- setdiff(food_items, names(food_category_map))
  if (length(unmapped) > 0) {
    stop(sprintf(
      "food_category_map missing entries for: %s",
      paste(unmapped, collapse = ", ")
    ))
  }
  
  # Sum MC samples across supply chain stages -> total impact per food item
  # mc_results$samples is a list keyed by stage; each element is matrix [n x food_items]
  produced_samples <- matrix(0, nrow = n, ncol = length(food_items),
                             dimnames = list(NULL, food_items))
  
  for (stage in stages_mc) {
    if (stage %in% names(mc_results$samples)) {
      produced_samples <- produced_samples + mc_results$samples[[stage]]
    }
  }
  
  # For each food item: draw Beta waste rates per stage, compute consumed impact
  consumed_samples    <- produced_samples          # initialise; overwritten below
  waste_rates_applied <- list()
  
  for (f in food_items) {
    
    category <- food_category_map[f]
    wd       <- waste_df[waste_df$food_category == category, ]
    
    if (nrow(wd) == 0) {
      warning(sprintf(
        "No waste data for category '%s' (food item: '%s'). No adjustment applied.",
        category, f
      ))
      waste_rates_applied[[f]] <- list()
      next
    }
    
    # Cumulative retention across all stages: prod(1 - w_s)
    retention   <- rep(1, n)
    stage_rates <- list()
    
    for (s in STAGES) {
      stage_row <- wd[wd$stage == s, ]
      if (nrow(stage_row) == 0) next   # no waste data for this stage
      
      m  <- stage_row$waste_mean[1]
      sd <- stage_row$waste_sd[1]
      
      # Fit Beta distribution; fall back to fixed rate if parameters invalid
      params <- tryCatch(
        beta_params(m, sd),
        error = function(e) {
          warning(sprintf(
            "Beta params failed for %s / %s: %s. Using fixed rate.",
            f, s, e$message
          ))
          NULL
        }
      )
      
      w <- if (is.null(params)) rep(m, n) else rbeta(n, params["alpha"], params["beta"])
      
      retention   <- retention * (1 - w)
      stage_rates[[s]] <- w
    }
    
    retention <- pmax(retention, 1e-6)   # guard against division by zero
    consumed_samples[, f]   <- produced_samples[, f] / retention
    waste_rates_applied[[f]] <- stage_rates
  }
  
  # Helper: summarise a sample matrix into a data frame
  .summarise_samples <- function(mat) {
    do.call(rbind, lapply(food_items, function(f) {
      x <- mat[, f]
      data.frame(
        food_item = f,
        mean      = mean(x),
        sd        = sd(x),
        median    = median(x),
        q2.5      = quantile(x, 0.025),
        q97.5     = quantile(x, 0.975),
        cv        = sd(x) / mean(x),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }))
  }
  
  list(
    produced          = .summarise_samples(produced_samples),
    consumed          = .summarise_samples(consumed_samples),
    samples           = list(produced = produced_samples,
                             consumed = consumed_samples),
    waste_rates       = waste_rates_applied,
    waste_df          = waste_df,
    food_category_map = food_category_map,
    meta              = mc_results$meta
  )
}


# ------------------------------------------------------------
# waste_summary_table()
# Prints a formatted console summary of waste rates applied
# ------------------------------------------------------------

waste_summary_table <- function(adjusted) {
  
  food_items <- names(adjusted$waste_rates)
  
  cat("\n── Waste Adjustment Summary ──────────────────────────────────────\n\n")
  
  for (f in food_items) {
    
    cat(sprintf("Food item : %s\n", f))
    cat(sprintf("Category  : %s\n", adjusted$food_category_map[f]))
    
    stage_rates <- adjusted$waste_rates[[f]]
    
    if (length(stage_rates) == 0) {
      cat("  (no waste data — impact unchanged)\n\n")
      next
    }
    
    for (s in names(stage_rates)) {
      w <- stage_rates[[s]]
      cat(sprintf("  %-15s  mean = %5.1f%%   95%% CI: [%4.1f%% – %4.1f%%]\n",
                  s,
                  mean(w) * 100,
                  quantile(w, 0.025) * 100,
                  quantile(w, 0.975) * 100))
    }
    
    prod_mean   <- adjusted$produced$mean[adjusted$produced$food_item == f]
    cons_mean   <- adjusted$consumed$mean[adjusted$consumed$food_item == f]
    pct_change  <- (cons_mean / prod_mean - 1) * 100
    
    cat(sprintf("  %-15s  %.4f\n", "Produced mean", prod_mean))
    cat(sprintf("  %-15s  %.4f   (+%.1f%%)\n\n", "Consumed mean", cons_mean, pct_change))
  }
  
  invisible(adjusted)
}


# ------------------------------------------------------------
# plot_waste_comparison()
# Side-by-side produced vs consumed point-range plot
# ------------------------------------------------------------

plot_waste_comparison <- function(adjusted,
                                  x_label = "Environmental impact",
                                  title   = "Produced vs Consumed Impact") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required. Install with: install.packages('ggplot2')")
  }
  
  library(ggplot2)
  
  prod       <- adjusted$produced;  prod$basis <- "Produced"
  cons       <- adjusted$consumed;  cons$basis <- "Consumed"
  plot_df    <- rbind(prod, cons)
  plot_df$basis <- factor(plot_df$basis, levels = c("Produced", "Consumed"))
  
  ggplot(plot_df,
         aes(x = mean, y = food_item, colour = basis, shape = basis)) +
    geom_pointrange(
      aes(xmin = q2.5, xmax = q97.5),
      position = position_dodge(width = 0.5),
      linewidth = 0.7,
      size      = 0.6
    ) +
    scale_colour_manual(
      values = c("Produced" = "#2166ac", "Consumed" = "#d73027")
    ) +
    scale_shape_manual(
      values = c("Produced" = 16, "Consumed" = 17)
    ) +
    labs(
      x      = x_label,
      y      = NULL,
      title  = title,
      colour = NULL,
      shape  = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "top",
      panel.grid.minor = element_blank()
    )
}


# ------------------------------------------------------------
# export_waste_results()
# Exports results to a multi-sheet Excel workbook
# ------------------------------------------------------------

export_waste_results <- function(adjusted,
                                 file = "waste_results.xlsx") {
  
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("openxlsx is required. Install with: install.packages('openxlsx')")
  }
  
  library(openxlsx)
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "Produced");  writeData(wb, "Produced", adjusted$produced)
  addWorksheet(wb, "Consumed");  writeData(wb, "Consumed", adjusted$consumed)
  
  # Waste rates applied — one row per food item × stage
  waste_summary <- do.call(rbind, lapply(names(adjusted$waste_rates), function(f) {
    rates <- adjusted$waste_rates[[f]]
    if (length(rates) == 0) return(NULL)
    do.call(rbind, lapply(names(rates), function(s) {
      w <- rates[[s]]
      data.frame(
        food_item   = f,
        stage       = s,
        waste_mean  = mean(w),
        waste_sd    = sd(w),
        waste_q2.5  = quantile(w, 0.025),
        waste_q97.5 = quantile(w, 0.975),
        stringsAsFactors = FALSE,
        row.names   = NULL
      )
    }))
  }))
  
  addWorksheet(wb, "Waste rates applied")
  writeData(wb, "Waste rates applied", waste_summary)
  
  addWorksheet(wb, "Waste source data")
  writeData(wb, "Waste source data", adjusted$waste_df)
  
  saveWorkbook(wb, file, overwrite = TRUE)
  message(sprintf("Exported to: %s", file))
  invisible(file)
}

# ---------------------------------------------------------------------------
# adjust_for_waste_simple() — standalone entry point (no lcaStats-R needed)
# ---------------------------------------------------------------------------

#' Adjust point-estimate LCA impacts for food loss and waste
#'
#' A standalone entry point that does NOT require \code{run_monte_carlo()}
#' output from lcaStats-R. Use this when you have impact figures from an
#' external source (e.g. Agribalyse, Ecoinvent, a published LCA study) and
#' want to adjust them from per-kg-produced to per-kg-consumed.
#'
#' The production impact is treated as a fixed point estimate. Only waste
#' rate uncertainty is propagated through the Monte Carlo simulation. This
#' is appropriate when:
#' \itemize{
#'   \item You are checking or exploring figures from a waste database
#'   \item You do not have pedigree scores or GSD² values for the impacts
#'   \item Your interest is specifically in the sensitivity of consumed-unit
#'         impacts to waste rate uncertainty
#' }
#'
#' If you also have uncertainty information for the production impacts
#' (e.g. a coefficient of variation), supply it via the \code{impact_cv}
#' column in \code{impacts} and set \code{propagate_production_uncertainty
#' = TRUE}.
#'
#' @param impacts Data frame with at minimum two columns:
#'   \describe{
#'     \item{food_item}{Character. Name of each food item.}
#'     \item{impact_mean}{Numeric. Mean LCA impact per kg produced.
#'       Any impact category (GWP, water use, land use, etc.) and any
#'       functional unit, as long as it is expressed per kg \emph{produced}.}
#'     \item{impact_cv}{Numeric (optional). Coefficient of variation for
#'       the production impact (SD / mean). If supplied and
#'       \code{propagate_production_uncertainty = TRUE}, production
#'       uncertainty is modelled as Lognormal and propagated jointly with
#'       waste rate uncertainty. Ignored otherwise.}
#'   }
#' @param waste_df Data frame. Output of \code{\link{load_waste_data}}.
#' @param food_category_map Named character vector mapping food items to
#'   waste categories. Names are food items (matching \code{impacts$food_item});
#'   values are category names matching \code{waste_df$food_category}.
#'   Example: \code{c("Wheat bread" = "Wheat", "Beef burger" = "Meat")}.
#' @param n Integer. Number of Monte Carlo iterations. Default 10,000.
#' @param seed Integer or NULL. Random seed for reproducibility. Default 42.
#' @param propagate_production_uncertainty Logical. If TRUE and
#'   \code{impact_cv} is present in \code{impacts}, production impacts are
#'   drawn from a Lognormal distribution rather than treated as fixed.
#'   Default FALSE.
#'
#' @return A list with the same structure as \code{\link{adjust_for_waste}},
#'   so all downstream functions (\code{waste_summary_table},
#'   \code{plot_waste_comparison}, \code{export_waste_results}) work
#'   identically:
#'   \describe{
#'     \item{\code{$produced}}{Summary statistics for impact per kg produced.}
#'     \item{\code{$consumed}}{Summary statistics for impact per kg consumed.}
#'     \item{\code{$samples}}{List with \code{produced} and \code{consumed}
#'       matrices (n rows × food items).}
#'     \item{\code{$waste_rates}}{Drawn waste rate samples per food item
#'       and stage.}
#'     \item{\code{$waste_df}}{The waste data used as input.}
#'     \item{\code{$meta}}{Run metadata.}
#'   }
#'
#' @examples
#' \dontrun{
#' # --- Minimal example: point estimates from Agribalyse ---
#'
#' impacts <- data.frame(
#'   food_item   = c("Wheat bread", "Beef burger", "Whole milk"),
#'   impact_mean = c(1.124, 17.4, 1.87)   # kg CO2-eq per kg produced
#' )
#'
#' food_category_map <- c(
#'   "Wheat bread" = "Cereals & bakery",
#'   "Beef burger" = "Meat",
#'   "Whole milk"  = "Dairy"
#' )
#'
#' waste_df <- load_waste_data("gustavsson_2011",
#'                             food_category = "Cereals & bakery")
#'
#' results <- adjust_for_waste_simple(
#'   impacts           = impacts,
#'   waste_df          = waste_df,
#'   food_category_map = food_category_map,
#'   n                 = 10000,
#'   seed              = 42
#' )
#'
#' waste_summary_table(results)
#' plot_waste_comparison(results, x_label = "GWP100 (kg CO2-eq per kg food)")
#' export_waste_results(results, file = "simple_waste_results.xlsx")
#'
#'
#' # --- With production uncertainty (if you have a CV) ---
#'
#' impacts_with_cv <- data.frame(
#'   food_item   = c("Wheat bread", "Beef burger"),
#'   impact_mean = c(1.124, 17.4),
#'   impact_cv   = c(0.12, 0.25)   # 12% and 25% coefficient of variation
#' )
#'
#' results_full <- adjust_for_waste_simple(
#'   impacts                          = impacts_with_cv,
#'   waste_df                         = waste_df,
#'   food_category_map                = food_category_map,
#'   propagate_production_uncertainty = TRUE
#' )
#' }
#'
#' @seealso \code{\link{adjust_for_waste}} for the full pipeline using
#'   \code{run_monte_carlo()} output from lcaStats-R.
#'
#' @export
adjust_for_waste_simple <- function(impacts,
                                    waste_df,
                                    food_category_map,
                                    n    = 10000,
                                    seed = 42,
                                    propagate_production_uncertainty = FALSE) {
  
  # --- Input validation -------------------------------------------------------
  
  if (!is.data.frame(impacts)) {
    stop("'impacts' must be a data frame.")
  }
  required_cols <- c("food_item", "impact_mean")
  missing_cols  <- setdiff(required_cols, names(impacts))
  if (length(missing_cols) > 0) {
    stop(sprintf("'impacts' is missing required column(s): %s",
                 paste(missing_cols, collapse = ", ")))
  }
  if (any(impacts$impact_mean <= 0)) {
    stop("All values in 'impact_mean' must be positive.")
  }
  if (!is.character(food_category_map) || is.null(names(food_category_map))) {
    stop("'food_category_map' must be a named character vector.")
  }
  unmapped <- setdiff(impacts$food_item, names(food_category_map))
  if (length(unmapped) > 0) {
    stop(sprintf(
      "The following food items have no entry in food_category_map: %s",
      paste(unmapped, collapse = ", ")))
  }
  
  if (propagate_production_uncertainty) {
    if (!"impact_cv" %in% names(impacts)) {
      warning(paste0(
        "propagate_production_uncertainty = TRUE but 'impact_cv' column not ",
        "found in impacts. Falling back to fixed point estimates."))
      propagate_production_uncertainty <- FALSE
    } else if (any(impacts$impact_cv <= 0, na.rm = TRUE)) {
      stop("All values in 'impact_cv' must be positive.")
    }
  }
  
  # --- Set seed ---------------------------------------------------------------
  
  if (!is.null(seed)) set.seed(seed)
  
  food_items <- impacts$food_item
  n_items    <- length(food_items)
  
  # --- Step 1: build produced_samples matrix (n × n_items) -------------------
  #
  # If propagate_production_uncertainty = FALSE  → replicate the point estimate
  #   n times (constant column). Only waste uncertainty contributes.
  #
  # If propagate_production_uncertainty = TRUE   → draw from Lognormal using
  #   mean and CV. Both sources of uncertainty contribute.
  
  produced_samples <- matrix(NA_real_, nrow = n, ncol = n_items,
                             dimnames = list(NULL, food_items))
  
  for (j in seq_along(food_items)) {
    mu <- impacts$impact_mean[j]
    
    if (propagate_production_uncertainty && "impact_cv" %in% names(impacts)) {
      cv       <- impacts$impact_cv[j]
      sigma_ln <- sqrt(log(cv^2 + 1))
      mu_ln    <- log(mu) - sigma_ln^2 / 2
      produced_samples[, j] <- rlnorm(n, meanlog = mu_ln, sdlog = sigma_ln)
    } else {
      produced_samples[, j] <- rep(mu, n)
    }
  }
  
  # --- Step 2: draw waste rates and compute consumed samples ------------------
  #
  # Identical logic to adjust_for_waste() — Beta draws per stage per iteration.
  
  consumed_samples <- matrix(NA_real_, nrow = n, ncol = n_items,
                             dimnames = list(NULL, food_items))
  
  waste_rate_draws <- vector("list", n_items)
  names(waste_rate_draws) <- food_items
  
  for (j in seq_along(food_items)) {
    
    category <- food_category_map[food_items[j]]
    item_waste <- waste_df[waste_df$food_category == category, ]
    
    if (nrow(item_waste) == 0) {
      stop(sprintf(
        "No waste data found for category '%s' (mapped from '%s'). ",
        "Check food_category_map and waste_df$food_category.",
        category, food_items[j]))
    }
    
    # Draw waste rates for each stage — Beta(alpha, beta) per iteration
    survival <- matrix(1, nrow = n, ncol = 1)   # cumulative survival fraction
    stage_draws <- list()
    
    for (s in STAGES) {
      row <- item_waste[item_waste$stage == s, ]
      
      if (nrow(row) == 0 || is.na(row$waste_mean) || is.na(row$waste_sd)) {
        # Stage not available for this category — assume zero loss
        stage_draws[[s]] <- rep(0, n)
        next
      }
      
      params        <- beta_params(row$waste_mean, row$waste_sd)
      draws         <- rbeta(n, params["alpha"], params["beta"])
      stage_draws[[s]] <- draws
      survival      <- survival * (1 - draws)
    }
    
    waste_rate_draws[[j]] <- as.data.frame(stage_draws)
    
    # Consumed impact = produced impact / survival fraction
    consumed_samples[, j] <- produced_samples[, j] / as.numeric(survival)
  }
  
  # --- Step 3: summarise ------------------------------------------------------
  
  summarise_samples <- function(mat, ci = 0.95) {
    alpha <- (1 - ci) / 2
    do.call(rbind, lapply(colnames(mat), function(item) {
      x <- mat[, item]
      data.frame(
        food_item = item,
        mean      = mean(x),
        sd        = sd(x),
        median    = median(x),
        ci_lower  = quantile(x, alpha),
        ci_upper  = quantile(x, 1 - alpha),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }))
  }
  
  produced_summary <- summarise_samples(produced_samples)
  consumed_summary <- summarise_samples(consumed_samples)
  
  # --- Step 4: waste rate summary (mirrors adjust_for_waste output) -----------
  
  waste_rates_summary <- do.call(rbind, lapply(seq_along(food_items), function(j) {
    draws <- waste_rate_draws[[j]]
    do.call(rbind, lapply(names(draws), function(s) {
      x <- draws[[s]]
      data.frame(
        food_item = food_items[j],
        stage     = s,
        mean      = mean(x),
        sd        = sd(x),
        ci_lower  = quantile(x, 0.025),
        ci_upper  = quantile(x, 0.975),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }))
  }))
  
  # --- Return -----------------------------------------------------------------
  
  list(
    produced    = produced_summary,
    consumed    = consumed_summary,
    samples     = list(
      produced  = produced_samples,
      consumed  = consumed_samples
    ),
    waste_rates = waste_rates_summary,
    waste_df    = waste_df,
    meta        = list(
      n                                = n,
      seed                             = seed,
      food_items                       = food_items,
      food_category_map                = food_category_map,
      propagate_production_uncertainty = propagate_production_uncertainty,
      entry_point                      = "adjust_for_waste_simple",
      timestamp                        = Sys.time()
    )
  )
}