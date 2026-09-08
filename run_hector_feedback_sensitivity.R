#!/usr/bin/env Rscript

# Hector 3.5.0 referee sensitivity for the V59.1/V60 carbon pathway.
# Purpose: test the two programme-contingent feedback objections raised after
# the first strict paired experiment:
#   (1) stored-carbon reversal must be zero in the removal-off world;
#   (2) exported permafrost must respond to each arm's own temperature.
#
# The script preserves the exact common Hector state at end-2025 and produces:
#   - the legacy strict pair for reference;
#   - a reversal-zero-off sensitivity with prescribed permafrost;
#   - a dynamic-permafrost sensitivity with reversal retained in both arms;
#   - the combined programme-effect pair: reversal=0 in off and permafrost
#     generated from each arm's own Hector TAS trajectory.
#
# Dynamic permafrost uses the project's pre-existing V52 causal formulation:
# 18 PgC per degC above 1.20 C, running-peak commitment, 50-year e-folding
# release, and no reinjection of historical commitment already implicit in the
# common 2025 atmospheric state.  The climate/permafrost coupling is solved by
# damped fixed-point iteration separately for each arm.

suppressPackageStartupMessages(library(hector))
suppressPackageStartupMessages(library(jsonlite))

options(error = function() {
  cat("\n=== R TRACEBACK ===\n", file = stderr())
  traceback(20)
  quit(save = "no", status = 1, runLast = FALSE)
})

say <- function(...) {
  cat(sprintf(...), "\n")
  flush.console()
}
fail <- function(...) stop(sprintf(...), call. = FALSE)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
traj_path <- file.path(root, "external_validation_net_co2_trajectory.csv")
if (!file.exists(traj_path)) fail("Trajectory not found: %s", traj_path)
if (as.character(packageVersion("hector")) != "3.5.0") {
  fail("Expected Hector 3.5.0, found %s", as.character(packageVersion("hector")))
}

x0 <- read.csv(traj_path, check.names = FALSE, stringsAsFactors = FALSE)
need <- c("year", "gross_co2_gtco2", "permafrost_co2_gtco2",
          "stored_carbon_reversal_gtco2", "cdr_gtco2", "net_co2_gtco2")
missing_cols <- setdiff(need, names(x0))
if (length(missing_cols)) fail("Missing columns: %s", paste(missing_cols, collapse = ", "))
x0 <- x0[order(x0$year), ]
if (min(x0$year) != 2026 || max(x0$year) != 2183) fail("Expected 2026-2183 trajectory")
if (anyDuplicated(x0$year)) fail("Duplicate trajectory years")
mb <- x0$gross_co2_gtco2 + x0$permafrost_co2_gtco2 + x0$stored_carbon_reversal_gtco2 - x0$cdr_gtco2
if (max(abs(mb - x0$net_co2_gtco2)) >= 1e-9) fail("Mass-balance check failed")

# Extend only the exogenous physical programme components.  Dynamic permafrost
# is NOT extended from the terminal exported value; it is recalculated below.
last <- x0[nrow(x0), ]
ext_year <- 2184:2300
ext <- data.frame(
  year = ext_year,
  gross_co2_gtco2 = last$gross_co2_gtco2,
  permafrost_co2_gtco2 = last$permafrost_co2_gtco2,
  stored_carbon_reversal_gtco2 = last$stored_carbon_reversal_gtco2,
  cdr_gtco2 = last$cdr_gtco2,
  net_co2_gtco2 = last$net_co2_gtco2
)
x <- rbind(x0[, names(ext)], ext)
yrs <- as.numeric(x$year)
stopifnot(identical(yrs, as.numeric(2026:2300)))

co2_to_c <- 12.011 / 44.009
gtc_to_gtco2 <- 44.009 / 12.011
gtco2_per_ppm <- 2.124 * (44.009 / 12.011)

ini <- system.file("input", "hector_ssp245.ini", package = "hector")
if (!nzchar(ini) || !file.exists(ini)) fail("hector_ssp245.ini not found")
unit_ffi <- getunits(FFI_EMISSIONS())
unit_luc <- getunits(LUC_EMISSIONS())
unit_luc_uptake <- getunits(LUC_UPTAKE())
unit_daccs <- getunits(DACCS_UPTAKE())

# Permafrost assumptions copied from the project's existing V52 implementation.
PF_CARBON_PGC_PER_C <- 18.0
PF_ACTIVATION_C <- 1.20
PF_TAU_YR <- 50.0
PF_DAMPING <- 0.50
PF_TOL_GTCO2_YR <- 1e-5
PF_MAX_ITER <- 30L

set_checked <- function(core, dates, var, values, unit, label) {
  tryCatch(setvar(core, dates, var, values, unit),
           error = function(e) fail("setvar failed for %s: %s", label, conditionMessage(e)))
}

fetch_scalar <- function(df, variable_name, year) {
  z <- df[df$variable == variable_name & df$year == year, "value"]
  if (length(z) != 1 || !is.finite(z)) fail("Bad %s at %d", variable_name, year)
  as.numeric(z)
}

# Return annual GtCO2/yr release for 2026:2300 using a common 2025 baseline.
permafrost_from_tas <- function(tas_year, tas_value, common_tas_2025) {
  if (!identical(as.numeric(tas_year), yrs)) fail("TAS years do not match 2026-2300")
  running_peak <- common_tas_2025
  prev_commit <- PF_CARBON_PGC_PER_C * max(0, common_tas_2025 - PF_ACTIVATION_C)
  pool <- 0.0
  release_fraction <- 1.0 - exp(-1.0 / PF_TAU_YR)
  out <- numeric(length(yrs))
  for (i in seq_along(yrs)) {
    t <- as.numeric(tas_value[i])
    running_peak <- max(running_peak, t)
    target <- PF_CARBON_PGC_PER_C * max(0, running_peak - PF_ACTIVATION_C)
    new_commit <- max(0, target - prev_commit)
    pool <- pool + new_commit
    release_pgc <- pool * release_fraction
    pool <- pool - release_pgc
    out[i] <- release_pgc * gtc_to_gtco2
    prev_commit <- target
  }
  out
}

# One Hector execution from the exact end-2025 common state.
run_once <- function(label, removal_on, reversal_off_zero,
                     permafrost_gtco2, reversal_override = NULL) {
  core <- newcore(ini, suppresslogging = TRUE, name = label)
  on.exit(try(shutdown(core), silent = TRUE), add = TRUE)
  invisible(run(core, 2025))
  hist <- fetchvars(core, 2025, vars = c(CONCENTRATIONS_CO2(), GLOBAL_TAS()))
  hist_co2 <- fetch_scalar(hist, CONCENTRATIONS_CO2(), 2025)
  hist_tas <- fetch_scalar(hist, GLOBAL_TAS(), 2025)

  reversal <- x$stored_carbon_reversal_gtco2
  if (!removal_on && reversal_off_zero) reversal <- rep(0, length(yrs))
  if (!is.null(reversal_override)) reversal <- reversal_override

  source_gtco2 <- x$gross_co2_gtco2 + permafrost_gtco2 + reversal
  if (any(source_gtco2 < 0) || any(!is.finite(source_gtco2))) fail("Bad source pathway in %s", label)
  ffi_c <- source_gtco2 * co2_to_c
  cdr_c <- if (removal_on) x$cdr_gtco2 * co2_to_c else rep(0, length(yrs))
  zero_c <- rep(0, length(yrs))

  set_checked(core, yrs, FFI_EMISSIONS(), ffi_c, unit_ffi, paste0(label, " FFI"))
  set_checked(core, yrs, LUC_EMISSIONS(), zero_c, unit_luc, paste0(label, " LUC emissions"))
  set_checked(core, yrs, LUC_UPTAKE(), zero_c, unit_luc_uptake, paste0(label, " LUC uptake"))
  set_checked(core, yrs, DACCS_UPTAKE(), cdr_c, unit_daccs, paste0(label, " DACCS"))
  invisible(run(core, 2300))

  out <- fetchvars(core, yrs, vars = c(CONCENTRATIONS_CO2(), GLOBAL_TAS()))
  tas <- out[out$variable == GLOBAL_TAS(), c("year", "value")]
  tas <- tas[order(tas$year), ]
  co2 <- out[out$variable == CONCENTRATIONS_CO2(), c("year", "value")]
  co2 <- co2[order(co2$year), ]

  audit_years <- c(2026, 2040, 2050, 2075, 2100, 2150, 2156, 2183, 2200, 2250, 2300)
  audit <- fetchvars(core, audit_years,
                     vars = c(FFI_EMISSIONS(), LUC_EMISSIONS(), LUC_UPTAKE(), DACCS_UPTAKE()))
  audit$case <- label

  list(output = out, tas = tas, co2 = co2, audit = audit,
       hist_co2 = hist_co2, hist_tas = hist_tas,
       reversal = reversal, permafrost = permafrost_gtco2)
}

# Iteratively close temperature -> permafrost -> temperature for one arm.
solve_dynamic <- function(label, removal_on, reversal_off_zero) {
  pf <- rep(0, length(yrs))
  conv <- data.frame(iteration = integer(), max_abs_flux_change_gtco2_yr = numeric(),
                     cumulative_pf_gtco2 = numeric(), peak_pf_gtco2_yr = numeric())
  common_tas <- NA_real_
  final <- NULL
  for (it in seq_len(PF_MAX_ITER)) {
    rr <- run_once(sprintf("%s_iter_%02d", label, it), removal_on,
                   reversal_off_zero, pf)
    if (it == 1L) common_tas <- rr$hist_tas
    target_pf <- permafrost_from_tas(rr$tas$year, rr$tas$value, common_tas)
    err <- max(abs(target_pf - pf))
    conv <- rbind(conv, data.frame(iteration = it,
                                   max_abs_flux_change_gtco2_yr = err,
                                   cumulative_pf_gtco2 = sum(target_pf),
                                   peak_pf_gtco2_yr = max(target_pf)))
    say("%s iteration %d: max dPF=%.8g GtCO2/yr; cumulative PF=%.6f GtCO2",
        label, it, err, sum(target_pf))
    if (err < PF_TOL_GTCO2_YR) {
      pf <- target_pf
      final <- run_once(paste0(label, "_final"), removal_on,
                        reversal_off_zero, pf)
      break
    }
    pf <- PF_DAMPING * target_pf + (1 - PF_DAMPING) * pf
  }
  if (is.null(final)) fail("Dynamic permafrost did not converge for %s", label)
  # Verify final fixed point against final TAS.
  check_pf <- permafrost_from_tas(final$tas$year, final$tas$value, final$hist_tas)
  final_err <- max(abs(check_pf - pf))
  if (final_err > 5 * PF_TOL_GTCO2_YR) {
    fail("Final dynamic permafrost closure failed for %s: %.9g", label, final_err)
  }
  final$permafrost <- check_pf
  final$convergence <- conv
  final$final_pf_closure_error <- final_err
  final
}

say("Hector feedback-sensitivity experiment starting")
say("Hector version %s", as.character(packageVersion("hector")))
say("Canonical CDR through 2183: %.6f GtCO2", sum(x$cdr_gtco2[x$year <= 2183]))
say("Exported stored-carbon reversal through 2300: %.6f GtCO2", sum(x$stored_carbon_reversal_gtco2))
say("Exported permafrost through 2300: %.6f GtCO2", sum(x$permafrost_co2_gtco2))

# Strict legacy pair, preserved for direct comparison with the first experiment.
legacy_on <- run_once("legacy_on", TRUE, FALSE, x$permafrost_co2_gtco2)
legacy_off <- run_once("legacy_off", FALSE, FALSE, x$permafrost_co2_gtco2)

# Isolate the stored-reversal counterfactual correction only.
no_rev_off <- run_once("no_reversal_off", FALSE, TRUE, x$permafrost_co2_gtco2)

# Dynamic permafrost sensitivity.  On arm is common to both dynamic comparisons.
dyn_on <- solve_dynamic("dynamic_on", TRUE, FALSE)
dyn_off_same_rev <- solve_dynamic("dynamic_off_same_reversal", FALSE, FALSE)
dyn_off_no_rev <- solve_dynamic("dynamic_off_no_reversal", FALSE, TRUE)

# Common-state check across all final arms.
all_cases <- list(legacy_on = legacy_on, legacy_off = legacy_off,
                  no_reversal_off = no_rev_off, dynamic_on = dyn_on,
                  dynamic_off_same_reversal = dyn_off_same_rev,
                  dynamic_off_no_reversal = dyn_off_no_rev)
hist_co2 <- sapply(all_cases, function(z) z$hist_co2)
hist_tas <- sapply(all_cases, function(z) z$hist_tas)
if (max(hist_co2) - min(hist_co2) > 1e-12 || max(hist_tas) - min(hist_tas) > 1e-12) {
  fail("Common-state audit failed across sensitivity arms")
}

co2_series <- function(z) {
  d <- z$co2
  setNames(as.numeric(d$value), as.character(d$year))
}

make_pair <- function(label, on_case, off_case) {
  onv <- co2_series(on_case); offv <- co2_series(off_case)
  delta <- offv - onv
  cum_cdr <- cumsum(x$cdr_gtco2)
  frac <- delta * gtco2_per_ppm / cum_cdr
  data.frame(experiment = label, year = yrs,
             co2_on_ppm = as.numeric(onv), co2_off_ppm = as.numeric(offv),
             delta_co2_ppm = as.numeric(delta),
             cumulative_cdr_gtco2 = cum_cdr,
             response_fraction = as.numeric(frac))
}

pairs <- rbind(
  make_pair("legacy_strict_pair", legacy_on, legacy_off),
  make_pair("reversal_zero_off_only", legacy_on, no_rev_off),
  make_pair("dynamic_permafrost_only", dyn_on, dyn_off_same_rev),
  make_pair("combined_programme_effect", dyn_on, dyn_off_no_rev)
)
write.csv(pairs, "hector_feedback_sensitivity_pair.csv", row.names = FALSE)

sel_years <- c(2040, 2050, 2075, 2100, 2150, 2156, 2200, 2250, 2300)
selected <- pairs[pairs$year %in% sel_years, ]
write.csv(selected, "hector_feedback_sensitivity_selected.csv", row.names = FALSE)

# Long outputs and forcing audits for independent reproduction/audit.
long <- do.call(rbind, lapply(names(all_cases), function(nm) {
  d <- all_cases[[nm]]$output
  d$case <- nm
  d
}))
write.csv(long, "hector_feedback_sensitivity_long.csv", row.names = FALSE)

audits <- do.call(rbind, lapply(names(all_cases), function(nm) all_cases[[nm]]$audit))
write.csv(audits, "hector_feedback_sensitivity_forcing_audit.csv", row.names = FALSE)

pf <- data.frame(
  year = yrs,
  exported_permafrost_gtco2 = x$permafrost_co2_gtco2,
  dynamic_on_gtco2 = dyn_on$permafrost,
  dynamic_off_same_reversal_gtco2 = dyn_off_same_rev$permafrost,
  dynamic_off_no_reversal_gtco2 = dyn_off_no_rev$permafrost
)
write.csv(pf, "hector_feedback_sensitivity_permafrost.csv", row.names = FALSE)

conv <- rbind(
  transform(dyn_on$convergence, case = "dynamic_on"),
  transform(dyn_off_same_rev$convergence, case = "dynamic_off_same_reversal"),
  transform(dyn_off_no_rev$convergence, case = "dynamic_off_no_reversal")
)
write.csv(conv, "hector_feedback_sensitivity_convergence.csv", row.names = FALSE)

getrow <- function(exp, year) {
  z <- pairs[pairs$experiment == exp & pairs$year == year, ]
  if (nrow(z) != 1) fail("Missing pair row %s %d", exp, year)
  z
}
legacy2300 <- getrow("legacy_strict_pair", 2300)
rev2300 <- getrow("reversal_zero_off_only", 2300)
dyn2300 <- getrow("dynamic_permafrost_only", 2300)
comb2300 <- getrow("combined_programme_effect", 2300)
comb2156 <- getrow("combined_programme_effect", 2156)
comb2200 <- getrow("combined_programme_effect", 2200)

cdr_2156 <- comb2156$cumulative_cdr_gtco2
cdr_2200 <- comb2200$cumulative_cdr_gtco2
attrib_burden_2156 <- comb2156$delta_co2_ppm * gtco2_per_ppm
attrib_burden_2200 <- comb2200$delta_co2_ppm * gtco2_per_ppm

summary <- list(
  model = "Hector",
  version = "3.5.0",
  experiment = "programme-contingent feedback sensitivity after strict paired attribution",
  common_state_year = 2025,
  common_state_co2_ppm = unname(hist_co2[1]),
  common_state_tas_degc = unname(hist_tas[1]),
  common_state_max_co2_error_ppm = max(hist_co2) - min(hist_co2),
  common_state_max_tas_error_degc = max(hist_tas) - min(hist_tas),
  reporting_gtco2_per_ppm = gtco2_per_ppm,
  permafrost_method = list(
    source = "project V52 causal running-peak commitment formulation",
    carbon_feedback_pgc_per_c = PF_CARBON_PGC_PER_C,
    activation_warming_c = PF_ACTIVATION_C,
    thaw_release_tau_yr = PF_TAU_YR,
    fixed_point_damping = PF_DAMPING,
    convergence_tolerance_gtco2_yr = PF_TOL_GTCO2_YR,
    historical_commitment_treatment = "excluded because implicit in the exact common 2025 state"
  ),
  cumulative_exported_reversal_2026_2300_gtco2 = sum(x$stored_carbon_reversal_gtco2),
  cumulative_exported_permafrost_2026_2300_gtco2 = sum(x$permafrost_co2_gtco2),
  cumulative_dynamic_permafrost_on_2026_2300_gtco2 = sum(dyn_on$permafrost),
  cumulative_dynamic_permafrost_off_same_reversal_2026_2300_gtco2 = sum(dyn_off_same_rev$permafrost),
  cumulative_dynamic_permafrost_off_no_reversal_2026_2300_gtco2 = sum(dyn_off_no_rev$permafrost),
  response_fraction_2300 = list(
    legacy_strict_pair = legacy2300$response_fraction,
    reversal_zero_off_only = rev2300$response_fraction,
    dynamic_permafrost_only = dyn2300$response_fraction,
    combined_programme_effect = comb2300$response_fraction
  ),
  delta_co2_2300_ppm = list(
    legacy_strict_pair = legacy2300$delta_co2_ppm,
    reversal_zero_off_only = rev2300$delta_co2_ppm,
    dynamic_permafrost_only = dyn2300$delta_co2_ppm,
    combined_programme_effect = comb2300$delta_co2_ppm
  ),
  reversal_window_combined = list(
    cumulative_cdr_added_2156_to_2200_gtco2 = cdr_2200 - cdr_2156,
    attributable_burden_2156_gtco2 = attrib_burden_2156,
    attributable_burden_2200_gtco2 = attrib_burden_2200,
    attributable_burden_change_2156_to_2200_gtco2 = attrib_burden_2200 - attrib_burden_2156,
    co2_on_2156_ppm = comb2156$co2_on_ppm,
    co2_on_2200_ppm = comb2200$co2_on_ppm,
    co2_off_2156_ppm = comb2156$co2_off_ppm,
    co2_off_2200_ppm = comb2200$co2_off_ppm
  ),
  fixed_point_closure_error_gtco2_yr = list(
    dynamic_on = dyn_on$final_pf_closure_error,
    dynamic_off_same_reversal = dyn_off_same_rev$final_pf_closure_error,
    dynamic_off_no_reversal = dyn_off_no_rev$final_pf_closure_error
  ),
  interpretation_boundary = paste0(
    "legacy_strict_pair is a direct CDR counterfactual conditional on exported feedbacks; ",
    "combined_programme_effect is a programme-effect sensitivity in which storage reversal is contingent on prior CDR ",
    "and permafrost responds to each arm's own Hector temperature."
  )
)
write_json(summary, "hector_feedback_sensitivity_summary.json", pretty = TRUE,
           auto_unbox = TRUE, digits = NA)

say("\n=== HECTOR FEEDBACK SENSITIVITY COMPLETE ===")
say("Legacy 2300 fraction: %.9f", legacy2300$response_fraction)
say("Reversal-zero-off 2300 fraction: %.9f", rev2300$response_fraction)
say("Dynamic-PF-only 2300 fraction: %.9f", dyn2300$response_fraction)
say("Combined programme-effect 2300 fraction: %.9f", comb2300$response_fraction)
say("Combined 2156->2200 CDR added: %.6f GtCO2", cdr_2200 - cdr_2156)
say("Combined 2156->2200 attributable burden change: %.6f GtCO2", attrib_burden_2200 - attrib_burden_2156)
say("Outputs written: hector_feedback_sensitivity_*.csv/json")
