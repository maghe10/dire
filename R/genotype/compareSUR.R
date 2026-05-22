# ============================================================
# Compare phenotype SUR vs genotype SUR
# and include breakpoint-distance information from disk diffusion
#
# Input files expected in:
#   paste(processedRootRassembly, "phenotype", sep = "/")
#
# Files:
#   - phenotype_SUR_wide.csv
#   - genotype_SUR_wide.csv
#   - genotype_SUR_combined_per_sample.csv
#   - phenotype_millimeters.csv
#   - phenotype_breakpoints.csv
#
# Main idea:
#   closer to the breakpoint generally means a more uncertain phenotype
#   so we quantify distance-to-breakpoint and include it in the analysis
# ============================================================
source("genotype/genotypeCommon.R")

comparePhenotypeAndGenotypeSUR <- function(
    input_dir = paste(processedRootRassembly, "phenotype", sep = "/"),
    out_dir   = input_dir
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(purrr)
    library(tibble)
  })
  
  # ----------------------------------------------------------
  # Helpers
  # ----------------------------------------------------------
  normalize_sample_id <- function(x) {
    x |>
      as.character() |>
      str_trim() |>
      str_replace_all("\\.0$", "") |>
      str_replace_all("^sample[_\\- ]*", "") |>
      str_replace_all("[^A-Za-z0-9]", "") |>
      toupper()
  }
  
  normalize_antibiotic_name <- function(x) {
    x |>
      as.character() |>
      str_trim() |>
      toupper()
  }
  
  detect_sample_column <- function(df) {
    nms <- names(df)
    hit <- nms[tolower(nms) %in% c("sample", "sampleid", "sample_id", "isolate", "isolateid", "id")]
    if (length(hit) == 0) stop("Could not find sample identifier column.")
    hit[1]
  }
  
  read_sur_csv <- function(path) {
    readr::read_delim(
      file = path,
      delim = ";",
      show_col_types = FALSE,
      trim_ws = TRUE,
      na = c("", "NA", "NaN")
    )
  }
  
  write_semicolon_csv <- function(df, path) {
    readr::write_delim(df, file = path, delim = ";", na = "")
  }
  
  standardize_sur_wide <- function(df, antibiotics = NULL) {
    sample_col <- detect_sample_column(df)
    
    df <- df |>
      rename(sample = all_of(sample_col)) |>
      mutate(sample = normalize_sample_id(sample))
    
    names(df) <- c("sample", normalize_antibiotic_name(names(df)[-1]))
    
    if (is.null(antibiotics)) {
      antibiotics <- setdiff(names(df), "sample")
    } else {
      antibiotics <- normalize_antibiotic_name(antibiotics)
      antibiotics <- intersect(antibiotics, names(df))
    }
    
    df |>
      select(sample, all_of(antibiotics)) |>
      mutate(across(-sample, ~ as.character(.x) |> str_trim() |> toupper()))
  }
  
  # ----------------------------------------------------------
  # Parse breakpoint table into a standard long format
  # Expected output columns:
  #   antibiotic, s_min_mm, r_max_mm
  #
  # Interpretation:
  #   >= s_min_mm  => S
  #   <= r_max_mm  => R
  #   between      => U
  #
  # This is EUCAST-style for zone diameters.
  # ----------------------------------------------------------
  parse_breakpoints_long <- function(bp_raw, antibiotics = NULL) {
    names(bp_raw) <- normalize_antibiotic_name(names(bp_raw))
    
    if (!is.null(antibiotics)) {
      antibiotics <- normalize_antibiotic_name(antibiotics)
    }
    
    # Case 1: already long-ish with antibiotic and breakpoint columns
    lower_map <- setNames(names(bp_raw), tolower(names(bp_raw)))
    
    if ("antibiotic" %in% names(lower_map)) {
      out <- bp_raw |>
        rename(
          antibiotic = all_of(lower_map[["antibiotic"]])
        )
      
      names(out) <- normalize_antibiotic_name(names(out))
      if (!"ANTIBIOTIC" %in% names(out)) {
        stop("Breakpoint table parsing failed.")
      }
      
      out <- out |>
        rename(antibiotic = ANTIBIOTIC) |>
        mutate(antibiotic = normalize_antibiotic_name(antibiotic))
      
      s_candidates <- names(out)[tolower(names(out)) %in% c("s_min_mm", "sbreakpoint", "s_breakpoint", "s", "s_mm", "smin", "sminmm")]
      r_candidates <- names(out)[tolower(names(out)) %in% c("r_max_mm", "rbreakpoint", "r_breakpoint", "r", "r_mm", "rmax", "rmaxmm")]
      
      if (length(s_candidates) >= 1 && length(r_candidates) >= 1) {
        out <- out |>
          transmute(
            antibiotic = antibiotic,
            s_min_mm = suppressWarnings(as.numeric(.data[[s_candidates[1]]])),
            r_max_mm = suppressWarnings(as.numeric(.data[[r_candidates[1]]]))
          )
        
        if (!is.null(antibiotics)) {
          out <- out |>
            filter(antibiotic %in% antibiotics)
        }
        
        return(out)
      }
    }
    
    # Case 2: wide table with one row per breakpoint type and one column per antibiotic
    first_col <- names(bp_raw)[1]
    first_values <- bp_raw[[first_col]] |> as.character() |> str_trim() |> tolower()
    
    if (any(first_values %in% c("s", "susceptible", "s_min_mm", "r", "resistant", "r_max_mm"))) {
      out <- bp_raw |>
        rename(kind = all_of(first_col)) |>
        mutate(kind = str_trim(tolower(as.character(kind)))) |>
        pivot_longer(
          cols = -kind,
          names_to = "antibiotic",
          values_to = "value"
        ) |>
        mutate(
          antibiotic = normalize_antibiotic_name(antibiotic),
          value = suppressWarnings(as.numeric(value)),
          kind = case_when(
            kind %in% c("s", "susceptible", "s_min_mm", "s_breakpoint", "sbreakpoint") ~ "s_min_mm",
            kind %in% c("r", "resistant", "r_max_mm", "r_breakpoint", "rbreakpoint") ~ "r_max_mm",
            TRUE ~ kind
          )
        ) |>
        filter(kind %in% c("s_min_mm", "r_max_mm")) |>
        pivot_wider(
          names_from = kind,
          values_from = value
        )
      
      if (!is.null(antibiotics)) {
        out <- out |>
          filter(antibiotic %in% antibiotics)
      }
      
      return(out)
    }
    
    stop("Could not parse phenotype_breakpoints.csv into columns antibiotic, s_min_mm, r_max_mm.")
  }
  
  # ----------------------------------------------------------
  # Parse millimeter table into long format:
  #   sample, antibiotic, mm
  # ----------------------------------------------------------
  parse_mm_long <- function(mm_raw, antibiotics = NULL) {
    mm_raw <- as_tibble(mm_raw)
    sample_col <- detect_sample_column(mm_raw)
    
    mm_raw <- mm_raw |>
      rename(sample = all_of(sample_col)) |>
      mutate(sample = normalize_sample_id(sample))
    
    names(mm_raw) <- c("sample", normalize_antibiotic_name(names(mm_raw)[-1]))
    
    if (!is.null(antibiotics)) {
      antibiotics <- normalize_antibiotic_name(antibiotics)
      keep_abs <- intersect(antibiotics, names(mm_raw))
    } else {
      keep_abs <- setdiff(names(mm_raw), "sample")
    }
    
    mm_raw |>
      select(sample, all_of(keep_abs)) |>
      pivot_longer(
        cols = -sample,
        names_to = "antibiotic",
        values_to = "mm"
      ) |>
      mutate(
        antibiotic = normalize_antibiotic_name(antibiotic),
        mm = suppressWarnings(as.numeric(mm))
      )
  }
  
  # ----------------------------------------------------------
  # Compute distance-to-breakpoint features
  #
  # nearest_breakpoint_distance:
  #   absolute distance in mm to nearest S/R breakpoint
  #   smaller => closer to decision boundary => more uncertain
  #
  # interval_width:
  #   size of U interval, if present
  #
  # margin_from_called_region:
  #   positive margin inside currently assigned phenotype region
  #   smaller => closer to border
  # ----------------------------------------------------------
  add_breakpoint_distance_features <- function(df) {
    df |>
      mutate(
        nearest_breakpoint_distance = pmin(
          abs(mm - s_min_mm),
          abs(mm - r_max_mm),
          na.rm = FALSE
        ),
        interval_width = s_min_mm - r_max_mm - 1,
        margin_from_called_region = case_when(
          phenotype_SUR == "S" ~ mm - s_min_mm,
          phenotype_SUR == "R" ~ r_max_mm - mm,
          phenotype_SUR == "U" ~ pmin(mm - r_max_mm, s_min_mm - mm, na.rm = FALSE),
          TRUE ~ NA_real_
        ),
        is_near_breakpoint_0mm = !is.na(nearest_breakpoint_distance) & nearest_breakpoint_distance <= 0,
        is_near_breakpoint_1mm = !is.na(nearest_breakpoint_distance) & nearest_breakpoint_distance <= 1,
        is_near_breakpoint_2mm = !is.na(nearest_breakpoint_distance) & nearest_breakpoint_distance <= 2,
        is_near_breakpoint_3mm = !is.na(nearest_breakpoint_distance) & nearest_breakpoint_distance <= 3
      )
  }
  
  # ----------------------------------------------------------
  # File paths
  # ----------------------------------------------------------
  phenotype_file   <- file.path(input_dir, "phenotype_SUR_wide.csv")
  genotype_file    <- file.path(input_dir, "genotype_SUR_wide.csv")
  combined_file    <- file.path(input_dir, "genotype_SUR_combined_per_sample.csv")
  millimeter_file  <- file.path(input_dir, "phenotype_millimeters.csv")
  breakpoint_file  <- file.path(input_dir, "phenotype_breakpoints.csv")
  
  stopifnot(file.exists(phenotype_file))
  stopifnot(file.exists(genotype_file))
  stopifnot(file.exists(combined_file))
  stopifnot(file.exists(millimeter_file))
  stopifnot(file.exists(breakpoint_file))
  
  # ----------------------------------------------------------
  # Read input
  # ----------------------------------------------------------
  phenotype_raw  <- read_sur_csv(phenotype_file)
  genotype_raw   <- read_sur_csv(genotype_file)
  combined_raw   <- read_sur_csv(combined_file)
  millimeter_raw <- read_sur_csv(millimeter_file)
  breakpoint_raw <- read_sur_csv(breakpoint_file)
  
  # ----------------------------------------------------------
  # Antibiotic order
  # Prefer global ANTIBIOTICS if already sourced
  # ----------------------------------------------------------
  if (exists("ANTIBIOTICS", inherits = TRUE)) {
    antibiotics <- get("ANTIBIOTICS", inherits = TRUE)
  } else {
    g_sample_col <- detect_sample_column(genotype_raw)
    p_sample_col <- detect_sample_column(phenotype_raw)
    antibiotics <- intersect(
      normalize_antibiotic_name(setdiff(names(genotype_raw), g_sample_col)),
      normalize_antibiotic_name(setdiff(names(phenotype_raw), p_sample_col))
    )
  }
  antibiotics <- normalize_antibiotic_name(antibiotics)
  
  # ----------------------------------------------------------
  # Standardize wide SUR tables
  # ----------------------------------------------------------
  phenotype <- standardize_sur_wide(phenotype_raw, antibiotics = antibiotics)
  genotype  <- standardize_sur_wide(genotype_raw,  antibiotics = antibiotics)
  
  combined_sample_col <- detect_sample_column(combined_raw)
  combined <- combined_raw |>
    rename(sample = all_of(combined_sample_col)) |>
    mutate(sample = normalize_sample_id(sample))
  
  # ----------------------------------------------------------
  # Standardize millimeters and breakpoints
  # ----------------------------------------------------------
  millimeters_long <- parse_mm_long(millimeter_raw, antibiotics = antibiotics)
  breakpoints_long <- parse_breakpoints_long(breakpoint_raw, antibiotics = antibiotics)
  
  # ----------------------------------------------------------
  # Basic checks
  # ----------------------------------------------------------
  phenotype_samples <- sort(unique(phenotype$sample))
  genotype_samples  <- sort(unique(genotype$sample))
  shared_samples    <- intersect(phenotype_samples, genotype_samples)
  
  sample_check <- tibble(
    dataset = c("phenotype", "genotype"),
    n_samples = c(length(phenotype_samples), length(genotype_samples))
  )
  
  phenotype_only <- setdiff(phenotype_samples, genotype_samples)
  genotype_only  <- setdiff(genotype_samples, phenotype_samples)
  
  antibiotic_check <- tibble(
    antibiotic = antibiotics,
    in_phenotype = antibiotics %in% names(phenotype),
    in_genotype  = antibiotics %in% names(genotype),
    in_mm        = antibiotics %in% unique(millimeters_long$antibiotic),
    in_bp        = antibiotics %in% unique(breakpoints_long$antibiotic)
  )
  
  # ----------------------------------------------------------
  # Long comparison
  # ----------------------------------------------------------
  phenotype_long <- phenotype |>
    filter(sample %in% shared_samples) |>
    arrange(sample) |>
    pivot_longer(
      cols = all_of(antibiotics),
      names_to = "antibiotic",
      values_to = "phenotype_SUR"
    )
  
  genotype_long <- genotype |>
    filter(sample %in% shared_samples) |>
    arrange(sample) |>
    pivot_longer(
      cols = all_of(antibiotics),
      names_to = "antibiotic",
      values_to = "genotype_SUR"
    )
  
  comparison_long <- genotype_long |>
    inner_join(
      phenotype_long,
      by = c("sample", "antibiotic")
    ) |>
    mutate(
      agree = genotype_SUR == phenotype_SUR,
      transition = paste(genotype_SUR, phenotype_SUR, sep = " -> ")
    )
  
  # ----------------------------------------------------------
  # Add millimeters and breakpoint-distance info
  # ----------------------------------------------------------
  comparison_with_mm <- comparison_long |>
    left_join(millimeters_long, by = c("sample", "antibiotic")) |>
    left_join(breakpoints_long, by = "antibiotic") |>
    add_breakpoint_distance_features()
  
  differences_long <- comparison_with_mm |>
    filter(!agree)
  
  # ----------------------------------------------------------
  # Core summaries
  # ----------------------------------------------------------
  overall_summary <- comparison_with_mm |>
    summarise(
      n_samples = n_distinct(sample),
      n_antibiotics = n_distinct(antibiotic),
      n_total_calls = n(),
      n_agree = sum(agree, na.rm = TRUE),
      n_differ = sum(!agree, na.rm = TRUE),
      agreement = n_agree / n_total_calls,
      n_with_mm = sum(!is.na(mm)),
      n_with_breakpoints = sum(!is.na(s_min_mm) & !is.na(r_max_mm)),
      n_with_distance = sum(!is.na(nearest_breakpoint_distance)),
      median_distance_to_breakpoint = median(nearest_breakpoint_distance, na.rm = TRUE)
    )
  
  summary_by_antibiotic <- comparison_with_mm |>
    group_by(antibiotic) |>
    summarise(
      n = n(),
      n_agree = sum(agree, na.rm = TRUE),
      n_differ = sum(!agree, na.rm = TRUE),
      agreement = n_agree / n,
      median_mm = median(mm, na.rm = TRUE),
      median_distance_to_breakpoint = median(nearest_breakpoint_distance, na.rm = TRUE),
      median_distance_disagree = median(nearest_breakpoint_distance[!agree], na.rm = TRUE),
      median_distance_agree = median(nearest_breakpoint_distance[agree], na.rm = TRUE),
      frac_within_1mm = mean(is_near_breakpoint_1mm, na.rm = TRUE),
      frac_within_2mm = mean(is_near_breakpoint_2mm, na.rm = TRUE),
      frac_within_3mm = mean(is_near_breakpoint_3mm, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(n_differ), antibiotic)
  
  summary_by_transition <- differences_long |>
    count(transition, sort = TRUE, name = "n")
  
  summary_by_antibiotic_transition <- differences_long |>
    count(antibiotic, transition, sort = TRUE, name = "n") |>
    arrange(antibiotic, desc(n))
  
  summary_by_sample <- differences_long |>
    count(sample, sort = TRUE, name = "n_differences")
  
  # ----------------------------------------------------------
  # Breakpoint-distance analysis
  # Main question:
  #   are disagreements enriched near breakpoints?
  # ----------------------------------------------------------
  breakpoint_distance_summary <- comparison_with_mm |>
    summarise(
      n = n(),
      n_with_distance = sum(!is.na(nearest_breakpoint_distance)),
      median_distance_all = median(nearest_breakpoint_distance, na.rm = TRUE),
      median_distance_agree = median(nearest_breakpoint_distance[agree], na.rm = TRUE),
      median_distance_disagree = median(nearest_breakpoint_distance[!agree], na.rm = TRUE),
      mean_distance_all = mean(nearest_breakpoint_distance, na.rm = TRUE),
      mean_distance_agree = mean(nearest_breakpoint_distance[agree], na.rm = TRUE),
      mean_distance_disagree = mean(nearest_breakpoint_distance[!agree], na.rm = TRUE),
      frac_disagree_all = mean(!agree, na.rm = TRUE),
      frac_disagree_within_1mm = mean(!agree[is_near_breakpoint_1mm], na.rm = TRUE),
      frac_disagree_within_2mm = mean(!agree[is_near_breakpoint_2mm], na.rm = TRUE),
      frac_disagree_within_3mm = mean(!agree[is_near_breakpoint_3mm], na.rm = TRUE),
      frac_disagree_beyond_3mm = mean(!agree[!is.na(nearest_breakpoint_distance) & nearest_breakpoint_distance > 3], na.rm = TRUE)
    )
  
  breakpoint_distance_by_antibiotic <- comparison_with_mm |>
    group_by(antibiotic) |>
    summarise(
      n = n(),
      n_with_distance = sum(!is.na(nearest_breakpoint_distance)),
      median_distance_all = median(nearest_breakpoint_distance, na.rm = TRUE),
      median_distance_agree = median(nearest_breakpoint_distance[agree], na.rm = TRUE),
      median_distance_disagree = median(nearest_breakpoint_distance[!agree], na.rm = TRUE),
      frac_disagree_all = mean(!agree, na.rm = TRUE),
      frac_disagree_within_1mm = mean(!agree[is_near_breakpoint_1mm], na.rm = TRUE),
      frac_disagree_within_2mm = mean(!agree[is_near_breakpoint_2mm], na.rm = TRUE),
      frac_disagree_within_3mm = mean(!agree[is_near_breakpoint_3mm], na.rm = TRUE),
      frac_disagree_beyond_3mm = mean(!agree[!is.na(nearest_breakpoint_distance) & nearest_breakpoint_distance > 3], na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(frac_disagree_within_2mm), desc(frac_disagree_all), antibiotic)
  
  breakpoint_distance_bins <- comparison_with_mm |>
    mutate(
      distance_bin = case_when(
        is.na(nearest_breakpoint_distance) ~ NA_character_,
        nearest_breakpoint_distance <= 0 ~ "0 mm",
        nearest_breakpoint_distance <= 1 ~ "1 mm",
        nearest_breakpoint_distance <= 2 ~ "2 mm",
        nearest_breakpoint_distance <= 3 ~ "3 mm",
        nearest_breakpoint_distance <= 5 ~ "4-5 mm",
        TRUE ~ ">5 mm"
      )
    ) |>
    filter(!is.na(distance_bin)) |>
    group_by(distance_bin) |>
    summarise(
      n = n(),
      n_differ = sum(!agree, na.rm = TRUE),
      disagreement = n_differ / n,
      .groups = "drop"
    ) |>
    mutate(
      distance_bin = factor(distance_bin, levels = c("0 mm", "1 mm", "2 mm", "3 mm", "4-5 mm", ">5 mm"))
    ) |>
    arrange(distance_bin)
  
  breakpoint_distance_bins_by_antibiotic <- comparison_with_mm |>
    mutate(
      distance_bin = case_when(
        is.na(nearest_breakpoint_distance) ~ NA_character_,
        nearest_breakpoint_distance <= 0 ~ "0 mm",
        nearest_breakpoint_distance <= 1 ~ "1 mm",
        nearest_breakpoint_distance <= 2 ~ "2 mm",
        nearest_breakpoint_distance <= 3 ~ "3 mm",
        nearest_breakpoint_distance <= 5 ~ "4-5 mm",
        TRUE ~ ">5 mm"
      )
    ) |>
    filter(!is.na(distance_bin)) |>
    group_by(antibiotic, distance_bin) |>
    summarise(
      n = n(),
      n_differ = sum(!agree, na.rm = TRUE),
      disagreement = n_differ / n,
      .groups = "drop"
    ) |>
    mutate(
      distance_bin = factor(distance_bin, levels = c("0 mm", "1 mm", "2 mm", "3 mm", "4-5 mm", ">5 mm"))
    ) |>
    arrange(antibiotic, distance_bin)
  
  # ----------------------------------------------------------
  # Phenotype U analysis
  # Useful because U includes I/ATU and should often sit near limits
  # ----------------------------------------------------------
  phenotype_u_breakpoint_summary <- comparison_with_mm |>
    group_by(phenotype_SUR) |>
    summarise(
      n = n(),
      median_mm = median(mm, na.rm = TRUE),
      median_distance_to_breakpoint = median(nearest_breakpoint_distance, na.rm = TRUE),
      frac_within_1mm = mean(is_near_breakpoint_1mm, na.rm = TRUE),
      frac_within_2mm = mean(is_near_breakpoint_2mm, na.rm = TRUE),
      frac_disagree = mean(!agree, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(match(phenotype_SUR, c("S", "U", "R")))
  
  phenotype_u_by_antibiotic <- comparison_with_mm |>
    group_by(antibiotic, phenotype_SUR) |>
    summarise(
      n = n(),
      median_distance_to_breakpoint = median(nearest_breakpoint_distance, na.rm = TRUE),
      frac_within_1mm = mean(is_near_breakpoint_1mm, na.rm = TRUE),
      frac_within_2mm = mean(is_near_breakpoint_2mm, na.rm = TRUE),
      frac_disagree = mean(!agree, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(antibiotic, match(phenotype_SUR, c("S", "U", "R")))
  
  # ----------------------------------------------------------
  # Confusion-style tables
  # ----------------------------------------------------------
  sur_levels <- c("S", "U", "R")
  
  confusion_by_antibiotic <- comparison_with_mm |>
    mutate(
      genotype_SUR  = factor(genotype_SUR, levels = sur_levels),
      phenotype_SUR = factor(phenotype_SUR, levels = sur_levels)
    ) |>
    count(antibiotic, genotype_SUR, phenotype_SUR, name = "n") |>
    complete(
      antibiotic,
      genotype_SUR = sur_levels,
      phenotype_SUR = sur_levels,
      fill = list(n = 0)
    ) |>
    arrange(antibiotic, genotype_SUR, phenotype_SUR)
  
  # ----------------------------------------------------------
  # Wide differences table
  # ----------------------------------------------------------
  comparison_wide <- comparison_with_mm |>
    mutate(value = if_else(agree, "", paste0(genotype_SUR, "|", phenotype_SUR))) |>
    select(sample, antibiotic, value) |>
    pivot_wider(
      names_from = antibiotic,
      values_from = value
    ) |>
    arrange(sample)
  
  # ----------------------------------------------------------
  # Add combined genotype details for debugging
  # ----------------------------------------------------------
  differences_with_combined <- differences_long |>
    left_join(combined, by = "sample")
  
  # ----------------------------------------------------------
  # Write outputs
  # ----------------------------------------------------------
  write_semicolon_csv(sample_check, file.path(out_dir, "SUR_compare_sample_check.csv"))
  write_semicolon_csv(
    tibble(sample = phenotype_only),
    file.path(out_dir, "SUR_compare_samples_only_in_phenotype.csv")
  )
  write_semicolon_csv(
    tibble(sample = genotype_only),
    file.path(out_dir, "SUR_compare_samples_only_in_genotype.csv")
  )
  write_semicolon_csv(antibiotic_check, file.path(out_dir, "SUR_compare_antibiotic_check.csv"))
  
  write_semicolon_csv(overall_summary, file.path(out_dir, "SUR_compare_overall_summary.csv"))
  write_semicolon_csv(summary_by_antibiotic, file.path(out_dir, "SUR_compare_summary_by_antibiotic.csv"))
  write_semicolon_csv(summary_by_transition, file.path(out_dir, "SUR_compare_summary_by_transition.csv"))
  write_semicolon_csv(summary_by_antibiotic_transition, file.path(out_dir, "SUR_compare_summary_by_antibiotic_transition.csv"))
  write_semicolon_csv(summary_by_sample, file.path(out_dir, "SUR_compare_summary_by_sample.csv"))
  
  write_semicolon_csv(comparison_with_mm, file.path(out_dir, "SUR_compare_long_with_mm.csv"))
  write_semicolon_csv(differences_long, file.path(out_dir, "SUR_compare_differences_long_with_mm.csv"))
  write_semicolon_csv(confusion_by_antibiotic, file.path(out_dir, "SUR_compare_confusion_by_antibiotic.csv"))
  write_semicolon_csv(comparison_wide, file.path(out_dir, "SUR_compare_differences_wide.csv"))
  write_semicolon_csv(differences_with_combined, file.path(out_dir, "SUR_compare_differences_with_combined.csv"))
  
  write_semicolon_csv(breakpoint_distance_summary, file.path(out_dir, "SUR_compare_breakpoint_distance_summary.csv"))
  write_semicolon_csv(breakpoint_distance_by_antibiotic, file.path(out_dir, "SUR_compare_breakpoint_distance_by_antibiotic.csv"))
  write_semicolon_csv(breakpoint_distance_bins, file.path(out_dir, "SUR_compare_breakpoint_distance_bins.csv"))
  write_semicolon_csv(breakpoint_distance_bins_by_antibiotic, file.path(out_dir, "SUR_compare_breakpoint_distance_bins_by_antibiotic.csv"))
  write_semicolon_csv(phenotype_u_breakpoint_summary, file.path(out_dir, "SUR_compare_phenotype_u_breakpoint_summary.csv"))
  write_semicolon_csv(phenotype_u_by_antibiotic, file.path(out_dir, "SUR_compare_phenotype_u_by_antibiotic.csv"))
  
  invisible(list(
    phenotype = phenotype,
    genotype = genotype,
    combined = combined,
    millimeters_long = millimeters_long,
    breakpoints_long = breakpoints_long,
    comparison_with_mm = comparison_with_mm,
    differences_long = differences_long,
    overall_summary = overall_summary,
    summary_by_antibiotic = summary_by_antibiotic,
    summary_by_transition = summary_by_transition,
    summary_by_antibiotic_transition = summary_by_antibiotic_transition,
    summary_by_sample = summary_by_sample,
    confusion_by_antibiotic = confusion_by_antibiotic,
    breakpoint_distance_summary = breakpoint_distance_summary,
    breakpoint_distance_by_antibiotic = breakpoint_distance_by_antibiotic,
    breakpoint_distance_bins = breakpoint_distance_bins,
    phenotype_u_breakpoint_summary = phenotype_u_breakpoint_summary
  ))
}

RUN <- function()
{
  res <- comparePhenotypeAndGenotypeSUR()
  
  res$overall_summary
  res$summary_by_antibiotic
  res$breakpoint_distance_summary
  res$breakpoint_distance_by_antibiotic
  res$phenotype_u_breakpoint_summary
}

RUN()