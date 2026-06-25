# ============================================================
# Plot mismatch diagnostics for phenotype vs genotype SUR
# Only considers hard mismatches: S -> R and R -> S
#
# Expected input files in:
#   paste(processedRootRassembly, "phenotype", sep = "/")
#
# Reads:
#   - SUR_compare_differences_long_with_mm.csv
#   - SUR_compare_summary_by_antibiotic_transition.csv
#   - genotype_SUR_combined_per_sample.csv (optional)
#
# Writes plots to:
#   paste(processedRootRassembly, "phenotype", "mismatch_plots", sep = "/")
# ============================================================

source("genotype/genotypeCommon.R")

plotHardMismatchSUR <- function(
    input_dir = paste(processedRootRassembly, "phenotype", sep = "/"),
    out_dir   = file.path(input_dir, "mismatch_plots")
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(ggplot2)
    library(forcats)
    library(scales)
  })
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  read_semicolon_csv <- function(path) {
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
  
  normalize_sample_id <- function(x) {
    x <- as.character(x)
    x <- stringr::str_trim(x)
    x <- stringr::str_replace_all(x, "\\.0$", "")
    x <- stringr::str_replace_all(x, "^sample[_\\- ]*", "")
    x <- stringr::str_replace_all(x, "[^A-Za-z0-9]", "")
    is_num <- grepl("^[0-9]+$", x)
    x[is_num] <- as.character(as.integer(x[is_num]))
    toupper(x)
  }
  
  normalize_antibiotic <- function(x) {
    x |>
      as.character() |>
      stringr::str_trim() |>
      toupper()
  }
  
  safe_frac <- function(num, den) {
    dplyr::if_else(den > 0, num / den, NA_real_)
  }
  
  save_plot <- function(plot, filename, width = 8, height = 5) {
    ggplot2::ggsave(
      filename = file.path(out_dir, filename),
      plot = plot,
      width = width,
      height = height,
      dpi = 600
    )
  }
  
  transition_palette <- c("S -> R" = "#d73027", "R -> S" = "#4575b4")
  
  diff_file     <- file.path(input_dir, "SUR_compare_differences_long_with_mm.csv")
  by_ab_tr_file <- file.path(input_dir, "SUR_compare_summary_by_antibiotic_transition.csv")
  combined_file <- file.path(input_dir, "genotype_SUR_combined_per_sample.csv")
  
  if (!file.exists(diff_file)) stop("Missing file: ", diff_file)
  if (!file.exists(by_ab_tr_file)) stop("Missing file: ", by_ab_tr_file)
  
  differences_long <- read_semicolon_csv(diff_file)
  summary_by_ab_transition <- read_semicolon_csv(by_ab_tr_file)
  
  combined <- NULL
  if (file.exists(combined_file)) {
    combined <- read_semicolon_csv(combined_file)
    if (ncol(combined) >= 1) {
      names(combined)[1] <- "sample"
      combined <- combined |>
        mutate(sample = normalize_sample_id(sample))
    }
  }
  
  req_cols <- c("sample", "antibiotic", "genotype_SUR", "phenotype_SUR", "transition")
  missing_req <- setdiff(req_cols, names(differences_long))
  if (length(missing_req) > 0) {
    stop("SUR_compare_differences_long_with_mm.csv is missing required columns: ",
         paste(missing_req, collapse = ", "))
  }
  
  differences_long <- differences_long |>
    mutate(
      sample = normalize_sample_id(sample),
      antibiotic = normalize_antibiotic(antibiotic),
      genotype_SUR = toupper(trimws(as.character(genotype_SUR))),
      phenotype_SUR = toupper(trimws(as.character(phenotype_SUR))),
      transition = paste(genotype_SUR, phenotype_SUR, sep = " -> ")
    )
  
  hard_mismatches <- differences_long |>
    filter(transition %in% c("S -> R", "R -> S"))
  
  if (nrow(hard_mismatches) == 0) {
    stop("No S -> R or R -> S mismatches found.")
  }
  
  if (exists("ANTIBIOTICS", inherits = TRUE)) {
    ab_levels <- normalize_antibiotic(get("ANTIBIOTICS", inherits = TRUE))
    ab_levels <- intersect(ab_levels, unique(hard_mismatches$antibiotic))
    ab_levels <- c(ab_levels, setdiff(sort(unique(hard_mismatches$antibiotic)), ab_levels))
  } else {
    ab_levels <- sort(unique(hard_mismatches$antibiotic))
  }
  
  hard_mismatches <- hard_mismatches |>
    mutate(
      antibiotic = factor(antibiotic, levels = ab_levels),
      transition = factor(transition, levels = c("S -> R", "R -> S"))
    )
  
  summary_by_ab_transition <- summary_by_ab_transition |>
    mutate(
      antibiotic = normalize_antibiotic(antibiotic),
      transition = trimws(as.character(transition))
    ) |>
    filter(transition %in% c("S -> R", "R -> S")) |>
    mutate(
      antibiotic = factor(antibiotic, levels = ab_levels),
      transition = factor(transition, levels = c("S -> R", "R -> S"))
    )
  
  hard_by_antibiotic <- hard_mismatches |>
    count(antibiotic, transition, name = "n") |>
    tidyr::complete(
      antibiotic = factor(ab_levels, levels = ab_levels),
      transition = factor(c("S -> R", "R -> S"), levels = c("S -> R", "R -> S")),
      fill = list(n = 0)
    ) |>
    group_by(antibiotic) |>
    mutate(total = sum(n), frac = safe_frac(n, total)) |>
    ungroup()
  
  hard_by_sample <- hard_mismatches |>
    count(sample, transition, name = "n") |>
    tidyr::complete(
      sample,
      transition = factor(c("S -> R", "R -> S"), levels = c("S -> R", "R -> S")),
      fill = list(n = 0)
    ) |>
    group_by(sample) |>
    mutate(total = sum(n)) |>
    ungroup()
  
  sample_order <- hard_by_sample |>
    distinct(sample, total) |>
    arrange(desc(total), sample) |>
    pull(sample)
  
  hard_by_sample <- hard_by_sample |>
    mutate(sample = factor(sample, levels = sample_order))
  
  mismatch_matrix <- hard_mismatches |>
    select(sample, antibiotic, transition) |>
    distinct() |>
    mutate(
      sample = factor(sample, levels = sample_order),
      antibiotic = factor(antibiotic, levels = ab_levels)
    )
  
  if ("nearest_breakpoint_distance" %in% names(hard_mismatches)) {
    hard_distance <- hard_mismatches |>
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
      mutate(distance_bin = factor(distance_bin, levels = c("0 mm", "1 mm", "2 mm", "3 mm", "4-5 mm", ">5 mm")))
    
    hard_distance_summary <- hard_distance |>
      count(distance_bin, transition, name = "n") |>
      group_by(distance_bin) |>
      mutate(total = sum(n), frac = safe_frac(n, total)) |>
      ungroup()
    
    hard_distance_by_ab <- hard_distance |>
      count(antibiotic, distance_bin, transition, name = "n") |>
      group_by(antibiotic, distance_bin) |>
      mutate(total = sum(n), frac = safe_frac(n, total)) |>
      ungroup()
  } else {
    hard_distance_summary <- NULL
    hard_distance_by_ab <- NULL
  }
  
  if (!is.null(combined) && "sample" %in% names(combined)) {
    mismatch_with_combined <- hard_mismatches |>
      left_join(combined, by = "sample")
    write_semicolon_csv(mismatch_with_combined,
                        file.path(out_dir, "SUR_compare_hard_mismatches_with_combined.csv"))
  }
  
  p1 <- hard_by_antibiotic |>
    ggplot(aes(x = fct_reorder(as.character(antibiotic), total, .desc = TRUE), y = n, fill = transition)) +
    geom_col(position = "stack") +
    coord_flip() +
    scale_fill_manual(values = transition_palette, drop = FALSE) +
    labs(
      x = NULL,
      y = "Number of hard mismatches",
      fill = "Transition",
      title = "Hard mismatches by antibiotic",
      subtitle = "Only S -> R and R -> S mismatches"
    ) +
    theme_minimal(base_size = 12)
  
  save_plot(p1, "hard_mismatches_by_antibiotic.png", width = 8, height = 5.5)
  
  p2 <- hard_by_antibiotic |>
    ggplot(aes(x = transition, y = antibiotic, fill = n)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n), size = 3) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Count",
      title = "Hard-mismatch transition heatmap"
    ) +
    theme_minimal(base_size = 12)
  
  save_plot(p2, "hard_mismatch_transition_heatmap.png", width = 6, height = 6.5)
  
  top_n <- 30
  sample_totals <- hard_by_sample |>
    distinct(sample, total) |>
    arrange(desc(total), sample)
  
  n_top <- min(top_n, nrow(sample_totals))
  
  top_samples <- sample_totals |>
    slice_head(n = n_top) |>
    pull(sample) |>
    as.character()
  
  p3 <- hard_by_sample |>
    filter(as.character(sample) %in% top_samples) |>
    ggplot(aes(x = sample, y = n, fill = transition)) +
    geom_col(position = "stack") +
    coord_flip() +
    scale_fill_manual(values = transition_palette, drop = FALSE) +
    labs(
      x = NULL,
      y = "Number of hard mismatches",
      fill = "Transition",
      title = "Top samples by hard mismatches",
      subtitle = paste0("Top ", n_top, " samples")
    ) +
    theme_minimal(base_size = 12)
  
  save_plot(p3, "hard_mismatches_top_samples.png", width = 8, height = 7)
  
  p4 <- mismatch_matrix |>
    filter(as.character(sample) %in% top_samples) |>
    ggplot(aes(x = antibiotic, y = sample, fill = transition)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_manual(values = transition_palette, drop = FALSE) +
    labs(
      x = NULL,
      y = NULL,
      fill = "Transition",
      title = "Hard-mismatch map for top samples"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot(p4, "hard_mismatch_map_top_samples.png", width = 8.5, height = 7)
  
  if (!is.null(hard_distance_summary) && nrow(hard_distance_summary) > 0) {
    p5 <- hard_distance_summary |>
      ggplot(aes(x = distance_bin, y = n, fill = transition)) +
      geom_col(position = "stack") +
      scale_fill_manual(values = transition_palette, drop = FALSE) +
      labs(
        x = "Nearest breakpoint distance",
        y = "Number of hard mismatches",
        fill = "Transition",
        title = "Hard mismatches by breakpoint distance"
      ) +
      theme_minimal(base_size = 12)
    
    save_plot(p5, "hard_mismatches_by_breakpoint_distance.png", width = 7.5, height = 5)
    
    p6 <- hard_distance_by_ab |>
      ggplot(aes(x = distance_bin, y = antibiotic, fill = n)) +
      geom_tile(color = "white") +
      labs(
        x = "Nearest breakpoint distance",
        y = NULL,
        fill = "Count",
        title = "Hard mismatches by antibiotic and breakpoint distance"
      ) +
      theme_minimal(base_size = 12)
    
    save_plot(p6, "hard_mismatches_by_antibiotic_and_distance.png", width = 8.5, height = 6.5)
  }
  
  write_semicolon_csv(hard_mismatches,
                      file.path(out_dir, "SUR_compare_hard_mismatches_long.csv"))
  write_semicolon_csv(hard_by_antibiotic,
                      file.path(out_dir, "SUR_compare_hard_mismatches_by_antibiotic.csv"))
  write_semicolon_csv(hard_by_sample,
                      file.path(out_dir, "SUR_compare_hard_mismatches_by_sample.csv"))
  if (!is.null(hard_distance_summary)) {
    write_semicolon_csv(hard_distance_summary,
                        file.path(out_dir, "SUR_compare_hard_mismatches_by_distance.csv"))
  }
  if (!is.null(hard_distance_by_ab)) {
    write_semicolon_csv(hard_distance_by_ab,
                        file.path(out_dir, "SUR_compare_hard_mismatches_by_antibiotic_and_distance.csv"))
  }
  
  invisible(list(
    hard_mismatches = hard_mismatches,
    hard_by_antibiotic = hard_by_antibiotic,
    hard_by_sample = hard_by_sample,
    mismatch_matrix = mismatch_matrix,
    hard_distance_summary = hard_distance_summary,
    hard_distance_by_ab = hard_distance_by_ab,
    out_dir = out_dir
  ))
}

RUN <- function()
{
  plotHardMismatchSUR()
}

RUN()