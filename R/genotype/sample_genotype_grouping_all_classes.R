# ============================================================================
# sampleGenotypeGrouping_all_classes.R
# ----------------------------------------------------------------------------
# Purpose
#   Build per-sample genotype grouping summaries for:
#     1) beta-lactamases
#     2) quinolone-associated genes
#     3) aminoglycoside-associated genes
#
# Inputs
#   - resolved_amr_gene_calls_long_core.csv
#   - allSamples() function returning all samples
#
# Output
#   - combined per-sample summary table
#   - class-specific tables and group maps
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
})

# ----------------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------------

source(file = "common.R")
source(file = "genotype/genotypeCommon.R")

indir <- file.path(processedRootRassembly, "genotype")
outdir <- file.path(processedRootRassembly, "genotype")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

resolved_file <- file.path(indir, "resolved_amr_gene_calls_long_core.csv")
stopifnot(file.exists(resolved_file))

# ----------------------------------------------------------------------------
# HELPERS
# ----------------------------------------------------------------------------

collapse_unique <- function(x, sep = "; ") {
  x <- as.character(x)
  x <- unique(x)
  x <- x[!is.na(x)]
  x <- x[x != ""]
  if (length(x) == 0) return("")
  paste(sort(x), collapse = sep)
}

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA")] <- NA_character_
  
  out <- x |>
    stringr::str_remove("\\.tsv$") |>
    stringr::str_remove("_amrfinderplus$") |>
    stringr::str_remove("_amrfinder$") |>
    stringr::str_extract("[0-9]+")
  
  out_num <- suppressWarnings(as.integer(out))
  
  out_chr <- rep(NA_character_, length(out_num))
  ok <- !is.na(out_num)
  out_chr[ok] <- sprintf("%03d", out_num[ok])
  
  out_chr
}

sample_order_key <- function(x) {
  suppressWarnings(as.integer(normalize_sample_id(x)))
}

make_group_labels <- function(values, prefix) {
  vals <- as.character(values)
  vals[is.na(vals)] <- ""
  
  uniq <- sort(unique(vals))
  map <- tibble(
    value = uniq,
    group = sprintf("%s%02d", prefix, seq_along(uniq))
  )
  
  groups <- tibble(value = vals) %>%
    left_join(map, by = "value") %>%
    pull(group)
  
  list(map = map, groups = groups)
}

# ----------------------------------------------------------------------------
# GENE CLASSIFICATION HELPERS
# ----------------------------------------------------------------------------

is_beta_lactamase <- function(gene) {
  gene <- as.character(gene)
  stringr::str_detect(gene, regex("^(bla|ampC)", ignore_case = TRUE))
}

infer_quinolone_family <- function(gene) {
  gene <- as.character(gene)
  
  dplyr::case_when(
    stringr::str_detect(gene, regex("^qnrA", ignore_case = TRUE)) ~ "QnrA",
    stringr::str_detect(gene, regex("^qnrB", ignore_case = TRUE)) ~ "QnrB",
    stringr::str_detect(gene, regex("^qnrS", ignore_case = TRUE)) ~ "QnrS",
    stringr::str_detect(gene, regex("^qnrC", ignore_case = TRUE)) ~ "QnrC",
    stringr::str_detect(gene, regex("^qnrD", ignore_case = TRUE)) ~ "QnrD",
    stringr::str_detect(gene, regex("^qepA", ignore_case = TRUE)) ~ "QepA",
    stringr::str_detect(gene, regex("^oqxA$", ignore_case = TRUE)) ~ "OqxA",
    stringr::str_detect(gene, regex("^oqxB$", ignore_case = TRUE)) ~ "OqxB",
    stringr::str_detect(gene, regex("^(gyrA|gyrB|parC|parE)", ignore_case = TRUE)) ~ "QRDR",
    TRUE ~ NA_character_
  )
}

infer_aminoglycoside_family <- function(gene) {
  gene <- as.character(gene)
  
  dplyr::case_when(
    stringr::str_detect(gene, regex("^aac\\(6", ignore_case = TRUE)) ~ "AAC(6)",
    stringr::str_detect(gene, regex("^aac\\(3", ignore_case = TRUE)) ~ "AAC(3)",
    stringr::str_detect(gene, regex("^aac\\(2", ignore_case = TRUE)) ~ "AAC(2)",
    stringr::str_detect(gene, regex("^aph\\(3''\\)", ignore_case = TRUE)) ~ "APH(3'')",
    stringr::str_detect(gene, regex("^aph\\(3'\\)", ignore_case = TRUE)) ~ "APH(3')",
    stringr::str_detect(gene, regex("^aph\\(4\\)", ignore_case = TRUE)) ~ "APH(4)",
    stringr::str_detect(gene, regex("^aph\\(6\\)", ignore_case = TRUE)) ~ "APH(6)",
    stringr::str_detect(gene, regex("^ant\\(", ignore_case = TRUE)) ~ "ANT",
    stringr::str_detect(gene, regex("^aadA", ignore_case = TRUE)) ~ "AAD",
    stringr::str_detect(gene, regex("^armA$", ignore_case = TRUE)) ~ "ArmA",
    stringr::str_detect(gene, regex("^rmt", ignore_case = TRUE)) ~ "Rmt",
    TRUE ~ NA_character_
  )
}

# ----------------------------------------------------------------------------
# READ RESOLVED LONG TABLE
# ----------------------------------------------------------------------------

resolved_long <- read.csv2(
  resolved_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(all(c("sample_id", "Gene_symbol") %in% names(resolved_long)))

if (!"Enzyme_family" %in% names(resolved_long)) {
  resolved_long$Enzyme_family <- NA_character_
}
if (!"Functional_group" %in% names(resolved_long)) {
  resolved_long$Functional_group <- NA_character_
}

resolved_long <- resolved_long %>%
  mutate(
    sample_id = normalize_sample_id(sample_id),
    Gene_symbol = as.character(Gene_symbol),
    Enzyme_family = as.character(Enzyme_family),
    Functional_group = as.character(Functional_group)
  )

# ----------------------------------------------------------------------------
# ALL SAMPLES
# ----------------------------------------------------------------------------

all_samples <- tibble(sample = normalize_sample_id(allSamples())) %>%
  distinct() %>%
  filter(!is.na(sample))

# ----------------------------------------------------------------------------
# 1) BETA-LACTAMASE SUMMARY
# ----------------------------------------------------------------------------

beta_annot <- resolved_long %>%
  filter(is_beta_lactamase(Gene_symbol))

beta_summary_per_sample <- all_samples %>%
  left_join(
    beta_annot %>%
      rename(sample = sample_id) %>%
      group_by(sample) %>%
      summarise(
        Beta_genes = collapse_unique(Gene_symbol),
        Functional_groups = collapse_unique(Functional_group),
        Enzyme_families = collapse_unique(Enzyme_family),
        .groups = "drop"
      ),
    by = "sample"
  ) %>%
  mutate(
    Beta_genes = if_else(is.na(Beta_genes), "", Beta_genes),
    Functional_groups = if_else(is.na(Functional_groups), "", Functional_groups),
    Enzyme_families = if_else(is.na(Enzyme_families), "", Enzyme_families)
  )

beta_group_obj <- make_group_labels(
  beta_summary_per_sample$Functional_groups,
  prefix = "FG"
)

beta_summary_per_sample <- beta_summary_per_sample %>%
  mutate(
    samplegroup_functional = beta_group_obj$groups
  ) %>%
  select(sample, samplegroup_functional, Functional_groups, Beta_genes, Enzyme_families) %>%
  arrange(samplegroup_functional, sample_order_key(sample), sample)

# ----------------------------------------------------------------------------
# 2) QUINOLONE SUMMARY
# ----------------------------------------------------------------------------

quinolone_annot <- resolved_long %>%
  transmute(
    sample_id,
    Gene_symbol,
    Quinolone_family = infer_quinolone_family(Gene_symbol)
  ) %>%
  filter(!is.na(Quinolone_family))

quinolone_summary_per_sample <- all_samples %>%
  left_join(
    quinolone_annot %>%
      rename(sample = sample_id) %>%
      group_by(sample) %>%
      summarise(
        Quinolone_genes = collapse_unique(Gene_symbol),
        Quinolone_families = collapse_unique(Quinolone_family),
        .groups = "drop"
      ),
    by = "sample"
  ) %>%
  mutate(
    Quinolone_genes = if_else(is.na(Quinolone_genes), "", Quinolone_genes),
    Quinolone_families = if_else(is.na(Quinolone_families), "", Quinolone_families)
  )

quinolone_group_obj <- make_group_labels(
  quinolone_summary_per_sample$Quinolone_families,
  prefix = "QG"
)

quinolone_summary_per_sample <- quinolone_summary_per_sample %>%
  mutate(
    samplegroup_quinolone = quinolone_group_obj$groups
  ) %>%
  select(sample, samplegroup_quinolone, Quinolone_genes, Quinolone_families) %>%
  arrange(samplegroup_quinolone, sample_order_key(sample), sample)

# ----------------------------------------------------------------------------
# 3) AMINOGLYCOSIDE SUMMARY
# ----------------------------------------------------------------------------

aminoglycoside_annot <- resolved_long %>%
  transmute(
    sample_id,
    Gene_symbol,
    Aminoglycoside_family = infer_aminoglycoside_family(Gene_symbol)
  ) %>%
  filter(!is.na(Aminoglycoside_family))

aminoglycoside_summary_per_sample <- all_samples %>%
  left_join(
    aminoglycoside_annot %>%
      rename(sample = sample_id) %>%
      group_by(sample) %>%
      summarise(
        Aminoglycoside_genes = collapse_unique(Gene_symbol),
        Aminoglycoside_families = collapse_unique(Aminoglycoside_family),
        .groups = "drop"
      ),
    by = "sample"
  ) %>%
  mutate(
    Aminoglycoside_genes = if_else(is.na(Aminoglycoside_genes), "", Aminoglycoside_genes),
    Aminoglycoside_families = if_else(is.na(Aminoglycoside_families), "", Aminoglycoside_families)
  )

aminoglycoside_group_obj <- make_group_labels(
  aminoglycoside_summary_per_sample$Aminoglycoside_families,
  prefix = "AG"
)

aminoglycoside_summary_per_sample <- aminoglycoside_summary_per_sample %>%
  mutate(
    samplegroup_aminoglycoside = aminoglycoside_group_obj$groups
  ) %>%
  select(sample, samplegroup_aminoglycoside, Aminoglycoside_genes, Aminoglycoside_families) %>%
  arrange(samplegroup_aminoglycoside, sample_order_key(sample), sample)

# ----------------------------------------------------------------------------
# 4) FINAL COMBINED SAMPLE SUMMARY
# ----------------------------------------------------------------------------

sample_genotype_summary <- beta_summary_per_sample %>%
  left_join(
    quinolone_summary_per_sample,
    by = "sample"
  ) %>%
  left_join(
    aminoglycoside_summary_per_sample,
    by = "sample"
  ) %>%
  mutate(
    samplegroup_functional = if_else(is.na(samplegroup_functional), "", samplegroup_functional),
    samplegroup_quinolone = if_else(is.na(samplegroup_quinolone), "", samplegroup_quinolone),
    samplegroup_aminoglycoside = if_else(is.na(samplegroup_aminoglycoside), "", samplegroup_aminoglycoside),
    
    Functional_groups = if_else(is.na(Functional_groups), "", Functional_groups),
    Beta_genes = if_else(is.na(Beta_genes), "", Beta_genes),
    Enzyme_families = if_else(is.na(Enzyme_families), "", Enzyme_families),
    
    Quinolone_genes = if_else(is.na(Quinolone_genes), "", Quinolone_genes),
    Quinolone_families = if_else(is.na(Quinolone_families), "", Quinolone_families),
    
    Aminoglycoside_genes = if_else(is.na(Aminoglycoside_genes), "", Aminoglycoside_genes),
    Aminoglycoside_families = if_else(is.na(Aminoglycoside_families), "", Aminoglycoside_families)
  ) %>%
  arrange(
    samplegroup_functional,
    samplegroup_quinolone,
    samplegroup_aminoglycoside,
    sample_order_key(sample),
    sample
  )

# ----------------------------------------------------------------------------
# 5) SIMPLE SUMMARY TABLE
# ----------------------------------------------------------------------------

sample_genotype_summary_simple <- sample_genotype_summary %>%
  select(
    sample,
    samplegroup_functional,
    Functional_groups,
    Beta_genes,
    Enzyme_families,
    samplegroup_quinolone,
    Quinolone_families,
    samplegroup_aminoglycoside,
    Aminoglycoside_families
  ) %>%
  arrange(
    samplegroup_functional,
    samplegroup_quinolone,
    samplegroup_aminoglycoside,
    sample_order_key(sample),
    sample
  )

# ----------------------------------------------------------------------------
# CHECKS
# ----------------------------------------------------------------------------

stopifnot(nrow(sample_genotype_summary) == nrow(all_samples))
stopifnot(length(unique(sample_genotype_summary$sample)) == nrow(all_samples))

# ----------------------------------------------------------------------------
# WRITE OUTPUT
# ----------------------------------------------------------------------------

write.csv2(
  beta_summary_per_sample,
  file.path(outdir, "sample_beta_grouping.csv"),
  row.names = FALSE
)

write.csv2(
  quinolone_summary_per_sample,
  file.path(outdir, "sample_quinolone_grouping.csv"),
  row.names = FALSE
)

write.csv2(
  aminoglycoside_summary_per_sample,
  file.path(outdir, "sample_aminoglycoside_grouping.csv"),
  row.names = FALSE
)

write.csv2(
  sample_genotype_summary,
  file.path(outdir, "sample_genotype_grouping_all_classes.csv"),
  row.names = FALSE
)

write.csv2(
  sample_genotype_summary_simple,
  file.path(outdir, "sample_genotype_grouping_simple.csv"),
  row.names = FALSE
)

write.csv2(
  beta_group_obj$map,
  file.path(outdir, "samplegroup_functional_map.csv"),
  row.names = FALSE
)

write.csv2(
  quinolone_group_obj$map,
  file.path(outdir, "samplegroup_quinolone_map.csv"),
  row.names = FALSE
)

write.csv2(
  aminoglycoside_group_obj$map,
  file.path(outdir, "samplegroup_aminoglycoside_map.csv"),
  row.names = FALSE
)

message("Done.")
message("Wrote:")
message(" - sample_beta_grouping.csv")
message(" - sample_quinolone_grouping.csv")
message(" - sample_aminoglycoside_grouping.csv")
message(" - sample_genotype_grouping_all_classes.csv")
message(" - sample_genotype_grouping_simple.csv")
message(" - samplegroup_functional_map.csv")
message(" - samplegroup_quinolone_map.csv")
message(" - samplegroup_aminoglycoside_map.csv")