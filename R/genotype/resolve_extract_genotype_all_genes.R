# ============================================================================
# resolve_extract_genotype_all_genes.R
# ----------------------------------------------------------------------------
# Purpose
#   Build resolved AMR gene calls per sample using:
#     - AMRFinder as primary source
#     - ARIBA as resolver for partial/generic AMRFinder calls
#   Keep only CORE genes.
#
# Outputs
#   1) Long resolved sample-gene table for all AMR genes
#   2) Manuscript table with one row per gene present in any sample
#   3) Beta-lactamase-only subtables for both
#
# Requested conventions
#   - class and subclass come from AMRFinder
#   - Molecular_class and Functional_group are only for beta-lactamases
#   - reference is only for beta-lactamases
#   - tool says whether resolved gene symbol came from AMRFinder or ariba
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(janitor)
})

# ----------------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------------

source(file = "common.R")

outdir <- file.path(processedRootRassembly, "genotype")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

amrf_file  <- file.path(outdir, "amrfinder_amr_long.csv")
ariba_file <- file.path(outdir, "ariba_genes_long.csv")
bldb_file  <- file.path(outdir, "bldb_lookup.csv")
bldb_alias_file <- file.path(outdir, "bldb_aliases.csv")

stopifnot(file.exists(amrf_file))
stopifnot(file.exists(ariba_file))
stopifnot(file.exists(bldb_file))

# ----------------------------------------------------------------------------
# HELPERS
# ----------------------------------------------------------------------------

na_if_empty <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

first_non_missing <- function(x) {
  x <- na_if_empty(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  x[[1]]
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(na_if_empty(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = sep)
}

normalize_sample_id2 <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace(x, "^sample", "")
  x <- str_replace(x, "\\.tsv$", "")
  x <- str_replace(x, "\\.csv$", "")
  x[x == ""] <- NA_character_
  x
}

normalize_gene_symbol <- function(x) {
  x <- na_if_empty(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\s+", "")
  x
}

normalize_gene_for_lookup <- function(x) {
  x |>
    na_if_empty() |>
    replace_na("") |>
    str_replace("^bla", "") |>
    str_replace_all("\\s+", "") |>
    str_to_upper()
}

is_beta_lactamase <- function(x) {
  x <- normalize_gene_symbol(x)
  str_detect(x, regex("^(bla|ampC)", ignore_case = TRUE))
}

is_generic_family_symbol <- function(x) {
  x <- normalize_gene_symbol(x)
  case_when(
    is.na(x) ~ FALSE,
    str_detect(x, regex("^blaTEM$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^blaSHV$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^blaCTX-M$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^blaOXA$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^blaCMY$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^blaDHA$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^ampC$", ignore_case = TRUE)) ~ TRUE,
    TRUE ~ FALSE
  )
}

is_specific_allele_symbol <- function(x) {
  x <- normalize_gene_symbol(x)
  case_when(
    is.na(x) ~ FALSE,
    str_detect(x, regex("^blaCTX-M-[0-9A-Za-z._-]+$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^blaTEM-[0-9A-Za-z._-]+$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^blaSHV-[0-9A-Za-z._-]+$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^blaOXA-[0-9A-Za-z._-]+$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^bla(KPC|NDM|VIM|IMP|GES|CMY|DHA|FOX|ACC|ACT|MIR|LAT|MOX|LAP)-[0-9A-Za-z._-]+$", ignore_case = TRUE)) ~ TRUE,
    TRUE ~ FALSE
  )
}

# ----------------------------------------------------------------------------
# BETA FAMILY INFERENCE
# ----------------------------------------------------------------------------

beta_family_from_symbol <- function(x) {
  x <- normalize_gene_symbol(x)
  
  case_when(
    is.na(x) ~ NA_character_,
    str_detect(x, regex("^blaCTX-M", ignore_case = TRUE)) ~ "CTX-M",
    str_detect(x, regex("^blaTEM", ignore_case = TRUE)) ~ "TEM",
    str_detect(x, regex("^blaSHV", ignore_case = TRUE)) ~ "SHV",
    str_detect(x, regex("^blaOXA", ignore_case = TRUE)) ~ "OXA",
    str_detect(x, regex("^blaNDM", ignore_case = TRUE)) ~ "NDM",
    str_detect(x, regex("^blaKPC", ignore_case = TRUE)) ~ "KPC",
    str_detect(x, regex("^blaVIM", ignore_case = TRUE)) ~ "VIM",
    str_detect(x, regex("^blaIMP", ignore_case = TRUE)) ~ "IMP",
    str_detect(x, regex("^blaLAP", ignore_case = TRUE)) ~ "LAP",
    str_detect(x, regex("^blaCMY", ignore_case = TRUE)) ~ "CMY",
    str_detect(x, regex("^blaDHA", ignore_case = TRUE)) ~ "DHA",
    str_detect(x, regex("^blaFOX", ignore_case = TRUE)) ~ "FOX",
    str_detect(x, regex("^blaACC", ignore_case = TRUE)) ~ "ACC",
    str_detect(x, regex("^blaACT", ignore_case = TRUE)) ~ "ACT",
    str_detect(x, regex("^blaMIR", ignore_case = TRUE)) ~ "MIR",
    str_detect(x, regex("^blaLAT", ignore_case = TRUE)) ~ "LAT",
    str_detect(x, regex("^blaMOX", ignore_case = TRUE)) ~ "MOX",
    str_detect(x, regex("^blaGES", ignore_case = TRUE)) ~ "GES",
    str_detect(x, regex("^blaVEB", ignore_case = TRUE)) ~ "VEB",
    str_detect(x, regex("^blaPER", ignore_case = TRUE)) ~ "PER",
    str_detect(x, regex("^blaEC", ignore_case = TRUE)) ~ "EC",
    str_detect(x, regex("^ampC", ignore_case = TRUE)) ~ "AMPC",
    TRUE ~ NA_character_
  )
}

infer_beta_annotation_from_family <- function(gene_symbol) {
  fam <- beta_family_from_symbol(gene_symbol)
  
  tibble(
    Enzyme_family = case_when(
      fam == "CTX-M" ~ "CTX-M",
      fam == "TEM" ~ "TEM",
      fam == "SHV" ~ "SHV",
      fam == "OXA" ~ "OXA",
      fam == "NDM" ~ "NDM",
      fam == "KPC" ~ "KPC",
      fam == "VIM" ~ "VIM",
      fam == "IMP" ~ "IMP",
      fam == "LAP" ~ "LAP",
      fam == "GES" ~ "GES",
      fam == "VEB" ~ "VEB",
      fam == "PER" ~ "PER",
      fam %in% c("CMY", "DHA", "FOX", "ACC", "ACT", "MIR", "LAT", "MOX", "EC", "AMPC") ~ "AMPC",
      TRUE ~ NA_character_
    ),
    Molecular_class = case_when(
      fam %in% c("CTX-M", "TEM", "SHV", "KPC", "LAP", "GES", "VEB", "PER",
                 "CMY", "DHA", "FOX", "ACC", "ACT", "MIR", "LAT", "MOX") ~ "A",
      fam %in% c("NDM", "VIM", "IMP") ~ "B",
      fam %in% c("EC", "AMPC") ~ "C",
      fam == "OXA" ~ "D",
      TRUE ~ NA_character_
    ),
    Functional_group = case_when(
      fam == "CTX-M" ~ "2be",
      fam == "TEM" ~ "2b",
      fam == "SHV" ~ "2b",
      fam == "KPC" ~ "2f",
      fam == "LAP" ~ "2b",
      fam %in% c("NDM", "VIM", "IMP") ~ "3a",
      fam == "OXA" ~ "2d",
      fam == "GES" ~ "2f",
      fam == "VEB" ~ "2be",
      fam == "PER" ~ "2be",
      fam %in% c("CMY", "DHA", "FOX", "ACC", "ACT", "MIR", "LAT", "MOX", "EC", "AMPC") ~ "1",
      TRUE ~ NA_character_
    ),
    reference = case_when(
      !is.na(fam) ~ "inferred from beta-lactamase family",
      TRUE ~ NA_character_
    )
  )
}

# ----------------------------------------------------------------------------
# READ INPUTS
# ----------------------------------------------------------------------------

amrf <- read.csv2(amrf_file, stringsAsFactors = FALSE, check.names = FALSE) |>
  janitor::clean_names()

ariba <- read.csv2(ariba_file, stringsAsFactors = FALSE, check.names = FALSE) |>
  janitor::clean_names()

bldb <- read.csv2(bldb_file, stringsAsFactors = FALSE, check.names = FALSE) |>
  janitor::clean_names()

bldb_aliases <- if (file.exists(bldb_alias_file)) {
  read.csv2(bldb_alias_file, stringsAsFactors = FALSE, check.names = FALSE) |>
    janitor::clean_names()
} else {
  tibble()
}

if (!"gene_symbol" %in% names(ariba)) ariba[["gene_symbol"]] <- NA_character_
if (!"original_name" %in% names(ariba)) ariba[["original_name"]] <- NA_character_
if (!"ref_name" %in% names(ariba)) ariba[["ref_name"]] <- NA_character_
if (!"gene" %in% names(ariba)) ariba[["gene"]] <- NA_character_
if (!"hit_completeness" %in% names(ariba)) ariba[["hit_completeness"]] <- NA_character_

# ----------------------------------------------------------------------------
# PREPARE AMRFINDER
# ----------------------------------------------------------------------------

amrf_core <- amrf |>
  filter(type == "AMR") |>
  mutate(
    sample_id = normalize_sample_id2(sample_id),
    scope = na_if_empty(scope),
    amrf_gene = normalize_gene_symbol(element_symbol),
    amrf_cov = suppressWarnings(as.numeric(str_replace(percent_coverage_of_reference, ",", "."))),
    amrf_ident = suppressWarnings(as.numeric(str_replace(percent_identity_to_reference, ",", "."))),
    amrf_method = na_if_empty(method),
    amrf_partial = str_detect(replace_na(amrf_method, ""), "PARTIAL") |
      (!is.na(amrf_cov) & amrf_cov < 90),
    amrf_class = na_if_empty(class),
    amrf_subclass = na_if_empty(subclass)
  ) |>
  filter(scope == "core") |>
  mutate(
    gene_family = case_when(
      is_beta_lactamase(amrf_gene) ~ coalesce(beta_family_from_symbol(amrf_gene), amrf_gene),
      TRUE ~ amrf_gene
    )
  ) |>
  group_by(sample_id, gene_family) |>
  arrange(desc(amrf_cov), desc(amrf_ident)) |>
  slice(1) |>
  ungroup()

# ----------------------------------------------------------------------------
# PREPARE ARIBA
# ----------------------------------------------------------------------------

ariba_core <- ariba |>
  mutate(
    sample_id = normalize_sample_id2(sample_id),
    ariba_gene = case_when(
      !is.na(gene_symbol) & gene_symbol != "" ~ gene_symbol,
      !is.na(original_name) & original_name != "" ~ original_name,
      !is.na(ref_name) & ref_name != "" ~ ref_name,
      !is.na(gene) & gene != "" ~ as.character(gene),
      TRUE ~ NA_character_
    ),
    ariba_gene = normalize_gene_symbol(ariba_gene),
    assembly_fraction = suppressWarnings(as.numeric(str_replace(assembly_fraction, ",", "."))),
    pc_ident = suppressWarnings(as.numeric(str_replace(pc_ident, ",", "."))),
    hit_completeness = replace_na(hit_completeness, "fragment"),
    ariba_complete = hit_completeness %in% c("complete", "near_complete"),
    gene_family = case_when(
      is_beta_lactamase(ariba_gene) ~ coalesce(beta_family_from_symbol(ariba_gene), ariba_gene),
      TRUE ~ ariba_gene
    )
  ) |>
  filter(!is.na(ariba_gene)) |>
  group_by(sample_id, gene_family) |>
  arrange(desc(assembly_fraction), desc(pc_ident)) |>
  slice(1) |>
  ungroup()

# ----------------------------------------------------------------------------
# RESOLUTION LOGIC
# ----------------------------------------------------------------------------

resolve_best_gene <- function(amrf_gene, ariba_gene, amrf_partial = FALSE, ariba_complete = FALSE) {
  amrf_gene <- normalize_gene_symbol(amrf_gene)
  ariba_gene <- normalize_gene_symbol(ariba_gene)
  
  # Treat incomplete ARIBA as unavailable for resolution
  if (!isTRUE(ariba_complete)) {
    ariba_gene <- NA_character_
  }
  
  if (is.na(amrf_gene) && is.na(ariba_gene)) return(NA_character_)
  if (is.na(amrf_gene)) return(ariba_gene)
  if (is.na(ariba_gene)) {
    if (isTRUE(amrf_partial)) return(NA_character_)
    return(amrf_gene)
  }
  
  if (is_beta_lactamase(amrf_gene) || is_beta_lactamase(ariba_gene)) {
    fam_amrf <- beta_family_from_symbol(amrf_gene)
    fam_ariba <- beta_family_from_symbol(ariba_gene)
    
    if ((isTRUE(amrf_partial) || is_generic_family_symbol(amrf_gene)) &&
        is_specific_allele_symbol(ariba_gene)) {
      return(ariba_gene)
    }
    
    if (!is.na(fam_amrf) && !is.na(fam_ariba) &&
        fam_amrf == fam_ariba &&
        is_specific_allele_symbol(ariba_gene) &&
        !is_specific_allele_symbol(amrf_gene)) {
      return(ariba_gene)
    }
    
    if (is_specific_allele_symbol(ariba_gene) && !is_specific_allele_symbol(amrf_gene)) {
      return(ariba_gene)
    }
    
    if (is_specific_allele_symbol(amrf_gene)) {
      return(amrf_gene)
    }
    
    if (!is.na(fam_amrf) && !is.na(fam_ariba) && fam_amrf == fam_ariba) {
      return(ariba_gene)
    }
    
    return(amrf_gene)
  }
  
  # Non-beta genes: AMRFinder primary; ARIBA only resolves partial AMRFinder
  if (isTRUE(amrf_partial) && !is.na(ariba_gene)) {
    return(ariba_gene)
  }
  
  amrf_gene
}

resolved_tool <- function(amrf_gene, ariba_gene, resolved_gene) {
  amrf_gene <- normalize_gene_symbol(amrf_gene)
  ariba_gene <- normalize_gene_symbol(ariba_gene)
  resolved_gene <- normalize_gene_symbol(resolved_gene)
  
  case_when(
    !is.na(resolved_gene) & !is.na(amrf_gene) & resolved_gene == amrf_gene ~ "AMRFinder",
    !is.na(resolved_gene) & !is.na(ariba_gene) & resolved_gene == ariba_gene ~ "ariba",
    TRUE ~ NA_character_
  )
}

comparison <- full_join(
  amrf_core |> mutate(sample_id = normalize_sample_id(sample_id)),
  ariba_core |> mutate(sample_id = normalize_sample_id(sample_id)),
  by = c("sample_id", "gene_family"),
  suffix = c(".amrf", ".ariba")
) |>
  mutate(
    comparison_status = case_when(
      !is.na(amrf_gene) & !is.na(ariba_gene) & amrf_gene == ariba_gene ~ "exact_match",
      !is.na(amrf_gene) & !is.na(ariba_gene) ~ "family_match",
      !is.na(amrf_gene) & is.na(ariba_gene) ~ "amrfinder_only",
      is.na(amrf_gene) & !is.na(ariba_gene) ~ "ariba_only",
      TRUE ~ "discordant"
    ),
    resolved_gene = purrr::pmap_chr(
      list(amrf_gene, ariba_gene, amrf_partial, ariba_complete),
      resolve_best_gene
    ),
    tool = purrr::pmap_chr(
      list(amrf_gene, ariba_gene, resolved_gene),
      resolved_tool
    ),
    resolved_gene_norm = normalize_gene_for_lookup(resolved_gene)
  )

# ----------------------------------------------------------------------------
# BLDB PREP
# ----------------------------------------------------------------------------

bldb_main <- bldb |>
  mutate(
    gene_norm = normalize_gene_for_lookup(gene_symbol_guess),
    Gene_symbol = normalize_gene_symbol(gene_symbol_guess),
    Enzyme_family_bldb = na_if_empty(subfamily),
    Molecular_class_bldb = na_if_empty(molecular_class),
    Functional_group_bldb = na_if_empty(bush_jacoby_class),
    reference_bldb = na_if_empty(source_url)
  ) |>
  transmute(
    gene_norm,
    Gene_symbol,
    Enzyme_family_bldb,
    Molecular_class_bldb,
    Functional_group_bldb,
    reference_bldb
  ) |>
  group_by(gene_norm) |>
  summarise(
    Gene_symbol_bldb = first_non_missing(Gene_symbol),
    Enzyme_family_bldb = first_non_missing(Enzyme_family_bldb),
    Molecular_class_bldb = first_non_missing(Molecular_class_bldb),
    Functional_group_bldb = first_non_missing(Functional_group_bldb),
    reference_bldb = first_non_missing(reference_bldb),
    .groups = "drop"
  )

bldb_alias <- if (nrow(bldb_aliases) > 0) {
  bldb_aliases |>
    mutate(
      alias_norm = normalize_gene_for_lookup(alias_bla),
      Enzyme_family_alias = na_if_empty(subfamily),
      Molecular_class_alias = na_if_empty(molecular_class),
      Functional_group_alias = na_if_empty(bush_jacoby_class)
    ) |>
    transmute(
      alias_norm,
      Enzyme_family_alias,
      Molecular_class_alias,
      Functional_group_alias
    ) |>
    group_by(alias_norm) |>
    summarise(
      Enzyme_family_alias = first_non_missing(Enzyme_family_alias),
      Molecular_class_alias = first_non_missing(Molecular_class_alias),
      Functional_group_alias = first_non_missing(Functional_group_alias),
      .groups = "drop"
    )
} else {
  tibble()
}

# ----------------------------------------------------------------------------
# ANNOTATE RESOLVED LONG TABLE
# ----------------------------------------------------------------------------

resolved_long <- comparison |>
  filter(!is.na(resolved_gene)) |>
  mutate(
    resolved_is_beta = is_beta_lactamase(resolved_gene)
  ) |>
  left_join(
    bldb_main,
    by = c("resolved_gene_norm" = "gene_norm")
  )

if (nrow(bldb_alias) > 0) {
  resolved_long <- resolved_long |>
    left_join(
      bldb_alias,
      by = c("resolved_gene_norm" = "alias_norm")
    )
} else {
  resolved_long <- resolved_long |>
    mutate(
      Enzyme_family_alias = NA_character_,
      Molecular_class_alias = NA_character_,
      Functional_group_alias = NA_character_
    )
}

beta_fallback <- purrr::map_dfr(resolved_long$resolved_gene, infer_beta_annotation_from_family)

resolved_long <- resolved_long |>
  bind_cols(
    beta_fallback |>
      rename(
        Enzyme_family_beta_fallback = Enzyme_family,
        Molecular_class_beta_fallback = Molecular_class,
        Functional_group_beta_fallback = Functional_group,
        reference_beta_fallback = reference
      )
  ) |>
  mutate(
    Gene_symbol = resolved_gene,
    
    # requested AMRFinder class/subclass columns
    class = amrf_class,
    subclass = amrf_subclass,
    
    # beta-lactamases only get BLDB/fallback annotation
    Enzyme_family = case_when(
      resolved_is_beta ~ coalesce(Enzyme_family_bldb, Enzyme_family_alias, Enzyme_family_beta_fallback),
      TRUE ~ NA_character_
    ),
    Molecular_class = case_when(
      resolved_is_beta ~ coalesce(Molecular_class_bldb, Molecular_class_alias, Molecular_class_beta_fallback),
      TRUE ~ NA_character_
    ),
    Functional_group = case_when(
      resolved_is_beta ~ coalesce(Functional_group_bldb, Functional_group_alias, Functional_group_beta_fallback),
      TRUE ~ NA_character_
    ),
    reference = case_when(
      resolved_is_beta ~ coalesce(reference_bldb, reference_beta_fallback),
      TRUE ~ NA_character_
    )
  ) |>
  select(
    sample_id,
    tool,
    scope,
    type,
    class,
    subclass,
    amrf_gene,
    ariba_gene,
    amrf_partial,
    ariba_complete,
    comparison_status,
    gene_family,
    resolved_gene,
    Gene_symbol,
    Enzyme_family,
    Molecular_class,
    Functional_group,
    reference,
    everything()
  )

resolved_long <- resolved_long |>
  mutate(scope = na_if_empty(scope)) |>
  filter(scope == "core")

# ----------------------------------------------------------------------------
# MANUSCRIPT TABLE
# ----------------------------------------------------------------------------

manuscript_table <- resolved_long |>
  transmute(
    Gene_symbol = resolved_gene,
    tool,
    class,
    subclass,
    Enzyme_family,
    Molecular_class,
    Functional_group,
    reference
  ) |>
  distinct() |>
  group_by(Gene_symbol) |>
  summarise(
    tool = first_non_missing(tool),
    class = first_non_missing(class),
    subclass = first_non_missing(subclass),
    Enzyme_family = first_non_missing(Enzyme_family),
    Molecular_class = first_non_missing(Molecular_class),
    Functional_group = first_non_missing(Functional_group),
    reference = first_non_missing(reference),
    .groups = "drop"
  ) |>
  arrange(Gene_symbol)

# ----------------------------------------------------------------------------
# BETA-LACTAMASE SUBTABLES
# ----------------------------------------------------------------------------

resolved_long_beta <- resolved_long |>
  filter(is_beta_lactamase(resolved_gene)) |>
  arrange(sample_id, resolved_gene)

manuscript_table_beta <- manuscript_table |>
  filter(is_beta_lactamase(Gene_symbol)) |>
  arrange(Gene_symbol)

# ----------------------------------------------------------------------------
# CHECKS
# ----------------------------------------------------------------------------

beta_unmapped <- manuscript_table_beta |>
  filter(
    is.na(Enzyme_family) |
      is.na(Molecular_class) |
      is.na(Functional_group)
  )

if (nrow(beta_unmapped) > 0) {
  print(beta_unmapped)
  stop("There are still unmapped beta-lactamase entries in the manuscript beta table.")
}

# ----------------------------------------------------------------------------
# WRITE OUTPUT
# ----------------------------------------------------------------------------

write.csv2(
  comparison,
  file.path(outdir, "resolved_amrfinder_ariba_comparison_all_genes.csv"),
  row.names = FALSE
)

write.csv2(
  resolved_long,
  file.path(outdir, "resolved_amr_gene_calls_long_core.csv"),
  row.names = FALSE
)

write.csv2(
  manuscript_table,
  file.path(outdir, "resolved_amr_genes_manuscript_table.csv"),
  row.names = FALSE
)

write.csv2(
  resolved_long_beta,
  file.path(outdir, "resolved_beta_lactamase_calls_long_core.csv"),
  row.names = FALSE
)

write.csv2(
  manuscript_table_beta,
  file.path(outdir, "resolved_beta_lactamases_manuscript_table.csv"),
  row.names = FALSE
)

message("Done.")
message("Wrote:")
message(" - amrfinder_ariba_comparison_all_genes.csv")
message(" - resolved_amr_gene_calls_long_core.csv")
message(" - resolved_amr_genes_manuscript_table.csv")
message(" - resolved_beta_lactamase_calls_long_core.csv")
message(" - resolved_beta_lactamases_manuscript_table.csv")

