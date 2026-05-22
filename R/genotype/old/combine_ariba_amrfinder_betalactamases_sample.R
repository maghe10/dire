# ============================================================================
# Resolve and annotate beta-lactamase calls with ARIBA-preferred specificity
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(readr)
})

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

na_if_empty <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace(x, "^sample", "")
  x <- str_replace(x, "\\.tsv$", "")
  x <- str_replace(x, "\\.csv$", "")
  x[x == ""] <- NA_character_
  x
}

normalize_beta_symbol <- function(x) {
  x <- na_if_empty(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\s+", "")
  x
}

first_non_missing <- function(x) {
  x <- na_if_empty(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  x[[1]]
}

is_generic_family_symbol <- function(x) {
  x <- normalize_beta_symbol(x)
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
  x <- normalize_beta_symbol(x)
  case_when(
    is.na(x) ~ FALSE,
    str_detect(x, regex("^blaCTX-M-[0-9A-Za-z._-]+$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^blaTEM-[0-9A-Za-z._-]+$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^blaSHV-[0-9A-Za-z._-]+$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^blaOXA-[0-9A-Za-z._-]+$", ignore_case = TRUE)) ~ TRUE,
    str_detect(x, regex("^bla(CMY|DHA|FOX|ACC|ACT|MIR|LAT|MOX)-[0-9A-Za-z._-]+$", ignore_case = TRUE)) ~ TRUE,
    TRUE ~ FALSE
  )
}

beta_family_from_symbol <- function(x) {
  x <- normalize_beta_symbol(x)
  
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
    str_detect(x, regex("^blaEC", ignore_case = TRUE)) ~ "EC",
    str_detect(x, regex("^ampC", ignore_case = TRUE)) ~ "AMPC",
    TRUE ~ NA_character_
  )
}


infer_from_family <- function(gene_symbol) {
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
      fam %in% c("CMY", "DHA", "FOX", "ACC", "ACT", "MIR", "LAT", "MOX", "EC", "AMPC") ~ "AMPC",
      TRUE ~ NA_character_
    ),
    Molecular_class = case_when(
      fam %in% c("CTX-M", "TEM", "SHV", "KPC", "LAP", "CMY", "DHA", "FOX", "ACC", "ACT", "MIR", "LAT", "MOX") ~ "A",
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
      fam %in% c("CMY", "DHA", "FOX", "ACC", "ACT", "MIR", "LAT", "MOX", "EC", "AMPC") ~ "1",
      TRUE ~ NA_character_
    ),
    reference = case_when(
      !is.na(fam) ~ "inferred from gene family",
      TRUE ~ NA_character_
    )
  )
}


# ----------------------------------------------------------------------------
# Resolution logic
# ----------------------------------------------------------------------------

resolve_best_beta_symbol <- function(amrf_gene, ariba_gene, amrf_partial = FALSE) {
  amrf_gene <- normalize_beta_symbol(amrf_gene)
  ariba_gene <- normalize_beta_symbol(ariba_gene)
  
  if (is.na(amrf_gene) && is.na(ariba_gene)) return(NA_character_)
  if (is.na(amrf_gene)) return(ariba_gene)
  if (is.na(ariba_gene)) return(amrf_gene)
  
  fam_amrf <- beta_family_from_symbol(amrf_gene)
  fam_ariba <- beta_family_from_symbol(ariba_gene)
  
  # Prefer ARIBA if AMRFinder is generic or partial and ARIBA is allele-level
  if ((isTRUE(amrf_partial) || is_generic_family_symbol(amrf_gene)) &&
      is_specific_allele_symbol(ariba_gene)) {
    return(ariba_gene)
  }
  
  # Prefer ARIBA if same family and ARIBA is more specific
  if (!is.na(fam_amrf) && !is.na(fam_ariba) &&
      fam_amrf == fam_ariba &&
      is_specific_allele_symbol(ariba_gene) &&
      !is_specific_allele_symbol(amrf_gene)) {
    return(ariba_gene)
  }
  
  # Prefer a specific allele over a generic family symbol
  if (is_specific_allele_symbol(ariba_gene) && !is_specific_allele_symbol(amrf_gene)) {
    return(ariba_gene)
  }
  
  if (is_specific_allele_symbol(amrf_gene)) {
    return(amrf_gene)
  }
  
  # If both generic and same family, prefer ARIBA
  if (!is.na(fam_amrf) && !is.na(fam_ariba) && fam_amrf == fam_ariba) {
    return(ariba_gene)
  }
  
  # Default fallback
  amrf_gene
}

# ----------------------------------------------------------------------------
# Build resolved calls from comparison table
# ----------------------------------------------------------------------------
# Expected input columns:
#   sample_id
#   amrf_gene
#   ariba_gene
#   amrf_partial
#
# If your file uses different names, adjust them below.

build_resolved_beta_calls <- function(comparison_df) {
  required <- c("sample_id", "amrf_gene", "ariba_gene", "amrf_partial")
  missing <- setdiff(required, names(comparison_df))
  if (length(missing) > 0) {
    stop("Missing required columns in comparison_df: ", paste(missing, collapse = ", "))
  }
  
  comparison_df %>%
    mutate(
      sample_id = normalize_sample_id(sample_id),
      amrf_gene = normalize_beta_symbol(amrf_gene),
      ariba_gene = normalize_beta_symbol(ariba_gene),
      amrf_partial = as.logical(amrf_partial),
      resolved_gene = purrr::pmap_chr(
        list(amrf_gene, ariba_gene, amrf_partial),
        resolve_best_beta_symbol
      )
    )
}

# ----------------------------------------------------------------------------
# BLDB annotation
# ----------------------------------------------------------------------------

prepare_bldb_main <- function(bldb_lookup) {
  bldb_lookup %>%
    mutate(
      Gene_symbol = normalize_beta_symbol(gene_symbol_guess),
      Enzyme_family = na_if_empty(subfamily),
      Molecular_class = na_if_empty(molecular_class),
      Functional_group = na_if_empty(bush_jacoby_class),
      reference = na_if_empty(source_url)
    ) %>%
    transmute(
      Gene_symbol,
      Enzyme_family,
      Molecular_class,
      Functional_group,
      reference
    ) %>%
    filter(!is.na(Gene_symbol)) %>%
    group_by(Gene_symbol) %>%
    summarise(
      Enzyme_family = first_non_missing(Enzyme_family),
      Molecular_class = first_non_missing(Molecular_class),
      Functional_group = first_non_missing(Functional_group),
      reference = first_non_missing(reference),
      .groups = "drop"
    )
}

prepare_bldb_alias <- function(bldb_aliases) {
  bldb_aliases %>%
    mutate(
      alias_symbol = normalize_beta_symbol(alias_bla),
      Enzyme_family_alias = na_if_empty(subfamily),
      Molecular_class_alias = na_if_empty(molecular_class),
      Functional_group_alias = na_if_empty(bush_jacoby_class)
    ) %>%
    transmute(
      alias_symbol,
      Enzyme_family_alias,
      Molecular_class_alias,
      Functional_group_alias
    ) %>%
    filter(!is.na(alias_symbol)) %>%
    group_by(alias_symbol) %>%
    summarise(
      Enzyme_family_alias = first_non_missing(Enzyme_family_alias),
      Molecular_class_alias = first_non_missing(Molecular_class_alias),
      Functional_group_alias = first_non_missing(Functional_group_alias),
      .groups = "drop"
    )
}

annotate_beta_calls <- function(beta_calls, bldb_lookup, bldb_aliases = NULL) {
  bldb_main <- prepare_bldb_main(bldb_lookup)
  
  out <- beta_calls %>%
    mutate(
      sample_id = normalize_sample_id(sample_id),
      Gene_symbol = normalize_beta_symbol(Gene_symbol)
    ) %>%
    left_join(bldb_main, by = "Gene_symbol")
  
  if (!is.null(bldb_aliases) && nrow(bldb_aliases) > 0) {
    bldb_alias <- prepare_bldb_alias(bldb_aliases)
    
    out <- out %>%
      left_join(bldb_alias, by = c("Gene_symbol" = "alias_symbol")) %>%
      mutate(
        Enzyme_family = coalesce(Enzyme_family, Enzyme_family_alias),
        Molecular_class = coalesce(Molecular_class, Molecular_class_alias),
        Functional_group = coalesce(Functional_group, Functional_group_alias)
      ) %>%
      select(-Enzyme_family_alias, -Molecular_class_alias, -Functional_group_alias)
  }
  
  # Family fallback for any still-unmapped entries
  fallback_tbl <- purrr::map_dfr(out$Gene_symbol, infer_from_family)
  
  out <- bind_cols(out, fallback_tbl %>%
                     rename(
                       Enzyme_family_fallback = Enzyme_family,
                       Molecular_class_fallback = Molecular_class,
                       Functional_group_fallback = Functional_group,
                       reference_fallback = reference
                     )) %>%
    mutate(
      Enzyme_family = coalesce(Enzyme_family, Enzyme_family_fallback),
      Molecular_class = coalesce(Molecular_class, Molecular_class_fallback),
      Functional_group = coalesce(Functional_group, Functional_group_fallback),
      reference = coalesce(reference, reference_fallback)
    ) %>%
    select(
      -Enzyme_family_fallback,
      -Molecular_class_fallback,
      -Functional_group_fallback,
      -reference_fallback
    )
  
  out
}

# ----------------------------------------------------------------------------
# Wrapper
# ----------------------------------------------------------------------------
# Example expected files:
#   amrfinder_ariba_comparison.csv
#   bldb_lookup.csv
#   bldb_aliases.csv
#
# Output:
#   resolved_and_annotated_beta_calls.csv
resolve_and_annotate_beta_calls <- function(
    comparison_file,
    bldb_lookup_file,
    bldb_aliases_file = NULL,
    out_file = "resolved_and_annotated_beta_calls.csv"
) {
  comparison_df <- read.csv2(comparison_file, stringsAsFactors = FALSE, check.names = FALSE)
  bldb_lookup <- read.csv2(bldb_lookup_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  bldb_aliases <- NULL
  if (!is.null(bldb_aliases_file) && file.exists(bldb_aliases_file)) {
    bldb_aliases <- read.csv2(bldb_aliases_file, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  resolved_calls <- build_resolved_beta_calls(comparison_df)
  
  beta_calls <- resolved_calls %>%
    transmute(
      sample_id = normalize_sample_id(sample_id),
      Gene_symbol = normalize_beta_symbol(resolved_gene)
    ) %>%
    filter(!is.na(Gene_symbol)) %>%
    filter(
      str_detect(Gene_symbol, regex("^(bla|ampC)", ignore_case = TRUE))
    ) %>%
    distinct()
  
  annotated <- annotate_beta_calls(
    beta_calls = beta_calls,
    bldb_lookup = bldb_lookup,
    bldb_aliases = bldb_aliases
  )
  
  beta_unmapped <- annotated %>%
    filter(
      is.na(Enzyme_family) |
        is.na(Molecular_class) |
        is.na(Functional_group)
    )
  
  if (nrow(beta_unmapped) > 0) {
    print(beta_unmapped)
    stop("There are still unmapped beta-lactamase entries.")
  }
  
  write.csv2(annotated, out_file, row.names = FALSE)
  
  invisible(list(
    resolved_calls = resolved_calls,
    beta_calls = beta_calls,
    annotated = annotated
  ))
}


# ----------------------------------------------------------------------------
# Example usage
# ----------------------------------------------------------------------------

EXAMPLE <- function()
{
  source("model/modelcommon.R")
  indir <- file.path(processedRootRassembly, "genotype")
  # file.path(indir,
  res <- resolve_and_annotate_beta_calls(
    comparison_file = file.path(indir,"amrfinder_ariba_comparison.csv"),
    bldb_lookup_file = file.path(indir,"bldb_lookup.csv"),
    bldb_aliases_file = file.path(indir,"bldb_aliases.csv"),
    out_file = file.path(indir,"resolved_and_annotated_beta_calls.csv")
  )
}
