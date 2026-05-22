library(dplyr)
library(tidyr)
library(stringr)

source(file='model/modelCommon.R')
outdir <- paste(processedRootRassembly, "genotype", sep="/")


extractBetalactamaseTable <- function(
    resolved_amr_gene_calls_annotated = "resolved_amr_gene_calls_annotated.csv",
    bldb_lookup = NULL
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(tibble)
  })
  
  # =========================================================
  # 1) Read input
  # =========================================================
  x <- if (is.character(resolved_amr_gene_calls_annotated)) {
    read.csv2(resolved_amr_gene_calls_annotated, stringsAsFactors = FALSE)
  } else {
    as.data.frame(resolved_amr_gene_calls_annotated, stringsAsFactors = FALSE)
  }
  
  bldb <- NULL
  if (!is.null(bldb_lookup)) {
    bldb <- if (is.character(bldb_lookup)) {
      read.csv2(bldb_lookup, stringsAsFactors = FALSE)
    } else {
      as.data.frame(bldb_lookup, stringsAsFactors = FALSE)
    }
  }
  
  required_cols <- c(
    "gene_family",
    "resolved_gene",
    "amrf_gene",
    "amrf_partial",
    "ariba_gene",
    "proteinname",
    "molecular_class",
    "bush_jacoby_class"
  )
  
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0) {
    stop(
      "extractBetalactamaseTable(): missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # Keep only CORE genes when scope exists
  if ("scope" %in% names(x)) {
    x$scope <- as.character(x$scope)
    x$scope[x$scope %in% c("", "NA")] <- NA_character_
    x <- x %>% filter(scope == "core")
  }
  
  # =========================================================
  # 2) Helper functions
  # =========================================================
  na_if_empty <- function(z) {
    z <- as.character(z)
    z[z %in% c("", "NA")] <- NA_character_
    z
  }
  
  first_non_missing <- function(z) {
    z <- na_if_empty(z)
    z <- z[!is.na(z)]
    if (length(z) == 0) NA_character_ else z[[1]]
  }
  
  is_beta_lactamase <- function(gene) {
    str_detect(gene, regex("^(bla|ampC)", ignore_case = TRUE))
  }
  
  is_ampc_promoter <- function(gene) {
    str_detect(gene, regex("^ampC_[ACGT]-[0-9]+[ACGT]$", ignore_case = TRUE)) |
      str_detect(gene, regex("^ampC_[0-9A-Za-z._-]+$", ignore_case = TRUE))
  }
  
  derive_enzyme_family <- function(gene_symbol) {
    case_when(
      is.na(gene_symbol) ~ NA_character_,
      str_detect(gene_symbol, regex("^ampC", ignore_case = TRUE)) ~ "AMPC",
      str_detect(gene_symbol, regex("^blaCTX-M", ignore_case = TRUE)) ~ "CTX-M",
      str_detect(gene_symbol, regex("^blaTEM", ignore_case = TRUE)) ~ "TEM",
      str_detect(gene_symbol, regex("^blaSHV", ignore_case = TRUE)) ~ "SHV",
      str_detect(gene_symbol, regex("^blaOXA", ignore_case = TRUE)) ~ "OXA",
      str_detect(gene_symbol, regex("^blaNDM", ignore_case = TRUE)) ~ "NDM",
      str_detect(gene_symbol, regex("^blaDHA", ignore_case = TRUE)) ~ "AMPC",
      str_detect(gene_symbol, regex("^blaLAP", ignore_case = TRUE)) ~ "LAP",
      str_detect(gene_symbol, regex("^blaEC", ignore_case = TRUE)) ~ "AMPC",
      TRUE ~ NA_character_
    )
  }
  
  make_fallback_symbol <- function(gene_symbol) {
    case_when(
      is.na(gene_symbol) ~ NA_character_,
      is_ampc_promoter(gene_symbol) ~ "ampC_promoter",
      str_detect(gene_symbol, regex("^[A-Za-z0-9()-]+-[0-9]+[A-Za-z]$", ignore_case = TRUE)) ~
        str_replace(gene_symbol, "([0-9]+)[A-Za-z]$", "\\1"),
      TRUE ~ gene_symbol
    )
  }
  
  family_prefix <- function(gene_symbol) {
    case_when(
      str_detect(gene_symbol, regex("^blaTEM", ignore_case = TRUE)) ~ "blaTEM",
      str_detect(gene_symbol, regex("^blaCTX-M", ignore_case = TRUE)) ~ "blaCTX-M",
      str_detect(gene_symbol, regex("^blaSHV", ignore_case = TRUE)) ~ "blaSHV",
      str_detect(gene_symbol, regex("^blaOXA", ignore_case = TRUE)) ~ "blaOXA",
      str_detect(gene_symbol, regex("^blaNDM", ignore_case = TRUE)) ~ "blaNDM",
      str_detect(gene_symbol, regex("^blaDHA", ignore_case = TRUE)) ~ "blaDHA",
      str_detect(gene_symbol, regex("^blaLAP", ignore_case = TRUE)) ~ "blaLAP",
      str_detect(gene_symbol, regex("^blaEC", ignore_case = TRUE)) ~ "blaEC",
      str_detect(gene_symbol, regex("^ampC", ignore_case = TRUE)) ~ "ampC",
      TRUE ~ gene_symbol
    )
  }
  
  is_family_only_symbol <- function(gene_symbol) {
    gene_symbol %in% c(
      "blaTEM", "blaCTX-M", "blaSHV", "blaOXA", "blaNDM", "blaDHA", "blaLAP", "blaEC", "ampC"
    )
  }
  
  infer_functional_group_from_subfamily <- function(subfamily) {
    case_when(
      is.na(subfamily) ~ NA_character_,
      str_detect(subfamily, regex("CTX-M", ignore_case = TRUE)) ~ "2be",
      str_detect(subfamily, regex("^OXA", ignore_case = TRUE)) ~ "2d",
      str_detect(subfamily, regex("^LAP", ignore_case = TRUE)) ~ "2b",
      str_detect(subfamily, regex("^NDM", ignore_case = TRUE)) ~ "3a",
      TRUE ~ NA_character_
    )
  }
  
  infer_functional_group_from_gene <- function(gene_symbol, enzyme_family, molecular_class) {
    case_when(
      is.na(gene_symbol) & is.na(enzyme_family) & is.na(molecular_class) ~ NA_character_,
      str_detect(gene_symbol, regex("^blaCTX-M", ignore_case = TRUE)) ~ "2be",
      str_detect(gene_symbol, regex("^blaSHV", ignore_case = TRUE)) ~ "2be",
      str_detect(gene_symbol, regex("^blaOXA", ignore_case = TRUE)) ~ "2d",
      str_detect(gene_symbol, regex("^blaNDM", ignore_case = TRUE)) ~ "3a",
      str_detect(gene_symbol, regex("^blaLAP", ignore_case = TRUE)) ~ "2b",
      str_detect(gene_symbol, regex("^blaDHA", ignore_case = TRUE)) ~ "1",
      str_detect(gene_symbol, regex("^blaEC", ignore_case = TRUE)) ~ "1",
      str_detect(gene_symbol, regex("^ampC", ignore_case = TRUE)) ~ "1",
      enzyme_family == "CTX-M" ~ "2be",
      enzyme_family == "SHV" ~ "2be",
      enzyme_family == "OXA" ~ "2d",
      enzyme_family == "NDM" ~ "3a",
      enzyme_family == "LAP" ~ "2b",
      enzyme_family == "AMPC" ~ "1",
      molecular_class == "C" ~ "1",
      TRUE ~ NA_character_
    )
  }
  
  # =========================================================
  # 3) Resolve calls: ARIBA overrides only when AMRFinder is partial
  # =========================================================
  x2 <- x %>%
    mutate(
      across(
        any_of(c(
          "resolved_gene", "amrf_gene", "ariba_gene", "gene_family",
          "proteinname", "molecular_class", "bush_jacoby_class",
          "source_url"
        )),
        na_if_empty
      ),
      amrf_partial = as.logical(amrf_partial),
      has_amrf = !is.na(amrf_gene),
      has_ariba = !is.na(ariba_gene),
      use_ariba_override = has_amrf & amrf_partial & has_ariba,
      resolved_gene = case_when(
        use_ariba_override ~ ariba_gene,
        has_amrf ~ amrf_gene,
        has_ariba ~ ariba_gene,
        TRUE ~ resolved_gene
      ),
      resolved_source = case_when(
        use_ariba_override ~ "ARIBA_override_partial_AMRFinder",
        has_amrf ~ "AMRFinder",
        has_ariba ~ "ARIBA_only",
        TRUE ~ NA_character_
      )
    ) %>%
    select(-has_amrf, -has_ariba, -use_ariba_override)
  
  # =========================================================
  # 4) Keep only beta-lactamases
  # =========================================================
  beta_long <- x2 %>%
    filter(!is.na(resolved_gene), resolved_gene != "") %>%
    filter(is_beta_lactamase(resolved_gene))
  
  # =========================================================
  # 5) Direct annotation from annotated file
  # =========================================================
  direct_annot <- beta_long %>%
    transmute(
      Gene_symbol = resolved_gene,
      Enzyme_family = gene_family,
      Molecular_class = molecular_class,
      Functional_group = bush_jacoby_class,
      reference = if ("source_url" %in% names(beta_long)) source_url else NA_character_
    )
  
  # =========================================================
  # 6) Optional BLDB lookup by gene symbol
  # =========================================================
  if (!is.null(bldb)) {
    bldb2 <- bldb %>%
      transmute(
        Gene_symbol = na_if_empty(if ("gene_symbol_guess" %in% names(bldb)) gene_symbol_guess else NA),
        Molecular_class_bldb = na_if_empty(if ("molecular_class" %in% names(bldb)) molecular_class else NA),
        Functional_group_bldb = na_if_empty(if ("bush_jacoby_class" %in% names(bldb)) bush_jacoby_class else NA),
        subfamily_bldb = na_if_empty(if ("subfamily" %in% names(bldb)) subfamily else NA),
        reference_bldb = na_if_empty(if ("source_url" %in% names(bldb)) source_url else NA)
      ) %>%
      filter(!is.na(Gene_symbol), Gene_symbol != "") %>%
      group_by(Gene_symbol) %>%
      summarise(
        Molecular_class_bldb = first_non_missing(Molecular_class_bldb),
        Functional_group_bldb = first_non_missing(Functional_group_bldb),
        subfamily_bldb = first_non_missing(subfamily_bldb),
        reference_bldb = first_non_missing(reference_bldb),
        .groups = "drop"
      )
    
    direct_annot <- direct_annot %>%
      left_join(bldb2, by = "Gene_symbol") %>%
      mutate(
        Molecular_class = coalesce(Molecular_class, Molecular_class_bldb),
        Functional_group = coalesce(Functional_group, Functional_group_bldb),
        reference = coalesce(reference, reference_bldb)
      )
  } else {
    direct_annot <- direct_annot %>%
      mutate(
        Molecular_class_bldb = NA_character_,
        Functional_group_bldb = NA_character_,
        subfamily_bldb = NA_character_,
        reference_bldb = NA_character_
      )
  }
  
  # =========================================================
  # 7) Fallback annotation map
  # =========================================================
  fallback_map <- direct_annot %>%
    mutate(Fallback_symbol = make_fallback_symbol(Gene_symbol)) %>%
    group_by(Fallback_symbol) %>%
    summarise(
      Molecular_class_fallback = first_non_missing(Molecular_class),
      Functional_group_fallback = first_non_missing(Functional_group),
      reference_fallback = first_non_missing(reference),
      .groups = "drop"
    )
  
  ampc_promoter_fallback <- tibble(
    Fallback_symbol = "ampC_promoter",
    Molecular_class_fallback = "C",
    Functional_group_fallback = "1",
    reference_fallback = NA_character_
  )
  
  fallback_map <- bind_rows(fallback_map, ampc_promoter_fallback) %>%
    group_by(Fallback_symbol) %>%
    summarise(
      Molecular_class_fallback = first_non_missing(Molecular_class_fallback),
      Functional_group_fallback = first_non_missing(Functional_group_fallback),
      reference_fallback = first_non_missing(reference_fallback),
      .groups = "drop"
    )
  
  # =========================================================
  # 8) Final output with inference
  # =========================================================
  out <- direct_annot %>%
    mutate(
      Fallback_symbol = make_fallback_symbol(Gene_symbol),
      Enzyme_family = derive_enzyme_family(Gene_symbol)
    ) %>%
    left_join(fallback_map, by = "Fallback_symbol") %>%
    mutate(
      Molecular_class = coalesce(Molecular_class, Molecular_class_fallback),
      Functional_group = coalesce(Functional_group, Functional_group_fallback),
      reference = coalesce(reference, reference_fallback),
      
      inferred_from_subfamily = infer_functional_group_from_subfamily(subfamily_bldb),
      inferred_from_family = infer_functional_group_from_gene(
        gene_symbol = Gene_symbol,
        enzyme_family = Enzyme_family,
        molecular_class = Molecular_class
      ),
      
      inference_label = case_when(
        is.na(Functional_group) & !is.na(inferred_from_subfamily) & !is.na(subfamily_bldb) ~
          paste0("infered by subfamily ", subfamily_bldb),
        is.na(Functional_group) & is.na(inferred_from_subfamily) & !is.na(inferred_from_family) & !is.na(Enzyme_family) ~
          paste0("infered by family ", Enzyme_family),
        is.na(Functional_group) & is.na(inferred_from_subfamily) & !is.na(inferred_from_family) & !is.na(Molecular_class) ~
          paste0("infered by molecular class ", Molecular_class),
        TRUE ~ NA_character_
      ),
      
      Functional_group = case_when(
        !is.na(Functional_group) ~ Functional_group,
        is.na(Functional_group) & !is.na(inferred_from_subfamily) ~ inferred_from_subfamily,
        is.na(Functional_group) & is.na(inferred_from_subfamily) & !is.na(inferred_from_family) ~ inferred_from_family,
        TRUE ~ Functional_group
      ),
      
      reference = case_when(
        !is.na(inference_label) ~ inference_label,
        is.na(reference) & !is.na(Enzyme_family) ~ paste0("infered by family ", Enzyme_family),
        TRUE ~ reference
      ),
      
      Molecular_class = case_when(
        is_ampc_promoter(Gene_symbol) ~ coalesce(Molecular_class, "C"),
        TRUE ~ Molecular_class
      ),
      Functional_group = case_when(
        is_ampc_promoter(Gene_symbol) ~ coalesce(Functional_group, "1"),
        TRUE ~ Functional_group
      ),
      
    ) %>%
    select(
      Gene_symbol,
      Enzyme_family,
      Molecular_class,
      Functional_group,
      reference
    ) %>%
    distinct()
  
  # =========================================================
  # 9) Remove generic family-only rows when allele-level rows exist
  # =========================================================
  out <- out %>%
    mutate(
      family_prefix = family_prefix(Gene_symbol),
      is_family_only = is_family_only_symbol(Gene_symbol)
    ) %>%
    group_by(family_prefix) %>%
    filter(!(is_family_only & any(Gene_symbol != family_prefix))) %>%
    ungroup() %>%
    select(-family_prefix, -is_family_only)
  
  # =========================================================
  # 10) Order
  # =========================================================
  out %>%
    arrange(Gene_symbol)
}


betalactamasesTableBLDB <- extractBetalactamaseTable(
  read.csv2(paste(outdir, "resolved_amr_gene_calls_annotated.csv", sep="/")),
  read.csv2(paste(outdir, "bldb_lookup.csv", sep="/"))
)

print(n = 50, betalactamasesTableBLDB)

write.csv2(
  betalactamasesTableBLDB,
  paste(outdir, "betalactamasesTable.csv", sep="/"),
  row.names = FALSE
)