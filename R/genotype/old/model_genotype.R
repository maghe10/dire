# ============================================================================
# prediction_vs_genotype_onefile.R
# ----------------------------------------------------------------------------
# Purpose:
# Analyse whether AI prediction performance is associated with genotype and
# phenotype complexity, both overall and per antibiotic.
#
# Expected inputs:
#   1) A function readCountSampleFrameLong() available in the environment
#      returning a long table like:
#         sample, antibiotic, metric, count, ...
#   2) resolved_amr_gene_calls_annotated.csv
#   3) optionally phenotype_SUR_wide.csv
#
# Main outputs:
#   - sample-level performance tables
#   - sample-antibiotic-level performance tables
#   - genotype feature tables
#   - phenotype feature tables
#   - merged analysis datasets
#   - sample-level models
#   - antibiotic-specific models
#   - plots
#   - heatmaps
#   - phenotype-pattern clustering
#
# Notes:
#   - Designed to be robust to slightly varying column names
#   - Uses only common CRAN packages
#   - Edit paths in the CONFIG section if needed
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(forcats)
  library(broom)
  library(scales)
})

# ----------------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------------

source(file = 'manuscript/manuscriptcommon.R')

indir <- file.path(processedRootRassembly, "genotype")

OUTDIR <- file.path(processedRootRassembly, "prediction_vs_genotype_output")
GENOTYPE_FILE <- file.path(indir,"resolved_amr_gene_calls_annotated.csv")
PHENOTYPE_FILE <-file.path(indir, "phenotype_SUR_wide.csv")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "models"), recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# HELPERS
# ----------------------------------------------------------------------------

message2 <- function(...) {
  cat(paste0(..., "\n"))
}

resolve_col <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    if (required) {
      stop("Could not find any of these columns: ",
           paste(candidates, collapse = ", "))
    }
    return(NULL)
  }
  hit[[1]]
}


safe_ratio <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

make_binary01 <- function(x) {
  case_when(
    is.logical(x) ~ as.integer(x),
    is.numeric(x) ~ as.integer(x > 0),
    is.character(x) ~ as.integer(toupper(x) %in% c("1", "TRUE", "T", "YES", "Y")),
    TRUE ~ NA_integer_
  )
}

save_csv <- function(x, filename) {
  readr::write_csv(x, file.path(OUTDIR, "tables", filename))
}

save_plot <- function(plot_obj, filename, width = 7, height = 5) {
  ggsave(
    filename = file.path(OUTDIR, "plots", filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 600
  )
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(x)
  x <- x[!is.na(x)]
  x <- x[x != ""]
  paste(sort(x), collapse = sep)
}

mode_string <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

# ----------------------------------------------------------------------------
# STEP 1. READ AI PREDICTION COUNTS
# ----------------------------------------------------------------------------

read_prediction_counts <- function() {
  if (!exists("readCountSampleFrameLong", mode = "function")) {
    stop("Function readCountSampleFrameLong() is not available in the environment.")
  }
  
  x <- readCountSampleFrameLong() %>% filter_and_drop(cpmode,"normal") %>% filter_and_drop(mode,MODE_A) %>% filter_and_drop(significanceLevel,NA)
  

  sample_col <- resolve_col(x, c("sample", "sample_id"))
  ab_col <- resolve_col(x, c("antibiotic", "ab", "abx"))
  metric_col <- resolve_col(x, c("metric"))
  count_col <- resolve_col(x, c("count"))
  
  x %>%
    rename(
      sample = !!sample_col,
      antibiotic = !!ab_col,
      metric = !!metric_col,
      count = !!count_col
    ) %>%
    mutate(
      sample = normalize_sample_id(sample),
      antibiotic = as.character(antibiotic),
      metric = as.character(metric),
      count = as.numeric(count)
    )
}

# ----------------------------------------------------------------------------
# STEP 2. DERIVE PERFORMANCE TABLES
# ----------------------------------------------------------------------------

derive_prediction_tables <- function(pred_long) {
  required_metrics <- c(
    "correctS", "correctR", "falseS", "falseR",
    "zerolabelS", "zerolabelR", "twolabelS", "twolabelR",
    "S", "R", "total"
  )
  
  perf_long <- pred_long %>%
    filter(metric %in% required_metrics) %>%
    group_by(sample, antibiotic, metric) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = metric,
      values_from = count,
      values_fill = 0
    ) %>%
    mutate(
      correct = correctS + correctR,
      errors = falseS + falseR,
      ambiguous = zerolabelS + zerolabelR + twolabelS + twolabelR,
      called = correct + errors,
      n_true_S = S,
      n_true_R = R,
      acc_called = safe_ratio(correct, called),
      acc_total = safe_ratio(correct, total),
      error_total = safe_ratio(errors, total),
      ambiguity_total = safe_ratio(ambiguous, total),
      ME_called = safe_ratio(falseS, called),
      VME_called = safe_ratio(falseR, called),
      ME_total = safe_ratio(falseS, total),
      VME_total = safe_ratio(falseR, total),
      correctS_rate_total = safe_ratio(correctS, total),
      correctR_rate_total = safe_ratio(correctR, total)
    )
  
  perf_sample <- perf_long %>%
    group_by(sample) %>%
    summarise(
      n_antibiotics = n_distinct(antibiotic),
      total = sum(total, na.rm = TRUE),
      called = sum(called, na.rm = TRUE),
      correct = sum(correct, na.rm = TRUE),
      errors = sum(errors, na.rm = TRUE),
      ambiguous = sum(ambiguous, na.rm = TRUE),
      correctS = sum(correctS, na.rm = TRUE),
      correctR = sum(correctR, na.rm = TRUE),
      falseS = sum(falseS, na.rm = TRUE),
      falseR = sum(falseR, na.rm = TRUE),
      n_true_S = sum(n_true_S, na.rm = TRUE),
      n_true_R = sum(n_true_R, na.rm = TRUE),
      acc_called = safe_ratio(correct, called),
      acc_total = safe_ratio(correct, total),
      error_total = safe_ratio(errors, total),
      ambiguity_total = safe_ratio(ambiguous, total),
      ME_called = safe_ratio(falseS, called),
      VME_called = safe_ratio(falseR, called),
      ME_total = safe_ratio(falseS, total),
      VME_total = safe_ratio(falseR, total),
      S_fraction = safe_ratio(n_true_S, total),
      R_fraction = safe_ratio(n_true_R, total),
      phenotype_balance = abs(S_fraction - 0.5),
      .groups = "drop"
    )
  
  list(
    perf_sample_ab = perf_long,
    perf_sample = perf_sample
  )
}

# ----------------------------------------------------------------------------
# STEP 3. READ AND FEATURE-ENGINEER GENOTYPE
# ----------------------------------------------------------------------------

read_genotype_table <- function(file = GENOTYPE_FILE) {
  if (!file.exists(file)) {
    stop("Genotype file not found: ", file)
  }
  
  x <- read.csv2(file, stringsAsFactors = FALSE, check.names = FALSE)
  
  sample_col <- resolve_col(x, c("sample", "sample_id"))
  gene_col <- resolve_col(x, c("resolved_gene", "Gene_symbol", "gene", "gene_symbol"))
  family_col <- resolve_col(
    x,
    c("gene_family", "Enzyme_family", "enzyme_family"),
    required = FALSE
  )
  class_col <- resolve_col(
    x,
    c("molecular_class", "Molecular_class"),
    required = FALSE
  )
  bush_col <- resolve_col(
    x,
    c("bush_jacoby_class", "Functional_group"),
    required = FALSE
  )
  
  x2 <- x
  
  # build canonical sample column without rename collisions
  x2$sample <- x[[sample_col]]
  x2$sample <- normalize_sample_id(x2$sample)
  
  # build canonical gene column without rename collisions
  x2$gene <- x[[gene_col]]
  x2$gene <- as.character(x2$gene)
  
  # canonical optional columns
  x2$gene_family <- if (!is.null(family_col)) as.character(x[[family_col]]) else NA_character_
  x2$molecular_class <- if (!is.null(class_col)) as.character(x[[class_col]]) else NA_character_
  x2$bush_jacoby_class <- if (!is.null(bush_col)) as.character(x[[bush_col]]) else NA_character_
  
  x2
}
infer_beta_group <- function(gene) {
  case_when(
    is.na(gene) ~ NA_character_,
    str_detect(gene, regex("^blaCTX-M", ignore_case = TRUE)) ~ "ESBL_CTXM",
    str_detect(gene, regex("^bla(CMY|DHA|FOX|ACC|ACT|MIR|LAT|MOX)", ignore_case = TRUE)) ~ "AmpC_acquired",
    str_detect(gene, regex("^blaOXA-1$", ignore_case = TRUE)) ~ "OXA1",
    str_detect(gene, regex("^bla(KPC|NDM|VIM|IMP|OXA-48|OXA-181|OXA-232|GES|SME|IMI|NMC)", ignore_case = TRUE)) ~ "Carbapenemase",
    str_detect(gene, regex("^blaTEM-1([A-Z]|$)|^blaSHV-1$", ignore_case = TRUE)) ~ "Penicillinase",
    str_detect(gene, regex("^blaTEM", ignore_case = TRUE)) ~ "TEM_other",
    str_detect(gene, regex("^blaSHV", ignore_case = TRUE)) ~ "SHV_other",
    str_detect(gene, regex("^blaOXA", ignore_case = TRUE)) ~ "OXA_other",
    str_detect(gene, regex("^blaEC", ignore_case = TRUE)) ~ "Intrinsic_EC_AmpC",
    str_detect(gene, regex("^ampC", ignore_case = TRUE)) ~ "AmpC_related",
    str_detect(gene, regex("^blaLAP", ignore_case = TRUE)) ~ "LAP",
    str_detect(gene, regex("^bla", ignore_case = TRUE)) ~ "Other_beta_lactamase",
    TRUE ~ "Other_AMR"
  )
}

derive_genotype_features <- function(geno_long, all_samples = NULL) {
  geno2 <- geno_long %>%
    mutate(
      gene = as.character(gene),
      gene_family = as.character(gene_family),
      molecular_class = as.character(molecular_class),
      bush_jacoby_class = as.character(bush_jacoby_class),
      beta_group = infer_beta_group(gene),
      
      has_beta_lactamase = str_detect(gene, regex("^(bla|ampC)", ignore_case = TRUE)),
      has_CTXM = str_detect(gene, regex("CTX-M", ignore_case = TRUE)),
      has_AmpC_gene = str_detect(gene, regex("^(bla(CMY|DHA|FOX|ACC|ACT|MIR|LAT|MOX)|ampC)", ignore_case = TRUE)),
      has_TEM = str_detect(gene, regex("TEM", ignore_case = TRUE)),
      has_SHV = str_detect(gene, regex("SHV", ignore_case = TRUE)),
      has_OXA = str_detect(gene, regex("OXA", ignore_case = TRUE)),
      has_OXA1 = str_detect(gene, regex("^blaOXA-1$", ignore_case = TRUE)),
      has_carbapenemase_gene = str_detect(gene, regex("^(bla(KPC|NDM|VIM|IMP|OXA-48|OXA-181|OXA-232|GES|SME|IMI|NMC))", ignore_case = TRUE)),
      has_qnr = str_detect(gene, regex("^qnr", ignore_case = TRUE)),
      has_dfr = str_detect(gene, regex("^dfr", ignore_case = TRUE)),
      has_sul = str_detect(gene, regex("^sul[123]$", ignore_case = TRUE)),
      has_aac6 = str_detect(gene, regex("aac\\(6", ignore_case = TRUE)),
      has_aac3 = str_detect(gene, regex("aac\\(3", ignore_case = TRUE)),
      has_aph2 = str_detect(gene, regex("aph\\(2", ignore_case = TRUE)),
      has_armA = str_detect(gene, regex("^armA$", ignore_case = TRUE)),
      has_rmt = str_detect(gene, regex("^rmt", ignore_case = TRUE)),
      has_oqx = str_detect(gene, regex("^oqx(A|B)$", ignore_case = TRUE)),
      has_gyr_par = str_detect(gene, regex("^(gyrA|gyrB|parC|parE)", ignore_case = TRUE))
    )
  
  geno_sample_detected <- geno2 %>%
    group_by(sample) %>%
    summarise(
      n_genes = n(),
      n_unique_genes = n_distinct(gene),
      n_beta_genes = sum(has_beta_lactamase, na.rm = TRUE),
      n_beta_groups = n_distinct(beta_group[!is.na(beta_group)]),
      
      genes = collapse_unique(gene),
      gene_families = collapse_unique(gene_family),
      beta_groups = collapse_unique(beta_group),
      
      ESBL = any(has_CTXM, na.rm = TRUE),
      AmpC = any(has_AmpC_gene, na.rm = TRUE),
      TEM = any(has_TEM, na.rm = TRUE),
      SHV = any(has_SHV, na.rm = TRUE),
      OXA = any(has_OXA, na.rm = TRUE),
      OXA1 = any(has_OXA1, na.rm = TRUE),
      carbapenemase = any(has_carbapenemase_gene, na.rm = TRUE),
      
      qnr = any(has_qnr, na.rm = TRUE),
      gyr_par = any(has_gyr_par, na.rm = TRUE),
      dfr = any(has_dfr, na.rm = TRUE),
      sul = any(has_sul, na.rm = TRUE),
      oqx = any(has_oqx, na.rm = TRUE),
      aminoglycoside_mod = any(has_aac6 | has_aac3 | has_aph2 | has_armA | has_rmt, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (!is.null(all_samples)) {
    geno_sample <- tibble(sample = unique(as.character(all_samples))) %>%
      left_join(geno_sample_detected, by = "sample") %>%
      mutate(
        n_genes = coalesce(n_genes, 0L),
        n_unique_genes = coalesce(n_unique_genes, 0L),
        n_beta_genes = coalesce(n_beta_genes, 0L),
        n_beta_groups = coalesce(n_beta_groups, 0L),
        genes = if_else(is.na(genes), "", genes),
        gene_families = if_else(is.na(gene_families), "", gene_families),
        beta_groups = if_else(is.na(beta_groups), "", beta_groups),
        ESBL = coalesce(ESBL, FALSE),
        AmpC = coalesce(AmpC, FALSE),
        TEM = coalesce(TEM, FALSE),
        SHV = coalesce(SHV, FALSE),
        OXA = coalesce(OXA, FALSE),
        OXA1 = coalesce(OXA1, FALSE),
        carbapenemase = coalesce(carbapenemase, FALSE),
        qnr = coalesce(qnr, FALSE),
        gyr_par = coalesce(gyr_par, FALSE),
        dfr = coalesce(dfr, FALSE),
        sul = coalesce(sul, FALSE),
        oqx = coalesce(oqx, FALSE),
        aminoglycoside_mod = coalesce(aminoglycoside_mod, FALSE)
      )
  } else {
    geno_sample <- geno_sample_detected
  }
  
  geno_sample <- geno_sample %>%
    mutate(
      beta_complexity = case_when(
        carbapenemase ~ "Carbapenemase",
        ESBL & AmpC ~ "ESBL_plus_AmpC",
        ESBL ~ "ESBL_only",
        AmpC ~ "AmpC_only",
        OXA1 ~ "OXA1_only_or_with_others",
        TEM & !ESBL & !AmpC & !OXA1 ~ "TEM_SHV_like_only",
        n_beta_genes > 0 ~ "Other_beta",
        TRUE ~ "No_beta_gene"
      ),
      resistance_gene_burden = case_when(
        n_unique_genes <= 1 ~ "0_1",
        n_unique_genes <= 3 ~ "2_3",
        n_unique_genes <= 6 ~ "4_6",
        TRUE ~ "7_plus"
      ),
      beta_complexity = factor(
        beta_complexity,
        levels = c(
          "No_beta_gene",
          "TEM_SHV_like_only",
          "OXA1_only_or_with_others",
          "AmpC_only",
          "ESBL_only",
          "ESBL_plus_AmpC",
          "Carbapenemase",
          "Other_beta"
        )
      ),
      resistance_gene_burden = factor(
        resistance_gene_burden,
        levels = c("0_1", "2_3", "4_6", "7_plus")
      )
    )
  
  list(
    genotype_long = geno2,
    genotype_sample = geno_sample
  )
}

# ----------------------------------------------------------------------------
# STEP 4. READ AND FEATURE-ENGINEER PHENOTYPE (OPTIONAL)
# ----------------------------------------------------------------------------

read_phenotype_wide <- function(file = PHENOTYPE_FILE) {
  if (!file.exists(file)) {
    return(NULL)
  }
  
  x <- read.csv2(file, stringsAsFactors = FALSE, check.names = FALSE)
  sample_col <- resolve_col(x, c("sample", "sample_id"))
  x %>%
    rename(sample = !!sample_col) %>%
    mutate(sample = normalize_sample_id(sample))
}

derive_phenotype_features <- function(pheno_wide) {
  if (is.null(pheno_wide)) return(NULL)
  
  pheno_long <- pheno_wide %>%
    pivot_longer(
      cols = -sample,
      names_to = "antibiotic",
      values_to = "SUR"
    ) %>%
    mutate(
      SUR = na_if(as.character(SUR), ""),
      is_R = SUR == "R",
      is_S = SUR == "S",
      is_U = SUR == "U"
    )
  
  pheno_sample <- pheno_long %>%
    group_by(sample) %>%
    summarise(
      n_pheno_tested = sum(!is.na(SUR)),
      n_R = sum(is_R, na.rm = TRUE),
      n_S = sum(is_S, na.rm = TRUE),
      n_U = sum(is_U, na.rm = TRUE),
      R_fraction_pheno = safe_ratio(n_R, n_pheno_tested),
      U_fraction_pheno = safe_ratio(n_U, n_pheno_tested),
      phenotype_complexity = case_when(
        n_R <= 1 ~ "low",
        n_R <= 3 ~ "medium",
        TRUE ~ "high"
      ),
      phenotype_pattern = collapse_unique(
        pheno_long$antibiotic[pheno_long$sample == first(sample) &
                                pheno_long$SUR[pheno_long$sample == first(sample)] == "R"]
      ),
      .groups = "drop"
    ) %>%
    mutate(
      phenotype_complexity = factor(
        phenotype_complexity,
        levels = c("low", "medium", "high")
      )
    )
  
  list(
    phenotype_long = pheno_long,
    phenotype_sample = pheno_sample
  )
}

# ----------------------------------------------------------------------------
# STEP 5. MERGE ANALYSIS DATASETS
# ----------------------------------------------------------------------------

merge_analysis_data <- function(perf, geno, pheno = NULL) {
  sample_df <- perf$perf_sample %>%
    left_join(geno$genotype_sample, by = "sample") %>%
    mutate(
      n_genes = coalesce(n_genes, 0L),
      n_unique_genes = coalesce(n_unique_genes, 0L),
      n_beta_genes = coalesce(n_beta_genes, 0L),
      n_beta_groups = coalesce(n_beta_groups, 0L),
      
      genes = coalesce(genes, ""),
      gene_families = coalesce(gene_families, ""),
      beta_groups = coalesce(beta_groups, ""),
      
      ESBL = coalesce(ESBL, FALSE),
      AmpC = coalesce(AmpC, FALSE),
      TEM = coalesce(TEM, FALSE),
      SHV = coalesce(SHV, FALSE),
      OXA = coalesce(OXA, FALSE),
      OXA1 = coalesce(OXA1, FALSE),
      carbapenemase = coalesce(carbapenemase, FALSE),
      
      qnr = coalesce(qnr, FALSE),
      gyr_par = coalesce(gyr_par, FALSE),
      dfr = coalesce(dfr, FALSE),
      sul = coalesce(sul, FALSE),
      oqx = coalesce(oqx, FALSE),
      aminoglycoside_mod = coalesce(aminoglycoside_mod, FALSE)
    ) %>%
    mutate(
      beta_complexity = case_when(
        carbapenemase ~ "Carbapenemase",
        ESBL & AmpC ~ "ESBL_plus_AmpC",
        ESBL ~ "ESBL_only",
        AmpC ~ "AmpC_only",
        OXA1 ~ "OXA1_only_or_with_others",
        TEM & !ESBL & !AmpC & !OXA1 ~ "TEM_SHV_like_only",
        n_beta_genes > 0 ~ "Other_beta",
        TRUE ~ "No_beta_gene"
      ),
      resistance_gene_burden = case_when(
        n_unique_genes <= 1 ~ "0_1",
        n_unique_genes <= 3 ~ "2_3",
        n_unique_genes <= 6 ~ "4_6",
        TRUE ~ "7_plus"
      ),
      beta_complexity = factor(
        beta_complexity,
        levels = c(
          "No_beta_gene",
          "TEM_SHV_like_only",
          "OXA1_only_or_with_others",
          "AmpC_only",
          "ESBL_only",
          "ESBL_plus_AmpC",
          "Carbapenemase",
          "Other_beta"
        )
      ),
      resistance_gene_burden = factor(
        resistance_gene_burden,
        levels = c("0_1", "2_3", "4_6", "7_plus")
      )
    )
  
  if (!is.null(pheno) && !is.null(pheno$phenotype_sample)) {
    sample_df <- sample_df %>%
      left_join(pheno$phenotype_sample, by = "sample")
  }
  
  sample_ab_df <- perf$perf_sample_ab %>%
    left_join(
      sample_df %>%
        select(
          sample,
          n_genes, n_unique_genes, n_beta_genes, n_beta_groups,
          genes, gene_families, beta_groups,
          ESBL, AmpC, TEM, SHV, OXA, OXA1, carbapenemase,
          qnr, gyr_par, dfr, sul, oqx, aminoglycoside_mod,
          beta_complexity, resistance_gene_burden
        ),
      by = "sample"
    )
  
  if (!is.null(pheno) && !is.null(pheno$phenotype_long)) {
    sample_ab_df <- sample_ab_df %>%
      left_join(
        pheno$phenotype_long %>% select(sample, antibiotic, SUR),
        by = c("sample", "antibiotic")
      )
  }
  
  list(
    sample_df = sample_df,
    sample_ab_df = sample_ab_df
  )
}


# ----------------------------------------------------------------------------
# STEP 6. SAMPLE-LEVEL MODELS
# ----------------------------------------------------------------------------

fit_sample_level_models <- function(sample_df) {
  df <- sample_df %>%
    mutate(
      ESBL = make_binary01(ESBL),
      AmpC = make_binary01(AmpC),
      TEM = make_binary01(TEM),
      OXA1 = make_binary01(OXA1),
      carbapenemase = make_binary01(carbapenemase),
      gene_burden_num = as.numeric(resistance_gene_burden),
      phenotype_complexity_num = if ("phenotype_complexity" %in% names(.)) {
        as.numeric(phenotype_complexity)
      } else {
        NA_real_
      }
    )
  
  model_specs <- list(
    acc_total = acc_total ~ ESBL + AmpC + TEM + OXA1 + carbapenemase + n_unique_genes,
    ambiguity_total = ambiguity_total ~ ESBL + AmpC + TEM + OXA1 + carbapenemase + n_unique_genes,
    ME_total = ME_total ~ ESBL + AmpC + TEM + OXA1 + carbapenemase + n_unique_genes,
    VME_total = VME_total ~ ESBL + AmpC + TEM + OXA1 + carbapenemase + n_unique_genes
  )
  
  fits <- purrr::imap(model_specs, function(form, name) {
    dat <- df %>% filter(complete.cases(model.frame(form, data = df)))
    if (nrow(dat) < 10) return(NULL)
    fit <- lm(form, data = dat)
    list(name = name, fit = fit)
  })
  
  fits <- fits[!vapply(fits, is.null, logical(1))]
  
  tidy_tbl <- purrr::map_dfr(fits, function(x) {
    broom::tidy(x$fit, conf.int = TRUE) %>%
      mutate(model = x$name, .before = 1)
  })
  
  glance_tbl <- purrr::map_dfr(fits, function(x) {
    broom::glance(x$fit) %>%
      mutate(model = x$name, .before = 1)
  })
  
  list(
    fits = fits,
    tidy = tidy_tbl,
    glance = glance_tbl
  )
}

# ----------------------------------------------------------------------------
# STEP 7. ANTIBIOTIC-SPECIFIC MODELS
# ----------------------------------------------------------------------------

fit_antibiotic_specific_models <- function(sample_ab_df, min_n = 10) {
  df <- sample_ab_df %>%
    mutate(
      ESBL = make_binary01(ESBL),
      AmpC = make_binary01(AmpC),
      TEM = make_binary01(TEM),
      OXA1 = make_binary01(OXA1),
      carbapenemase = make_binary01(carbapenemase)
    ) %>%
    filter(!is.na(acc_total))
  
  by_ab <- split(df, df$antibiotic)
  
  fits <- purrr::imap(by_ab, function(dat, ab) {
    dat2 <- dat %>%
      filter(complete.cases(
        acc_total, ESBL, AmpC, TEM, OXA1, carbapenemase, n_unique_genes
      ))
    
    if (nrow(dat2) < min_n) return(NULL)
    
    fit <- lm(
      acc_total ~ ESBL + AmpC + TEM + OXA1 + carbapenemase + n_unique_genes,
      data = dat2
    )
    
    list(antibiotic = ab, fit = fit, n = nrow(dat2))
  })
  
  fits <- fits[!vapply(fits, is.null, logical(1))]
  
  tidy_tbl <- purrr::map_dfr(fits, function(x) {
    broom::tidy(x$fit, conf.int = TRUE) %>%
      mutate(antibiotic = x$antibiotic, n = x$n, .before = 1)
  })
  
  glance_tbl <- purrr::map_dfr(fits, function(x) {
    broom::glance(x$fit) %>%
      mutate(antibiotic = x$antibiotic, n = x$n, .before = 1)
  })
  
  list(
    fits = fits,
    tidy = tidy_tbl,
    glance = glance_tbl
  )
}

# ----------------------------------------------------------------------------
# STEP 8. PLOTS
# ----------------------------------------------------------------------------

plot_sample_level <- function(sample_df) {
  df <- sample_df
  
  p1 <- ggplot(df, aes(x = beta_complexity, y = acc_total)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Genotype group",
      y = "Overall accuracy (correct / total)",
      title = "AI prediction accuracy by genotype group"
    )
  
  p2 <- ggplot(df, aes(x = n_unique_genes, y = acc_total)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    theme_bw() +
    labs(
      x = "Number of unique AMR genes",
      y = "Overall accuracy (correct / total)",
      title = "AI prediction accuracy vs gene burden"
    )
  
  p3 <- ggplot(df, aes(x = beta_complexity, y = ambiguity_total)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Genotype group",
      y = "Ambiguity rate",
      title = "AI ambiguity by genotype group"
    )
  
  p4 <- ggplot(df, aes(x = beta_complexity, y = VME_total)) +
    geom_violin(trim = FALSE) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Genotype group",
      y = "VME rate",
      title = "AI VME rate by genotype group"
    )
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}

plot_antibiotic_group_summary <- function(sample_ab_df) {
  summary_tbl <- sample_ab_df %>%
    group_by(antibiotic, beta_complexity) %>%
    summarise(
      mean_acc = mean(acc_total, na.rm = TRUE),
      n = sum(!is.na(acc_total)),
      .groups = "drop"
    )
  
  p <- ggplot(summary_tbl, aes(x = beta_complexity, y = mean_acc)) +
    geom_col() +
    facet_wrap(~ antibiotic, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Genotype group",
      y = "Mean overall accuracy",
      title = "Antibiotic-specific accuracy by genotype group"
    )
  
  list(summary_tbl = summary_tbl, plot = p)
}

plot_heatmaps <- function(sample_ab_df) {
  heat1 <- sample_ab_df %>%
    group_by(antibiotic, beta_complexity) %>%
    summarise(
      mean_acc = mean(acc_total, na.rm = TRUE),
      .groups = "drop"
    )
  
  p1 <- ggplot(heat1, aes(x = antibiotic, y = beta_complexity, fill = mean_acc)) +
    geom_tile() +
    scale_fill_viridis_c(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = "Genotype group",
      fill = "Mean accuracy",
      title = "Heatmap: genotype group vs antibiotic accuracy"
    )
  
  heat2 <- sample_ab_df %>%
    group_by(sample, antibiotic) %>%
    summarise(acc_total = mean(acc_total, na.rm = TRUE), .groups = "drop")
  
  p2 <- ggplot(heat2, aes(x = antibiotic, y = fct_reorder(sample, acc_total, .fun = mean, na.rm = TRUE), fill = acc_total)) +
    geom_tile() +
    scale_fill_viridis_c(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = "Sample",
      fill = "Accuracy",
      title = "Heatmap: sample-specific prediction accuracy"
    )
  
  list(
    heat_beta_ab_tbl = heat1,
    heat_sample_ab_tbl = heat2,
    p1 = p1,
    p2 = p2
  )
}

plot_model_coefficients <- function(model_tidy, title = "Model coefficients") {
  df <- model_tidy %>%
    filter(term != "(Intercept)")
  
  ggplot(df, aes(x = estimate, y = fct_rev(term))) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point() +
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      width = 0.2,
      orientation = "y"
    ) +
    facet_wrap(~ model, scales = "free_x") +
    theme_bw() +
    labs(
      x = "Estimate",
      y = NULL,
      title = title
    )
}

plot_ab_coefficients <- function(ab_model_tidy) {
  df <- ab_model_tidy %>%
    filter(term %in% c("ESBL", "AmpC", "TEM", "OXA1", "carbapenemase", "n_unique_genes"))
  
  ggplot(df, aes(x = antibiotic, y = estimate, ymin = conf.low, ymax = conf.high)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_pointrange() +
    facet_wrap(~ term, scales = "free_y") +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = "Coefficient estimate",
      title = "Antibiotic-specific model coefficients"
    )
}

# ----------------------------------------------------------------------------
# STEP 9. PHENOTYPE-PATTERN CLUSTERING
# ----------------------------------------------------------------------------

cluster_isolates_by_prediction_pattern <- function(sample_ab_df) {
  wide <- sample_ab_df %>%
    select(sample, antibiotic, acc_total) %>%
    pivot_wider(
      names_from = antibiotic,
      values_from = acc_total
    )
  
  if (nrow(wide) < 3) return(NULL)
  
  mat_df <- as.data.frame(wide, stringsAsFactors = FALSE)
  rownames(mat_df) <- mat_df$sample
  mat_df$sample <- NULL
  mat <- as.matrix(mat_df)
  
  # Impute missing values with column means
  for (j in seq_len(ncol(mat))) {
    miss <- is.na(mat[, j])
    if (any(miss)) {
      col_mean <- mean(mat[, j], na.rm = TRUE)
      if (is.nan(col_mean)) col_mean <- 0
      mat[miss, j] <- col_mean
    }
  }
  
  d <- dist(mat)
  hc <- hclust(d, method = "ward.D2")
  
  cluster_tbl <- tibble(
    sample = rownames(mat),
    cluster_k3 = factor(cutree(hc, k = min(3, nrow(mat)))),
    cluster_k4 = factor(cutree(hc, k = min(4, nrow(mat))))
  )
  
  list(
    matrix = mat,
    hclust = hc,
    clusters = cluster_tbl
  )
}
plot_cluster_heatmap <- function(sample_ab_df, cluster_tbl) {
  df <- sample_ab_df %>%
    left_join(cluster_tbl, by = "sample") %>%
    group_by(sample, cluster_k3) %>%
    mutate(sample_mean_acc = mean(acc_total, na.rm = TRUE)) %>%
    ungroup()
  
  sample_order <- df %>%
    distinct(sample, cluster_k3, sample_mean_acc) %>%
    arrange(cluster_k3, sample_mean_acc) %>%
    pull(sample)
  
  df <- df %>%
    mutate(sample = factor(sample, levels = sample_order))
  
  ggplot(df, aes(x = antibiotic, y = sample, fill = acc_total)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    facet_grid(cluster_k3 ~ ., scales = "free_y", space = "free_y") +
    theme_bw() +
    labs(
      x = "Antibiotic",
      y = "Sample",
      fill = "Accuracy",
      title = "Clustering of isolates by AI prediction pattern"
    )
}

# ----------------------------------------------------------------------------
# STEP 10. SUMMARIES
# ----------------------------------------------------------------------------

make_summary_tables <- function(sample_df, sample_ab_df) {
  overall <- sample_df %>%
    summarise(
      n_samples = n(),
      mean_acc_total = mean(acc_total, na.rm = TRUE),
      median_acc_total = median(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_ME_total = mean(ME_total, na.rm = TRUE),
      mean_VME_total = mean(VME_total, na.rm = TRUE)
    )
  
  by_group <- sample_df %>%
    group_by(beta_complexity) %>%
    summarise(
      n_samples = n(),
      mean_acc_total = mean(acc_total, na.rm = TRUE),
      median_acc_total = median(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_ME_total = mean(ME_total, na.rm = TRUE),
      mean_VME_total = mean(VME_total, na.rm = TRUE),
      .groups = "drop"
    )
  
  by_antibiotic <- sample_ab_df %>%
    group_by(antibiotic) %>%
    summarise(
      n_samples = sum(!is.na(acc_total)),
      mean_acc_total = mean(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_ME_total = mean(ME_total, na.rm = TRUE),
      mean_VME_total = mean(VME_total, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    overall = overall,
    by_group = by_group,
    by_antibiotic = by_antibiotic
  )
}

# ----------------------------------------------------------------------------
# STEP 11. MAIN RUNNER
# ----------------------------------------------------------------------------

run_prediction_vs_genotype_analysis <- function() {
  message2("Reading prediction counts ...")
  pred_long <- read_prediction_counts()
  
  message2("Deriving prediction performance tables ...")
  perf <- derive_prediction_tables(pred_long)
  
  message2("Reading genotype data ...")
  geno_long <- read_genotype_table(GENOTYPE_FILE)
  
  message2("Deriving genotype features ...")
  all_samples <- unique(perf$perf_sample$sample)
  geno <- derive_genotype_features(geno_long, all_samples = all_samples)
  
  message2("Checking sample overlap ...")
  pred_samples <- tibble(sample = sort(unique(perf$perf_sample$sample)))
  geno_samples <- tibble(sample = sort(unique(geno$genotype_sample$sample)))
  
  only_in_pred <- pred_samples %>% anti_join(geno_samples, by = "sample")
  only_in_geno <- geno_samples %>% anti_join(pred_samples, by = "sample")
  
  message2("Samples only in prediction table: ", nrow(only_in_pred))
  message2("Samples only in genotype table: ", nrow(only_in_geno))
  
  if (nrow(only_in_pred) > 0) {
    print(head(only_in_pred, 20))
  }
  if (nrow(only_in_geno) > 0) {
    print(head(only_in_geno, 20))
  }
  
  
  message2("Reading phenotype data (optional) ...")
  pheno_wide <- read_phenotype_wide(PHENOTYPE_FILE)
  pheno <- derive_phenotype_features(pheno_wide)
  
  message2("Merging analysis datasets ...")
  merged <- merge_analysis_data(perf, geno, pheno)
  
  message2("Fitting sample-level models ...")
  sample_models <- fit_sample_level_models(merged$sample_df)
  
  message2("Fitting antibiotic-specific models ...")
  ab_models <- fit_antibiotic_specific_models(merged$sample_ab_df)
  
  message2("Creating summary tables ...")
  summaries <- make_summary_tables(merged$sample_df, merged$sample_ab_df)
  
  message2("Saving tables ...")
  save_csv(pred_long, "prediction_counts_long.csv")
  save_csv(perf$perf_sample_ab, "prediction_performance_sample_antibiotic.csv")
  save_csv(perf$perf_sample, "prediction_performance_sample.csv")
  
  save_csv(geno$genotype_long, "genotype_long_features.csv")
  save_csv(geno$genotype_sample, "genotype_sample_features.csv")
  
  if (!is.null(pheno) && !is.null(pheno$phenotype_long)) {
    save_csv(pheno$phenotype_long, "phenotype_long.csv")
  }
  if (!is.null(pheno) && !is.null(pheno$phenotype_sample)) {
    save_csv(pheno$phenotype_sample, "phenotype_sample_features.csv")
  }
  
  save_csv(merged$sample_df, "analysis_dataset_sample.csv")
  save_csv(merged$sample_ab_df, "analysis_dataset_sample_antibiotic.csv")
  
  save_csv(summaries$overall, "summary_overall.csv")
  save_csv(summaries$by_group, "summary_by_genotype_group.csv")
  save_csv(summaries$by_antibiotic, "summary_by_antibiotic.csv")
  
  if (!is.null(sample_models$tidy)) {
    save_csv(sample_models$tidy, "sample_level_model_tidy.csv")
  }
  if (!is.null(sample_models$glance)) {
    save_csv(sample_models$glance, "sample_level_model_glance.csv")
  }
  
  if (!is.null(ab_models$tidy)) {
    save_csv(ab_models$tidy, "antibiotic_specific_model_tidy.csv")
  }
  if (!is.null(ab_models$glance)) {
    save_csv(ab_models$glance, "antibiotic_specific_model_glance.csv")
  }
  
  message2("Creating plots ...")
  sample_plots <- plot_sample_level(merged$sample_df)
  save_plot(sample_plots$p1, "sample_accuracy_by_genotype_group.png", 8, 5)
  save_plot(sample_plots$p2, "sample_accuracy_vs_gene_burden.png", 7, 5)
  save_plot(sample_plots$p3, "sample_ambiguity_by_genotype_group.png", 8, 5)
  save_plot(sample_plots$p4, "sample_VME_by_genotype_group.png", 8, 5)
  
  ab_plot_obj <- plot_antibiotic_group_summary(merged$sample_ab_df)
  save_csv(ab_plot_obj$summary_tbl, "antibiotic_group_summary_table.csv")
  save_plot(ab_plot_obj$plot, "antibiotic_specific_accuracy_by_genotype_group.png", 10, 7)
  
  heat_obj <- plot_heatmaps(merged$sample_ab_df)
  save_csv(heat_obj$heat_beta_ab_tbl, "heatmap_genotype_group_vs_antibiotic.csv")
  save_csv(heat_obj$heat_sample_ab_tbl, "heatmap_sample_vs_antibiotic.csv")
  save_plot(heat_obj$p1, "heatmap_genotype_group_vs_antibiotic.png", 9, 5)
  save_plot(heat_obj$p2, "heatmap_sample_vs_antibiotic.png", 8, 9)
  
  if (!is.null(sample_models$tidy) && nrow(sample_models$tidy) > 0) {
    p_coef <- plot_model_coefficients(sample_models$tidy, "Sample-level model coefficients")
    save_plot(p_coef, "sample_level_model_coefficients.png", 10, 6)
  }
  
  if (!is.null(ab_models$tidy) && nrow(ab_models$tidy) > 0) {
    p_ab_coef <- plot_ab_coefficients(ab_models$tidy)
    save_plot(p_ab_coef, "antibiotic_specific_model_coefficients.png", 10, 7)
  }
  
  message2("Clustering isolates by prediction pattern ...")
  clust <- cluster_isolates_by_prediction_pattern(merged$sample_ab_df)
  if (!is.null(clust)) {
    save_csv(clust$clusters, "prediction_pattern_clusters.csv")
    p_cluster <- plot_cluster_heatmap(merged$sample_ab_df, clust$clusters)
    save_plot(p_cluster, "clustered_prediction_pattern_heatmap.png", 9, 9)
    
    png(file.path(OUTDIR, "plots", "prediction_pattern_dendrogram.png"), width = 1200, height = 800, res = 150)
    plot(clust$hclust, main = "Dendrogram: isolates clustered by prediction pattern")
    dev.off()
  }
  
  message2("Writing model summaries ...")
  if (!is.null(sample_models$fits) && length(sample_models$fits) > 0) {
    sink(file.path(OUTDIR, "models", "sample_level_models_summary.txt"))
    for (x in sample_models$fits) {
      cat("\n============================================================\n")
      cat("MODEL:", x$name, "\n")
      cat("============================================================\n")
      print(summary(x$fit))
      cat("\n")
    }
    sink()
  }
  
  if (!is.null(ab_models$fits) && length(ab_models$fits) > 0) {
    sink(file.path(OUTDIR, "models", "antibiotic_specific_models_summary.txt"))
    for (x in ab_models$fits) {
      cat("\n============================================================\n")
      cat("ANTIBIOTIC:", x$antibiotic, "\n")
      cat("N:", x$n, "\n")
      cat("============================================================\n")
      print(summary(x$fit))
      cat("\n")
    }
    sink()
  }
  
  message2("Done.")
  invisible(list(
    prediction_long = pred_long,
    perf = perf,
    geno = geno,
    pheno = pheno,
    merged = merged,
    sample_models = sample_models,
    ab_models = ab_models,
    summaries = summaries
  ))
}

# ----------------------------------------------------------------------------
# OPTIONAL QUICK INTERPRETATION TABLES
# ----------------------------------------------------------------------------

make_quick_interpretation_tables <- function(sample_df, sample_model_tidy = NULL, ab_model_tidy = NULL) {
  top_hard_groups <- sample_df %>%
    group_by(beta_complexity) %>%
    summarise(
      n = n(),
      mean_acc = mean(acc_total, na.rm = TRUE),
      mean_ambiguity = mean(ambiguity_total, na.rm = TRUE),
      mean_VME = mean(VME_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_acc)
  
  save_csv(top_hard_groups, "quick_interpretation_hardest_groups.csv")
  
  if (!is.null(sample_model_tidy)) {
    top_effects <- sample_model_tidy %>%
      filter(term != "(Intercept)") %>%
      arrange(p.value)
    
    save_csv(top_effects, "quick_interpretation_sample_model_effects.csv")
  }
  
  if (!is.null(ab_model_tidy)) {
    top_ab_effects <- ab_model_tidy %>%
      filter(term != "(Intercept)") %>%
      arrange(p.value)
    
    save_csv(top_ab_effects, "quick_interpretation_antibiotic_model_effects.csv")
  }
  
  invisible(TRUE)
}

# ----------------------------------------------------------------------------
# RUN
# ----------------------------------------------------------------------------
results <- run_prediction_vs_genotype_analysis()
make_quick_interpretation_tables(
  results$merged$sample_df,
  results$sample_models$tidy,
  results$ab_models$tidy
)