`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(tibble)
})

source(file = "model/modelCommon.R")
source(file = "genotype/genotypeCommon.R")

indir  <- file.path(processedRootRassembly, "genotype")
outdir <- file.path(processedRootRassembly, "phenotype")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 1) Inputs
# =========================================================
gene_calls_file <- file.path(indir, "resolved_amr_gene_calls_long_core.csv")
gene_table_file <- file.path(indir, "resolved_amr_genes_manuscript_table.csv")
sample_group_file <- file.path(indir, "sample_genotype_grouping_all_classes.csv")

stopifnot(file.exists(gene_calls_file))
stopifnot(file.exists(gene_table_file))

has_sample_grouping <- file.exists(sample_group_file)

# =========================================================
# 2) Helpers
# =========================================================
na_if_empty <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA")] <- NA_character_
  x
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(stats::na.omit(as.character(x)))
  x <- x[x != ""]
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = sep)
}

normalize_sample <- function(x) {
  normalize_sample_id(as.character(x))
}

sample_order_key <- function(x) {
  suppressWarnings(as.integer(normalize_sample(x)))
}

severity_to_sur <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x <= 0 ~ "S",
    x == 1 ~ "U",
    x >= 2 ~ "R"
  )
}

coalesce_chr <- function(x, y) {
  dplyr::coalesce(as.character(x), as.character(y))
}

normalize_sample_column <- function(df, candidates = c("sample", "sample_id", "Sample", "SampleID")) {
  nm <- intersect(candidates, names(df))
  if (length(nm) == 0) {
    stop("Could not find a sample identifier column.")
  }
  df %>% mutate(sample = normalize_sample(.data[[nm[1]]]))
}

write_csv2_safe <- function(x, path) {
  write.csv2(x, file = path, row.names = FALSE, na = "")
}

add_subclass_rule <- function(pattern, antibiotics, severity, note, confidence = "specific") {
  tibble(
    subclass_pattern = pattern,
    Antibiotic = antibiotics,
    severity = severity,
    rule_note = note,
    confidence = confidence
  )
}

# =========================================================
# 3) Antibiotic universe
# =========================================================
MODEL_ANTIBIOTICS_ORDER <- if (exists("ANTIBIOTICS", inherits = TRUE)) {
  get("ANTIBIOTICS", inherits = TRUE)
} else if (exists("ALL_ANTIBIOTICS_IN_MODEL", inherits = TRUE)) {
  get("ALL_ANTIBIOTICS_IN_MODEL", inherits = TRUE)
} else {
  stop("Could not find ANTIBIOTICS or ALL_ANTIBIOTICS_IN_MODEL in sourced files.")
}

TARGET_ANTIBIOTICS <- MODEL_ANTIBIOTICS_ORDER[
  MODEL_ANTIBIOTICS_ORDER %in% c(
    "AMP", "AMC", "PIP", "TZP", "CAZ", "CRO", "CTX", "FEP",
    "CIP", "OFX", "LVX", "MFX", "GEN", "TOB"
  )
]

quinolone_antibiotics <- c("CIP", "OFX", "LVX", "MFX")

# =========================================================
# 4) Load and normalize data
# =========================================================
gene_calls <- read.csv2(gene_calls_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    sample_id = normalize_sample(sample_id),
    scope = if ("scope" %in% names(.)) na_if_empty(scope) else NA_character_,
    resolved_gene = dplyr::coalesce(
      if ("resolved_gene" %in% names(.)) na_if_empty(resolved_gene) else NA_character_,
      if ("Gene_symbol" %in% names(.)) na_if_empty(Gene_symbol) else NA_character_
    )
  )

if ("scope" %in% names(gene_calls)) {
  gene_calls_core <- gene_calls %>% filter(is.na(scope) | scope == "core")
} else {
  gene_calls_core <- gene_calls
}

gross_gene_table <- read.csv2(gene_table_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    Gene_symbol = na_if_empty(Gene_symbol),
    class = na_if_empty(class),
    subclass = na_if_empty(subclass),
    Enzyme_family = na_if_empty(Enzyme_family),
    Molecular_class = na_if_empty(Molecular_class),
    Functional_group = na_if_empty(Functional_group),
    reference = na_if_empty(reference)
  )

all_samples <- tibble(sample = normalize_sample(allSamples())) %>%
  distinct() %>%
  filter(!is.na(sample), sample != "")

# =========================================================
# 5) Distinct detected genes per sample + annotation
# =========================================================
gene_detected <- gene_calls_core %>%
  transmute(
    sample = normalize_sample(sample_id),
    Gene_symbol = resolved_gene
  ) %>%
  filter(!is.na(sample), sample != "", !is.na(Gene_symbol)) %>%
  distinct()

gene_annot <- gene_detected %>%
  left_join(gross_gene_table, by = "Gene_symbol") %>%
  mutate(
    class = case_when(
      !is.na(class) ~ class,
      str_detect(Gene_symbol, regex("^(bla|ampC)", ignore_case = TRUE)) ~ "BETA-LACTAM",
      TRUE ~ NA_character_
    ),
    Enzyme_family = case_when(
      !is.na(Enzyme_family) ~ Enzyme_family,
      str_detect(Gene_symbol, regex("^blaCTX-M", ignore_case = TRUE)) ~ "CTX-M",
      str_detect(Gene_symbol, regex("^blaOXA", ignore_case = TRUE)) ~ "OXA",
      str_detect(Gene_symbol, regex("^blaTEM", ignore_case = TRUE)) ~ "TEM",
      str_detect(Gene_symbol, regex("^blaSHV", ignore_case = TRUE)) ~ "SHV",
      str_detect(Gene_symbol, regex("^blaNDM", ignore_case = TRUE)) ~ "NDM",
      str_detect(Gene_symbol, regex("^blaDHA", ignore_case = TRUE)) ~ "AMPC",
      str_detect(Gene_symbol, regex("^blaCMY", ignore_case = TRUE)) ~ "AMPC",
      str_detect(Gene_symbol, regex("^blaACC", ignore_case = TRUE)) ~ "AMPC",
      str_detect(Gene_symbol, regex("^blaFOX", ignore_case = TRUE)) ~ "AMPC",
      str_detect(Gene_symbol, regex("^blaMOX", ignore_case = TRUE)) ~ "AMPC",
      str_detect(Gene_symbol, regex("^blaACT", ignore_case = TRUE)) ~ "AMPC",
      str_detect(Gene_symbol, regex("^blaMIR", ignore_case = TRUE)) ~ "AMPC",
      str_detect(Gene_symbol, regex("^ampC", ignore_case = TRUE)) ~ "AMPC",
      TRUE ~ Enzyme_family
    ),
    Functional_group = case_when(
      !is.na(Functional_group) ~ Functional_group,
      str_detect(Gene_symbol, regex("^blaCTX-M", ignore_case = TRUE)) ~ "2be",
      str_detect(Gene_symbol, regex("^blaOXA", ignore_case = TRUE)) ~ "2d",
      str_detect(Gene_symbol, regex("^blaDHA|^blaCMY|^blaACC|^blaFOX|^blaMOX|^blaACT|^blaMIR", ignore_case = TRUE)) ~ "ampc_beta_lactamase",
      str_detect(Gene_symbol, regex("^ampC", ignore_case = TRUE)) ~ "ampc_hyperproduction",
      str_detect(Gene_symbol, regex("^blaNDM", ignore_case = TRUE)) ~ "3a",
      TRUE ~ Functional_group
    ),
    beta_rule_key = case_when(
      str_detect(Gene_symbol, regex("^blaCTX-M-(15|16|28|55|3|57|82)", ignore_case = TRUE)) ~ "2be_CTXM15_like",
      str_detect(Gene_symbol, regex("^blaCTX-M-(14|27|24|17|65)", ignore_case = TRUE)) ~ "2be_CTXM14_like",
      str_detect(Gene_symbol, regex("^blaOXA-1$|^blaOXA-30$", ignore_case = TRUE)) ~ "2d_OXA1_like",
      str_detect(Gene_symbol, regex("^blaDHA|^blaCMY|^blaACC|^blaFOX|^blaMOX|^blaACT|^blaMIR", ignore_case = TRUE)) ~ "ampc_beta_lactamase",
      str_detect(Gene_symbol, regex("^ampC", ignore_case = TRUE)) ~ "ampc_hyperproduction",
      TRUE ~ Functional_group
    )
  )

# =========================================================
# 6) Rule tables
#    severity: 0 = S, 1 = U, 2 = R
# =========================================================
beta_rules <- tribble(
  ~beta_rule_key,      ~Antibiotic, ~severity, ~rule_note,
  "ampc_beta_lactamase",     "AMP", 2, "AmpC beta-lactamase",
  "ampc_beta_lactamase",     "AMC", 2, "AmpC beta-lactamase not inhibited by clavulanate",
  "ampc_beta_lactamase",     "PIP", 2, "AmpC beta-lactamase",
  "ampc_beta_lactamase",     "TZP", 1, "AmpC beta-lactamase with reduced inhibitor reliability",
  "ampc_beta_lactamase",     "CAZ", 2, "AmpC beta-lactamase",
  "ampc_beta_lactamase",     "CRO", 2, "AmpC beta-lactamase",
  "ampc_beta_lactamase",     "CTX", 2, "AmpC beta-lactamase",
  "ampc_beta_lactamase",     "FEP", 1, "AmpC beta-lactamase may spare cefepime partly; keep conservative",
  
  "ampc_hyperproduction",    "AMP", 2, "Chromosomal ampC hyperproduction",
  "ampc_hyperproduction",    "AMC", 1, "Chromosomal ampC hyperproduction not inhibited by clavulanate",
  "ampc_hyperproduction",    "PIP", 1, "Chromosomal ampC hyperproduction",
  "ampc_hyperproduction",    "TZP", 1, "Chromosomal ampC hyperproduction with reduced inhibitor reliability",
  "ampc_hyperproduction",    "CAZ", 1, "Chromosomal ampC hyperproduction",
  "ampc_hyperproduction",    "CRO", 1, "Chromosomal ampC hyperproduction",
  "ampc_hyperproduction",    "CTX", 1, "Chromosomal ampC hyperproduction",
  "ampc_hyperproduction",    "FEP", 1, "Chromosomal ampC hyperproduction may spare cefepime partly; keep conservative",
  
  "2b",               "AMP", 2, "Broad-spectrum penicillinase",
  "2b",               "AMC", 1, "Default to U for inhibitor combinations",
  "2b",               "PIP", 2, "Broad-spectrum penicillinase",
  "2b",               "TZP", 1, "Default to U for inhibitor combinations",
  "2b",               "CAZ", 0, "No extended-spectrum activity",
  "2b",               "CRO", 0, "No extended-spectrum activity",
  "2b",               "CTX", 0, "No extended-spectrum activity",
  "2b",               "FEP", 0, "No extended-spectrum activity",
  
  "2be",              "AMP", 2, "ESBL",
  "2be",              "AMC", 1, "ESBL with variable inhibitor effect",
  "2be",              "PIP", 2, "ESBL",
  "2be",              "TZP", 1, "ESBL with variable inhibitor effect",
  "2be",              "CAZ", 1, "Generic ESBL; CAZ depends on subtype/background",
  "2be",              "CRO", 2, "ESBL",
  "2be",              "CTX", 2, "ESBL",
  "2be",              "FEP", 1, "Keep cefepime conservative in generic ESBL background",
  
  "2be_CTXM15_like",  "AMP", 2, "CTX-M-15-like ESBL",
  "2be_CTXM15_like",  "AMC", 1, "CTX-M-15-like ESBL with variable inhibitor effect",
  "2be_CTXM15_like",  "PIP", 2, "CTX-M-15-like ESBL",
  "2be_CTXM15_like",  "TZP", 1, "CTX-M-15-like; keep TZP conservative without OXA/permeability data",
  "2be_CTXM15_like",  "CAZ", 2, "CTX-M-15-like may impact ceftazidime more strongly",
  "2be_CTXM15_like",  "CRO", 2, "CTX-M-15-like ESBL",
  "2be_CTXM15_like",  "CTX", 2, "CTX-M-15-like ESBL",
  "2be_CTXM15_like",  "FEP", 1, "Cefepime still kept conservative",
  
  "2be_CTXM14_like",  "AMP", 2, "CTX-M-14/27-like ESBL",
  "2be_CTXM14_like",  "AMC", 1, "CTX-M-14/27-like ESBL with variable inhibitor effect",
  "2be_CTXM14_like",  "PIP", 2, "CTX-M-14/27-like ESBL",
  "2be_CTXM14_like",  "TZP", 1, "Keep TZP conservative without OXA/permeability data",
  "2be_CTXM14_like",  "CAZ", 1, "CTX-M-14/27-like may spare CAZ more often than CTX-M-15-like",
  "2be_CTXM14_like",  "CRO", 2, "CTX-M-14/27-like ESBL",
  "2be_CTXM14_like",  "CTX", 2, "CTX-M-14/27-like ESBL",
  "2be_CTXM14_like",  "FEP", 1, "Cefepime kept conservative",
  
  "2br",              "AMP", 2, "Inhibitor-resistant penicillinase",
  "2br",              "AMC", 2, "Inhibitor-resistant",
  "2br",              "PIP", 2, "Inhibitor-resistant penicillinase",
  "2br",              "TZP", 2, "Inhibitor-resistant",
  "2br",              "CAZ", 0, "No ESBL activity",
  "2br",              "CRO", 0, "No ESBL activity",
  "2br",              "CTX", 0, "No ESBL activity",
  "2br",              "FEP", 0, "No ESBL activity",
  
  "2ber",             "AMP", 2, "ESBL plus inhibitor-resistant",
  "2ber",             "AMC", 2, "ESBL plus inhibitor-resistant",
  "2ber",             "PIP", 2, "ESBL plus inhibitor-resistant",
  "2ber",             "TZP", 2, "ESBL plus inhibitor-resistant",
  "2ber",             "CAZ", 1, "Generic ESBL plus inhibitor-resistant background; keep CAZ conservative",
  "2ber",             "CRO", 2, "ESBL",
  "2ber",             "CTX", 2, "ESBL",
  "2ber",             "FEP", 1, "Possible reduced susceptibility",
  
  "2d",               "AMP", 2, "OXA",
  "2d",               "AMC", 1, "Generic OXA; inhibitor effect variable",
  "2d",               "PIP", 2, "OXA",
  "2d",               "TZP", 1, "Generic OXA; keep inhibitor combinations conservative",
  "2d",               "CAZ", 0, "Usually not ESBL by default",
  "2d",               "CRO", 0, "Usually not ESBL by default",
  "2d",               "CTX", 0, "Usually not ESBL by default",
  "2d",               "FEP", 0, "Usually not ESBL by default",
  
  "2d_OXA1_like",     "AMP", 2, "OXA-1-like narrow-spectrum oxacillinase",
  "2d_OXA1_like",     "AMC", 2, "OXA-1-like often compromises clavulanate combinations",
  "2d_OXA1_like",     "PIP", 2, "OXA-1-like",
  "2d_OXA1_like",     "TZP", 2, "OXA-1-like often compromises TZP",
  "2d_OXA1_like",     "CAZ", 0, "No ESBL activity by default",
  "2d_OXA1_like",     "CRO", 0, "No ESBL activity by default",
  "2d_OXA1_like",     "CTX", 0, "No ESBL activity by default",
  "2d_OXA1_like",     "FEP", 0, "No ESBL activity by default",
  
  "3a",               "AMP", 2, "Metallo-beta-lactamase",
  "3a",               "AMC", 2, "Metallo-beta-lactamase",
  "3a",               "PIP", 2, "Metallo-beta-lactamase",
  "3a",               "TZP", 2, "Metallo-beta-lactamase",
  "3a",               "CAZ", 2, "Metallo-beta-lactamase",
  "3a",               "CRO", 2, "Metallo-beta-lactamase",
  "3a",               "CTX", 2, "Metallo-beta-lactamase",
  "3a",               "FEP", 2, "Metallo-beta-lactamase"
) %>%
  filter(Antibiotic %in% TARGET_ANTIBIOTICS)

beta_gene_effects <- gene_annot %>%
  filter(str_detect(class %||% "", regex("BETA-LACTAM", ignore_case = TRUE)) |
           str_detect(Gene_symbol, regex("^(bla|ampC)", ignore_case = TRUE))) %>%
  filter(!is.na(beta_rule_key)) %>%
  left_join(beta_rules, by = "beta_rule_key", relationship = "many-to-many") %>%
  mutate(
    rule_class = "beta_lactam",
    Functional_group_original = Functional_group,
    Functional_group = beta_rule_key
  )

beta_unmapped <- gene_annot %>%
  filter(str_detect(class %||% "", regex("BETA-LACTAM", ignore_case = TRUE)) |
           str_detect(Gene_symbol, regex("^(bla|ampC)", ignore_case = TRUE))) %>%
  filter(is.na(beta_rule_key)) %>%
  arrange(sample_order_key(sample), sample, Gene_symbol)

quinolone_gene_catalog <- tribble(
  ~gene_pattern,             ~gene_family, ~mechanism_group, ~strength, ~rule_note,
  "^qnr",                  "qnr",       "PMQR",          1L,        "Plasmid-mediated quinolone protection",
  "^qepA",                 "qepA",      "PMQR",          1L,        "QepA-mediated quinolone efflux",
  "^oqxA$",                "oqxA",      "PMQR",          1L,        "OqxAB-associated quinolone efflux",
  "^oqxB$",                "oqxB",      "PMQR",          1L,        "OqxAB-associated quinolone efflux",
  "^aac\\(6'\\)-Ib-cr", "aac6Ibcr",  "PMQR",          1L,        "aac(6')-Ib-cr may reduce fluoroquinolone susceptibility",
  "^gyrA",                 "gyrA",      "QRDR",          1L,        "QRDR determinant",
  "^gyrB",                 "gyrB",      "QRDR",          1L,        "QRDR determinant",
  "^parC",                 "parC",      "QRDR",          1L,        "QRDR determinant",
  "^parE",                 "parE",      "QRDR",          1L,        "QRDR determinant"
)

quinolone_detected <- purrr::map_dfr(seq_len(nrow(quinolone_gene_catalog)), function(i) {
  rule <- quinolone_gene_catalog[i, ]
  gene_annot %>%
    filter(str_detect(Gene_symbol, regex(rule$gene_pattern, ignore_case = TRUE))) %>%
    transmute(
      sample,
      Gene_symbol,
      class,
      subclass,
      Enzyme_family,
      Functional_group,
      gene_family = rule$gene_family,
      mechanism_group = rule$mechanism_group,
      strength = rule$strength,
      gene_rule_note = rule$rule_note
    )
}) %>%
  distinct()

quinolone_combination_rules <- tribble(
  ~combination_rule,          ~severity, ~rule_note,
  "PMQR_only",              1,         "PMQR only -> U",
  "single_QRDR",            1,         "Single QRDR determinant -> U",
  "multiple_QRDR",          2,         "Multiple QRDR determinants -> R",
  "PMQR_plus_single_QRDR",  2,         "PMQR plus QRDR -> lean R",
  "PMQR_plus_multiple_QRDR",2,         "PMQR plus multiple QRDR determinants -> R"
)

quinolone_sample_summary <- quinolone_detected %>%
  group_by(sample) %>%
  summarise(
    pmqr_genes = collapse_unique(Gene_symbol[mechanism_group == "PMQR"]),
    qrdr_genes = collapse_unique(Gene_symbol[mechanism_group == "QRDR"]),
    pmqr_families = collapse_unique(gene_family[mechanism_group == "PMQR"]),
    qrdr_families = collapse_unique(gene_family[mechanism_group == "QRDR"]),
    n_pmqr_genes = n_distinct(Gene_symbol[mechanism_group == "PMQR"]),
    n_qrdr_genes = n_distinct(Gene_symbol[mechanism_group == "QRDR"]),
    n_qrdr_families = n_distinct(gene_family[mechanism_group == "QRDR"]),
    has_efflux_pmqr = any(gene_family %in% c("qepA", "oqxA", "oqxB")),
    .groups = "drop"
  ) %>%
  mutate(
    has_pmqr = n_pmqr_genes > 0,
    has_qrdr = n_qrdr_genes > 0,
    combination_rule = case_when(
      has_pmqr & (n_qrdr_genes >= 2 | n_qrdr_families >= 2) ~ "PMQR_plus_multiple_QRDR",
      has_pmqr & has_qrdr ~ "PMQR_plus_single_QRDR",
      n_qrdr_genes >= 2 | n_qrdr_families >= 2 ~ "multiple_QRDR",
      has_qrdr ~ "single_QRDR",
      has_pmqr ~ "PMQR_only",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(quinolone_combination_rules, by = "combination_rule")

quinolone_gene_rules <- quinolone_combination_rules %>%
  mutate(
    rule_family = "quinolone_combination",
    Antibiotic = list(quinolone_antibiotics)
  ) %>%
  tidyr::unnest(Antibiotic) %>%
  filter(Antibiotic %in% TARGET_ANTIBIOTICS)

quinolone_gene_effects <- quinolone_sample_summary %>%
  filter(!is.na(combination_rule)) %>%
  tidyr::expand_grid(Antibiotic = quinolone_antibiotics) %>%
  filter(Antibiotic %in% TARGET_ANTIBIOTICS) %>%
  transmute(
    sample,
    Gene_symbol = case_when(
      !is.na(pmqr_genes) & !is.na(qrdr_genes) ~ paste(pmqr_genes, qrdr_genes, sep = "; "),
      !is.na(pmqr_genes) ~ pmqr_genes,
      !is.na(qrdr_genes) ~ qrdr_genes,
      TRUE ~ NA_character_
    ),
    class = "QUINOLONE",
    subclass = combination_rule,
    Enzyme_family = qrdr_families,
    Functional_group = pmqr_families,
    Antibiotic,
    severity,
    rule_note,
    rule_class = "quinolone"
  )

aminoglycoside_subclass_rules <- bind_rows(
  add_subclass_rule("GENTAMICIN", c("GEN"), 2, "Gentamicin-specific AME subclass"),
  add_subclass_rule("GENTAMICIN|TOBRAMYCIN", c("TOB"), 1, "Tobramycin overlap uncertain without exact AME substrate profile", confidence = "conservative"),
  add_subclass_rule("TOBRAMYCIN", c("TOB"), 1, "Tobramycin-specific call kept conservative as U", confidence = "conservative")
) %>%
  filter(Antibiotic %in% TARGET_ANTIBIOTICS)

aminoglycoside_gene_effects <- purrr::map_dfr(seq_len(nrow(aminoglycoside_subclass_rules)), function(i) {
  rule <- aminoglycoside_subclass_rules[i, ]
  gene_annot %>%
    filter(str_detect(class %||% "", regex("AMINOGLYCOSIDE", ignore_case = TRUE))) %>%
    filter(str_detect(subclass %||% "", regex(rule$subclass_pattern, ignore_case = TRUE))) %>%
    transmute(
      sample,
      Gene_symbol,
      class,
      subclass,
      Enzyme_family,
      Functional_group,
      Antibiotic = rule$Antibiotic,
      severity = rule$severity,
      rule_note = rule$rule_note,
      rule_class = "aminoglycoside"
    )
}) %>%
  distinct()

# =========================================================
# 7) Union all effects and aggregate to sample-level SUR
# =========================================================
all_gene_effects <- bind_rows(
  beta_gene_effects %>%
    select(sample, Gene_symbol, class, subclass, Enzyme_family, Functional_group,
           Antibiotic, severity, rule_note, rule_class),
  quinolone_gene_effects,
  aminoglycoside_gene_effects
) %>%
  mutate(
    sample = normalize_sample(sample),
    Antibiotic = as.character(Antibiotic)
  ) %>%
  filter(!is.na(sample), sample != "", Antibiotic %in% TARGET_ANTIBIOTICS)

predictions_long <- all_gene_effects %>%
  group_by(sample, Antibiotic) %>%
  summarise(
    max_severity = max(severity, na.rm = TRUE),
    Predicted_SUR = severity_to_sur(max_severity),
    Contributing_genes = collapse_unique(Gene_symbol[severity == max_severity]),
    Contributing_classes = collapse_unique(rule_class[severity == max_severity]),
    Contributing_groups = collapse_unique(Functional_group[severity == max_severity]),
    Rule_notes = collapse_unique(rule_note[severity == max_severity]),
    .groups = "drop"
  ) %>%
  mutate(
    max_severity = ifelse(is.infinite(max_severity), NA_real_, max_severity),
    Predicted_SUR = ifelse(is.na(max_severity), NA_character_, Predicted_SUR)
  )

predictions_long_complete <- tidyr::expand_grid(
  sample = all_samples$sample,
  Antibiotic = factor(TARGET_ANTIBIOTICS, levels = TARGET_ANTIBIOTICS)
) %>%
  left_join(predictions_long, by = c("sample", "Antibiotic")) %>%
  mutate(
    max_severity = ifelse(is.na(max_severity), 0, max_severity),
    Predicted_SUR = coalesce(Predicted_SUR, "S"),
    Antibiotic = factor(as.character(Antibiotic), levels = TARGET_ANTIBIOTICS)
  ) %>%
  arrange(sample_order_key(sample), sample, Antibiotic)

predictions_wide <- predictions_long_complete %>%
  select(sample, Antibiotic, Predicted_SUR) %>%
  pivot_wider(names_from = Antibiotic, values_from = Predicted_SUR) %>%
  select(sample, all_of(TARGET_ANTIBIOTICS)) %>%
  arrange(sample_order_key(sample), sample)

sur_cols <- TARGET_ANTIBIOTICS[TARGET_ANTIBIOTICS %in% names(predictions_wide)]

predictions_grouped <- predictions_wide %>%
  unite("sur_pattern", all_of(sur_cols), sep = "|", remove = FALSE) %>%
  mutate(samplegroup_sur = paste0("SUR", sprintf("%02d", dense_rank(sur_pattern)))) %>%
  arrange(sample_order_key(sample), sample)

predictions_groups_summary <- predictions_grouped %>%
  group_by(samplegroup_sur, sur_pattern) %>%
  summarise(
    n_samples = n(),
    samples = paste(sort(sample), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(samplegroup_sur)

# =========================================================
# 8) Join with sample genotype grouping from all classes
# =========================================================
if (has_sample_grouping) {
  sample_grouping <- read.csv2(sample_group_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
    normalize_sample_column() %>%
    distinct(sample, .keep_all = TRUE)
} else {
  sample_grouping <- all_samples
}

combined_per_sample <- all_samples %>%
  left_join(sample_grouping, by = "sample") %>%
  left_join(predictions_grouped, by = "sample") %>%
  arrange(sample_order_key(sample), sample)

# =========================================================
# 9) Crosswalks and summaries
# =========================================================
sur_group_summary <- combined_per_sample %>%
  group_by(samplegroup_sur, sur_pattern) %>%
  summarise(
    n_samples = n(),
    n_beta_groups = if ("samplegroup_functional" %in% names(cur_data())) n_distinct(samplegroup_functional, na.rm = TRUE) else NA_integer_,
    n_quinolone_groups = if ("samplegroup_quinolone" %in% names(cur_data())) n_distinct(samplegroup_quinolone, na.rm = TRUE) else NA_integer_,
    n_aminoglycoside_groups = if ("samplegroup_aminoglycoside" %in% names(cur_data())) n_distinct(samplegroup_aminoglycoside, na.rm = TRUE) else NA_integer_,
    functional_groupsets = if ("Functional_groups" %in% names(cur_data())) collapse_unique(Functional_groups) else NA_character_,
    quinolone_family_sets = if ("Quinolone_families" %in% names(cur_data())) collapse_unique(Quinolone_families) else NA_character_,
    aminoglycoside_family_sets = if ("Aminoglycoside_families" %in% names(cur_data())) collapse_unique(Aminoglycoside_families) else NA_character_,
    gene_examples = collapse_unique(c(
      if ("Beta_genes" %in% names(cur_data())) Beta_genes else NULL,
      if ("Quinolone_genes" %in% names(cur_data())) Quinolone_genes else NULL,
      if ("Aminoglycoside_genes" %in% names(cur_data())) Aminoglycoside_genes else NULL
    )),
    samples = paste(sort(sample), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(samplegroup_sur)

group_crosswalk <- combined_per_sample %>%
  select(any_of(c(
    "sample", "samplegroup_sur", "sur_pattern",
    "samplegroup_functional", "samplegroup_quinolone", "samplegroup_aminoglycoside",
    "Functional_groups", "Quinolone_families", "Aminoglycoside_families"
  ))) %>%
  count(across(-sample), name = "n_samples") %>%
  arrange(samplegroup_sur)

# =========================================================
# 10) Dump rule tables and realized rule hits
# =========================================================
beta_rules_export <- beta_rules %>%
  mutate(
    rule_family = "beta_lactam",
    Predicted_SUR = severity_to_sur(severity)
  ) %>%
  select(rule_family, beta_rule_key, Antibiotic, severity, Predicted_SUR, rule_note) %>%
  arrange(beta_rule_key, Antibiotic)

quinolone_rules_export <- quinolone_gene_rules %>%
  mutate(
    Predicted_SUR = severity_to_sur(severity)
  ) %>%
  select(rule_family, combination_rule, Antibiotic, severity, Predicted_SUR, rule_note) %>%
  arrange(combination_rule, Antibiotic)

aminoglycoside_rules_export <- aminoglycoside_subclass_rules %>%
  mutate(
    rule_family = "aminoglycoside",
    Predicted_SUR = severity_to_sur(severity)
  ) %>%
  select(rule_family, subclass_pattern, confidence, Antibiotic, severity, Predicted_SUR, rule_note) %>%
  arrange(subclass_pattern, Antibiotic)

all_rules_export <- bind_rows(
  beta_rules_export %>% rename(rule_key = beta_rule_key),
  quinolone_rules_export %>% rename(rule_key = combination_rule),
  aminoglycoside_rules_export %>% rename(rule_key = subclass_pattern)
) %>%
  arrange(rule_family, rule_key, Antibiotic)

beta_gene_effects_export <- beta_gene_effects %>%
  mutate(Predicted_SUR = severity_to_sur(severity)) %>%
  select(sample, Gene_symbol, class, subclass, Enzyme_family,
         Functional_group, Antibiotic, severity, Predicted_SUR, rule_note, rule_class) %>%
  arrange(sample_order_key(sample), sample, Antibiotic, Gene_symbol)

quinolone_gene_effects_export <- quinolone_gene_effects %>%
  mutate(Predicted_SUR = severity_to_sur(severity)) %>%
  select(sample, Gene_symbol, class, subclass, Enzyme_family,
         Functional_group, Antibiotic, severity, Predicted_SUR, rule_note, rule_class) %>%
  arrange(sample_order_key(sample), sample, Antibiotic)

aminoglycoside_gene_effects_export <- aminoglycoside_gene_effects %>%
  mutate(Predicted_SUR = severity_to_sur(severity)) %>%
  select(sample, Gene_symbol, class, subclass, Enzyme_family,
         Functional_group, Antibiotic, severity, Predicted_SUR, rule_note, rule_class) %>%
  arrange(sample_order_key(sample), sample, Antibiotic, Gene_symbol)

all_gene_effects_export <- bind_rows(
  beta_gene_effects_export,
  quinolone_gene_effects_export,
  aminoglycoside_gene_effects_export
) %>%
  arrange(rule_class, sample_order_key(sample), sample, Antibiotic, Gene_symbol)

# =========================================================
# 11) Write outputs
# =========================================================
write_csv2_safe(predictions_long_complete, file.path(outdir, "genotype_SUR_long.csv"))
write_csv2_safe(predictions_wide, file.path(outdir, "genotype_SUR_wide.csv"))
write_csv2_safe(predictions_grouped, file.path(outdir, "genotype_SUR_grouped.csv"))
write_csv2_safe(predictions_groups_summary, file.path(outdir, "genotype_SUR_groups_summary.csv"))
write_csv2_safe(combined_per_sample, file.path(outdir, "genotype_SUR_combined_per_sample.csv"))
write_csv2_safe(sur_group_summary, file.path(outdir, "genotype_SUR_group_summary.csv"))
write_csv2_safe(group_crosswalk, file.path(outdir, "genotype_SUR_group_crosswalk.csv"))
write_csv2_safe(beta_unmapped, file.path(outdir, "genotype_SUR_beta_lactamases_missing_functional_group.csv"))

write_csv2_safe(beta_rules_export, file.path(outdir, "genotype_SUR_rules_beta_lactam.csv"))
write_csv2_safe(quinolone_rules_export, file.path(outdir, "genotype_SUR_rules_quinolone.csv"))
write_csv2_safe(aminoglycoside_rules_export, file.path(outdir, "genotype_SUR_rules_aminoglycoside.csv"))
write_csv2_safe(all_rules_export, file.path(outdir, "genotype_SUR_rules_all.csv"))

write_csv2_safe(beta_gene_effects_export, file.path(outdir, "genotype_SUR_rule_hits_beta_lactam.csv"))
write_csv2_safe(quinolone_gene_effects_export, file.path(outdir, "genotype_SUR_rule_hits_quinolone.csv"))
write_csv2_safe(aminoglycoside_gene_effects_export, file.path(outdir, "genotype_SUR_rule_hits_aminoglycoside.csv"))
write_csv2_safe(all_gene_effects_export, file.path(outdir, "genotype_SUR_rule_hits_all.csv"))

message("Done.")
message("Wrote:")
message(" - genotype_SUR_long.csv")
message(" - genotype_SUR_wide.csv")
message(" - genotype_SUR_grouped.csv")
message(" - genotype_SUR_groups_summary.csv")
message(" - genotype_SUR_combined_per_sample.csv")
message(" - genotype_SUR_group_summary.csv")
message(" - genotype_SUR_group_crosswalk.csv")
message(" - genotype_SUR_beta_lactamases_missing_functional_group.csv")
message(" - genotype_SUR_rules_beta_lactam.csv")
message(" - genotype_SUR_rules_quinolone.csv")
message(" - genotype_SUR_rules_aminoglycoside.csv")
message(" - genotype_SUR_rules_all.csv")
message(" - genotype_SUR_rule_hits_beta_lactam.csv")
message(" - genotype_SUR_rule_hits_quinolone.csv")
message(" - genotype_SUR_rule_hits_aminoglycoside.csv")
message(" - genotype_SUR_rule_hits_all.csv")
