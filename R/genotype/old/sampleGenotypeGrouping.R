library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)

source(file = "model/modelCommon.R")
source(file = "genotype/genotypeCommon.R")
outdir <- paste(processedRootRassembly, "genotype", sep = "/")

ALL_ANTIBIOTICS_IN_MODEL <- c(
   "AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP",
   "CIP","OFX","LVX","MFX","GEN","TOB"
 )

# =========================================================
# 1) Load data
# =========================================================
gene_calls <- read.csv2(
  file.path(outdir, "resolved_amr_gene_calls_long_core.csv"),
  stringsAsFactors = FALSE
)

betalactamases <- read.csv2(
  file.path(outdir, "resolved_amr_genes_manuscript_table.csv"),
  stringsAsFactors = FALSE
)

if ("sample_id" %in% names(gene_calls)) {
  gene_calls <- gene_calls %>%
    mutate(sample_id = normalize_sample_id(sample_id))
}

if ("scope" %in% names(gene_calls)) {
  gene_calls$scope <- as.character(gene_calls$scope)
  gene_calls$scope[gene_calls$scope %in% c("", "NA")] <- NA_character_
}

# =========================================================
# 2) Helpers
# =========================================================
na_if_empty <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA")] <- NA_character_
  x
}

normalize_sample <- function(x) {
  normalize_sample_id(x)
}


collapse_unique <- function(x, sep = "; ") {
  x <- unique(stats::na.omit(as.character(x)))
  if (length(x) == 0) NA_character_ else paste(sort(x), collapse = sep)
}


sample_order_key <- function(x) {
  suppressWarnings(as.numeric(str_extract(x, "^[0-9]+")))
}

# normalize sample names once, upstream
if ("sample_id" %in% names(gene_calls)) {
  gene_calls <- gene_calls %>%
    mutate(sample_id = normalize_sample(sample_id))
}


# =========================================================
# 4) Build per-sample beta-lactamase annotation
# =========================================================
beta_calls <- gene_calls %>%
  mutate(
    sample_id = normalize_sample(sample_id),
    resolved_gene = na_if_empty(resolved_gene)
  ) %>%
  filter(!is.na(sample_id), !is.na(resolved_gene)) %>%
  filter(str_detect(resolved_gene, regex("^(bla|ampC)", ignore_case = TRUE))) %>%
  distinct(sample_id, resolved_gene)

beta_annot <- beta_calls %>%
  left_join(
    betalactamases %>%
      transmute(
        Gene_symbol = na_if_empty(Gene_symbol),
        Enzyme_family = na_if_empty(Enzyme_family),
        Molecular_class = na_if_empty(Molecular_class),
        Functional_group = na_if_empty(Functional_group),
        reference = na_if_empty(reference)
      ),
    by = c("resolved_gene" = "Gene_symbol")
  ) %>%
  rename(Gene_symbol = resolved_gene)



# =========================================================
# 5) Rule table based on Bush-Jacoby functional group
# 0 = S, 1 = U, 2 = R
# =========================================================
beta_rules <- tribble(
  ~Functional_group, ~Antibiotic, ~severity, ~rule_note,
  
  "1",    "AMP", 2, "AmpC",
  "1",    "AMC", 2, "AmpC not inhibited by clavulanate",
  "1",    "PIP", 2, "AmpC",
  "1",    "TZP", 1, "AmpC may reduce piperacillin-tazobactam susceptibility",
  "1",    "CAZ", 2, "AmpC",
  "1",    "CRO", 2, "AmpC",
  "1",    "CTX", 2, "AmpC",
  "1",    "FEP", 1, "AmpC may spare cefepime partly",
  
  "2b",   "AMP", 2, "Broad-spectrum penicillinase",
  "2b",   "AMC", 0, "Usually inhibitor-susceptible",
  "2b",   "PIP", 2, "Broad-spectrum penicillinase",
  "2b",   "TZP", 0, "Usually inhibitor-susceptible",
  "2b",   "CAZ", 0, "No extended-spectrum activity",
  "2b",   "CRO", 0, "No extended-spectrum activity",
  "2b",   "CTX", 0, "No extended-spectrum activity",
  "2b",   "FEP", 0, "No extended-spectrum activity",
  
  "2be",  "AMP", 2, "ESBL",
  "2be",  "AMC", 1, "ESBL with variable inhibitor effect",
  "2be",  "PIP", 2, "ESBL",
  "2be",  "TZP", 1, "ESBL with variable inhibitor effect",
  "2be",  "CAZ", 2, "ESBL",
  "2be",  "CRO", 2, "ESBL",
  "2be",  "CTX", 2, "ESBL",
  "2be",  "FEP", 1, "Possible reduced susceptibility",
  
  "2br",  "AMP", 2, "Inhibitor-resistant penicillinase",
  "2br",  "AMC", 2, "Inhibitor-resistant",
  "2br",  "PIP", 2, "Inhibitor-resistant penicillinase",
  "2br",  "TZP", 2, "Inhibitor-resistant",
  "2br",  "CAZ", 0, "No ESBL activity",
  "2br",  "CRO", 0, "No ESBL activity",
  "2br",  "CTX", 0, "No ESBL activity",
  "2br",  "FEP", 0, "No ESBL activity",
  
  "2ber", "AMP", 2, "ESBL plus inhibitor-resistant",
  "2ber", "AMC", 2, "ESBL plus inhibitor-resistant",
  "2ber", "PIP", 2, "ESBL plus inhibitor-resistant",
  "2ber", "TZP", 2, "ESBL plus inhibitor-resistant",
  "2ber", "CAZ", 2, "ESBL",
  "2ber", "CRO", 2, "ESBL",
  "2ber", "CTX", 2, "ESBL",
  "2ber", "FEP", 1, "Possible reduced susceptibility",
  
  "2d",   "AMP", 2, "OXA",
  "2d",   "AMC", 2, "OXA with poor inhibitor effect",
  "2d",   "PIP", 2, "OXA",
  "2d",   "TZP", 1, "Variable OXA effect",
  "2d",   "CAZ", 0, "Usually not ESBL by default",
  "2d",   "CRO", 0, "Usually not ESBL by default",
  "2d",   "CTX", 0, "Usually not ESBL by default",
  "2d",   "FEP", 0, "Usually not ESBL by default",
  
  "3a",   "AMP", 2, "Metallo-beta-lactamase",
  "3a",   "AMC", 2, "Metallo-beta-lactamase",
  "3a",   "PIP", 2, "Metallo-beta-lactamase",
  "3a",   "TZP", 2, "Metallo-beta-lactamase",
  "3a",   "CAZ", 2, "Metallo-beta-lactamase",
  "3a",   "CRO", 2, "Metallo-beta-lactamase",
  "3a",   "CTX", 2, "Metallo-beta-lactamase",
  "3a",   "FEP", 2, "Metallo-beta-lactamase"
) 

# =========================================================
# 6) Join rules to each detected beta-lactamase
# =========================================================
beta_gene_effects <- beta_annot %>%
  filter(!is.na(Functional_group)) %>%
  left_join(beta_rules, by = "Functional_group", relationship = "many-to-many")

beta_unmapped <- beta_annot %>%
  filter(is.na(Functional_group)) %>%
  arrange(sample_id, Gene_symbol)


# =========================================================
# 8) Complete all sample x beta-lactam combinations
# default to S when no beta-lactamase rule contributes
# =========================================================
all_samples <- tibble(sample = allSamples())


# =========================================================
# 10) Per-sample summary
#     Include also samples without beta-lactamases
#     and create samplegroup based on Functional_groups
# =========================================================

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
  ) %>%
  arrange(sample_order_key(sample), sample) %>%
  select(sample, Beta_genes, Functional_groups, Enzyme_families)




# =========================================================
# 10b) Gruppindelning baserat på identiska Functional_groups
# =========================================================
beta_functional_groups_per_sample <- beta_summary_per_sample %>%
  mutate(
    Functional_groups = if_else(is.na(Functional_groups), "", Functional_groups),
    functional_pattern = Functional_groups,
    samplegroup_functional = paste0("FG", sprintf("%02d", dense_rank(functional_pattern)))
  ) %>%
  arrange(sample_order_key(sample), sample) %>%
  select(sample, samplegroup_functional, Functional_groups, Beta_genes, Enzyme_families)

beta_functional_groups_summary <- beta_functional_groups_per_sample %>%
  group_by(samplegroup_functional, Functional_groups) %>%
  summarise(
    n_samples = n(),
    samples = paste(sort(sample), collapse = "; "),
    Functional_groups = first(Functional_groups),
    Beta_genes_examples = paste(unique(stats::na.omit(Beta_genes))[1:min(length(unique(stats::na.omit(Beta_genes))), 5)], collapse = " | "),
    Enzyme_families_examples = paste(unique(stats::na.omit(Enzyme_families))[1:min(length(unique(stats::na.omit(Enzyme_families))), 5)], collapse = " | "),
    .groups = "drop"
  ) %>%
  arrange(samplegroup_functional)



# =========================================================
# 11) Save outputs
# =========================================================
write.csv2(
  beta_summary_per_sample,
  file.path(outdir, "genotype_beta_summary_per_sample.csv"),
  row.names = FALSE
)


write.csv2(
  beta_functional_groups_per_sample,
  file.path(outdir, "genotype_beta_functional_groups_per_sample.csv"),
  row.names = FALSE
)


write.csv2(
  beta_functional_groups_summary,
  file.path(outdir, "genotype_beta_functional_groups_summary.csv"),
  row.names = FALSE
)

print(beta_predictions_wide)
