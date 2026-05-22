source(file = "common.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

# =========================================================
# 0) Paths and input files
# =========================================================
outdir <- file.path(processedRootRassembly, "genotype")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

amrf_file  <- file.path(outdir, "amrfinder_amr_long.csv")
ariba_file <- file.path(outdir, "ariba_genes_long.csv")
bldb_file  <- file.path(outdir, "bldb_lookup.tsv")
bldb_alias_file <- file.path(outdir, "bldb_aliases.tsv")

stopifnot(file.exists(amrf_file))
stopifnot(file.exists(ariba_file))
stopifnot(file.exists(bldb_file))

# =========================================================
# 1) Load outputs from prior scripts
# =========================================================
amrf <- readr::read_csv(amrf_file, show_col_types = FALSE) %>%
  janitor::clean_names()

ariba <- readr::read_csv(ariba_file, show_col_types = FALSE) %>%
  janitor::clean_names()

bldb <- readr::read_tsv(bldb_file, show_col_types = FALSE) %>%
  janitor::clean_names()

bldb_aliases <- if (file.exists(bldb_alias_file)) {
  readr::read_tsv(bldb_alias_file, show_col_types = FALSE) %>%
    janitor::clean_names()
} else {
  tibble()
}

# =========================================================
# 2) Helper functions
# =========================================================
normalize_gene <- function(x) {
  x %>%
    replace_na("") %>%
    str_replace("^bla", "") %>%
    str_replace_all("\\s+", "") %>%
    str_to_upper()
}

get_gene_family <- function(x) {
  case_when(
    str_detect(x, "^blaCTX-M") ~ "blaCTX-M",
    str_detect(x, "^blaTEM")   ~ "blaTEM",
    str_detect(x, "^blaSHV")   ~ "blaSHV",
    str_detect(x, "^blaOXA")   ~ "blaOXA",
    str_detect(x, "^blaKPC")   ~ "blaKPC",
    str_detect(x, "^blaNDM")   ~ "blaNDM",
    str_detect(x, "^blaVIM")   ~ "blaVIM",
    str_detect(x, "^blaIMP")   ~ "blaIMP",
    str_detect(x, "^blaCMY")   ~ "blaCMY",
    str_detect(x, "^blaDHA")   ~ "blaDHA",
    str_detect(x, "^blaFOX")   ~ "blaFOX",
    str_detect(x, "^blaACC")   ~ "blaACC",
    str_detect(x, "^blaACT")   ~ "blaACT",
    str_detect(x, "^blaMIR")   ~ "blaMIR",
    str_detect(x, "^blaLAT")   ~ "blaLAT",
    str_detect(x, "^blaMOX")   ~ "blaMOX",
    str_detect(x, "^blaGES")   ~ "blaGES",
    str_detect(x, "^blaVEB")   ~ "blaVEB",
    str_detect(x, "^blaPER")   ~ "blaPER",
    str_detect(x, "^blaEC")    ~ "blaEC",
    TRUE ~ x
  )
}

infer_beta_group <- function(gene) {
  case_when(
    gene == "blaEC" ~ "Intrinsik blaEC (E. coli klass C)",
    
    str_detect(gene, "^bla(KPC|NDM|VIM|IMP|OXA-48|OXA-181|OXA-232|GES|SME|IMI|NMC)") ~
      "Carbapenemase",
    
    str_detect(gene, "^blaCTX-M") ~
      "ESBL (CTX-M)",
    
    str_detect(gene, "^bla(CMY|DHA|FOX|ACC|ACT|MIR|LAT|MOX)") ~
      "AmpC (förvärvad/plasmid)",
    
    gene == "blaOXA-1" ~
      "OXA-1",
    
    gene %in% c("blaTEM-1", "blaTEM-1A", "blaTEM-1B", "blaTEM-1C", "blaSHV-1") ~
      "Penicillinas (TEM-1/SHV-1)",
    
    str_detect(gene, "^blaTEM") ~ "TEM-familj (övrig)",
    str_detect(gene, "^blaSHV") ~ "SHV-familj (övrig)",
    str_detect(gene, "^blaOXA") ~ "OXA-familj (övrig)",
    str_detect(gene, "^bla")    ~ "Övrig β-laktamas",
    TRUE ~ NA_character_
  )
}

# =========================================================
# 3) Prepare AMRFinder best hits
#    Only core genes
# =========================================================
amrf_best <- amrf %>%
  filter(type == "AMR", scope == "core") %>%
  mutate(
    amrf_gene   = element_symbol,
    amrf_family = get_gene_family(amrf_gene),
    amrf_cov    = suppressWarnings(as.numeric(percent_coverage_of_reference)),
    amrf_ident  = suppressWarnings(as.numeric(percent_identity_to_reference)),
    amrf_method = replace_na(method, ""),
    amrf_partial = str_detect(amrf_method, "PARTIAL") |
      (!is.na(amrf_cov) & amrf_cov < 90),
    
    amrf_rank = case_when(
      !amrf_partial & amrf_method == "ALLELEX" ~ 1,
      !amrf_partial                            ~ 2,
      str_detect(amrf_method, "PARTIAL_CONTIG_END") ~ 3,
      amrf_partial                             ~ 4,
      TRUE                                     ~ 5
    )
  ) %>%
  arrange(sample_id, amrf_family, amrf_rank, desc(amrf_cov), desc(amrf_ident)) %>%
  group_by(sample_id, amrf_family) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    sample_id,
    gene_family = amrf_family,
    amrf_gene,
    amrf_method,
    amrf_cov,
    amrf_ident,
    amrf_partial,
    amrf_element_name = element_name,
    amrf_class = class,
    amrf_subclass = subclass,
    amrf_scope = scope
  )

# =========================================================
# 4) Prepare ARIBA best hits
# =========================================================
ariba_best <- ariba %>%
  mutate(
    ariba_gene   = gene_symbol,
    ariba_family = get_gene_family(ariba_gene),
    assembly_fraction = suppressWarnings(as.numeric(assembly_fraction)),
    pc_ident = suppressWarnings(as.numeric(pc_ident)),
    hit_completeness = replace_na(hit_completeness, "fragment"),
    hit_confidence   = replace_na(hit_confidence, "weak"),
    is_partial       = replace_na(is_partial, TRUE),
    
    ariba_complete = hit_completeness %in% c("complete", "near_complete"),
    
    ariba_rank = case_when(
      hit_completeness == "complete"      ~ 1,
      hit_completeness == "near_complete" ~ 2,
      hit_completeness == "partial"       ~ 3,
      hit_completeness == "fragment"      ~ 4,
      TRUE                                ~ 5
    )
  ) %>%
  arrange(sample_id, ariba_family, ariba_rank, desc(assembly_fraction), desc(pc_ident)) %>%
  group_by(sample_id, ariba_family) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    sample_id,
    gene_family = ariba_family,
    ariba_gene,
    ref_name,
    original_name,
    assembly_fraction,
    pc_ident,
    hit_completeness,
    hit_confidence,
    ariba_partial = is_partial,
    ariba_complete
  )

# =========================================================
# 5) Merge AMRFinder + ARIBA and resolve partial AMRFinder calls
# =========================================================
comparison <- full_join(amrf_best, ariba_best, by = c("sample_id", "gene_family")) %>%
  mutate(
    comparison_status = case_when(
      !is.na(amrf_gene) & !is.na(ariba_gene) & amrf_gene == ariba_gene ~ "exact_match",
      !is.na(amrf_gene) & !is.na(ariba_gene) &
        get_gene_family(amrf_gene) == get_gene_family(ariba_gene) ~ "family_match",
      !is.na(amrf_gene) & is.na(ariba_gene) ~ "amrfinder_only",
      is.na(amrf_gene) & !is.na(ariba_gene) ~ "ariba_only",
      TRUE ~ "discordant"
    ),
    
    ariba_resolves_amrf_partial = case_when(
      isTRUE(amrf_partial) &
        !is.na(ariba_gene) &
        isTRUE(ariba_complete) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    resolved_gene = case_when(
      !is.na(amrf_gene) & !isTRUE(amrf_partial) ~ amrf_gene,
      isTRUE(ariba_resolves_amrf_partial) ~ ariba_gene,
      !is.na(amrf_gene) ~ amrf_gene,
      is.na(amrf_gene) & !is.na(ariba_gene) ~ ariba_gene,
      TRUE ~ NA_character_
    ),
    
    resolved_source = case_when(
      !is.na(amrf_gene) & !isTRUE(amrf_partial) ~ "AMRFinder",
      isTRUE(ariba_resolves_amrf_partial) ~ "ARIBA_resolved_partial",
      !is.na(amrf_gene) ~ "AMRFinder_partial_unresolved",
      is.na(amrf_gene) & !is.na(ariba_gene) ~ "ARIBA_only",
      TRUE ~ NA_character_
    ),
    
    needs_manual_review = case_when(
      comparison_status == "discordant" ~ TRUE,
      isTRUE(amrf_partial) & !isTRUE(ariba_resolves_amrf_partial) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  relocate(sample_id, gene_family, amrf_gene, ariba_gene, resolved_gene, resolved_source)

# =========================================================
# 6) Annotate resolved genes with local BLDB incl. Bush-Jacoby
# =========================================================
bldb_main <- bldb %>%
  mutate(
    gene_norm = normalize_gene(gene_symbol_guess),
    protein_norm = normalize_gene(proteinname)
  )

bldb_alias_tbl <- if (nrow(bldb_aliases) > 0) {
  bldb_aliases %>%
    mutate(
      alias_norm = normalize_gene(alias_bla)
    ) %>%
    select(
      alias_norm,
      proteinname,
      bldb_group,
      molecular_class,
      metallo_subclass,
      subfamily,
      phenotype,
      functionalinformation,
      natural_n_or_acquired_a,
      bush_jacoby_raw,
      bush_jacoby_class,
      bush_jacoby_uncertain
    ) %>%
    distinct()
} else {
  tibble(
    alias_norm = character(),
    proteinname = character(),
    bldb_group = character(),
    molecular_class = character(),
    metallo_subclass = character(),
    subfamily = character(),
    phenotype = character(),
    functionalinformation = character(),
    natural_n_or_acquired_a = character(),
    bush_jacoby_raw = character(),
    bush_jacoby_class = character(),
    bush_jacoby_uncertain = logical()
  )
}

resolved_calls <- comparison %>%
  filter(!is.na(resolved_gene)) %>%
  mutate(
    resolved_gene_norm = normalize_gene(resolved_gene)
  ) %>%
  left_join(
    bldb_main %>%
      select(
        gene_norm,
        proteinname,
        gene_symbol_guess,
        bldb_group,
        molecular_class,
        metallo_subclass,
        subfamily,
        phenotype,
        functionalinformation,
        natural_n_or_acquired_a,
        bush_jacoby_raw,
        bush_jacoby_class,
        bush_jacoby_uncertain
      ),
    by = c("resolved_gene_norm" = "gene_norm")
  ) %>%
  left_join(
    bldb_alias_tbl,
    by = c("resolved_gene_norm" = "alias_norm"),
    suffix = c("", "_alias")
  ) %>%
  mutate(
    bldb_proteinname = coalesce(proteinname, proteinname_alias),
    bldb_group_final = coalesce(bldb_group, bldb_group_alias),
    molecular_class_final = coalesce(molecular_class, molecular_class_alias),
    metallo_subclass_final = coalesce(metallo_subclass, metallo_subclass_alias),
    subfamily_final = coalesce(subfamily, subfamily_alias),
    phenotype_final = coalesce(phenotype, phenotype_alias),
    functionalinformation_final = coalesce(functionalinformation, functionalinformation_alias),
    natural_n_or_acquired_a_final = coalesce(natural_n_or_acquired_a, natural_n_or_acquired_a_alias),
    bush_jacoby_raw_final = coalesce(bush_jacoby_raw, bush_jacoby_raw_alias),
    bush_jacoby_class_final = coalesce(bush_jacoby_class, bush_jacoby_class_alias),
    bush_jacoby_uncertain_final = coalesce(bush_jacoby_uncertain, bush_jacoby_uncertain_alias),
    
    resolved_beta_group = coalesce(
      infer_beta_group(resolved_gene),
      bldb_group_final
    )
  ) %>%
  select(
    sample_id,
    gene_family,
    resolved_gene,
    resolved_source,
    amrf_gene,
    amrf_method,
    amrf_cov,
    amrf_ident,
    amrf_partial,
    ariba_gene,
    ref_name,
    assembly_fraction,
    pc_ident,
    hit_completeness,
    hit_confidence,
    comparison_status,
    needs_manual_review,
    bldb_proteinname,
    bldb_group_final,
    molecular_class_final,
    metallo_subclass_final,
    subfamily_final,
    phenotype_final,
    functionalinformation_final,
    natural_n_or_acquired_a_final,
    bush_jacoby_raw_final,
    bush_jacoby_class_final,
    bush_jacoby_uncertain_final,
    resolved_beta_group
  ) %>%
  arrange(sample_id, gene_family)

# =========================================================
# 7) File focused on AMRFinder partials and their ARIBA rescue
# =========================================================
partials_resolved <- comparison %>%
  filter(amrf_partial %in% TRUE) %>%
  mutate(
    resolution = case_when(
      ariba_resolves_amrf_partial ~ "resolved_by_ariba",
      !is.na(ariba_gene) & !isTRUE(ariba_complete) ~ "ariba_supports_family_but_not_fully_resolved",
      is.na(ariba_gene) ~ "no_ariba_support",
      TRUE ~ "unresolved"
    )
  ) %>%
  arrange(sample_id, gene_family)

# =========================================================
# 8) Summaries from resolved calls
# =========================================================
resolved_beta_hits <- resolved_calls %>%
  filter(str_detect(resolved_gene, "^bla")) %>%
  mutate(
    beta_group = coalesce(resolved_beta_group, infer_beta_group(resolved_gene))
  )

resolved_beta_flags <- resolved_beta_hits %>%
  group_by(sample_id) %>%
  summarise(
    beta_genes_all = paste(sort(unique(resolved_gene)), collapse = "; "),
    beta_groups_all = paste(sort(unique(beta_group)), collapse = "; "),
    bush_jacoby_all = paste(sort(unique(na.omit(bush_jacoby_class_final))), collapse = "; "),
    any_partial_beta_hit = any(amrf_partial %in% TRUE | hit_completeness %in% c("partial", "fragment"), na.rm = TRUE),
    
    has_carbapenemase = any(beta_group == "Carbapenemase" & resolved_gene != "blaEC"),
    has_esbl          = any(str_detect(beta_group, "^ESBL") & resolved_gene != "blaEC"),
    has_ampc          = any(beta_group == "AmpC (förvärvad/plasmid)" & resolved_gene != "blaEC"),
    has_oxa1          = any(beta_group == "OXA-1" & resolved_gene != "blaEC"),
    has_penicillinase = any(beta_group == "Penicillinas \\(TEM-1/SHV-1\\)" & resolved_gene != "blaEC"),
    has_intrinsic_blaEC = any(resolved_gene == "blaEC"),
    
    primary_beta_category = case_when(
      has_carbapenemase ~ "Carbapenemase",
      has_esbl ~ "ESBL",
      has_ampc ~ "AmpC (förvärvad/plasmid)",
      has_oxa1 ~ "OXA-1",
      has_penicillinase ~ "Penicillinas (TEM-1/SHV-1)",
      has_intrinsic_blaEC ~ "Endast intrinsik blaEC",
      TRUE ~ "Inga β-laktamaser hittade"
    ),
    .groups = "drop"
  )

# Non-beta features from AMRFinder retained as primary source
# Only core genes
other_flags <- amrf %>%
  filter(type == "AMR", scope == "core") %>%
  group_by(sample_id) %>%
  summarise(
    has_dfr = any(str_detect(element_symbol, "^dfr")),
    has_sul = any(element_symbol %in% c("sul1", "sul2", "sul3")),
    has_qnr = any(str_detect(element_symbol, "^qnr")),
    has_gyr_par = any(str_detect(element_symbol, "^(gyrA|parC|parE)_")),
    has_gen_tob = any(
      str_detect(element_symbol, "aac\\(3") |
        str_detect(element_symbol, "aac\\(6") |
        str_detect(element_symbol, "aph\\(2") |
        str_detect(element_symbol, "^rmt") |
        str_detect(element_symbol, "^armA$")
    ),
    has_oqx = any(str_detect(element_symbol, "^oqx(A|B)$")),
    .groups = "drop"
  )

resolved_features <- resolved_beta_flags %>%
  left_join(other_flags, by = "sample_id") %>%
  mutate(across(where(is.logical), ~replace_na(.x, FALSE)))

# =========================================================
# 9) Final phenotype heuristics from resolved calls
# =========================================================
drug_key <- tribble(
  ~abx, ~drug,
  "F",   "nitrofurantoin",
  "MEL", "mecillinam",
  "CFR", "cefadroxil",
  "W",   "trimethoprim",
  "SXT", "trimethoprim/sulfamethoxazole",
  "MEM", "meropenem",
  "AMP", "ampicillin",
  "AMC", "amoxicillin/clavulanic acid",
  "PIP", "piperacillin",
  "TZP", "piperacillin/tazobactam",
  "CAZ", "ceftazidime",
  "CRO", "ceftriaxone",
  "CTX", "cefotaxime",
  "FEP", "cefepime",
  "CIP", "ciprofloxacin",
  "OFX", "ofloxacin",
  "LVX", "levofloxacin",
  "MFX", "moxifloxacin",
  "GEN", "gentamicin",
  "TOB", "tobramycin"
)

resolved_predictions_long <- resolved_features %>%
  select(
    sample_id,
    has_carbapenemase, has_esbl, has_ampc, has_oxa1, has_penicillinase,
    has_dfr, has_sul, has_qnr, has_gyr_par, has_gen_tob, has_oqx
  ) %>%
  tidyr::crossing(drug_key) %>%
  mutate(
    has_any_acquired_beta = has_carbapenemase | has_esbl | has_ampc | has_oxa1 | has_penicillinase,
    
    predicted = case_when(
      abx == "MEM" ~ if_else(has_carbapenemase, "R", "S"),
      
      abx %in% c("AMP", "PIP") ~ case_when(
        has_carbapenemase     ~ "R",
        has_any_acquired_beta ~ "R",
        TRUE                  ~ "S?"
      ),
      
      abx %in% c("AMC", "TZP") ~ case_when(
        has_carbapenemase ~ "R",
        has_ampc          ~ "R?",
        has_oxa1          ~ "R?",
        has_esbl          ~ "I/R?",
        has_penicillinase ~ "S?",
        TRUE              ~ "S?"
      ),
      
      abx %in% c("CAZ", "CRO", "CTX") ~ case_when(
        has_carbapenemase ~ "R",
        has_esbl | has_ampc ~ "R",
        TRUE ~ "S?"
      ),
      
      abx == "FEP" ~ case_when(
        has_carbapenemase ~ "R",
        has_esbl | has_ampc ~ "I/R?",
        TRUE ~ "S?"
      ),
      
      abx == "CFR" ~ case_when(
        has_carbapenemase ~ "R",
        has_esbl | has_ampc ~ "R",
        TRUE ~ "S?"
      ),
      
      abx == "MEL" ~ "?",
      abx == "F"   ~ if_else(has_oqx, "I/R?", "?"),
      
      abx == "W"   ~ if_else(has_dfr, "R", "S?"),
      abx == "SXT" ~ if_else(has_dfr | has_sul, "R", "S?"),
      
      abx %in% c("CIP", "OFX", "LVX", "MFX") ~ case_when(
        has_gyr_par ~ "R?",
        has_qnr     ~ "I/R?",
        TRUE        ~ "S?"
      ),
      
      abx %in% c("GEN", "TOB") ~ if_else(has_gen_tob, "R?", "S?"),
      TRUE ~ "?"
    )
  ) %>%
  select(sample_id, abx, drug, predicted)

resolved_predictions_wide <- resolved_predictions_long %>%
  select(sample_id, abx, predicted) %>%
  pivot_wider(names_from = abx, values_from = predicted) %>%
  arrange(sample_id)

resolved_predictions_summary <- resolved_predictions_long %>%
  count(abx, drug, predicted) %>%
  group_by(abx, drug) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(abx, desc(n))

# =========================================================
# 10) Write outputs
# =========================================================
readr::write_csv(comparison, file.path(outdir, "amrfinder_ariba_comparison.csv"))
readr::write_csv(partials_resolved, file.path(outdir, "amrfinder_partials_resolved_by_ariba.csv"))
readr::write_csv(resolved_calls, file.path(outdir, "resolved_amr_gene_calls_annotated.csv"))

readr::write_csv(resolved_beta_hits, file.path(outdir, "resolved_beta_hits_long.csv"))
readr::write_csv(resolved_beta_flags, file.path(outdir, "resolved_beta_flags_per_sample.csv"))

readr::write_csv(resolved_predictions_long, file.path(outdir, "resolved_predicted_SIR_long.csv"))
readr::write_csv(resolved_predictions_wide, file.path(outdir, "resolved_predicted_SIR_wide.csv"))
readr::write_csv(resolved_predictions_summary, file.path(outdir, "resolved_predicted_SIR_summary.csv"))
