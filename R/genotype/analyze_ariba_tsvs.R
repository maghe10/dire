source(file='common.R')

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

# =========================
# 0) Inställningar
# =========================
folder <- aribaDirectory
outdir <- paste(processedRootRassembly, "genotype", sep="/")

files <- list.files(folder, pattern = "^sample.*\\.tsv$", full.names = TRUE)
stopifnot(length(files) > 0)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 1) Läs in ARIBA report.tsv-filer
# =========================
read_ariba_tsv <- function(f) {
  readr::read_tsv(f, col_types = cols(.default = "c"), progress = FALSE) %>%
    janitor::clean_names() %>%
    mutate(
      source_file = basename(f),
      sample_id = normalize_sample_id(source_file)
    ) %>%
    filter(!is.na(sample_id))
}

all_hits <- purrr::map_dfr(files, read_ariba_tsv) %>%
  relocate(sample_id, source_file) %>%
  filter(!sample_id %in% c("014", "038"))





# gör numeriska kolumner numeriska där det passar
all_hits <- all_hits %>%
  mutate(
    ref_len = suppressWarnings(as.numeric(ref_len)),
    ref_base_assembled = suppressWarnings(as.numeric(ref_base_assembled)),
    pc_ident = suppressWarnings(as.numeric(pc_ident)),
    ctg_len = suppressWarnings(as.numeric(ctg_len)),
    ctg_cov = suppressWarnings(as.numeric(ctg_cov)),
    reads = suppressWarnings(as.numeric(reads)),
    flag = suppressWarnings(as.numeric(flag)),
    known_var = suppressWarnings(as.numeric(known_var)),
    var_only = suppressWarnings(as.numeric(var_only))
  )

# =========================
# 2) Extrahera mer lättläst gennamn
# =========================
# free_text brukar se ut ungefär:
# "Original name: blaTEM-1B_1_AY458016"
# eller "Original name: aadA1_5_JX185132"
all_hits <- all_hits %>%
  mutate(
    original_name = str_match(free_text, "Original name:\\s*(.+)$")[, 2],
    gene_symbol = case_when(
      !is.na(original_name) ~ original_name %>%
        str_remove("_[0-9]+_[A-Z0-9]+$") %>%   # tar bort suffix som _1_AY458016
        str_remove("_[0-9]+$"),
      TRUE ~ ref_name
    )
  )

# =========================
# 3) Klassificera hit-kvalitet / partial
# =========================
all_hits <- all_hits %>%
  mutate(
    assembly_fraction = if_else(!is.na(ref_len) & ref_len > 0,
                                ref_base_assembled / ref_len,
                                NA_real_),
    
    hit_completeness = case_when(
      !is.na(assembly_fraction) & assembly_fraction >= 0.95 ~ "complete",
      !is.na(assembly_fraction) & assembly_fraction >= 0.80 ~ "near_complete",
      !is.na(assembly_fraction) & assembly_fraction >= 0.50 ~ "partial",
      !is.na(assembly_fraction) ~ "fragment",
      TRUE ~ NA_character_
    ),
    
    hit_confidence = case_when(
      !is.na(pc_ident) & !is.na(assembly_fraction) &
        pc_ident >= 99 & assembly_fraction >= 0.95 ~ "high",
      !is.na(pc_ident) & !is.na(assembly_fraction) &
        pc_ident >= 95 & assembly_fraction >= 0.80 ~ "moderate",
      !is.na(pc_ident) & !is.na(assembly_fraction) &
        pc_ident >= 90 & assembly_fraction >= 0.50 ~ "low",
      TRUE ~ "weak"
    ),
    
    is_partial = case_when(
      !is.na(assembly_fraction) & assembly_fraction < 0.95 ~ TRUE,
      TRUE ~ FALSE
    )
  )

# =========================
# 4) Gen-tabell: en rad per prov + gen/hit
#    Kollapsa ARIBA-variant-rader till en gen-nivå
# =========================
# Vi grupperar på sample + ref_name + cluster + gene_symbol
# eftersom samma gene/hit ofta ligger på flera rader.
ariba_genes <- all_hits %>%
  group_by(sample_id, source_file, cluster, ref_name, gene_symbol) %>%
  summarise(
    ref_len = max(ref_len, na.rm = TRUE),
    ref_base_assembled = max(ref_base_assembled, na.rm = TRUE),
    assembly_fraction = max(assembly_fraction, na.rm = TRUE),
    pc_ident = max(pc_ident, na.rm = TRUE),
    ctg_cov = max(ctg_cov, na.rm = TRUE),
    reads = max(reads, na.rm = TRUE),
    n_rows = n(),
    
    # sammanfatta variantsignaler
    n_variant_rows = sum(!is.na(ref_ctg_change) & ref_ctg_change != ".", na.rm = TRUE),
    ref_ctg_changes = paste(sort(unique(ref_ctg_change[!is.na(ref_ctg_change) & ref_ctg_change != "."])),
                            collapse = "; "),
    ref_ctg_effects = paste(sort(unique(ref_ctg_effect[!is.na(ref_ctg_effect) & ref_ctg_effect != "."])),
                            collapse = "; "),
    has_known_var = any(has_known_var %in% c("1", "TRUE", "True"), na.rm = TRUE),
    
    original_name = paste(sort(unique(original_name[!is.na(original_name)])), collapse = "; "),
    hit_completeness = case_when(
      max(assembly_fraction, na.rm = TRUE) >= 0.95 ~ "complete",
      max(assembly_fraction, na.rm = TRUE) >= 0.80 ~ "near_complete",
      max(assembly_fraction, na.rm = TRUE) >= 0.50 ~ "partial",
      TRUE ~ "fragment"
    ),
    hit_confidence = case_when(
      max(pc_ident, na.rm = TRUE) >= 99 & max(assembly_fraction, na.rm = TRUE) >= 0.95 ~ "high",
      max(pc_ident, na.rm = TRUE) >= 95 & max(assembly_fraction, na.rm = TRUE) >= 0.80 ~ "moderate",
      max(pc_ident, na.rm = TRUE) >= 90 & max(assembly_fraction, na.rm = TRUE) >= 0.50 ~ "low",
      TRUE ~ "weak"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    is_partial = hit_completeness != "complete"
  ) %>%
  relocate(sample_id, source_file, gene_symbol, cluster, ref_name)

# fix för grupper där alla värden var NA
ariba_genes <- ariba_genes %>%
  mutate(
    ref_ctg_changes = na_if(ref_ctg_changes, ""),
    ref_ctg_effects = na_if(ref_ctg_effects, ""),
    original_name = na_if(original_name, "")
  )

# =========================
# 5) β-laktamfilter + kategorisering
# =========================
beta_hits <- ariba_genes %>%
  filter(str_detect(gene_symbol, "^bla") | str_detect(cluster, "^bla")) %>%
  mutate(
    gene_symbol_lc = str_to_lower(gene_symbol),
    
    beta_group = case_when(
      gene_symbol == "blaEC" ~ "Intrinsik blaEC (E. coli klass C)",
      
      str_detect(gene_symbol, "^bla(KPC|NDM|VIM|IMP|OXA-48|OXA-181|OXA-232|GES|SME|IMI|NMC)") ~ "Carbapenemase",
      
      str_detect(gene_symbol, "^blaCTX-M") ~ "ESBL (CTX-M)",
      
      str_detect(gene_symbol, "^bla(CMY|DHA|FOX|ACC|ACT|MIR|LAT|MOX)") ~ "AmpC (förvärvad/plasmid)",
      
      gene_symbol == "blaOXA-1" ~ "OXA-1",
      
      gene_symbol %in% c("blaTEM-1", "blaTEM-1A", "blaTEM-1B", "blaTEM-1C", "blaSHV-1") ~
        "Penicillinas (TEM-1/SHV-1)",
      
      str_detect(gene_symbol, "^blaTEM") ~ "TEM-familj (övrig)",
      str_detect(gene_symbol, "^blaSHV") ~ "SHV-familj (övrig)",
      str_detect(gene_symbol, "^blaOXA") ~ "OXA-familj (övrig)",
      TRUE ~ "Övrig β-laktamas"
    )
  ) %>%
  relocate(sample_id, gene_symbol, beta_group, everything())

# =========================
# 6) Summeringar
# =========================
# per sample: alla AMR-gener
amr_by_sample <- ariba_genes %>%
  group_by(sample_id) %>%
  summarise(
    n_genes = n(),
    genes = paste(sort(unique(gene_symbol)), collapse = "; "),
    n_partial = sum(is_partial, na.rm = TRUE),
    .groups = "drop"
  )

# genfrekvens över alla prov
gene_frequency <- ariba_genes %>%
  distinct(sample_id, gene_symbol, cluster, hit_completeness, hit_confidence) %>%
  count(gene_symbol, cluster, sort = TRUE, name = "n_samples")

# partial hits
partial_hits <- ariba_genes %>%
  filter(is_partial) %>%
  arrange(sample_id, desc(assembly_fraction), gene_symbol)

# beta-laktam per sample
beta_flags <- beta_hits %>%
  group_by(sample_id) %>%
  summarise(
    beta_genes_all = paste(sort(unique(gene_symbol)), collapse = "; "),
    beta_groups_all = paste(sort(unique(beta_group)), collapse = "; "),
    any_partial_beta_hit = any(is_partial, na.rm = TRUE),
    
    has_carbapenemase = any(beta_group == "Carbapenemase" & gene_symbol != "blaEC"),
    has_esbl          = any(beta_group == "ESBL (CTX-M)" & gene_symbol != "blaEC"),
    has_ampc          = any(beta_group == "AmpC (förvärvad/plasmid)" & gene_symbol != "blaEC"),
    has_oxa1          = any(beta_group == "OXA-1" & gene_symbol != "blaEC"),
    has_penicillinase = any(beta_group == "Penicillinas (TEM-1/SHV-1)" & gene_symbol != "blaEC"),
    has_intrinsic_blaEC = any(gene_symbol == "blaEC"),
    
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

beta_gene_frequency <- beta_hits %>%
  distinct(sample_id, gene_symbol, beta_group, hit_completeness, hit_confidence) %>%
  count(beta_group, gene_symbol, sort = TRUE, name = "n_samples")

# presence/absence-matris
presence_matrix <- ariba_genes %>%
  distinct(sample_id, gene_symbol) %>%
  mutate(present = 1L) %>%
  pivot_wider(
    names_from = gene_symbol,
    values_from = present,
    values_fill = 0L
  ) %>%
  arrange(sample_id)

# =========================
# 7) Spara utdata
# =========================
readr::write_csv2(all_hits,        file.path(outdir, "ariba_all_hits_long.csv"))
readr::write_csv2(ariba_genes,     file.path(outdir, "ariba_genes_long.csv"))
readr::write_csv2(amr_by_sample,   file.path(outdir, "ariba_genes_per_sample.csv"))
readr::write_csv2(gene_frequency,  file.path(outdir, "ariba_gene_frequency.csv"))
readr::write_csv2(partial_hits,    file.path(outdir, "ariba_partial_hits.csv"))

readr::write_csv2(beta_hits,           file.path(outdir, "ariba_beta_hits_long.csv"))
readr::write_csv2(beta_flags,          file.path(outdir, "ariba_beta_flags_per_sample.csv"))
readr::write_csv2(beta_gene_frequency, file.path(outdir, "ariba_beta_gene_frequency.csv"))
readr::write_csv2(presence_matrix,     file.path(outdir, "ariba_presence_absence_matrix.csv"))

