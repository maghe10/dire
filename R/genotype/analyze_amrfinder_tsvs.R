source(file='common.R')

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

# =========================
# 0) Inställningar
# =========================
folder <- amrfinderDirectory
files <- list.files(folder, pattern = "\\.tsv$", full.names = TRUE)

files <- files[!grepl(paste0("sample14.tsv", "$"), files, perl = TRUE) ]
files <- files[!grepl(paste0("sample38.tsv", "$"), files, perl = TRUE) ]

stopifnot(length(files) > 0)

outdir <- paste(processedRootRassembly, "genotype", sep="/")


# Antibiotikalista (din)
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

# =========================
# 1) Läs in AMRFinder-filer
# =========================
files <- files[
  !normalize_sample_id(basename(files)) %in% c("014", "038")
]

read_amrfinder_tsv <- function(f) {
  readr::read_tsv(f, col_types = cols(.default = "c"), progress = FALSE) %>%
    janitor::clean_names() %>%
    mutate(
      source_file = basename(f),
      sample_id = normalize_sample_id(source_file)
    ) %>%
    filter(!is.na(sample_id))
}

all_hits <- purrr::map_dfr(files, read_amrfinder_tsv) %>%
  relocate(sample_id, source_file)


# Behåll AMR-rader (AMRFinder har ofta virulens/stress också)
amr <- all_hits %>%
  filter(type == "AMR") %>%
  mutate(
    element_symbol = replace_na(element_symbol, ""),
    element_name   = replace_na(element_name, ""),
    class          = replace_na(class, ""),
    subclass       = replace_na(subclass, ""),
    method         = replace_na(method, ""),
    # gör coverage/identity numeriska om de finns
    percent_coverage_of_reference = suppressWarnings(as.numeric(percent_coverage_of_reference)),
    percent_identity_to_reference = suppressWarnings(as.numeric(percent_identity_to_reference))
  ) %>%
  relocate(sample_id, source_file)

# =========================
# 2) Plocka ut β-laktamrelevanta träffar + kategorisera
# =========================
beta_hits <- amr %>%
  filter(class == "BETA-LACTAM" | str_detect(element_symbol, "^bla")) %>%
  mutate(
    element_name_lc = str_to_lower(element_name),
    
    # Flagga partiella träffar (contig-ände etc.)
    hit_is_partial = str_detect(method, "PARTIAL") |
      (!is.na(percent_coverage_of_reference) & percent_coverage_of_reference < 90),
    
    # Kategorisera β-laktamaser
    beta_group = case_when(
      # Intrinsik (vanlig i E. coli; ofta "bakgrund")
      element_symbol == "blaEC" ~ "Intrinsik blaEC (E. coli klass C)",
      
      # Karbapenemaser (utöka gärna listan om du vill)
      str_detect(element_symbol, "^bla(KPC|NDM|VIM|IMP|OXA-48|OXA-181|OXA-232|GES|SME|IMI|NMC)") ~ "Carbapenemase",
      str_detect(element_name_lc, "carbapenemase") ~ "Carbapenemase",
      
      # ESBL
      str_detect(element_symbol, "^blaCTX") ~ "ESBL (CTX-M)",
      str_detect(element_name_lc, "extended-spectrum class a beta-lactamase|extended-spectrum beta-lactamase") ~ "ESBL",
      
      # AmpC (förvärvad/plasmid)
      str_detect(element_symbol, "^bla(CMY|DHA|FOX|ACC|ACT|MIR|LAT|MOX)") ~ "AmpC (förvärvad/plasmid)",
      str_detect(element_name_lc, "extended-spectrum class c beta-lactamase") ~ "AmpC (förvärvad/plasmid)",
      
      # OXA-1 (kan påverka AMC/TZP mer än “vanlig” TEM-1)
      element_symbol == "blaOXA-1" ~ "OXA-1",
      
      # “Penicillinas”
      element_symbol %in% c("blaTEM-1", "blaSHV-1") ~ "Penicillinas (TEM-1/SHV-1)",
      str_detect(element_name_lc, "broad-spectrum class a beta-lactamase tem-1") ~ "Penicillinas (TEM-1/SHV-1)",
      
      # Övriga familjer
      str_detect(element_symbol, "^blaTEM") ~ "TEM-familj (övrig)",
      str_detect(element_symbol, "^blaSHV") ~ "SHV-familj (övrig)",
      str_detect(element_symbol, "^blaOXA") ~ "OXA-familj (övrig)",
      str_detect(element_symbol, "^bla") ~ "Övrig β-laktamas",
      TRUE ~ NA_character_
    )
  )

# =========================
# 3) Bygg provvisa “flags” (ESBL/AmpC/Carba osv.)
#    OBS: Vi exkluderar blaEC när vi avgör "förvärvad β-laktamas"
# =========================
beta_flags <- beta_hits %>%
  group_by(sample_id) %>%
  summarise(
    beta_genes_all = paste(sort(unique(element_symbol)), collapse = "; "),
    beta_groups_all = paste(sort(unique(na.omit(beta_group))), collapse = "; "),
    any_partial_beta_hit = any(hit_is_partial, na.rm = TRUE),
    
    # Förvärvade (allt utom blaEC)
    has_carbapenemase = any(beta_group == "Carbapenemase" & element_symbol != "blaEC"),
    has_esbl          = any(str_detect(beta_group, "^ESBL") & element_symbol != "blaEC"),
    has_ampc          = any(beta_group == "AmpC (förvärvad/plasmid)" & element_symbol != "blaEC"),
    has_oxa1          = any(beta_group == "OXA-1" & element_symbol != "blaEC"),
    has_penicillinase = any(beta_group == "Penicillinas (TEM-1/SHV-1)" & element_symbol != "blaEC"),
    
    has_intrinsic_blaEC = any(element_symbol == "blaEC"),
    
    primary_beta_category = case_when(
      has_carbapenemase ~ "Carbapenemase",
      has_esbl          ~ "ESBL",
      has_ampc          ~ "AmpC (förvärvad/plasmid)",
      has_oxa1          ~ "OXA-1",
      has_penicillinase ~ "Penicillinas (TEM-1/SHV-1)",
      has_intrinsic_blaEC ~ "Endast intrinsik blaEC",
      TRUE              ~ "Inga β-laktamaser hittade"
    ),
    
    .groups = "drop"
  )

# =========================
# 4) Summeringar över alla prov
# =========================
beta_category_counts <- beta_flags %>%
  count(primary_beta_category, sort = TRUE)

beta_gene_frequency <- beta_hits %>%
  distinct(sample_id, element_symbol, element_name, beta_group) %>%
  count(beta_group, element_symbol, element_name, sort = TRUE, name = "n_samples")

# =========================
# 5) (Valfritt men ofta nyttigt) Features för förmodad S/R
# =========================
other_flags <- amr %>%
  group_by(sample_id) %>%
  summarise(
    has_dfr = any(str_detect(element_symbol, "^dfr")),
    has_sul = any(element_symbol %in% c("sul1", "sul2", "sul3")),
    
    has_qnr = any(str_detect(element_symbol, "^qnr")),
    has_gyr_par = any(str_detect(element_symbol, "^(gyrA|parC|parE)_")),  # t.ex. gyrA_S83L
    
    # GEN/TOB-relevant (grov heuristik)
    has_gen_tob = any(
      str_detect(element_symbol, "aac\\(3") |
        str_detect(element_symbol, "aac\\(6") |
        str_detect(element_symbol, "aph\\(2") |
        str_detect(element_symbol, "^rmt") |
        str_detect(element_symbol, "^armA$")
    ),
    
    # Nitrofurantoin (svag heuristik)
    has_oqx = any(str_detect(element_symbol, "^oqx(A|B)$")),
    
    .groups = "drop"
  )

features <- beta_flags %>%
  left_join(other_flags, by = "sample_id") %>%
  mutate(across(everything(), ~replace_na(.x, FALSE)))

# =========================
# 6) Regelbaserad “förmodad S/R” (genotyp -> grov fenotyp)
#    OBS: heuristik! Anpassa gärna efter era lokala tolkningar.
# =========================
predictions_long <- features %>%
  select(sample_id,
         has_carbapenemase, has_esbl, has_ampc, has_oxa1, has_penicillinase,
         has_dfr, has_sul, has_qnr, has_gyr_par, has_gen_tob, has_oqx) %>%
  tidyr::crossing(drug_key) %>%
  mutate(
    has_any_acquired_beta = has_carbapenemase | has_esbl | has_ampc | has_oxa1 | has_penicillinase,
    
    predicted = case_when(
      # Karbapenem
      abx == "MEM" ~ if_else(has_carbapenemase, "R", "S"),
      
      # Penicilliner
      abx %in% c("AMP","PIP") ~ case_when(
        has_carbapenemase    ~ "R",
        has_any_acquired_beta ~ "R",
        TRUE                 ~ "S?"
      ),
      
      # Hämmarkombinationer
      abx %in% c("AMC","TZP") ~ case_when(
        has_carbapenemase ~ "R",
        has_ampc          ~ "R?",
        has_oxa1          ~ "R?",
        has_esbl          ~ "I/R?",
        has_penicillinase ~ "S?",
        TRUE              ~ "S?"
      ),
      
      # Cefalosporiner (3:e gen)
      abx %in% c("CAZ","CRO","CTX") ~ case_when(
        has_carbapenemase ~ "R",
        has_esbl | has_ampc ~ "R",
        TRUE ~ "S?"
      ),
      
      # Cefepim
      abx == "FEP" ~ case_when(
        has_carbapenemase ~ "R",
        has_esbl | has_ampc ~ "I/R?",
        TRUE ~ "S?"
      ),
      
      # Cefadroxil
      abx == "CFR" ~ case_when(
        has_carbapenemase ~ "R",
        has_esbl | has_ampc ~ "R",
        TRUE ~ "S?"
      ),
      
      abx == "MEL" ~ "?",
      abx == "F"   ~ if_else(has_oqx, "I/R?", "?"),
      
      abx == "W"   ~ if_else(has_dfr, "R", "S?"),
      abx == "SXT" ~ if_else(has_dfr | has_sul, "R", "S?"),
      
      abx %in% c("CIP","OFX","LVX","MFX") ~ case_when(
        has_gyr_par ~ "R?",
        has_qnr ~ "I/R?",
        TRUE ~ "S?"
      ),
      
      abx %in% c("GEN","TOB") ~ if_else(has_gen_tob, "R?", "S?"),
      
      TRUE ~ "?"
    )
  ) %>%
  select(sample_id, abx, drug, predicted)

predictions_wide <- predictions_long %>%
  select(sample_id, abx, predicted) %>%
  tidyr::pivot_wider(names_from = abx, values_from = predicted) %>%
  arrange(sample_id)

predictions_summary <- predictions_long %>%
  count(abx, drug, predicted) %>%
  group_by(abx, drug) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(abx, desc(n))

# =========================
# 7) Spara utdata
# =========================
write_csv2(all_hits, file.path(outdir, "amrfinder_all_hits_long.csv"))
write_csv2(amr, file.path(outdir, "amrfinder_amr_long.csv"))

write_csv2(beta_hits, file.path(outdir, "amrfinder_beta_hits_long.csv"))
write_csv2(beta_flags, file.path(outdir, "amrfinder_beta_flags_per_sample.csv"))
write_csv2(beta_category_counts, file.path(outdir, "amrfinder_beta_primary_category_counts.csv"))
write_csv2(beta_gene_frequency, file.path(outdir, "amrfinder_beta_gene_frequency_n_samples.csv"))

write_csv2(predictions_long, file.path(outdir, "amrfinder_predicted_SIR_long.csv"))
write_csv2(predictions_wide, file.path(outdir, "amrfinder_predicted_SIR_wide.csv"))
write_csv2(predictions_summary, file.path(outdir, "amrfinder_predicted_SIR_summary.csv"))