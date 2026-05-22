source(file = "common.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(janitor)
})

outdir <- file.path(processedRootRassembly, "genotype")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

amrf_file  <- file.path(outdir, "amrfinder_amr_long.csv")
ariba_file <- file.path(outdir, "ariba_genes_long.csv")
bldb_file  <- file.path(outdir, "bldb_lookup.csv")
bldb_alias_file <- file.path(outdir, "bldb_aliases.csv")

stopifnot(file.exists(amrf_file))
stopifnot(file.exists(ariba_file))
stopifnot(file.exists(bldb_file))
stopifnot(file.exists(bldb_alias_file))

amrf <- read.csv2(
  amrf_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  janitor::clean_names()

ariba <- read.csv2(
  ariba_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  janitor::clean_names()

bldb <- read.csv2(
  bldb_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  janitor::clean_names()

bldb_aliases <- if (file.exists(bldb_alias_file)) {
  read.csv2(
    bldb_alias_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) |>
    janitor::clean_names()
} else {
  tibble()
}

normalize_gene <- function(x) {
  
  x |>
    replace_na("") |>
    str_replace("^bla","") |>
    str_replace_all("\\s+","") |>
    str_to_upper()
  
}

get_gene_family <- function(x){
  
  case_when(
    str_detect(x,"^blaCTX-M") ~ "blaCTX-M",
    str_detect(x,"^blaTEM") ~ "blaTEM",
    str_detect(x,"^blaSHV") ~ "blaSHV",
    str_detect(x,"^blaOXA") ~ "blaOXA",
    str_detect(x,"^blaKPC") ~ "blaKPC",
    str_detect(x,"^blaNDM") ~ "blaNDM",
    str_detect(x,"^blaVIM") ~ "blaVIM",
    str_detect(x,"^blaIMP") ~ "blaIMP",
    str_detect(x,"^blaCMY") ~ "blaCMY",
    str_detect(x,"^blaDHA") ~ "blaDHA",
    str_detect(x,"^blaFOX") ~ "blaFOX",
    str_detect(x,"^blaACC") ~ "blaACC",
    str_detect(x,"^blaACT") ~ "blaACT",
    str_detect(x,"^blaMIR") ~ "blaMIR",
    str_detect(x,"^blaLAT") ~ "blaLAT",
    str_detect(x,"^blaMOX") ~ "blaMOX",
    str_detect(x,"^blaGES") ~ "blaGES",
    str_detect(x,"^blaVEB") ~ "blaVEB",
    str_detect(x,"^blaPER") ~ "blaPER",
    str_detect(x,"^blaEC") ~ "blaEC",
    TRUE ~ x
  )
  
}

get_gene_family_plain <- function(x){
  
  x2 <- normalize_gene(x)
  
  case_when(
    str_detect(x2,"^CTX-M") ~ "CTX-M",
    str_detect(x2,"^TEM") ~ "TEM",
    str_detect(x2,"^SHV") ~ "SHV",
    str_detect(x2,"^OXA") ~ "OXA",
    str_detect(x2,"^KPC") ~ "KPC",
    str_detect(x2,"^NDM") ~ "NDM",
    str_detect(x2,"^VIM") ~ "VIM",
    str_detect(x2,"^IMP") ~ "IMP",
    str_detect(x2,"^CMY") ~ "CMY",
    str_detect(x2,"^DHA") ~ "DHA",
    str_detect(x2,"^FOX") ~ "FOX",
    str_detect(x2,"^ACC") ~ "ACC",
    str_detect(x2,"^ACT") ~ "ACT",
    str_detect(x2,"^MIR") ~ "MIR",
    str_detect(x2,"^LAT") ~ "LAT",
    str_detect(x2,"^MOX") ~ "MOX",
    str_detect(x2,"^GES") ~ "GES",
    str_detect(x2,"^VEB") ~ "VEB",
    str_detect(x2,"^PER") ~ "PER",
    str_detect(x2,"^EC") ~ "EC",
    TRUE ~ x2
  )
  
}

infer_beta_group <- function(gene){
  
  case_when(
    gene=="blaEC" ~ "Intrinsik blaEC (E. coli klass C)",
    str_detect(gene,"^bla(KPC|NDM|VIM|IMP|OXA-48|OXA-181|OXA-232|GES|SME|IMI|NMC)") ~ "Carbapenemase",
    str_detect(gene,"^blaCTX-M") ~ "ESBL (CTX-M)",
    str_detect(gene,"^bla(CMY|DHA|FOX|ACC|ACT|MIR|LAT|MOX)") ~ "AmpC (förvärvad/plasmid)",
    gene=="blaOXA-1" ~ "OXA-1",
    gene %in% c("blaTEM-1","blaTEM-1A","blaTEM-1B","blaTEM-1C","blaSHV-1") ~ "Penicillinas (TEM-1/SHV-1)",
    str_detect(gene,"^blaTEM") ~ "TEM-familj (övrig)",
    str_detect(gene,"^blaSHV") ~ "SHV-familj (övrig)",
    str_detect(gene,"^blaOXA") ~ "OXA-familj (övrig)",
    str_detect(gene,"^bla") ~ "Övrig β-laktamas",
    TRUE ~ NA_character_
  )
  
}

amrf_best <- amrf |>
  filter(type=="AMR") |>
  mutate(
    amrf_gene = element_symbol,
    amrf_family = get_gene_family(amrf_gene),
    amrf_cov = suppressWarnings(as.numeric(percent_coverage_of_reference)),
    amrf_ident = suppressWarnings(as.numeric(percent_identity_to_reference)),
    amrf_method = replace_na(method,""),
    amrf_partial = str_detect(amrf_method,"PARTIAL") |
      (!is.na(amrf_cov) & amrf_cov < 90)
  ) |>
  group_by(sample_id,amrf_family) |>
  slice_max(order_by=amrf_cov,with_ties=FALSE) |>
  ungroup()


if (!"gene_symbol" %in% names(ariba)) ariba[["gene_symbol"]] <- NA_character_
if (!"original_name" %in% names(ariba)) ariba[["original_name"]] <- NA_character_
if (!"ref_name" %in% names(ariba)) ariba[["ref_name"]] <- NA_character_
if (!"gene" %in% names(ariba)) ariba[["gene"]] <- NA_character_
if (!"hit_completeness" %in% names(ariba)) ariba[["hit_completeness"]] <- NA_character_

ariba_best <- ariba |>
  mutate(
    ariba_gene = case_when(
      !is.na(gene_symbol) & gene_symbol != "" ~ gene_symbol,
      !is.na(original_name) & original_name != "" ~ original_name,
      !is.na(ref_name) & ref_name != "" ~ ref_name,
      !is.na(gene) & gene != "" ~ as.character(gene),
      TRUE ~ NA_character_
    ),
    ariba_family = get_gene_family(ariba_gene),
    assembly_fraction = suppressWarnings(as.numeric(assembly_fraction)),
    pc_ident = suppressWarnings(as.numeric(pc_ident)),
    hit_completeness = tidyr::replace_na(hit_completeness, "fragment"),
    ariba_complete = hit_completeness %in% c("complete", "near_complete")
  ) |>
  filter(!is.na(ariba_gene)) |>
  group_by(sample_id, ariba_family) |>
  slice_max(order_by = assembly_fraction, with_ties = FALSE) |>
  ungroup()



comparison <- full_join(
  amrf_best,
  ariba_best,
  by = c("sample_id", "amrf_family" = "ariba_family")
)  |>
  rename(gene_family = amrf_family)



comparison <- comparison |>
  mutate(
    
    comparison_status = case_when(
      !is.na(amrf_gene) & !is.na(ariba_gene) & amrf_gene==ariba_gene ~ "exact_match",
      !is.na(amrf_gene) & !is.na(ariba_gene) ~ "family_match",
      !is.na(amrf_gene) & is.na(ariba_gene) ~ "amrfinder_only",
      is.na(amrf_gene) & !is.na(ariba_gene) ~ "ariba_only",
      TRUE ~ "discordant"
    ),
    
    resolved_gene = case_when(
      !is.na(amrf_gene) & !amrf_partial ~ amrf_gene,
      amrf_partial & !is.na(ariba_gene) & ariba_complete ~ ariba_gene,
      !is.na(amrf_gene) ~ amrf_gene,
      is.na(amrf_gene) & !is.na(ariba_gene) ~ ariba_gene,
      TRUE ~ NA_character_
    )
    
  )

bldb_main <- bldb |>
  mutate(
    gene_norm = normalize_gene(gene_symbol_guess)
  )

resolved_calls <- comparison |>
  filter(!is.na(resolved_gene)) |>
  mutate(
    resolved_gene_norm = normalize_gene(resolved_gene)
  ) |>
  left_join(
    bldb_main,
    by=c("resolved_gene_norm"="gene_norm")
  )

resolved_calls <- resolved_calls |>
  mutate(
    resolved_beta_group = infer_beta_group(resolved_gene)
  )

write.csv2(
  comparison,
  file.path(outdir,"amrfinder_ariba_comparison.csv"),
  row.names=FALSE
)

write.csv2(
  resolved_calls,
  file.path(outdir,"resolved_amr_gene_calls_annotated.csv"),
  row.names=FALSE
)