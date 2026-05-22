# localBLBD.R (fixad, vektoriserad)
suppressPackageStartupMessages({
  library(rvest)
  library(xml2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(janitor)
  library(tidyr)
})

# ---- konfigurera ----
# se till att common.R finns och sätter processedRootRassembly
if (file.exists("common.R")) source("common.R") else {
  stop("Missing common.R — skapa eller kopiera den innan du kör.")
}
outdir <- file.path(processedRootRassembly, "genotype")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pages <- tibble::tribble(
  ~prot, ~bldb_group,   ~url,
  "A",   "Class A",     "http://bldb.eu/BLDB.php?prot=A#AA1-1",
  "B1",  "Subclass B1", "http://bldb.eu/BLDB.php?prot=B1#AFM-1",
  "B2",  "Subclass B2", "http://bldb.eu/BLDB.php?prot=B2#CphA-1",
  "B3",  "Subclass B3", "http://bldb.eu/BLDB.php?prot=B3#AHM-1",
  "C",   "Class C",     "http://bldb.eu/BLDB.php?prot=C#ACC-1",
  "D",   "Class D",     "http://bldb.eu/BLDB.php?prot=D#AFD-1"
)

# ---- vektoriserade helperfunktioner ----

safe_read_html <- function(url) {
  tryCatch(read_html(url), error = function(e) {
    warning("Failed to read ", url, " : ", e$message)
    return(NULL)
  })
}

clean_table_names <- function(tab) {
  tab <- tab %>% janitor::clean_names()
  tab <- tab %>% rename_with(~str_replace(.x, "^number_of_pd_bstructures$", "number_of_pdbstructures"))
  # tvinga till character och trimma
  tab <- tab %>% mutate(across(everything(), ~ as.character(.x)))
  tab <- tab %>% mutate(across(where(is.character), ~ str_squish(.x)))
  tab
}

# Är en hel rad en 'header/section' (t.ex. "Sequence alignment for ACC")?
is_header_row_vec <- function(amblerclass, proteinname, phenotype, natural_col) {
  # vektoriserad logik: returnerar logisk vektor
  amblerclass <- replace_na(amblerclass, "")
  proteinname <- replace_na(proteinname, "")
  phenotype <- replace_na(phenotype, "")
  natural_col <- replace_na(natural_col, "")
  
  cond1 <- str_detect(phenotype, regex("^Sequence alignment for\\b", ignore_case = TRUE))
  cond2 <- str_detect(natural_col, regex("^Sequence alignment for\\b", ignore_case = TRUE))
  cond3 <- (proteinname == amblerclass) & !(amblerclass %in% c("A","B1","B2","B3","C","D",""))
  cond1 | cond2 | cond3
}

# Bush-Jacoby extraktion (vektoriserad via map_chr där nödigt)
extract_bush_jacoby_one <- function(x) {
  x <- replace_na(x, "") %>% str_squish()
  if (x == "" || x %in% c("—", "-")) return(NA_character_)
  x2 <- str_replace_all(x, "\\?", "") %>%
    str_replace_all("\\band\\b", ";") %>%
    str_replace_all(",", ";") %>%
    str_replace_all("\\+", " + ") %>%
    str_squish()
  hits <- str_extract_all(x2, "\\b[123](?:[a-z]{1,3})?\\b")[[1]]
  hits <- unique(hits)
  if (length(hits) == 0) return(NA_character_)
  paste(hits, collapse = ";")
}
is_bush_jacoby_uncertain_one <- function(x) {
  str_detect(replace_na(x, ""), "\\?")
}

# Vektoriserad to_bla_name (säker för vektorer)
to_bla_name_vec <- function(x) {
  x0 <- replace_na(x, "")
  x0 <- str_squish(x0)
  # tomt -> NA
  res <- ifelse(x0 == "", NA_character_, x0)
  # om redan börjar med 'bla' behåll, annars prefixa 'bla'
  has_bla <- ifelse(is.na(res), FALSE, str_detect(res, regex("^bla", ignore_case = TRUE)))
  res2 <- ifelse(has_bla, res, ifelse(is.na(res), NA_character_, paste0("bla", str_replace_all(res, "\\s+|\\(|\\)", ""))))
  # rena tomma strängar till NA
  res2[res2 == ""] <- NA_character_
  res2
}

normalize_name_vec <- function(x) {
  replace_na(x, "") %>% str_squish()
}

# ---- läs och rengör tabeller från varje sida ----
read_bldb_table <- function(url, prot, bldb_group) {
  message("Reading: ", url)
  page <- safe_read_html(url)
  if (is.null(page)) return(tibble())
  
  tables <- html_table(page, fill = TRUE)
  if (length(tables) < 2) {
    warning("No enzyme table found at ", url)
    return(tibble())
  }
  
  tab <- tables[[2]] %>% clean_table_names()
  
  # säkerställ att nödvändiga kolumnnamn finns (lägg till NA om saknas)
  required_cols <- c("amblerclass","proteinname","alternativeprotein_names",
                     "subfamily","gen_pept_id","gen_bank_id","pub_med_id_doi",
                     "sequence","number_of_pdbstructures","mutants",
                     "phenotype","functionalinformation","natural_n_or_acquired_a")
  missing <- setdiff(required_cols, names(tab))
  if (length(missing) > 0) tab[missing] <- NA_character_
  
  tab <- tab %>%
    mutate(
      prot = prot,
      bldb_group = bldb_group,
      source_url = url,
      molecular_class = case_when(
        prot == "A" ~ "A",
        prot %in% c("B1","B2","B3") ~ "B",
        prot == "C" ~ "C",
        prot == "D" ~ "D",
        TRUE ~ NA_character_
      ),
      metallo_subclass = case_when(
        prot == "B1" ~ "B1",
        prot == "B2" ~ "B2",
        prot == "B3" ~ "B3",
        TRUE ~ NA_character_
      )
    )
  
  # identifiera header-rows med vektoriserad funktion
  tab <- tab %>%
    mutate(.is_header = is_header_row_vec(amblerclass, proteinname, phenotype, natural_n_or_acquired_a)) %>%
    filter(!.is_header) %>%
    select(-.is_header)
  
  # drop tomma proteinname-rader
  tab <- tab %>% filter(!is.na(proteinname) & proteinname != "")
  
  # slutstädning
  tab %>% mutate(across(where(is.character), ~ str_squish(.x)))
}

# ---- bygg bldb_raw genom att slå ihop sidor ----
bldb_raw <- purrr::pmap_dfr(
  pages,
  function(prot, bldb_group, url) {
    read_bldb_table(url, prot, bldb_group)
  }
)

# säkerhetsfilter: inga tomma proteinname
bldb_raw <- bldb_raw %>% filter(!is.na(proteinname) & proteinname != "")

# ---- joinable: standardisera och skapa gene_symbol_guess ----
bldb_joinable <- bldb_raw %>%
  mutate(
    proteinname = normalize_name_vec(proteinname),
    gene_symbol_guess = to_bla_name_vec(proteinname),
    gene_symbol_guess = str_replace_all(gene_symbol_guess, "\\s+", ""),
    enzyme_name = proteinname
  ) %>%
  distinct()

bldb_aliases <- bldb_joinable %>%
  mutate(
    bush_jacoby_raw = phenotype,
    bush_jacoby_class = map_chr(phenotype, ~ extract_bush_jacoby_one(.x)),
    bush_jacoby_uncertain = map_lgl(phenotype, ~ is_bush_jacoby_uncertain_one(.x))
  ) %>%
  select(
    proteinname,
    alternativeprotein_names,
    bldb_group,
    molecular_class,
    metallo_subclass,
    subfamily,
    phenotype,
    functionalinformation,
    natural_n_or_acquired_a,
    bush_jacoby_raw,
    bush_jacoby_class,
    bush_jacoby_uncertain,
    source_url,
    gene_symbol_guess
  ) %>%
  mutate(alternativeprotein_names = replace_na(alternativeprotein_names, "")) %>%
  separate_rows(alternativeprotein_names, sep = "\\s*;\\s*|\\s*,\\s*") %>%
  mutate(
    alternativeprotein_names = na_if(str_squish(alternativeprotein_names), ""),
    alias_bla = if_else(
      is.na(alternativeprotein_names),
      NA_character_,
      if_else(
        str_detect(alternativeprotein_names, regex("^bla", ignore_case = TRUE)),
        alternativeprotein_names,
        paste0("bla", str_replace_all(alternativeprotein_names, "\\s+", ""))
      )
    )
  ) %>%
  filter(!is.na(alternativeprotein_names)) %>%
  distinct()

# ---- lookup: kompakt tabell + bush-jacoby ----
bldb_lookup <- bldb_joinable %>%
  mutate(
    bush_jacoby_raw = phenotype,
    bush_jacoby_class = map_chr(phenotype, ~ extract_bush_jacoby_one(.x)),
    bush_jacoby_uncertain = map_lgl(phenotype, ~ is_bush_jacoby_uncertain_one(.x))
  ) %>%
  select(
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
    bush_jacoby_uncertain,
    source_url
  ) %>%
  distinct()

# ---- skriv ut ----
write_tsv(bldb_raw, file.path(outdir, "bldb_raw.tsv"))
write_csv(bldb_raw, file.path(outdir, "bldb_raw.csv"))

write_tsv(bldb_joinable, file.path(outdir, "bldb_joinable.tsv"))
write_csv(bldb_joinable, file.path(outdir, "bldb_joinable.csv"))

write_tsv(bldb_aliases, file.path(outdir, "bldb_aliases.tsv"))
write_csv(bldb_aliases, file.path(outdir, "bldb_aliases.csv"))

write_tsv(bldb_lookup, file.path(outdir, "bldb_lookup.tsv"))
write_csv(bldb_lookup, file.path(outdir, "bldb_lookup.csv"))

message("Wrote outputs to: ", outdir)
message("Rows: raw=", nrow(bldb_raw), " joinable=", nrow(bldb_joinable),
        " aliases=", nrow(bldb_aliases), " lookup=", nrow(bldb_lookup))