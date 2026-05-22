source(file = "common.R")

suppressPackageStartupMessages({
  library(rvest)
  library(xml2)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(janitor)
  library(tidyr)
})

outdir <- file.path(processedRootRassembly, "genotype")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pages <- tibble::tribble(
  ~prot, ~bldb_group, ~url,
  "A","Class A","http://bldb.eu/BLDB.php?prot=A",
  "B1","Subclass B1","http://bldb.eu/BLDB.php?prot=B1",
  "B2","Subclass B2","http://bldb.eu/BLDB.php?prot=B2",
  "B3","Subclass B3","http://bldb.eu/BLDB.php?prot=B3",
  "C","Class C","http://bldb.eu/BLDB.php?prot=C",
  "D","Class D","http://bldb.eu/BLDB.php?prot=D"
)

safe_read_html <- function(url) {
  
  tryCatch(
    read_html(url),
    error=function(e){
      warning("Could not read ",url)
      NULL
    }
  )
  
}

clean_table_names <- function(tab){
  
  tab |>
    janitor::clean_names() |>
    mutate(across(everything(), as.character)) |>
    mutate(across(where(is.character), str_squish))
  
}

extract_bush_jacoby <- function(x){
  
  x <- replace_na(x,"")
  
  if(x==""){
    return(NA_character_)
  }
  
  hits <- str_extract_all(x,"\\b[123][a-z]*\\b")[[1]]
  
  if(length(hits)==0){
    return(NA_character_)
  }
  
  paste(unique(hits), collapse=";")
  
}

read_bldb_table <- function(url,prot,bldb_group){
  
  page <- safe_read_html(url)
  
  if(is.null(page)){
    return(tibble())
  }
  
  tables <- html_table(page, fill=TRUE)
  
  if(length(tables)<2){
    return(tibble())
  }
  
  tab <- tables[[2]] |>
    clean_table_names()
  
  tab |>
    mutate(
      
      prot=prot,
      bldb_group=bldb_group,
      source_url=url,
      
      molecular_class=case_when(
        prot=="A" ~ "A",
        prot %in% c("B1","B2","B3") ~ "B",
        prot=="C" ~ "C",
        prot=="D" ~ "D",
        TRUE ~ NA_character_
      ),
      
      metallo_subclass=case_when(
        prot=="B1" ~ "B1",
        prot=="B2" ~ "B2",
        prot=="B3" ~ "B3",
        TRUE ~ NA_character_
      ),
      
      bush_jacoby_raw=phenotype,
      
      bush_jacoby_class=
        purrr::map_chr(
          phenotype,
          extract_bush_jacoby
        )
      
    ) |>
    filter(!is.na(proteinname))
  
}

bldb_raw <- purrr::pmap_dfr(
  pages,
  function(prot,bldb_group,url){
    read_bldb_table(url,prot,bldb_group)
  }
)

bldb_joinable <- bldb_raw |>
  mutate(
    
    proteinname=str_squish(proteinname),
    
    gene_symbol_guess=
      paste0(
        "bla",
        str_replace_all(proteinname,"\\s+","")
      ),
    
    gene_symbol_guess_no_prefix=
      str_replace(
        proteinname,
        regex("^bla",ignore_case=TRUE),
        ""
      )
    
  ) |>
  distinct()

bldb_aliases <- bldb_joinable |>
  select(
    proteinname,
    alternativeprotein_names,
    bldb_group,
    molecular_class,
    metallo_subclass,
    subfamily,
    phenotype,
    bush_jacoby_raw,
    bush_jacoby_class,
    gene_symbol_guess
  ) |>
  mutate(
    alternativeprotein_names=
      replace_na(alternativeprotein_names,"")
  ) |>
  separate_rows(
    alternativeprotein_names,
    sep="\\s*;\\s*|\\s*,\\s*"
  ) |>
  mutate(
    
    alternativeprotein_names=
      na_if(str_squish(alternativeprotein_names),""),
    
    alias_bla=
      ifelse(
        is.na(alternativeprotein_names),
        NA_character_,
        paste0(
          "bla",
          str_replace_all(
            alternativeprotein_names,
            "\\s+",""
          )
        )
      )
    
  ) |>
  filter(!is.na(alternativeprotein_names)) |>
  distinct()

bldb_lookup <- bldb_joinable |>
  select(
    proteinname,
    gene_symbol_guess,
    gene_symbol_guess_no_prefix,
    bldb_group,
    molecular_class,
    metallo_subclass,
    subfamily,
    phenotype,
    bush_jacoby_raw,
    bush_jacoby_class,
    source_url
  ) |>
  distinct()

write.csv2(
  bldb_raw,
  file.path(outdir,"bldb_raw.csv"),
  row.names=FALSE
)

write.csv2(
  bldb_joinable,
  file.path(outdir,"bldb_joinable.csv"),
  row.names=FALSE
)

write.csv2(
  bldb_aliases,
  file.path(outdir,"bldb_aliases.csv"),
  row.names=FALSE
)

write.csv2(
  bldb_lookup,
  file.path(outdir,"bldb_lookup.csv"),
  row.names=FALSE
)