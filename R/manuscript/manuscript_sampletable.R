library(readr)
library(dplyr)
library(tibble)
library(openxlsx)

source(file = "common.R")
source(file = "genotype/genotypecommon.R")
library(readr)
library(dplyr)
library(tibble)
library(openxlsx)
library(jsonlite)

source(file = "common.R")


normalize_sample_id <- function(x) {
  x <- as.character(x)
  
  x <- gsub("^sample", "", x)
  x <- gsub("^DIRE_EC_", "", x)
  x <- gsub("[^0-9]", "", x)
  
  sprintf("%03d", as.integer(x))
}


clean_empty <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x[x %in% c("NA", "NaN", "NULL", "None")] <- ""
  x
}


split_determinants <- function(x) {
  x <- clean_empty(x)
  
  if (length(x) == 0 || all(x == "")) {
    return(character(0))
  }
  
  y <- unlist(strsplit(x, ";", fixed = TRUE))
  y <- trimws(y)
  y <- y[y != ""]
  
  unique(y)
}


combine_determinants <- function(...) {
  x <- c(...)
  
  y <- unlist(lapply(x, split_determinants))
  y <- unique(y)
  y <- y[y != ""]
  
  if (length(y) == 0) {
    return("")
  }
  
  paste(y, collapse = "; ")
}


read_mlst_summary <- function(
    mlst_file = file.path(processedRoot, "python", "mlst", "mlst_ecoli_achtman.tsv")
) {
  if (!file.exists(mlst_file)) {
    stop("MLST file not found: ", mlst_file)
  }
  
  read_tsv(
    mlst_file,
    col_names = TRUE,
    col_types = cols(.default = "c")
  ) %>%
    transmute(
      sample = normalize_sample_id(sample),
      `Sequence type` = clean_empty(sequence_type)
    ) %>%
    mutate(
      `Sequence type` = if_else(
        `Sequence type` == "" | `Sequence type` == "-",
        "Not assigned",
        `Sequence type`
      )
    )
}


read_genotype_summary <- function() {
  genotype_table <- getGenotypeGroupTable() %>%
    select(
      sample,
      Functional_groups,
      Beta_genes,
      Quinolone_genes,
      Aminoglycoside_genes
    ) %>%
    mutate(
      sample = normalize_sample_id(sample),
      Functional_groups = clean_empty(Functional_groups),
      Beta_genes = clean_empty(Beta_genes),
      Quinolone_genes = clean_empty(Quinolone_genes),
      Aminoglycoside_genes = clean_empty(Aminoglycoside_genes)
    ) %>%
    rowwise() %>%
    mutate(
      `Resistance determinants` = combine_determinants(
        Beta_genes,
        Quinolone_genes,
        Aminoglycoside_genes
      )
    ) %>%
    ungroup() %>%
    transmute(
      sample,
      `Beta-lactamase functional class` = Functional_groups,
      `Beta-lactamase genes` = Beta_genes,
      `Quinolone resistance determinants` = Quinolone_genes,
      `Aminoglycoside resistance genes` = Aminoglycoside_genes,
      `Resistance determinants`
    )
  
  genotype_table
}


read_ena_sample_summary <- function(
    sample_json_file = file.path(processedRoot, "ENA", "samples.json")
) {
  if (!file.exists(sample_json_file)) {
    stop("ENA sample JSON file not found: ", sample_json_file)
  }
  
  ena_samples <- fromJSON(sample_json_file, flatten = TRUE)
  
  as_tibble(ena_samples) %>%
    transmute(
      sample = normalize_sample_id(name),
      `ENA sample alias` = name,
      `ENA sample accession` = accession
    )
}

read_sample_dates <- function(
    demographic_file = file.path(
      processedRootRcommon,
      "demographicstable.csv"
    )
) {
  if (!file.exists(demographic_file)) {
    stop("Demographics file not found: ", demographic_file)
  }
  
  demographic_table <- read_delim(
    demographic_file,
    delim = ";",
    col_types = cols(.default = "c"),
    trim_ws = TRUE
  )
  
  names(demographic_table) <- trimws(names(demographic_table))
  
  required_columns <- c("sample", "date")
  missing_columns <- setdiff(required_columns, names(demographic_table))
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing column(s) in demographics file: ",
      paste(missing_columns, collapse = ", "),
      ". Columns found: ",
      paste(names(demographic_table), collapse = ", ")
    )
  }
  
  demographic_table %>%
    transmute(
      sample = normalize_sample_id(.data[["sample"]]),
      `Collection date` = as.Date(.data[["date"]])
    )
}


read_ena_run_summary <- function(
    run_report_file = file.path(processedRoot, "ENA", "filereport_read_run_PRJEB115520.tsv")
) {
  if (!file.exists(run_report_file)) {
    warning("ENA run report not found: ", run_report_file)
    return(
      tibble(
        `ENA sample accession` = character(),
        `ENA experiment accession` = character(),
        `ENA run accession` = character()
      )
    )
  }
  
  read_tsv(
    run_report_file,
    col_types = cols(.default = "c")
  ) %>%
    transmute(
      `ENA sample accession` = sample_accession,
      `ENA experiment accession` = experiment_accession,
      `ENA run accession` = run_accession
    )
}


create_sample_summary <- function(
    output_csv = file.path(processedRootRcommon, "sample_summary.csv"),
    mlst_file = file.path(processedRoot, "python", "mlst", "mlst_ecoli_achtman.tsv"),
    sample_json_file = file.path(processedRoot, "ENA", "samples.json"),
    run_report_file = file.path(processedRoot, "ENA", "filereport_read_run_PRJEB115520.tsv"),
    demographic_file = file.path(processedRootRcommon, "demographicstable.csv"),
    exclude_samples = c("014", "038")
) {
  mlst_table <- read_mlst_summary(mlst_file)
  genotype_table <- read_genotype_summary()
  ena_samples <- read_ena_sample_summary(sample_json_file)
  ena_runs <- read_ena_run_summary(run_report_file)
  sample_dates <- read_sample_dates(demographic_file)
  
  sample_summary <- ena_samples %>%
    left_join(
      ena_runs,
      by = "ENA sample accession"
    ) %>%
    left_join(
      sample_dates,
      by = "sample"
    ) %>%
    left_join(
      mlst_table,
      by = "sample"
    ) %>%
    left_join(
      genotype_table,
      by = "sample"
    ) %>%
    mutate(
      `Sequence type` = if_else(
        is.na(`Sequence type`) | `Sequence type` == "",
        "Not assigned",
        `Sequence type`
      ),
      across(
        c(
          `ENA experiment accession`,
          `ENA run accession`,
          `Beta-lactamase functional class`,
          `Beta-lactamase genes`,
          `Quinolone resistance determinants`,
          `Aminoglycoside resistance genes`,
          `Resistance determinants`
        ),
        ~ clean_empty(.x)
      )
    ) %>%
    filter(!sample %in% exclude_samples) %>%
    arrange(sample) %>%
    transmute(
      `Sample ID` = sample,
      `Collection date`,
      `ENA sample alias`,
      `ENA sample accession`,
      `ENA experiment accession`,
      `ENA run accession`,
      `Sequence type`,
      `Beta-lactamase functional class`,
      `Beta-lactamase genes`,
      `Quinolone resistance determinants`,
      `Aminoglycoside resistance genes`,
      `Resistance determinants`
    )
  
  if (any(is.na(sample_summary$`Collection date`))) {
    missing_dates <- sample_summary %>%
      filter(is.na(`Collection date`)) %>%
      pull(`Sample ID`)
    
    stop(
      "Collection date missing for sample(s): ",
      paste(missing_dates, collapse = ", ")
    )
  }
  
  write_delim(
    sample_summary,
    output_csv,
    delim = ";",
    na = ""
  )
  
  message("Sample summary written to: ", output_csv)
  message("Number of isolates: ", nrow(sample_summary))
  
  invisible(sample_summary)
}


make_supplementary_table_s3 <- function(
    input_csv = file.path(processedRootRcommon, "sample_summary.csv"),
    output_xlsx = file.path(manuscriptDirectory, "Supplementary_Table_S3_isolate_collection_accessions_ST_resistance_determinants.xlsx"),
    output_csv = file.path(manuscriptDirectory, "Supplementary_Table_S3_isolate_collection_accessions_ST_resistance_determinants.csv")
) {
  supp_table_s3 <- read_delim(
    input_csv,
    delim = ";",
    col_types = cols(.default = "c")
  )
  
  supp_table_s3 <- supp_table_s3 %>%
    arrange(`Sample ID`)
  
  expected_columns <- c(
    "Sample ID",
    "Collection date",
    "ENA sample alias",
    "ENA sample accession",
    "ENA experiment accession",
    "ENA run accession",
    "Sequence type",
    "Beta-lactamase functional class",
    "Beta-lactamase genes",
    "Quinolone resistance determinants",
    "Aminoglycoside resistance genes",
    "Resistance determinants"
  )
  
  missing_columns <- setdiff(expected_columns, names(supp_table_s3))
  if (length(missing_columns) > 0) {
    stop(
      "Missing expected columns in input CSV: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  supp_table_s3 <- supp_table_s3 %>%
    select(all_of(expected_columns))
  
  column_descriptions_s3 <- tribble(
    ~Column, ~Description,
    "Sample ID", "Study-specific isolate identifier used in the manuscript.",
    "Collection date", "Date of collection of the clinical urine sample.",
    "ENA sample alias", "Study-specific ENA sample alias.",
    "ENA sample accession", "BioSample accession assigned by ENA.",
    "ENA experiment accession", "ENA experiment accession linked to the sequencing library and instrument metadata.",
    "ENA run accession", "ENA run accession linked to the deposited paired-end FASTQ files.",
    "Sequence type", "Multilocus sequence type according to the Achtman Escherichia coli scheme; \"Not assigned\" indicates that no sequence type was assigned.",
    "Beta-lactamase functional class", "Functional beta-lactamase class detected from genome data.",
    "Beta-lactamase genes", "Detected beta-lactamase genes.",
    "Quinolone resistance determinants", "Detected plasmid-mediated quinolone resistance genes and/or quinolone resistance-determining region mutations.",
    "Aminoglycoside resistance genes", "Detected aminoglycoside resistance genes.",
    "Resistance determinants", "Combined list of detected resistance genes and resistance-associated point mutations, including QRDR mutations."
  )
  
  readme_s3 <- tribble(
    ~Item, ~Description,
    "Table name",
    "Supplementary Table S3. Isolate collection dates, accessions, sequence types, and resistance determinants.",
    "Table description",
    "The table links isolate identifiers used in the manuscript to their collection dates and corresponding ENA accessions and summarizes sequence types and detected resistance determinants.",
    "Column descriptions",
    "Sheet 2 named 'Column descriptions' provides descriptions for columns in sheet 3 'Data'.",
    "Data",
    "Sheet 3 named 'Data' contains all table data.",
    "Notes",
    "QRDR mutations are included among resistance determinants. Constants such as organism, isolation source, and country are not included because they are identical for all isolates. Group and family-level resistance summaries are not included."
  )
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "README")
  writeData(wb, "README", readme_s3)

  addWorksheet(wb, "Column descriptions")
  writeData(wb, "Column descriptions", column_descriptions_s3)

  addWorksheet(wb, "Data")
  writeData(wb, "Data", supp_table_s3)
  
  
  
  header_style <- createStyle(
    textDecoration = "bold",
    border = "Bottom"
  )
  
  addStyle(
    wb, "Data", header_style,
    rows = 1,
    cols = seq_len(ncol(supp_table_s3)),
    gridExpand = TRUE
  )
  
  addStyle(
    wb, "Column descriptions", header_style,
    rows = 1,
    cols = seq_len(ncol(column_descriptions_s3)),
    gridExpand = TRUE
  )
  
  addStyle(
    wb, "README", header_style,
    rows = 1,
    cols = seq_len(ncol(readme_s3)),
    gridExpand = TRUE
  )
  
  freezePane(wb, "Data", firstRow = TRUE)
  freezePane(wb, "Column descriptions", firstRow = TRUE)
  freezePane(wb, "README", firstRow = TRUE)
  
  setColWidths(
    wb,
    "Data",
    cols = seq_len(ncol(supp_table_s3)),
    widths = "auto"
  )
  
  setColWidths(
    wb,
    "Column descriptions",
    cols = 1:2,
    widths = c(35, 120)
  )
  
  setColWidths(
    wb,
    "README",
    cols = 1:2,
    widths = c(30, 120)
  )
  
  saveWorkbook(
    wb,
    output_xlsx,
    overwrite = TRUE
  )
  
  write_delim(
    supp_table_s3,
    output_csv,
    delim = ";",
    na = ""
  )
  
  invisible(
    list(
      data = supp_table_s3,
      column_descriptions = column_descriptions_s3,
      readme = readme_s3,
      output_xlsx = output_xlsx,
      output_csv = output_csv
    )
  )
}


ALL <- function() {
  create_sample_summary()
  make_supplementary_table_s3()
  
  message("Supplementary Table s3 generated.")
  
  invisible(TRUE)
}

# ALL()