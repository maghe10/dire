library(dplyr)
library(readr)
library(tidyr)
library(pheatmap)
library(factoextra)
library(ggplot2)
library(patchwork)
library(pvclust)
library(gtable) 
source("common.R")
source("genotype/genotypeCommon.R")
source("manuscript/generic_plot_helpers.R")
source("manuscript/manuscriptcommon.R")
# ============================================================
# Configuration
# ============================================================

PHENOTYPE_CLUSTER_LABEL <- "Phenotype cluster"
PHENOTYPE_ESBL_LABEL <- "Phenotype ESBL"

MODEL_ANTIBIOTICS <- ANTIBIOTICS
BETA_LACTAMS <- c("AMP", "AMC", "PIP", "TZP", "CAZ", "CRO", "CTX", "FEP")

cluster_palette <- function(n) {
  base <- c("blue", "green", "purple",
            "pink", "brown", "cyan","black")
  setNames(base[seq_len(n)], as.character(seq_len(n)))
}

# lower rank = placed higher up
ESBL_PRIORITY <- c(
  "C"   = 1,
  "A+M" = 2,
  "M"   = 3,
  "A"   = 4,
  "no"  = 5
)



gradient_red <- function(n = 100) colorRampPalette(c("white", "red"))(n)
gradient_orange <- function(n = 100) colorRampPalette(c("white", "orange"))(n)
gradient_green_yellow_red <- function(n = 199) colorRampPalette(c("red", "yellow", "green"))(n)
gradient_red_white_orange <- function(n = 199) {
  stopifnot(n %% 2 == 1)  # must be odd
  
  half <- (n - 1) / 2
  
  c(
    colorRampPalette(c("red", "white"))(half + 1)[- (half + 1)],
    "#FFFFFF",
    colorRampPalette(c("white", "orange"))(half + 1)[-1]
  )
}

# ============================================================
# IO
# ============================================================

read_semicolon_csv <- function(path) {
  read.csv2(path, check.names = FALSE, stringsAsFactors = FALSE)
}

load_cluster_data <- function(data_dir = ".") {
  list(
    atu                  = read_semicolon_csv(file.path(data_dir, "atu.csv")),
    breakpoints          = read_semicolon_csv(file.path(data_dir, "breakpoints.csv")),
    demographics         = read_semicolon_csv(file.path(data_dir, "demographicstable.csv")),
    esbl                 = read_semicolon_csv(file.path(data_dir, "esbl.csv")),
    mm_full              = read_semicolon_csv(file.path(data_dir, "millimetertable.csv")),
    mm_model             = read_semicolon_csv(file.path(data_dir, "modelmillimetertable.csv")),
    mm_full_rescaled     = read_semicolon_csv(file.path(data_dir, "millimetertableRescaled.csv")),
    mm_model_rescaled    = read_semicolon_csv(file.path(data_dir, "modelmillimetertableRescaled.csv")),
    sir                  = read_semicolon_csv(file.path(data_dir, "sirAntibioticsModel.csv")),
    sir_mode_a           = read_semicolon_csv(file.path(data_dir, "sirAntibioticsModel_Mode-A.csv")),
    sir_mode_b           = read_semicolon_csv(file.path(data_dir, "sirAntibioticsModel_Mode-B.csv")),
    sir_mode_c           = read_semicolon_csv(file.path(data_dir, "sirAntibioticsModel_Mode-C.csv"))
  )
}

# ============================================================
# Small utilities
# ============================================================

as_sample_rownames <- function(df) {
  stopifnot("sample" %in% colnames(df))
  rn <- df$sample
  df <- as.data.frame(df)
  rownames(df) <- rn
  df
}

select_model_antibiotics <- function(df) {
  df %>%
    as.data.frame() %>%
    select(sample, all_of(MODEL_ANTIBIOTICS))
}

drop_sample_column <- function(df) {
  df %>% select(-sample)
}

safe_rescale_01 <- function(df_numeric) {
  out <- df_numeric
  for (j in seq_len(ncol(out))) {
    x <- out[[j]]
    rng <- range(x, na.rm = TRUE)
    if (is.infinite(rng[1]) || is.infinite(rng[2]) || isTRUE(all.equal(rng[1], rng[2]))) {
      out[[j]] <- 0
    } else {
      out[[j]] <- (x - rng[1]) / (rng[2] - rng[1])
    }
  }
  out
}

sir_to_numeric <- function(sir_df) {
  out <- sir_df
  out[out == "S"] <- 1
  out[out == "I"] <- 0.5
  out[out == "R"] <- 0
  out[] <- lapply(out, as.numeric)
  out
}

wrap_pheatmap <- function(ph) {
  wrap_pheatmap_tight(ph)
  #  patchwork::wrap_elements(full = ph$gtable, clip = FALSE)
}

wrap_pheatmap_tight <- function(ph) {
  g <- ph$gtable
  
  # remove empty outer rows/cols as much as possible
  g <- gtable::gtable_trim(g)
  
  patchwork::wrap_elements(full = g, clip = FALSE)
}



panel_box_theme <- theme(
  plot.background = element_rect(fill = NA, colour = NA),
  plot.margin = margin(4, 4, 4, 4)
)

# ============================================================
# Data extractors
# ============================================================

read_millimeter_table <- function(dat) {
  dat$mm_model %>%
    select_model_antibiotics() %>%
    as_sample_rownames()
}

read_rescaled_millimeter_table <- function(dat, use_existing = TRUE) {
  if (use_existing) {
    dat$mm_model_rescaled %>%
      select_model_antibiotics() %>%
      as_sample_rownames()
  } else {
    mm <- read_millimeter_table(dat)
    mm_only <- drop_sample_column(mm)
    mm_scaled <- safe_rescale_01(mm_only)
    out <- cbind(sample = rownames(mm), mm_scaled)
    as_sample_rownames(out)
  }
}

read_sir_antibiotics_model <- function(dat, mode = NA) {
  sir_tbl <- if (is.na(mode)) {
    dat$sir
  } else {
    switch(
      mode,
      "Mode-A" = dat$sir_mode_a,
      "Mode-B" = dat$sir_mode_b,
      "Mode-C" = dat$sir_mode_c,
      stop("Unknown mode: ", mode)
    )
  }
  
  sir_tbl %>%
    select_model_antibiotics() %>%
    as_sample_rownames()
}

read_atu_antibiotics_model <- function(dat) {
  dat$atu %>%
    select_model_antibiotics() %>%
    as_sample_rownames()
}

read_esbl_annotation <- function(dat) {
  esbl <- dat$esbl %>%
    mutate(
      ESBL = case_when(
        ESBL_CARBA ~ "C",
        ESBL_A & ESBL_M ~ "A+M",
        ESBL_M ~ "M",
        ESBL_A ~ "A",
        TRUE ~ "no"
      )
    ) %>%
    select(sample, ESBL)
  
  rownames(esbl) <- esbl$sample
  esbl$sample <- NULL
  colnames(esbl)[1] <- PHENOTYPE_ESBL_LABEL
  esbl
}

fetch_tables <- function(dat, use_existing_rescaled = TRUE) {
  mm <- read_millimeter_table(dat)
  mm_rescaled <- read_rescaled_millimeter_table(dat, use_existing = use_existing_rescaled)
  sir <- read_sir_antibiotics_model(dat)
  sir_num <- sir_to_numeric(drop_sample_column(sir))
  prediction_errors <- fetchPredictionErrors()
  
  list(
    SIR = sir_num,
    RESCALEDMM = drop_sample_column(mm_rescaled),
    MM = drop_sample_column(mm),
    ERRORS = prediction_errors
  )
}

# ============================================================
# Clustering + annotations
# ============================================================
order_rows_by_cluster_esbl_resistance <- function(
    caa,
    mat,
    cluster_order = NULL,
    esbl_priority = c("C" = 1, "A+M" = 2, "M" = 3, "A" = 4, "no" = 5)
) {
  stopifnot(all(rownames(caa$annotation) %in% rownames(mat)))
  
  ann <- caa$annotation
  mat <- mat[rownames(ann), , drop = FALSE]
  
  # lower mean == more resistant
  isolate_mean <- rowMeans(mat, na.rm = TRUE)
  
  # original tree order as final tie-breaker
  order_names <- rownames(caa$data)[caa$hc$order]
  rank <- match(rownames(ann), order_names)
  
  ann2 <- ann %>%
    mutate(
      sample = rownames(.),
      .cluster = as.character(.data[[PHENOTYPE_CLUSTER_LABEL]]),
      .esbl = as.character(.data[[PHENOTYPE_ESBL_LABEL]]),
      .esbl_rank = unname(esbl_priority[.esbl]),
      .isolate_mean = isolate_mean,
      .rank = rank
    )
  
  # fallback for unknown ESBL labels
  ann2$.esbl_rank[is.na(ann2$.esbl_rank)] <- max(esbl_priority, na.rm = TRUE) + 1
  
  # if cluster order not supplied, keep numeric/alphabetic cluster order
  if (is.null(cluster_order)) {
    cl <- unique(ann2$.cluster)
    suppressWarnings(cl_num <- as.numeric(cl))
    if (all(!is.na(cl_num))) {
      cluster_order <- as.character(sort(cl_num))
    } else {
      cluster_order <- sort(cl)
    }
  } else {
    cluster_order <- as.character(cluster_order)
  }
  
  # any clusters not listed go last
  missing_clusters <- setdiff(unique(ann2$.cluster), cluster_order)
  cluster_order <- c(cluster_order, sort(missing_clusters))
  
  cluster_rank <- setNames(seq_along(cluster_order), cluster_order)
  ann2$.cluster_rank <- unname(cluster_rank[ann2$.cluster])
  
  ann2 <- ann2 %>%
    arrange(
      .cluster_rank,       # user-defined cluster order
      .esbl_rank,          # C first within cluster
      .isolate_mean,    # lowest mm first within ESBL
      .rank                # tie-breaker
    )
  
  ann2$sample
}
  

# order_rows_by_cluster_esbl_group_error <- function(caa,
#                                                    error_tbl,
#                                                    cols,
#                                                    cluster_order = NULL,
#                                                    esbl_priority = ESBL_PRIORITY) {
#   ann <- caa$annotation
#   error_tbl <- error_tbl[rownames(ann), cols, drop = FALSE]
#   
#   mean_abs_error <- rowMeans(abs(error_tbl), na.rm = TRUE)
#   
#   order_names <- rownames(caa$data)[caa$hc$order]
#   tree_rank <- match(rownames(ann), order_names)
#   
#   ann2 <- ann %>%
#     mutate(
#       sample = rownames(.),
#       .cluster = as.character(Cluster),
#       .esbl = as.character(ESBL),
#       .esbl_rank = unname(esbl_priority[.esbl]),
#       .mean_abs_error = mean_abs_error,
#       .tree_rank = tree_rank
#     )
#   
#   ann2$.esbl_rank[is.na(ann2$.esbl_rank)] <- max(esbl_priority) + 1
#   
#   if (is.null(cluster_order)) {
#     cluster_order <- sort(unique(ann2$.cluster))
#   }
#   
#   cluster_rank <- setNames(seq_along(cluster_order), as.character(cluster_order))
#   
#   ann2 %>%
#     mutate(.cluster_rank = unname(cluster_rank[.cluster])) %>%
#     arrange(
#       .cluster_rank,
#       .esbl_rank,
#       desc(.mean_abs_error),
#       .tree_rank
#     ) %>%
#     pull(sample)
# }

order_rows_by_group_error <- function(error_tbl, cols, decreasing = TRUE) {
  stopifnot(all(cols %in% colnames(error_tbl)))
  
  score <- rowMeans(abs(error_tbl[, cols, drop = FALSE]), na.rm = TRUE)
  
  names(sort(score, decreasing = decreasing))
}

order_rows_by_error <- function(error_tbl, decreasing = TRUE) {

  score <- rowMeans(abs(error_tbl[,, drop = FALSE]), na.rm = TRUE)
  
  names(sort(score, decreasing = decreasing))
}



renumber_clusters_by_order <- function(caa, cluster_order) {
  ann <- caa$annotation
  
  cluster_order <- as.character(cluster_order)
  current_clusters <- as.character(unique(ann[[PHENOTYPE_CLUSTER_LABEL]]))
  
  # add any missing clusters at the end
  missing <- setdiff(current_clusters, cluster_order)
  cluster_order <- c(cluster_order, sort(missing))
  
  # mapping old → new labels
  new_labels <- seq_along(cluster_order)
  names(new_labels) <- cluster_order
  
  ann[[PHENOTYPE_CLUSTER_LABEL]] <- factor(
    new_labels[as.character(ann[[PHENOTYPE_CLUSTER_LABEL]])],
    levels = seq_along(cluster_order)
  )
  
  caa$annotation <- ann
  caa$clusters <- as.integer(as.character(ann[[PHENOTYPE_CLUSTER_LABEL]]))
  
  # update colors to match new numbering
  caa$annotation_colors[[PHENOTYPE_CLUSTER_LABEL]] <- cluster_palette(length(cluster_order))
  
  caa
}


clusters_and_annotation_umap <- function(data_frame, k) {
  data_mat <- as.matrix(data_frame)
  print(head(data_frame))
  set.seed(123)
  um <- uwot::umap(
    data_mat,
    n_components = 2
  )
  rownames(um) <- rownames(data_mat)
  
  hc_umap <- hclust(dist(um), method = "ward.D2")
  clusters_umap <- cutree(hc_umap, k = k)
  
  annotation <- data.frame(
    Cluster = factor(as.character(clusters_umap), levels = as.character(seq_len(k))),
    row.names = rownames(data_frame)
  )
  
  annotation_colors <- list(
    Cluster = cluster_palette(k)
  )
  
  list(
    data = data_mat,
    umap = um,
    hc = hc_umap,
    clusters = clusters_umap,
    annotation = annotation,
    annotation_colors = annotation_colors
  )
}



clusters_and_annotation_pca <- function(data_frame, k) {
  data_mat <- as.matrix(data_frame)

  hc <- hclust(dist(data_mat), method = "ward.D2")
  clusters <- cutree(hc, k = k)

  annotation <- data.frame(
    row.names = rownames(data_frame)
  )
  
  annotation[[PHENOTYPE_CLUSTER_LABEL]] <- factor(
    clusters,
    levels = as.character(seq_len(k))
  )
  
  annotation_colors <- list()
  annotation_colors[[PHENOTYPE_CLUSTER_LABEL]] <- cluster_palette(k)
  

  list(
    data = data_mat,
    hc = hc,
    clusters = clusters,
    annotation = annotation,
    annotation_colors = annotation_colors
  )
}

clusters_and_annotation <- function(data_frame, k) {
  caa <- clusters_and_annotation_pca(data_frame, k)
  
  caa
}



decorate_with_esbl <- function(caa, dat) {
  esbl <- read_esbl_annotation(dat)
  
  esbl <- esbl[rownames(caa$annotation), , drop = FALSE]
  caa$annotation[[PHENOTYPE_ESBL_LABEL]] <- factor(
    esbl[[PHENOTYPE_ESBL_LABEL]],
    levels = c("A", "M", "A+M", "C", "no")
  )
  
  caa$annotation_colors[[PHENOTYPE_ESBL_LABEL]] <- c(
    "C" = "darkred",
    "A+M" = "red",
    "M" = "orange",
    "A" = "yellow",
    "no" = "white"
  )
  caa
}



decorate_with_quinolone_family_presence <- function(caa, genotypeGroupTable) {
  
  qf <- genotypeGroupTable %>%
    dplyr::select(sample, Quinolone_families) %>%
    dplyr::mutate(
      sample = as.character(sample),
      Quinolone_families = dplyr::if_else(
        is.na(Quinolone_families), "", Quinolone_families
      )
    )
  
  rownames(qf) <- qf$sample
  qf <- qf[rownames(caa$annotation), , drop = FALSE]
  
  all_families <- qf$Quinolone_families %>%
    strsplit(";", fixed = TRUE) %>%
    unlist() %>%
    trimws() %>%
    unique()
  
  all_families <- sort(all_families[all_families != ""])
  added_cols <- character(0)
  
  for (fam in all_families) {
    present <- grepl(fam, qf$Quinolone_families, fixed = TRUE)
    col_name <- fam
    
    caa$annotation[[col_name]] <- factor(
      ifelse(present, "yes", "no"),
      levels = c("yes", "no")
    )
    
    yes_color <- if (fam == "QRDR") "orange" else "yellow"
    
    caa$annotation_colors[[col_name]] <- c(
      "yes" = yes_color,
      "no" = "white"
    )
    
    added_cols <- c(added_cols, col_name)
  }
  
  caa <- reorder_presence_annotations_by_color(caa, added_cols)
  caa
}


decorate_with_aminoglycoside_family_presence <- function(caa, genotypeGroupTable) {
  
  agf <- genotypeGroupTable %>%
    dplyr::select(sample, Aminoglycoside_families) %>%
    dplyr::mutate(
      sample = as.character(sample),
      Aminoglycoside_families = dplyr::if_else(
        is.na(Aminoglycoside_families), "", Aminoglycoside_families
      )
    )
  
  rownames(agf) <- agf$sample
  agf <- agf[rownames(caa$annotation), , drop = FALSE]
  
  all_families <- agf$Aminoglycoside_families %>%
    strsplit(";", fixed = TRUE) %>%
    unlist() %>%
    trimws() %>%
    unique()
  
  all_families <- sort(all_families[all_families != ""])
  added_cols <- character(0)
  
  for (fam in all_families) {
    present <- grepl(fam, agf$Aminoglycoside_families, fixed = TRUE)
    col_name <- fam
    
    caa$annotation[[col_name]] <- factor(
      ifelse(present, "yes", "no"),
      levels = c("yes", "no")
    )
    
    yes_color <- dplyr::case_when(
      fam %in% c("AAC(3)", "AAC(6)") ~ "red",
      grepl("^APH", fam) ~ "orange",
      fam %in% c("AAD", "ANT") ~ "yellow",
      TRUE ~ "grey70"
    )
    
    caa$annotation_colors[[col_name]] <- c(
      "yes" = yes_color,
      "no" = "white"
    )
    
    added_cols <- c(added_cols, col_name)
  }
  
  caa <- reorder_presence_annotations_by_color(caa, added_cols)
  caa
}

decorate_with_functional_group_presence <- function(caa, genotypeGroupTable) {
  
  fg <- genotypeGroupTable %>%
    dplyr::select(sample, Functional_groups) %>%
    dplyr::mutate(
      sample = as.character(sample),
      Functional_groups = dplyr::if_else(
        is.na(Functional_groups), "", Functional_groups
      )
    )
  
  rownames(fg) <- fg$sample
  fg <- fg[rownames(caa$annotation), , drop = FALSE]
  
  all_groups <- fg$Functional_groups %>%
    strsplit(";", fixed = TRUE) %>%
    unlist() %>%
    trimws() %>%
    unique()
  
  all_groups <- sort(all_groups[all_groups != ""])
  added_cols <- character(0)
  
  get_group_color <- function(grp) {
    if (grepl("^3", grp)) return("darkred")
    if (grp == "1") return("orange")
    if (grp == "2be") return("yellow")
    if (grp %in% c("2b", "2br", "2d")) return("grey70")
    "white"
  }
  
  for (grp in all_groups) {
    present <- grepl(grp, fg$Functional_groups, fixed = TRUE)
    col_name <- grp
    
    caa$annotation[[col_name]] <- factor(
      ifelse(present, "yes", "no"),
      levels = c("yes", "no")
    )
    
    caa$annotation_colors[[col_name]] <- c(
      "yes" = get_group_color(grp),
      "no" = "white"
    )
    
    added_cols <- c(added_cols, col_name)
  }
  
  caa <- reorder_presence_annotations_by_color(caa, added_cols)
  caa
}

decorate_with_beta_gene_presence <- function(caa,
                                             genotypeGroupTable,
                                             mapping_table,
                                             gene_col = "Gene_symbol",
                                             functional_group_col = "Functional_group") {
  
  bg <- genotypeGroupTable %>%
    dplyr::select(sample, Beta_genes) %>%
    dplyr::mutate(
      sample = as.character(sample),
      Beta_genes = dplyr::if_else(
        is.na(Beta_genes), "", Beta_genes
      )
    )
  
  rownames(bg) <- bg$sample
  bg <- bg[rownames(caa$annotation), , drop = FALSE]
  
  # Build gene -> functional group lookup
  gene_map <- mapping_table %>%
    dplyr::select(
      gene = dplyr::all_of(gene_col),
      functional_group = dplyr::all_of(functional_group_col)
    ) %>%
    dplyr::mutate(
      gene = trimws(as.character(gene)),
      functional_group = trimws(as.character(functional_group))
    )
  
  gene_to_fg <- stats::setNames(gene_map$functional_group, gene_map$gene)
  
  all_genes <- bg$Beta_genes %>%
    strsplit(";", fixed = TRUE) %>%
    unlist() %>%
    trimws() %>%
    unique()
  
  all_genes <- sort(all_genes[all_genes != ""])
  added_cols <- character(0)
  
  get_fg_color <- function(fg) {
    if (is.na(fg) || fg == "" || fg == "NA") return("white")
    if (grepl("^3", fg)) return("darkred")
    if (fg == "1") return("orange")
    if (fg == "2be") return("yellow")
    if (fg %in% c("2b", "2br", "2d")) return("grey70")
    "white"
  }
  
  for (gene in all_genes) {
    present <- grepl(gene, bg$Beta_genes, fixed = TRUE)
    col_name <- gene
    
    fg <- unname(gene_to_fg[gene])
    yes_color <- get_fg_color(fg)
    
    caa$annotation[[col_name]] <- factor(
      ifelse(present, "yes", "no"),
      levels = c("yes", "no")
    )
    
    caa$annotation_colors[[col_name]] <- c(
      "yes" = yes_color,
      "no" = "white"
    )
    
    added_cols <- c(added_cols, col_name)
  }
  
  caa <- reorder_presence_annotations_by_color(caa, added_cols)
  caa
}

decorate_with_aminoglycoside_gene_presence <- function(caa,
                                                       genotypeGroupTable,
                                                       geneMappingTable) {
  
  ag <- genotypeGroupTable %>%
    dplyr::select(sample, Aminoglycoside_genes) %>%
    dplyr::mutate(
      sample = as.character(sample),
      Aminoglycoside_genes = dplyr::if_else(
        is.na(Aminoglycoside_genes), "", Aminoglycoside_genes
      )
    )
  
  rownames(ag) <- ag$sample
  ag <- ag[rownames(caa$annotation), , drop = FALSE]
  
  gene_map <- geneMappingTable %>%
    dplyr::select(Gene_symbol, class, subclass) %>%
    dplyr::mutate(
      Gene_symbol = trimws(as.character(Gene_symbol)),
      class = trimws(as.character(class)),
      subclass = trimws(as.character(subclass))
    )
  
  gene_to_class <- stats::setNames(gene_map$class, gene_map$Gene_symbol)
  gene_to_subclass <- stats::setNames(gene_map$subclass, gene_map$Gene_symbol)
  
  all_genes <- ag$Aminoglycoside_genes %>%
    strsplit(";", fixed = TRUE) %>%
    unlist() %>%
    trimws() %>%
    unique()
  
  all_genes <- sort(all_genes[all_genes != ""])
  added_cols <- character(0)
  
  get_ag_color <- function(gene) {
    g <- trimws(gene)
    g_lower <- tolower(g)
    
    mapped_class <- toupper(unname(gene_to_class[g]))
    mapped_subclass <- toupper(unname(gene_to_subclass[g]))
    
    if (grepl("^aac\\(", g_lower)) return("red")
    if (grepl("^aph\\(", g_lower)) return("orange")
    if (grepl("^aad", g_lower)) return("yellow")
    if (grepl("^ant\\(", g_lower)) return("yellow")
    
    # fallback from mapped text if raw naming is odd
    if (grepl("AMINOGLYCOSIDE", mapped_class, fixed = TRUE)) {
      if (grepl("AAC", toupper(g), fixed = TRUE)) return("red")
      if (grepl("APH", toupper(g), fixed = TRUE)) return("orange")
      if (grepl("AAD", toupper(g), fixed = TRUE)) return("yellow")
      if (grepl("ANT", toupper(g), fixed = TRUE)) return("yellow")
    }
    
    if (grepl("AAC", mapped_subclass, fixed = TRUE)) return("red")
    if (grepl("APH", mapped_subclass, fixed = TRUE)) return("orange")
    if (grepl("AAD", mapped_subclass, fixed = TRUE)) return("yellow")
    if (grepl("ANT", mapped_subclass, fixed = TRUE)) return("yellow")
    
    "grey70"
  }
  
  for (gene in all_genes) {
    present <- grepl(gene, ag$Aminoglycoside_genes, fixed = TRUE)
    col_name <- gene
    
    caa$annotation[[col_name]] <- factor(
      ifelse(present, "yes", "no"),
      levels = c("yes", "no")
    )
    
    caa$annotation_colors[[col_name]] <- c(
      "yes" = get_ag_color(gene),
      "no" = "white"
    )
    
    added_cols <- c(added_cols, col_name)
  }
  
  caa <- reorder_presence_annotations_by_color(caa, added_cols)
  caa
}

decorate_with_quinolone_gene_presence <- function(caa,
                                                  genotypeGroupTable,
                                                  geneMappingTable) {
  
  qg <- genotypeGroupTable %>%
    dplyr::select(sample, Quinolone_genes) %>%
    dplyr::mutate(
      sample = as.character(sample),
      Quinolone_genes = dplyr::if_else(
        is.na(Quinolone_genes), "", Quinolone_genes
      )
    )
  
  rownames(qg) <- qg$sample
  qg <- qg[rownames(caa$annotation), , drop = FALSE]
  
  gene_map <- geneMappingTable %>%
    dplyr::select(Gene_symbol, class, subclass) %>%
    dplyr::mutate(
      Gene_symbol = trimws(as.character(Gene_symbol)),
      class = trimws(as.character(class)),
      subclass = trimws(as.character(subclass))
    )
  
  gene_to_class <- stats::setNames(gene_map$class, gene_map$Gene_symbol)
  gene_to_subclass <- stats::setNames(gene_map$subclass, gene_map$Gene_symbol)
  
  all_genes <- qg$Quinolone_genes %>%
    strsplit(";", fixed = TRUE) %>%
    unlist() %>%
    trimws() %>%
    unique()
  
  all_genes <- sort(all_genes[all_genes != ""])
  added_cols <- character(0)
  
  get_q_color <- function(gene) {
    g <- trimws(gene)
    g_upper <- toupper(g)
    
    mapped_class <- toupper(unname(gene_to_class[g]))
    mapped_subclass <- toupper(unname(gene_to_subclass[g]))
    
    # QRDR mutations often won't be in the gene map as AMR genes
    if (grepl("GYRA_|PARC_|PARE_|GYRB_", g_upper)) return("orange")
    if (g_upper == "QRDR") return("orange")
    
    # PMQR genes / mixed AMINOGLYCOSIDE-QUINOLONE genes
    if (grepl("^QNR", g_upper)) return("yellow")
    if (grepl("^QEPA", g_upper)) return("yellow")
    if (grepl("QUINOLONE", mapped_class, fixed = TRUE)) return("yellow")
    if (grepl("QUINOLONE", mapped_subclass, fixed = TRUE)) return("yellow")
    
    "grey70"
  }
  
  for (gene in all_genes) {
    present <- grepl(gene, qg$Quinolone_genes, fixed = TRUE)
    col_name <- gene
    
    caa$annotation[[col_name]] <- factor(
      ifelse(present, "yes", "no"),
      levels = c("yes", "no")
    )
    
    caa$annotation_colors[[col_name]] <- c(
      "yes" = get_q_color(gene),
      "no" = "white"
    )
    
    added_cols <- c(added_cols, col_name)
  }
  
  caa <- reorder_presence_annotations_by_color(caa, added_cols)
  caa
}


decorate_with_mean_absolute_error <- function(
    caa,
    error_tbl,
    cols = colnames(error_tbl),
    label = "Mean error rate"
) {
  x <- error_tbl[rownames(caa$annotation), cols, drop = FALSE]
  
  mean_error <- rowMeans(abs(as.matrix(x)), na.rm = TRUE)
  mean_error[is.nan(mean_error)] <- 0
  mean_error <- pmin(pmax(mean_error, 0), 1)
  
  caa$annotation[[label]] <- mean_error
  
  caa$annotation <- caa$annotation[
    , c(label, setdiff(colnames(caa$annotation), label)),
    drop = FALSE
  ]
  
  full_bw <- colorRampPalette(c("white", "black"))(200)
  
  min_idx <- max(1, round(min(mean_error, na.rm = TRUE) * 199) + 1)
  max_idx <- min(200, round(max(mean_error, na.rm = TRUE) * 199) + 1)
  
  caa$annotation_colors[[label]] <- full_bw[min_idx:max_idx]
  
  caa
}

keep_annotation_legend <- function(caa, keep) {
  caa$annotation <- caa$annotation[, keep, drop = FALSE]
  caa$annotation_colors <- caa$annotation_colors[keep]
  caa
}

reorder_presence_annotations_by_color <- function(caa,
                                                  new_cols,
                                                  color_order = c("grey70", "yellow", "orange", "red", "darkred")) {
  new_cols <- intersect(new_cols, colnames(caa$annotation))
  if (length(new_cols) == 0) return(caa)
  
  yes_cols <- sapply(new_cols, function(col) {
    col_map <- caa$annotation_colors[[col]]
    if (!is.null(col_map) && "yes" %in% names(col_map)) unname(col_map[["yes"]]) else NA_character_
  }, USE.NAMES = TRUE)
  
  color_rank <- match(yes_cols, color_order)
  color_rank[is.na(color_rank)] <- length(color_order) + 1
  
  ord <- order(color_rank, new_cols)
  new_cols_ordered <- new_cols[ord]
  
  old_cols <- setdiff(colnames(caa$annotation), new_cols)
  caa$annotation <- caa$annotation[, c(old_cols, new_cols_ordered), drop = FALSE]
  
  caa
}


aggregate_presence_by_color <- function(caa,
                                        decorate_fun,
                                        ...,
                                        label_map = NULL,
                                        prefix = NULL,
                                        drop_original = TRUE,
                                        no_color = "white") {
  stopifnot(is.function(decorate_fun))
  
  old_ann_cols <- colnames(caa$annotation)
  old_color_cols <- names(caa$annotation_colors)
  
  # Apply decorator
  caa2 <- decorate_fun(caa, ...)
  
  new_ann_cols <- setdiff(colnames(caa2$annotation), old_ann_cols)
  new_color_cols <- setdiff(names(caa2$annotation_colors), old_color_cols)
  target_cols <- intersect(new_ann_cols, new_color_cols)
  
  if (length(target_cols) == 0) {
    warning("Decorator did not add any new colored annotation columns.")
    return(caa2)
  }
  
  # yes-color for each new presence column
  yes_colors <- sapply(target_cols, function(col) {
    col_map <- caa2$annotation_colors[[col]]
    if (!("yes" %in% names(col_map))) return(NA_character_)
    unname(col_map[["yes"]])
  }, USE.NAMES = TRUE)
  
  keep <- !is.na(yes_colors) & yes_colors != no_color
  target_cols <- target_cols[keep]
  yes_colors <- yes_colors[keep]
  
  if (length(target_cols) == 0) {
    warning("No aggregatable colored presence columns found after excluding no_color.")
    return(caa2)
  }
  
  color_groups <- split(target_cols, yes_colors)
  rename_after_drop <- character(0)
  aggregated_cols <- character(0)
  
  make_unique_name <- function(x, existing) {
    if (!(x %in% existing)) return(x)
    i <- 2
    cand <- paste0(x, "_", i)
    while (cand %in% existing) {
      i <- i + 1
      cand <- paste0(x, "_", i)
    }
    cand
  }
  
  for (clr in names(color_groups)) {
    cols <- color_groups[[clr]]
    
    desired_name <- if (!is.null(label_map) && clr %in% names(label_map)) {
      label_map[[clr]]
    } else {
      clr
    }
    
    if (!is.null(prefix)) {
      desired_name <- paste0(prefix, desired_name)
    }
    
    if (isTRUE(drop_original) && desired_name %in% cols) {
      out_name <- make_unique_name(
        paste0(".__tmp__.", desired_name),
        colnames(caa2$annotation)
      )
      rename_after_drop[out_name] <- desired_name
    } else {
      out_name <- make_unique_name(desired_name, colnames(caa2$annotation))
    }
    
    any_yes <- apply(
      caa2$annotation[, cols, drop = FALSE],
      1,
      function(x) any(as.character(x) == "yes", na.rm = TRUE)
    )
    
    caa2$annotation[[out_name]] <- factor(
      ifelse(any_yes, "yes", "no"),
      levels = c("yes", "no")
    )
    
    caa2$annotation_colors[[out_name]] <- c(
      "yes" = clr,
      "no" = no_color
    )
    
    aggregated_cols <- c(aggregated_cols, out_name)
  }
  
  if (isTRUE(drop_original)) {
    caa2$annotation <- caa2$annotation[, setdiff(colnames(caa2$annotation), target_cols), drop = FALSE]
    caa2$annotation_colors[target_cols] <- NULL
    
    if (length(rename_after_drop) > 0) {
      for (tmp_name in names(rename_after_drop)) {
        final_name <- rename_after_drop[[tmp_name]]
        
        colnames(caa2$annotation)[colnames(caa2$annotation) == tmp_name] <- final_name
        caa2$annotation_colors[[final_name]] <- caa2$annotation_colors[[tmp_name]]
        caa2$annotation_colors[[tmp_name]] <- NULL
        
        aggregated_cols[aggregated_cols == tmp_name] <- final_name
      }
    }
  }
  
  # order aggregated columns by yes-color: grey70, yellow, orange, red, darkred
  caa2 <- reorder_presence_annotations_by_color(caa2, aggregated_cols)
  
  caa2
}

# ============================================================
# Plot helpers
# make_annotation_name_strip <- function(labels, angle = 90, base_size = 10) {
#   df <- data.frame(
#     x = seq_along(labels),
#     y = 1,
#     lab = labels,
#     stringsAsFactors = FALSE
#   )
#   
#   ggplot2::ggplot(df, ggplot2::aes(x, y, label = lab)) +
#     ggplot2::geom_text(
#       angle = angle,
#       hjust = 1,
#       vjust = 0.5,
#       size = base_size / 3
#     ) +
#     ggplot2::scale_x_continuous(
#       limits = c(0.5, length(labels) + 0.5),
#       expand = c(0, 0)
#     ) +
#     ggplot2::coord_cartesian(clip = "off") +
#     ggplot2::theme_void() +
#     ggplot2::theme(
#       plot.margin = ggplot2::margin(0, 2, 0, 2)
#     )
# }

# get_annotation_track_names <- function(annotation_row, exclude = "Cluster") {
#   labs <- colnames(annotation_row)
#   labs <- setdiff(labs, exclude)
#   c(labs, if (exclude %in% colnames(annotation_row)) exclude)
# }


# wrap_pheatmap_with_annotation_names <- function(ph, annotation_labels, tag = NULL, strip_height = 0.10) {
#   body <- patchwork::wrap_elements(full = ph$gtable, clip = FALSE)
#   
#   if (!is.null(tag)) {
#     body <- body +
#       ggplot2::labs(title = tag) +
#       ggplot2::theme(
#         plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 18),
#         plot.margin = ggplot2::margin(2, 2, 2, 2)
#       )
#   }
#   
#   strip <- make_annotation_name_strip(annotation_labels)
#   
#   body / strip + patchwork::plot_layout(heights = c(1, strip_height))
# }





make_cluster_plot_pca <- function(caa, title = NULL) {
  p <- fviz_cluster(
    list(data = caa$data, cluster = caa$clusters),
    palette = unlist(caa$annotation_colors[[PHENOTYPE_CLUSTER_LABEL]]),
    ellipse.type = "convex",
    geom = "point",
    show.clust.cent = TRUE
  ) +
    labs(
      colour = PHENOTYPE_CLUSTER_LABEL,
      shape = PHENOTYPE_CLUSTER_LABEL,
      fill = PHENOTYPE_CLUSTER_LABEL
    ) +
    ggtitle(title %||% "")
  
  # extract % variance from PCA
  pca <- prcomp(caa$data, scale. = TRUE)
  var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  
  p + 
    labs(
      x = sprintf("PC1 (%.1f%%)", var[1]),
      y = sprintf("PC2 (%.1f%%)", var[2])
    )
  
}

make_cluster_plot_umap <- function(caa, title = NULL, seed = 123) {
  set.seed(seed)
  
  um <- uwot::umap(
    as.matrix(caa$data),
    n_components = 2
  )
  
  colnames(um) <- c("UMAP1", "UMAP2")
  
  fviz_cluster(
    list(data = um, cluster = caa$clusters),
    palette = unlist(caa$annotation_colors$Cluster),
    ellipse.type = "convex",
    geom = "point",
    show.clust.cent = TRUE
  ) +
    labs(
      title = title,
      x = "UMAP1",
      y = "UMAP2"
    )
}


make_cluster_plot <- function(caa, title = NULL) {
  make_cluster_plot_pca(caa,title)
}


`%||%` <- function(x, y) if (is.null(x)) y else x




reorder_hclust_by_weights <- function(hc, weights) {
  dend <- as.dendrogram(hc)
  dend <- reorder(dend, wts = weights, agglo.FUN = mean)
  as.hclust(dend)
}

make_row_order_callback <- function(esbl_named_vector, esbl_priority) {
  force(esbl_named_vector)
  force(esbl_priority)
  
  function(hc, mat) {
    rn <- rownames(mat)
    
    esbl <- esbl_named_vector[rn]
    esbl_rank <- unname(esbl_priority[as.character(esbl)])
    
    # fallback if something is missing
    esbl_rank[is.na(esbl_rank)] <- max(esbl_priority, na.rm = TRUE) + 1
    
    # low mean resistance => higher up
    mm_mean <- rowMeans(mat, na.rm = TRUE)
    
    # ESBL should dominate, mm only breaks ties
    weights <- esbl_rank * 1000 + mm_mean
    
    reorder_hclust_by_weights(hc, weights = weights)
  }
}

make_annotation_legend <- function(title, labels, colors) {
  
  df <- data.frame(
    label = factor(labels, levels = rev(labels)),
    y = seq_along(labels)
  )
  
  ggplot(df, aes(x = 1, y = y, fill = label)) +
    geom_tile(width = 0.7, height = 0.7) +
    geom_text(
      aes(label = as.character(label)),
      x = 1.6,
      hjust = 0,
      size = 3
    ) +
    scale_fill_manual(
      values = setNames(colors, labels),
      guide = "none"
    ) +
    scale_y_reverse() +
    coord_cartesian(clip = "off") +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.margin = margin(2, 20, 2, 2)
    )
}


make_pheatmap <- function(mat,
                          annotation_row = NULL,
                          annotation_colors = NULL,
                          color = gradient_red(),
                          show_rownames = FALSE,
                          show_colnames = TRUE,
                          cluster_rows = TRUE,
                          cluster_cols = TRUE,
                          breaks = NULL,
                          legend = FALSE,
                          legend_breaks = NULL,
                          legend_labels = NULL,
                          display_numbers = NULL,
                          fontsize_row = 4,
                          fontsize_col = 7,
                          fontsize_number = 6,
                          fontsize = 7,
                          clustering_callback = NULL,
                          annotation_names_row = TRUE,
                          annotation_legend = TRUE) {
  

  args <- list(
    mat = mat,
    clustering_method = "ward.D2",
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    color = color,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    breaks = breaks,
    legend = legend,
    legend_labels = legend_labels,
    legend_breaks = legend_breaks,
    border_color = NA,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    fontsize_number = fontsize_number,
    fontsize = fontsize,
    main = "",
    annotation_names_row = annotation_names_row,
    annotation_legend = annotation_legend,
    cellwidth = 7
  )
  
  if (!is.null(display_numbers)) {
    args$display_numbers <- display_numbers
  }
  
  if (!is.null(clustering_callback)) {
    args$clustering_callback <- clustering_callback
  }
  
  do.call(pheatmap::pheatmap, args)
}


# ============================================================
# Choosing k
# ============================================================

estimate_optimal_number_of_clusters <- function(data_matrix, max_k = 10, seed = 123) {
  set.seed(seed)
  
  pv = pvclust(
    t(data_matrix),
    method.hclust = "ward.D2",
    method.dist = "euclidean",
    nboot = 10)
  plot(pv)    
  
  library(pvclust)
  list(
    silhouette = fviz_nbclust(data_matrix, hcut, method = "silhouette", hc_method = "ward.D2") +
      labs(subtitle = "Silhouette method"),
    gap = fviz_nbclust(data_matrix, hcut, method = "gap_stat", nboot = 50, hc_method = "ward.D2") +
      labs(subtitle = "Gap statistic method"),
    wss = fviz_nbclust(data_matrix, hcut, method = "wss", hc_method = "ward.D2") +
      labs(subtitle = "Elbow method"),
    pv
    
  )
}

fetchPredictionErrors <- function()
{

  metrics_long_ab <- derivedMetricsFrame(countFrameWide = readCountSampleFrameWide(),vars = SAMPLE_GROUPS_VAR)
  
  errors <- metrics_long_ab %>% 
    filter_and_drop(mode,MODE_A) %>% 
    filter_and_drop(cpmode,"normal") %>%
    filter_and_drop(significanceLevel,"STD") %>%
    filter_and_drop(noinputab,6) %>%
    filter(metric %in% c("ME","VME"))
  
  error_matrix <- errors %>%
    pivot_wider(
      names_from = metric,
      values_from = value
    ) %>%
    mutate(
      n_na = rowSums(is.na(across(c(ME, VME))))
    ) %>%
    {
      if (any(.$n_na != 1)) {
        stop("Sanity check failed: each sample-antibiotic pair must have exactly one NA among ME and VME")
      }
      .
    } %>%
    mutate(
      combined = if_else(!is.na(ME), ME, -VME)
    ) %>%
    select(sample, antibiotic, combined) %>%
    mutate(
      antibiotic = factor(antibiotic, levels = ANTIBIOTICS)
    ) %>%
    complete(sample, antibiotic = factor(ANTIBIOTICS, levels = ANTIBIOTICS)) %>%
    pivot_wider(
      names_from = antibiotic,
      values_from = combined
    ) %>%
    arrange(sample) %>%
    as.data.frame()
  
  rownames(error_matrix) <- error_matrix$sample
  error_matrix$sample <- NULL
  
  error_matrix <- error_matrix[, ANTIBIOTICS]
  
  stopifnot(all(colnames(error_matrix) == ANTIBIOTICS))
  error_matrix
}

# ============================================================
# Main workflow
# ============================================================

make_cluster_and_heatmaps <- function(data_dir = processedRootRcommon,
                                      primary_table = "MM",
                                      k = 6,
                                      add_esbl = TRUE) {
  dat <- load_cluster_data(data_dir)
  tables <- fetch_tables(dat)
  
  stopifnot(primary_table %in% names(tables))
  
  primary <- tables[[primary_table]]
  caa <- clusters_and_annotation(primary, k = k)

  if(k==4){
    cluster_order <- c(3, 4, 1, 2)
    caa <- renumber_clusters_by_order(caa, cluster_order)
  }

  if(k==5){
    cluster_order <- c(3, 4, 1, 2, 5)
    caa <- renumber_clusters_by_order(caa, cluster_order)
  }
  if(k==6){
    cluster_order <- c(4,5,1,3,2, 6)
    caa <- renumber_clusters_by_order(caa, cluster_order)
  }
  
  
  
  cluster_plot <- make_cluster_plot(caa)
  #cluster_plot <- make_cluster_plot_ggplot(caa,label_esbl = TRUE)
  cluster_legend <- make_annotation_legend(
    title = PHENOTYPE_CLUSTER_LABEL,
    labels = as.character(1:k),
    colors = unname(cluster_palette(k))
  )
  
  #Save caa with only clusters 
  caaBasic <- caa
  
  #
  if (add_esbl) {
    caa <- decorate_with_esbl(caa, dat)
  }
  esbl_vec <- caa$annotation[[PHENOTYPE_ESBL_LABEL]]
  names(esbl_vec) <- rownames(caa$annotation)
  
  row_callback <- make_row_order_callback(
    esbl_named_vector = esbl_vec,
    esbl_priority = ESBL_PRIORITY
  )
  
  row_names_ordered <- order_rows_by_cluster_esbl_resistance(
    caa = caa,
    mat = primary
  )

  # # rename cluster annotation display name
  # caaBasic$annotation <- caaBasic$annotation |>
  #   dplyr::rename(`Phenotype cluster` = Cluster)
  # 
  # names(caaBasic$annotation_colors)[1] <- "Phenotype cluster"
  
  
  
  
  
  primary_heatmap <- make_pheatmap(
    mat = primary[row_names_ordered, , drop = FALSE],
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = rev(gradient_red()),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE,
    clustering_callback = row_callback
  )
  
  row_order <- row_names_ordered

  # SIR heatmap with ATU stars
  sir_tbl <- tables$SIR
  atu_tbl <- read_atu_antibiotics_model(dat) %>% drop_sample_column()
  error_tbl <- tables$ERRORS

  
  display_numbers <- matrix(
    ifelse(as.matrix(atu_tbl), "\u2217", ""),
    nrow = nrow(atu_tbl),
    dimnames = dimnames(atu_tbl)
  )
  
  sir_heatmap <- make_pheatmap(
    mat = sir_tbl[row_order, , drop = FALSE],
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = gradient_green_yellow_red(),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = FALSE,
    display_numbers = display_numbers[row_order, , drop = FALSE]
  )


  
  
  combined_1 <- (
    (cluster_plot + panel_box_theme) +
      (wrap_pheatmap(primary_heatmap) & panel_box_theme) + 
      (wrap_pheatmap(sir_heatmap) & panel_box_theme)
  ) +
    plot_layout(ncol = 3) +
    plot_annotation(tag_levels = "A")
  
  genotypeGroupTable <- getGenotypeGroupTable()
  caa <- caaBasic
  
  caa <- aggregate_presence_by_color(
    caa = caa,
    decorate_fun = decorate_with_functional_group_presence,
    genotypeGroupTable = genotypeGroupTable,
    label_map = c("darkred" = "ESBL carba", "orange" = "AMPC", "yellow" = "ESBL classic", "grey70" = "Non-ESBL")
  )
  print(colnames(caa$annotation))
  
  beta_cols <- intersect(BETA_LACTAMS, colnames(primary))
  
  primary_heatmap_betalactams <- make_pheatmap(
    mat = primary[row_order, beta_cols, drop = FALSE], 
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = rev(gradient_red()),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE
  )

  caa <- caaBasic
  caa <- aggregate_presence_by_color(
    caa = caa,
    decorate_fun = decorate_with_aminoglycoside_family_presence,
    genotypeGroupTable = genotypeGroupTable,
    label_map = c("red" = "AAC", "orange" = "APH", "yellow" = "AAD/ANT")
  )

  ag_cols <- c("GEN","TOB")
  
  print(colnames(caa$annotation))
  primary_heatmap_aminoglycosides <- make_pheatmap(
    mat = primary[row_order, ag_cols, drop = FALSE],
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = rev(gradient_red()),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE
  )

  caa <- caaBasic

  caa <- aggregate_presence_by_color(
    caa = caa,
    decorate_fun = decorate_with_quinolone_family_presence,
    genotypeGroupTable = genotypeGroupTable,
    label_map = c(
      "orange" = "QRDR",
      "yellow" = "PMQR"
    ),
    drop_original = TRUE
  )

  print(colnames(caa$annotation))
  
  quinolones_cols <- c("CIP" ,"OFX" ,"LVX", "MFX")
  
  primary_heatmap_quinolones <- make_pheatmap(
    mat = primary[row_order, quinolones_cols, drop = FALSE], 
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = rev(gradient_red()),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE,
  )
  
  ph_list <- list(primary_heatmap_quinolones,
                  primary_heatmap_aminoglycosides,
                  primary_heatmap_betalactams,
                  sir_heatmap)

  
  combined_2A <- (
    (wrap_pheatmap(ph_list[[1]]) & panel_box_theme) +
      (wrap_pheatmap(ph_list[[2]]) & panel_box_theme) +
      (wrap_pheatmap(ph_list[[3]]) & panel_box_theme) +
      (wrap_pheatmap(ph_list[[4]]) & panel_box_theme)
  ) +
    plot_layout(ncol = 4) +
    plot_annotation(tag_levels = "A")
     
  
  
  
  caa <- caaBasic
  caa <- decorate_with_quinolone_family_presence(genotypeGroupTable = genotypeGroupTable,caa = caa) 
  sir_heatmap_quinolones_family <- make_pheatmap(
    mat = sir_tbl[row_order, quinolones_cols, drop = FALSE], 
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = gradient_green_yellow_red(),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = FALSE,
    display_numbers = display_numbers[row_order, quinolones_cols, drop = FALSE],
  )
  
  caa <- caaBasic
  caa <- decorate_with_aminoglycoside_family_presence(genotypeGroupTable = genotypeGroupTable,caa = caa) 
  sir_heatmap_aminoglycosides_family <- make_pheatmap(
    mat = sir_tbl[row_order, ag_cols, drop = FALSE], 
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = gradient_green_yellow_red(),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = FALSE,
    display_numbers = display_numbers[row_order, ag_cols, drop = FALSE]
  )
  
  caa <- caaBasic
  caa <- decorate_with_esbl(caa,dat) 
  caa <- decorate_with_functional_group_presence(genotypeGroupTable = genotypeGroupTable,caa = caa) 
  sir_heatmap_betalactams_bj <- make_pheatmap(
    mat = sir_tbl[row_order, beta_cols, drop = FALSE], 
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = gradient_green_yellow_red(),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = FALSE,
    display_numbers = display_numbers[row_order, beta_cols, drop = FALSE]
  )
  
  
  
  geneMappingTable <- getGeneMappingTable()
  
  caa <- caaBasic
  caa <- decorate_with_esbl(caa,dat)
  caa <- decorate_with_beta_gene_presence(caa = caa,
                                   genotypeGroupTable = genotypeGroupTable,
                                   mapping_table = geneMappingTable)
  
  primary_heatmap_betalactam_gene <- make_pheatmap(
    mat = primary[row_order, beta_cols, drop = FALSE], 
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = gradient_green_yellow_red(),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE,
    annotation_legend = FALSE
  )
  
  caa <- caaBasic
  caa <- decorate_with_aminoglycoside_gene_presence(caa = caa,
                                          genotypeGroupTable = genotypeGroupTable,
                                          geneMappingTable = geneMappingTable)

  primary_heatmap_aminoglycoside_gene <- make_pheatmap(
    mat = primary[row_order, ag_cols, drop = FALSE],
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = gradient_green_yellow_red(),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE,
    annotation_legend = FALSE
  )
  caa <- caaBasic
  caa <- decorate_with_quinolone_gene_presence(caa = caa,
                                                    genotypeGroupTable = genotypeGroupTable,
                                               geneMappingTable = geneMappingTable)

  # caa_legend <- keep_annotation_legend(
  #   caa,
  #   keep = PHENOTYPE_CLUSTER_LABEL
  # )
  # 
  primary_heatmap_quinolone_gene <- make_pheatmap(
    mat = primary[row_order, quinolones_cols, drop = FALSE],
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = gradient_green_yellow_red(),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE,
    annotation_legend = FALSE
  )
  

  ph_list <- list(primary_heatmap_quinolone_gene,
                  primary_heatmap_aminoglycoside_gene,
                  primary_heatmap_betalactam_gene,
                  sir_heatmap)
  combined_2B <- (
    (wrap_pheatmap(ph_list[[1]]) & panel_box_theme) +
      (wrap_pheatmap(ph_list[[2]]) & panel_box_theme) +
      (wrap_pheatmap(ph_list[[3]]) & panel_box_theme)+
      (wrap_pheatmap(ph_list[[4]]) & panel_box_theme)
  ) +
    plot_layout(ncol = 4) +
    plot_annotation(tag_levels = "A")


  
  
  
  
  caa <- caaBasic
  
  
  cols <- gradient_red_white_orange(199)
  
  error_heatmap <- make_pheatmap(
    mat = error_tbl[row_order, colnames(primary), drop = FALSE],
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = cols,
    breaks = seq(-1, 1, length.out = length(cols) + 1),
    legend_breaks = c(-1, 0, 1),
    legend_labels = c("VME", "", "ME"),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE,
    annotation_legend = FALSE
  )
  

  combined_3 <- (
      (wrap_pheatmap(error_heatmap) & panel_box_theme) + 
      (wrap_pheatmap(sir_heatmap) & panel_box_theme)
  ) +
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = "A")
    
  
  
  

  
  error_colors <- cols
  error_breaks <- seq(-1, 1, length.out = length(error_colors) + 1)

  
  
  row_order_error_betalactams <- order_rows_by_group_error(
    error_tbl = error_tbl,
    cols = beta_cols
  )
  
  row_order_error_aminoglycosides <- order_rows_by_group_error(
    error_tbl = error_tbl,
    cols = ag_cols
  )
  
  row_order_error_quinolones <- order_rows_by_group_error(
    error_tbl = error_tbl,
    cols = quinolones_cols
  )
  
  row_order_error <- order_rows_by_error(
    error_tbl = error_tbl
  )
  
  row_order_error_betalactams <- row_order_error
  row_order_error_aminoglycosides <- row_order_error
  row_order_error_quinolones <- row_order_error
  
  caa <- caaBasic
  caa <- decorate_with_mean_absolute_error(caa = caa,
                                           error_tbl = error_tbl)
  
  
  
  caa <- aggregate_presence_by_color(
    caa = caa,
    decorate_fun = decorate_with_functional_group_presence,
    genotypeGroupTable = genotypeGroupTable,
    label_map = c("darkred" = "ESBL carba", "orange" = "AMPC", "yellow" = "ESBL classic", "grey70" = "Non-ESBL")
  )

    
  error_heatmap_betalactams <- make_pheatmap(
    mat = error_tbl[row_order_error_betalactams, beta_cols, drop = FALSE],
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = error_colors,
    breaks = error_breaks,
    legend_breaks = c(-1, 0, 1),
    legend_labels = c("VME", "", "ME"),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE
  )
  
  caa <- caaBasic
  caa <- decorate_with_mean_absolute_error(caa = caa,
                                           error_tbl = error_tbl)
  caa <- aggregate_presence_by_color(
    caa = caa,
    decorate_fun = decorate_with_aminoglycoside_family_presence,
    genotypeGroupTable = genotypeGroupTable,
    label_map = c("red" = "AAC", "orange" = "APH", "yellow" = "AAD/ANT")
  )
  error_heatmap_aminoglycosides <- make_pheatmap(
    mat = error_tbl[row_order_error_aminoglycosides, ag_cols, drop = FALSE],
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = error_colors,
    breaks = error_breaks,
    legend_breaks = c(-1, 0, 1),
    legend_labels = c("VME", "", "ME"),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE
  )
  caa <- caaBasic
  caa <- decorate_with_mean_absolute_error(caa = caa,
                                           error_tbl = error_tbl)
  
  
  caa <- aggregate_presence_by_color(
    caa = caa,
    decorate_fun = decorate_with_quinolone_family_presence,
    genotypeGroupTable = genotypeGroupTable,
    label_map = c(
      "orange" = "QRDR",
      "yellow" = "PMQR"
    ),
    drop_original = TRUE
  )
  
  error_heatmap_quinolones <- make_pheatmap(
    mat = error_tbl[row_order_error_quinolones, quinolones_cols, drop = FALSE],
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = error_colors,
    breaks = error_breaks,
    legend_breaks = c(-1, 0, 1),
    legend_labels = c("VME", "", "ME"),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE
  )
  
  caa <- caaBasic
  caa <- decorate_with_mean_absolute_error(caa = caa,
                                           error_tbl = error_tbl)
  sir_heatmap_ordered_by_error <- make_pheatmap(
      mat = sir_tbl[row_order_error, , drop = FALSE],
      annotation_row = caa$annotation,
      annotation_colors = caa$annotation_colors,
      color = gradient_green_yellow_red(),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      legend = FALSE,
      display_numbers = display_numbers[row_order_error, , drop = FALSE]
    )
  
  
  
  ph_list <- list(error_heatmap_quinolones,
                  error_heatmap_aminoglycosides,
                  error_heatmap_betalactams,
                  sir_heatmap_ordered_by_error)
  
  
  combined_4 <- (
    (wrap_pheatmap(ph_list[[1]]) & panel_box_theme) +
      (wrap_pheatmap(ph_list[[2]]) & panel_box_theme) +
      (wrap_pheatmap(ph_list[[3]]) & panel_box_theme) +
      (wrap_pheatmap(ph_list[[4]]) & panel_box_theme)
  ) +
    plot_layout(ncol = 4) +
    plot_annotation(tag_levels = "A")
  
  
  
  caa <- caaBasic
  caa <- decorate_with_mean_absolute_error(caa = caa,
                                           error_tbl = error_tbl)
  
  cols <- gradient_red_white_orange(199)
  caa <- decorate_with_beta_gene_presence(caa = caa,
                                          genotypeGroupTable = genotypeGroupTable,
                                          mapping_table = geneMappingTable)

  caa <- decorate_with_aminoglycoside_gene_presence(caa = caa,
                                                    genotypeGroupTable = genotypeGroupTable,
                                                    geneMappingTable = geneMappingTable)
  
  
  caa <- decorate_with_quinolone_gene_presence(caa = caa,
                                               genotypeGroupTable = genotypeGroupTable,
                                               geneMappingTable = geneMappingTable)
  
  
  error_heatmap_allgenes <- make_pheatmap(
    mat = error_tbl[row_order_error, colnames(primary), drop = FALSE],
    annotation_row = caa$annotation,
    annotation_colors = caa$annotation_colors,
    color = cols,
    breaks = seq(-1, 1, length.out = length(cols) + 1),
    legend_breaks = c(-1, 0, 1),
    legend_labels = c("VME", "", "ME"),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    legend = TRUE,
    annotation_legend = FALSE
  )
  
  
  ph_list <- list(error_heatmap_allgenes,
                  sir_heatmap_ordered_by_error)
  
  combined_5 <- (
    (wrap_pheatmap(ph_list[[1]]) & panel_box_theme) +
      (wrap_pheatmap(ph_list[[2]]) & panel_box_theme) 
  ) +
    plot_layout(ncol = 2) +
    plot_layout(widths = c(3, 1)) +
    plot_annotation(tag_levels = "A")
  
  
  
  
  list(
    data = dat,
    tables = tables,
    cluster = caa,
    k = k,
    row_order = row_order,
    cluster_plot = cluster_plot,
    primary_heatmap = primary_heatmap,
    sir_heatmap = sir_heatmap,
    combined_phenotype_cluster_primary_sir = combined_1,
    combined_genotype_aggregated_sir = combined_2A,
    combined_genotype_sir = combined_2B,
    combined_predictionerror_sir = combined_3,
    combined_predictionerror_genotype_aggregated_sir = combined_4,
    combined_predictionerror_genotype_sir = combined_5
  )
}




ESTIMATE <- function()
{
  dat <- load_cluster_data(processedRootRcommon)
  tables <- fetch_tables(dat)

#  estimate_optimal_number_of_clusters(tables$MM)
#  estimate_optimal_number_of_clusters(tables$RESCALEDMM)
  estimate_optimal_number_of_clusters(tables$MM)
}
  

errorSampleOrder <- function()
{
  data <- load_cluster_data(processedRootRcommon)
  tables <- fetch_tables(data)
  errors <- tables$ERRORS
  absErrors <- abs(errors)
  zzz <- estimate_optimal_number_of_clusters(errors)
  caa <- clusters_and_annotation(tables$MM,6)
  colnames(caa$annotation)
  colnames(caa$annotation_colors)
  
}



ALL <- function()
{
  # create clustering using hclust
  # rotating dendrogram to get clusters
  #
  result <- make_cluster_and_heatmaps(
    data_dir = processedRootRcommon,
    primary_table = "MM",
    k = 6,
    add_esbl = TRUE
  )
  
  
  result$tables$ERRORS[result$row_order,]
  
  clusterInformation <- tibble(
    sample = result$data$mm_full_rescaled$sample,
    roworder = result$row_order,
    cluster = result$cluster$annotation$Cluster,
    ESBL = result$cluster$annotation[[PHENOTYPE_ESBL_LABEL]]
  )
  clusterInformation %>% 
    write.csv2(
      file.path(processedRootRcluster,
                sprintf("phenotype_cluster%d.csv",result$k))
                ,row.names = FALSE)

  result$tables$ERRORS[result$row_order,] %>% write.csv2(
    file.path(processedRootRcluster,
              sprintf("errorrate_cluster%d.csv",result$k))
    ,row.names = TRUE)
  
  
  ggplot2::ggsave(
    filename = file.path(processedRootRcluster,
                         sprintf("cluster_phenotype%d.png",result$k)),
    plot = result$combined_phenotype_cluster_primary_sir,
    width = 13,
    height = 8,
    dpi = 300
  )
  print(result$combined_phenotype_cluster_primary_sir)
  
  ggplot2::ggsave(
    filename = file.path(processedRootRcluster,
                         sprintf("genotype_aggregated_sir%d.png",result$k)),
    plot = result$combined_genotype_aggregated_sir,
    width = 13,
    height = 8,
    dpi = 300
  )
  
  ggplot2::ggsave(
    filename = file.path(processedRootRcluster,
                         sprintf("genotype_sir%d.png",result$k)),
    plot = result$combined_genotype_sir,
    width = 15,
    height = 8,
    dpi = 300
  )
  
  ggplot2::ggsave(
    filename = file.path(processedRootRcluster,
                         sprintf("predictionerror_sir%d.png",result$k)),
    plot = result$combined_predictionerror_sir,
    width = 8,
    height = 8,
    dpi = 300
  ) 
  
  ggplot2::ggsave(
    filename = file.path(processedRootRcluster,
                         sprintf("predictionerror_genotype_aggregated_sir%d.png",result$k)),
    plot = result$combined_predictionerror_genotype_aggregated_sir,
    width = 13,
    height = 8,
    dpi = 300
  ) 
  ggplot2::ggsave(
    filename = file.path(processedRootRcluster,
                         sprintf("predictionerror_genotype_sir%d.png",result$k)),
    plot = result$combined_predictionerror_genotype_sir,
    width = 13,
    height = 8,
    dpi = 300
  ) 
  
}



