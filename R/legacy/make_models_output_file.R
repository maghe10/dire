suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
})

make_model_output_files <- function(input_dir,
                                    output_dir,
                                    k,
                                    file_prefix     = "sirAntibioticsModelWordsJuan_",
                                    sample_range    = 1:99,
                                    true_sample_ids = sample_range,
                                    template_file   = NULL) {
  ## ---- Validate sample mapping -----------------------------------------
  if (length(sample_range) != length(true_sample_ids)) {
    stop("sample_range and true_sample_ids must be the same length")
  }
  sample_map <- setNames(true_sample_ids, sample_range)
  
  ## ---- Paths ------------------------------------------------------------
  input_dir  <- normalizePath(input_dir, mustWork = TRUE)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  
  ## ---- Logical antibiotic order ----------------------------------------
  antibiotic_order <- c(
    "AMP","AMC","PIP","TZP",
    "CAZ","CRO","CTX","FEP",
    "CIP","OFX","LVX","MFX",
    "GEN","TOB"
  )
  
  ### Token sorting helper (logical order)
  order_tokens <- function(tokens) {
    tokens <- tokens[!is.na(tokens)]
    if (!length(tokens)) return(character(0))
    parts <- str_split(tokens, "_", n=2)
    abs   <- vapply(parts, `[`, character(1), 1)
    idx   <- match(abs, antibiotic_order)
    # Unknown antibiotics → place last
    idx[is.na(idx)] <- max(idx, na.rm=TRUE) + seq_len(sum(is.na(idx)))
    tokens[order(idx)]
  }
  
  ### Predictor parsing with logical AB ordering
  parse_predictors <- function(x) {
    if (is.na(x) || x == "") return(character(0))
    toks <- str_split(x, "\\s+", simplify = TRUE)
    if (length(toks) <= 1) return(character(0))
    toks <- toks[-1]
    ab   <- unique(str_replace(toks, "_.*$", ""))
    idx  <- match(ab, antibiotic_order)
    idx[is.na(idx)] <- length(antibiotic_order) + seq_len(sum(is.na(idx)))
    ab[order(idx)]
  }
  
  ### cp → token
  cp_to_token <- function(ab, cp_val) {
    if (is.na(cp_val) || cp_val %in% c("", "NA")) return(NA)
    cp_val <- toupper(cp_val)
    if (cp_val == "S/R") return(paste0(ab,"_SR"))
    if (cp_val %in% c("S","R")) return(paste0(ab,"_",cp_val))
    NA
  }
  
  ## ---- Optional template -------------------------------------------------
  template <- NULL
  if (!is.null(template_file)) {
    template_file <- normalizePath(template_file, mustWork = TRUE)
    template <- read.csv2(template_file, stringsAsFactors=FALSE)
  }
  
  ## ---- Read all sample files --------------------------------------------
  df_list <- list()
  
  for (sample_idx in sample_range) {
    in_file <- file.path(input_dir, sprintf("%s%d_%d_cp.csv", file_prefix, k, sample_idx))
    if (!file.exists(in_file)) {
      message("Skipping missing file: ", in_file)
      next
    }
    
    tmp <- read.csv(in_file, stringsAsFactors=FALSE)
    if ("index" %in% names(tmp)) tmp <- tmp %>% filter(index != "id")
    
    true_sample <- sample_map[as.character(sample_idx)]
    
    tmp <- tmp %>%
      mutate(
        sample                = as.integer(true_sample),
        AST_true              = as.character(AST_true),
        AST_pred              = as.character(AST_pred),
        cp_85                 = as.character(cp_85),
        cp_90                 = as.character(cp_90),
        cp_95                 = as.character(cp_95),
        cp_975                = as.character(cp_975),
        antibiotic            = as.character(antibiotic),
        Antibiotic_predictors = as.character(Antibiotic_predictors)
      )
    
    df_list[[length(df_list)+1]] <- tmp
  }
  
  if (!length(df_list)) stop("No usable files.")
  
  df <- bind_rows(df_list)
  
  ## ---- Compute combo_name -----------------------------------------------
  df <- df %>%
    rowwise() %>%
    mutate(
      predictors_vec = list(parse_predictors(Antibiotic_predictors)),
      combo_name     = paste(predictors_vec, collapse="_")
    ) %>%
    ungroup()
  
  ## ---- Template-controlled or free ordering -----------------------------
  if (!is.null(template)) {
    sample_ids   <- template$sample
    combo_levels <- setdiff(names(template),"sample")
  } else {
    sample_ids   <- sort(unique(df$sample))
    combo_levels <- sort(unique(df$combo_name))
  }
  
  ## ---- Initialize WIDE tables ------------------------------------------
  make_empty <- function() {
    out <- data.frame(
      sample = sample_ids,
      matrix("<empty>", length(sample_ids), length(combo_levels)),
      stringsAsFactors=FALSE
    )
    colnames(out)[-1] <- combo_levels
    out
  }
  
  answer_wide   <- make_empty()
  preds_wide    <- make_empty()
  conf01_wide   <- make_empty()
  conf005_wide  <- make_empty()
  conf0025_wide <- make_empty()
  
  ## ---- Build LONG dump (NEW) -------------------------------------------
  long_rows <- list()   ### NEW
  
  ## ---- Fill both WIDE and LONG -----------------------------------------
  df_grouped <- df %>% group_by(sample, combo_name)
  
  groups <- split(df_grouped, interaction(df_grouped$sample, df_grouped$combo_name))
  
  for (g in groups) {
    smp   <- unique(g$sample)
    combo <- unique(g$combo_name)
    
    if (!(combo %in% combo_levels)) next
    
    row <- match(smp, sample_ids)
    col <- match(combo, combo_levels)
    
    # All ANTIBIOTICS in this group
    ans_tokens  <- sprintf("%s_%s", g$antibiotic, g$AST_true)
    pred_tokens <- sprintf("%s_%s", g$antibiotic, g$AST_pred)
    
    conf01_tokens   <- mapply(cp_to_token, g$antibiotic, g$cp_90)
    conf005_tokens  <- mapply(cp_to_token, g$antibiotic, g$cp_95)
    conf0025_tokens <- mapply(cp_to_token, g$antibiotic, g$cp_975)
    
    # Sort tokens logically
    answer_wide[row, col+1]   <- paste(order_tokens(ans_tokens),  collapse=" ")
    preds_wide[row, col+1]    <- paste(order_tokens(pred_tokens), collapse=" ")
    conf01_wide[row, col+1]   <- paste(order_tokens(conf01_tokens), collapse=" ")
    conf005_wide[row, col+1]  <- paste(order_tokens(conf005_tokens), collapse=" ")
    conf0025_wide[row, col+1] <- paste(order_tokens(conf0025_tokens), collapse=" ")
    
    ### Add LONG rows (NEW)
    long_rows[[length(long_rows)+1]] <- data.frame(
      sample   = smp,
      combo    = combo,
      antibiotic = g$antibiotic,
      AST_true   = g$AST_true,
      AST_pred   = g$AST_pred,
      cp_01      = g$cp_90,
      cp_005     = g$cp_95,
      cp_0025    = g$cp_975,
      stringsAsFactors = FALSE
    )
  }
  
  ### Combine LONG dump (NEW)
  long_df <- bind_rows(long_rows)
  
  ## ---- Write output files -----------------------------------------------
  write.csv2(answer_wide,   file.path(output_dir, sprintf("modelOutput_answer_%d.csv", k)), row.names=FALSE)
  write.csv2(preds_wide,    file.path(output_dir, sprintf("modelOutput_preds_%d.csv",  k)), row.names=FALSE)
  write.csv2(conf01_wide,   file.path(output_dir, sprintf("modelOutput_confpreds_%d_01.csv",   k)), row.names=FALSE)
  write.csv2(conf005_wide,  file.path(output_dir, sprintf("modelOutput_confpreds_%d_005.csv",  k)), row.names=FALSE)
  write.csv2(conf0025_wide, file.path(output_dir, sprintf("modelOutput_confpreds_%d_0025.csv", k)), row.names=FALSE)
  
  ### Write LONG file (NEW)
  colnames(long_df) <- c("sample","word",	"antibiotic" ,"answer","pred","confpred_01",	"confpred_005"	,"confpred_0025")
  write.csv2(long_df, file.path(output_dir, sprintf("modelOutput_long_%d.csv", k)), row.names=FALSE)
  
  invisible(list(
    answer = answer_wide,
    preds  = preds_wide,
    conf01 = conf01_wide,
    conf005= conf005_wide,
    conf0025 = conf0025_wide,
    long   = long_df   ### NEW
  ))
}
