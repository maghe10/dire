# ============================================================================
# Word / input-antibiotic analysis versus MMC
# ============================================================================
# Purpose
#   1) Analyze whether presence of a given input antibiotic is associated with
#      higher MMC, stratified by number of input antibiotics (4, 6, 8).
#   2) Find the best-performing input antibiotic words for 4, 6, and 8 inputs.
#   3) Evaluate whether adding a specific antibiotic tends to improve MMC.
#   4) Create illustrations and export tabular summaries.
#
# Data assumptions
#   - readCountWordFrameLong() returns a long tibble with columns:
#       cpmode, mode, word, noinputab, antibiotic, significanceLevel, metric, count
#   - The metric column contains confusion-count metrics named:
#       correctR, falseR, correctS, falseS
#     where
#       TP = correctR, FP = falseR, TN = correctS, FN = falseS
#
# Analysis scope fixed by request
#   - cpmode == "normal"
#   - mode == "Mode-A"
#   - significanceLevel is NA
#   - noinputab in c(4, 6, 8)
#
# Notes
#   - MMC is interpreted here as Matthews correlation coefficient (MCC).
#   - MCC is computed from confusion counts after aggregating counts across all
#     predicted antibiotics for each input word.
#   - To avoid NA, MCC is set to 0 when the denominator is 0.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
  library(readr)
  library(tibble)
})

source('model/modelcommon.R')
source('manuscript/manuscriptcommon.R')

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
TARGET_NOINPUTAB <- c(4, 6, 8)
TARGET_CPMODE <- "normal"
TARGET_MODE <- "Mode-A"
OUT_DIR <- file.path(processedRootRassembly, "model")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
# safe_mcc <- function(tp, fp, tn, fn) {
#   num <- tp * tn - fp * fn
#   den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
#   ifelse(is.na(den) | den == 0, 0, num / den)
# }
# MCC defined for all confusion matrices using the convention discussed
safe_mcc <- function(tp, tn, fp, fn) {
  tp <- as.numeric(tp)
  tn <- as.numeric(tn)
  fp <- as.numeric(fp)
  fn <- as.numeric(fn)
  
  n <- length(tp)
  out <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    TPi <- tp[i]
    TNi <- tn[i]
    FPi <- fp[i]
    FNi <- fn[i]
    
    vals <- c(TPi, TNi, FPi, FNi)
    
    if (any(is.na(vals))) {
      out[i] <- NA_real_
      next
    }
    
    num <- TPi * TNi - FPi * FNi
    den <- sqrt((TPi + FPi) * (TPi + FNi) * (TNi + FPi) * (TNi + FNi))
    
    if (!is.na(den) && den > 0) {
      out[i] <- num / den
      next
    }
    
    nonzero_n <- sum(vals > 0)
    
    if (nonzero_n == 1) {
      if (TPi > 0 || TNi > 0) {
        out[i] <- 1
      } else if (FPi > 0 || FNi > 0) {
        out[i] <- -1
      } else {
        out[i] <- 0
      }
      next
    }
    
    out[i] <- 0
  }
  
  out
}


split_word <- function(word) {
  strsplit(word, "_", fixed = TRUE)[[1]]
}

word_contains_ab <- function(word, ab) {
  ab %in% split_word(word)
}

compute_word_level_mmc <- function(word_df) {
  required_metrics <- c("correctR", "falseR", "correctS", "falseS")

  wide <- word_df %>%
    filter(metric %in% required_metrics) %>%
    group_by(word, noinputab, metric) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = metric,
      values_from = count,
      values_fill = 0
    )

  missing_metrics <- setdiff(required_metrics, colnames(wide))
  if (length(missing_metrics) > 0) {
    stop(
      "Missing required metric columns after pivot: ",
      paste(missing_metrics, collapse = ", ")
    )
  }

  wide %>%
    mutate(
      TP = correctR,
      FP = falseR,
      TN = correctS,
      FN = falseS,
      MMC = safe_mcc(tp=TP, fp=FP, tn=TN, fn=FN)
    ) %>%
    select(word, noinputab, TP, FP, TN, FN, MMC)
}

make_presence_table <- function(word_metrics) {
  all_input_abs <- sort(unique(unlist(lapply(word_metrics$word, split_word))))

  tidyr::expand_grid(
    word = unique(word_metrics$word),
    input_ab = all_input_abs
  ) %>%
    mutate(present = map2_lgl(word, input_ab, word_contains_ab)) %>%
    left_join(word_metrics, by = "word") %>%
    relocate(noinputab, input_ab, present, .after = word)
}

compute_presence_effects <- function(presence_df) {
  presence_df %>%
    group_by(noinputab, input_ab) %>%
    summarise(
      n_present = sum(present),
      n_absent = sum(!present),
      mean_MMC_present = mean(MMC[present], na.rm = TRUE),
      mean_MMC_absent = mean(MMC[!present], na.rm = TRUE),
      median_MMC_present = median(MMC[present], na.rm = TRUE),
      median_MMC_absent = median(MMC[!present], na.rm = TRUE),
      delta_mean_MMC = mean_MMC_present - mean_MMC_absent,
      p_value = tryCatch(
        wilcox.test(MMC[present], MMC[!present], exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    group_by(noinputab) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      rank_delta = rank(-delta_mean_MMC, ties.method = "min")
    ) %>%
    ungroup() %>%
    arrange(noinputab, desc(delta_mean_MMC), input_ab)
}

compute_best_words <- function(word_metrics, top_n = 15) {
  word_metrics %>%
    group_by(noinputab) %>%
    arrange(desc(MMC), desc(TP + TN), word, .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    filter(rank <= top_n)
}

compute_addition_effects <- function(word_metrics, target_sizes = c(4, 6, 8)) {
  words_by_size <- split(word_metrics, word_metrics$noinputab)
  out <- list()

  for (k in target_sizes) {
    parent_k <- k - 1
    if (!as.character(parent_k) %in% names(words_by_size) ||
        !as.character(k) %in% names(words_by_size)) {
      next
    }

    parents <- words_by_size[[as.character(parent_k)]] %>%
      transmute(
        parent_word = word,
        parent_MMC = MMC,
        parent_set = map(word, ~ sort(split_word(.x)))
      )

    children <- words_by_size[[as.character(k)]] %>%
      transmute(
        child_word = word,
        noinputab = k,
        child_MMC = MMC,
        child_set = map(word, ~ sort(split_word(.x)))
      )

    child_parent_rows <- purrr::pmap_dfr(
      children,
      function(child_word, noinputab, child_MMC, child_set) {
        child_vec <- child_set[[1]]
        tibble(
          child_word = child_word,
          noinputab = noinputab,
          child_MMC = child_MMC,
          parent_word = map_chr(seq_along(child_vec), function(i) {
            paste(child_vec[-i], collapse = "_")
          }),
          added_ab = child_vec
        )
      }
    )

    out[[as.character(k)]] <- child_parent_rows %>%
      left_join(parents %>% select(parent_word, parent_MMC), by = "parent_word") %>%
      mutate(delta_MMC = child_MMC - parent_MMC)
  }

  bind_rows(out)
}

summarise_addition_effects <- function(addition_df) {
  addition_df %>%
    group_by(noinputab, added_ab) %>%
    summarise(
      n = sum(!is.na(delta_MMC)),
      mean_delta_MMC = mean(delta_MMC, na.rm = TRUE),
      median_delta_MMC = median(delta_MMC, na.rm = TRUE),
      sd_delta_MMC = sd(delta_MMC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(noinputab, desc(mean_delta_MMC), added_ab)
}

save_plot <- function(plot_obj, filename, width = 10, height = 7) {
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 320
  )
}

# ---------------------------------------------------------------------------
# Load once and filter once
# ---------------------------------------------------------------------------
message("Reading word-level count data...")
word_df_raw <- readCountWordFrameLong()

message("Filtering analysis frame...")
word_df <- word_df_raw %>%
  filter(
    cpmode == TARGET_CPMODE,
    mode == TARGET_MODE,
    is.na(significanceLevel),
    noinputab %in% c(TARGET_NOINPUTAB, TARGET_NOINPUTAB - 1)
  )

# ---------------------------------------------------------------------------
# Compute word-level MMC
# ---------------------------------------------------------------------------
message("Computing word-level MMC...")
word_metrics_all <- compute_word_level_mmc(word_df)
word_metrics <- word_metrics_all %>%
  filter(noinputab %in% TARGET_NOINPUTAB)

# ---------------------------------------------------------------------------
# Presence/absence analysis for each input antibiotic
# ---------------------------------------------------------------------------
message("Computing presence/absence effects...")
presence_df <- make_presence_table(word_metrics)
presence_effects <- compute_presence_effects(presence_df)

write_csv(presence_effects, file.path(OUT_DIR, "presence_effects_vs_MMC.csv"))

presence_plot <- presence_effects %>%
  mutate(input_ab = fct_reorder(input_ab, delta_mean_MMC)) %>%
  ggplot(aes(x = input_ab, y = delta_mean_MMC)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ noinputab, scales = "free_y") +
  labs(
    title = "Association between input antibiotic presence and MMC",
    subtitle = "Delta mean MMC = mean MMC(words with antibiotic) - mean MMC(words without antibiotic)",
    x = "Input antibiotic",
    y = "Delta mean MMC"
  ) +
  theme_bw(base_size = 12)

save_plot(presence_plot, "presence_effects_vs_MMC.png", width = 11, height = 8)

presence_sig_plot <- presence_effects %>%
  mutate(
    neg_log10_padj = -log10(p_adj),
    input_ab = fct_reorder(input_ab, delta_mean_MMC)
  ) %>%
  ggplot(aes(x = delta_mean_MMC, y = neg_log10_padj)) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~ noinputab) +
  labs(
    title = "Effect size versus significance for input antibiotic presence",
    x = "Delta mean MMC",
    y = "-log10 adjusted p-value"
  ) +
  theme_bw(base_size = 12)

save_plot(presence_sig_plot, "presence_effects_volcano_like.png", width = 10, height = 7)

# ---------------------------------------------------------------------------
# Best words for 4, 6, and 8 input antibiotics
# ---------------------------------------------------------------------------
message("Finding best-performing words...")
best_words <- compute_best_words(word_metrics, top_n = 15)
write_csv(best_words, file.path(OUT_DIR, "best_words_by_MMC.csv"))

best_words_plot <- best_words %>%
  group_by(noinputab) %>%
  mutate(word_ranked = forcats::fct_reorder(word, MMC, .desc = FALSE)) %>%
  ungroup() %>%
  ggplot(aes(x = word_ranked, y = MMC)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ noinputab, scales = "free_y") +
  labs(
    title = "Top input antibiotic words by MMC",
    x = "Input antibiotic word",
    y = "MMC"
  ) +
  theme_bw(base_size = 11)

save_plot(best_words_plot, "best_words_by_MMC.png", width = 12, height = 10)

# ---------------------------------------------------------------------------
# Adding one antibiotic: change in MMC from parent word to child word
# ---------------------------------------------------------------------------
message("Computing addition effects...")
addition_df <- compute_addition_effects(word_metrics_all, target_sizes = TARGET_NOINPUTAB)
write_csv(addition_df, file.path(OUT_DIR, "word_addition_delta_MMC_raw.csv"))

addition_effects <- summarise_addition_effects(addition_df)
write_csv(addition_effects, file.path(OUT_DIR, "word_addition_delta_MMC_summary.csv"))

# addition_plot <- addition_effects %>%
#   mutate(added_ab = fct_reorder(added_ab, mean_delta_MMC)) %>%
#   ggplot(aes(x = added_ab, y = mean_delta_MMC)) +
#   geom_col() +
#   coord_flip() +
#   facet_wrap(~ noinputab, scales = "free_y") +
#   labs(
#     title = "Mean change in MMC when adding one input antibiotic",
#     subtitle = "Child word MMC - parent word MMC; panels correspond to resulting word size",
#     x = "Added input antibiotic",
#     y = "Mean delta MMC"
#   ) +
#   theme_bw(base_size = 12)
# 
# save_plot(addition_plot, "word_addition_delta_MMC.png", width = 11, height = 8)
# 
# addition_boxplot <- addition_df %>%
#   filter(!is.na(delta_MMC)) %>%
#   ggplot(aes(x = fct_reorder(added_ab, delta_MMC, .fun = median, na.rm = TRUE), y = delta_MMC)) +
#   geom_boxplot(outlier.size = 0.4) +
#   coord_flip() +
#   facet_wrap(~ noinputab, scales = "free_y") +
#   labs(
#     title = "Distribution of MMC changes after adding one input antibiotic",
#     subtitle = "Parent size is one less than the panel label",
#     x = "Added input antibiotic",
#     y = "Delta MMC"
#   ) +
#   theme_bw(base_size = 12)
# 
# save_plot(addition_boxplot, "word_addition_delta_MMC_boxplot.png", width = 11, height = 8)

# ---------------------------------------------------------------------------
# Compact textual summaries
# ---------------------------------------------------------------------------
summary_best_presence <- presence_effects %>%
  group_by(noinputab) %>%
  slice_max(order_by = delta_mean_MMC, n = 10, with_ties = FALSE) %>%
  ungroup()
write_csv(summary_best_presence, file.path(OUT_DIR, "top_input_antibiotics_by_presence_effect.csv"))

# summary_best_added <- addition_effects %>%
#   group_by(noinputab) %>%
#   slice_max(order_by = mean_delta_MMC, n = 10, with_ties = FALSE) %>%
#   ungroup()
# write_csv(summary_best_added, file.path(OUT_DIR, "top_added_antibiotics_by_delta_MMC.csv"))

message("Done. Outputs written to: ", OUT_DIR)
