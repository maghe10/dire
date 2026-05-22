source('model/modelcommon.R')

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

analyse_long_merged <- function()
{
  # read merged 
  inFile <- "modelOutput_long_merged.rds"
  frame <- readRDS(file.path(getPredictionAltFolder(MODE_A),inFile))
  nrow(frame)
  head(frame)
  
  subframe <- frame %>% filter(noinputab==6) %>% filter(confpred_01 == "") %>% select(antibiotic)
  abCounts <- subframe %>%
    count(antibiotic, sort = TRUE)                                                      
  hist(abCounts$n)
}

readMergedCpRaw <- function()
{
  readRDS(file.path(getCommonModelFolder(),"merged_cp_raw_with_softmax.rds"))
}
  
analyze_labels <-function(frame)
{
  head(frame)
  
  cp_long <- frame %>%
    pivot_longer(
      cols = c(cp_85, cp_90, cp_95, cp_975),
      names_to = "cp_level",
      values_to = "cp_value"
    ) %>%
    mutate(
      label_type = case_when(
        cp_value == ""    ~ "zero",
        cp_value == "S/R" ~ "two",
        TRUE              ~ "one"
      )
    )
  label_counts <- cp_long %>%
    count(antibiotic, cp_level, label_type)
  
  label_presence <- label_counts %>%
    filter(label_type %in% c("zero", "two")) %>%
    mutate(present = n > 0)
  
  label_pct <- cp_long %>%
    count(antibiotic, cp_level, label_type) %>%
    group_by(antibiotic, cp_level) %>%
    mutate(
      pct = 100 * n / sum(n)
    ) %>%
    ungroup()
  ggplot(label_pct %>% filter(label_type %in% c("zero", "two")),
         aes(x = cp_level, y = antibiotic, fill = pct)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.1f%%", pct)), size = 3) +
    facet_wrap(~ label_type, ncol = 1) +
    scale_fill_viridis_c(limits = c(0, 100)) +
    theme_bw()
  
}
  
analyze_cp_metrics <- function(frame)
{
  dens_frame <- frame %>%
    select(antibiotic, ps, pr) %>%
    pivot_longer(cols = c(ps, pr),
                 names_to = "metric",
                 values_to = "value")

  ggplot(dens_frame, aes(x = value, linetype = metric)) +
    geom_density(linewidth = 1) +
    facet_wrap(~ antibiotic, scales = "free_y") +
    labs(x = NULL, y = "Density", linetype = NULL) +
    theme_bw()
}

analyze_correctness_stratified <-function(frame)
{
  df3 <- frame %>%
    mutate(
      outcome = case_when(
        AST_true == "R" & AST_pred == "R" ~ "True R",
        AST_true == "S" & AST_pred == "S" ~ "True S",
        AST_true == "S" & AST_pred == "R" ~ "False R",
        AST_true == "R" & AST_pred == "S" ~ "False S",
        TRUE ~ NA_character_
      ),
      outcome = factor(
        outcome,
        levels = c("True R", "False R", "False S", "True S")
      )
    )
  ggplot(df3, aes(x = psSoftmax, fill = outcome)) +
    geom_density(alpha = 0.4) +
    labs(
      x = "psSoftmax",
      y = "Density",
      fill = "Outcome",
      title = "psSoftmax by prediction outcome (all antibiotics)"
    ) +
    theme_bw()
  
  ggplot(df3, aes(x = psSoftmax, fill = outcome)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~ antibiotic, scales = "free_y") +
    labs(
      x = "psSoftmax",
      y = "Density",
      fill = "Outcome",
      title = "psSoftmax by outcome and antibiotic"
    ) +
    theme_bw()
 
  ggplot(df3, aes(x = outcome, y = psSoftmax, fill = outcome)) +
#    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ antibiotic) +
    labs(
      x = NULL,
      y = "psSoftmax",
      title = "Confidence by outcome class and antibiotic"
    ) +
    theme_bw() +
    theme(legend.position = "none") 
  
  
  agg_outcome <- df3 %>%
    group_by(antibiotic, outcome) %>%
    summarise(
      mean_ps = mean(psSoftmax, na.rm = TRUE),
      median_ps = median(psSoftmax, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  ggplot(agg_outcome, aes(x = antibiotic, y = mean_ps, fill = outcome)) +
    geom_col(position = position_dodge()) +
    labs(
      x = "Antibiotic",
      y = "Mean psSoftmax",
      fill = "Outcome",
      title = "Mean confidence by outcome and antibiotic"
    ) +
    theme_bw()
  
  calibration <- df3 %>%
    mutate(
      correct = outcome %in% c("True R", "True S"),
      bin = cut(psSoftmax, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
    ) %>%
    group_by(antibiotic, bin) %>%
    summarise(
      accuracy = mean(correct),
      n = n(),
      .groups = "drop"
    )
  
  ggplot(calibration, aes(x = bin, y = accuracy, group = antibiotic)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ antibiotic) +
    labs(
      x = "psSoftmax bin",
      y = "Empirical accuracy",
      title = "Calibration of psSoftmax by antibiotic"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
}


analyze_correctness_conformal_pr_ps <- function(frame)
{
  df4 <- frame %>%
    mutate(
      correct = AST_true == AST_pred,
      correct = factor(correct, levels = c(FALSE, TRUE),
                       labels = c("Incorrect", "Correct"))
    )
  
  ggplot(df4, aes(x = ps, fill = correct)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~ antibiotic, scales = "free_y") +
    labs(
      x = "ps",
      y = "Density",
      fill = "Prediction",
      title = "ps vs correctness by antibiotic"
    ) +
    theme_bw()
  
  ggplot(df4, aes(x = pr, fill = correct)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~ antibiotic, scales = "free_y") +
    labs(
      x = "pr",
      y = "Density",
      fill = "Prediction",
      title = "pr vs correctness by antibiotic"
    ) +
    theme_bw()
  
  ggplot(df4, aes(x = correct, y = ps, fill = correct)) +
#    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ antibiotic) +
    labs(
      x = NULL,
      y = "ps",
      title = "ps by correctness and antibiotic"
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  ggplot(df4, aes(x = correct, y = pr, fill = correct)) +
    #    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ antibiotic) +
    labs(
      x = NULL,
      y = "pr",
      title = "pr by correctness and antibiotic"
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  agg_pspr <- df4 %>%
    group_by(antibiotic, correct) %>%
    summarise(
      mean_ps = mean(ps, na.rm = TRUE),
      median_ps = median(ps, na.rm = TRUE),
      mean_pr = mean(pr, na.rm = TRUE),
      median_pr = median(pr, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  agg_long <- agg_pspr %>%
    tidyr::pivot_longer(
      cols = c(mean_ps, mean_pr),
      names_to = "metric",
      values_to = "value"
    )
  
  ggplot(agg_long, aes(x = antibiotic, y = value, fill = correct)) +
    geom_col(position = position_dodge()) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(
      x = "Antibiotic",
      y = "Mean value",
      fill = "Prediction",
      title = "Mean ps and pr by correctness and antibiotic"
    ) +
    theme_bw()
  
  cnt <- df4 %>%
    mutate(
      ps_high = ps > 0.1,
      pr_high = pr > 0.1
    ) %>%
    count(ps_high, pr_high, correct)
  
  cnt_pct <- df4 %>%
    mutate(
      ps_high = ps > 0.1,
      pr_high = pr > 0.1
    ) %>%
    count(ps_high, pr_high, correct) %>%
    group_by(correct) %>%
    mutate(
      pct = 100 * n / sum(n)
    ) %>%
    ungroup()
  
  ggplot(cnt_pct, aes(x = ps_high, y = pr_high, fill = pct)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f%%", pct)), color = "black") +
    facet_wrap(~ correct) +
    scale_fill_viridis_c(
      name = "Percentage",
      limits = c(0, 100)
    ) +
    labs(
      x = "ps > 0.1",
      y = "pr > 0.1",
      title = "Correctness by ps/pr threshold regions (percentages)"
    ) +
    theme_bw()
  
  
}


analyze_correctness <-function(frame)
{
  df2 <- frame %>%
    mutate(
      correct = (AST_true == AST_pred),
      correct = factor(correct, levels = c(FALSE, TRUE),
                       labels = c("Incorrect", "Correct"))
    )
    
  
  ggplot(df2, aes(x = psSoftmax, fill = correct)) +
    geom_density(alpha = 0.4) +
    labs(
      x = "psSoftmax",
      y = "Density",
      fill = "Prediction",
      title = "Confidence vs correctness (all antibiotics)"
    ) +
    theme_bw()
  
  ggplot(df2, aes(x = psSoftmax, fill = correct)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~ antibiotic, scales = "free_y") +
    labs(
      x = "psSoftmax",
      y = "Density",
      fill = "Prediction",
      title = "psSoftmax vs correctness by antibiotic"
    ) +
    theme_bw()
  
  ggplot(df2, aes(x = correct, y = psSoftmax, fill = correct)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    facet_wrap(~ antibiotic) +
    labs(
      x = NULL,
      y = "psSoftmax",
      title = "Distribution of psSoftmax by correctness and antibiotic"
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  
  agg_ps <- df2 %>%
    group_by(antibiotic, correct) %>%
    summarise(
      mean_ps = mean(psSoftmax, na.rm = TRUE),
      median_ps = median(psSoftmax, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  ggplot(agg_ps, aes(x = antibiotic, y = mean_ps, fill = correct)) +
    geom_col(position = position_dodge()) +
    labs(
      x = "Antibiotic",
      y = "Mean psSoftmax",
      fill = "Prediction",
      title = "Mean confidence by antibiotic and correctness"
    ) +
    theme_bw()
  
  # ggplot(df2, aes(x = factor(sample), y = psSoftmax, color = correct)) +
  #   geom_point(alpha = 0.6) +
  #   facet_wrap(~ antibiotic, scales = "free_x") +
  #   labs(
  #     x = "Sample",
  #     y = "psSoftmax",
  #     color = "Prediction",
  #     title = "Per-sample confidence vs correctness"
  #   ) +
  #   theme_bw() +
  #   theme(axis.text.x = element_blank(),
  #         axis.ticks.x = element_blank())
  calibration <- df2 %>%
    mutate(bin = cut(psSoftmax, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
    group_by(antibiotic, bin) %>%
    summarise(
      accuracy = mean(correct == "Correct"),
      n = n(),
      .groups = "drop"
    )
  
  ggplot(calibration, aes(x = bin, y = accuracy, group = antibiotic)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ antibiotic) +
    labs(
      x = "psSoftmax bin",
      y = "Empirical accuracy",
      title = "Calibration of psSoftmax by antibiotic"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
}


analyze_diff_pr_ps <- function(frame)
{
  df4 <- frame %>%
    mutate(
      correct = AST_true == AST_pred,
      correct = factor(correct, levels = c(FALSE, TRUE),
                       labels = c("Incorrect", "Correct"))
    )
  
  df5 <- df4 %>%
    mutate(
      delta_pr_ps = pr - ps,
      abs_delta   = abs(delta_pr_ps)
    )
  
  
  ggplot(df5, aes(x = delta_pr_ps, fill = correct)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      x = "pr − ps",
      y = "Density",
      fill = "Prediction",
      title = "Decision margin (pr − ps) by correctness"
    ) +
    theme_bw()
  
  df_ab <- df4 %>%
    mutate(
      delta_pr_ps = pr - ps,
      abs_delta   = abs(delta_pr_ps),
      correct = factor(correct, levels = c("Incorrect", "Correct"))
    )

    ggplot(df_ab, aes(x = delta_pr_ps, fill = correct)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~ antibiotic, scales = "free_y") +
    labs(
      x = "pr − ps",
      y = "Density",
      fill = "Prediction",
      title = "Decision margin (pr − ps) by antibiotic"
    ) +
    theme_bw()
  
    
    limit <- 0.025
    
    ggplot(df4 %>% filter( pr> limit ), aes(x = ps)) +
      geom_density(alpha = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = limit, linetype = "dashed") +
      facet_wrap(~ antibiotic, scales = "free_y") +
      labs(
        x = "ps",
        y = "Density",
        fill = "Prediction",
        title = sprintf("ps by antibiotic where pr is over %.3f",limit)
      ) +
      theme_bw()

    ggplot(df4 %>% filter( ps> limit), aes(x = pr)) +
      geom_density(alpha = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = limit, linetype = "dashed") +
      facet_wrap(~ antibiotic, scales = "free_y") +
      labs(
        x = "pr",
        y = "Density",
        fill = "Prediction",
        title = sprintf("pr by antibiotic where ps is over %.3f",limit)
      ) +
      theme_bw()
    
}


RUN <- function()
{
  df <- readMergedCpRaw()

  frame <- df %>% filter(mode==MODE_A & noinputab == 6)
  
  head(df)
  analyze_labels(frame)
  analyze_cp_metrics(frame)
  analyze_correctness(frame)
  analyze_correctness_conformal_pr_ps(frame)
  analyze_correctness_stratified(frame)
  analyze_diff_pr_ps(frame)
  
}