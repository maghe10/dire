source(file = 'manuscript/plotStatisticsTables_common.R')
library(magick)



dir.create(manuscriptPlotDirectory, recursive = TRUE, showWarnings = FALSE)


CLUSTER_PNG_OUTPUT_DIR <- paste(processedRootR,"cluster","png",sep="/")





metric_colors <- c(correct = "#70AD47", ME = "#FFFF00", VME = "#FF0000")

input_antibiotics_legend_title <- "Number of input antibiotics"
significance_level_legend_title <- "Significance level"

unambiguous_correct_predictions_y_legend <- "Unambiguous correct predictions"
correct_predictions_y_legend <- "Correct predictions"

percentages_results_y_legend <- "Prediction results"
percentages_unambiguous_results_y_legend <- "Unambiguous prediction results"


metric_level_labels <- c("Correct",
                      "Major error",
                      "Very major error")






plot_one_metric_grouped <- function(metrics_filtered,
                                     metric_name,
                                     tag,
                                     ylab,
                                     id_col,         # string, e.g. "antibiotic" or "ab_group"
                                     y_as_percent,
                                     column_var,     # string, e.g. "noinputab" or "significanceLevel"
                                     id_levels = NULL,
                                     column_levels = NULL,
                                     column_labels = column_levels,
                                     legend_title = NULL) {
  
  df_wide <- metrics_filtered %>%
    dplyr::filter(metric == metric_name) %>%
    dplyr::select(dplyr::all_of(c(id_col, column_var, "value"))) %>%
    dplyr::mutate(
      .col = as.character(.data[[column_var]])
    ) %>%
    dplyr::select(dplyr::all_of(id_col), .col, value) %>%
    tidyr::pivot_wider(names_from = .col, values_from = value)
  
  if (is.null(column_levels)) {
    column_levels <- setdiff(names(df_wide), id_col)
  }
  
  p <- plot_grouped_bars(
    df_wide,
    id_col = id_col,
    measure_cols = column_levels,
    id_levels = id_levels,
    group_levels = column_levels,
    group_labels = column_labels,
    legend_title = ifelse(is.null(legend_title), column_var, legend_title),
    y_label = ylab,
    y_as_percent = y_as_percent,
    y_accuracy = ifelse(y_as_percent, 1,0.05),
    tag = tag,
    export = FALSE
  )
  
  # rotate only if it's the antibiotic axis (optional tweak)
  if (id_col == "antibiotic") {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  }
  
  p
}

plot_metric_panels_grouped <- function(metrics_long,
                                        metric_defs,        # your PERFORMANCE_METRICS_DEFS: metric, tag, ylab
                                        id_col,             # string
                                        column_var,         # string
                                        filter_fun,         # function(df) -> filtered df
                                        fig_stub,
                                        ncol = 2,
                                        out_dir = manuscriptPlotDirectory,
                                        export = TRUE,
                                        id_levels = NULL,
                                        column_levels = NULL,
                                        column_labels = column_levels,
                                        legend_title = NULL) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork is required.", call. = FALSE)
  }
  
  metrics_filtered <- filter_fun(metrics_long)
  
  plots <- vector("list", nrow(metric_defs))
  for (i in seq_len(nrow(metric_defs))) {
    plots[[i]] <- plot_one_metric_grouped(
      metrics_filtered = metrics_filtered,
      metric_name = metric_defs$metric[[i]],
      tag = metric_defs$tag[[i]],
      ylab = metric_defs$ylab[[i]],
      y_as_percent = metric_defs$y_as_percent[[i]],
      id_col = id_col,
      column_var = column_var,
      id_levels = id_levels,
      column_levels = column_levels,
      column_labels = column_labels,
      legend_title = legend_title
    )
  }
  
  combined <- patchwork::wrap_plots(plots, ncol = ncol)
  print(combined)
  
  if (isTRUE(export)) {
    export_patchwork(
      patch = combined,
      out_dir = out_dir,
      file_stub = fig_stub,
      width = 13,
      height = 8,
      dpi = 600,
      export = TRUE
    )
  }
  
  invisible(combined)
}

plotAntibioticsMetricsVsNoinputAb <- function(tag)
{
  metrics_long_ab <- derivedMetricsFrame()
  
  aPlot <- plot_metric_panels_grouped(
    metrics_long = metrics_long_ab,
    metric_defs  = PERFORMANCE_METRICS_DEFS,
    id_col       = "antibiotic",
    column_var   = "noinputab",
    filter_fun   = function(df) {
      df %>% dplyr::filter(
        significanceLevel == "STD",
        mode == "Mode-A",
        cpmode == "normal"
      )
    },
    fig_stub     = sprintf("Figure_%s", tag),
    ncol         = 2,
    export       = TRUE,
    legend_title = input_antibiotics_legend_title ,
    id_levels = ANTIBIOTICS_LEVELS
  )
  
  aPlot
}


plotAntibioticsMetricsVsNoinputAb_Abgroups <- function(tag) {
  
  metrics_long_ab <- derivedMetricsFrame(
    prefix = NA,
    countFrameWide = readCountFrameWideSumByGroup(),
    vars = AB_GROUP_GROUPS_VAR
  )
  
  plot_metric_panels_grouped(
    metrics_long = metrics_long_ab,
    metric_defs  = PERFORMANCE_METRICS_DEFS,
    id_col       = "ab_group",
    column_var   = "noinputab",
    filter_fun   = function(df) {
      df %>%
        dplyr::filter(
          significanceLevel == "STD",
          mode == "Mode-A",
          cpmode == "normal"
        )
    },
    
    fig_stub     = sprintf("Figure_%s", tag),
    ncol         = 2,
    export       = TRUE,
    legend_title = input_antibiotics_legend_title ,
    id_levels = AB_GROUPS_LEVELS
  )
}

plotAntibioticsMetricsVsSignificanceLevel_Abgroups <- function(tag) {
  
  metrics_long_ab <- derivedMetricsFrame(
    prefix = NA,
    countFrameWide = readCountFrameWideSumByGroup(),
    vars = AB_GROUP_GROUPS_VAR
  )
  column_labels = setNames(CONFIDENCE_LEVELS, SIGNIFICANCE_LEVELS)
  
  plot_metric_panels_grouped(
    metrics_long = metrics_long_ab,
    metric_defs  = PERFORMANCE_METRICS_DEFS,
    id_col       = "ab_group",
    column_var   = "significanceLevel",
    filter_fun = function(df) {
      df %>%
        dplyr::mutate(
          significanceLevel = factor(as.character(significanceLevel), levels = SIGNIFICANCE_LEVELS)
        ) %>%
        dplyr::filter(
          noinputab == 6,
          mode == "Mode-A",
          cpmode == "normal"
        )
    },
    fig_stub     = sprintf("Figure_%s", tag),
    ncol         = 2,
    export       = TRUE,
    legend_title = significance_level_legend_title,
    column_levels = SIGNIFICANCE_LEVELS,
    column_labels = column_labels,
    id_levels = AB_GROUPS_LEVELS
  )
}


plotAntibioticsMetricsVsMode <- function(tag) {
  
  metrics_long_ab <- derivedMetricsFrame()
  
  plot_metric_panels_grouped(
    metrics_long = metrics_long_ab,
    metric_defs  = PERFORMANCE_METRICS_DEFS,
    id_col       = "antibiotic",
    column_var   = "mode",
    filter_fun   = function(df) {
      df %>% dplyr::filter(significanceLevel == "STD", noinputab == 6, cpmode == "normal")
    },
    fig_stub     = sprintf("Figure_%s", tag),
    ncol         = 2,
    export       = TRUE,
    legend_title = "Mode",
    id_levels = ANTIBIOTICS_LEVELS 
  )
}

plotAntibioticsMetricsVsNoinputAbOnePerAbGroup <- function(tag) {
  
  metrics_long_ab <- derivedMetricsFrame("oneperabgroup")
  
  plot_metric_panels_grouped(
    metrics_long = metrics_long_ab,
    metric_defs  = PERFORMANCE_METRICS_DEFS,
    id_col       = "antibiotic",
    column_var   = "noinputab",
    filter_fun   = function(df) {
      df %>% dplyr::filter(significanceLevel == "STD", mode == "Mode-A", cpmode == "normal")
    },
    fig_stub     = sprintf("Figure_%s", tag),
    ncol         = 2,
    export       = TRUE,
    legend_title = input_antibiotics_legend_title,
    id_levels = ANTIBIOTICS_LEVELS
  )
}


plotAntibioticsMetricsVsSignificanceLevel <- function(tag) {
  
  metrics_long_ab <- derivedMetricsFrame()
  column_labels = setNames(CONFIDENCE_LEVELS, SIGNIFICANCE_LEVELS)
  
  plot_metric_panels_grouped(
    metrics_long = metrics_long_ab,
    metric_defs  = PERFORMANCE_METRICS_DEFS,
    id_col       = "antibiotic",
    column_var   = "significanceLevel",
    filter_fun = function(df) {
      df %>%
        dplyr::mutate(
          significanceLevel = factor(as.character(significanceLevel), levels = SIGNIFICANCE_LEVELS)
        ) %>%
        dplyr::filter(
          noinputab == 6,
          mode == "Mode-A",
          cpmode == "normal"
        )
    },
    fig_stub     = sprintf("Figure_%s", tag),
    ncol         = 2,
    export       = TRUE,
    legend_title = significance_level_legend_title,
    column_levels = SIGNIFICANCE_LEVELS,   # if your helper supports this, set it too
    column_labels = column_labels,
    id_levels = ANTIBIOTICS_LEVELS
  )
}

plotAntibioticsMetricsVsSignificanceLevelSixBest <- function(tag) {
  
  countFrameWide = readCountWordFrameWide() %>% filter_and_drop(word,"AMC_TZP_CTX_CIP_OFX_TOB") 
  metrics_long_ab <- derivedMetricsFrame(countFrameWide=countFrameWide)
  #metrics_long_ab <- derivedMetricsFrame("sixbest")
  column_labels = setNames(CONFIDENCE_LEVELS, SIGNIFICANCE_LEVELS)
  
  predicted_antibiotics <- unique(metrics_long_ab$antibiotic)
  
  predicted_antibiotics_ordered <- ANTIBIOTICS_LEVELS[ANTIBIOTICS_LEVELS %in% predicted_antibiotics]
  
  #mean(metrics_long_ab[metrics_long_ab$metric == "MCC",]$value,na.rm = TRUE)
  plot_metric_panels_grouped(
    metrics_long = metrics_long_ab,
    metric_defs  = PERFORMANCE_METRICS_DEFS,
    id_col       = "antibiotic",
    column_var   = "significanceLevel",
    filter_fun = function(df) {
      df %>%
        dplyr::mutate(
          significanceLevel = factor(as.character(significanceLevel), levels = SIGNIFICANCE_LEVELS)
        ) %>%
        dplyr::filter(
          noinputab == 6,
          mode == "Mode-A",
          cpmode == "normal"
        )
    },
    fig_stub     = sprintf("Figure_%s", tag),
    ncol         = 2,
    export       = TRUE,
    legend_title = significance_level_legend_title,
    column_levels = SIGNIFICANCE_LEVELS,   # if your helper supports this, set it too
    column_labels =  column_labels,
    id_levels = predicted_antibiotics_ordered
  )
}

plotF1AbgroupsVsSet <- function(tag) {
  # F1, 6 ab, non-conformal
  metrics_long_evaluation <- derivedMetricsFrame(
    prefix = NA,
    countFrameWide = readCountFrameWideSumByGroup(),
    vars = AB_GROUP_GROUPS_VAR
  ) %>% 
    filter_and_drop(mode, "Mode-A") %>% 
    filter_and_drop(cpmode ,"normal") %>% 
    filter_and_drop(noinputab ,6) %>% 
    filter_and_drop(metric,"F1") %>% 
    filter_and_drop(significanceLevel ,"STD")

  # > metrics_long_evaluation
  # # A tibble: 4 × 2
  # ab_group         value
  # <chr>            <dbl>
  #   1 Aminoglycosides  0.574
  # 2 Cephalosporins   0.850
  # 3 Fluoroquinolones 0.867
  # 4 Penicillins      0.835  
  
  metrics_long_reference <- tibble::tribble(
    ~ab_group,          ~value,
    "Aminoglycosides",   0.65,
    "Cephalosporins",    0.91,
    "Fluoroquinolones",  0.86,
    "Penicillins",       0.74
  ) |>
    dplyr::mutate(ab_group = factor(ab_group, levels = AB_GROUPS_LEVELS))
  
  merged_long <- dplyr::bind_rows(
    metrics_long_evaluation %>% dplyr::mutate(source = "evaluation"),
    metrics_long_reference  %>% dplyr::mutate(source = "reference")
  )
  
  merged_long <- merged_long %>%
    dplyr::rename(set = source) %>%
    dplyr::mutate(set = factor(set, levels = c("evaluation", "reference")))
  
  set_levels <- c("evaluation", "reference")

  merged_long <- merged_long %>%
    dplyr::mutate(
      ab_group = factor(as.character(ab_group), levels = AB_GROUPS_LEVELS),
      set = factor(as.character(set), levels = set_levels)
    )
  
  df_wide <- merged_long %>%
    mutate(set = as.character(set)) %>%
    pivot_wider(names_from = set, values_from = value)
  
  p <- plot_grouped_bars(
    df_wide,
    id_col       = "ab_group",
    file_stub =  sprintf("Figure_%s",tag),
    measure_cols = set_levels,
    group_levels = set_levels,
    legend_title = "Set",
    y_label      = "F1 score",
    export       = TRUE,
    id_levels = AB_GROUPS_LEVELS,
    y_as_percent = FALSE,
    y_accuracy = 0.05,   # e.g. show 0.00, 0.05, 0.10 ...
    y_limits = c(0, 1)
  )
  
  print(p)
  

}

plotMetricsVsSet <- function(tag) {
  
  metrics_long_evaluation <- derivedMetricsFrame(prefix = NA,countFrameWide = readCountFrameWideSum()) %>% 
    filter_and_drop(mode, "Mode-A") %>% 
    filter_and_drop(cpmode ,"normal") %>% 
    filter_and_drop(noinputab ,6) %>% 
    dplyr::filter(metric %in% c("ME","VME")) %>% 
    dplyr::filter(significanceLevel %in% c("5%","10%"))
  
  metrics_long_reference <- tibble(
    significanceLevel = factor(
      c("5%", "5%",
        "10%", "10%"),
      levels = SIGNIFICANCE_LEVELS
    ),
    metric = c("ME", "VME",
               "ME", "VME"),
    value = c(0.051, 0.05,
              0.101, 0.097)
  )
  
  merged_long <- dplyr::bind_rows(
    metrics_long_evaluation %>% dplyr::mutate(source = "evaluation"),
    metrics_long_reference  %>% dplyr::mutate(source = "reference")
  )
  
  merged_long <- merged_long %>%
    dplyr::rename(set = source) %>%
    dplyr::mutate(set = factor(set, levels = c("evaluation", "reference")))
  

  
  # (Optional but recommended) make sure ordering is stable
  # If you have SIGNIFICANCE_LEVELS constant, use it; otherwise specify manually.
  sig_levels <- c("5%", "10%")        # or SIGNIFICANCE_LEVELS
  set_levels <- c("evaluation", "reference")
  
  # Ensure factors follow desired order
  merged_long <- merged_long %>%
    dplyr::mutate(
      significanceLevel = factor(as.character(significanceLevel), levels = sig_levels),
      set = factor(as.character(set), levels = set_levels),
      metric = as.character(metric)
    )
  
  make_one_metric_plot <- function(metric_name, tag, ylab) {
    df_wide <- merged_long %>%
      dplyr::filter(metric == metric_name) %>%
      dplyr::select(significanceLevel, set, value) %>%
      dplyr::mutate(set = as.character(set)) %>%
      tidyr::pivot_wider(names_from = set, values_from = value)
    
    plot_grouped_bars(
      df_wide,
      id_col       = "significanceLevel",
      id_levels    = sig_levels,
      measure_cols = set_levels,
      group_levels = set_levels,
      legend_title = "Set",
      y_label      = ylab,
      tag          = tag,
      export       = FALSE
    ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5)
      )
  }
  
  pA <- make_one_metric_plot(
    metric_name = VME_ME_METRICS_DEFS$metric[[1]],
    tag         = VME_ME_METRICS_DEFS$tag[[1]],
    ylab        = VME_ME_METRICS_DEFS$ylab[[1]]
  )
  
  pB <- make_one_metric_plot(
    metric_name = VME_ME_METRICS_DEFS$metric[[2]],
    tag         = VME_ME_METRICS_DEFS$tag[[2]],
    ylab        = VME_ME_METRICS_DEFS$ylab[[2]]
  )
  
  pAB <- patchwork::wrap_plots(pA, pB, ncol = 2)
  print(pAB)
  export_patchwork(patch = pAB,out_dir = manuscriptPlotDirectory,sprintf("Figure_%s",tag))
  
  
}


plotUnambiguousCorrectPredictionsAndErrors <- function(tag)
{
  metrics_long <- derivedMetricsFrameUnambiguous(prefix = NA,countFrameWide = readCountFrameWideSum()) %>% 
    filter_and_drop(mode, "Mode-A") %>% filter_and_drop(cpmode ,"normal") 
  
  frameCorrectSignificanceLevelVsIndexModeA <- metrics_long %>% filter_and_drop(metric, "correct") %>% 
    pivot_wider(names_from = noinputab, values_from = value) %>%
    slice(n():1) 
  
  frameCorrectSignificanceLevelVsIndexModeA$significanceLevel <- SIGNIFICANCE_LEVELS
#  frameCorrectSignificanceLevelVsIndexModeA$confidenceLevel <- CONFIDENCE_LEVELS
  
  id_labels = setNames(CONFIDENCE_LEVELS, SIGNIFICANCE_LEVELS)
  
  
  pA <- plot_grouped_bars(
    frameCorrectSignificanceLevelVsIndexModeA,
    id_col = "significanceLevel",
    id_levels = SIGNIFICANCE_LEVELS,
    id_labels = id_labels,
    group_levels = NO_INPUT_LEVELS,
    legend_title = input_antibiotics_legend_title,
    export = FALSE
  )
  
  
  ncol <- 2

  long <- metrics_long %>%
    dplyr::mutate(
      noinputab = factor(noinputab, NO_INPUT_LEVELS),
      metric = factor(metric, levels = METRIC_LEVELS),
      significanceLevel = factor(significanceLevel, levels = SIGNIFICANCE_LEVELS)
    )
  # 4) Plot
  
  confidence_labels <- setNames(CONFIDENCE_LEVELS, SIGNIFICANCE_LEVELS)
  
  
  pB <- plot_stacked_bars(
    long,
    x_col = "noinputab",
    fill_col = "metric",
    values_to = "value",
    x_levels = NO_INPUT_LEVELS,
    fill_levels = METRIC_LEVELS,
    fill_colors = metric_colors,
    fill_labels = metric_level_labels,
    y_label = percentages_unambiguous_results_y_legend,
    bar_width = 0.8,
    clamp = TRUE,
    stack_reverse = TRUE,
    export = FALSE
  ) +
#    ggplot2::facet_wrap(~ significanceLevel, ncol = ncol, drop = FALSE) +
    ggplot2::facet_wrap(
      ~ significanceLevel,
      ncol = ncol,
      drop = FALSE,
      labeller = ggplot2::as_labeller(confidence_labels)
    ) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(0.75, "lines")
    )
  
  pAB <- combine_tagged(pA, pB, ncol = 2, tag_levels = "A", base_size = 12)
  print(pAB)
  export_patchwork(patch = pAB,out_dir = manuscriptPlotDirectory,sprintf("Figure_%s",tag))
  
}

plotUnambiguousCorrectPredictionsAndErrorsVersusIndexAntibioticsLevel <- function(tag,significanceLevel = "STD")
{

  metrics_long <- derivedMetricsFrameUnambiguous(prefix = NA,countFrameWide = readCountFrameWide()) %>% 
    filter_and_drop(mode, "Mode-A") %>% filter_and_drop(cpmode ,"normal") 
  
  frameCorrectAntibioticsVsIndex <- metrics_long %>% filter_and_drop(metric, "correct") %>% filter_and_drop(significanceLevel, significanceLevel) %>%  
    pivot_wider(names_from = noinputab, values_from = value) 
  
  long_for_stacked <- metrics_long %>% filter_and_drop(significanceLevel, significanceLevel)

  # --- Panel A: grouped bars (antibiotic x #inputs)
  pA <- plot_grouped_bars(
    frameCorrectAntibioticsVsIndex,
    id_col = "antibiotic",
    id_levels = ANTIBIOTICS_LEVELS,
    measure_cols = NO_INPUT_LEVELS,
    group_levels = NO_INPUT_LEVELS,
    dodge_width = 0.78,
    bar_width = 0.66,
    legend_title = input_antibiotics_legend_title,
    y_label = ifelse(
      significanceLevel == "STD",
      correct_predictions_y_legend,
      unambiguous_correct_predictions_y_legend
    ),
    
    y_limits = c(0, 1),
    tag = "A",
    export = FALSE
  ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7, angle = 45, hjust = 1, vjust = 1))
  
  # --- Panel B: STACKED needs long form for our generic function (self-contained)
  ncol <- 5
  long <- pad_facets_to_full_row(long = long_for_stacked, facet_col = "antibiotic", ncol = ncol,
                                 x_col = "noinputab", fill_col = "metric", values_col = "value")
  
  pad_labeller <- ggplot2::labeller(
    antibiotic = function(x) ifelse(grepl("^__PAD\\d+__$", x), "", x)
  )
  
  
  pB <- plot_stacked_bars(
    long,
    x_col = "noinputab",
    fill_col = "metric",
    values_to = "value",
    x_levels = NO_INPUT_LEVELS,
    fill_levels = METRIC_LEVELS,
    fill_colors = metric_colors,
    fill_labels = metric_level_labels,
    y_label = ifelse(
      significanceLevel == "STD",
      percentages_results_y_legend,
      percentages_unambiguous_results_y_legend
    ),
    y_limits = c(0, 1),
    clamp = TRUE,
    stack_reverse = TRUE,
    tag = "B",
    export = FALSE,
    bar_width = 0.8
  ) +
    ggplot2::facet_wrap(~ antibiotic, ncol = ncol, drop = FALSE, labeller = pad_labeller) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7, angle = 0, hjust = 0.5, vjust = 0.5),
      strip.text  = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(0.75, "lines")
    )
  
  # --- Combine
  pAB <- combine_tagged(pA, pB, ncol = 2, tag_levels = "A", base_size = 12)
  print(pAB)

  export_patchwork(patch = pAB,out_dir = manuscriptPlotDirectory,sprintf("Figure_%s",tag))
}


plotUnambiguousCorrectPredictionsAndErrorsSignificanceLevelAntibioticsLevel <- function(tag)
{

  metrics_long <- derivedMetricsFrameUnambiguous(prefix = NA,countFrameWide = readCountFrameWide()) %>% 
    filter_and_drop(mode, "Mode-A") %>% filter_and_drop(cpmode ,"normal") 
  
  frameCorrectAntibioticsVsSignificanceLevel <- metrics_long %>% filter_and_drop(metric, "correct") %>% filter_and_drop(noinputab, 6) %>% 
    pivot_wider(names_from = significanceLevel, values_from = value) 
  
  long_for_stacked <- metrics_long %>% filter_and_drop(noinputab, 6)
  
  confidence_labels <- setNames(CONFIDENCE_LEVELS, SIGNIFICANCE_LEVELS)
  
  pA <- plot_grouped_bars(
    frameCorrectAntibioticsVsSignificanceLevel,
    id_col = "antibiotic",
    id_levels = ANTIBIOTICS_LEVELS,
    group_levels = SIGNIFICANCE_LEVELS,
    group_labels = confidence_labels,
    legend_title = significance_level_legend_title,
    y_label = unambiguous_correct_predictions_y_legend,
    dodge_width = 0.78,
    bar_width = 0.66,
    y_limits = c(0, 1),
    export = FALSE
  )
  # --- Panel B: STACKED needs long form for our generic function (self-contained)
  ncol <- 5

  long <- pad_facets_to_full_row(long = long_for_stacked, facet_col = "antibiotic", ncol = ncol,
                                         x_col = "significanceLevel", fill_col = "metric", values_col = "value")
  
  pad_labeller <- ggplot2::labeller(
    antibiotic = function(x) ifelse(grepl("^__PAD\\d+__$", x), "", x)
  )
  
  
  pB <- plot_stacked_bars(
    long,
    x_col = "significanceLevel",
    fill_col = "metric",
    values_to = "value",
    x_levels = SIGNIFICANCE_LEVELS,
    fill_levels = METRIC_LEVELS,
    fill_colors = metric_colors,
    fill_labels = metric_level_labels,
    y_label = percentages_unambiguous_results_y_legend,
    y_limits = c(0, 1),
    clamp = TRUE,
    stack_reverse = TRUE,
    tag = "B",
    export = FALSE,
    bar_width = 0.8
  ) +
    ggplot2::facet_wrap(~ antibiotic, ncol = ncol, drop = FALSE, labeller = pad_labeller) +
    ggplot2::scale_x_discrete(
      labels = confidence_labels,
      expand = c(0, 0)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7, angle = 0, hjust = 0.5, vjust = 0.5),
      strip.text  = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(0.75, "lines")
    )
  
  # --- Combine
  pAB <- combine_tagged(pA, pB, ncol = 2, tag_levels = "A", base_size = 12)
  print(pAB)
  
  export_patchwork(patch = pAB,out_dir = manuscriptPlotDirectory,sprintf("Figure_%s",tag))
}


plot_sir_distribution <- function(sir_wide,
                                  sample_col = "sample",
                                  antibiotic_order = NULL,
                                  out_file = NULL,
                                  width = 12,
                                  height = 6) {
  long <- sir_wide %>%
    pivot_longer(
      cols = -all_of(sample_col),
      names_to = "antibiotic",
      values_to = "sir"
    ) %>%
    filter(!is.na(sir), sir %in% c("S", "I", "R")) %>%
    mutate(
      sir = factor(sir, levels = c("S", "I", "R"))
    )
  
  if (is.null(antibiotic_order)) {
    antibiotic_order <- long %>%
      distinct(antibiotic) %>%
      pull(antibiotic)
  }
  
  plot_df <- long %>%
    mutate(antibiotic = factor(antibiotic, levels = antibiotic_order))
  
  p <- ggplot(plot_df, aes(x = antibiotic, fill = sir)) +
    geom_bar(position = "fill", width = 0.72, colour = NA) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(
      values = c(
        "S" = "#70AD47",  # green
        "I" = "#FFFF00",  # yellow
        "R" = "#FF0000"   # red
      ),
      breaks = c("S", "I", "R")
    ) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "bottom"
    ) + theme_minimal() +
    theme(
      panel.grid = element_blank()
    )
  
  if (!is.null(out_file)) {
    ggsave(out_file, p, width = width, height = height, dpi = 600)
  }
  
  p
}


plotSirAllAntibiotics <- function(tag)
{
  sir_tbl <- read.csv2(
    file.path(processedRootRcommon, "sirAntibiotics.csv"), stringsAsFactors = FALSE)
  
  ab_order <- c(
    "F", "MEL", "CFR", "W", "SXT", "MEM",
    "AMP", "AMC", "PIP", "TZP",
    "CAZ", "CRO", "CTX", "FEP",
    "CIP", "OFX", "LVX", "MFX",
    "TOB", "GEN"
  )
  
  plot_sir_distribution(
    sir_wide = sir_tbl,
    sample_col = "sample",
    antibiotic_order = ab_order,
    out_file =   file.path(manuscriptPlotDirectory,sprintf("Figure_%s.png",tag))
  )
}

make_included_vs_estimated_plot <-function()
{
  patients <- read.csv2(file.path(processedRootRcommon, "patients.csv"))
  
  patients_long <- patients %>%
    mutate(
      Age = factor(Age, levels = c("<18", "18-64", ">=65")),
      sex = factor(sex, levels = c("Female", "Male")),
      id = interaction(Age, sex, sep = "_")
    ) %>%
    select(id, source, prop) %>%
    pivot_wider(names_from = source, values_from = prop)

#   pA <- plot_grouped_bars(
#     patients_long,
#     id_col = "id",
#     measure_cols = c("Estimated", "Included"),
#     group_levels = c("Estimated", "Included"),
#     y_limits = c(0, 0.6),
#     y_label = NULL,
#     legend_title = NULL,
#     tag = "A"
#   ) +
# #    scale_fill_manual(values = c("Estimated" = "#00B0F0", "Included" = "#FF0000")) +
#     scale_x_discrete(
#       labels = c(
#         "<18\nFemale", "<18\nMale",
#         "18-64\nFemale", "18-64\nMale",
#         "≥65\nFemale", "≥65\nMale"
#       )
#     ) +
#     theme(
#       axis.text.x = element_text(size = 12, lineheight = 0.9),
#       legend.position = "bottom",
#       )
#   
  pA <- plot_grouped_bars(
    patients_long,
    id_col = "id",
    measure_cols = c("Estimated", "Included"),
    group_levels = c("Estimated", "Included"),
    y_limits = c(0, 0.6),
    y_label = NULL,
    legend_title = NULL,
    tag = "A",
    id_labels = c(
      "<18\nFemale", "<18\nMale",
      "18-64\nFemale", "18-64\nMale",
      "≥65\nFemale", "≥65\nMale"
    )
  ) +
    theme(
      axis.text.x = element_text(size = 12, lineheight = 0.9),
      legend.position = "bottom"
    )
  
  
  samples <- read.csv2(file.path(processedRootRcommon, "samples.csv"))
  
  pB <- plot_grouped_bars(
    samples,
    id_col = "Category",
    measure_cols = c("Estimated", "Included"),
    group_levels = c("Estimated", "Included"),
    y_limits = c(0, 0.6),
    y_label = NULL,
    legend_title = NULL,
    tag = "B",
    id_labels = samples$Category
  )
  
  pAB <- combine_tagged(
    pA, pB,
    ncol = 2
  )
  
  pAB  
}

plotPatientsAndSamples <- function(tag)
{
  p <- make_included_vs_estimated_plot()
  ggsave(file.path(manuscriptPlotDirectory,sprintf("Figure_%s.png",tag))
                       , p, width = 13, height = 8, dpi = 600)
}


ALL_PEK_FIGURES <- function()
{
  
  pekString <- function(n){
    sprintf("%d",n)
  }
  
  pekFigure <-function(n){
    sprintf("Figure_%d.png",n)
  }
  
  
  #start with counter
  n <- 1
  
  #Figures
  #Overview => From Biorender
  file <- paste(processedRoot,"biorender","Clinical evaluation setup.png",sep="/")
  stopifnot(file.exists(file))
  file.copy(from = file, to = paste(manuscriptPlotDirectory,pekFigure(n),sep = "/"),overwrite = TRUE)

  n <- n + 1

  # SIR
  plotSirAllAntibiotics(pekString(n))
  
  n <- n + 1
  
  # Cluster on phenotype
  file <- paste(processedRootRcluster,"cluster_phenotype6.png",sep="/")
  stopifnot(file.exists(file))
  file.copy(from = file, to = paste(manuscriptPlotDirectory,pekFigure(n),sep = "/"),overwrite = TRUE)
  
  n <- n + 1
  
  #Unambiguous correct predictions vs input antibiotic and confidence level. Fractions.
  plotUnambiguousCorrectPredictionsAndErrors(pekString(n))
  
  n <- n + 1
  
  
  # Antibiotic level VME/ME/MCC/F1 rates. Non-conformal. 4 to 8 antibiotics.
  plotAntibioticsMetricsVsNoinputAb(pekString(n))
  n <- n + 1
  
  
  # Figure 6, panel A: HeatMapErrorSir
  file_a <- paste(processedRootRcluster, "predictionerror_sir6.png", sep = "/")
  stopifnot(file.exists(file_a))
  
  # Figure 6, panel B: Abs mean errors
  file_b <- paste(
    manuscriptDirectory,
    "error_vs_genotype",
    "plot",
    "ttest_difference_forest_plot.png",
    sep = "/"
  )
  stopifnot(file.exists(file_b))
  
  output_file <- paste(manuscriptPlotDirectory, pekFigure(n), sep = "/")
  
  make_two_panel_png(
    file_a = file_a,
    file_b = file_b,
    output_file = output_file,
    label_a = "A",
    label_b = "B",
    stack = FALSE
  )
  
# 
#   # HeatMapErrorSir 
#   file <- paste(processedRootRcluster,"predictionerror_sir6.png",sep="/")
#   stopifnot(file.exists(file))
#   file.copy(from = file, to = paste(manuscriptPlotDirectory,pekFigure(n),sep = "/"),overwrite = TRUE)
# 
#   n <- n + 1
#   
#   # Abs mean errors 
#   file <- paste(manuscriptDirectory,"error_vs_genotype","plot", "ttest_difference_forest_plot.png",sep="/")
#   stopifnot(file.exists(file))
#   file.copy(from = file, to = paste(manuscriptPlotDirectory,suppFigure(n),sep = "/"),overwrite = TRUE)
#   n <- n + 1
#   
  

  suppString <- function(n){
    sprintf("S%d",n)
  }
  
  suppFigure <-function(n){
    sprintf("Figure_S%d.png",n)
  }
  
  n <- 1
  

  #Patients and samples
  plotPatientsAndSamples(suppString(n))
  n <- n + 1
  
  
  #Antibiotic level: Unambiguous correct predictions vs input antibiotic. Fractions.
  plotUnambiguousCorrectPredictionsAndErrorsVersusIndexAntibioticsLevel(suppString(n),"STD")
  n <- n + 1
  
  #Antibiotic level: Unambiguous correct predictions vs signi. Fractions.
  plotUnambiguousCorrectPredictionsAndErrorsSignificanceLevelAntibioticsLevel(suppString(n))
  n <- n + 1

  # Antibiotic level VME/ME/MCC/F1 rates. Conformal, 2,5, 5, and 10%.
  plotAntibioticsMetricsVsSignificanceLevel(suppString(n))
  n <- n + 1
  
  
  #Performance antibiotics Group - compare with pek. Non-conformal.
  plotAntibioticsMetricsVsNoinputAb_Abgroups(suppString(n))
  n <- n + 1
  #Performance antibiotics Group - compare with pek. significance level.
  plotAntibioticsMetricsVsSignificanceLevel_Abgroups(suppString(n))
  n <- n + 1
  
  #Antibiotic level VME/ME/MCC/F1 rates. Non-conformal. 4 to 8 antibiotics. One per group.
  plotAntibioticsMetricsVsNoinputAbOnePerAbGroup(suppString(n))
  n <- n + 1
  
  #Fig. S2 - Antibiotic level VME/ME/MCC/F1 rates. Non-conformal. Best choice. 
  plotAntibioticsMetricsVsSignificanceLevelSixBest(suppString(n))
  n <- n + 1
  
  
  #Fig. S5 - Antibiotic level VME/ME/MCC/F1 rates. Non-conformal. Modes A B C.
  plotAntibioticsMetricsVsMode(suppString(n))
  n <- n + 1

  ################### Genotype ####################
  
  file <- paste(processedRootRcluster,"genotype_aggregated_sir6.png",sep="/")
  stopifnot(file.exists(file))
  file.copy(from = file, to = paste(manuscriptPlotDirectory,suppFigure(n),sep = "/"),overwrite = TRUE)
  n <- n + 1
  
  file <- paste(processedRootRcluster,"genotype_sir6.png",sep="/")
  stopifnot(file.exists(file))
  file.copy(from = file, to = paste(manuscriptPlotDirectory,suppFigure(n),sep = "/"),overwrite = TRUE)
  n <- n + 1

  ################### prediction error vs genotype ####################
  # file <- paste(processedRootRcluster,"predictionerror_genotype_aggregated_sir6.png",sep="/")
  # stopifnot(file.exists(file))
  # file.copy(from = file, to = paste(manuscriptPlotDirectory,suppFigure(n),sep = "/"),overwrite = TRUE)
  # n <- n + 1
  # 
  # file <- paste(processedRootRcluster,"predictionerror_genotype_sir6.png",sep="/")
  # stopifnot(file.exists(file))
  # file.copy(from = file, to = paste(manuscriptPlotDirectory,suppFigure(n),sep = "/"),overwrite = TRUE)
  # n <- n + 1

  # Supplementary Figure: prediction errors ordered by error + mean absolute error by property
  file_a <- paste(processedRootRcluster, "sir_heatmap_ordered_by_error6.png", sep = "/")
  stopifnot(file.exists(file_a))
  
  file_b <- paste(
    manuscriptDirectory,
    "error_vs_genotype",
    "plot",
    "boxplot_mean_abs_error_by_property.png",
    sep = "/"
  )
  stopifnot(file.exists(file_b))
  
  output_file <- paste(manuscriptPlotDirectory, suppFigure(n), sep = "/")
  
  make_two_panel_png(
    file_a = file_a,
    file_b = file_b,
    output_file = output_file,
    label_a = "A",
    label_b = "B",
    stack = FALSE
  )
  
  n <- n + 1
  

  # file <- paste(manuscriptDirectory,"error_vs_genotype","plot", "ttest_difference_forest_plot.png",sep="/")
  # stopifnot(file.exists(file))
  # file.copy(from = file, to = paste(manuscriptPlotDirectory,suppFigure(n),sep = "/"),overwrite = TRUE)
  # n <- n + 1
  
    
  ########## Extra figure that might be reintroduces
  xString <- function(n){
    sprintf("X%d",n)
  }
  
  xFigure <-function(n){
    sprintf("Figure_X%d.png",n)
  }
  
  n <- 1
  
  
  #Fig. QC
  file <- file.path(processedRootR,"qc","qc_bands.png")
  stopifnot(file.exists(file))
  file.copy(from = file, to = paste(manuscriptPlotDirectory,xFigure(n),sep = "/"),overwrite = TRUE)
  n <- n + 1
  
  #Fig. Compare to PEK. ME and VME Conformal 10%,5%.
  plotMetricsVsSet(xString(n))
  n <- n + 1
  #Fig. Compare to PEK. F1 Conformal 10%,5%.
  plotF1AbgroupsVsSet(xString(n))

  
  
}



ALL <- function()
{
  ENSURE_READS()
  ALL_PEK_FIGURES()
  
}
