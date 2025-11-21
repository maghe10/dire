plot_siglevel_grouped <- function(df,
                                  id_col = significanceLevel,
                                  measure_cols = NULL,
                                  names_to = "n_inputs",
                                  values_to = "share",
                                  tag = "A",
                                  out_dir = NULL,
                                  file_stub = "figure_5A",
                                  export = TRUE) {
  # --- Resolve output directory (uses manuscriptPlotDirectory if it exists)
  if (is.null(out_dir)) {
    od <- get0("manuscriptPlotDirectory", inherits = TRUE)
    out_dir <- if (!is.null(od)) od else "manuscriptPlotDirectory"
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # --- Input normalization: wide or long
  nm <- names(df)
  id_name <- as.character(rlang::ensym(id_col))
  is_long <- all(c(id_name, names_to, values_to) %in% nm)
  
  if (!is_long) {
    stopifnot(id_name %in% nm)
    if (is.null(measure_cols)) {
      measure_cols <- setdiff(nm, id_name)
    }
    if (length(measure_cols) == 0) stop("No measure columns found.")
    
    # Clean decimals & pivot
    df[measure_cols] <- lapply(df[measure_cols], function(x) {
      x <- as.character(x)
      x <- gsub(",", ".", x, fixed = TRUE)
      suppressWarnings(as.numeric(x))
    })
    
    long <- tidyr::pivot_longer(df,
                                cols = tidyselect::all_of(measure_cols),
                                names_to = names_to,
                                values_to = values_to)
  } else {
    long <- df
  }
  
  # --- Preserve order
  long[[id_name]] <- factor(long[[id_name]], levels = unique(long[[id_name]]))
  long[[names_to]] <- factor(long[[names_to]], levels = unique(long[[names_to]]))
  
  # --- Dynamic colors
  group_levels <- levels(long[[names_to]])
  greens <- grDevices::colorRampPalette(c("#C3D8BB", "#70AD47", "#538233"))(length(group_levels))
  names(greens) <- group_levels
  
  # --- Plot (no aes_string!)
  p <- ggplot2::ggplot(long, ggplot2::aes(x = !!rlang::sym(id_name),
                                          y = !!rlang::sym(values_to),
                                          fill = !!rlang::sym(names_to))) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.68) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::scale_fill_manual(values = greens, name = "Number of input antibiotics") +
    ggplot2::labs(x = NULL, y = "Unambiguous correct predictions", tag = tag) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.title = ggplot2::element_text(size = 10),
                   legend.text  = ggplot2::element_text(size = 10),
                   plot.tag = ggplot2::element_text(face = "bold", size = 14),
                   plot.tag.position = c(0.01, 0.99))
  
  # --- Export
  if (isTRUE(export)) {
    
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".png")),
                    p, width = 6.5, height = 4, dpi = 300, bg = "white")
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".svg")),
                    p, width = 6.5, height = 4, bg = "white")
    file_xlsx <- file.path(out_dir, paste0(file_stub, ".xlsx"))
    write.xlsx(x = df, file = file_xlsx)
  }
  p
}


plot_prediction_metrics_stacked <- function(metric_list,
                                            pick = "6",                 # e.g. "4", "6", "8" or numeric index
                                            tag = "B",
                                            out_dir = NULL,
                                            file_stub = "figure_5B",
                                            export = TRUE) {
  # --- Resolve output directory (uses manuscriptPlotDirectory if available)
  if (is.null(out_dir)) {
    od <- get0("manuscriptPlotDirectory", inherits = TRUE)
    out_dir <- if (!is.null(od)) od else "manuscriptPlotDirectory"
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # --- Select dataset
  df_wide <- if (is.character(pick)) metric_list[[pick]] else metric_list[[pick]]
  
  # --- Pivot to long
  long <- tidyr::pivot_longer(
    df_wide,
    cols = -metric,
    names_to = "significanceLevel",
    values_to = "share"
  )
  
  # --- Preserve and order factors
  # order: "correct" (bottom), "ME" (middle), "VME" (top)
  long$metric <- factor(long$metric, levels = c("VME", "ME", "correct"))
  long$significanceLevel <- factor(
    long$significanceLevel,
    levels = c("non-conformal", "10%", "5%", "2.5%")
  )
  
  # --- Colors
  metric_colors <- c(
    "correct" = "#70AD47",   # green
    "ME"      = "#FFFF00",   # yellow
    "VME"     = "#FF0000"    # red
  )
  
  # --- Plot
  p <- ggplot2::ggplot(long,
                       ggplot2::aes(x = significanceLevel,
                                    y = share,
                                    fill = metric)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                limits = c(0, 1)) +
    ggplot2::scale_fill_manual(values = metric_colors, name = NULL) +
    ggplot2::labs(x = NULL,
                  y = "Percentages of prediction results",
                  tag = tag) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 14),
      plot.tag.position = c(0.01, 0.99)
    )
  
  # --- Export
  if (isTRUE(export)) {
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".png")),
                    p, width = 6.5, height = 4, dpi = 300, bg = "white")
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".svg")),
                    p, width = 6.5, height = 4, bg = "white")
    file_xlsx <- file.path(out_dir, paste0(file_stub, ".xlsx"))
    write.xlsx(x = df_wide, file = file_xlsx)
  }
  
  p
}