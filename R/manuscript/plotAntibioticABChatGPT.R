plot_fig7A <- function(df,
                       id_col = antibiotic,
                       measure_cols = c("4","5","6","7","8"),
                       tag = "A",
                       out_dir = NULL,
                       file_stub = "figure_7A",
                       export = TRUE) {
  # Output dir (defaults to manuscriptPlotDirectory if available)
  if (is.null(out_dir)) {
    od <- get0("manuscriptPlotDirectory", inherits = TRUE)
    out_dir <- if (!is.null(od)) od else "manuscriptPlotDirectory"
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  id_name <- as.character(rlang::ensym(id_col))
  
  # Wide → long
  long <- tidyr::pivot_longer(df,
                              cols = tidyselect::all_of(measure_cols),
                              names_to = "n_inputs",
                              values_to = "share")
  
  # Keep given order of antibiotics and n
  long[[id_name]] <- factor(long[[id_name]], levels = unique(long[[id_name]]))
  long$n_inputs   <- factor(long$n_inputs, levels = measure_cols)
  
  # Greens like 5A (light → dark)
  greens <- grDevices::colorRampPalette(c("#D9EAD3", "#6AA84F", "#38761D"))(length(measure_cols))
  names(greens) <- levels(long$n_inputs)
  
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = !!rlang::sym(id_name), y = share, fill = n_inputs)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78), width = 0.66) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    ggplot2::scale_fill_manual(values = greens, name = "Number of input antibiotics") +
    ggplot2::labs(x = NULL, y = "Unambiguous correct predictions", tag = tag) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 10),
      legend.text  = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      plot.tag = ggplot2::element_text(face = "bold", size = 14),
      plot.tag.position = c(0.01, 0.99)
    )
  
  if (isTRUE(export)) {
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".png")),
                    p, width = 8.5, height = 5, dpi = 300, bg = "white")
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".svg")),
                    p, width = 8.5, height = 5, bg = "white")
    file_xlsx <- file.path(out_dir, paste0(file_stub, ".xlsx"))
    write.xlsx(x = df, file = file_xlsx)
  }
  p
}


plot_fig7B <- function(metric_list,
                       pick = "6",                 # "4", "6", or "8" (or numeric index)
                       tag = "B",
                       out_dir = NULL,
                       file_stub = "figure_7B",
                       export = TRUE) {
  # Output dir (uses manuscriptPlotDirectory if available)
  if (is.null(out_dir)) {
    od <- get0("manuscriptPlotDirectory", inherits = TRUE)
    out_dir <- if (!is.null(od)) od else "manuscriptPlotDirectory"
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Select and reshape data
  df <- if (is.character(pick)) metric_list[[pick]] else metric_list[[pick]]
  long <- tidyr::pivot_longer(df, cols = c(correct, ME, VME),
                              names_to = "metric", values_to = "share")
  
  # Order & factors
  long$antibiotic <- factor(long$antibiotic, levels = df$antibiotic) # keep given order
  long$metric     <- factor(long$metric, levels = c("correct","ME","VME"))  # bottom→top
  
  # Colors (same idea as 5B; adjust if you want your vivid ramps instead)
  metric_colors <- c(correct = "#70AD47",  # green
                     ME      = "#FFFF00",  # yellow
                     VME     = "#FF0000")  # red
  
  # Plot
  p <- ggplot2::ggplot(long, ggplot2::aes(x = antibiotic, y = share, fill = metric)) +
    ggplot2::geom_col(width = 0.75, position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_fill_manual(values = metric_colors, name = NULL) +
    ggplot2::labs(x = NULL, y = "Percentages of prediction results", tag = tag) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      plot.tag = ggplot2::element_text(face = "bold", size = 14),
      plot.tag.position = c(0.01, 0.99)
    )
  
  # Export
  if (isTRUE(export)) {
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".png")),
                    p, width = 8.8, height = 5.2, dpi = 300, bg = "white")
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".svg")),
                    p, width = 8.8, height = 5.2, bg = "white")
    file_xlsx <- file.path(out_dir, paste0(file_stub, ".xlsx"))
    write.xlsx(x = df, file = file_xlsx)
  }
  p
}
