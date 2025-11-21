## metric_plot_helpers.R
## Helpers to:
##  - optionally collapse antibiotics -> ab_group
##  - compute metrics from correctR/correctS/falseR/falseS
##  - produce bar plots & stacked error plots (S4/S5/S6-style)
##  - support dual-line titles: bold main + wrapped subcaption
##  - standardized export dimensions (default width=8, height=4, dpi=300)
##  - centralized export helper (PNG + XLSX + CSV only)

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(rlang)
library(scales)
library(patchwork)
library(ggtext)
library(stringr)
## NOTE: write.xlsx() comes from e.g. openxlsx or xlsx package.

# ---- Utility: safe_div -------------------------------------------------------

safe_div <- function(num, den) {
  ifelse(is.na(den) | den == 0, NA_real_, num / den)
}

# ---- Title helper (bold main + wrapped subcaption) --------------------------

make_fig_title <- function(main = NULL, sub = NULL, width = 85) {
  if (is.null(main) && is.null(sub)) return(NULL)
  
  if (!is.null(main) && nchar(main) > 0) {
    main_wrapped <- stringr::str_wrap(main, width = as.integer(width/1.1))
    main_html <- gsub("\n", "<br>", main_wrapped, fixed = TRUE)
    main_html <- paste0("<b>", main_html, "</b>")
    
  } else {
    main_html <- NULL
  }
  
  sub_html <- NULL
  if (!is.null(sub) && nchar(sub) > 0) {
    sub_wrapped <- stringr::str_wrap(sub, width = width)
    sub_html <- gsub("\n", "<br>", sub_wrapped, fixed = TRUE)
  }
  
  if (!is.null(main_html) && !is.null(sub_html)) {
    paste0(main_html, "<br>", sub_html)
  } else if (!is.null(main_html)) {
    main_html
  } else {
    sub_html
  }
}

# ---- Centralized export helper (PNG + XLSX + CSV) ---------------------------

export_figure <- function(plot,
                          data,
                          out_dir,
                          file_stub,
                          width = 8,
                          height = 4,
                          dpi = 300) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  file_png  <- file.path(out_dir, paste0(file_stub, ".png"))
  file_xlsx <- file.path(out_dir, paste0(file_stub, ".xlsx"))
  file_csv  <- file.path(out_dir, paste0(file_stub, ".csv"))
  # # ---- Add a thin frame INSIDE the plot ----
  # plot <- plot +
  #   theme(
  #     panel.background = element_rect(fill = "white", colour = NA),
  #     panel.border     = element_rect(color = "black", linewidth = 0.4, fill = NA)
  #   )
  
  
  # PNG
  ggplot2::ggsave(
    filename = file_png,
    plot     = plot,
    width    = width,
    height   = height,
    units    = "in",
    dpi      = dpi,
    bg       = "white"
  )
  
  # Table exports (if data is supplied)
  if (!is.null(data)) {
    write.xlsx(x = data, file = file_xlsx)
    write.csv2(data, file = file_csv, row.names = FALSE)
  }
}

# ---- Antibiotic group definition ---------------------------------------------

make_ab_groups <- function() {
  tibble::tribble(
    ~antibiotic, ~ab_group,
    # Penicillins / β-lactam + inhibitor
    "AMP", "Penicillins",
    "AMC", "Penicillins",
    "PIP", "Penicillins",
    "TZP", "Penicillins",
    
    # Cephalosporins
    "CAZ", "Cephalosporins",
    "CRO", "Cephalosporins",
    "CTX", "Cephalosporins",
    "FEP", "Cephalosporins",
    
    # Fluoroquinolones
    "CIP", "Fluoroquinolones",
    "OFX", "Fluoroquinolones",
    "LVX", "Fluoroquinolones",
    "MFX", "Fluoroquinolones",
    
    # Aminoglycosides
    "GEN", "Aminoglycosides",
    "TOB", "Aminoglycosides"
  )
}

# ---- Collapse counts to antibiotic groups ------------------------------------

collapse_to_ab_groups <- function(
    df,
    ab_groups  = make_ab_groups(),
    group_vars = c("noinputab", "significanceLevel", "mode")
) {
  df %>%
    left_join(ab_groups, by = "antibiotic") %>%
    mutate(ab_group = if_else(is.na(ab_group), antibiotic, ab_group)) %>%
    {
      gv <- intersect(c(group_vars, "ab_group"), names(.))
      cols_to_sum <- setdiff(names(.), c(gv, "antibiotic"))
      
      group_by(., across(all_of(gv))) %>%
        summarise(
          across(all_of(cols_to_sum), ~ sum(.x, na.rm = TRUE)),
          .groups = "drop"
        )
    } %>%
    rename(antibiotic = ab_group)
}

# ---- Compute metrics from count data (TP/FP/TN/FN) ---------------------------

compute_metrics_long_counts <- function(
    df,
    group_vars = c("noinputab", "antibiotic", "significanceLevel", "mode"),
    mapping    = c(TP = "correctR",
                   FP = "falseR",
                   TN = "correctS",
                   FN = "falseS")
) {
  missing_cols <- setdiff(unname(mapping), names(df))
  if (length(missing_cols) > 0) {
    stop(
      "These count columns from 'mapping' are missing in df: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  group_vars <- intersect(group_vars, names(df))
  
  summarised <- df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      across(all_of(unname(mapping)), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  TP_col <- mapping["TP"]
  FP_col <- mapping["FP"]
  TN_col <- mapping["TN"]
  FN_col <- mapping["FN"]
  
  metrics_wide <- summarised %>%
    mutate(
      precision = safe_div(.data[[TP_col]], .data[[TP_col]] + .data[[FP_col]]),
      recall    = safe_div(.data[[TP_col]], .data[[TP_col]] + .data[[FN_col]]),
      F1        = safe_div(
        2 * .data[[TP_col]],
        2 * .data[[TP_col]] + .data[[FP_col]] + .data[[FN_col]]
      ),
      ME        = safe_div(.data[[FP_col]], .data[[TN_col]] + .data[[FP_col]]),
      VME       = safe_div(.data[[FN_col]], .data[[TP_col]] + .data[[FN_col]]),
      correct   = safe_div(
        .data[[TP_col]] + .data[[TN_col]],
        .data[[TP_col]] + .data[[FP_col]] + .data[[TN_col]] + .data[[FN_col]]
      )
    )
  
  metrics_wide %>%
    pivot_longer(
      cols      = c(precision, recall, F1, ME, VME, correct),
      names_to  = "metric",
      values_to = "value"
    )
}

# ---- Helper: metric_wide_for_plot (id x column_var) --------------------------

metric_wide_for_plot <- function(
    metrics_long,
    metric,
    id_col,
    column_var,
    filter_expr
) {
  id_col     <- enquo(id_col)
  column_var <- enquo(column_var)
  
  df <- metrics_long %>%
    filter(.data$metric == !!metric)
  
  if (!missing(filter_expr)) {
    df <- df %>% filter({{ filter_expr }})
  }
  
  df %>%
    select(!!id_col, !!column_var, value) %>%
    distinct() %>%
    pivot_wider(
      names_from  = !!column_var,
      values_from = value
    ) %>%
    arrange(!!id_col)
}

# ---- Plot metric bars (S4-style, grouped by column_var) ----------------------

plot_metric_bars <- function(
    metrics_long,
    metric,
    id_col,
    column_var,
    tag            = "A",
    ylab           = NULL,
    filter_expr,
    out_dir        = NULL,
    file_stub      = NULL,
    export         = FALSE,
    fig_title_main = NULL,
    fig_title_sub  = NULL,
    width          = 8,
    height         = 4,
    dpi            = 300
) {
  if (is.null(ylab)) {
    ylab <- metric
  }
  
  if (is.null(out_dir)) {
    od <- get0("manuscriptPlotDirectory", inherits = TRUE)
    out_dir <- if (!is.null(od)) od else "manuscriptPlotDirectory"
  }
  
  id_quo     <- enquo(id_col)
  column_quo <- enquo(column_var)
  
  df_wide <- metric_wide_for_plot(
    metrics_long  = metrics_long,
    metric        = metric,
    id_col        = !!id_quo,
    column_var    = !!column_quo,
    filter_expr   = {{ filter_expr }}
  )
  
  id_name <- as_name(id_quo)
  measure_cols <- setdiff(names(df_wide), id_name)
  
  long <- df_wide %>%
    pivot_longer(
      cols      = tidyselect::all_of(measure_cols),
      names_to  = "group_var",
      values_to = "value"
    )
  
  long[[id_name]] <- factor(long[[id_name]], levels = unique(long[[id_name]]))
  long$group_var  <- factor(long$group_var, levels = measure_cols)
  
  pal <- grDevices::colorRampPalette(
    c("#D9EAD3", "#6AA84F", "#38761D")
  )(length(measure_cols))
  names(pal) <- levels(long$group_var)
  
  p <- ggplot(
    long,
    aes(x = !!id_quo, y = value, fill = group_var)
  ) +
    geom_col(
      position = position_dodge(width = 0.78),
      width    = 0.66,
      na.rm    = TRUE
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(
      values = pal,
      name   = as_name(column_quo)
    ) +
    labs(
      x     = NULL,
      y     = ylab,
      tag   = tag,
      title = make_fig_title(fig_title_main, fig_title_sub)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom",
      legend.title      = element_text(size = 10),
      legend.text       = element_text(size = 10),
      axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.tag          = element_text(face = "bold", size = 14),
      plot.tag.position = c(0.01, 0.99),
      plot.title        = ggtext::element_markdown(hjust = 0, size = 12)
      
    )
  
  if (isTRUE(export)) {
    if (is.null(file_stub)) {
      file_stub <- paste0("figure_", tag, "_", metric)
    }
    export_figure(
      plot     = p,
      data     = df_wide,
      out_dir  = out_dir,
      file_stub = file_stub,
      width    = width,
      height   = height,
      dpi      = dpi
    )
  }
  
  p
}

# ---- Multi-panel metric figure (S4-style, patchwork) -------------------------

plot_metric_panels <- function(
    metrics_long,
    metric_defs,                 # tibble: metric, tag, ylab
    id_col,
    column_var,
    filter_expr,
    fig_stub       = "figure",
    out_dir        = NULL,
    nrow           = NULL,
    width          = 8,
    height         = 4,
    dpi            = 300,
    export         = TRUE,
    fig_title_main = NULL,
    fig_title_sub  = NULL
) {
  if (is.null(out_dir)) {
    od <- get0("manuscriptPlotDirectory", inherits = TRUE)
    out_dir <- if (!is.null(od)) od else "manuscriptPlotDirectory"
  }
  
  if (is.null(nrow)) {
    nrow <- nrow(metric_defs)
  }
  
  plots <- pmap(
    metric_defs,
    function(metric, tag, ylab) {
      plot_metric_bars(
        metrics_long   = metrics_long,
        metric         = metric,
        id_col         = {{ id_col }},
        column_var     = {{ column_var }},
        tag            = tag,
        ylab           = ylab,
        filter_expr    = {{ filter_expr }},
        out_dir        = out_dir,
        file_stub      = paste0(fig_stub, "_", tag),
        export         = FALSE,
        fig_title_main = NULL,
        fig_title_sub  = NULL,
        width          = width,
        height         = height,
        dpi            = dpi
      )
    }
  )
  
  combined <- patchwork::wrap_plots(plots, nrow = nrow)
  
  if (!is.null(fig_title_main) || !is.null(fig_title_sub)) {
    combined <- combined +
      patchwork::plot_annotation(
        title = make_fig_title(fig_title_main, fig_title_sub),
        theme = ggplot2::theme(
          plot.title = ggtext::element_markdown(
            hjust = 0,
            size  = 12
          )
        )
      )
  }
  
  if (isTRUE(export)) {
    metrics_wide <- metrics_long %>%
      tidyr::pivot_wider(names_from = metric, values_from = value)
    
    export_figure(
      plot      = combined,
      data      = metrics_wide,
      out_dir   = out_dir,
      file_stub = fig_stub,
      width     = width,
      height    = height,
      dpi       = dpi
    )
  }
  
  combined
}

# ---- Simple metric bar plot (sorted, single panel) ---------------------------

plot_metric_bars_simple <- function(
    metrics_long,
    metric,
    id_col,
    tag            = "",
    ylab           = NULL,
    export         = FALSE,
    show_x_labels  = FALSE,
    fig_title_main = NULL,
    fig_title_sub  = NULL,
    width          = 8,
    height         = 4,
    dpi            = 300
) {
  id_quo <- enquo(id_col)
  
  df_plot <- metrics_long %>%
    dplyr::filter(.data$metric == !!metric) %>%
    dplyr::filter(!is.na(value), is.finite(value)) %>%
    dplyr::arrange(dplyr::desc(value)) %>%
    dplyr::mutate(!!id_quo := factor(!!id_quo, levels = !!id_quo))
  
  if (is.null(ylab)) ylab <- metric
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = !!id_quo, y = value)) +
    ggplot2::geom_col() +
    ggplot2::labs(
      x     = if (show_x_labels) rlang::as_label(id_quo) else NULL,
      y     = ylab,
      title = make_fig_title(fig_title_main, fig_title_sub)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggtext::element_markdown(
        hjust = 0,
        size  = 12
      )
    )
  
  if (!show_x_labels) {
    p <- p +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  } else {
    p <- p +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  }
  
  if (nzchar(tag)) {
    p <- p +
      ggplot2::annotate(
        "text",
        x = -Inf, y = Inf,
        label = tag,
        fontface = "bold",
        size = 6,
        hjust = -0.1,
        vjust = 1.3
      ) +
      ggplot2::coord_cartesian(clip = "off")
  }
  
  if (export) {
    fname_base <- paste0(
      "metric_", metric,
      if (nzchar(tag)) paste0("_", tag) else ""
    )
    
    export_figure(
      plot      = p,
      data      = df_plot,
      out_dir   = getwd(),
      file_stub = fname_base,
      width     = width,
      height    = height,
      dpi       = dpi
    )
  }
  
  p
}

# ---- Simple multi-panel metric figure (sorted bars) --------------------------

plot_metrics_panels_simple <- function(
    metrics_long,
    metric_defs,
    id_col,
    fig_stub       = "figure",
    out_dir        = NULL,
    nrow           = NULL,
    width          = 8,
    height         = 4,
    dpi            = 300,
    export         = TRUE,
    show_x_labels_last = TRUE,
    fig_title_main = NULL,
    fig_title_sub  = NULL
) {
  id_quo <- rlang::enquo(id_col)
  
  if (is.null(out_dir)) {
    od <- get0("manuscriptPlotDirectory", inherits = TRUE)
    out_dir <- if (!is.null(od)) od else "manuscriptPlotDirectory"
  }
  
  n_panels <- nrow(metric_defs)
  if (is.null(nrow)) {
    nrow <- n_panels
  }
  
  plots <- vector("list", n_panels)
  
  for (i in seq_len(n_panels)) {
    metric <- metric_defs$metric[i]
    tag    <- metric_defs$tag[i]
    ylab   <- metric_defs$ylab[i]
    
    show_x <- show_x_labels_last && (i == n_panels)
    
    plots[[i]] <- plot_metric_bars_simple(
      metrics_long   = metrics_long,
      metric         = metric,
      id_col         = !!id_quo,
      tag            = tag,
      ylab           = ylab,
      export         = FALSE,
      show_x_labels  = show_x,
      fig_title_main = NULL,
      fig_title_sub  = NULL,
      width          = width,
      height         = height,
      dpi            = dpi
    )
  }
  
  combined <- patchwork::wrap_plots(plots, nrow = nrow)
  
  if (!is.null(fig_title_main) || !is.null(fig_title_sub)) {
    combined <- combined +
      patchwork::plot_annotation(
        title = make_fig_title(fig_title_main, fig_title_sub),
        theme = ggplot2::theme(
          plot.title = ggtext::element_markdown(
            hjust = 0,
            size  = 13
          )
        )
      )
  }
  
  if (isTRUE(export)) {
    metrics_wide <- metrics_long %>%
      tidyr::pivot_wider(names_from = metric, values_from = value)
    
    export_figure(
      plot      = combined,
      data      = metrics_wide,
      out_dir   = out_dir,
      file_stub = fig_stub,
      width     = width,
      height    = height,
      dpi       = dpi
    )
  }
  
  combined
}

# ---- Stacked errors by antibiotic & column_var (correct/ME/VME) -------------

plot_stacked_errors <- function(metrics_long,
                                id_col,
                                column_var,
                                value_col,
                                filter_expr    = TRUE,
                                ylab           = "Percentage of prediction results",
                                export         = FALSE,
                                out_dir        = NULL,
                                file_stub      = "figure_S4A_stacked",
                                width          = 8,
                                height         = 4,
                                dpi            = 300,
                                fig_title_main = NULL,
                                fig_title_sub  = NULL) {
  
  id_col      <- rlang::ensym(id_col)
  column_var  <- rlang::ensym(column_var)
  value_col   <- rlang::ensym(value_col)
  filter_expr <- rlang::enquo(filter_expr)
  
  dat <- metrics_long %>%
    dplyr::filter(
      !!filter_expr,
      metric %in% c("correct", "ME", "VME")
    ) %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("correct", "ME", "VME"))
    ) %>%
    dplyr::filter(!is.na(!!value_col))
  
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = !!column_var, y = !!value_col, fill = metric)
  ) +
    ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::facet_grid(cols = vars(!!id_col), switch = "x") +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1),
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        correct = "#4daf4a",
        ME      = "yellow",
        VME     = "red"
      ),
      name = NULL
    ) +
    ggplot2::labs(
      x     = NULL,
      y     = ylab,
      title = make_fig_title(fig_title_main, fig_title_sub)
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      strip.text.x     = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 0, hjust = 0.5, size = 6),
      plot.title       = ggtext::element_markdown(
        hjust = 0, size = 11
      )
    )
  
  if (export) {
    if (is.null(out_dir)) {
      od <- get0("manuscriptPlotDirectory", inherits = TRUE)
      out_dir <- if (!is.null(od)) od else "manuscriptPlotDirectory"
    }
    
    metrics_wide <- metrics_long %>%
      tidyr::pivot_wider(names_from = metric, values_from = value)
    
    export_figure(
      plot      = p,
      data      = metrics_wide,
      out_dir   = out_dir,
      file_stub = file_stub,
      width     = width,
      height    = height,
      dpi       = dpi
    )
  }
  
  p
}

# ---- Stacked errors by id_col only (single x-variable) ----------------------

plot_stacked_errors_single <- function(metrics_long,
                                       id_col,
                                       value_col,
                                       filter_expr    = TRUE,
                                       ylab           = "Percentage of prediction results",
                                       export         = FALSE,
                                       out_dir        = NULL,
                                       file_stub      = "figure_stacked_single",
                                       width          = 8,
                                       height         = 4,
                                       dpi            = 300,
                                       fig_title_main = NULL,
                                       fig_title_sub  = NULL) {
  
  id_col      <- rlang::ensym(id_col)
  value_col   <- rlang::ensym(value_col)
  filter_expr <- rlang::enquo(filter_expr)
  
  dat <- metrics_long %>%
    dplyr::filter(
      !!filter_expr,
      metric %in% c("correct", "ME", "VME")
    ) %>%
    dplyr::mutate(
      metric = factor(metric, levels = c("correct", "ME", "VME"))
    ) %>%
    dplyr::filter(!is.na(!!value_col))
  
  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(x = !!id_col, y = !!value_col, fill = metric)
  ) +
    ggplot2::geom_col(position = ggplot2::position_stack(reverse = TRUE)) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1),
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        correct = "#4daf4a",
        ME      = "yellow",
        VME     = "red"
      ),
      name = NULL
    ) +
    ggplot2::labs(
      x     = NULL,
      y     = ylab,
      title = make_fig_title(fig_title_main, fig_title_sub)
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(),
      plot.title       = ggtext::element_markdown(
        hjust = 0, size = 11
      )
    )
  
  if (export) {
    if (is.null(out_dir)) {
      od <- get0("manuscriptPlotDirectory", inherits = TRUE)
      out_dir <- if (!is.null(od)) od else "manuscriptPlotDirectory"
    }
    
    metrics_wide <- metrics_long %>%
      tidyr::pivot_wider(names_from = metric, values_from = value)
    
    export_figure(
      plot      = p,
      data      = metrics_wide,
      out_dir   = out_dir,
      file_stub = file_stub,
      width     = width,
      height    = height,
      dpi       = dpi
    )
  }
  
  p
}


# ---- Call plotting via machinery with callback ----------------------

plot_fig_via_callback <- function(df,
                       id_col = antibiotic,
                       measure_cols = NULL,
                       tag = "A",
                       ylab = "F1 score",
                       out_dir = NULL,
                       file_stub = paste("figure_S4", tag),
                       export = TRUE) {
  # Output dir (defaults to manuscriptPlotDirectory if available)
  if (is.null(out_dir)) {
    od <- get0("manuscriptPlotDirectory", inherits = TRUE)
    out_dir <- if (!is.null(od)) od else "manuscriptPlotDirectory"
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  id_name <- as.character(rlang::ensym(id_col))
  
  # Auto-detect measure columns if not supplied:
  # everything except the id column
  if (is.null(measure_cols)) {
    measure_cols <- setdiff(names(df), id_name)
  }
  
  # Wide → long
  long <- tidyr::pivot_longer(
    df,
    cols      = tidyselect::all_of(measure_cols),
    names_to  = "n_inputs",
    values_to = "f1"
  )
  
  # Keep the given order
  long[[id_name]] <- factor(long[[id_name]], levels = unique(long[[id_name]]))
  long$n_inputs   <- factor(long$n_inputs, levels = measure_cols)
  
  # Greens (light → dark), one per n_inputs
  greens <- grDevices::colorRampPalette(c("#D9EAD3", "#6AA84F", "#38761D"))(length(measure_cols))
  names(greens) <- levels(long$n_inputs)
  
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = !!rlang::sym(id_name), y = f1, fill = n_inputs)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.78),
      width = 0.66,
      na.rm = TRUE
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::scale_fill_manual(
      values = greens,
      name   = "Number of input antibiotics"
    ) +
    ggplot2::labs(x = NULL, y = ylab, tag = tag) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      legend.position   = "bottom",  # <-- FIXED HERE
      legend.title      = ggplot2::element_text(size = 10),
      legend.text       = ggplot2::element_text(size = 10),
      axis.text.x       = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      plot.tag          = ggplot2::element_text(face = "bold", size = 14),
      plot.tag.position = c(0.01, 0.99)
    )
  
  if (isTRUE(export)) {
    file_png  <- file.path(out_dir, paste0(file_stub, ".png"))
    file_svg  <- file.path(out_dir, paste0(file_stub, ".svg"))
    file_xlsx <- file.path(out_dir, paste0(file_stub, ".xlsx"))
    
    ggplot2::ggsave(file_png, p, width = 8.5, height = 5, dpi = 300, bg = "white")
    ggplot2::ggsave(file_svg, p, width = 8.5, height = 5, bg = "white")
    
    # underlying wide frame used for the plot
    write.xlsx(x = df, file = file_xlsx)
  }
  
  p
}

# ---- Plotting function with metric defs and callback function -----

plot_figure_panels_with_callback <- function(
    selectedMode,
    metric_defs,
    fetch_fun  = fetchfractionMetricAntibioticsVsIndex,
    fig_stub   = "figure",
    out_dir    = NULL,
    nrow       = NULL,
    width      = 13,
    height     = 5.2,
    dpi        = 600,
    export     = TRUE
) {
  library(dplyr)
  library(purrr)
  library(patchwork)
  library(ggplot2)
  
  # Output directory (defaults to manuscriptPlotDirectory if present)
  if (is.null(out_dir)) {
    od <- get0("manuscriptPlotDirectory", inherits = TRUE)
    out_dir <- if (!is.null(od)) od else "manuscriptPlotDirectory"
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Default: one row per metric (vertical stacking)
  if (is.null(nrow)) {
    nrow <- nrow(metric_defs)
  }
  
  # 1) Compute all metrics in long format
  metrics_long <- computeIndexMetricsLong(selectedMode, fetch_fun = fetch_fun)
  metrics_wide <- metrics_long %>%
    tidyr::pivot_wider(names_from = metric, values_from = value)
  
  # 2) Build a single panel for each metric definition
  plots <- purrr::pmap(
    metric_defs,
    function(metric, tag, ylab) {
      df_wide <- makeWide(metrics_long, metric)
      plot_fig_via_callback(df_wide, tag = tag, ylab = ylab, export = FALSE)
    }
  )
  
  # 3) Combine into a multi-panel figure
  combined <- patchwork::wrap_plots(plots, nrow = nrow)
  
  # 4) Save final combined figure + underlying data
  if (isTRUE(export)) {
    file_png  <- file.path(out_dir, paste0(fig_stub, ".png"))
    file_svg  <- file.path(out_dir, paste0(fig_stub, ".svg"))
    file_xlsx <- file.path(out_dir, paste0(fig_stub, ".xlsx"))
    
    ggplot2::ggsave(file_png, combined, width = width, height = height, dpi = dpi, bg = "white")
    ggplot2::ggsave(file_svg, combined, width = width, height = height, bg = "white")
    
    # underlying wide frame for all metrics
    write.xlsx(x = metrics_wide, file = file_xlsx)
  }
  
  combined
}
