# ======================================================================
# generic_plot_helpers.R
# Self-contained ggplot + patchwork helpers for manuscript figures
# ======================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(rlang)
  library(scales)
})

# -------------------------------
# Column resolution (robust NSE)
# -------------------------------

resolve_col <- function(df, col) {
  if (is.character(col) && length(col) == 1) return(col)
  rlang::as_name(rlang::ensym(col))
}

# -------------------------------
# Utilities: output + export
# -------------------------------

resolve_out_dir <- function(out_dir = NULL, default = "manuscriptPlotDirectory") {
  if (!is.null(out_dir)) return(out_dir)
  od <- get0("manuscriptPlotDirectory", inherits = TRUE)
  if (!is.null(od)) return(od)
  default
}

export_plot_bundle <- function(plot,
                               df_raw,
                               out_dir,
                               file_stub,
                               width = 6.5,
                               height = 4,
                               dpi = 300,
                               export = TRUE) {
  if (!isTRUE(export)) return(invisible(NULL))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".png")),
                  plot, width = width, height = height, dpi = dpi, bg = "white")
  invisible(NULL)
}

# -------------------------------
# Utilities: data checks + coercion
# -------------------------------

coerce_numeric_01 <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

assert_shares_01 <- function(x, tol = 1e-8, allow_na = TRUE) {
  if (!is.numeric(x)) stop("Shares must be numeric after coercion.", call. = FALSE)
  if (!allow_na && anyNA(x)) stop("Missing values found (allow_na = FALSE).", call. = FALSE)
  if (any(x < -tol, na.rm = TRUE) || any(x > 1 + tol, na.rm = TRUE)) {
    stop("Shares must be in [0,1] (within tolerance). Found values outside range.", call. = FALSE)
  }
  invisible(TRUE)
}

clamp_01 <- function(x) pmax(0, pmin(1, x))

as_long_01 <- function(df,
                       id_col,
                       measure_cols = NULL,
                       names_to = "group",
                       values_to = "share",
                       check_01 = TRUE,
                       allow_na = TRUE) {
  
  id_name <- resolve_col(df, id_col)
  nm <- names(df)
  
  is_long <- all(c(id_name, names_to, values_to) %in% nm)
  
  if (is_long) {
    long <- df
    long[[values_to]] <- coerce_numeric_01(long[[values_to]])
  } else {
    if (!(id_name %in% nm)) {
      stop(
        sprintf("Column '%s' not found. Available columns: %s",
                id_name, paste(nm, collapse = ", ")),
        call. = FALSE
      )
    }
    
    if (is.null(measure_cols)) {
      measure_cols <- setdiff(nm, id_name)
    } else if (!is.character(measure_cols)) {
      stop("measure_cols must be a character vector of column names.", call. = FALSE)
    }
    
    if (length(measure_cols) == 0) stop("No measure columns found.", call. = FALSE)
    
    df[measure_cols] <- lapply(df[measure_cols], coerce_numeric_01)
    
    long <- tidyr::pivot_longer(
      df,
      cols = tidyselect::all_of(measure_cols),
      names_to = names_to,
      values_to = values_to
    )
  }
  
  if (check_01) assert_shares_01(long[[values_to]], allow_na = allow_na)
  
  long[[id_name]]  <- factor(long[[id_name]], levels = unique(long[[id_name]]))
  long[[names_to]] <- factor(long[[names_to]], levels = unique(long[[names_to]]))
  
  long
}

# -------------------------------
# Utilities: colors + theme + tags
# -------------------------------

greens_ramp <- function(n) {
  grDevices::colorRampPalette(c("#C3D8BB", "#70AD47", "#538233"))(n)
}

# UPDATED: remove ALL support lines (grid lines)
theme_manuscript <- function(base_size = 12,
                             tag_pos = c(0.01, 0.99),
                             legend_pos = "bottom") {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = legend_pos,
      legend.title = ggplot2::element_text(size = base_size - 2),
      legend.text  = ggplot2::element_text(size = base_size - 2),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.tag.position = tag_pos
    )
}

combine_tagged <- function(...,
                           ncol = 2,
                           tag_levels = "A",
                           base_size = 12) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for combine_tagged(). Please install it.", call. = FALSE)
  }
  
  patchwork::wrap_plots(..., ncol = ncol) +
    patchwork::plot_annotation(tag_levels = tag_levels) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = base_size + 2,
        hjust = 0,
        vjust = 1
      ),
      plot.tag.position = "topleft",
      plot.tag.location = "margin"
    )
}


# # ======================================================================
# # Plot: grouped bars (dodged)
# # ======================================================================
# 
# plot_grouped_bars <- function(df,
#                               id_col,
#                               measure_cols = NULL,
#                               names_to = "group",
#                               values_to = "share",
#                               id_levels = NULL,
#                               group_levels = NULL,
#                               y_label = "Unambiguous correct predictions",
#                               legend_title = "Number of input antibiotics",
#                               y_limits = c(0, 1),
#                               percent_accuracy = 1,
#                               dodge_width = 0.75,
#                               bar_width = 0.68,
#                               clamp = FALSE,
#                               tag = NULL,
#                               base_size = 12,
#                               out_dir = NULL,
#                               file_stub = "figure_A",
#                               export = FALSE,
#                               width = 6.5,
#                               height = 4,
#                               check_01 = TRUE,
#                               allow_na = TRUE) {
#   
#   out_dir <- resolve_out_dir(out_dir)
#   id_name <- resolve_col(df, id_col)
#   
#   long <- as_long_01(df, id_col = id_name, measure_cols = measure_cols,
#                      names_to = names_to, values_to = values_to,
#                      check_01 = check_01, allow_na = allow_na)
#   
#   if (isTRUE(clamp)) long[[values_to]] <- clamp_01(long[[values_to]])
#   
#   if (!is.null(id_levels))    long[[id_name]]  <- factor(long[[id_name]], levels = id_levels)
#   if (!is.null(group_levels)) long[[names_to]] <- factor(long[[names_to]], levels = group_levels)
#   
#   group_levs <- levels(long[[names_to]])
#   cols <- setNames(greens_ramp(length(group_levs)), group_levs)
#   
#   p <- ggplot2::ggplot(long, ggplot2::aes(
#     x = !!rlang::sym(id_name),
#     y = !!rlang::sym(values_to),
#     fill = !!rlang::sym(names_to)
#   )) +
#     ggplot2::geom_col(
#       position = ggplot2::position_dodge(width = dodge_width),
#       width = bar_width
#     ) +
#     ggplot2::scale_y_continuous(
#       labels = scales::percent_format(accuracy = percent_accuracy),
#       expand = ggplot2::expansion(mult = c(0, 0.02))
#     ) +
#     ggplot2::coord_cartesian(ylim = y_limits) +
#     ggplot2::scale_fill_manual(values = cols, name = legend_title) +
#     ggplot2::labs(x = NULL, y = y_label, tag = tag) +
#     theme_manuscript(base_size = base_size)
#   
#   export_plot_bundle(p, df, out_dir, file_stub, width, height, export = export)
#   p
# }


plot_grouped_bars <- function(df,
                              id_col,
                              measure_cols = NULL,
                              names_to = "group",
                              values_to = "share",
                              id_levels = NULL,
                              group_levels = NULL,
                              y_label = "Unambiguous correct predictions",
                              legend_title = "Number of input antibiotics",
                              y_limits = c(0, 1),
                              y_as_percent = TRUE,          # NEW
                              y_accuracy = 1,               # NEW (replaces percent_accuracy)
                              dodge_width = 0.75,
                              bar_width = 0.68,
                              clamp = FALSE,
                              tag = NULL,
                              base_size = 12,
                              out_dir = NULL,
                              file_stub = "figure_A",
                              export = FALSE,
                              width = 6.5,
                              height = 4,
                              check_01 = TRUE,
                              allow_na = TRUE) {
  
  out_dir <- resolve_out_dir(out_dir)
  id_name <- resolve_col(df, id_col)
  
  long <- as_long_01(
    df,
    id_col = id_name,
    measure_cols = measure_cols,
    names_to = names_to,
    values_to = values_to,
    check_01 = check_01,
    allow_na = allow_na
  )
  
  if (isTRUE(clamp)) long[[values_to]] <- clamp_01(long[[values_to]])
  
  if (!is.null(id_levels))    long[[id_name]]  <- factor(long[[id_name]], levels = id_levels)
  if (!is.null(group_levels)) long[[names_to]] <- factor(long[[names_to]], levels = group_levels)
  
  group_levs <- levels(long[[names_to]])
  cols <- setNames(greens_ramp(length(group_levs)), group_levs)
  
  y_labeller <- if (isTRUE(y_as_percent)) {
    scales::percent_format(accuracy = y_accuracy)
  } else {
    scales::number_format(accuracy = y_accuracy)
  }
  
  p <- ggplot2::ggplot(long, ggplot2::aes(
    x = !!rlang::sym(id_name),
    y = !!rlang::sym(values_to),
    fill = !!rlang::sym(names_to)
  )) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = dodge_width),
      width = bar_width
    ) +
    ggplot2::scale_y_continuous(
      labels = y_labeller,
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::scale_fill_manual(values = cols, name = legend_title) +
    ggplot2::labs(x = NULL, y = y_label, tag = tag) +
    theme_manuscript(base_size = base_size)
  
  export_plot_bundle(p, df, out_dir, file_stub, width, height, export = export)
  p
}


# ======================================================================
# Plot: stacked bars
# ======================================================================

plot_stacked_bars <- function(df,
                              x_col,
                              fill_col,
                              measure_cols = NULL,
                              values_to = "share",
                              x_levels = NULL,
                              fill_levels = NULL,      # bottom -> top
                              fill_colors = NULL,      # named vector
                              fill_labels = fill_levels,
                              y_label = "Percentages of prediction results",
                              y_limits = c(0, 1),
                              percent_accuracy = 1,
                              bar_width = 0.72,
                              clamp = TRUE,
                              stack_reverse = TRUE,    # <-- key fix: makes fill_levels behave bottom->top
                              tag = NULL,
                              base_size = 12,
                              out_dir = NULL,
                              file_stub = "figure_B",
                              export = FALSE,
                              width = 6.5,
                              height = 4,
                              check_01 = TRUE,
                              allow_na = TRUE) {
  
  out_dir <- resolve_out_dir(out_dir)
  
  x_name <- resolve_col(df, x_col)
  fill_name <- resolve_col(df, fill_col)
  nm <- names(df)
  
  is_long <- all(c(x_name, fill_name, values_to) %in% nm)
  
  if (!is_long) {
    if (!(fill_name %in% nm)) {
      stop(
        sprintf("Column '%s' not found. Available columns: %s",
                fill_name, paste(nm, collapse = ", ")),
        call. = FALSE
      )
    }
    if (is.null(measure_cols)) {
      measure_cols <- setdiff(nm, fill_name)
    } else if (!is.character(measure_cols)) {
      stop("measure_cols must be a character vector of column names.", call. = FALSE)
    }
    if (length(measure_cols) == 0) stop("No measure columns found for wide stacked plot.", call. = FALSE)
    
    df2 <- df
    df2[measure_cols] <- lapply(df2[measure_cols], coerce_numeric_01)
    
    long <- tidyr::pivot_longer(df2,
                                cols = tidyselect::all_of(measure_cols),
                                names_to = x_name,
                                values_to = values_to)
  } else {
    long <- df
    long[[values_to]] <- coerce_numeric_01(long[[values_to]])
  }
  
  if (check_01) assert_shares_01(long[[values_to]], allow_na = allow_na)
  if (isTRUE(clamp)) long[[values_to]] <- clamp_01(long[[values_to]])
  
  if (!is.null(x_levels)) {
    long[[x_name]] <- factor(long[[x_name]], levels = x_levels)
  } else {
    long[[x_name]] <- factor(long[[x_name]], levels = unique(long[[x_name]]))
  }
  
  if (!is.null(fill_levels)) {
    long[[fill_name]] <- factor(long[[fill_name]], levels = fill_levels)
  } else {
    long[[fill_name]] <- factor(long[[fill_name]], levels = unique(long[[fill_name]]))
  }
  
  # Colors + legend order consistent with fill_levels
  if (is.null(fill_colors)) {
    levs <- levels(long[[fill_name]])
    fill_colors <- setNames(greens_ramp(length(levs)), levs)
  } else if (is.null(names(fill_colors))) {
    stop("fill_colors must be a named vector.", call. = FALSE)
  }
  
  breaks_use <- if (!is.null(fill_levels)) fill_levels else levels(long[[fill_name]])
  labels_use <- if (!is.null(fill_labels)) fill_labels else breaks_use
  
  p <- ggplot2::ggplot(long, ggplot2::aes(
    x = !!rlang::sym(x_name),
    y = !!rlang::sym(values_to),
    fill = !!rlang::sym(fill_name)
  )) +
    ggplot2::geom_col(width = bar_width,
                      position = ggplot2::position_stack(reverse = stack_reverse)) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = percent_accuracy),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::scale_fill_manual(values = fill_colors, breaks = breaks_use, labels = labels_use, name = NULL) +
    ggplot2::labs(x = NULL, y = y_label, tag = tag) +
    theme_manuscript(base_size = base_size)
  
  export_plot_bundle(p, df, out_dir, file_stub, width, height, export = export)
  p
}

plot_grouped_stacked_bars <- function(df,
                                      x_col,
                                      group_col,
                                      stack_col,
                                      values_to = "share",
                                      x_levels = NULL,
                                      group_levels = NULL,
                                      stack_levels = NULL,       # bottom -> top
                                      group_colors = NULL,       # optional: colors for groups (outline/legend not used by default)
                                      stack_colors = NULL,       # named vector: metric colors
                                      y_label = "Percentages",
                                      legend_title = NULL,
                                      y_limits = c(0, 1),
                                      percent_accuracy = 1,
                                      dodge_width = 0.78,
                                      bar_width = 0.66,
                                      clamp = TRUE,
                                      stack_reverse = TRUE,      # keeps stack_levels bottom->top visually
                                      tag = NULL,
                                      base_size = 12,
                                      check_01 = TRUE,
                                      allow_na = TRUE) {
  # LONG ONLY contract
  x_name     <- resolve_col(df, x_col)
  group_name <- resolve_col(df, group_col)
  stack_name <- resolve_col(df, stack_col)
  
  required <- c(x_name, group_name, stack_name, values_to)
  nm <- names(df)
  if (!all(required %in% nm)) {
    stop(
      sprintf(
        "plot_grouped_stacked_bars() expects LONG data with columns: %s. Available columns: %s",
        paste(required, collapse = ", "),
        paste(nm, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  long <- df
  long[[values_to]] <- coerce_numeric_01(long[[values_to]])
  
  if (check_01) assert_shares_01(long[[values_to]], allow_na = allow_na)
  if (isTRUE(clamp)) long[[values_to]] <- clamp_01(long[[values_to]])
  
  if (!is.null(x_levels)) {
    long[[x_name]] <- factor(long[[x_name]], levels = x_levels)
  } else {
    long[[x_name]] <- factor(long[[x_name]], levels = unique(long[[x_name]]))
  }
  
  if (!is.null(group_levels)) {
    long[[group_name]] <- factor(long[[group_name]], levels = group_levels)
  } else {
    long[[group_name]] <- factor(long[[group_name]], levels = unique(long[[group_name]]))
  }
  
  if (!is.null(stack_levels)) {
    long[[stack_name]] <- factor(long[[stack_name]], levels = stack_levels)
  } else {
    long[[stack_name]] <- factor(long[[stack_name]], levels = unique(long[[stack_name]]))
  }
  
  # Stack colors (required for clear metric legend)
  if (is.null(stack_colors)) {
    levs <- levels(long[[stack_name]])
    stack_colors <- setNames(greens_ramp(length(levs)), levs)
  } else if (is.null(names(stack_colors))) {
    stop("stack_colors must be a named vector.", call. = FALSE)
  }
  breaks_use <- if (!is.null(stack_levels)) stack_levels else levels(long[[stack_name]])
  
  # Key trick: x = interaction(x, group) and dodge on group
  # This yields grouped stacks with a clean group legend if desired.
  long$.xg <- interaction(long[[x_name]], long[[group_name]], lex.order = TRUE, drop = TRUE)
  
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = .xg,
      y = !!rlang::sym(values_to),
      fill = !!rlang::sym(stack_name),
      group = !!rlang::sym(group_name)
    )
  ) +
    ggplot2::geom_col(
      width = bar_width,
      position = ggplot2::position_stack(reverse = stack_reverse)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = percent_accuracy),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::scale_fill_manual(values = stack_colors, breaks = breaks_use, name = legend_title) +
    ggplot2::labs(x = NULL, y = y_label, tag = tag) +
    theme_manuscript(base_size = base_size) +
    # Show only the main x labels (antibiotic), not the interaction labels
    ggplot2::scale_x_discrete(labels = function(labs) sub("\\..*$", "", labs)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  
  p
}

export_patchwork <- function(patch,
                             out_dir = NULL,
                             file_stub = "figure_AB",
                             width = 13,
                             height = 4,
                             dpi = 300,
                             bg = "white",
                             export = TRUE) {
  if (!isTRUE(export)) return(invisible(NULL))
  
  out_dir <- resolve_out_dir(out_dir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(file_stub, ".png")),
    plot = patch,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg
  )
  
  invisible(NULL)
}


# ======================================================================
# Examples (self-contained)
# ======================================================================

example_grouped_wide_basic <- function() {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for examples. Please install it.", call. = FALSE)
  }
  
  df <- tibble::tibble(
    significanceLevel = c("non-conformal", "10%", "5%", "2.5%"),
    `4` = c(0.806, 0.713, 0.699, 0.644),
    `5` = c(0.826, 0.747, 0.729, 0.677),
    `6` = c(0.842, 0.773, 0.754, 0.703),
    `7` = c(0.854, 0.794, 0.774, 0.723),
    `8` = c(0.864, 0.811, 0.790, 0.738)
  )
  
  plot_grouped_bars(
    df,
    id_col = "significanceLevel",
    id_levels = c("non-conformal", "10%", "5%", "2.5%"),
    group_levels = c("4", "5", "6", "7", "8"),
    legend_title = "Number of input antibiotics",
    export = FALSE
  )
}

example_grouped_wide_basic_not_percent <- function() {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for examples. Please install it.", call. = FALSE)
  }
  
  df <- tibble::tibble(
    significanceLevel = c("non-conformal", "10%", "5%", "2.5%"),
    `4` = c(0.806, 0.713, 0.699, 0.644),
    `5` = c(0.826, 0.747, 0.729, 0.677),
    `6` = c(0.842, 0.773, 0.754, 0.703),
    `7` = c(0.854, 0.794, 0.774, 0.723),
    `8` = c(0.864, 0.811, 0.790, 0.738)
  )
  
  plot_grouped_bars(
    df,
    id_col = "significanceLevel",
    id_levels = c("non-conformal", "10%", "5%", "2.5%"),
    group_levels = c("4", "5", "6", "7", "8"),
    legend_title = "Number of input antibiotics",
    export = FALSE,
    y_as_percent = FALSE,
    y_accuracy = 0.05,   # e.g. show 0.00, 0.05, 0.10 ...
    y_limits = c(0, 1)
  )
}



example_stacked_wide_metrics <- function() {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for examples. Please install it.", call. = FALSE)
  }
  
  df6 <- tibble::tibble(
    metric = c("correct", "ME", "VME"),
    `non-conformal` = c(0.842, 0.106, 0.0513),
    `10%` = c(0.773, 0.141, 0.0454),
    `5%` = c(0.754, 0.0919, 0.0276),
    `2.5%` = c(0.703, 0.0612, 0.0170)
  )
  
  metric_colors <- c(correct = "#70AD47", ME = "#FFFF00", VME = "#FF0000")
  
  # fill_levels is bottom -> top, and stack_reverse=TRUE makes it actually draw that way
  plot_stacked_bars(
    df6,
    x_col = "significanceLevel",
    fill_col = "metric",
    measure_cols = c("non-conformal", "10%", "5%", "2.5%"),
    x_levels = c("non-conformal", "10%", "5%", "2.5%"),
    fill_levels = c("correct", "ME", "VME"),
    fill_colors = metric_colors,
    stack_reverse = TRUE,
    clamp = TRUE,
    export = FALSE
  )
}

example_grouped_stacked_metrics_by_index <- function() {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for examples. Please install it.", call. = FALSE)
  }
  
  antibiotics <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")
  idx_levels <- c("4","5","6","7","8")
  
  # Self-contained demo data (using your provided numbers where available; rest are plausible)
  metrics4 <- tibble::tibble(
    antibiotic = antibiotics,
    correct = c(0.883,0.655,0.906,0.728,0.770,0.837,0.831,0.833,0.842,0.854,0.817,0.852,0.734,0.743),
    ME      = c(0.0692,0.283,0.0341,0.234,0.0981,0.0689,0.0366,0.0604,0.107,0.0373,0.154,0.0254,0.238,0.233),
    VME     = c(0.0480,0.0621,0.0604,0.0384,0.132,0.0940,0.132,0.106,0.0515,0.109,0.0291,0.123,0.0288,0.0236)
  )
  
  metrics6 <- tibble::tibble(
    antibiotic = antibiotics,
    correct = c(0.935,0.668,0.941,0.724,0.804,0.888,0.890,0.885,0.877,0.876,0.854,0.874,0.797,0.778),
    ME      = c(0.0333,0.293,0.0256,0.242,0.102,0.0634,0.0262,0.0568,0.0885,0.0220,0.136,0.0128,0.182,0.206),
    VME     = c(0.0315,0.0381,0.0337,0.0340,0.0939,0.0484,0.0834,0.0581,0.0349,0.102,0.0104,0.114,0.0211,0.0155)
  )
  
  metrics8 <- tibble::tibble(
    antibiotic = antibiotics,
    correct = c(0.965,0.677,0.958,0.728,0.821,0.917,0.920,0.915,0.893,0.881,0.881,0.882,0.854,0.807),
    ME      = c(0.0163,0.297,0.0235,0.240,0.105,0.0604,0.0183,0.0534,0.0789,0.0108,0.116,0.00538,0.129,0.182),
    VME     = c(0.0190,0.0261,0.0180,0.0313,0.0739,0.0229,0.0616,0.0316,0.0285,0.108,0.00319,0.113,0.0167,0.0103)
  )
  
  # For a complete 4..8 demo, synthesize 5 and 7 by interpolation
  metrics5 <- metrics4
  metrics5$correct <- (metrics4$correct + metrics6$correct) / 2
  metrics5$ME      <- (metrics4$ME      + metrics6$ME)      / 2
  metrics5$VME     <- (metrics4$VME     + metrics6$VME)     / 2
  
  metrics7 <- metrics6
  metrics7$correct <- (metrics6$correct + metrics8$correct) / 2
  metrics7$ME      <- (metrics6$ME      + metrics8$ME)      / 2
  metrics7$VME     <- (metrics6$VME     + metrics8$VME)     / 2
  
  # Long format with both grouping and stacking
  long4 <- tidyr::pivot_longer(metrics4, cols = c("correct","ME","VME"),
                               names_to = "metric", values_to = "share") %>%
    mutate(index = "4")
  long5 <- tidyr::pivot_longer(metrics5, cols = c("correct","ME","VME"),
                               names_to = "metric", values_to = "share") %>%
    mutate(index = "5")
  long6 <- tidyr::pivot_longer(metrics6, cols = c("correct","ME","VME"),
                               names_to = "metric", values_to = "share") %>%
    mutate(index = "6")
  long7 <- tidyr::pivot_longer(metrics7, cols = c("correct","ME","VME"),
                               names_to = "metric", values_to = "share") %>%
    mutate(index = "7")
  long8 <- tidyr::pivot_longer(metrics8, cols = c("correct","ME","VME"),
                               names_to = "metric", values_to = "share") %>%
    mutate(index = "8")
  
  long <- dplyr::bind_rows(long4, long5, long6, long7, long8) %>%
    mutate(
      index = factor(index, levels = idx_levels),
      metric = factor(metric, levels = c("correct","ME","VME"))
    )
  
  metric_colors <- c(correct = "#70AD47", ME = "#FFFF00", VME = "#FF0000")
  
  plot_grouped_stacked_bars(
    long,
    x_col = "antibiotic",
    group_col = "index",
    stack_col = "metric",
    values_to = "share",
    x_levels = antibiotics,
    group_levels = idx_levels,
    stack_levels = c("correct","ME","VME"),
    stack_colors = metric_colors,
    legend_title = NULL,
    y_label = "Percentages of prediction results",
    y_limits = c(0, 1),
    clamp = TRUE,
    stack_reverse = TRUE,
    tag = "A"
  )
}



pad_facets_to_full_row <- function(long,
                                   facet_col = "antibiotic",
                                   ncol = 4,
                                   x_col = "index",
                                   fill_col = "metric",
                                   values_col = "share") {
  
  facet_col  <- rlang::ensym(facet_col)
  x_col      <- rlang::ensym(x_col)
  fill_col   <- rlang::ensym(fill_col)
  values_col <- rlang::ensym(values_col)
  
  facet_vals <- long %>% dplyr::distinct(!!facet_col) %>% dplyr::pull(!!facet_col)
  n <- length(facet_vals)
  rem <- n %% ncol
  if (rem == 0) return(long)
  
  n_pad <- ncol - rem
  pad_levels <- paste0("__PAD", seq_len(n_pad), "__")
  
  # IMPORTANT: keep original types (integer stays integer, etc.)
  x_vals <- long %>% dplyr::distinct(!!x_col) %>% dplyr::pull(!!x_col)
  fill_vals <- long %>% dplyr::distinct(!!fill_col) %>% dplyr::pull(!!fill_col)
  
  pad <- tidyr::expand_grid(
    !!facet_col := pad_levels,
    !!x_col := x_vals,
    !!fill_col := fill_vals
  ) %>%
    dplyr::mutate(!!values_col := 0)
  
  long2 <- dplyr::bind_rows(long, pad)
  
  # Force facet column to include pad levels, but keep existing order
  all_levels <- c(as.character(facet_vals), pad_levels)
  long2[[rlang::as_name(facet_col)]] <- factor(long2[[rlang::as_name(facet_col)]],
                                               levels = all_levels)
  
  long2
}

example_faceted_by_antibiotic_metrics_vs_index <- function(export = FALSE,
                                                           out_dir = NULL,
                                                           file_stub = "figure_faceted_by_antibiotic",
                                                           ncol = 4) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for examples. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    # facet_wrap does not require patchwork, but we keep checks consistent with other examples
    message("patchwork not installed (not needed for this example).")
  }
  
  antibiotics <- c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB")
  idx_levels <- c("4","5","6","7","8")
  
  # --- Self-contained demo data (your provided values for 4/6/8; 5 and 7 interpolated)
  metrics4 <- tibble::tibble(
    antibiotic = antibiotics,
    correct = c(0.883,0.655,0.906,0.728,0.770,0.837,0.831,0.833,0.842,0.854,0.817,0.852,0.734,0.743),
    ME      = c(0.0692,0.283,0.0341,0.234,0.0981,0.0689,0.0366,0.0604,0.107,0.0373,0.154,0.0254,0.238,0.233),
    VME     = c(0.0480,0.0621,0.0604,0.0384,0.132,0.0940,0.132,0.106,0.0515,0.109,0.0291,0.123,0.0288,0.0236)
  )
  
  metrics6 <- tibble::tibble(
    antibiotic = antibiotics,
    correct = c(0.935,0.668,0.941,0.724,0.804,0.888,0.890,0.885,0.877,0.876,0.854,0.874,0.797,0.778),
    ME      = c(0.0333,0.293,0.0256,0.242,0.102,0.0634,0.0262,0.0568,0.0885,0.0220,0.136,0.0128,0.182,0.206),
    VME     = c(0.0315,0.0381,0.0337,0.0340,0.0939,0.0484,0.0834,0.0581,0.0349,0.102,0.0104,0.114,0.0211,0.0155)
  )
  
  metrics8 <- tibble::tibble(
    antibiotic = antibiotics,
    correct = c(0.965,0.677,0.958,0.728,0.821,0.917,0.920,0.915,0.893,0.881,0.881,0.882,0.854,0.807),
    ME      = c(0.0163,0.297,0.0235,0.240,0.105,0.0604,0.0183,0.0534,0.0789,0.0108,0.116,0.00538,0.129,0.182),
    VME     = c(0.0190,0.0261,0.0180,0.0313,0.0739,0.0229,0.0616,0.0316,0.0285,0.108,0.00319,0.113,0.0167,0.0103)
  )
  
  # Interpolate 5 between 4 and 6; interpolate 7 between 6 and 8
  metrics5 <- metrics4
  metrics5$correct <- (metrics4$correct + metrics6$correct) / 2
  metrics5$ME      <- (metrics4$ME      + metrics6$ME)      / 2
  metrics5$VME     <- (metrics4$VME     + metrics6$VME)     / 2
  
  metrics7 <- metrics6
  metrics7$correct <- (metrics6$correct + metrics8$correct) / 2
  metrics7$ME      <- (metrics6$ME      + metrics8$ME)      / 2
  metrics7$VME     <- (metrics6$VME     + metrics8$VME)     / 2
  
  make_long <- function(df, idx) {
    tidyr::pivot_longer(
      df,
      cols = c("correct","ME","VME"),
      names_to = "metric",
      values_to = "share"
    ) %>%
      dplyr::mutate(index = idx)
  }
  
  long <- dplyr::bind_rows(
    make_long(metrics4, "4"),
    make_long(metrics5, "5"),
    make_long(metrics6, "6"),
    make_long(metrics7, "7"),
    make_long(metrics8, "8")
  ) %>%
    dplyr::mutate(
      antibiotic = factor(antibiotic, levels = antibiotics),
      index      = factor(index, levels = idx_levels),
      metric     = factor(metric, levels = c("correct","ME","VME"))
    )
  
  long <- pad_facets_to_full_row(long, facet_col = "antibiotic", ncol = ncol,
                                 x_col = "index", fill_col = "metric", values_col = "share")
  
  pad_labeller <- ggplot2::labeller(
    antibiotic = function(x) ifelse(grepl("^__PAD\\d+__$", x), "", x)
  )
  

  
  metric_colors <- c(correct = "#70AD47", ME = "#FFFF00", VME = "#FF0000")
  
  p <- plot_stacked_bars(
    long,
    x_col = "index",
    fill_col = "metric",
    values_to = "share",
    x_levels = idx_levels,
    fill_levels = c("correct","ME","VME"),
    fill_colors = metric_colors,
    y_label = "Percentages of prediction results",
    y_limits = c(0, 1),
    clamp = TRUE,
    stack_reverse = TRUE,
    tag = "A",
    export = FALSE,
    bar_width = 0.95
  ) +
    ggplot2::facet_wrap(~ antibiotic, ncol = ncol, drop = FALSE, labeller = pad_labeller) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      strip.text  = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(0.35, "lines")
    )
  
  print(p)
  
  if (isTRUE(export)) {
    out_dir <- resolve_out_dir(out_dir)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    # Size suggestion for 14 facets; adjust if you change ncol
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".png")),
                    p, width = 12, height = 8, dpi = 300, bg = "white")
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".svg")),
                    p, width = 12, height = 8, bg = "white")
  }
  
  invisible(p)
}

example_patchwork_AB <- function(export = FALSE,
                                 out_dir = NULL,
                                 file_stub = "figure_AB") {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for patchwork examples. Please install it.", call. = FALSE)
  }
  
  pA <- example_grouped_wide_basic()
  pB <- example_stacked_wide_metrics()
  
  pAB <- combine_tagged(pA, pB, ncol = 2, tag_levels = "A", base_size = 12)
  print(pAB)
  
  export_patchwork(
    pAB,
    out_dir = out_dir,
    file_stub = file_stub,
    width = 13,   # two 6.5-inch panels side-by-side
    height = 4,
    export = export
  )
  
  invisible(pAB)
}


example_patchwork_fig7AB <- function(export = FALSE,
                                     out_dir = NULL,
                                     file_stub = "figure_7AB",
                                     pick = "6") {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for examples. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for this example. Please install it.", call. = FALSE)
  }
  
  # --- Data (self-contained; mirrors your Fig7A/Fig7B structure)
  frameCorrectAntibioticsVsIndex <- tibble::tibble(
    antibiotic = c("AMP","AMC","PIP","TZP","CAZ","CRO","CTX","FEP","CIP","OFX","LVX","MFX","GEN","TOB"),
    `4` = c(0.883,0.655,0.906,0.728,0.770,0.837,0.831,0.833,0.842,0.854,0.817,0.852,0.734,0.743),
    `5` = c(0.913,0.662,0.927,0.726,0.789,0.866,0.865,0.862,0.862,0.868,0.837,0.865,0.766,0.759),
    `6` = c(0.935,0.668,0.941,0.724,0.804,0.888,0.890,0.885,0.877,0.876,0.854,0.874,0.797,0.778),
    `7` = c(0.952,0.674,0.951,0.725,0.814,0.904,0.907,0.902,0.886,0.880,0.868,0.879,0.826,0.793),
    `8` = c(0.965,0.677,0.958,0.728,0.821,0.917,0.920,0.915,0.893,0.881,0.881,0.882,0.854,0.807)
  )
  
  frameAntibioticsVsMetric <- list(
    `4` = tibble::tibble(
      antibiotic = frameCorrectAntibioticsVsIndex$antibiotic,
      correct = frameCorrectAntibioticsVsIndex$`4`,
      ME  = c(0.0692,0.283,0.0341,0.234,0.0981,0.0689,0.0366,0.0604,0.107,0.0373,0.154,0.0254,0.238,0.233),
      VME = c(0.0480,0.0621,0.0604,0.0384,0.132,0.0940,0.132,0.106,0.0515,0.109,0.0291,0.123,0.0288,0.0236)
    ),
    `6` = tibble::tibble(
      antibiotic = frameCorrectAntibioticsVsIndex$antibiotic,
      correct = frameCorrectAntibioticsVsIndex$`6`,
      ME  = c(0.0333,0.293,0.0256,0.242,0.102,0.0634,0.0262,0.0568,0.0885,0.0220,0.136,0.0128,0.182,0.206),
      VME = c(0.0315,0.0381,0.0337,0.0340,0.0939,0.0484,0.0834,0.0581,0.0349,0.102,0.0104,0.114,0.0211,0.0155)
    ),
    `8` = tibble::tibble(
      antibiotic = frameCorrectAntibioticsVsIndex$antibiotic,
      correct = frameCorrectAntibioticsVsIndex$`8`,
      ME  = c(0.0163,0.297,0.0235,0.240,0.105,0.0604,0.0183,0.0534,0.0789,0.0108,0.116,0.00538,0.129,0.182),
      VME = c(0.0190,0.0261,0.0180,0.0313,0.0739,0.0229,0.0616,0.0316,0.0285,0.108,0.00319,0.113,0.0167,0.0103)
    )
  )
  
  # --- Panel A: grouped bars (antibiotic x #inputs)
  pA <- plot_grouped_bars(
    frameCorrectAntibioticsVsIndex,
    id_col = "antibiotic",
    measure_cols = c("4","5","6","7","8"),
    group_levels = c("4","5","6","7","8"),
    dodge_width = 0.78,
    bar_width = 0.66,
    legend_title = "Number of input antibiotics",
    y_label = "Unambiguous correct predictions",
    y_limits = c(0, 1),
    tag = "A",
    export = FALSE
  ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  
  # --- Panel B: STACKED needs long form for our generic function (self-contained)
  metric_colors <- c(correct = "#70AD47", ME = "#FFFF00", VME = "#FF0000")
  df_pick <- frameAntibioticsVsMetric[[as.character(pick)]]
  
  df_pick_long <- tidyr::pivot_longer(
    df_pick,
    cols = c("correct", "ME", "VME"),
    names_to = "metric",
    values_to = "share"
  )
  
  pB <- plot_stacked_bars(
    df_pick_long,
    x_col = "antibiotic",
    fill_col = "metric",
    values_to = "share",
    x_levels = df_pick$antibiotic,              # keep given order
    fill_levels = c("correct", "ME", "VME"),    # bottom -> top
    fill_colors = metric_colors,
    y_label = "Percentages of prediction results",
    y_limits = c(0, 1),
    tag = "B",
    clamp = TRUE,
    stack_reverse = TRUE,
    export = FALSE
  ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  
  # --- Combine
  pAB <- combine_tagged(pA, pB, ncol = 2, tag_levels = "A", base_size = 12)
  print(pAB)
  
  # --- Optional export of combined figure (patchwork)
  if (isTRUE(export)) {
    out_dir <- resolve_out_dir(out_dir)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".png")),
                    pAB, width = 17.3, height = 5.2, dpi = 300, bg = "white")
    ggplot2::ggsave(file.path(out_dir, paste0(file_stub, ".svg")),
                    pAB, width = 17.3, height = 5.2, bg = "white")
  }
  
  invisible(pAB)
}




EXAMPLES <- function() {
  p1 <- example_grouped_wide_basic();  print(p1)
  p2 <- example_stacked_wide_metrics();print(p2)
  p3 <- example_grouped_wide_basic_not_percent();print(p3)
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    example_patchwork_AB()
    example_patchwork_fig7AB()
  } else {
    message("patchwork not installed; skipping patchwork examples.")
  }
  
  invisible(list(grouped = p1, stacked = p2))
}
