library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# ------------------------------------------------------------
# Apply fixed antibiotic order
# ------------------------------------------------------------
qc_apply_antibiotic_order <- function(qc_tbl, ab_order) {
  qc_tbl %>%
    mutate(
      antibiotic = factor(antibiotic, levels = ab_order)
    )
}

# ------------------------------------------------------------
# Convert low/high rows to wide
# ------------------------------------------------------------
qc_limits_wide <- function(limits_tbl, prefix) {
  limits_tbl %>%
    pivot_wider(
      names_from = kind,
      values_from = limit,
      names_prefix = paste0(prefix, "_")
    )
}

# ------------------------------------------------------------
# Signed offset to interval [low, high]
# ------------------------------------------------------------
qc_signed_offset <- function(x, low, high) {
  case_when(
    is.na(x) | is.na(low) | is.na(high) ~ NA_real_,
    x < low  ~ x - low,
    x > high ~ x - high,
    TRUE     ~ 0
  )
}

# ------------------------------------------------------------
# Main QC comparison
# ------------------------------------------------------------
qc_compare <- function(meas_tbl, target_tbl, range_tbl, ab_order = NULL) {
  
  target_w <- qc_limits_wide(target_tbl, "target")
  range_w  <- qc_limits_wide(range_tbl,  "range")
  
  qc <- meas_tbl %>%
    left_join(target_w, by = c("Studienummer", "antibiotic")) %>%
    left_join(range_w,  by = c("Studienummer", "antibiotic")) %>%
    mutate(
      on_target = if_else(
        !is.na(value) & !is.na(target_low) & !is.na(target_high),
        value >= target_low & value <= target_high,
        NA
      ),
      in_range = if_else(
        !is.na(value) & !is.na(range_low) & !is.na(range_high),
        value >= range_low & value <= range_high,
        NA
      ),
      target_offset = qc_signed_offset(value, target_low, target_high),
      range_offset  = qc_signed_offset(value, range_low,  range_high),
      qc_status = case_when(
        is.na(target_low) | is.na(target_high) |
          is.na(range_low)  | is.na(range_high)  ~ "missing limits",
        on_target ~ "on target",
        !on_target & in_range ~ "off target, in range",
        TRUE ~ "out of range"
      )
    )
  
  if (!is.null(ab_order)) {
    qc <- qc_apply_antibiotic_order(qc, ab_order)
  }
  
  qc
}

# ------------------------------------------------------------
# Detailed QC table
# ------------------------------------------------------------
qc_table_detail <- function(qc_tbl) {
  qc_tbl %>%
    arrange(Studienummer, antibiotic) %>%
    select(
      Studienummer, antibiotic, value,
      target_low, target_high, on_target, target_offset,
      range_low, range_high, in_range, range_offset,
      qc_status
    )
}

# ------------------------------------------------------------
# Summary QC table
# ------------------------------------------------------------
qc_table_summary <- function(qc_tbl) {
  qc_tbl %>%
    group_by(Studienummer) %>%
    summarise(
      n = n(),
      on_target_n = sum(on_target %in% TRUE, na.rm = TRUE),
      off_target_in_range_n = sum(qc_status == "off target, in range", na.rm = TRUE),
      out_of_range_n = sum(qc_status == "out of range", na.rm = TRUE),
      missing_limits_n = sum(qc_status == "missing limits", na.rm = TRUE),
      on_target_pct = 100 * on_target_n / n,
      in_range_pct  = 100 * sum(in_range %in% TRUE, na.rm = TRUE) / n,
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# Helper: save plots to manuscript directory
# ------------------------------------------------------------
qc_save_plot <- function(plot, filename, width = 7, height = 9) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )
  invisible(plot)
}

# ------------------------------------------------------------
# Plot: values vs target + range
# ------------------------------------------------------------
qc_plot_bands <- function(qc_tbl, filename) {
  
  p <- ggplot(
    qc_tbl %>% mutate(antibiotic = fct_rev(antibiotic)),
    aes(x = antibiotic, y = value)
  ) +
    geom_linerange(
      aes(ymin = range_low, ymax = range_high),
      linewidth = 0.6, alpha = 0.25, na.rm = TRUE
    ) +
    geom_linerange(
      aes(ymin = target_low, ymax = target_high),
      linewidth = 2.0, alpha = 0.9, na.rm = TRUE
    ) +
    geom_point(
      aes(y = (target_low + target_high) / 2),
      shape = 124, size = 7, stroke = 1.4, na.rm = TRUE
    ) +
    geom_point(aes(shape = qc_status), size = 2.3, na.rm = TRUE) +
    facet_wrap(~Studienummer) +
    coord_flip() +
    scale_y_continuous(limits = c(5, 40)) +   # <-- fixed scale
    labs(
      x = NULL,
      y = "Measured value",
      shape = "QC status"
    ) +
    theme_bw()
  
  qc_save_plot(p, filename)
}


