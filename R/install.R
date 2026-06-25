# Main packages you explicitly work with
main_pkgs <- c(
  "tidyverse",
  "cluster",
  "factoextra",
  "readxl",
  "writexl",
  "openxlsx",
  "pheatmap",
  "png",
  "reticulate",
  "gridpattern",
  "R.utils",
  "svglite",
  "patchwork",
  "fpc",
  "fst",
  "arrow",
  "ggtext",
  "furrr",
  "future",
  "janitor",
  "httr2",
  "uwot",
  "pvclust",
  "vegan",
  "ape",
  "magick"
)

# Dependencies (for reference only)
dep_pkgs <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "scales",
  "forcats",
  "ggplot2"
)

# Install only missing main packages
to_install <- setdiff(main_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) {
  install.packages(to_install)
}

# ---- Dependency Check ----
missing_deps <- setdiff(dep_pkgs, rownames(installed.packages()))

if (length(missing_deps) == 0) {
  message("All dependency packages are installed.")
} else {
  message("Missing dependency packages:\n  - ",
          paste(missing_deps, collapse = "\n  - "))
}

