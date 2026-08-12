library(readr)
library(dplyr)
library(grid)
library(png)
library(jpeg)
library(magick)

make_figure_pdf <- function(
    caption_csv,
    output_pdf,
    page_width = 8.27,
    page_height = 11.69
) {
  captions <- read_delim(
    caption_csv,
    delim = ";",
    col_types = cols(.default = "c")
  )
  
  pdf(output_pdf, width = page_width, height = page_height, onefile = TRUE)
  
  on.exit(dev.off(), add = TRUE)
  
  for (i in seq_len(nrow(captions))) {
    fig_label <- captions$figure[i]
    fig_file <- captions$file[i]
    fig_caption <- captions$caption[i]
    
    if (!file.exists(fig_file)) {
      warning("Figure file not found: ", fig_file)
      next
    }
    
    grid.newpage()
    
    # Page title / caption text
    grid.text(
      paste0(captions$label[i], ". ", captions$caption[i]),
      x = unit(0.7, "in"),
      y = unit(page_height - 0.6, "in"),
      just = c("left", "top"),
      gp = gpar(fontsize = 10)
    )
    
    # Read image
    img <- magick::image_read(fig_file)
    img_raster <- as.raster(img)
    
    # Draw image below caption
    grid.raster(
      img_raster,
      x = unit(0.5, "npc"),
      y = unit(0.45, "npc"),
      width = unit(0.9, "npc"),
      height = unit(0.78, "npc"),
      just = "center",
      interpolate = TRUE
    )
  }
  
  invisible(output_pdf)
}

