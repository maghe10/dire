# flatten_calibration.R
library(dplyr)
library(purrr)

getwd()
setwd(paste(getwd(),"../Confidence-based-Prediction-of-Antibiotic-Resistance",sep="/"))
dir.exists(getwd())
calibration_list <- readRDS("CICP/calibration.rds")

# calibration_list[[model]][[pathogen]][[antibiotic]] has
#   $susceptible and $resistant data.frames with columns: score, n, id

flat_df <- imap_dfr(calibration_list, function(model_list, model_name) {
  imap_dfr(model_list, function(path_list, pathogen_name) {
    imap_dfr(path_list, function(ab_list, antibiotic_name) {
      
      sus <- ab_list$susceptible %>%
        mutate(
          model      = model_name,
          pathogen   = pathogen_name,
          antibiotic = antibiotic_name,
          pheno      = "susceptible"
        )
      
      res <- ab_list$resistant %>%
        mutate(
          model      = model_name,
          pathogen   = pathogen_name,
          antibiotic = antibiotic_name,
          pheno      = "resistant"
        )
      
      bind_rows(sus, res)
    })
  })
})

# Reorder columns a bit
flat_df <- flat_df %>%
  select(model, pathogen, antibiotic, pheno, score, n, id)

write.csv(flat_df, "CICP/calibration_flat.csv", row.names = FALSE)
