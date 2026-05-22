source(file = 'manuscript/manuscriptcommon.R')
library(dplyr)
library(tidyr)


correlationZonemillimeters <- function()
{
  zones <- readModelMillimeterTable()
  sampleCountFrame <- readCountSampleFrameLong()
  head(zones)
  head(sampleCountFrame)
  unique(sampleCountFrame$metric)
  
  frame <- sampleCountFrame %>% filter(cpmode == "normal" & mode == MODE_A & is.na(significanceLevel))
  

  zones_long <- zones %>%
    pivot_longer(
      cols = -sample,
      names_to = "antibiotic",
      values_to = "zone_mm"
    )
  
  df_zone <- frame %>%
    left_join(zones_long, by = c("sample", "antibiotic"))
  

  df_zone <- df_zone %>%
    mutate(
      ambiguity_type = case_when(
        metric %in% c("twolabelS", "twolabelR") ~ "two-label",
        metric %in% c("zerolabelS", "zerolabelR") ~ "zero-label",
        TRUE ~ NA_character_
      )
    )

    
}