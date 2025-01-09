#' -----------------------------------------
#' Standardized helper functions and styles
#' -----------------------------------------

library(ggplot2)
theme_set(theme_classic(base_size = 14)+
            theme(axis.text.x = 
                    element_text(angle = 40, 
                                 vjust = 1, 
                                 hjust=1),
                  legend.position = "bottom"))

berm_control_colors <- c("#af8dc3", "#7fbf7b")

fix_data_for_plots <- function(dat){
  dat |>
    # make duxbury not have a _
    mutate(site = gsub("_", " ", site)) |>
    
    # make zone a factor and order so it's low, mid, high
    mutate( height = factor(height, levels = c("Low", "Mid", "High")))
}
