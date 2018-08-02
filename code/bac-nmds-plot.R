bacord <- readRDS("data/bac-ordination.Rds")
bacord
bacord_plots <- bacord$points %>% data.frame
bacord_family <- bacord$species %>% data.frame
bacord_family$names <- rownames(bacord$species)
library(tidyverse)

theme_gsk <- function() {
  theme_minimal()+
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.tag = element_text(face = "bold")
    ) 
}
ggplot(bacord_plots) + 
  geom_point(aes(x = MDS1, y = MDS2), 
             color = "tomato4", size = 4) +
  geom_point(data = bacord_family, 
             aes(x = MDS1, y = MDS2),
             size = .5, color = alpha("grey25", 0.75)) + 
  geom_text_repel(data = bacord_family, 
             aes(x = MDS1, y = MDS2, label = names),
             size = 2.5, color = alpha("grey25", 0.75),
             segment.alpha = 0, fontface = "italic") + 
  theme_gsk() +
  theme(axis.title = element_text(size = 16),
        axis.line = element_line(colour = 'grey25', size = 0.8))
  
