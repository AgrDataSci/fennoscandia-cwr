# Make a boxplot of AUCs from model outputs
# this produces Fig. S7

library("ggplot2")
library("svglite")


spp_dt <- read.csv("data/species_names.csv")

spp <- list.files("processing/enm")

# read AUC files
auc <- data.frame()

for(i in seq_along(spp)){
  sp_i <- spp[[i]]
  f <- paste0("processing/enm/", sp_i, "/outputs/auc_testing.csv")
  f <- read.csv(f)
  f <- data.frame(acronym = sp_i,
                  AUC = as.vector(na.omit(f$auc)))
  
  auc <- rbind(auc, f)
  
}

auc <- merge(auc, spp_dt[,c("acronym","genus","species","authority")], by = "acronym")

auc$taxa <- with(auc, paste(genus, species, authority))

auc$taxa <- factor(auc$taxa, levels = sort(unique(auc$taxa)))

p <- ggplot(auc, aes(x = taxa, y = AUC)) +
  geom_boxplot() + 
  ylim(0.6, 1) +
  labs(x = NULL,
       y = "AUC") + 
  theme(axis.text.y = element_text(size = 9, angle = 0, 
                                   hjust = 1, vjust = 0.5, colour = "black"),
        axis.text.x = element_text(size = 9, angle = 70, 
                                   hjust = 1, vjust = 1, colour = "black"),
        axis.title = element_text(size = 9, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 9, colour = "black"),
        strip.text.y = element_text(size=9, colour = "black"),
        strip.background = element_rect(colour="black", fill="#FFFFFF"))

p

output <- "output/auc_models"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

ggsave(paste0(output, "/auc_models.png"),
       plot = p,
       width = 10,
       height = 9,
       dpi = 500,
       units = "cm")



