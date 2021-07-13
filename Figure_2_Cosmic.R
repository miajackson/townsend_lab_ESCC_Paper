library(cancereffectsizeR)
library(ggplot2)
library(cowplot)

setwd("/Users/miajackson/Documents/Townsend Lab/ESCC_Paper")
#load("Analysis.12.17.ep.results.Rdata")
load("Analysis.5.7.results.Rdata")

#####################################################################################################################
twostage_final <- ep_results
twostage_results <- snv_results(twostage_final)
twostage_results <- twostage_results$selection.1

#seperating into pre and pri
tumor_names_twostage <- unique(twostage_final@maf$Unique_Patient_Identifier)

tumor_pre <- c()
tumor_pri <- c()

samples_twostage <- twostage_final@samples
signature_table_twostage <- twostage_final@trinucleotide_mutation_weights$signature_weight_table

for (x in tumor_names_twostage){
  if (samples_twostage$group[which(samples_twostage$Unique_Patient_Identifier == x)] == "Pre"){
    tumor_pre <- c(tumor_pre, x)
  }
  else if (samples_twostage$group[which(samples_twostage$Unique_Patient_Identifier == x)] == "Pri"){
    tumor_pri <- c(tumor_pri, x)
  }
}

tumor_pre_signature <- data.frame(matrix(ncol = 72, nrow = 0))
colnames(tumor_pre_signature) <- colnames(signature_table_twostage)[5:76]

tumor_pri_signature <- data.frame(matrix(ncol = 72, nrow = 0))
colnames(tumor_pri_signature) <- colnames(signature_table_twostage)[5:76]

#Separating primary and tumors and extracting their cosmic signatures
for(x in tumor_pre){
  tumor_pre_signature <- rbind(tumor_pre_signature, signature_table_twostage[which(signature_table_twostage$Unique_Patient_Identifier == x),5:76])
}

for(x in tumor_pri){
  tumor_pri_signature <- rbind(tumor_pri_signature, signature_table_twostage[which(signature_table_twostage$Unique_Patient_Identifier == x),5:76])
}

#Remove columns that do not have any weight
remove_col <- which(colnames(tumor_pre_signature) == "SBS4"|
                      colnames(tumor_pre_signature) == "SBS7a"|
                      colnames(tumor_pre_signature) == "SBS7b"|
                      colnames(tumor_pre_signature) == "SBS7c"|
                      colnames(tumor_pre_signature) == "SBS7d"|
                      colnames(tumor_pre_signature) == "SBS9"|
                      colnames(tumor_pre_signature) == "SBS10a"|
                      colnames(tumor_pre_signature) == "SBS10b"|
                      colnames(tumor_pre_signature) == "SBS11"|
                      colnames(tumor_pre_signature) == "SBS14"|
                      colnames(tumor_pre_signature) == "SBS15"|
                      colnames(tumor_pre_signature) == "SBS16"|
                      colnames(tumor_pre_signature) == "SBS17a"|
                      colnames(tumor_pre_signature) == "SBS17b"|
                      colnames(tumor_pre_signature) == "SBS19"|
                      colnames(tumor_pre_signature) == "SBS20"|
                      colnames(tumor_pre_signature) == "SBS21"|
                      colnames(tumor_pre_signature) == "SBS22"|
                      colnames(tumor_pre_signature) == "SBS23"|
                      colnames(tumor_pre_signature) == "SBS24"|
                      colnames(tumor_pre_signature) == "SBS25"|
                      colnames(tumor_pre_signature) == "SBS26"|
                      colnames(tumor_pre_signature) == "SBS27"|
                      colnames(tumor_pre_signature) == "SBS28"|
                      colnames(tumor_pre_signature) == "SBS29"|
                      colnames(tumor_pre_signature) == "SBS30"|
                      colnames(tumor_pre_signature) == "SBS31"|
                      colnames(tumor_pre_signature) == "SBS32"|
                      colnames(tumor_pre_signature) == "SBS34"|
                      colnames(tumor_pre_signature) == "SBS35"|
                      colnames(tumor_pre_signature) == "SBS36"|
                      colnames(tumor_pre_signature) == "SBS38"|
                      colnames(tumor_pre_signature) == "SBS42"|
                      colnames(tumor_pre_signature) == "SBS43"|
                      colnames(tumor_pre_signature) == "SBS44"|
                      colnames(tumor_pre_signature) == "SBS45"|
                      colnames(tumor_pre_signature) == "SBS46"|
                      colnames(tumor_pre_signature) == "SBS47"|
                      colnames(tumor_pre_signature) == "SBS48"|
                      colnames(tumor_pre_signature) == "SBS49"|
                      colnames(tumor_pre_signature) == "SBS50"|
                      colnames(tumor_pre_signature) == "SBS51"|
                      colnames(tumor_pre_signature) == "SBS52"|
                      colnames(tumor_pre_signature) == "SBS53"|
                      colnames(tumor_pre_signature) == "SBS54"|
                      colnames(tumor_pre_signature) == "SBS55"|
                      colnames(tumor_pre_signature) == "SBS56"|
                      colnames(tumor_pre_signature) == "SBS57"|
                      colnames(tumor_pre_signature) == "SBS58"|
                      colnames(tumor_pre_signature) == "SBS59"|
                      colnames(tumor_pre_signature) == "SBS60"|
                      colnames(tumor_pre_signature) == "SBS84"|
                      colnames(tumor_pre_signature) == "SBS85"|
                      colnames(tumor_pre_signature) == "SBS88")


tumor_pre_signature <- subset(tumor_pre_signature,
                                select = -remove_col)

tumor_pri_signature <- subset(tumor_pri_signature,
                               select = -remove_col)

library(stringr)
colnames(tumor_pre_signature) <- str_replace(colnames(tumor_pre_signature), "SBS", "")
colnames(tumor_pri_signature) <- str_replace(colnames(tumor_pri_signature), "SBS", "")

library(reshape2)
#Need to convert data frame from wide to long format
tumor_pre_signature <- as.data.frame(t(as.matrix(tumor_pre_signature)))
tumor_pre_signature$group <- row.names(tumor_pre_signature)
tumor_pre_signature.m <- melt(tumor_pre_signature, id.vars = "group")

tumor_pri_signature <- as.data.frame(t(as.matrix(tumor_pri_signature)))
tumor_pri_signature$group <- row.names(tumor_pri_signature)
tumor_pri_signature.m <- melt(tumor_pri_signature, id.vars = "group")

#Ordering the signatures for ggplot2
tumor_pre_signature.m$group <- factor(tumor_pre_signature.m$group, levels=unique(as.character(tumor_pre_signature.m$group)))
tumor_pri_signature.m$group <- factor(tumor_pri_signature.m$group, levels=unique(as.character(tumor_pri_signature.m$group)))


save(tumor_pre_signature.m, file="ESCC_tumor_pre_signature.RData")
save(tumor_pri_signature.m, file="ESCC_tumor_pri_signature.RData")

tumor_pre_signature_plot <- ggplot(data=tumor_pre_signature.m, aes(x=group, y=value, fill=group)) +
  geom_violin(scale="width") +
  geom_violin(scale="width") +
  theme(plot.title = element_text(hjust = 0), axis.text.x = element_text (hjust = 0.5)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("COSMIC Signature") +
  ylab("Signature Weight") + ylim(0, 1) + scale_x_discrete() +
  ggtitle("Normal")
tumor_pre_signature_plot


tumor_pri_signature_plot <- ggplot(data=tumor_pri_signature.m, aes(x=group, y=value, fill=group)) +
  geom_violin(scale="width") +
  theme(plot.title = element_text(hjust = 0), axis.text.x = element_text (hjust = 0.5)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  xlab("COSMIC Signature") +
  ylab("Signature Weight") + ylim(0, .7) +
  ggtitle("Primary")
tumor_pri_signature_plot
#combining all three plots
ggsave("ESCC_figures/Pri_ESCC_signature.png", width = 4.5, height = 3)

combined_ESCC_signature <- plot_grid(tumor_pre_signature_plot, tumor_pri_signature_plot,
                                     ncol=1, nrow=2, labels = c("A", "B"), label_size = 12)
combined_ESCC_signature
ggsave("ESCC_figures/combined_ESCC_signature.png", width = 9, height = 6)
