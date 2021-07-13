library(cancereffectsizeR)
library(cowplot)
library(dplyr)
library(ggplot2)

setwd("/Users/miajackson/Documents/Townsend Lab/ESCC_Paper")
#load("analysis_escc_10-21.Rdata")
load("Analysis.5.7.results.Rdata")

ESCC_analysis <- analysis_escc


ESCC_trinuc <- ESCC_analysis@trinucleotide_mutation_weights$trinuc_proportion_matrix

#Seperating trinuc data into pre and pri
tumor_names <- unique(ESCC_analysis@maf$Unique_Patient_Identifier)

tumor_pre <- c()
tumor_pri <- c()

samples <- ESCC_analysis@samples

for (x in tumor_names){
  if (samples$group[which(samples$Unique_Patient_Identifier == x)] == "Pre"){
    tumor_pre <- c(tumor_pre, x)
  }
  else if (samples$group[which(samples$Unique_Patient_Identifier == x)] == "Pri"){
    tumor_pri <- c(tumor_pri, x)
  }
}

#---------------Added Stuff-------------#


tumor_stage_pre <- ESCC_trinuc[tumor_pre, ]
tumor_stage_pri <- ESCC_trinuc[tumor_pri, ]

#Averaging trinuc data
tumor_stage_pre_avg <- 
  apply(tumor_stage_pre, 2, mean)
tumor_stage_pri_avg <- 
  apply(tumor_stage_pri, 2, mean)

# Ordering trinuc average data
tumor_stage_pre_avg_ordered <- 
  data.frame(average = tumor_stage_pre_avg,
             mutation = names(tumor_stage_pre_avg)) %>%
  mutate(Upstream = substr(mutation, 1, 1),
         Downstream = substr(mutation, 7, 7),
         mutated_from = substr(mutation, 3, 3),
         mutated_to = substr(mutation, 5, 5),
         mutation = paste0(mutated_from, "\u2192", mutated_to)) %>%
  arrange(Downstream, Upstream)

tumor_stage_pri_avg_ordered <- 
  data.frame(average = tumor_stage_pri_avg,
             mutation = names(tumor_stage_pr_avg)) %>% ###there was a mistake there!!
  mutate(Upstream = substr(mutation, 1, 1),
         Downstream = substr(mutation, 7, 7),
         mutated_from = substr(mutation, 3, 3),
         mutated_to = substr(mutation, 5, 5),
         mutation = paste0(mutated_from, "\u2192", mutated_to)) %>%
  arrange(Downstream, Upstream)

trinuc.mutation_data <- tumor_stage_pre_avg_ordered[, -1]



# ---- Generate heatmap ----




tumor_trinuc_pre <- data.frame(matrix(ncol = 96, nrow = 0))
tumor_trinuc_pri <- data.frame(matrix(ncol = 96, nrow = 0))
colnames(tumor_trinuc_pri) <- colnames(ESCC_trinuc)

for (x in tumor_pre){
  tumor_trinuc_pre <- rbind(tumor_trinuc_pre, ESCC_trinuc[x,])
}
colnames(tumor_trinuc_pre) <- colnames(ESCC_trinuc)

for(x in tumor_pri){
  tumor_trinuc_pri <- rbind(tumor_trinuc_pri, ESCC_trinuc[x,])
}
colnames(tumor_trinuc_pri) <- colnames(ESCC_trinuc)

#Averaging trinuc data
tumor_trinuc_pre_avg <- apply(tumor_trinuc_pre, 2, mean)
tumor_trinuc_pri_avg <- apply(tumor_trinuc_pri, 2, mean)

tumor_trinuc_pre_avg_ordered <- c()
tumor_trinuc_pri_avg_ordered <- c()

#Reordering tumor_trinuc_pre_avg
#N_A
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
#N_C
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
#N_G
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
#N_T
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_pre)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_pre)[x], 7, 7)) == TRUE){
    tumor_trinuc_pre_avg_ordered <- c(tumor_trinuc_pre_avg_ordered, tumor_trinuc_pre_avg[x])
    
  }
}

#Reordering tumor_trinuc_pri_avg
#N_A
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("A", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
#N_C
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("C", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
#N_G
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("G", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
#N_T
for(x in 1:96){
  if (grepl("A", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("C", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("G", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
  }
}
for(x in 1:96){
  if (grepl("T", substr(colnames(tumor_trinuc_pri)[x], 1, 1)) & grepl("T", substr(colnames(tumor_trinuc_pri)[x], 7, 7)) == TRUE){
    tumor_trinuc_pri_avg_ordered <- c(tumor_trinuc_pri_avg_ordered, tumor_trinuc_pri_avg[x])
    
  }
}

ESCC_trinuc_heatmap_data <- data.frame(mutation=trinuc.mutation_data$mutation,
                                       Upstream=trinuc.mutation_data$Upstream,
                                       Downstream=trinuc.mutation_data$Downstream,
                                       mutated_from=trinuc.mutation_data$mutated_from,
                                       mutated_to=trinuc.mutation_data$mutated_to,
                                       trinuc_context=rep(NA, 96),
                                       trinuc_pre=tumor_trinuc_pre_avg_ordered,
                                       trinuc_pri=tumor_trinuc_pri_avg_ordered)


library(stringr)
ESCC_trinuc_heatmap_data$mutation <- str_replace(ESCC_trinuc_heatmap_data$mutation, "%->%", "\u2192")

for (x in 1:96){
  ESCC_trinuc_heatmap_data$trinuc_context[x] <- paste(ESCC_trinuc_heatmap_data$Upstream[x], 
                                                      ESCC_trinuc_heatmap_data$mutated_from[x],
                                                      ESCC_trinuc_heatmap_data$Downstream[x],
                                                      sep = "")  
}

levels(ESCC_trinuc_heatmap_data$mutation) <- c("C\u2192A", "C\u2192G", "C\u2192T", "T\u2192A", "T\u2192C", "T\u2192G")

# save(ESCC_trinuc_heatmap_data, file="ESCC_trinuc_heatmap_data.RData")
# 
# load("ESCC_trinuc_heatmap_data.RData")

ESCC_trinuc_bar_pre <- ggplot(data=ESCC_trinuc_heatmap_data, aes(x = trinuc_context, y= trinuc_pre*100, fill = mutation)) +
  geom_bar(stat = "identity") + theme_bw()+
  facet_wrap(.~mutation, nrow = 1, scales = "free_x")+
  theme(axis.line = element_line(color = 'white'), axis.ticks.y = element_blank(), panel.border = element_blank(),
        axis.text.y = element_text(), axis.text.x = element_text(angle=90, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())+ 
  theme(legend.position="none")+
  xlab("") + ylab("Mutation type probability")+ ggtitle("Normal")+
  scale_y_continuous(labels = function(x) paste0(x, "%"))

ESCC_trinuc_bar_pre
ggsave("ESCC_Figures/ESCC_trinuc_bar_pre.5.8.png", width=12, height=3)


ESCC_trinuc_bar_pri <- ggplot(data=ESCC_trinuc_heatmap_data, aes(x = trinuc_context, y= trinuc_pri*100, fill = mutation)) +
  geom_bar(stat = "identity") + theme_bw()+
  facet_wrap(.~mutation, nrow = 1, scales = "free_x")+
  theme(axis.line = element_line(color = 'white'), axis.ticks.y = element_blank(), panel.border = element_blank(),
        axis.text.y = element_text(), axis.text.x = element_text(angle=90, vjust = 0.5))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())+ 
  theme(legend.position="none")+
  xlab("") + ylab("Mutation type probability")+ ggtitle("Primary")+
  scale_y_continuous(labels = function(x) paste0(x, "%"))

ESCC_trinuc_bar_pri
ggsave("ESCC_Figures/ESCC_trinuc_bar_pri.5.8.png", width=12, height=3)


combined_trinuc_bar <-plot_grid(ESCC_trinuc_bar_pre, ESCC_trinuc_bar_pri,
                                labels = c("A", "B"), label_size = 12,
                                align="h", axis="t", nrow=2, ncol=1, rel_heights = c(1,1,1))

combined_trinuc_bar

ggsave("ESCC_Figures/combined_trinuc_bar.5.8.png",width=12, height=6)

