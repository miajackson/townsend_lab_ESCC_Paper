library(cancereffectsizeR)
library(scales)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)

setwd("/Users/miajackson/Documents/Townsend Lab/ESCC_Paper")
#load("analysis_escc_10-21.Rdata")
load("Analysis.5.7.results.Rdata")
selection_results <- analysis_escc@selection_results$`stage-specific`

notch1_selections <- selection_results[selection_results$variant_id %like% "NOTCH1",]


twostage_results <- snv_results(analysis_escc)
tmp = twostage_results$`stage-specific`

pre = tmp[, .(variant_name, ci_low_95_si_1, ci_high_95_si_1, variant_type, gene, si_1,maf_freq_in_Pre, progression = "Pre")]
pri = tmp[, .(variant_name, ci_low_95_si_2, ci_high_95_si_2, variant_type, gene, si_2,maf_freq_in_Pri, progression = "Pri")]

names(pre)[names(pre) == "si_1"] <- "selection_intensity"
names(pri)[names(pri) == "si_2"] <- "selection_intensity"
names(pre)[names(pre) == "maf_freq_in_Pre"] <- "tumors_with_variant"
names(pri)[names(pri) == "maf_freq_in_Pri"] <- "tumors_with_variant"
names(pre)[names(pre) == "ci_low_95_si_1"] <- "ci_low_95"
names(pri)[names(pri) == "ci_low_95_si_2"] <- "ci_low_95"
names(pre)[names(pre) == "ci_high_95_si_1"] <- "ci_high_95"
names(pri)[names(pri) == "ci_high_95_si_2"] <- "ci_high_95"

twostage_results <- rbind(pre, pri)

names(twostage_results)[names(twostage_results) == "variant_name"] <- "variant"

variant <- gsub("_", " ", twostage_results$variant)
twostage_results$variant <- variant


#only looks at amino acid changes
aac <- twostage_results$variant_type == "aac"
twostage_results <- twostage_results[aac,]

twostage_results <- twostage_results[order(-selection_intensity),]


all <- twostage_results$tumors_with_variant >= 1
twostage_results <- twostage_results[all,]
recurrent <- twostage_results$tumors_with_variant > 1
twostage_results_recur <- twostage_results[recurrent,]

normal <- twostage_results$progression == "Pre"
primary <- twostage_results$progression == "Pri"

twostage_results_normal <- twostage_results[normal,]
twostage_results_primary <- twostage_results[primary,]

normal <- twostage_results_recur$progression == "Pre"
primary <- twostage_results_recur$progression == "Pri"

twostage_results_normal_recur <- twostage_results_recur[normal,]
twostage_results_primary_recur <- twostage_results_recur[primary,]


##########################################
#Summary function

summary_multivar <- function(data, group_name) {
  data_clean <- data %>% 
    filter(group == group_name) %>%
    arrange(desc(selection_intensity)) %>%
    filter(selection_intensity > 1)
  
  # Summarise information of gene with multiple variants
  info1 <- data_clean %>% group_by(gene) %>%
    summarise(cum_si = sum(selection_intensity), # change sum to mean and sd
              mean_si = mean(selection_intensity),
              sd = sd(selection_intensity),
              max_si = max(selection_intensity),
              n_variant = n_distinct(variant)) %>%
    filter(n_variant > 1)
  
  top_variant <- data_clean %>%
    group_by(gene) %>% filter(row_number() == 1)
  
  merge_info <- merge(info1, top_variant[, -3], by.x = "gene") %>%
    arrange(desc(cum_si), desc(n_variant))
  return(merge_info)
}

###############################################################################################
#Gene level SI
################################################################################################
library(dplyr)

stage_data <- data.frame(variant = twostage_results$variant,
                         gene = twostage_results$gene,
                         selection_intensity = twostage_results$selection_intensity,
                         group = twostage_results$progression)


pre_info <- summary_multivar(stage_data, group = "Pre")
pri_info <- summary_multivar(stage_data, group = "Pri")



fill_na <- function(x, fill = 0) {
  x = ifelse(is.na(x), fill, x)
  return(x)
}

sd_info <- stage_data %>% 
  filter(selection_intensity > 1) %>%
  group_by(gene) %>% summarise(sd = sd(selection_intensity))


prim_info <- merge(pre_info, pri_info, by = "gene", all = T, 
                   suffixes = c(".e", ".l")) %>%
  mutate_at(c("cum_si.e", "cum_si.l", 
              "mean_si.e", "mean_si.l", 
              "sd.e", "sd.l",
              "n_variant.e", "n_variant.l"), fill_na) %>%
  mutate(n_variant_prim = n_variant.e + n_variant.l, 
         mean_si_prim = (cum_si.e + cum_si.l) / n_variant_prim) %>%
  arrange(desc(n_variant_prim))

stage_merge <- prim_info

########################################################################################
# Normal

stage_merge_normal_ordered <- stage_merge[order(-stage_merge$mean_si.e),]
selected_normal_genes <- stage_merge_normal_ordered$gene[1:10]

stage_merge_normal_ordered[stage_merge_normal_ordered$gene == "NOTCH1",]
head(stage_merge_normal_ordered, 10)

normal_list <- twostage_results %>% filter(gene %in% selected_normal_genes)

#  normal_list <- twostage_results %>%
#  filter(gene %in% selected_normal_genes) %>%
#  select(1,3,5,6,10)

#set order of genes for plot
normal_list$gene <- normal_list$gene %>%
  factor(levels = selected_normal_genes)

#Dummy points

dummy_PRKDC.l <- list("PRKDC Variant.l", as.double(0.001), "primary", "1", "PRKDC")
dummy_COL5A1.l <- list("COL5A1 Variant.l", as.double(0.001), "primary", "1", "COL5A1")
dummy_PTEN.l <- list("PTEN Variant.l", as.double(0.001), "primary", "1", "PTEN")
dummy_CNOT3.l <- list("CNOT3 Variant.l", as.double(0.001), "primary", "1", "CNOT3")

normal_list <- normal_list %>%
  rbind(dummy_PRKDC.l) %>%
  rbind(dummy_COL5A1.l) %>%
  rbind(dummy_PTEN.l) %>%
  rbind(dummy_CNOT3.l)

library(ggplot2)

normal_boxplot <- ggplot(normal_list, aes(x=gene, y=selection_intensity, fill=progression)) +
  #stat_boxplot(geom ='errorbar') +
  #geom_boxplot(outlier.size = 0) +
  geom_point(pch = 18, position = position_jitterdodge(jitter.width = .15),aes(color=progression)) +
  #geom_boxplot(position = "dodge2", outlier.shape = 1) + 
  xlab("Gene or protein") + ylab("Scaled selection coefficient") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        legend.position = c(0.92, .8))+
  #legend.position = c(0.8, 0.3)) 
  #scale_fill_discrete(name = "Stage", labels = c("Normal", "Primary")) +
  geom_vline(xintercept=seq(1.5, length(unique(normal_list$gene))-0.5, 1), lwd=.5, colour="lightgrey") +
  scale_y_continuous(labels=scientific, trans = 'log10') +
  scale_y_log10(labels=scientific)
normal_boxplot


########################################################################################
# Primary

stage_merge_primary_ordered <- stage_merge[order(-stage_merge$mean_si.l),]
selected_primary_genes <- stage_merge_primary_ordered$gene[1:10]


#select all variants within gene list
#primary_list <- twostage_results %>%
#  filter(gene %in% selected_primary_genes) %>%
#  select(1,3,5,6,10)
primary_list <- twostage_results %>% filter(gene %in% selected_primary_genes)

#set order of genes for plot
primary_list$gene <- primary_list$gene %>%
  factor(levels = selected_primary_genes)

#Dummy points

dummy_HRAS.e <- list("HRAS Variant.e", as.double(0.001), "normal", "1", "HRAS")
dummy_KMT2C.e <- list("KMT2C Variant.e", as.double(0.001), "normal", "1", "KMT2C")
dummy_CHRNA10.e <- list("CHRNA10 Variant.e", as.double(0.001), "normal", "1", "CHRNA10")
dummy_CACNA1E.e <- list("CACNA1E Variant.e", as.double(0.001), "normal", "1", "CACNA1E")

primary_list <- primary_list %>%
  rbind(dummy_HRAS.e) %>%
  rbind(dummy_KMT2C.e) %>%
  rbind(dummy_CHRNA10.e) %>%
  rbind(dummy_CACNA1E.e)

primary_list
library(ggplot2)
primary_boxplot <- ggplot(primary_list, aes(x=gene, y=selection_intensity, fill=progression)) +
  #stat_boxplot(geom ='errorbar') +
  geom_point(pch = 18, position = position_jitterdodge(jitter.width = .15),aes(color=progression)) +
  #geom_boxplot(position = "dodge2", outlier.shape = 1) + 
  xlab("Gene or protein") + ylab("Scaled selection coefficient") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        legend.position = c(0.92, .8))+
  #legend.position = c(.8,.3)) +
  #scale_fill_discrete(name = "Stage", labels = c("Normal", "Primary")) +
  geom_vline(xintercept=seq(1.5, length(unique(primary_list$gene))-0.5, 1), 
             lwd=.5, colour="lightgrey") +
  #scale_y_continuous(labels = scientific, breaks = trans_breaks("log10", function(x) 10^x))
  scale_y_log10(labels=scientific) 

primary_boxplot

##########################################################################
library(cowplot)

combined_boxplot <- plot_grid(normal_boxplot, primary_boxplot,
                              align = "h", axis = "t", nrow = 2, ncol = 1, scale = 1,
                              labels = c("A", "B"), label_size = 10)

combined_boxplot

ggsave("ESCC_Figures/epistasis/boxplot_twostage.4.19.png", width = 10, height = 12)

##############################################################################
# Variants
twostage_results_normal_recur <- twostage_results_normal_recur[order(-twostage_results_normal_recur$selection_intensity),]
twostage_results_primary_recur <- twostage_results_primary_recur[order(-twostage_results_primary_recur$selection_intensity),]

normal_top10 <- twostage_results_normal_recur[1:10,]
primary_top10 <- twostage_results_primary_recur[1:10,]
primary_top10[is.na(primary_top10$ci_low_95)]$ci_low_95 <- 0


normal_colors <- rep("grey85", 10)
normal_colors[c(3,5:10)] <- "#F8766D"


bargraph_normal <- ggplot(data=normal_top10, aes(x=reorder(variant, selection_intensity), y=selection_intensity, fill=reorder(variant, selection_intensity)))+
  #geom_bar(stat="identity") +
  geom_point(pch = 22, position = position_jitterdodge(jitter.width = .15), aes(color=progression)) +
  #geom_errorbar(aes(ymin=ci_low_95, ymax=ci_high_95), width=0.25) +
  theme(axis.line = element_line(color = 'black'), axis.ticks.y = element_blank(), panel.border = element_blank(),
        axis.text.y = element_text()) +
  theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.45)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin = margin(0,0,0,13, "pt")) +
  xlab("") + ylab("Scaled selection coefficient") +
  scale_y_continuous(labels=scientific, breaks = c(0,  1e6, 2e6), limits = c(-1e4, 2e6)) + coord_flip() + 
  scale_y_log10(labels=scientific) + coord_flip() +
  scale_fill_manual(values = normal_colors) +
  geom_text(aes(label=tumors_with_variant), position = position_stack(vjust=0), hjust = 1.4) 


bargraph_normal

primary_colors <- rep("grey85", 10)
primary_colors[c(1,2,4,5,8,9,10)] <- "#00BFC5"

#00BA38

bargraph_primary <- ggplot(data=primary_top10, aes(x=reorder(variant, selection_intensity), y=selection_intensity, fill=reorder(variant, selection_intensity)))+
  #geom_bar(stat="identity") + 
  geom_point(pch = 18, position = position_jitterdodge(jitter.width = .15), aes(color=progression)) +
  #geom_errorbar(aes(ymin=ci_low_95, ymax=ci_high_95), width=0.25) +
  theme(axis.line = element_line(color = 'black'), axis.ticks.y = element_blank(), panel.border = element_blank()) +
  theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.45)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(legend.position="none") +
  theme(plot.margin = margin(0,0,0,11, "pt")) +
  xlab("") + ylab("Scaled selection coefficient") +
  #scale_y_continuous(labels=scientific, breaks = c(0,  1e6, 2e6), limits = c(-1e4,2e6)) + coord_flip() + 
  scale_y_log10(labels=scientific) + coord_flip() +
  scale_fill_manual(values = primary_colors) +
  geom_text(aes(label=tumors_with_variant), position = position_stack(vjust=0), hjust = 1.4)

bargraph_primary


bargraph_combined <- plot_grid(bargraph_normal, bargraph_primary,
                               align = "h", axis = "t", nrow = 2, ncol = 1, scale = 1,
                               labels = c(), label_size = 10)

bargraph_combined
ggsave("New_Epistasis_Graphs_2_11/bargraph_twostage.12.4.png", width = 5, height = 10)


normal_title <- ggdraw() + 
  draw_label(
    "Normal Skin",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 49)
  )

primary_title <- ggdraw() + 
  draw_label(
    "Primary Tumor",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 49)
  )


normal_combined <- plot_grid(normal_boxplot, bargraph_normal, 
                             align = "h", axis = "t", nrow = 1, ncol = 2, scale = 1, rel_widths = c(2,1),
                             labels = c("A", "B"), label_size = 10)
normal_combined_title <- plot_grid(normal_title, normal_combined, 
                                   align = "h", axis = "t", nrow = 2, ncol = 1, scale = 1, rel_heights = c(0.1,1))

normal_combined_title
ggsave("New_Epistasis_Graphs_2_11/normal_combined_title.12.4.png", width = 12, height = 4)

primary_combined <- plot_grid(primary_boxplot, bargraph_primary, 
                              align = "h", axis = "t", nrow = 1, ncol = 2, scale = 1, rel_widths = c(2,1),
                              labels = c("C", "D"), label_size = 10)
primary_combined_title <- plot_grid(primary_title, primary_combined, 
                                    align = "h", axis = "t", nrow = 2, ncol = 1, scale = 1, rel_heights = c(0.1,1))

primary_combined_title
ggsave("New_Epistasis_Graphs_2_11/primary_combined_title.12.4.png", width = 12, height = 4)


bar_box_combined <- plot_grid(normal_combined_title, primary_combined_title,
                              align = "h", axis = "t", nrow = 2, ncol = 1, scale = 1)

bar_box_combined
ggsave("New_Epistasis_Graphs_2_11/bar_box_combined.12.4.png", width = 12, height = 9)


