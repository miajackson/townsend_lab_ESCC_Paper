library(cancereffectsizeR)
library(scales)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)

setwd("/Users/miajackson/Documents/Townsend Lab/ESCC_Paper")
#load("analysis_escc_10-21.Rdata")
load("Analysis.Pri.12.2.results.Rdata")

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

normal_jitter <- ggplot(normal_list, aes(x=gene, y=selection_intensity, color=progression)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1))+ 
  xlab("Gene or protein") + ylab("Scaled selection coefficient") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        legend.position = c(0.95, 0.85))+
  scale_color_discrete(name = "Stage", labels = c("Normal", "Primary")) +
  geom_vline(xintercept=seq(1.5, length(unique(normal_list$gene))-0.5, 1), 
             lwd=.5, colour="lightgrey") +
  scale_y_continuous(labels=scientific, #limits = c(-0.8e4, 1e5)
                     )

normal_jitter



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


primary_jitter<- ggplot(primary_list, aes(x=gene, y=selection_intensity, color=progression)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1))+ 
  xlab("Gene or protein") + ylab("Scaled selection coefficient") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        panel.border = element_blank(),
        legend.position = c(0.95, 0.85))+
  scale_color_discrete(name = "Stage", labels = c("Normal", "Primary")) +
  geom_vline(xintercept=seq(1.5, length(unique(primary_list$gene))-0.5, 1), 
             lwd=.5, colour="lightgrey") +
  scale_y_continuous(labels=scientific, #limits = c(-.24e4, 3e4), 
                     #breaks = c(0, 1e4, 1.5e4, 2e4, 3e4)
                     )

primary_jitter


########################################################################################
#Box plots

twostage_results_normal_recur <- twostage_results_normal_recur[order(-twostage_results_normal_recur$selection_intensity),]
twostage_results_primary_recur <- twostage_results_primary_recur[order(-twostage_results_primary_recur$selection_intensity),]

normal_top10 <- twostage_results_normal_recur[1:10,]
primary_top10 <- twostage_results_primary_recur[1:10,]
primary_top10[is.na(primary_top10$ci_low_95)]$ci_low_95 <- 0


unique_jitter_normal_all <- ggplot(data = normal_top10, aes(x = 1, y = selection_intensity, color=progression)) +
  geom_jitter(position=position_jitter(0.0, seed = 5))+
  #geom_text(aes(label=variant),hjust=-.3, vjust=.3, color = "black") +
  geom_text_repel(data = normal_top10, aes(label = variant), color = "black", direction="y", nudge_x=0.03) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none")+
  xlab("Top 10 variants") + ylab("") +
  scale_color_manual(values = "#F8766D") + 
  scale_x_continuous(limits = c(0.99, 1.1)) +
  #scale_y_continuous(labels=scientific, )
  scale_y_log10(labels=scientific) 

unique_jitter_normal_all



unique_jitter_primary_all<- ggplot(data = primary_top10, aes(x = 1, y = selection_intensity, color=progression)) +
  geom_jitter(position=position_jitter(0.0, seed = 5)) +
  #geom_text_repel(data = primary_top10,mapping=aes(x=ILE2, y=TE,label=primary_top10$CA),
  #  size=4, size=6, box.padding = unit(0.5, "lines")
  #)+
  geom_text_repel(data = primary_top10, aes(label = variant), color = "black", direction="y", nudge_x=0.03) +
  #geom_text(aes(label=variant), hjust=-.1, color = "black", size=3) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none")+
  xlab("Top 10 variants") + ylab("") +
  scale_color_manual(values = "#00BFC5") +
  
  scale_x_continuous(limits = c(0.99, 1.1))+
  #scale_y_continuous(labels=scientific, #limits = c(-.24e4, 3e4), 
                     #breaks = c(0, 1e4, 1.5e4, 2e4, 3.0e4)
                     #)
  scale_y_log10(labels=scientific) 

unique_jitter_primary_all



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


normal_combined <- plot_grid(normal_jitter, unique_jitter_normal_all, 
                             align = "h", axis = "t", nrow = 1, ncol = 2, scale = 1, rel_widths = c(2,1),
                             labels = c("A", "B"), label_size = 10)
normal_combined_title <- plot_grid(normal_title, normal_combined, 
                                   align = "h", axis = "t", nrow = 2, ncol = 1, scale = 1, rel_heights = c(0.1,1))

normal_combined_title
ggsave("ESCC_Figures/normal_combined_title.3.23.png", width = 12, height = 4)

primary_combined <- plot_grid(primary_jitter, unique_jitter_primary_all, 
                              align = "h", axis = "t", nrow = 1, ncol = 2, scale = 1, rel_widths = c(2,1),
                              labels = c("C", "D"), label_size = 10)
primary_combined_title <- plot_grid(primary_title, primary_combined, 
                                    align = "h", axis = "t", nrow = 2, ncol = 1, scale = 1, rel_heights = c(0.1,1))

primary_combined_title
ggsave("ESCC_Figures/primary_combined_title.3.23.png", width = 12, height = 4)


bar_box_combined <- plot_grid(normal_combined_title, primary_combined_title,
                              align = "h", axis = "t", nrow = 2, ncol = 1, scale = 1)

bar_box_combined
ggsave("ESCC_Figures/bar_box_combined.3.23.png", width = 10, height = 9)


