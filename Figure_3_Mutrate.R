library(cancereffectsizeR)
library(scales)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(RColorBrewer)

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}


load("Analysis.5.7.results.Rdata")
twostage_final <- analysis_escc

mut_rate_stageless <- data.frame(gene=twostage_final@mutrates$gene,
                       pri_mutation_rate=twostage_final@mutrates$rate_grp_1)

mut_rate_stageless[mut_rate_stageless$gene == "BRAF",]
#applicable to both stageless and prim vs met
highlight <- rep(FALSE, length(mutation_rate$gene))
highlight[c(11128, 17456, 5551,5487,11237,13449, 13004)] <- TRUE

#####################################################################################################################

#scatter plot of gene-level mutation rates, prim vs met
selected_colors <- brewer.pal(n = length(highlight[highlight==TRUE]), name = 'Dark2')
selected_colors[2] <- "#027cd9"


mutrates_stageless_plot <- ggplot() +
  geom_jitter(data = mut_rate_stageless[!highlight,], aes(x=1, y= pri_mutation_rate, color = pri_mutation_rate), shape=16, position=position_jitter(0.05), size = .75, alpha=0.6) +
  scale_colour_gradient(low="#D0ECCD", high="#D0ECCD") +
  ylab("") + xlab("") + ggtitle("Mutation rates") +
  geom_jitter(data = mut_rate_stageless[highlight,], aes(x=1, y= pri_mutation_rate), shape=16, position=position_jitter(0.05, seed = 5), size = 3, color=selected_colors, alpha=1) +
  geom_text_repel(data = mut_rate_stageless[highlight,], aes(x=1, y= pri_mutation_rate, label = gene), position=position_jitter(0.05, seed = 5)) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(color = "black"),
        legend.position="none") +
  theme(plot.margin = margin(0,0,0,45, "pt")) +
  coord_flip() + scale_y_continuous(labels=scientific, limits=c(0.0000004,0.0000015)) #+ expand_limits(y=c(0, 0.000002))
mutrates_stageless_plot

ggsave("ESCC_Figures/ESCC_mutrate.5.8.png", width = 9, height = 3)

?xlab
