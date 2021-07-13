library(cancereffectsizeR)
library(ggplot2)
library(stringr)
library(scales)


scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", label_scientific()(x)))))
}

setwd("/Users/miajackson/Documents/Townsend Lab/ESCC_Paper")
load("Analysis.5.7.ep.results.all.Rdata")
#Load in your epistasis data and change the object name

original <- ep_results
ESCC_epistasis <- epistasis_results(original)$gene_epistasis_1

ESCC_epistasis

######################################################################

# Combined Graphs

######################################################################

###NOTCH1###
ESCC_epistasis
#Select your gene pairs of interest
NOTCH1_list <- c(3,10,16, 22, 23, 24, 25, 26)

epistatic_change_NOTCH1 <- c()

#Decoupling the gene pairs
for(x in NOTCH1_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(ESCC_epistasis[x,5] - ESCC_epistasis[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(ESCC_epistasis[x,6] - ESCC_epistasis[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(ESCC_epistasis[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(ESCC_epistasis[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(ESCC_epistasis[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(ESCC_epistasis[x,2]))
  epistatic_change_NOTCH1 <- rbind(epistatic_change_NOTCH1, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_NOTCH1 <- data.frame(gene = epistatic_change_NOTCH1[,1], change = as.numeric(epistatic_change_NOTCH1[,2]))

#Separating "Before" and "After"
epistatic_change_NOTCH1_before <- epistatic_change_NOTCH1[grep("NOTCH1_", epistatic_change_NOTCH1[,1]),]

epistatic_change_NOTCH1_before$time <- rep("Before", 8)
epistatic_change_NOTCH1_before <- epistatic_change_NOTCH1_before[order(epistatic_change_NOTCH1_before$change),]
epistatic_change_NOTCH1_after <- epistatic_change_NOTCH1[grep("_NOTCH1", epistatic_change_NOTCH1[,1]),]
epistatic_change_NOTCH1_after$time <- rep("After", 8)
epistatic_change_NOTCH1_after <- epistatic_change_NOTCH1_after[order(epistatic_change_NOTCH1_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_NOTCH1_before[,1] <- sub("NOTCH1_after_", "", epistatic_change_NOTCH1_before[,1])
epistatic_change_NOTCH1_after[,1] <- sub("after_NOTCH1", "", epistatic_change_NOTCH1_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = 0, time = "Before")

epistatic_change_NOTCH1 <- rbind(epistatic_change_NOTCH1_before, blank, epistatic_change_NOTCH1_after)
epistatic_change_NOTCH1

epistatic_change_NOTCH1$gene <- factor(epistatic_change_NOTCH1$gene, levels = c("ERBB4", "TP53", "EP300",
                                                                                "BRAF", "PPM1D", "SETD2","PTEN",
                                                                                "NRAS", "BLANK", "EP300_",
                                                                                 "PPM1D_","SETD2_","ERBB4_","BRAF_",
                                                                                "TP53_","PTEN_", "NRAS_"))

gene_labels_NOTCH1 <- c(epistatic_change_NOTCH1_before$gene, "", epistatic_change_NOTCH1_after$gene)
gene_labels_NOTCH1 <- sub("_", "", gene_labels_NOTCH1)
epistatic_change_NOTCH1

rownames(epistatic_change_NOTCH1) <- 1:15
epistatic_change_NOTCH1$time_factor <- as.factor(epistatic_change_NOTCH1$time)
epistatic_change_NOTCH1

b1 <- bracketsGrob(0.33, 0.05, 0, 0.05, h=0.05, lwd=2, col="red")

waterfall_NOTCH1 <- ggplot(epistatic_change_NOTCH1, aes(x=gene, y=change, fill=time)) +
  geom_bar(stat = "identity") + theme_classic() + #scale_fill_discrete(breaks = c("Before", "After")) +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = .5),
        legend.position = "none", legend.title = element_blank()) +
  ggtitle("NOTCH1 gene pairs") + xlab("Gene") + ylab("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels_NOTCH1) +
  scale_y_continuous(labels = scientific, limits=c(-50000, NA))

waterfall_NOTCH1

ggsave("./ESCC_Figures/epistasis/waterfall_NOTCH1.png", width = 10, dpi=300, height = 7)

##ERBB4##
ERBB4_list <- c(8,9,11,12)

epistatic_change_ERBB4 <- c()

#Decoupling the gene pairs
for(x in ERBB4_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(ESCC_epistasis[x,5] - ESCC_epistasis[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(ESCC_epistasis[x,6] - ESCC_epistasis[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(ESCC_epistasis[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(ESCC_epistasis[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(ESCC_epistasis[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(ESCC_epistasis[x,2]))
  epistatic_change_ERBB4 <- rbind(epistatic_change_ERBB4, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_ERBB4 <- data.frame(gene = epistatic_change_ERBB4[,1], change = as.numeric(epistatic_change_ERBB4[,2]))


#Separating "Before" and "After"
epistatic_change_ERBB4_before <- epistatic_change_ERBB4[grep("ERBB4_", epistatic_change_ERBB4[,1]),]

epistatic_change_ERBB4_before$time <- rep("Before", 4)
epistatic_change_ERBB4_before <- epistatic_change_ERBB4_before[order(epistatic_change_ERBB4_before$change),]
epistatic_change_ERBB4_after <- epistatic_change_ERBB4[grep("_ERBB4", epistatic_change_ERBB4[,1]),]
epistatic_change_ERBB4_after$time <- rep("After", 4)
epistatic_change_ERBB4_after <- epistatic_change_ERBB4_after[order(epistatic_change_ERBB4_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_ERBB4_before[,1] <- sub("ERBB4_after_", "", epistatic_change_ERBB4_before[,1])
epistatic_change_ERBB4_after[,1] <- sub("after_ERBB4", "", epistatic_change_ERBB4_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = 0, time = "Before")

epistatic_change_ERBB4 <- rbind(epistatic_change_ERBB4_before, blank, epistatic_change_ERBB4_after)

epistatic_change_ERBB4$gene <- factor(epistatic_change_ERBB4$gene, levels = c("NRAS", "PTEN", "NOTCH1",
                                                                              "SETD2", "BLANK", "NRAS_",
                                                                              "PTEN_", "NOTCH1_", "SETD2_"))

gene_labels_ERBB4 <- c(epistatic_change_ERBB4_before$gene, "", epistatic_change_ERBB4_after$gene)
gene_labels_ERBB4 <- sub("_", "", gene_labels_ERBB4)

epistatic_change_ERBB4

waterfall_ERBB4 <- ggplot(epistatic_change_ERBB4, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = .5),legend.position = "none", legend.title = element_blank()) +
  ggtitle("ERBB4 gene pairs") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels_ERBB4) + 
  scale_y_continuous(labels = scientific, limits=c(-35000, 30000))

waterfall_ERBB4



ggsave("./ESCC_Figures/epistasis/waterfall_ERBB4.png", width = 10, dpi=300, height = 7)


##TP53
TP53_list <- c(8,15,21,26,33,35,36)
ESCC_epistasis

epistatic_change_TP53 <- c()
ESCC_epistasis
#Decoupling the gene pairs
for(x in TP53_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(ESCC_epistasis[x,5] - ESCC_epistasis[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(ESCC_epistasis[x,6] - ESCC_epistasis[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(ESCC_epistasis[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(ESCC_epistasis[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(ESCC_epistasis[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(ESCC_epistasis[x,2]))
  epistatic_change_TP53 <- rbind(epistatic_change_TP53, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_TP53 <- data.frame(gene = epistatic_change_TP53[,1], change = as.numeric(epistatic_change_TP53[,2]))


#Separating "Before" and "After"
epistatic_change_TP53_before <- epistatic_change_TP53[grep("TP53_", epistatic_change_TP53[,1]),]

epistatic_change_TP53_before$time <- rep("Before", 7)
epistatic_change_TP53_before <- epistatic_change_TP53_before[order(epistatic_change_TP53_before$change),]
epistatic_change_TP53_after <- epistatic_change_TP53[grep("_TP53", epistatic_change_TP53[,1]),]
epistatic_change_TP53_after$time <- rep("After", 7)
epistatic_change_TP53_after <- epistatic_change_TP53_after[order(epistatic_change_TP53_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_TP53_before[,1] <- sub("TP53_after_", "", epistatic_change_TP53_before[,1])
epistatic_change_TP53_after[,1] <- sub("after_TP53", "", epistatic_change_TP53_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = 0, time = "Before")

epistatic_change_TP53 <- rbind(epistatic_change_TP53_before, blank, epistatic_change_TP53_after)
epistatic_change_TP53
epistatic_change_TP53$gene <- factor(epistatic_change_TP53$gene, levels = c("ERBB4", "BRAF","EP300","PPM1D","SETD2","PTEN","NOTCH1","BLANK", "PPM1D_",
                                                                              "EP300_", "NOTCH1_", "SETD2_", "ERBB4_", "PTEN_", "BRAF_"))

gene_labels_TP53 <- c(epistatic_change_TP53_before$gene, "", epistatic_change_TP53_after$gene)
gene_labels_TP53 <- sub("_", "", gene_labels_TP53)

epistatic_change_TP53

waterfall_TP53 <- ggplot(epistatic_change_TP53, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = .5),legend.position = "none", legend.title = element_blank()) +
  ggtitle("TP53 gene pairs") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels_TP53) + 
  scale_y_continuous(labels = scientific, limits=c(-22000, 22000))

waterfall_TP53


ggsave("./ESCC_Figures/epistasis/waterfall_TP53.png", width = 10, dpi=300, height = 7)





library(cowplot)

epistasis_all <- plot_grid(
  waterfall_NOTCH1, waterfall_TP53,
  align="h", axis="t", nrow=2, ncol=1, scale = 1, labels = c("A","B"), label_size = 10)

epistasis_all

ggsave("./ESCC_Figures/epistasis/waterfall_combined.png", width = 6, dpi=300, height = 9)

######################################################################

# Precancerous Graphs

######################################################################


setwd("/Users/miajackson/Documents/Townsend Lab/ESCC_Paper")
load("Analysis.5.7.ep.results.pre.Rdata")
original_pre <- analysis_escc_ep_pre
ESCC_epistasis_pre <- epistasis_results(original_pre)$gene_epistasis_1

###NOTCH1###
ESCC_epistasis_pre
NOTCH1_list <- c(2,7,10,15,16,17,18,19)

epistatic_change_NOTCH1 <- c()

#Decoupling the gene pairs
for(x in NOTCH1_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(ESCC_epistasis_pre[x,5] - ESCC_epistasis_pre[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(ESCC_epistasis_pre[x,6] - ESCC_epistasis_pre[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(ESCC_epistasis_pre[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(ESCC_epistasis_pre[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(ESCC_epistasis_pre[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(ESCC_epistasis_pre[x,2]))
  epistatic_change_NOTCH1 <- rbind(epistatic_change_NOTCH1, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_NOTCH1 <- data.frame(gene = epistatic_change_NOTCH1[,1], change = as.numeric(epistatic_change_NOTCH1[,2]))


#Separating "Before" and "After"
epistatic_change_NOTCH1_before <- epistatic_change_NOTCH1[grep("NOTCH1_", epistatic_change_NOTCH1[,1]),]

epistatic_change_NOTCH1_before$time <- rep("Before", 8)
epistatic_change_NOTCH1_before <- epistatic_change_NOTCH1_before[order(epistatic_change_NOTCH1_before$change),]
epistatic_change_NOTCH1_after <- epistatic_change_NOTCH1[grep("_NOTCH1", epistatic_change_NOTCH1[,1]),]
epistatic_change_NOTCH1_after$time <- rep("After", 8)
epistatic_change_NOTCH1_after <- epistatic_change_NOTCH1_after[order(epistatic_change_NOTCH1_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_NOTCH1_before[,1] <- sub("NOTCH1_after_", "", epistatic_change_NOTCH1_before[,1])
epistatic_change_NOTCH1_after[,1] <- sub("after_NOTCH1", "", epistatic_change_NOTCH1_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = 0, time = "Before")

epistatic_change_NOTCH1 <- rbind(epistatic_change_NOTCH1_before, blank, epistatic_change_NOTCH1_after)
epistatic_change_NOTCH1
epistatic_change_NOTCH1$gene <- factor(epistatic_change_NOTCH1$gene, levels = c("ERBB4", "EP300", "BRAF",
                                                                              "SETD2", "PPM1D", "TP53",
                                                                              "PTEN", "NRAS", "BLANK", "EP300_",
                                                                              "PPM1D_", "ERBB4_", "SETD2_",
                                                                              "TP53_", "BRAF_", "PTEN_", "NRAS_"))


gene_labels_NOTCH1_pre <- c(epistatic_change_NOTCH1_after$gene, "", epistatic_change_NOTCH1_before$gene)
gene_labels_NOTCH1_pre <- sub("_", "", gene_labels_NOTCH1_pre)


waterfall_NOTCH1_pre <- ggplot(epistatic_change_NOTCH1, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = .5),
        legend.position = "none", legend.title = element_blank()) +
  ggtitle("NOTCH1 gene pairs (Normal)") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels_NOTCH1_pre) +
  scale_y_continuous(labels = scientific, limits=c(-380000, NA))

waterfall_NOTCH1_pre
ggsave("./ESCC_Figures/epistasis/waterfall_NOTCH1_pre.png", width = 10, dpi=300, height = 7)



###ERBB4###
#skipping for now
ERBB4_list <- c(1,2,3,4,5)

epistatic_change_ERBB4 <- c()

#Decoupling the gene pairs
for(x in ERBB4_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(ESCC_epistasis_pre[x,5] - ESCC_epistasis_pre[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(ESCC_epistasis_pre[x,6] - ESCC_epistasis_pre[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(ESCC_epistasis_pre[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(ESCC_epistasis_pre[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(ESCC_epistasis_pre[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(ESCC_epistasis_pre[x,2]))
  epistatic_change_ERBB4 <- rbind(epistatic_change_ERBB4, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_ERBB4 <- data.frame(gene = epistatic_change_ERBB4[,1], change = as.numeric(epistatic_change_ERBB4[,2]))


#Separating "Before" and "After"
epistatic_change_ERBB4_before <- epistatic_change_ERBB4[grep("ERBB4_", epistatic_change_ERBB4[,1]),]

epistatic_change_ERBB4_before$time <- rep("Before", 5)
epistatic_change_ERBB4_before <- epistatic_change_ERBB4_before[order(-epistatic_change_ERBB4_before$change),]
epistatic_change_ERBB4_after <- epistatic_change_ERBB4[grep("_ERBB4", epistatic_change_ERBB4[,1]),]
epistatic_change_ERBB4_after$time <- rep("After", 5)
epistatic_change_ERBB4_after <- epistatic_change_ERBB4_after[order(-epistatic_change_ERBB4_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_ERBB4_before[,1] <- sub("ERBB4_after_", "", epistatic_change_ERBB4_before[,1])
epistatic_change_ERBB4_after[,1] <- sub("after_ERBB4", "", epistatic_change_ERBB4_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = 0, time = "Before")

epistatic_change_ERBB4 <- rbind(epistatic_change_ERBB4_before, blank, epistatic_change_ERBB4_after)
epistatic_change_ERBB4
epistatic_change_ERBB4$gene <- factor(epistatic_change_ERBB4$gene, levels = c("SETD2", "NOTCH1", "TP53",
                                                                              "PTEN", "NRAS", "BLANK",
                                                                              "SETD2_", "NOTCH1_", "TP53_",
                                                                              "PTEN_", "NRAS_"))

gene_labels <- c(epistatic_change_ERBB4_after$gene, "", epistatic_change_ERBB4_before$gene)
gene_labels <- sub("_", "", gene_labels)


waterfall_ERBB4_pre <- ggplot(epistatic_change_ERBB4, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = .5),
        legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("ERBB4 gene pairs (Normal)") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels) + scale_fill_discrete(breaks = c("Before", "After"))+
  scale_y_continuous(labels = scientific) #+ scale_x_continuous(trans = pseudolog10_trans)

waterfall_ERBB4_pre

ggsave("./ESCC_Figures/epistasis/waterfall_ERBB4_pre.png", width = 10, dpi=300, height = 7)






###TP53###
ESCC_epistasis_pre
TP53_list <- c(6,9,14,19,22,23,25,26)

epistatic_change_TP53 <- c()

#Decoupling the gene pairs
for(x in TP53_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(ESCC_epistasis_pre[x,5] - ESCC_epistasis_pre[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(ESCC_epistasis_pre[x,6] - ESCC_epistasis_pre[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(ESCC_epistasis_pre[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(ESCC_epistasis_pre[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(ESCC_epistasis_pre[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(ESCC_epistasis_pre[x,2]))
  epistatic_change_TP53 <- rbind(epistatic_change_TP53, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_TP53 <- data.frame(gene = epistatic_change_TP53[,1], change = as.numeric(epistatic_change_TP53[,2]))


#Separating "Before" and "After"
epistatic_change_TP53_before <- epistatic_change_TP53[grep("TP53_", epistatic_change_TP53[,1]),]

epistatic_change_TP53_before$time <- rep("Before", 8)
epistatic_change_TP53_before <- epistatic_change_TP53_before[order(-epistatic_change_TP53_before$change),]
epistatic_change_TP53_after <- epistatic_change_TP53[grep("_TP53", epistatic_change_TP53[,1]),]
epistatic_change_TP53_after$time <- rep("After", 8)
epistatic_change_TP53_after <- epistatic_change_TP53_after[order(-epistatic_change_TP53_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_TP53_before[,1] <- sub("TP53_after_", "", epistatic_change_TP53_before[,1])
epistatic_change_TP53_after[,1] <- sub("after_TP53", "", epistatic_change_TP53_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = 0, time = "Before")

epistatic_change_TP53 <- rbind(epistatic_change_TP53_before, blank, epistatic_change_TP53_after)
epistatic_change_TP53
epistatic_change_TP53$gene <- factor(epistatic_change_TP53$gene, levels = c("PTEN", "NOTCH1", "EP300",
                                                                              "PPM1D", "SETD2", "BRAF",
                                                                              "ERBB4", "NRAS", "BLANK",
                                                                              "PTEN_", "BRAF_", "NOTCH1_",
                                                                            "SETD2_", "ERBB4_", "NRAS_",
                                                                            "PPM1D_", "EP300_"))
epistatic_change_TP53
gene_labels <- c(epistatic_change_TP53_after$gene, "", epistatic_change_TP53_before$gene)
gene_labels <- sub("_", "", gene_labels)

waterfall_TP53_pre <- ggplot(epistatic_change_TP53, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = .5),
        legend.position = "none", legend.title = element_blank()) +
  ggtitle("TP53 gene pairs (Normal)") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels) +
  scale_y_continuous(labels = scientific, limits=c(-30000, 30000))

waterfall_TP53_pre

ggsave("./ESCC_Figures/epistasis/waterfall_TP53_pre.5.8.png", width = 9, dpi=300, height = 6)


epistasis_pre <- plot_grid(
  waterfall_NOTCH1_pre, waterfall_TP53_pre,
  align="h", axis="t", nrow=2, ncol=1, scale = 1, labels = c("A","B"), label_size = 10)

epistasis_pre

ggsave("./ESCC_Figures/epistasis/waterfall_combined_pre.png", width = 6, dpi=300, height = 9)



######################################################################

# Primary Graphs

######################################################################

#primary tumor
load("Analysis.Pri.12.17.ep.results.Rdata")
original_pri <- ep_results_pri
ESCC_epistasis_pri <- epistasis_results(original_pri)$gene_epistasis_1
ESCC_epistasis_pri


##EP300##
EP300_list <- c(1,2)

epistatic_change_EP300 <- c()

#Decoupling the gene pairs
for(x in EP300_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(ESCC_epistasis_pri[x,5] - ESCC_epistasis_pri[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(ESCC_epistasis_pri[x,6] - ESCC_epistasis_pri[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(ESCC_epistasis_pri[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(ESCC_epistasis_pri[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(ESCC_epistasis_pri[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(ESCC_epistasis_pri[x,2]))
  epistatic_change_EP300 <- rbind(epistatic_change_EP300, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_EP300 <- data.frame(gene = epistatic_change_EP300[,1], change = as.numeric(epistatic_change_EP300[,2]))


#Separating "Before" and "After"
epistatic_change_EP300_before <- epistatic_change_EP300[grep("EP300_", epistatic_change_EP300[,1]),]

epistatic_change_EP300_before$time <- rep("Before", 2)
epistatic_change_EP300_before <- epistatic_change_EP300_before[order(-epistatic_change_EP300_before$change),]
epistatic_change_EP300_after <- epistatic_change_EP300[grep("_EP300", epistatic_change_EP300[,1]),]
epistatic_change_EP300_after$time <- rep("After", 2)
epistatic_change_EP300_after <- epistatic_change_EP300_after[order(-epistatic_change_EP300_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_EP300_before[,1] <- sub("EP300_after_", "", epistatic_change_EP300_before[,1])
epistatic_change_EP300_after[,1] <- sub("after_EP300", "", epistatic_change_EP300_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = 0, time = "Before")

epistatic_change_EP300 <- rbind(epistatic_change_EP300_before, blank, epistatic_change_EP300_after)
epistatic_change_EP300$gene <- factor(epistatic_change_EP300$gene, levels = c("TP53", "NOTCH1", "BLANK",
                                                                             "NOTCH1_", "TP53_"))
  
gene_labels <- c(epistatic_change_EP300_after$gene, "", epistatic_change_EP300_before$gene)
gene_labels <- sub("_", "", gene_labels)

waterfall_EP300_pri <- ggplot(epistatic_change_EP300, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom", legend.title = element_blank()) +
  ggtitle("EP300 gene pairs (Primary)") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels) + scale_fill_discrete(breaks = c("Before", "After"))+
  scale_y_continuous(labels = scientific #breaks = c(-1e4, 0, 1e4), #limits = c(-2e4, 2e4)
  ) #+ scale_x_continuous(trans = pseudolog10_trans)

waterfall_EP300_pri
ggsave("./ESCC_Figures/epistasis/waterfall_EP300_pri.png", width = 10, dpi=300, height = 7)



##NOTCH1##
NOTCH1_list <- c(1,3)

epistatic_change_NOTCH1 <- c()

#Decoupling the gene pairs
for(x in NOTCH1_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(ESCC_epistasis_pri[x,5] - ESCC_epistasis_pri[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(ESCC_epistasis_pri[x,6] - ESCC_epistasis_pri[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(ESCC_epistasis_pri[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(ESCC_epistasis_pri[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(ESCC_epistasis_pri[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(ESCC_epistasis_pri[x,2]))
  epistatic_change_NOTCH1 <- rbind(epistatic_change_NOTCH1, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_NOTCH1 <- data.frame(gene = epistatic_change_NOTCH1[,1], change = as.numeric(epistatic_change_NOTCH1[,2]))


#Separating "Before" and "After"
epistatic_change_NOTCH1_before <- epistatic_change_NOTCH1[grep("NOTCH1_", epistatic_change_NOTCH1[,1]),]

epistatic_change_NOTCH1_before$time <- rep("Before", 2)
epistatic_change_NOTCH1_before <- epistatic_change_NOTCH1_before[order(-epistatic_change_NOTCH1_before$change),]
epistatic_change_NOTCH1_after <- epistatic_change_NOTCH1[grep("_NOTCH1", epistatic_change_NOTCH1[,1]),]
epistatic_change_NOTCH1_after$time <- rep("After", 2)
epistatic_change_NOTCH1_after <- epistatic_change_NOTCH1_after[order(-epistatic_change_NOTCH1_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_NOTCH1_before[,1] <- sub("NOTCH1_after_", "", epistatic_change_NOTCH1_before[,1])
epistatic_change_NOTCH1_after[,1] <- sub("after_NOTCH1", "", epistatic_change_NOTCH1_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = 0, time = "Before")

epistatic_change_NOTCH1 <- rbind(epistatic_change_NOTCH1_before, blank, epistatic_change_NOTCH1_after)
epistatic_change_NOTCH1$gene <- factor(epistatic_change_NOTCH1$gene, levels = c("TP53", "EP300", "BLANK",
                                                                               "EP300_", "TP53_"))


gene_labels <- c(epistatic_change_NOTCH1_after$gene, "", epistatic_change_NOTCH1_before$gene)
gene_labels <- sub("_", "", gene_labels)


waterfall_NOTCH1_pri <- ggplot(epistatic_change_NOTCH1, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = .5),
        legend.position = "none", legend.title = element_blank()) +
  ggtitle("NOTCH1 gene pairs (Primary)") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels) +
  scale_y_continuous(labels = scientific, limits=c(NA, 10000)) #+ scale_x_continuous(trans = pseudolog10_trans)

waterfall_NOTCH1_pri
ggsave("./ESCC_Figures/epistasis/waterfall_NOTCH1_pri.png", width = 10, dpi=300, height = 7)




###TP53###


TP53_list <- c(2,3)

epistatic_change_TP53 <- c()

#Decoupling the gene pairs
for(x in TP53_list){
  gene1_after_gene2 <- unlist(c(as.character("gene1_after_gene2"), as.numeric(ESCC_epistasis_pri[x,5] - ESCC_epistasis_pri[x,3])))
  gene2_after_gene1 <- unlist(c(as.character("gene2_after_gene1"), as.numeric(ESCC_epistasis_pri[x,6] - ESCC_epistasis_pri[x,4])))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene1", as.character(ESCC_epistasis_pri[x,1]))
  gene1_after_gene2[1] <- str_replace(gene1_after_gene2[1], "gene2", as.character(ESCC_epistasis_pri[x,2]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene1", as.character(ESCC_epistasis_pri[x,1]))
  gene2_after_gene1[1] <- str_replace(gene2_after_gene1[1], "gene2", as.character(ESCC_epistasis_pri[x,2]))
  epistatic_change_TP53 <- rbind(epistatic_change_TP53, gene1_after_gene2, gene2_after_gene1)
}

epistatic_change_TP53 <- data.frame(gene = epistatic_change_TP53[,1], change = as.numeric(epistatic_change_TP53[,2]))


#Separating "Before" and "After"
epistatic_change_TP53_before <- epistatic_change_TP53[grep("TP53_", epistatic_change_TP53[,1]),]

epistatic_change_TP53_before$time <- rep("Before", 2)
epistatic_change_TP53_before <- epistatic_change_TP53_before[order(-epistatic_change_TP53_before$change),]
epistatic_change_TP53_after <- epistatic_change_TP53[grep("_TP53", epistatic_change_TP53[,1]),]
epistatic_change_TP53_after$time <- rep("After", 2)
epistatic_change_TP53_after <- epistatic_change_TP53_after[order(-epistatic_change_TP53_after$change),]

#Need to have extra underscore to have unique names
epistatic_change_TP53_before[,1] <- sub("TP53_after_", "", epistatic_change_TP53_before[,1])
epistatic_change_TP53_after[,1] <- sub("after_TP53", "", epistatic_change_TP53_after[,1])

#Blank spot
blank <- data.frame(gene = "BLANK", change = 0, time = "Before")

epistatic_change_TP53 <- rbind(epistatic_change_TP53_before, blank, epistatic_change_TP53_after)
epistatic_change_TP53$gene <- factor(epistatic_change_TP53$gene, levels = c("EP300", "NOTCH1", "BLANK",
                                                                                "EP300_", "NOTCH1_"))



gene_labels <- c(epistatic_change_TP53_after$gene, "", epistatic_change_TP53_before$gene)
gene_labels <- sub("_", "", gene_labels)

waterfall_TP53_pri <- ggplot(epistatic_change_TP53, aes(x= gene, y=change, fill=time)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() +
  scale_fill_manual(values = c("#F8766D","#00BFC4"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = .5),
        legend.position = "none", legend.title = element_blank()) +
  ggtitle("TP53 gene pairs (Primary)") + xlab("Gene") + ylab ("Epistatic change in selection")+
  scale_x_discrete(labels = gene_labels) +
  scale_y_continuous(labels = scientific, limits=c(-30000,NA)) #+ scale_x_continuous(trans = pseudolog10_trans)

waterfall_TP53_pri
ggsave("./ESCC_Figures/epistasis/waterfall_TP53_pri.png", width = 10, dpi=300, height = 7)



epistasis_pri <- plot_grid(
  waterfall_NOTCH1_pri,waterfall_TP53_pri,
  align="h", axis="t", nrow=2, ncol=1, scale = 1, labels = c("A","B"), label_size = 10)

epistasis_pri

ggsave("./ESCC_Figures/epistasis/waterfall_combined_pri.png", width = 6, dpi=300, height = 9)



