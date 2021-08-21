library(cancereffectsizeR)
setwd("/Users/miajackson/Documents/Townsend Lab/ESCC_Paper/ESCC_Data")

combined_data_wes = "combined_data_wes.maf"
combined_data_wgs = "combined_data_wgs.maf"
martincorena_TS_data = "martincorena_TS_data.maf"
ucla_TS_data = "ucla_TS_data.maf"
yokoyama_TS_data = "yokoyama_TS_data.maf"
covered_regions_mart = "covered_regions_mart.bed"
covered_regions_yoko = "covered_regions_yoko.bed"
covered_regions_ucla = "covered_regions_ucla.bed"

#yokoyama_TS_data_table <- read.table(yokoyama_TS_data, sep = "\t", header = T, stringsAsFactors = F)
#head(yokoyama_TS_data_table)


analysis = CESAnalysis(refset = "ces.refset.hg19", sample_groups = c("Pre", "Pri"))

analysis = load_maf(analysis, maf = combined_data_wes,coverage="exome", group_col = "Pre_or_Pri")
analysis = load_maf(analysis, maf = combined_data_wgs,coverage="exome", group_col = "Pre_or_Pri")

analysis = load_maf(analysis, maf = yokoyama_TS_data, group_col = "Pre_or_Pri", coverage = "targeted", 
                    covered_regions_name = "TGS1", covered_regions = covered_regions_yoko)

analysis = load_maf(analysis, maf = martincorena_TS_data, group_col = "Pre_or_Pri", coverage = "targeted", 
                    covered_regions_name = "TGS2", covered_regions = covered_regions_mart)

analysis = load_maf(analysis, maf = ucla_TS_data, group_col = "Pre_or_Pri", coverage = "targeted", 
                    covered_regions_name = "TGS3", covered_regions = covered_regions_ucla)


analysis = trinuc_mutation_rates(analysis,signature_set = "COSMIC_v3.1")

analysis = gene_mutation_rates(analysis, covariates = "ESCA")

analysis = ces_variant(cesa = analysis, run_name = "sswm")
analysis = ces_variant(cesa = analysis, model = "sswm_sequential", run_name = "stage-specific",
                       groups = list("Pre", "Pri"))

analysis_escc <- analysis

setwd("/Users/miajackson/Documents/Townsend Lab/ESCC_Paper")

save(analysis_escc,file = "Analysis.5.7.results.Rdata")

genes_selected <- c("NOTCH1", "TP53","NRAS","BRAF","EP300", "ERBB4", "NRAS", 
                    "SETD2", "PPM1D", "HRAS", "PTEN")
ep_results <- ces_gene_epistasis(analysis, conf = .95, genes = genes_selected)

save(ep_results,file = "Analysis.5.7.ep.results.all.Rdata")

###########Only Precancerous###########
#running epistasis of only primary
pre_analysis <- analysis_escc
pre_analysis@samples <- pre_analysis@samples[group == "Pre"]
pre_analysis@maf = pre_analysis@maf[Unique_Patient_Identifier %in% pre_analysis@samples$Unique_Patient_Identifier]

selected_genes <- c("NOTCH1", "TP53","NRAS","BRAF","EP300", "ERBB4", "NRAS", 
                    "SETD2", "PPM1D", "HRAS", "PTEN")
analysis_escc_ep_pre <- ces_gene_epistasis(pre_analysis, conf = .95, genes = selected_genes)
save(analysis_escc_ep_pre,file="Analysis.5.7.ep.results.pre.Rdata")


###########Primary###########
pri_analysis <- analysis_escc
pri_analysis@samples <- pri_analysis@samples[group == "Pri"]
pri_analysis@maf = pri_analysis@maf[Unique_Patient_Identifier %in% pri_analysis@samples$Unique_Patient_Identifier]

selected_genes <- c("NOTCH1", "TP53","NRAS","BRAF","EP300", "ERBB4", "NRAS", 
                    "SETD2", "PPM1D", "HRAS", "PTEN")
analysis_escc_ep_pri <- ces_gene_epistasis(pri_analysis, conf = .95, genes = selected_genes)

save(analysis_escc_ep_pri,file="Analysis.5.7.ep.results.pri.Rdata")




###ignore this###
selection <- analysis_escc$selection$`stage-specific`

selection_notch1 <- selection[selection$gene == "NOTCH1",]
selection_notch1
