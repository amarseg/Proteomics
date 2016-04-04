rm(list = ls())
library('betr')
library(gtools)
library(DESeq2)
library(edgeR)
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
library(gplots)
library(Biobase)
setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/betr_clustering_protein/')
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')
source('C:/Users/am4613/Documents/GitHub/Misc/gene_plotter.R')

prot_data <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/SQ_Results_PROTEIN.tsv', header = T, strings = F)
norm_data <- normalise_ProtDataset(prot_data, what = 'nada')
norm_data <- norm_data[,7:42]

heatmap_data <- normalise_ProtDataset(prot_data, what= 'heatmap')
heatmap_data <- heatmap_data[,7:42]

exp_design = data.frame(Time = rep(0:11,each = 3), 
												Replicate = rep(1:3,12), 
												group = rep(1,ncol(norm_data)),
												row.names = colnames(norm_data))

es_norm <- ExpressionSet(as.matrix(norm_data))

prob <- betr(eset = es_norm, timepoint = exp_design$Time, 
						 replicate = exp_design$Replicate, alpha = 0.05, twoCondition = FALSE)

prob_cutoff <- 0.99

sigs <- prob[prob > prob_cutoff]

write.table(sigs,'results_prot_betr.txt', sep = '\t')

tree_sig <- gene_plotter(names(sigs), what = 'Protein')

hc.rows <- hclust(dist(heatmap_data[which(row.names(heatmap_data) %in% names(sigs)),]))
plot(hc.rows)
k = 5
rect.hclust(hc.rows, k = k)


pval <- 0.05
clusters<- cutree(hc.rows, k)


pdf(file  = 'Results_betr_prot.pdf')
for(i in 1:max(clusters))
{
	cl <- names(clusters[clusters == i])
	write.table(cl, file = paste0('Cluster_',i,'.txt'), sep = '\t')
	
	go_results <- GOanalysis(cl, GOtable, all = 5123)
	go_results <- go_results[go_results[,2] > pval, ]
	write.table(go_results, file = paste0('GO_cluster',i,'.txt'), sep = '\t')
	
	gene_plotter(cl, what = 'Protein')
}
dev.off()

