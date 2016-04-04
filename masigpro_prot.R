rm(list = ls())
library(maSigPro)
library(gtools)
library(DESeq2)
library(edgeR)
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
source('C:/Users/am4613/Documents/GitHub/Misc/gene_plotter.R')
library(gplots)
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')
setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/Masigpro_clustering_protein/')

prot_data <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/SQ_Results_PROTEIN.tsv', header = T, strings = F)
norm_data <- normalise_ProtDataset(prot_data, what = 'nada')
norm_data <- norm_data[,7:42]

heatmap_data <- normalise_ProtDataset(prot_data, what= 'heatmap')
heatmap_data <- heatmap_data[,7:42]

norm_data <- norm_data[which(rowSums(norm_data) != 0),]

exp_design = data.frame(Time = rep(0:11,each = 3), 
												Replicate = rep(1:3,12), 
												group = rep(1,ncol(norm_data)),
												row.names = colnames(norm_data))

design <- make.design.matrix(exp_design, degree = 3, time.col = 1, repl.col = 2, group.col = 3)
fit <- p.vector(norm_data, design)
tstep <- T.fit(fit)
sigs <- get.siggenes(tstep, vars = 'all')
see.genes(sigs, edesign = exp_design , 
					group.cols = 3, show.fit = T, summary.mode = 'median',
					color.mode = 'gray')


masigpro_sig <- data.frame(ID = sigs$summary, pval = sigs$sig.genes$sig.pvalues$'p-value', beta0 =sigs$sig.genes$coefficients$beta0, beta1 = sigs$sig.genes$coefficients$betaTime, beta2 = sigs$sig.genes$coefficients$betaTime2,
													 beta3 = sigs$sig.genes$coefficients$betaTime3)
write.table(masigpro_sig, 'masigpro_proteins.txt', sep = '\t')

arbol <- gene_plotter(sigs$summary, what = 'Protein')

hc.rows <- hclust(dist(heatmap_data[which(row.names(heatmap_data) %in% sigs$summary),]))

k = 4

clusters<- cutree(hc.rows, k)

plot(hc.rows, cex = 0.5)
rect.hclust(hc.rows,k = k)
dev.copy(pdf, file  = 'Dendrogram.pdf')
dev.off()
pval <- 0.05

pdf(file  = 'C:/Users/am4613/Documents/Summaries_as_timecourses/Masigpro_clustering_protein/Results_masig.pdf')
for(i in 1:max(clusters))
{
	cl <- names(clusters[clusters == i])
	write.table(cl, file = paste0('C:/Users/am4613/Documents/Summaries_as_timecourses/Masigpro_clustering_protein/Cluster_',i,'.txt'), sep = '\t')
	
	go_results <- GOanalysis(cl, GOtable, all = 5123)
	go_results <- go_results[go_results[,2] > pval, ]
	write.table(go_results, file = paste0('GO_cluster',i,'.txt'), sep = '\t')
	
	gene_plotter(cl, what = 'Protein')
}
dev.off()

