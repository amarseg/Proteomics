rm(list = ls())
library(maSigPro)
library(gtools)
library(DESeq2)
library(edgeR)
library(pheatmap)
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
source('C:/Users/am4613/Documents/GitHub/Misc/gene_plotter.R')
library(gplots)
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')
setwd('C:/Users/am4613/OneDrive/Summaries_as_timecourses/Masigpro_clustering_protein/')

prot_data <- read.delim('C:/Users/am4613/OneDrive - Imperial College London/ondedriveBACK/Summaries_as_timecourses/analysis/SQ_Results_PROTEIN.tsv', header = T, strings = F)
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


###Intersection between masigpro and BETR
###I load the data again because I don't want to run the whole MaSigPro
masig <- read.delim('masigpro_proteins.txt', header = T, strings = F)
betr <- read.delim('../betr_clustering_protein/results_prot_betr.txt', header = T, strings = F)
both <- intersect(masig$ID, row.names(betr))

t<- gene_plotter(both, what = 'Protein')
both_heatmap <- heatmap_data[which(row.names(heatmap_data) %in% both),]
both_heatmap<- reorder_proteomics(both_heatmap)
t<- pheatmap(both_heatmap, cluster_cols = F, cutree_rows = 3)
m <- cutree(t$tree_row, k = 3)
write.table(m, sep = '\t', 'both_protein_clusters.txt')

mRNA <- read.delim('../both_clustering/Cluster_1.txt', header = T, strings = F)
mRNA2 <- read.delim('../both_clustering/Cluster_2.txt', header = T, strings = F)
all_rna <- rbind(mRNA, mRNA2)

venn(list(all_rna[,1], both))

not_in_prot <- all_rna[which(!(all_rna[,1] %in% both)),]

venn(list(not_in_prot, row.names(norm_data)))

rna_seq <- read.delim('../DESeq_norm_counts.txt', header = T, strings = F)

all_data <- merge(rna_seq, norm_data, by.x = 'row.names', by.y = 'row.names')
all_data2 <- t(all_data)
corr <- c(1:ncol(all_data2))

for(i in 1:ncol(all_data2))
{
	corr[i] <- cor(as.numeric(all_data2[2:37,i]), as.numeric(all_data2[38:73,i]), use = 'complete.obs')
}
hist(corr, col = 'darkgrey')
plot(density(corr), col = 'darkgrey')

##Create joint dataset with normalised!!! data
rna_seq_norm <- normalise_rna(reorder_proteomics(rna_seq))
heatmap_data <- reorder_proteomics(heatmap_data)

all_norm <- merge(rna_seq_norm, heatmap_data, by.x ='row.names', by.y = 'row.names')
row.names(all_norm) <- all_norm[,1]
post_trans <- all_data2[1,which(corr < -0.2)]

post_trans_data <- all_norm[which(all_norm$Row.names %in% post_trans),]

col <- colorRampPalette(c('blue','gray','yellow'))(100)
cl_post <- pheatmap(post_trans_data[,2:73], cluster_cols = F, clustering_method = 'ward.D2', color = col, cutree_rows = 2, labels_row = post_trans_data[,1])

cl_names <- cutree(cl_post$tree_row, k = 2)
write.table(cl_names, 'post_transcriptional_regulation.txt', sep = '\t')

cluster_1 <- names(cl_names[which(cl_names == 1)])
cluster_2 <- names(cl_names[which(cl_names == 2)])

source('C:/Users/am4613/Documents/GitHub/Misc/half_life_plotter.R')
half_life_plotter(cluster_1, data = 'protein')
half_life_plotter(cluster_1, data = 'rna')
half_life_plotter(cluster_2, data = 'protein')
half_life_plotter(cluster_2, data = 'rna')
