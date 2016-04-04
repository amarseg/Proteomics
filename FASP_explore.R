##Exploring FASP dataset
rm(list = ls())
library(plyr)
par(mfrow = c(1,1))
library(limma)
library(gplots)

MAplot <- function(d1,d2)
{
	d1 = d1 + 0.001 #Removing zeroes
	d2 = d2 + 0.001
	
	M <- log2(d1) - log2(d2)
	A <- 1/2*(log2(d1)+log2(d2))
	
	plot(x = A, y = M)
}

data <- read.delim('P:/max_quant_FASP_comparison.txt', header = T, strings = F)
data_pg <- read.delim('P:/progenesis_FASP.txt', header = T, strings = F, skip = 1)

ids <- strsplit(data[,1],split = '|', fixed = T)
ids_df <- data.frame(Reduce(rbind, ids))

row.names(data) <- ids_df[,2]
data <- data[,5:12]


design <- data.frame(row.names = colnames(data[,5:12]), 
										 Nitrogen_source = rep(c('Phenilanaline','NH4Cl','NH4Cl','Phenilalanine'),each =2),
										 Techinique = rep(c('FASP','Amalia'),each = 4))


col = colorRampPalette(c('blue','gray','yellow'))

