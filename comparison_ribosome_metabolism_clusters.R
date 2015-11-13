##load clusters proteomics
library('gplots')

load('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/clusters_prot.rda')

##See 301015_talk for referece to where these clusters are coming from 
clusters_metabolism <- c(f_clusters[[4]],f_clusters[[6]])
clusters_ribosomes <- c(f_clusters[[2]],f_clusters[[7]])

proteomics <- read.delim('C:/Users/am4613/Documents/normalised_proteomics.txt', header = T, strings = F)
proteomics <- proteomics[,7:42]

clusters_data <- proteomics[which(row.names(proteomics) %in% clusters_metabolism | row.names(proteomics) %in% clusters_ribosomes),]

separate <- list()

for(i in 1:3)
{
	separate[[i]] <- cbind(clusters_data[,seq(i, 33 + i,3)])
}


reodered <- cbind(separate[[1]], separate[[2]], separate[[3]])
norm <- reodered

norm[,1:12] <- log2(reodered[,1:12]/reodered[,1])
norm[,13:24] <- log2(reodered[,13:24]/reodered[,13])
norm[,25:36] <- log2(reodered[,25:36]/reodered[,25])

is.na(norm) <- sapply(norm, is.nan)
is.na(norm) <- sapply(norm, is.infinite)

norm[row.names(norm) %in% clusters_metabolism,37] <- 'red'
norm[row.names(norm) %in% clusters_ribosomes,37] <- 'blue'

norm <- norm[order(norm[,37]),]

heatmap.2(as.matrix(norm[,1:36]), trace = 'none', 
					RowSideColors = norm[,37], 
					col = colorRampPalette(c('blue','blue','grey','yellow','yellow')), Colv = F, Rowv = F, breaks = 100)

