##AS proteomics analysis
library('gplots')
library('gtools')
library('Mfuzz')
library('matrixStats')
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')

as_safequant <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/SQ_Results_PROTEIN.tsv', header = T, strings = F)


##Preprocessing, and averaging between technical replicates.Reordering

numbers <- as_safequant[,1:78]
numbers <- numbers[grep(numbers[,1], pattern = 'SP'),]

titulos <- colnames(numbers)[7:78]
titulos <- substr(titulos, 0, nchar(titulos) - 3)
titulos <- unique(titulos)

norm_sq <- matrix(ncol = 6+36, nrow = nrow(numbers), NA)
norm_sq <- as.data.frame(norm_sq)
colnames(norm_sq) <- c(colnames(numbers[,1:6]), titulos)

norm_sq[,1:6] <- numbers[,1:6]

##Average between technical replicates

j = 7

for(i in seq(7,78,2))
{
	  norm_sq[,j] <- rowMeans(cbind(numbers[,i],numbers[,i+1]), na.rm = T)
	  
	  j = j +1
}

##Averaging between biological replicates for fuzzy clustering

median_sq <- matrix(ncol = 6+12, nrow = nrow(numbers), NA)
median_sq <- as.data.frame(median_sq)
colnames(median_sq) <- c(colnames(norm_sq[,1:6]),titles)

j = 7

for(i in seq(7,42,3))
{
	median_sq[,j] <- rowMedians(cbind(norm_sq[,i],norm_sq[,i+1],norm_sq[,i+2]), na.rm = T)
	
	j = j +1
}

row.names(median_sq) <- norm_sq[,1]
a <- strsplit(row.names(median_sq), "\\|")
vector_with_ids <- sapply(a,"[[",1)
row.names(median_sq) <- vector_with_ids
##Create Expression set for fuzzy 
eset <- ExpressionSet(as.matrix(median_sq[,7:18]))

eset.r <- filter.NA(eset,0.5) #Threshold for discarding missing data 
eset.f <- fill.NA(eset.r, mode = 'mean') #Fill the gaps with the mean of the other values

tmp <- filter.std(eset.f, min.std = 0)

eset.s <- standardise(tmp) #Data with average zero and sd = 1, so samples with similar changes are close, despite differeces in expression levels
n_cl = 9
m_fuzzy = mestimate(eset.s)
cl <- mfuzz(eset.s, c = n_cl, m = m_fuzzy)
mfuzz.plot(eset.s, cl = cl, mfrow = c(3,3), min.mem = 0.7, colo = rainbow(n=10))

##Code to extract the members of the different clusters and do GO analysis on them 

f_clusters <- list()

for(i in 1:n_cl)
{
	f_clusters[[i]] <- names(cl$cluster[cl$cluster == i])
}

##GO analysis of the different clusters
p = 0.05

results_go <- list()

for(i in 1:n_cl)
{
	temp <- GOanalysis(li = f_clusters[[i]], go = GOtable, all = 5123, sort = T)
	temp[,6] <- i 
	results_go[[i]] <- temp[temp[,2] <= p,]
}

output_final <- data.frame(Reduce(rbind,results_go))
write.table(output_final, 'C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/tableGO_fuzzy.txt', sep = '\t')
save(f_clusters,file = 'C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/clusters_prot.rda')


median_trajectories <- as.data.frame(matrix(nrow = n_cl, ncol = 12, NA))
par(mfrow = c(3,3))
for(i in 1:n_cl)
{
	todo <- exprs(eset.s)[row.names(exprs(eset.s)) %in% f_clusters[[i]],]
	median_trajectories[i,] <- colMedians(as.matrix(todo), na.rm = T)
	plot(y = median_trajectories[i,], x = c(0:11), type = 'l', lwd = 3, xlab = 'Time', ylab = paste0('Cluster ',i))
}

