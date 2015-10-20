##AS proteomics analysis
library('gplots')
library('gtools')
library('Mfuzz')
library('matrixStats')
load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')

as_safequant <- read.delim('C:/Users/am4613/Desktop/SQ_Results_PROTEIN.tsv', header = T, strings = F)


##Preprocessing, and averaging between technical replicates.Reordering

numbers <- as_safequant[,1:78]

titulos <- colnames(numbers)[7:78]
titulos <- substr(titulos, 0, nchar(titulos) - 3)
titulos <- unique(titulos)

norm_sq <- matrix(ncol = 6+36, nrow = nrow(as_safequant), NA)
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

##Normalisation to tp0

renorm_sq <- norm_sq

for(i in 1:3)
{
	renorm_sq[,seq(i+6,39+i,3)] <- norm_sq[,seq(i+6,39+i,3)]/norm_sq[,i+6]
}

#Log2 normalisation. NAs are changed to zeroes beacouse heatmap does not deal very well with NAs

renorm_sq[,7:42] <- log2(renorm_sq[,7:42])
is.na(renorm_sq[,7:42]) <- sapply(renorm_sq[,7:42],is.nan)
is.na(renorm_sq[,7:42]) <- sapply(renorm_sq[,7:42],is.infinite)
renorm_sq[,is.na(renorm_sq)] <- 0 

##Separate biological replicates into different lists, so they appear together on the heatmap

separate_prot <- list()

for(i in 1:3)
{
	separate_prot[[i]] <- renorm_sq[,seq(i+6,39+i,3)]
}

new_order <- cbind(separate_prot[[1]],separate_prot[[2]],separate_prot[[3]])

heatmap.2(as.matrix(new_order), trace = 'none', col = colorRampPalette(c('blue','gray','yellow')), Colv = F)

##Try out fuzzy clustering


renorm_sq[renorm_sq == 0] <- NA
row.names(renorm_sq) <- renorm_sq[,1]

titles <- colnames(norm_sq[,seq(7,42,3)])

##Averaging between biological replicates for fuzzy clustering

median_sq <- matrix(ncol = 6+12, nrow = nrow(as_safequant), NA)
median_sq <- as.data.frame(median_sq)
colnames(median_sq) <- c(colnames(norm_sq[,1:6]),titles)

j = 7

for(i in seq(7,42,3))
{
	median_sq[,j] <- rowMedians(cbind(renorm_sq[,i],renorm_sq[,i+1],renorm_sq[,i+2]), na.rm = T)
	
	j = j +1
}

##Create Expression set for fuzzy 
eset <- ExpressionSet(as.matrix(median_sq[,7:18]))

eset.r <- filter.NA(eset,0.5) #Threshold for discarding missing data 
eset.f <- fill.NA(eset.r, mode = 'mean') #Fill the gaps with the mean of the other values

tmp <- filter.std(eset.f, min.std = 0)

eset.s <- standardise(tmp)
n_cl = 9
m_fuzzy = 1.25
cl <- mfuzz(eset.s, c = n_cl, m = m_fuzzy)
mfuzz.plot(eset.s, cl = cl, mfrow = c(3,3), min.mem = 0.7, colo = rainbow(n=10))

##Code to extract the members of the different clusters and do GO analysis on them 

f_clusters <- list()

for(i in 1:n_cl)
{
	temp <- renorm_sq[as.numeric(names(cl$cluster[cl$cluster == i])),1]
	a <- strsplit(temp, "\\|")
	vector_with_ids <- sapply(a,"[[",1)
	f_clusters[[i]] <- vector_with_ids
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
write.table(output_final, 'C:/Users/am4613/Desktop/tableGO_fuzzy.txt', sep = '\t')

