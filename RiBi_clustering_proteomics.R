##

library('gplots')
library('gtools')
library('Mfuzz')
library('matrixStats')

biogenesis <- read.delim('GO_ribosome_biogenesis', header = T, strings = F)
rps <- read.delim('rproteins.txt', header = T, strings = F)

load('C:/Users/am4613/Documents/GitHub/Misc/GO.analysis.110914.rda')

as_safequant <- read.delim('C:/Users/am4613/Desktop/SQ_Results_PROTEIN.tsv', header = T, strings = F)


##Preprocessing, and averaging between technical replicates.Reordering

numbers <- as_safequant[,1:78]
numbers <- numbers[grep(numbers[,1], pattern = 'SP'),]

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

median_sq <- median_sq[,7:18]
ids <- renorm_sq[,1]
a <- strsplit(ids, '\\|')
ids <- sapply(a,"[[",1)
row.names(median_sq) <- ids

prot_bio <- subset(median_sq,row.names(median_sq) %in% biogenesis[,1])
prot_rp <- subset(median_sq,row.names(median_sq) %in% rps[,3])

bioset <- ExpressionSet(as.matrix(prot_bio))
rpset <- ExpressionSet(as.matrix(prot_rp))

bioset.r <- filter.NA(bioset,0.5) #Threshold for discarding missing data 
bioset.f <- fill.NA(bioset.r, mode = 'mean') #Fill the gaps with the mean of the other values
rpset.r <- filter.NA(rpset,0.5) #Threshold for discarding missing data 
rpset.f <- fill.NA(rpset.r, mode = 'mean') #Fill the gaps with the mean of the other values

bio_tmp <- filter.std(bioset.f, min.std = 0)
rp_tmp <- filter.std(rpset.f, min.std = 0)

bioset.s <- standardise(bio_tmp)
rpset.s <- standardise(rp_tmp)
n_cl = 5
m_fuzzy = 1.25
biocl <- mfuzz(bioset.s, c = n_cl, m = m_fuzzy)
mfuzz.plot(bioset.s, cl = biocl, mfrow = c(3,3), min.mem = 0.7, colo = rainbow(n=10))

rpcl <- mfuzz(rpset.s, c = n_cl, m = m_fuzzy)
mfuzz.plot(rpset.s, cl = rpcl, mfrow = c(3,3), min.mem = 0.7, colo = rainbow(n=10))

##Code to extract the members of the different clusters and do GO analysis on them 

bio_clusters <- list()

for(i in 1:n_cl)
{
	bio_clusters[[i]] <- names(biocl$cluster[biocl$cluster == i])
}


rp_clusters <- list()

for(i in 1:n_cl)
{
	rp_clusters[[i]] <- names(rpcl$cluster[rpcl$cluster == i])
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

