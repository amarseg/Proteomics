###Combined heatmap transcriptomics and proteomics

library('gplots')

as_safequant <- read.delim('C:/Users/am4613/Desktop/SQ_Results_PROTEIN.tsv', header = T, strings = F)


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
ids <- renorm_sq[,1]
a <- strsplit(ids, '\\|')
ids <- sapply(a,"[[",1)
row.names(new_order) <- ids


load('P:/CDC2_RNA_051015.rda')

sam_rpkm <- list()

sam_rpkm[[1]] <- CDC2
sam_rpkm[[2]] <- CDC2_2
sam_rpkm[[3]] <- CDC2_3
sam_rpkm[[4]] <- CDC2_4

for(i in 2:4)
{
	sam_rpkm[[i]][,1:12] <- log2(sam_rpkm[[i]][,1:12]/sam_rpkm[[i]][,1])
	is.na(sam_rpkm[[i]]) <- sapply(sam_rpkm[[i]], is.nan)
	is.na(sam_rpkm[[i]]) <- sapply(sam_rpkm[[i]], is.infinite)
	sam_rpkm[[i]][is.na(sam_rpkm[[i]])] <- 0
	
}

norm_sam <- cbind(sam_rpkm[[2]][,1:12],sam_rpkm[[3]][,1:12],sam_rpkm[[4]][,1:12])

merg <- merge(norm_sam, new_order, by.x = 'row.names', by.y = 'row.names', all.y = T)

heatmap.2(as.matrix(merg[,2:73]), col = colorRampPalette(c('blue','gray','yellow')), trace = 'none', Colv = F)

##Identify rows with zeroes

