###Combined heatmap transcriptomics and proteomics

library('gplots')

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


##Separate biological replicates into different lists, so they appear together on the heatmap

separate_prot <- list()

for(i in 1:3)
{
	separate_prot[[i]] <- norm_sq[,seq(i+6,39+i,3)]
}

new_order <- cbind(separate_prot[[1]],separate_prot[[2]],separate_prot[[3]])
ids <- norm_sq[,1]
a <- strsplit(ids, '\\|')
ids <- sapply(a,"[[",1)
row.names(new_order) <- ids

load('P:/CDC2_RNA_051015.rda')

sam_rpkm <- list()

sam_rpkm[[1]] <- CDC2
sam_rpkm[[2]] <- CDC2_2
sam_rpkm[[3]] <- CDC2_3
sam_rpkm[[4]] <- CDC2_4

norm_sam <- cbind(sam_rpkm[[2]][,1:12],sam_rpkm[[3]][,1:12],sam_rpkm[[4]][,1:12])

merg <- merge(new_order, norm_sam, by.x = 'row.names', by.y = 'row.names', all.x = T)

cor_spearman <- as.data.frame(matrix(ncol = 36, nrow = 1, NA))
cor_pearson <- as.data.frame(matrix(ncol = 36, nrow = 1, NA))

for(i in 1:36)
{
	cor_spearman[,i] <- cor(merg[,i+1], merg[i+37], method = 'spearman')
	cor_pearson[,i] <- cor(merg[,i+1], merg[i+37], method = 'pearson')
}

time = rep(0:11,each = 3)

all <- rbind(cor_spearman, cor_pearson, time)
all = t(all)
all[,3] <- factor(all[,3])
colnames(all) <- c('spearman','pearson', 'time')

par(mfrow = c(1,2))
lmts <- range(all[,1], all[,2])
boxplot(pearson ~ time, data = all, names = c(0:11), ylim = lmts)
boxplot(spearman ~ time, data = all, names = c(0:11), ylim = lmts)
