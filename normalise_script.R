##Normalising proteomics dataset

cat('This function normalises the proteomics data in different ways
			depending on the what argument
			what = nada -> averaging between techinical replicates
			what = median -> median between biological replicates
			what = heatmap -> log2 fold change being tp0 the standard')

normalise_ProtDataset <- function(as_safequant, what)
{
	library('matrixStats')
	
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
	
	a <- strsplit(norm_sq[,1], "\\|")
	vector_with_ids <- sapply(a,"[[",1)
	row.names(norm_sq) <- vector_with_ids
	
	if(what == 'median')
	{
		median_sq <- matrix(ncol = 6+12, nrow = nrow(numbers), NA)
		median_sq <- as.data.frame(median_sq)
		colnames(median_sq) <- c(colnames(norm_sq[,1:6]),unique(sapply(strsplit(titulos, '_'),"[[",2)))
		row.names(median_sq) <- row.names(norm_sq)
		
		j = 7
		
		for(i in seq(7,40,3))
		{
			median_sq[,j] <- rowMedians(cbind(norm_sq[,i],norm_sq[,i+1],norm_sq[,i+2]), na.rm = T)
			j = j +1
		}
		
		return(median_sq)
		
	}else if(what == 'heatmap')
	{
		renorm_sq <- norm_sq
		
		for(i in 1:3)
		{
			renorm_sq[,seq(i+6,39+i,3)] <- norm_sq[,seq(i+6,39+i,3)]/norm_sq[,i+6]
		}
		
		#Log2 normalisation. NAs are changed to zeroes beacouse heatmap does not deal very well with NAs
		
		renorm_sq[,7:42] <- log2(renorm_sq[,7:42])
		is.na(renorm_sq[,7:42]) <- sapply(renorm_sq[,7:42],is.nan)
		is.na(renorm_sq[,7:42]) <- sapply(renorm_sq[,7:42],is.infinite) 
		return(renorm_sq)
	}else{
		return(norm_sq)
	}
}

reorder_proteomics <- function(dataset)
{
	reorder <- list()
	
	for(i in seq(1,3))
	{	
		reorder[[i]] <- dataset[,seq(i,33+i,3)]	
	}
	
	back <- cbind(reorder[[1]], reorder[[2]], reorder[[3]])
	
	return(back)
}

normalise_rna <- function(table)
{
	for(i in c(1,13,25))
	{
		table[,i:(i+11)] <- log2(table[,i:(i+11)]/table[,i])
	}
	
	is.na(table) <- sapply(table, is.nan)
	is.na(table) <- sapply(table, is.infinite)
	table[is.na(table)] <- 0 
	return(table)
}

df2list <- function(df)
{
	out <- list()
	for(i in seq(1,ncol(df)))
	{
		out[[i]] <- df[,i]
	}
	return(out)
}

ggplotRegression <- function (fit) {
	
	require(ggplot2)
	
	ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
		geom_point() +
		stat_smooth(method = "lm", col = "red") +
		labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
											 "Intercept =",signif(fit$coef[[1]],5 ),
											 " Slope =",signif(fit$coef[[2]], 5),
											 " P =",signif(summary(fit)$coef[2,4], 5)))
}

median_rna <- function(df)
{
	avg_rpkm <- df[,1:12]
	for(i in 1:12)
	{
		avg_rpkm[,i] <- apply(df[,seq(i,i+33,12)], 1, mean)
	}
	return(avg_rpkm)
	
}

nc_remover <- function(ids)
{
	new_ids <- ids[grep(ids, pattern = 'SPNC', invert = T)]
	return(new_ids)
}