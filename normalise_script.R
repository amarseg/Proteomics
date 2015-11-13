##Normalising proteomics dataset

normalise_ProtDataset <- function(as_safequant, what)
{
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
		colnames(median_sq) <- c(colnames(norm_sq[,1:6]),titles)
		
		j = 7
		
		for(i in seq(7,42,3))
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
		renorm_sq[,is.na(renorm_sq)] <- 0 
		return(renorm_sq)
	}else{
		return(norm_sq)
	}
}