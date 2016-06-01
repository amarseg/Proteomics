###Polynomial  regression to proteomics data
library('locfit')
library('rootSolve')
library('plyr')
library('cvTools')
library('polynom')

as_safequant <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/SQ_Results_PROTEIN.tsv', header = T, strings = F)
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

ids <- norm_sq[,1]
a <- strsplit(ids, '\\|')
ids <- sapply(a,"[[",1)
row.names(norm_sq) <- ids

##Fit polynomial models of degrees 1,2 or 3 to every gene and choosing the best model using K-fold crossvalidation
poly1 <- list()
poly2 <- list()
poly3 <- list()

best_model <- list()

time <- rep(0:11, each = 3)

for(i in 1:nrow(norm_sq))
{
	test <- cbind(t(norm_sq[i,7:42]),time)
	colnames(test) <- c('gene','time')
	test <- as.data.frame(test)
	poly1 <- lm(test$gene ~ poly(test$time,1,raw=TRUE))
	poly2 <- lm(test$gene ~ poly(test$time,2,raw=TRUE))
	poly3 <- lm(test$gene ~ poly(test$time,3,raw=TRUE))
	
	adj_r_squared <- c(summary(poly1)$adj.r.squared,summary(poly2)$adj.r.squared,summary(poly3)$adj.r.squared)
	pos <- which(adj_r_squared == max(adj_r_squared))
	if(pos == 1)
	{
		best_model[[i]] <- poly1
	}else if(pos == 2)
	{
		best_model[[i]] <- poly2
	}else if (pos == 3)
	{
		best_model[[i]] <- poly3
	}
	

}

model_degree <- c(1:nrow(norm_sq))
for(i in 1:nrow(norm_sq))
{
	model_degree[i] <-  length(best_model[[i]]$coefficients) -1
}

hist(model_degree)

first_derivative <- list()
second_derivative <- list()
sol_first <- list()
sol_second <- list()
for(i in 1:nrow(norm_sq))
{
	first_derivative[[i]] <- deriv(polynomial(best_model[[i]]$coefficients))
	sol_first[[i]] <- Re(polyroot(first_derivative[[i]]))
	second_derivative[[i]] <- deriv(polynomial(first_derivative[[i]]))
	sol_second[[i]] <- Re(polyroot(second_derivative[[i]]))
}

sol_firstdf <- ldply(sol_first, rbind)
sol_firstdf <- as.data.frame(sol_firstdf)
row.names(sol_firstdf) <- row.names(norm_sq)

hist(na.omit(rbind(sol_firstdf[,1], sol_firstdf[,2])), breaks = 100, freq = F)
lines(density(sol_firstdf[,1], na.rm = T, to = max(sol_firstdf[,1], na.rm = T)), col = 'red')
write.table(sol_firstdf, sep = '\t', 'C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/first_sol_prot.txt')

sol_seconddf <- ldply(sol_second, rbind)
sol_seconddf <- as.data.frame(sol_seconddf)
row.names(sol_seconddf) <- row.names(norm_sq)

hist(sol_seconddf[,1], breaks = 100, freq =F)
lines(density(sol_seconddf[,1], na.rm = T, to = max(sol_seconddf[,1], na.rm = T)), col = 'red')
write.table(sol_seconddf, sep = '\t', 'C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/second_sol_prot.txt')


