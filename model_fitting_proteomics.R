###Fit LOESS regression to proteomics data
library('locfit')
library('rootSolve')
library('plyr')
library('cvTools')

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

time <- rep(0:11,each = 3)

result_loess <- list()
result_spline <- list()
result_locfit <- list()
result_locfit2 <- list()
result_locfit1 <- list()
rev <- list()

for(i in 1:nrow(norm_sq))
{
	test <- cbind(t(norm_sq[i,7:42]),time)
	colnames(test) <- c('gene','time')
	test <- as.data.frame(test)
	result_loess[[i]] <- loess(gene~time, data = test)
	result_locfit[[i]] <- locfit(gene~time, data = test)
	result_locfit2[[i]] <- locfit(gene~time, data = test, deriv = c(1,1)) #That's how you express the second derivative
	result_locfit1[[i]] <- locfit(gene~time, data = test, deriv = 1)
	result_spline[[i]] <- smooth.spline(x = test$time, y= test$gene, all.knots = T)
}

###Code to plot the fit of any gene
plot_loess_fit <- function(gene_name)
{
	j = which(row.names(norm_sq) == gene_name)
	plot(y = norm_sq[j,7:42], x = time)
	lines(time,result_loess[[j]]$fitted, col = 'red', lwd = 3)
	lines(result_spline[[j]], col = 'blue', lwd = 3)
	derivada <- predict.smooth.spline(result[[j]])
	lines(derivada$y, col = 'green', lwd = 3)
	lines(result_locfit[[j]], col = 'purple', lwd = 3)
	plot(result_locfit2[[j]])
}

roots_locfit2 = list()
for(i in 1:length(result_locfit2))
{
	roots_locfit2[[i]] <- uniroot.all(function(t) predict(result_locfit2[[i]],t), interval = c(0,11))
}

roots_locfit2 <- ldply(roots_locfit2, rbind)

meh <- rbind(roots_locfit2[,1],roots_locfit2[,2],roots_locfit2[,3])
hist(meh, breaks = 100, main = 'Histogram of roots', freq = F)
lines(density(meh, na.rm = T, from = 0, to = max(meh, na.rm = T)), col = 'red')

roots_locfit1 = list()
for(i in 1:length(result_locfit1))
{
	roots_locfit1[[i]] <- uniroot.all(function(t) predict(result_locfit1[[i]],t), interval = c(0,11))
}

roots_locfit1 <- ldply(roots_locfit1, rbind)

meh <- rbind(roots_locfit1[,1],roots_locfit1[,2],roots_locfit1[,3], roots_locfit1[,4])
hist(meh, breaks = 100, main = 'Histogram of roots', freq = F)
lines(density(meh, na.rm = T, from = 0, to = max(meh, na.rm = T)), col = 'red')

myInterval <- rowSums(roots_locfit2 > 2 & roots_locfit2 < 5 , na.rm = T)
genesInterval <- row.names(norm_sq[myInterval > 0,])
write.table(genesInterval, sep = '\t', 'P_2&5_2ndDerivative.txt')


