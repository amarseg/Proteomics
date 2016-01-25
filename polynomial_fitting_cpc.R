###Polynomial  regression to proteomics data
rm(list = ls())
library('locfit')
library('rootSolve')
library('plyr')
library('cvTools')
library('polynom')
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')

use_replicates_data = T

if(use_replicates_data)
{
	rep_cpc <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/rep_cpc.txt', header = T, strings = F)
	rep_cpc <- 2^rep_cpc
	rep_cpc <- reorder_proteomics(rep_cpc)
	time <- rep(0:11,each = 3)
}else{
	rep_cpc <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/cpc_proteins.txt', header = T, strings = F)
	rep_cpc <- 2^rep_cpc
	rep_cpc <- reorder_proteomics(rep_cpc)
	time <- rep(0:11)
}

##Fit polynomial models of degrees 1,2 or 3 to every gene and choosing the best model using K-fold crossvalidation
poly1 <- list()
poly2 <- list()
poly3 <- list()

best_model <- list()



for(i in 1:nrow(rep_cpc))
{
	test <- cbind(rep_cpc[i,],time)
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

model_degree <- c(1:nrow(rep_cpc))
for(i in 1:nrow(rep_cpc))
{
	model_degree[i] <-  length(best_model[[i]]$coefficients) -1
}

hist(model_degree)

first_derivative <- list()
second_derivative <- list()
sol_first <- list()
sol_second <- list()
for(i in 1:nrow(rep_cpc))
{
	first_derivative[[i]] <- deriv(polynomial(best_model[[i]]$coefficients))
	sol_first[[i]] <- Re(polyroot(first_derivative[[i]]))
	second_derivative[[i]] <- deriv(polynomial(first_derivative[[i]]))
	sol_second[[i]] <- Re(polyroot(second_derivative[[i]]))
}

sol_firstdf <- ldply(sol_first, rbind)
sol_firstdf <- as.data.frame(sol_firstdf)
row.names(sol_firstdf) <- row.names(rep_cpc)

hist(na.omit(c(sol_firstdf[,1],sol_firstdf[,2])), breaks = 100, freq = F)
lines(density(c(sol_firstdf[,1],sol_firstdf[,2]), na.rm = T), col = 'red')
#write.table(sol_firstdf, sep = '\t', 'C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/first_sol_prot.txt')

sol_seconddf <- ldply(sol_second, rbind)
sol_seconddf <- as.data.frame(sol_seconddf)
row.names(sol_seconddf) <- row.names(rep_cpc)

hist(sol_seconddf[,1], breaks = 100, freq =F)
lines(density(sol_seconddf[,1], na.rm = T, to = max(sol_seconddf[,1], na.rm = T)), col = 'red')
#write.table(sol_seconddf, sep = '\t', 'C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/second_sol_prot.txt')

