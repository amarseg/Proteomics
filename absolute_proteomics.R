rm(list = ls())


library("dplyr")
library('matrixStats')
library('DAAG')
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')
source('C:/Users/am4613/Documents/GitHub/Misc/loov.r')
library('gplots')



sam_cpc <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/sam_cpc.txt', header = T, strings = F)

pep_amount <- 150 #femptomol per 1 ugr of protein extract 
micrograms_per_injection <- 3

cell_number <- read.delim("C:/Users/am4613/Documents/Summaries_as_timecourses/protein_per_cell_as.txt")
colnames(cell_number) <- c('Time','Replicate','Protein.per.cell')
cell_number$Cell.per.injection <- micrograms_per_injection/cell_number$Protein.per.cell
rep_cell_number <- cell_number

avg_cell_number <- aggregate(Cell.per.injection ~ Time, data = cell_number, median)

######
# Preparation PRM results
######

##Load PRM results

prm_results <- read.csv("C:/Users/am4613/Desktop/Peptide Ratio Results_cdc2as.csv", header = T, strings = F)
prm_results[prm_results == "#N/A"] <- NA
prm_results[,4:6] <- apply(prm_results[,4:6], c(1,2), as.numeric)

standard_key <- read.delim("C:/Users/am4613/Documents/Skyline/Standards_key.txt", header = T, strings = F)
standard_key <- standard_key[,1:2]

temp<-t(matrix(unlist(strsplit(as.character(prm_results$Replicate.Name), "_",fixed=TRUE)),nrow=3))
prm_results$timepoint <- temp[,1]

#Calculate median ratios between technical replicates
standard_quant <- as.data.frame(prm_results %>% group_by(timepoint, Peptide.Sequence) %>% summarize(RatioToStandard.median = median(Ratio.To.Standard, na.rm = T)))
standard_quant$ID <- standard_key[match(standard_quant[,2], standard_key[,1]),2]

standard_quant$AmountLight <- standard_quant$RatioToStandard * pep_amount
standard_quant$ProteinNumber <- (standard_quant$AmountLight*6.022*10^23)/10^15

standard_quant <- standard_quant[standard_quant$ID != 'SPCC191.02c',]
standard_quant <- standard_quant[standard_quant$ID != 'SPBC28F2.12',]
standard_quant <- standard_quant[standard_quant$ID != 'SPBP22H7.08',]



wo_list <- split(standard_quant, standard_quant[,1])
###########Preparation proteomics results (coming from SafeQuant)
###########

##IN SILICO DIGESTION OF THE PROTEOME NEEDS TO BE REDONE, so far i am sticking to the one that Erik made for me in Basel

digest <- read.delim("C:/Users/am4613/Documents/Amalia_Start_201113/Absolute/pompep_290811_plus_Pseudo.decoy.DIGEST.tsv",
										 strings = F, header = T)

#Remove all the REV hits (reverse sequences from the decoy database)
digest <- digest[grep("REV", digest[,1], value = F, invert = T),]
n_tryptic <- data.frame(digest$ac, digest$nbTrypticPeptides)


as_safequant <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/SQ_Results_PROTEIN.tsv', header = T, strings = F)


##Preprocessing, and averaging between technical replicates.Reordering

norm_sq <- normalise_ProtDataset(as_safequant, what = 'nada')
norm_sq <- norm_sq[,7:42]

#Match digest with protein_data

proteomics_tryptic <- merge(norm_sq, n_tryptic, by.x = "row.names", by.y = "digest.ac", all.x = T)

#Normalising to number of tryptic peptide

norm_prot <- proteomics_tryptic[,2:37]/proteomics_tryptic$digest.nbTrypticPeptides
row.names(norm_prot) <- proteomics_tryptic[,1]


##Combine norm_prot w/ standard quant for lineal regression

for(i in seq(1,36))
{
	toDo <- wo_list[[i]]
	x <- merge(toDo, norm_prot, by.x = "ID", by.y = "row.names", all.x = T)
	wo_list[[i]] <- cbind(x[,1:6], x[,i+6])
	colnames(wo_list[[i]]) <- c(colnames(wo_list[[i]])[1:6],'ExpIntensity')
	
}


##Linear regression for all lists

parameters <- data.frame(timepoint = rep(0:11, each = 3), intercept = NA, slope = NA)
corrs <- c(1:36)
for(i in seq(1,36))
{
	thing <- wo_list[[i]]$ProteinNumber/rep_cell_number$Cell.per.injection[i]
	fit <- lm(log2(wo_list[[i]]$ExpIntensity) ~ log2(thing) + 0)
	corrs[i] <- summary(fit)$r.squared
	parameters[i,2] <- fit$coefficients[[1]]
	#parameters[i,3] <- fit$coefficients[[2]]
}

###Plotting biolog2ical replicates together to see variance

par(mfrow = c(3,4))
for(i in seq(1,34, 3))
{
	plot(x=log2(wo_list[[i]]$ProteinNumber/rep_cell_number[i]), y=log2(wo_list[[i]]$ExpIntensity), pch = 2)
	points(x=log2(wo_list[[i+1]]$ProteinNumber/rep_cell_number[i]), y=log2(wo_list[[i+1]]$ExpIntensity),   pch = 2, col = 'red')
	points(x=log2(wo_list[[i+2]]$ProteinNumbe/rep_cell_number[i]), y=log2(wo_list[[i+2]]$ExpIntensity),   pch = 2, col = 'blue')
	
	abline(0, parameters[i,2])
	abline(0, parameters[i+1,2], col = 'red')
	abline(0, parameters[i+2,2], col = 'blue')
}

##Calculate copies per cell for the different replicates

rep_cpc <- norm_sq
for(i in 1:36)
{
	rep_cpc[,i] <- log2(norm_sq[,i])/parameters$intercept[i]
}

par(mfrow = c(1,3))
boxplot(reorder_proteomics(rep_cpc)[,1:12], outline = F)
boxplot(reorder_proteomics(rep_cpc)[,13:24], outline = F)
boxplot(reorder_proteomics(rep_cpc)[,25:36], outline = F)

write.table(rep_cpc, sep = '\t','C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/rep_cpc.txt')

##Averaging biolog2ical repeats

avg_st_quant <- list()

j = 1
for(i in seq(1,34,3))
{
	avg_st_quant[[j]] <- wo_list[[i]]
	toDo <- cbind(wo_list[[i]]$ProteinNumber, wo_list[[i+1]]$ProteinNumber, wo_list[[i+2]]$ProteinNumber)
	avg_st_quant[[j]]$ProteinNumber <- rowMedians(toDo, na.rm = T)
	toDo <- cbind(wo_list[[i]]$ExpIntensity, wo_list[[i+1]]$ExpIntensity, wo_list[[i+2]]$ExpIntensity)
	avg_st_quant[[j]]$ExpIntensity <- rowMedians(toDo, na.rm = T)
	j = j+1
}

##Calculate median of the whole proteomics dataset
median_norm <- normalise_ProtDataset(as_safequant, what = 'median')
median_norm <- median_norm[,7:18]

##Calculate linear regression for all the calibration curves, also plots them
par(mfrow = c(3,4))
bio_corrs <- c(1:12)
bio_parameters <- data.frame(timepoint = 0:11, intercept = NA, slope = NA)
r_squared <- data.frame(matrix(nrow = nrow(avg_st_quant[[1]]), ncol = 12, NA))
for(j in c(1:12))
{
	rm(fit)
	avg_st_quant[[j]]$ProteinPerCell <- avg_st_quant[[j]]$ProteinNumber/avg_cell_number$Cell.per.injection[j]
	fit <- lm(log2(avg_st_quant[[j]]$ExpIntensity) ~ log2(avg_st_quant[[j]]$ProteinPerCell) + 0)
	bio_corrs[j] <- summary(fit)$r.squared
	bio_parameters[j,2] <- fit$coefficients[[1]]
	#bio_parameters[j,3] <- fit$coefficients[[2]]
	r_squared[,j] <- loov(y = avg_st_quant[[j]]$ExpIntensity, x = avg_st_quant[[j]]$ProteinPerCell, log = T)
	#confidence95 <- confint(fit, level = 0.95)
	plot(x=log2(avg_st_quant[[j]]$ProteinPerCell), y=log2(avg_st_quant[[j]]$ExpIntensity),  col = rainbow(nrow(avg_st_quant[[j]])), pch = 2)
	abline(fit)
	#abline(confidence95[,1], col = 'red')
	#abline(confidence95[,2], col = 'red')

	legend('topleft', legend = paste('R squared: ',round(summary(fit)$r.squared,3)), bty = 'n')
}

##Plot to see which peptides affect most the r squared
par(mfrow = c(3,4))
for(i in 1:12)
{
	plot(r_squared[,i], type = 'h')
}

##Calculating copies per cell 
cpc <- median_norm
par(mfrow = c(1,1))
for(i in 1:12)
{
	cpc[,i] <- log2(median_norm[,i])/bio_parameters$intercept[i] 
}

write.table(cpc, sep = '\t','C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/cpc_proteins.txt')


f <- colorRampPalette(c('blue','red'))
col = f(12)

plot(density(cpc[,1]), col = col[1])
for(i in 2:12)
{
	lines(density(cpc[,i]), col = col [i])
}

plot(x = colMedians(as.matrix(cpc)), y = bio_parameters$slope, col = f(12))

norm_cpc <- cpc - cpc[,1]
real_cpc <- 2^cpc
heatmap.2(as.matrix(cpc), Colv = F, col = colorRampPalette(c('lightblue','darkblue')), trace = 'none')


