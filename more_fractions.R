setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/')
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')

avg_isx <- read.delim('isx_data_summary.txt')

aa_metabolic <- read.delim('aa_mtabolism', header = T, strings = F)
aa_metabolic_id <- aa_metabolic$ensembl_id
ribosomes <- read.delim('ribosomal_proteins', header = T, strings = F)
ribosomes_id <- ribosomes$ensembl_id
ribi <- read.delim('GO_ribosome_biogenesis', header = T, strings = F)
ribi_id <- ribi$ensembl_id
carb_metabolic <- read.delim('carb_metabolism', header = T, strings = F)
carb_id <- carb_metabolic$ensembl_id

ribi_id <- setdiff(ribi_id, ribosomes_id)

prot_data <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/rep_cpc.txt', header = T, strings = F)

##Preprocessing, and averaging between technical replicates.Reordering

# norm_prot <- normalise_ProtDataset(prot_data, what = 'nada')
# norm_prot <- norm_prot[,7:42]

norm_prot <- 2^prot_data

fractions <- as.data.frame(matrix(nrow = 5, ncol = ncol(norm_prot), NA))
row.names(fractions) <- c('ribosome','ribi','aa', 'carb', 'rest')
colnames(fractions) <- colnames(norm_prot)

ribi <- colSums(norm_prot[which(row.names(norm_prot) %in% ribi_id),])
ribo <- colSums(norm_prot[which(row.names(norm_prot) %in% ribosomes_id), ])
aa <- colSums(norm_prot[which(row.names(norm_prot) %in% aa_metabolic_id), ])
carb <- colSums(norm_prot[which(row.names(norm_prot) %in% carb_id), ])
total <- colSums(norm_prot)

fractions[1,] <- ribo/total
fractions[2,] <- ribi/total
fractions[3,] <- aa/total
fractions[4,] <- carb/total

fractions[5,] <- (total - ribo - ribi -aa-carb)/total

separate_frac <- list()

for(i in 1:3)
{
	separate_frac[[i]] <- fractions[,seq(i,33+i,3)]
}

par(mfrow = c(1,3))
for(i in 1:3)
{
	dfbar <- barplot(as.matrix(separate_frac[[i]]), col = c('darkgreen','lightgreen','darkorange','orange','lightgray'))
	par(new= T)
	length = avg_isx[avg_isx$rep == i, ]
	plot(x = dfbar, y = length$Length_Erode.M03..4., type = 'l', col = 'black', lwd = 3,xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	axis(4)
}

