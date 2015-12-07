setwd('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/')
source('C:/Users/am4613/Documents/GitHub/Proteomics/normalise_script.R')

avg_isx <- read.delim('isx_data_summary.txt')

pol1 <- read.delim('polymerase1', header = T, strings = F)
pol1 <- pol1$ensembl_id
pol2 <- read.delim('polymerase2', header = T, strings = F)
pol2 <- pol2$ensembl_id
pol3 <- read.delim('polymerase3', header = T, strings = F)
pol3 <- pol3$ensembl_id
chro <- read.delim('chromatin_remodelleres', header = T, strings = F)
chro <- chro$ensembl_id

common <- Reduce(intersect, list(pol1,pol2,pol3))
pol1and3 <- intersect(pol1,pol3)
pol1and3 <- setdiff(pol1and3, common)
	
prot_data <- read.delim('C:/Users/am4613/Documents/Summaries_as_timecourses/analysis/rep_cpc.txt', header = T, strings = F)

##Preprocessing, and averaging between technical replicates.Reordering

# norm_prot <- normalise_ProtDataset(prot_data, what = 'nada')
# norm_prot <- norm_prot[,7:42]

norm_prot <- 2^prot_data

fractions <- as.data.frame(matrix(nrow = 5, ncol = ncol(norm_prot), NA))
row.names(fractions) <- c('pol1','pol2','pol3','common', 'pol1&3')
colnames(fractions) <- colnames(norm_prot)

one <- colSums(norm_prot[which(row.names(norm_prot) %in% pol1),])
two <- colSums(norm_prot[which(row.names(norm_prot) %in% pol2), ])
three <- colSums(norm_prot[which(row.names(norm_prot) %in% pol3), ])
all <- colSums(norm_prot[which(row.names(norm_prot) %in% common), ])
oneAndThree <- colSums(norm_prot[which(row.names(norm_prot) %in% pol1and3), ])
remodellers <-  colSums(norm_prot[which(row.names(norm_prot) %in% chro), ])

total <- rbind(one, two, three, all, oneAndThree)
total <- colSums(total)

fractions[1,] <- one/total
fractions[2,] <- two/total
fractions[3,] <- three/total
fractions[4,] <- all/total
fractions[5,] <- oneAndThree/total

separate_frac <- list()

for(i in 1:3)
{
	separate_frac[[i]] <- fractions[,seq(i,33+i,3)]
}

par(mfrow = c(1,3))
for(i in 1:3)
{
	dfbar <- barplot(as.matrix(separate_frac[[i]]), col = c('lightblue','darkorange','red','white','purple'))
	par(new= T)
	length = avg_isx[avg_isx$rep == i, ]
	plot(x = dfbar, y = length$Length_Erode.M03..4., type = 'l', col = 'black', lwd = 3,xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	axis(4)
}

total_fractions <- fractions

all_prot_total <- colSums(norm_prot)

total_fractions[1,] <- one/all_prot_total
total_fractions[2,] <- two/all_prot_total
total_fractions[3,] <- three/all_prot_total
total_fractions[4,] <- all/all_prot_total
total_fractions[5,] <- oneAndThree/all_prot_total
total_fractions[6,] <- (all_prot_total - one - two - three - all - oneAndThree - remodellers)/all_prot_total
total_fractions[7,] <- remodellers/all_prot_total

separate_frac <- list()

for(i in 1:3)
{
	separate_frac[[i]] <- total_fractions[,seq(i,33+i,3)]
}

par(mfrow = c(1,3))
for(i in 1:3)
{
	dfbar <- barplot(as.matrix(separate_frac[[i]]), col = c('lightblue','darkorange','red','white','purple', 'grey','purple'))
	par(new= T)
	length = avg_isx[avg_isx$rep == i, ]
	plot(x = dfbar, y = length$Length_Erode.M03..4., type = 'l', col = 'black', lwd = 3,xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	axis(4)
}

par(mfrow = c(1,3))
for(i in 1:3)
{
	plot(x = c(0:11), y = separate_frac[[i]][7,], type = 'l')
}