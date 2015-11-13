library("dplyr")
library('matrixStats')

pep_amount <- 150 #femptomol per 1 ugr of protein extract 
micrograms_per_injection <- 3
cell_number <- read.delim("C:/Users/am4613/Documents/Summaries_as_timecourses/protein_per_cell_sho1.txt")

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

st_quant_list <- split(standard_quant, standard_quant[,1])
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
#Trim the IDs to get sistematic IDs

a <- strsplit(norm_sq[,1], "\\|")
row.names(norm_sq) <- sapply(a,"[[",1)
norm_sq <- norm_sq[,7:42]

#Match digest with protein_data

proteomics_tryptic <- merge(norm_sq, n_tryptic, by.x = "row.names", by.y = "digest.ac", all.x = T)

#Normalising to number of tryptic peptide

norm_prot <- proteomics_tryptic[,2:37]/proteomics_tryptic$digest.nbTrypticPeptides
row.names(norm_prot) <- proteomics_tryptic[,1]


##Combine norm_prot w/ standard quant for lineal regression

for(i in seq(1,36))
{
	toDo <- st_quant_list[[i]]
	x <- merge(toDo, norm_prot, by.x = "ID", by.y = "row.names", all.x = T)
	st_quant_list[[i]]$ExpIntensity <- x[,i+6]
	
}


##Linear regression for all lists

parameters <- data.frame(timepoint = rep(0:11, each = 3), intercept = NA, slope = NA)
corrs <- c(1:36)
for(i in seq(1,36))
{
	fit <- lm(st_quant_list[[i]]$ProteinNumber ~ st_quant_list[[i]]$ExpIntensity)
	corrs[i] <- cor(x = st_quant_list[[i]]$ProteinNumber,y = st_quant_list[[i]]$ExpIntensity, use = 'complete.obs')	
	parameters[i,2] <- fit$coefficients[[1]]
	parameters[i,3] <- fit$coefficients[[2]]
}

cpc <- norm_prot

cell_number$Cell.per.injection <- micrograms_per_injection/cell_number$Protein.per.cell
rep_cell_number <- rep(cell_number$Cell.per.injection, each = 3)

##Calibrate normalised proteomics data

for(i in seq(1,ncol(norm_prot)))
{
	cpc[,i] <- norm_prot[,i]*parameters[i,3]
	cpc[,i] <- cpc[,i]/rep_cell_number[i]
}

avg_st_quant <- list()

j = 0
for(i in seq(1,34,3))
{
	avg_st_quant[[j]] <- st_quant_list[[i]]
	avg_st_quant[[j]]$ProteinNumber <- median(st_quant_list[[i]]$ProteinNumber, st_quant_list[[i+1]]$ProteinNumber, st_quant_list[[i+2]]$ProteinNumber)
	avg_st_quant[[j]]$ExpIntensity <- median(st_quant_list[[i]]$ExpIntensity, st_quant_list[[i+1]]$ExpIntensity, st_quant_list[[i+2]]$ExpIntensity)
	j = j+1
}
