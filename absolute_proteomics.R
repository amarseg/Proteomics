library("dplyr")

pep_amount <- 50 #femptomol per 1 ugr of protein extract 
cell_number <- read.delim("C:/Users/am4613/Desktop/Skyline/cell_numbers_ts.txt")
cell_number <- as.data.frame(cbind(cell_number[,1], cell_number[,24]))

######
# Preparation PRM results
######

##Load PRM results

prm_results <- read.csv("C:/Users/am4613/Desktop/Skyline/Alex_files/Peptide Ratio Results.csv", header = T, strings = F)
prm_results[prm_results == "#N/A"] <- NA
prm_results[,4:6] <- apply(prm_results[,4:6], c(1,2), as.numeric)

standard_key <- read.delim("C:/Users/am4613/Desktop/Skyline/Standards_key.txt", header = T, strings = F)
standard_key <- standard_key[,1:2]

temp<-t(matrix(unlist(strsplit(as.character(prm_results$Replicate.Name), "_",fixed=TRUE)),nrow=3))
prm_results$timepoint <- temp[,1]

#Calculate median ratios between technical replicates
standard_quant <- as.data.frame(prm_results %>% group_by(timepoint, Peptide.Sequence) %>% summarize(RatioToStandard.median = median(RatioToStandard, na.rm = T)))
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


protein_data <- read.csv("C:/Users/am4613/Desktop/Top2rank Protein measurements_onlyNorm.csv", header = T, strings = F)
protein_data[protein_data == 0] <- NA 

avg_proteindata <- matrix(nrow = nrow(protein_data), ncol = 12,0)
avg_proteindata <- as.data.frame(avg_proteindata)

##Averaging between techinical replicates (still need to take care of cases where one is 0)
j = 1
for(i in seq(10,33,2))
{
	avg_proteindata[,j] <- rowMeans(cbind(protein_data[,i], protein_data[,i+1]), na.rm = T)
	j= j + 1
}

#Trim the IDs to get sistematic IDs

a <- strsplit(protein_data[,1], "\\|")
row.names(avg_proteindata) <- sapply(a,"[[",1)

#Get sample names from raw data
sample_names <- colnames(protein_data[,10:33])
sample_names <- substr(sample_names, start = 0, stop = nchar(sample_names) - 3 )
colnames(avg_proteindata) <- unique(sample_names)

#Match digest with protein_data

proteomics_tryptic <- merge(avg_proteindata, n_tryptic, by.x = "row.names", by.y = "digest.ac", all.x = T)

#Normalising to number of tryptic peptide

norm_prot <- proteomics_tryptic[,2:13]/proteomics_tryptic[,14]
row.names(norm_prot) <- proteomics_tryptic[,1]
colnames(norm_prot) <- as.character(seq(0,11))

##Combine norm_prot w/ standard quant for lineal regression

for(i in seq(1,12))
{
	toDo <- st_quant_list[[i]]
	x <- merge(toDo, norm_prot, by.x = "ID", by.y = "row.names", all.x = T)
	st_quant_list[[i]]$ExpIntensity <- x[,6 + i]
	
}


##Linear regression for all lists

parameters <- data.frame(timepoint = c(0:11), intercept = NA, slope = NA)

for(i in seq(1:12))
{
	fit <- lm(st_quant_list[[i]]$ProteinNumber ~ st_quant_list[[i]]$ExpIntensity)
	parameters[i,2] <- fit$coefficients[[1]]
	parameters[i,3] <- fit$coefficients[[2]]
}



##Calibrate normalised proteomics data

for(i in seq(1:12))
{
	norm_prot[,i] <- norm_prot[,i]*parameters[i,2] + parameters[i,1]
	norm_prot[,i] <- norm_prot[,i]/cell_number[i,2]
}


