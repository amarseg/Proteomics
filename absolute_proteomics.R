library("dplyr")

######
# Preparation PRM results
######

##Load PRM results

prm_results <- read.csv("C:/Users/am4613/Desktop/Skyline/PRM_pombe.csv", header = T, strings = F)
prm_results[prm_results == "#N/A"] <- NA



#########
#Preparation proteomics results (coming from SafeQuant)
###########

##IN SILICO DIGESTION OF THE PROTEOME NEEDS TO BE REDONE

digest <- read.delim("C:/Users/am4613/Documents/Amalia_Start_201113/Absolute/pompep_290811_plus_Pseudo.decoy.DIGEST.tsv",
										 strings = F, header = T)

#Remove all the REV hits (reverse sequences from the decoy database)

digest <- digest[grep("REV", digest[,1], value = F, invert = T),]

protein_data <- read.csv("C:/Users/am4613/Desktop/Top2rank Protein measurements_onlyNorm.csv", header = T, strings = F)
avg_proteindata <- matrix(nrow = nrow(protein_data), ncol = 13,0)
avg_proteindata <- as.data.frame(avg_proteindata)

##Averaging between techinical replicates (still need to take care of cases where one is 0)
j = 1
for(i in seq(10,33,2))
{
	avg_proteindata[,j+1] <- rowMeans(cbind(protein_data[,i], protein_data[,i+1]), na.rm = T)
	j= j + 1
}

avg_proteindata[,1] <- protein_data[,1]

#Trim the IDs to get sistematic IDs

a <- strsplit(avg_proteindata[,1], "\\|")
avg_proteindata[,1] <- sapply(a,"[[",1)

#Match digest with protein_data

proteomics_tryptic <- merge(avg_proteindata, digest, by.x = 1, by.y = 1, all.x = T)

#Normalising to number of tryptic peptide

proteomics_tryptic[,2:13] <- proteomics_tryptic[,2:13]/proteomics_tryptic$nbTrypticPeptides




