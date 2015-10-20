##Analysis haploinsufficiency with new analogue sensitive datset 
library('gplots')


insu <- read.delim("C:/Users/am4613/Desktop/am4613/proteomics0804/haploinsuficient.txt", strings = F, header = F)
pro <- read.delim("C:/Users/am4613/Desktop/am4613/proteomics0804/haploproficient.txt", strings = F, header = F)
ids <- rbind(insu[1], pro[1])
category <-  c(rep("insufficcient", nrow(insu[1])), rep("proficient", nrow(pro[1])))
haplo_ids <- cbind(ids, category)

proteomics <- read.delim('C:/Users/am4613/Desktop/norm_as_proteomics.txt', header = T, strings = F)
transcriptomics <- read.delim('C:/Users/am4613/Desktop/me_rpkm.txt', header = T, strings = F)

transcriptomics[,1:12] <- transcriptomics[,1:12]/transcriptomics[,1]
transcriptomics[,13:24] <- transcriptomics[,13:24]/transcriptomics[,13]
transcriptomics[,25:36] <- transcriptomics[,25:36]/transcriptomics[,25]

transcriptomics <- log2(transcriptomics)
is.na(transcriptomics) <- sapply(transcriptomics, is.nan)
is.na(transcriptomics) <- sapply(transcriptomics, is.infinite)
transcriptomics <- na.omit(transcriptomics)


all_data <- merge(transcriptomics, proteomics, by.x = 'row.names', by.y = 'row.names', all.x = T)
all_data <- merge(haplo_ids, all_data, by.x = 'V1', by.y = 'Row.names', all.y = T)
all_data[,2] <- as.character(all_data[,2])

only_haplo <- all_data[!is.na(all_data[,2]),]

heatmap.2(as.matrix(only_haplo[,3:73]), col = colorRampPalette(c('blue','gray','yellow')),Colv = F, trace = 'none')


