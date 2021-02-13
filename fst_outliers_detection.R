## R script for Fst scan. It takes as input the .fst file of Popoolation2.
install.packages("data.table")
library(data.table)
getwd()
  
setwd("C:/Users/User/Desktop/stickleback/santa_cruz/sc/fst_scan/")
list.files()
FST_all <- fread("fst_estuaries.fst", header = F, data.table = F) # Load FST data to use the first few columns
dim(FST_all)  
  
# Get info and metadata
FST_all_info <- FST_all[,1:5] # Keep the info columns
  
  
# Organise data frame
FST_all <- FST_all[,6:11] # Get rid of the info columns from the Fst file
names_all <- sapply(FST_all, function(i)unique(sub('=.*', '', i))) # extract names for columns
FST_all <- data.frame(lapply(FST_all, function(i)sub('.*=', '', i))) # separate values from "="
names(FST_all) <- names_all  # add column names
FST_all <- apply(as.matrix(FST_all), 2, as.numeric) # transform values as numeric

FST_all <- data.frame(t(FST_all)) # Transform into dataframe. Columns are SNPs and rows are pairwise comparisons
  
  
FST_pairs <- FST_all[,]
FST_pairs <- as.data.frame(FST_pairs)
row.names(FST_pairs) <- c("Waddell","Lombardi",
                            "OldDary","Younger"
                            ,"Scott","Laguna")

###OUTLIER DETECTION using all windows
FSToutliers_95 <- FST_pairs > apply(FST_all, 1, quantile, probs = 0.95, na.rm = T) # This identifies SNPs whose FST fall within the 5% estimate of being an outlier (SNPs in columns, pops in rows)
# FSToutliers_99 <- FST_pairs > apply(FST_all, 1, quantile, probs = 0.99, na.rm = T) # This identifies SNPs whose FST fall within the 5% estimate of being an outlier (SNPs in columns, pops in rows)


##chromosome info
colnames(FSToutliers_95) <- FST_all_info$V1
#colnames(FSToutliers_99) <- FST_all_info$V1

write.table(FSToutliers_95,"FST_outliers_95_v1.tsv", col.names = T, row.names = F, quote = F)
#write.table(FSToutliers_99,"FST_outliers_99.tsv", col.names = T, row.names = F, quote = F)


##########################################################
###get snp position for fst x3 and x4 outleirs
setwd("C:/Users/User/Desktop/stickleback/santa_cruz/sc/fst_scan")
list.files()
FST_all <- fread("fst_estuaries.fst", header = F, data.table = F) # Load FST data to use the first few columns
dim(FST_all)  

# Get info and metadata
FST_all_info <- FST_all[,1:5] # Keep the info columns


# Organise data frame
FST_all <- FST_all[,6:11] # Get rid of the info columns from the Fst file
names_all <- sapply(FST_all, function(i)unique(sub('=.*', '', i))) # extract names for columns
FST_all <- data.frame(lapply(FST_all, function(i)sub('.*=', '', i))) # separate values from "="
names(FST_all) <- names_all  # add column names
FST_all <- apply(as.matrix(FST_all), 2, as.numeric) # transform values as numeric

FST_pairs <- FST_all[,]
FST_pairs <- as.data.frame(FST_pairs)
col.names(FST_pairs) <- c("Waddell","Lombardi",
                          "OldDairy","Younger"
                          ,"Scott","Laguna")
col.names(FST_all) <- c("Waddell","Lombardi",
                          "OldDairy","Younger"
                          ,"Scott","Laguna")

par_3up <- read.table("FST_outliers_95_3up.tsv")
head(par_3up)
par_3up_index <- par_3up$V1

chromosome <- as.character(FST_all$V1[par_3up_index]) # chromosome position for outlier
snp_pos <- FST_all$V2[par_3up_index] # SNP position for outlier
num_est_outlier <- par_3up$V2 # number of estuaries (out of 6) that the snp is an outlier in
wad_fst <- 
lom_fst
oldd_fst
young_est
scott_est
lagun_est

fst_outlier_3up <- as.data.frame(cbind(chromosome, snp_pos, num_est_outlier))
head(fst_outlier_3up)
dim(fst_outlier_3up)
write.table(fst_outlier_3up,"fst_outlier_3up_snp_pos_top5.tsv", col.names = F, row.names = F,quote = F)


###############################################################
### get outliers for each estuary to get 'neutral allefreq change'
waddell.outlierindex <- which(FSToutliers_95[1,], arr.ind = TRUE)
lombardi.outlierindex <- which(FSToutliers_95[2,], arr.ind = TRUE)
olddairy.outlierindex <- which(FSToutliers_95[3,], arr.ind = TRUE)
younger.outlierindex <- which(FSToutliers_95[4,], arr.ind = TRUE)
scott.outlierindex <- which(FSToutliers_95[5,], arr.ind = TRUE)
laguna.outlierindex <- which(FSToutliers_95[6,], arr.ind = TRUE)


wadell_outliers <- as.numeric(as.matrix(FST_all[waddell.outlierindex]))
lombardi_outliers <- as.numeric(as.matrix(FST_all[lombardi.outlierindex]))
olddairy_outliers <- as.numeric(as.matrix(FST_all[olddairy.outlierindex]))
younger_outliers <- as.numeric(as.matrix(FST_all[younger.outlierindex]))
scott_outliers <- as.numeric(as.matrix(FST_all[scott.outlierindex]))
laguna_outliers <- as.numeric(as.matrix(FST_all[laguna.outlierindex]))


