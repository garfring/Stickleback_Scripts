## R script for Fst scan. It takes as input the .fst file of Popoolation2. Written by Antoine Paccard and modified by me (AGE). 
install.packages("data.table")
library(data.table)
getwd()
  
setwd("") # folder with Fst output from Popoolation2
list.files()
FST_all <- fread("output.fst", header = F, data.table = F) # Load FST data to use the first few columns
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
# using estuary names from which samples were taken
row.names(FST_pairs) <- c("Waddell","Lombardi",
                            "OldDairy","Younger",
                            "Scott","Laguna")

###OUTLIER DETECTION using all windows
FSToutliers_95 <- FST_pairs > apply(FST_all, 1, quantile, probs = 0.95, na.rm = T) # This identifies SNPs whose FST fall within the 5% estimate of being an outlier (SNPs in columns, pops in rows)
# FSToutliers_99 <- FST_pairs > apply(FST_all, 1, quantile, probs = 0.99, na.rm = T) # This identifies SNPs whose FST fall within the 1% estimate of being an outlier (SNPs in columns, pops in rows)


##chromosome info
colnames(FSToutliers_95) <- FST_all_info$V1
#colnames(FSToutliers_99) <- FST_all_info$V1

write.table(FSToutliers_95,"FST_outliers_95_v1.tsv", col.names = T, row.names = F, quote = F)
#write.table(FSToutliers_99,"FST_outliers_99.tsv", col.names = T, row.names = F, quote = F)
