###Code to generate IPASS and assign T ell infiltration status
#T-cell infiltration indicates a tumour based off the IPASS as either
#immune-inflamed or altered
#Reads in a gene expression data frame of TPM values
#code can be run on any number of samples

#Modified on May 12, 2023 to fit in the RICO container

#load("IPASS.RData")
library(tidyr)
library(dplyr)

# gather the arguments
args = commandArgs(trailingOnl = TRUE)

# test if there are 3 arguments: if not, return an error
if (length(args) != 3) { 
  stop("Three arguments are required", call. = FALSE)
}

TPM_file = args[1]
MATRIX_file = args[2]
IPASS_data = args[3]

##load IPASS related information
load(IPASS_data)

#load TPM expression data - columns will be patients and rows hugo gene symbols
TPMs <- read.delim(TPM_file, sep = " ",
            header = TRUE,
            stringsAsFactors = FALSE)

#load the comparison matrix data
M = read.csv(MATRIX_file, sep = " ", header = TRUE)

#add the sample to the big matrix
TPMs = TPMs %>% 
    inner_join(M)

## extract out ipass genes
# note gene_id is column name of the gene names
#check to ensure there are 15 rows in the ipassTPM dataframe
ipassTpm <- TPMs[which(TPMs$HGNC %in% ipassGenes$IPASS),]

##turn into a matrix and add 1 to handle for 0 TPM values
rownames(ipassTpm)<-ipassTpm$gene_id
ipassTpm <- ipassTpm[, -1]
ipassTpm <- as.matrix(ipassTpm)
ipassTpm <- ipassTpm + 1

#calculate score for each sample
ipassScore<-log(ipassTpm,base=10)
Score<-colMeans(ipassScore)
ipassScore<-rbind(ipassScore,Score)

#convert into data frame with IPASS value and infiltration status
#added in a column for the sample's percentile against the MATRIX cohort
tipassScore<-as.data.frame(t(ipassScore))
tipassScore$Sample<-rownames(tipassScore)
tipassScore<-tipassScore[,c("Sample","Score")]
tipassScore$Status<-ifelse(tipassScore$Score >= threshold,"T-cell_Infiltrated","Cold")

#Score percentile against the matrix
#Assumes TPM is the top sample in the list, which should be true when
#using inner_join as above
tipassScore = tipassScore %>%
  mutate(score_percentile = percent_rank(Score))

#only report results for the TPM sample (not the matrix samples)
tipassScore = tipassScore %>% dplyr::filter(!grepl("^X", Sample))

#Save results
write.table(tipassScore,"IPASS_Results.txt",sep="\t",row.names=F,quote=F)
