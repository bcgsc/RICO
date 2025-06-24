###Code to generate M1M2 scores
#Reads in a gene expression data frame of TPM values
#code can be run on any number of samples

#Modified on May 15, 2023 to fit in the RICO container

library(tidyr)
library(dplyr)

# gather the arguments
args = commandArgs(trailingOnl = TRUE)

# test if there are 3 arguments: if not, return an error
if (length(args) != 2) {
  stop("Two arguments are required", call. = FALSE)
}

TPM_file = args[1]
MATRIX_file = args[2]

#load TPM expression data - columns will be patients and rows hugo gene symbols
TPMs <- read.delim(TPM_file, sep = " ",
            header = TRUE,
            stringsAsFactors = FALSE)

#load the comparison matrix data
M = read.csv(MATRIX_file, sep = " ", header = TRUE)

#add the sample to the big matrix
TPMs = TPMs %>%
    inner_join(M)

# select the M1M2 genes from the matrix
# calculate the score per sample
M1M2_genes <- c("IDO1", "LAMP3", "CXCL10", "CXCL11", "APOL6","CXCL9", "PLA1A", "TNIP3", "CCL19", "CCR7")
TPMs = TPMs %>% filter(HGNC %in% M1M2_genes)
SCORES = TPMs %>% pivot_longer(cols=-HGNC, names_to = "anon", values_to = "estimate") %>%
    group_by(anon) %>%
    summarize(mean_exp = mean(estimate)) %>%
    mutate(mean_percentile = percent_rank(mean_exp)) %>%
    mutate(m1m2_score = if_else(mean_exp >= 41.29, 3,
                        if_else(mean_exp >= 22.84, 2,
                        if_else(mean_exp >= 6.96, 1, 0))))

#m1m2_score <- ifelse(lymph == "No" && m1m2 >= 41.29, 3, # 90th
#                ifelse(lymph == "Yes" && m1m2 >= 54.63, 2, # 80th
#                    ifelse(lymph == "No" && m1m2 >= 22.84, 2, # 80th
#                        ifelse(lymph == "No" && m1m2 >= 6.96, 1, # median
#                            0))))


#only report results for the TPM sample (not the matrix samples)
exp_scores = SCORES %>% dplyr::filter(!grepl("^X", anon))

#Save results
write.table(exp_scores,paste0(basename(TPM_file),"_M1M2.txt"), sep="\t", row.names=F, quote=F)
