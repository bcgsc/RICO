# Runs the immunedeconv code from the container downloaded by nextflow
# 2 params are expected:  
# 1) path to space separated 2 column file
# 2) CIBERSORT.R file path
# 3) Cibersort matrix file
# 4) Path to separate matrix of TPMs used to rank the results of 1)

library(immunedeconv)
library(tibble)
library(tidyr)
library(dplyr)

# gather the arguments
args = commandArgs(trailingOnly=TRUE)

# test if there are 4 arguments: if not, return an error
if (length(args)!=4) {
  stop("Four arguments are required", call.=FALSE)
}

TPM_file = args[1]
set_cibersort_binary(args[2])
set_cibersort_mat(args[3])
MATRIX_file = args[4]

#load the control matrix
M = read.csv(MATRIX_file, sep=" ", header=TRUE)

#load the sample to be processed and join with the control matrix
TPMs = read.csv(TPM_file, sep=" ", header=TRUE)
TPMs = TPMs %>% 
    inner_join(M)

#reformat with row_names
TPMs = TPMs %>%
    arrange(HGNC) %>%
    remove_rownames %>% 
    column_to_rownames(var="HGNC")
    
summary(TPMs)

#process with immunedeconv
res = deconvolute(TPMs, "cibersort_abs", QN=FALSE, abs_method="no.sumto1")
#show top of results
knitr::kable(res, digits=2)

#print the percentile of the new sample estimates in context of the BDAC samples
#from the provided matrix file
ptile_ests = res %>% 
  pivot_longer(cols=starts_with("X"), names_to = "anon", values_to = "estimate") %>%
  group_by(cell_type) %>%
  summarize(ptile = 100*sum(TPM>estimate)/n())

#remove cohort columns from matrix file before saving
res = res %>% select(-contains('X')) %>%
  inner_join(ptile_ests) %>%
  rename("Score"=TPM)
summary(res)
write.csv(res, 'immunedeconv_out.csv', row.names = FALSE)

