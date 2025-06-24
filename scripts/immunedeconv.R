# Runs the immunedeconv code from the container downloaded by nextflow
# 3 params are expected:
# 1) path to space separated 2 column file
# 2) name of tool to use
# 3) Matrix of sample data to use for comparison against 1)

library(immunedeconv)
library(tibble)
library(tidyr)
library(dplyr)

# gather the arguments
args = commandArgs(trailingOnly=TRUE)

# test if there are 3 arguments: if not, return an error
if (length(args)!=3) { 
  stop("Three arguments are required", call.=FALSE)
}

TPM_file = args[1]
method = args[2]
MATRIX_file = args[3]

#load the control matrix
M = read.csv(MATRIX_file, sep=" ", header=TRUE)

#load the matrix of samples and assign row names
TPMs = read.csv(TPM_file, sep=" ", header=TRUE)
TPMs = TPMs %>% 
    inner_join(M)

#uniqify the gene names.  THis is exactly what Cibersort does.
TPMs = TPMs %>% 
  arrange(HGNC) %>%
  mutate_if(is.factor, as.character) %>%
  mutate_if(is.character, make.unique)

TPMs = TPMs %>% 
    remove_rownames %>% 
    column_to_rownames(var="HGNC") #%>%
#summary(TPMs)

#process with immunedeconv
res = deconvolute(TPMs, method)

#Show top results
knitr::kable(res, digits=2)

#print the percentile of the new sample estimates in context of the BDAC samples
#from the provided matrix file
ptile_ests = res %>%
   pivot_longer(cols=starts_with("X"),
      names_to = "anon",
      values_to = "estimate") %>%
   group_by(cell_type) %>%
   summarize(ptile = 100 * sum(TPM > estimate) / n())

#remove cohort columns from matrix file before saving
res = res %>% select( - contains("X")) %>% 
   inner_join(ptile_ests) %>%
   rename("Score" = TPM)
 summary(res)

write.csv(res, 'immunedeconv_out.csv', row.names = FALSE)
