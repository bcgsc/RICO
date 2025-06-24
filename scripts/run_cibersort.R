library(base)
library(tidyr)
library(dplyr)
library(tibble)
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-s", "--sig_matrix"), default=NULL, type='character',
              help="Leukocyte signature matrix"),
  make_option(c("-m", "--mixture_file"), default=NULL, type='character',
              help="Expression matrix with hugo gene ids"),
  make_option(c("-p", "--perm"), default=0, type='character',
              help="Number of permutations; set to >=100 to calculate p-values (default = 0)"),
  make_option(c("-q", "--QN"), action="store_false", default=TRUE,
              help="Use -q if you do not want to quantile normalize the data"),
  make_option(c("-a", "--absolute"), action="store_true", default=FALSE,
              help="Run CIBERSORT in absolute mode (default = FALSE)\n- note that cell subsets will be scaled by their absolute levels and will not be represented as fractions (to derive the default output, normalize absolute levels such that they sum to 1 for each mixture sample)\n- the sum of all cell subsets in each mixture sample will be added to the ouput ('Absolute score'). If LM22 is used, this score will capture total immune content."),
  make_option(c("-c", "--abs_method"), default='sig.score', type='character',
              help="if absolute is set to TRUE, choose method: 'no.sumto1' or 'sig.score'\n- sig.score = for each mixture sample, define S as the median expression level of all genes in the signature matrix divided by the median expression level of all genes in the mixture. Multiple cell subset fractions by S.\n - no.sumto1 = remove sum to 1 constraint"),
  make_option(c("-o", "--out_dir"), default='.', type='character',
              help="The output directory"),
  make_option(c("-t", "--script_path"), type='character',
              help="The path to the CIBERSORT.R script"),
  make_option(c("-M", "--matrix_file"), type='character',
              help="The path to the control matrix for percentile calculations")
)

#parse the arguments
opt = parse_args(OptionParser(option_list=option_list))

#check that the mixture matrix uses hugo ids
signames = as.character(read.table(file=opt$sig_matrix, head=T, sep="\t")[, 1])
mixnames = as.character(read.table(file=opt$mixture_file, head=T, sep=" ")[, 1])
if(length(intersect(signames, mixnames))==0) {
	print("ERROR: The mixture file must use hugo id as the gene symbol in the first column. You are using:")
	print(head(mixnames))
	print("It needs to look like:")
	print(head(signames))
	q()
}

opt$QN

#Cibersort crashes when only 1 sample is provided.
TPMs = read.table(opt$mixture_file, sep = " ", header = TRUE)

#load the control matrix
M = read.table(opt$matrix_file, sep = " ", header = TRUE)

#add the control samples to the TPM matrix
#reformat with row_names
#uniqify the gene names.  THis is exactly what Cibersort does.
TPMs = TPMs %>%
  inner_join(M) %>%
  arrange(HGNC) %>%
  mutate_if(is.factor, as.character) %>%
  mutate_if(is.character, make.unique)

write.table(TPMs, "for_cibersort.csv", sep = "\t", quote = FALSE, row.names = FALSE )

source(opt$script_path)
results <- CIBERSORT(opt$sig_matrix, "for_cibersort.csv", opt$perm, opt$QN, opt$absolute, opt$abs_method)
write.csv(t(results), 'tmp.csv', row.names = TRUE)

myresults = tibble(read.csv('tmp.csv', sep = ","))
myresults = myresults %>%
  rename("cell_type" = "X")

#print the percentile of the new sample estimates in context of the BDAC samples
#from the provided matrix file
ptile_ests = myresults %>%
   pivot_longer(cols=starts_with("X"),
      names_to = "anon",
      values_to = "estimate") %>%
   group_by(cell_type) %>%
   summarize(ptile = 100 * sum(TPM > estimate) / n())

#remove cohort columns from matrix file before saving
res = myresults %>% select( - contains("X")) %>%
   inner_join(ptile_ests) %>%
   rename("Score" = TPM)

write.csv(res, 'immunedeconv_out.csv', row.names = FALSE)

#pretty formatted table to stdout
knitr::kable(res, digits=2)

#save results to a file
res = as.data.frame(res)
write.table(res, 'abundance_results.txt', quote = FALSE, sep = "\t", row.names = F)
