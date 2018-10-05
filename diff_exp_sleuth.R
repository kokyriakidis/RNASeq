#!/usr/bin/Rscript

#sleuth_object = diff_exp_sleuth("Control", "SOX15")

diff_exp_sleuth <- function(condition1, condition2) {

library(sleuth)
library(biomaRt)

tx2gene <- function(){
	mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
	t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id", "external_gene_name"), mart = mart)
	t2g <- dplyr::mutate(t2g, target_id = paste(ensembl_transcript_id, transcript_version, sep = "."))
	t2g <- dplyr::rename(t2g, target_id= target_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
	return(t2g)
	}

t2g <- tx2gene()

metadata <- read.delim("./Metadata/metadata.txt", comment.char="#")

# re-format metadata to meet sleuth's required format
# later on, the 'sleuth_prep' function of sleuth will assign 0=untreated and 1=treated condition
# importantly, by default such function will set '0' to first alphabetical match in the values
# found in the 'condition' field so here we force that 0=condition1
tab_metadata <- metadata
rownames(tab_metadata) <- tab_metadata$sample_id

conditions <- c()
for (i in 1:nrow(tab_metadata)) {
    tag <- tab_metadata[i, 'condition']
		if (tag == condition1) {
			conditions <- c(conditions, 0)
		}
		else if (tag == condition2) {
			conditions <- c(conditions, 1)
		}
		else {
			conditions <- c(conditions, -1)
		}
}
tab_metadata$condition <- conditions

# subset samples to include only those from the 2 conditions compared
cond1 <- tab_metadata$condition == 0
cond2 <- tab_metadata$condition == 1
tab_metadata <- tab_metadata[cond1 | cond2, ]

samples <- tab_metadata$sample_id

# get paths to data
tab_metadata$path <- as.character(tab_metadata$path)

# rename to meet sleuth's format requirements
tab_metadata <- tab_metadata[, c('sample_id', 'condition', 'path')]
names(tab_metadata)[1] <- 'sample'

# (1) load the kallisto processed data and make a regression model using 'condition' as the dependent variable
so <- sleuth_prep(tab_metadata, ~ condition, target_mapping = t2g)
# (2) estimate parameters for the sleuth response error measurement (full) model as responding to the 'condition' factor
so <- sleuth_fit(so)
# (3) Create another model where the gene expression is not dependent on any factor.
so <- sleuth_fit(so, ~1, 'reduced')

# (4.1) Run a likelihood ratio test (LRT) between the two models to see what transcripts appear
# to really be affected by the time factor value
so <- sleuth_lrt(so, 'reduced', 'full')

# (4.2) Run the Wald test (WT), a statistical tests which:
# - is somewhat related to the LRT and is also used to test for differential expression
# - LRT is considerd a better test than the WT but
# - WT is used becase it generates the beta statistic, which approximates to the fold change in expression between
# the 2 condition tested, which is typically reported in differential expression analysis
so <- sleuth_wt(so, paste0('condition'))

# export normalised abundance values
ANALYSIS <- "./Analysis"
condition1_name <- tolower(condition1)
condition2_name <- tolower(condition2)
otab = paste0(ANALYSIS, "/normalized_abundance_transcript_level_sleuth_", condition1_name, "_", condition2_name, ".tsv")
write.table(kallisto_table(so), otab, sep = "\t", quote = FALSE, row.names = FALSE)

# add beta (b), beta's standard error (se_b) and the mean expression in the samples (mean_obs)
res_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
res_wt <- sleuth_results(so, 'condition')
res <- merge(res_lrt, res_wt[, c('target_id', 'b', 'se_b', 'mean_obs')], on = 'target_id', sort = FALSE)

# export
condition1_name <- tolower(condition1)
condition2_name <- tolower(condition2)
otab = paste0(ANALYSIS, "/differential_expression_analysis_transcript_level_sleuth_", condition1_name, "_", condition2_name, ".tsv")
write.table(res, otab, sep = "\t", quote = FALSE, row.names = FALSE)

return(so)

}
