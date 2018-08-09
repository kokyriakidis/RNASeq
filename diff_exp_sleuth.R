#sleuth_object = diff_exp_sleuth("Control", "SOX15")


diff_exp_sleuth <- function(condition1, condition2) {

library(sleuth)

metadata <- read.delim("./Metadata/metadata.txt", comment.char="#")

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

cond1 <- tab_metadata$condition == 0
cond2 <- tab_metadata$condition == 1
tab_metadata <- tab_metadata[cond1 | cond2, ]

samples <- tab_metadata$sample_id

tab_metadata$path <- as.character(tab_metadata$path)

tab_metadata <- tab_metadata[, c('sample_id', 'condition', 'path')]
names(tab_metadata)[1] <- 'sample'

so <- sleuth_prep(tab_metadata, ~ condition)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
so <- sleuth_wt(so, paste0('condition'))

ANALYSIS <- "./Analysis"

condition1_name <- tolower(condition1)
condition2_name <- tolower(condition2)
otab = paste0(ANALYSIS, "/normalized_abundance_transcript_level_sleuth_", condition1_name, "_", condition2_name, ".tsv")
write.table(kallisto_table(so), otab, sep = "\t", quote = FALSE, row.names = FALSE)

res_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
res_wt <- sleuth_results(so, 'condition')
res <- merge(res_lrt, res_wt[, c('target_id', 'b', 'se_b', 'mean_obs')], on = 'target_id', sort = FALSE)

condition1_name <- tolower(condition1)
condition2_name <- tolower(condition2)
otab = paste0(ANALYSIS, "/differential_expression_analysis_transcript_level_sleuth_", condition1_name, "_", condition2_name, ".tsv")
write.table(res, otab, sep = "\t", quote = FALSE, row.names = FALSE)

return(so)

}
