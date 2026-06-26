library(TwoSampleMR)
library(ggplot2)
library(foreach)
library(data.table)
library(tidyr) 

args <- commandArgs(trailingOnly = TRUE)
outcome_file <- args[1]
result_dir <- args[2]

if (!dir.exists(result_dir)) {
  dir.create(result_dir)
}

exposure_datas <- read_exposure_data(
  filename = "exposures.csv", 
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  pval_col = "pval.exposure",
)

exposure_datas$phenotype <- "BMI"
exposure_datas$samplesize.exposure <- 1122049

outcome_name <- gsub(".*/|\\.h\\.tsv\\.gz", "", outcome_file)

temp_file <- tempfile()
system(paste("gunzip -c", outcome_file, ">", temp_file))
outcome_data <- fread(temp_file, select = c("rsid", "beta", "standard_error", "effect_allele", "other_allele", "p_value","effect_allele_frequency"))
colnames(outcome_data) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "pval.outcome","eaf.outcome")
outcome_data <- outcome_data[complete.cases(outcome_data), ]

outcome_csv <- file.path(result_dir, paste0(outcome_name, "_outcome.csv"))
write.csv(outcome_data, outcome_csv, row.names = FALSE)

outcome_datas <- read_outcome_data(
  filename = outcome_csv,
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  eaf_col = "eaf.outcome",
  pval_col = "pval.outcome"
)
outcome_datas$samplesize.outcome <- 136016
outcome_datas$phenotype <- outcome_name

dat <- harmonise_data(exposure_dat = exposure_datas, outcome_dat = outcome_datas, action = 2)
file.remove(outcome_csv)

if (nrow(dat) > 0) {
  
  mr_results <- mr(dat)
  
  beta <- mr_results$b[3]
  se <- mr_results$se[3]
  lower_ci <- beta - 1.96 * se
  upper_ci <- beta + 1.96 * se
  pval <- mr_results$pval[3]
  num_snps <- nrow(dat)
  
  result <- data.frame(
    Outcome = outcome_name,
    Beta = beta,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    P_value = pval,
    Num_SNPs = num_snps
  )
  
  filename <- file.path(result_dir, outcome_name)
  if (!dir.exists(filename)) {
    dir.create(filename, showWarnings = FALSE)
  }
  
  write.table(dat, file = file.path(filename, "harmonise.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(generate_odds_ratios(mr_results), file = file.path(filename, "5_Method_OR.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(mr_pleiotropy_test(dat), file = file.path(filename, "pleiotropy.txt"), sep = "\t", quote = FALSE)
  write.table(mr_heterogeneity(dat), file = file.path(filename, "heterogeneity.txt"), sep = "\t", quote = FALSE)
  write.table(result, file = file.path(filename, "ivw_results.txt"), sep = "\t", quote = FALSE)
  
  ggsave(mr_scatter_plot(mr_results, dat)[[1]], file = file.path(filename, "scatter.pdf"), width = 8, height = 8)
  ggsave(mr_forest_plot(mr_singlesnp(dat))[[1]], file = file.path(filename, "forest.pdf"), width = 8, height = 8)
  ggsave(mr_funnel_plot(mr_singlesnp(dat))[[1]], file = file.path(filename, "funnelplot.pdf"), width = 8, height = 8)
  
} else {
  message("No harmonized data available for outcome: ", outcome_name)
}
