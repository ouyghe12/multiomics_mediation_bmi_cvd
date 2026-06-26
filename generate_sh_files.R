total_subsets <- 55

target_dir <- "sh_scripts"

project_id <- "your project"

working_dir <- "your folder"

modules_to_load <- c(
  "bioinfo-tools",
  "python",
  "R/4.3.1",
  "R_packages/4.3.1"
)

if (!dir.exists(target_dir)) dir.create(target_dir)

for (i in 1:total_subsets) {
  sh_file <- file.path(target_dir, paste0("run_subset", i, ".sh"))

  writeLines(c(
    "#!/bin/bash",
    paste0("#SBATCH -A ", project_id),
    "#SBATCH -p core",
    "#SBATCH -n 8",
    "#SBATCH -t 45:00:00",
    "#SBATCH --mem=64G",
    "#SBATCH --tmp=40G",
    paste0("#SBATCH -J MR_subset", i),
    paste0("#SBATCH -o ", working_dir, "/logs/MR_subset", i, ".out"),
    paste0("#SBATCH -e ", working_dir, "/logs/MR_subset", i, ".err"),
    "",
    "module load " , paste(modules_to_load, collapse = "\nmodule load "), 
    "",
    "cd $SNIC_TMP",
    "",
    paste0("cp ", working_dir, "/all_proteins.csv ."),
    paste0("cp ", working_dir, "/analysis_functions.R ."),
    paste0("cp ", working_dir, "/exposures.csv ."),
    paste0("cp ", working_dir, "/snp_map.csv ."),
    paste0("cp ", working_dir, "/subset_analysis.R ."),
    paste0("cp ", working_dir, "/token.txt ."),
    "",
    paste0("Rscript subset_analysis.R ", i, " ", total_subsets),
    "",
    paste0("cp -r results/* ", working_dir, "/results/"),
    "",
    paste0("if [ -f skipped_proteins.txt ]; then cp skipped_proteins.txt ", working_dir, "/results/skipped_proteins_", i, ".txt; fi"),
    "",
    "echo \"Subset ", i, " completed!\""
  ), con = sh_file)

}

submit_all <- file.path(target_dir, "submit_all.sh")
writeLines(c(
  "#!/bin/bash",
  "",
  "mkdir -p logs",
  "",
  "for sh_file in run_subset*.sh; do",
  "  echo \"Submitting $sh_file...\"",
  "  sbatch \"$sh_file\"",
  "  sleep 10",
  "done"
), con = submit_all)

cat("All sh scripts and submit_all.sh generated under ", target_dir, "\n")
