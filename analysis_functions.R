#library(reticulate)
#use_virtualenv("/domus/h1/ouyghe/.virtualenvs/r-synapser", required = TRUE)
library(synapser)
library(synapserutils)
library(parallel)
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ggplot2)
library(tidyr)

your_analysis <- function(data_dir) {
  message("Analyzing data in directory: ", data_dir)
  result_dir <- "./results"
  if (!dir.exists(result_dir)) {
    dir.create(result_dir)
  }
  
  result <- data.frame(
    Outcome = character(),
    Beta = numeric(),
    Lower_CI = numeric(),
    Upper_CI = numeric(),
    P_value = numeric(),
    Num_SNPs = numeric(),
    stringsAsFactors = FALSE
  )

  file_list <- list.files(data_dir, pattern = "\\.gz$", full.names = TRUE)
  if (length(file_list) == 0) {
    message("No .gz files found in ", data_dir)
    return(NULL)
  }
  
  message("Found ", length(file_list), " .gz files")
  
  gwas_data <- tryCatch({
    data.table::rbindlist(lapply(file_list, function(f) {
      message("Reading file: ", basename(f))
      data.table::fread(f)
    }), fill = TRUE)
  }, error = function(e) {
    message("Error reading GWAS data: ", e$message)
    return(NULL)
  })
  
  if (is.null(gwas_data) || nrow(gwas_data) == 0) {
    message("No data could be read from files in ", data_dir)
    return(NULL)
  }
  
  message("GWAS data loaded with ", nrow(gwas_data), " rows and ", ncol(gwas_data), " columns")
  message("Available columns: ", paste(colnames(gwas_data), collapse = ", "))
  
  if (!"ID" %in% colnames(gwas_data)) {
    message("Required column 'ID' not found in GWAS data")
    return(NULL)
  }
  
  snp_map <- tryCatch({
    read.csv("./snp_map.csv", header = TRUE)
  }, error = function(e) {
    message("Error reading SNP map: ", e$message)
    return(NULL)
  })
  
  if (is.null(snp_map) || nrow(snp_map) == 0) {
    message("SNP map is empty or could not be read")
    return(NULL)
  }
  
  
  gwas_data <- merge(snp_map, gwas_data, by = "ID", all.x = TRUE)
  
  if (!"rsid" %in% colnames(gwas_data)) {
    message("rsid column not found after merge")
    return(NULL)
  }
  
  mr_data <- tryCatch({
    gwas_data %>%
      dplyr::select(rsid, CHROM, POS19, POS38, REF, ALT, A1FREQ, BETA, SE, LOG10P, N)
  }, error = function(e) {
    stop("Error selecting MR columns: ", e$message)
  })
  
  exposure_datas <- tryCatch({
    read_exposure_data(
      filename = "exposures.csv",
      sep = ",",
      snp_col = "SNP",
      beta_col = "beta.exposure",
      se_col = "se.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure",
      eaf_col = "eaf.exposure",
      pval_col = "pval.exposure"
    )
  }, error = function(e) {
    message("Error reading exposure data: ", e$message)
    return(NULL)
  })
  
  if (is.null(exposure_datas) || nrow(exposure_datas) == 0) {
    message("No exposure data could be read")
    return(NULL)
  }
  
  exposure_datas$phenotype <- "BMI"
  exposure_datas$samplesize.exposure <- 1122049
  
  exposure_snps <- exposure_datas$SNP
  mr_data_filtered <- mr_data %>%
    dplyr::filter(rsid %in% exposure_snps) %>%
    dplyr::distinct(rsid, .keep_all = TRUE) %>%
    dplyr::mutate(pval = 10^(-LOG10P))

  if (nrow(mr_data_filtered) == 0) {
    message("No matching SNPs found between GWAS and exposure data")
    return(NULL)
  }
  
  tryCatch({
    write.csv(mr_data_filtered, "outcome.csv", row.names = FALSE)
    
    outcome_datas <- read_outcome_data(
      filename = "outcome.csv",
      sep = ",",
      snp_col = "rsid",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "ALT",
      other_allele_col = "REF",
      eaf_col = "A1FREQ",
      pval_col = "pval"
    )
    
    file.remove("outcome.csv")
  }, error = function(e) {
    message("Error processing outcome data: ", e$message)
    if (file.exists("outcome.csv")) file.remove("outcome.csv")
    return(NULL)
  })
  
  if (is.null(outcome_datas) || nrow(outcome_datas) == 0) {
    message("No outcome data could be processed")
    return(NULL)
  }
  
  dat <- tryCatch({
    harmonise_data(exposure_dat = exposure_datas, outcome_dat = outcome_datas, action = 2)
  }, error = function(e) {
    message("Error harmonizing data: ", e$message)
    return(NULL)
  })
  
  if (is.null(dat) || nrow(dat) == 0) {
    message("No harmonized data available")
    return(NULL)
  }
  
  if (nrow(dat) > 0) {
    dat$R2 <- (2 * (dat$beta.exposure^2) * dat$eaf.exposure * (1 - dat$eaf.exposure)) /
      (2 * (dat$beta.exposure^2) * dat$eaf.exposure * (1 - dat$eaf.exposure) +
         2 * dat$samplesize.exposure * dat$eaf.exposure * (1 - dat$eaf.exposure) * dat$se.exposure^2)
    dat$f <- dat$R2 * (dat$samplesize.exposure - 2) / (1 - dat$R2)
  }
  
  mr_results <- tryCatch({
    mr(dat)
  }, error = function(e) {
    message("Error in MR analysis: ", e$message)
    return(NULL)
  })
  
  if (is.null(mr_results) || nrow(mr_results) < 3) {
    message("MR analysis yielded insufficient results")
    return(NULL)
  }
  
  beta <- mr_results$b[3]
  se <- mr_results$se[3]
  lower_ci <- beta - 1.96 * se
  upper_ci <- beta + 1.96 * se
  pval <- mr_results$pval[3]
  num_snps <- nrow(dat)
  
  outcome_name <- basename(data_dir)

  single_result <- data.frame(
    Outcome = outcome_name,
    Beta = beta,
    Lower_CI = lower_ci,
    Upper_CI = upper_ci,
    P_value = pval,
    Num_SNPs = num_snps,
    stringsAsFactors = FALSE
  )
  
  out_dir <- file.path(result_dir, outcome_name)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE)
  }
  
  tryCatch({
    write.table(dat, file = file.path(out_dir, "harmonise.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(generate_odds_ratios(mr_results), file = file.path(out_dir, "OR.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(mr_pleiotropy_test(dat), file = file.path(out_dir, "pleiotropy.txt"), sep = "\t", quote = FALSE)
    write.table(mr_heterogeneity(dat), file = file.path(out_dir, "heterogeneity.txt"), sep = "\t", quote = FALSE)
    write.table(single_result, file = file.path(out_dir, "results.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
    
    ggsave(mr_scatter_plot(mr_results, dat)[[1]], file = file.path(out_dir, "scatter.pdf"), width = 8, height = 8)
    ggsave(mr_forest_plot(mr_singlesnp(dat))[[1]], file = file.path(out_dir, "forest.pdf"), width = 8, height = 8)
    ggsave(mr_leaveoneout_plot(mr_leaveoneout(dat))[[1]], file = file.path(out_dir, "sensitivity-analysis.pdf"), width = 8, height = 8)
    ggsave(mr_funnel_plot(mr_singlesnp(dat))[[1]], file = file.path(out_dir, "funnelplot.pdf"), width = 8, height = 8)
  }, error = function(e) {
    message("Error saving results: ", e$message)
  })
  
  return(single_result)
}

token_file <- "token.txt" # need to download your own token

if (!file.exists(token_file)) {
  stop("Token file not found: ", token_file)
}

token <- readLines(token_file, warn = FALSE)
token <- trimws(token[1])

synLogin(authToken = token)


download_and_analyze <- function(protein_info) {
  name     <- protein_info$name
  synid    <- protein_info$synid
  proc_dir <- "./processed_files"
  if (!dir.exists(proc_dir)) dir.create(proc_dir, recursive = TRUE)
  
  tar_path <- file.path(proc_dir, paste0(name, ".tar"))
  
  if (!file.exists(tar_path)) {
    message("Tar file for ", name, " not found. Downloading...")
    
    download_success <- FALSE
    
    tryCatch({
      syn_obj <- synapserutils::syncFromSynapse(synid, path = proc_dir)
      if (!is.null(syn_obj) && "path" %in% names(syn_obj) && file.exists(syn_obj$path)) {
        tar_path <- as.character(syn_obj$path)
        download_success <- TRUE
        message("Download successful with syncFromSynapse")
      }
    }, error = function(e) {
      message("syncFromSynapse failed: ", e$message)
    })
    
    if (!download_success) {
      tryCatch({
        message("Trying alternative download with synGet...")

        if (!requireNamespace("synapser", quietly = TRUE)) {
          message("synapser package is required but not available")
        } else {
          syn_file <- synapser::synGet(synid, downloadLocation = proc_dir, ifcollision = "overwrite.local")
          if (!is.null(syn_file) && "path" %in% names(syn_file) && file.exists(syn_file$path)) {
            tar_path <- as.character(syn_file$path)
            download_success <- TRUE
            message("Download successful with synGet")
          }
        }
      }, error = function(e) {
        message("synGet failed: ", e$message)
      })
    }
    
    if (!download_success) {
      tryCatch({
        message("Trying download with direct curl command...")
        
        auth_token <- NULL
        if (exists("synGetSessionToken", envir = asNamespace("synapser"))) {
          auth_token <- get("synGetSessionToken", envir = asNamespace("synapser"))()
        }
        
        if (!is.null(auth_token)) {
          download_url <- paste0("https://repo-prod.prod.sagebase.org/repo/v1/entity/", synid, "/file")
          
          curl_cmd <- paste0(
            "curl -L -k -s ", 
            "-H \"Authorization: Bearer ", auth_token, "\" ",
            "-o ", shQuote(tar_path), " ",
            shQuote(download_url)
          )
          
          system(curl_cmd)
          
          if (file.exists(tar_path) && file.info(tar_path)$size > 0) {
            download_success <- TRUE
            message("Download successful with curl command")
          }
        }
      }, error = function(e) {
        message("Curl download failed: ", e$message)
      })
    }
    
    if (!download_success) {
      message("All download methods failed for ", name)
      cat(name, file = "skipped_proteins.txt", sep = "\n", append = TRUE)
      return(NULL)
    }
  } else {
    message("Found existing tar file for ", name, ", skipping download.")
  }
  
  if (!file.exists(tar_path)) {
    message("Tar file still not found for ", name)
    cat(name, file = "skipped_proteins.txt", sep = "\n", append = TRUE)
    return(NULL)
  }
  
  if (file.info(tar_path)$size == 0) {
    message("Tar file is empty for ", name)
    cat(name, file = "skipped_proteins.txt", sep = "\n", append = TRUE)
    return(NULL)
  }
  
  fname      <- basename(tar_path)
  base_name  <- sub("\\.tar(?:\\.gz)?$", "", fname)
  local_dir  <- file.path(proc_dir, base_name)
  
  if (!dir.exists(local_dir)) {
    dir.create(local_dir, recursive = TRUE)
    message("Extracting ", fname, " into ", local_dir)
    
    extract_status <- tryCatch({
      strip_flag <- "--strip-components=1"
      gz_flag    <- if (grepl("\\.tar\\.gz$", tar_path)) "-xzf" else "-xf"
      
      result <- system2("tar", args = c(gz_flag, shQuote(tar_path), strip_flag, "-C", shQuote(local_dir)), 
                        stdout = TRUE, stderr = TRUE)

      if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
        message("Extraction failed with status ", attr(result, "status"), " for ", name)
        message("Error message: ", paste(result, collapse = "\n"))
        FALSE
      } else {
        TRUE
      }
    }, error = function(e) {
      message("Error during extraction of ", name, ": ", e$message)
      FALSE
    })
    
    if (!extract_status) {
      message("Failed to extract ", name)
      cat(name, file = "skipped_proteins.txt", sep = "\n", append = TRUE)
      return(NULL)
    }
    
    if (length(list.files(local_dir)) == 0) {
      message("Extracted directory is empty for ", name)
      cat(name, file = "skipped_proteins.txt", sep = "\n", append = TRUE)
      return(NULL)
    }
  } else {
    message("Protein ", name, " already extracted, skipping extraction.")
  }

  required_files <- c("./snp_map.csv", "./exposures.csv")
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    message("Missing required files: ", paste(missing_files, collapse = ", "))
    cat(name, ": missing files: ", paste(missing_files, collapse = ", "), 
        file = "skipped_proteins.txt", sep = "", append = TRUE)
    return(NULL)
  }

  result <- tryCatch({
    your_analysis(local_dir)
  }, error = function(e) {
    message("Error in analysis of ", name, ": ", e$message)
    cat(name, ": analysis error: ", e$message, file = "skipped_proteins.txt", sep = "", append = TRUE)
    NULL
  })
  
  message("Cleaning up extracted files and tar for protein: ", name)
  unlink(local_dir, recursive = TRUE)
  if (file.exists(tar_path)) {
    file.remove(tar_path)
  }
  
  return(result)
}


run_analysis_proteins <- function(protein_infos, n_cores = 2) {
  message("Starting analysis of ", nrow(protein_infos), " proteins with ", n_cores, " cores")

  cat("", file = "skipped_proteins.txt", append = FALSE)
  
  results_list <- parallel::mclapply(
    seq_len(nrow(protein_infos)),
    function(i) {
      message("\n--- Processing protein ", i, " of ", nrow(protein_infos), ": ", protein_infos$name[i], " ---")
      tryCatch(
        download_and_analyze(protein_infos[i, , drop = FALSE]),
        error = function(e) {
          message("Critical error processing ", protein_infos$name[i], ": ", e$message)
          cat(protein_infos$name[i], ": critical error: ", e$message, "\n", 
              file = "skipped_proteins.txt", append = TRUE)
          NULL
        }
      )
    },
    mc.cores = n_cores
  )
  
  successful <- sum(!sapply(results_list, is.null))
  failed <- nrow(protein_infos) - successful
  message("\nAnalysis complete: ", successful, " successful, ", failed, " failed")
  
  valid_results <- Filter(Negate(is.null), results_list)
  if (length(valid_results) == 0) {
    message("No valid results were produced")
    return(NULL)
  }
  
  combined_results <- do.call(rbind, valid_results)
  message("Final results contain ", nrow(combined_results), " rows")
  return(combined_results)
}


