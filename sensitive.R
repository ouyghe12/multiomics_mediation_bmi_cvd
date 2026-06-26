library(data.table)
library(glmnet)
library(dplyr)
library(DescTools)
library(survival)
library(stats)
library(dplyr)
library(car)

TARGET_PRO_RDS  <- "target_data.rds"
TARGET_NMR_RDS  <- "target_data_nmr.rds"
TARGET_COMB_RDS <- "combined_data.rds"

MED_TAB <- "/proj/sens2017538/proj_15152/ukb671783.tab"
SBP_TAB <- "/proj/sens2017538/proj_15152/ukb671784.tab"

OUTDIR <- "/home/ouyghe/bmi_m_p/results_sensitivity"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
OUT_LOG <- file.path(OUTDIR, "sensitivity_baseline_lipid_exclusion_plus_SBP_new.log")

sink(OUT_LOG, split = TRUE)
cat("==== Sensitivity analysis (baseline lipid-lowering exclusion + SBP covariate) ====\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

print_step <- function(step, n_before, n_after) {
  removed <- n_before - n_after
  cat(sprintf("[%s] %s: %d -> %d  removed=%d\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              step, n_before, n_after, removed))
}

read_rds_dt <- function(path) {
  if (!file.exists(path)) stop("RDS not found: ", path)
  x <- readRDS(path)
  setDT(x)
  x
}

get_lipid_exclude_eids <- function(med_dt) {

  cols <- c("f.6153.0.0","f.6153.0.1","f.6153.0.2","f.6153.0.3",
            "f.6177.0.0","f.6177.0.1","f.6177.0.2")

  cols_present <- intersect(cols, names(med_dt))
  cols_missing <- setdiff(cols, cols_present)

  cat("Medication columns present:", paste(cols_present, collapse = ", "), "\n")
  if (length(cols_missing) > 0) {
    cat("Medication columns missing (ignored):", paste(cols_missing, collapse = ", "), "\n")
  }

  if (!("f.eid" %in% names(med_dt))) stop("med_dt missing f.eid")
  if (length(cols_present) == 0) stop("No medication columns found in med_dt")


  for (cc in cols_present) med_dt[, (cc) := as.integer(get(cc))]

  expr <- paste0("(", paste(paste0("`", cols_present, "` == 1L"), collapse = " | "), ")")
  idx <- med_dt[eval(parse(text = expr)), unique(f.eid)]
  idx
}


prepare_dataset <- function(dt, sbp_dt, exclude_eids, label, drop_sbp_missing = TRUE) {
  cat("\n=============================\n")
  cat("Dataset:", label, "\n")

  if (!("f.eid" %in% names(dt))) stop(label, ": missing f.eid")
  if (!("f.eid" %in% names(sbp_dt))) stop("SBP table missing f.eid")
  if (!("f.4080.0.0" %in% names(sbp_dt))) stop("SBP table missing f.4080.0.0")

  n0 <- uniqueN(dt$f.eid)
  cat("Unique participants at start:", n0, "\n")

  n_before <- uniqueN(dt$f.eid)
  dt1 <- dt[!(f.eid %in% exclude_eids)]
  n_after <- uniqueN(dt1$f.eid)
  print_step(paste0(label, " exclude baseline lipid meds"), n_before, n_after)

  n_before2 <- uniqueN(dt1$f.eid)
  dt2 <- merge(dt1, sbp_dt, by = "f.eid", all.x = TRUE)
  n_after2 <- uniqueN(dt2$f.eid)
  print_step(paste0(label, " merge SBP"), n_before2, n_after2)

  if (drop_sbp_missing) {
    n_before3 <- uniqueN(dt2$f.eid)
    dt2 <- dt2[!is.na(f.4080.0.0)]
    n_after3 <- uniqueN(dt2$f.eid)
    print_step(paste0(label, " drop missing SBP"), n_before3, n_after3)
  }

  cat("Final unique participants:", uniqueN(dt2$f.eid), "\n")
  dt2
}


check_required_cols <- function(dt, cols, context = "") {
  missing <- setdiff(cols, names(dt))
  if (length(missing) > 0) {
    stop("Missing columns (", context, "): ", paste(missing, collapse = ", "))
  }
}

run_pair <- function(dt, outcome, time_var, status_var, mediator_rhs, covars_rhs) {
  cat("\n---------------------------------\n")
  cat("Outcome:", outcome, "\n")
  cat("Time var:", time_var, " | Status var:", status_var, "\n")
  cat("Mediator(s):", mediator_rhs, "\n")

  rhs_terms <- c("f.21001.0.0")  # BMI exposure
  cov_terms <- trimws(unlist(strsplit(covars_rhs, "\\+")))
  med_terms <- trimws(unlist(strsplit(mediator_rhs, "\\+")))

  needed <- unique(c("f.eid", time_var, status_var, rhs_terms, cov_terms, med_terms))
  check_required_cols(dt, needed, context = paste0("run_pair(", outcome, ")"))

  te_formula <- as.formula(
    paste0("Surv(", time_var, ", ", status_var, ") ~ f.21001.0.0 + ", covars_rhs)
  )
  de_formula <- as.formula(
    paste0("Surv(", time_var, ", ", status_var, ") ~ f.21001.0.0 + ", mediator_rhs, " + ", covars_rhs)
  )

  cat("TE:", deparse(te_formula), "\n")
  cat("DE:", deparse(de_formula), "\n")

  mf_te <- model.frame(te_formula, data = dt, na.action = na.omit)
  mf_de <- model.frame(de_formula, data = dt, na.action = na.omit)

  n_te <- nrow(mf_te)
  n_de <- nrow(mf_de)

  s_te <- mf_te[[1]]
  s_de <- mf_de[[1]]
  ev_te <- sum(s_te[, "status"] == 1)
  ev_de <- sum(s_de[, "status"] == 1)

  cat(sprintf("Usable N (TE): %d | Events (TE): %d\n", n_te, ev_te))
  cat(sprintf("Usable N (DE): %d | Events (DE): %d\n", n_de, ev_de))

  fit_te <- coxph(te_formula, data = dt, ties = "efron")
  fit_de <- coxph(de_formula, data = dt, ties = "efron")

  cat("\n[TE summary]\n"); print(summary(fit_te))
  cat("\n[DE summary]\n"); print(summary(fit_de))

  beta_te <- as.numeric(coef(fit_te)["f.21001.0.0"])
  beta_de <- as.numeric(coef(fit_de)["f.21001.0.0"])

  prop_med <- NA_real_
  if (!is.na(beta_te) && abs(beta_te) >= 1e-8) {
    prop_med <- (beta_te - beta_de) / beta_te
  }

  cat(sprintf("\nMediated proportion (log-HR): %s\n",
              ifelse(is.na(prop_med), "NA (beta_te ~ 0 or NA)", sprintf("%.2f%%", 100 * prop_med))))
  invisible(prop_med)
}

cat("Reading RDS datasets\n")
target_data     <- read_rds_dt(TARGET_PRO_RDS)
target_data_nmr <- read_rds_dt(TARGET_NMR_RDS)
combined_data   <- read_rds_dt(TARGET_COMB_RDS)

cat("Reading medication + SBP tab files\n")

med_header <- names(fread(MED_TAB, sep = "\t", nrows = 0))
med_cols_want <- c("f.eid",
                   "f.6153.0.0","f.6153.0.1","f.6153.0.2","f.6153.0.3",
                   "f.6177.0.0","f.6177.0.1","f.6177.0.2")
med_cols_use <- intersect(med_cols_want, med_header)

if (!("f.eid" %in% med_cols_use)) stop("MED_TAB has no f.eid column")
cat("Reading MED_TAB columns:", paste(med_cols_use, collapse = ", "), "\n")

medication_data <- fread(MED_TAB, sep = "\t", select = med_cols_use)
setDT(medication_data)

sbp_header <- names(fread(SBP_TAB, sep = "\t", nrows = 0))
sbp_cols_want <- c("f.eid","f.4080.0.0")
sbp_cols_use <- intersect(sbp_cols_want, sbp_header)
if (!all(sbp_cols_want %in% sbp_cols_use)) {
  stop("SBP_TAB missing required columns: ", paste(setdiff(sbp_cols_want, sbp_cols_use), collapse = ", "))
}
cat("Reading SBP_TAB columns:", paste(sbp_cols_use, collapse = ", "), "\n")

blood_pressure_data <- fread(SBP_TAB, sep = "\t", select = sbp_cols_use)
setDT(blood_pressure_data)

exclude_eids <- get_lipid_exclude_eids(medication_data)
cat("\nBaseline lipid-lowering users identified:", length(exclude_eids), "\n")

sens_pro  <- prepare_dataset(target_data,     blood_pressure_data, exclude_eids, "PRO",      drop_sbp_missing = TRUE)
sens_nmr  <- prepare_dataset(target_data_nmr, blood_pressure_data, exclude_eids, "NMR",      drop_sbp_missing = TRUE)
sens_comb <- prepare_dataset(combined_data,   blood_pressure_data, exclude_eids, "COMBINED", drop_sbp_missing = TRUE)

covars_rhs <- paste(
  c("f.4080.0.0", "f.31.0.0", "f.21003.0.0", "f.20116.0.0",
    "f.21000.0.0", "f.22040.0.0", "f.22189.0.0"),
  collapse = " + "
)

cat("\n\n==== Running sensitivity ====\n")

cat("\n\n######## PRO (mediator pBMI) ########\n")
run_pair(sens_pro,  "T2D",    "T2D_followup_time",    "T2D_status",    "pBMI",        covars_rhs)
run_pair(sens_pro,  "Stroke", "stroke_followup_time", "stroke_status", "pBMI",        covars_rhs)
run_pair(sens_pro,  "CAD",    "CAD_followup_time",    "CAD_status",    "pBMI",        covars_rhs)

cat("\n\n######## NMR (mediator mBMI) ########\n")
run_pair(sens_nmr,  "T2D",    "T2D_followup_time",    "T2D_status",    "mBMI",        covars_rhs)
run_pair(sens_nmr,  "Stroke", "stroke_followup_time", "stroke_status", "mBMI",        covars_rhs)
run_pair(sens_nmr,  "CAD",    "CAD_followup_time",    "CAD_status",    "mBMI",        covars_rhs)

cat("\n\n######## COMBINED (mediators mBMI + pBMI) ########\n")
run_pair(sens_comb, "T2D",    "T2D_followup_time",    "T2D_status",    "mBMI + pBMI", covars_rhs)
run_pair(sens_comb, "Stroke", "stroke_followup_time", "stroke_status", "mBMI + pBMI", covars_rhs)
run_pair(sens_comb, "CAD",    "CAD_followup_time",    "CAD_status",    "mBMI + pBMI", covars_rhs)

cat("\n\n==== DONE ====\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Log written to:", OUT_LOG, "\n")
sink()

