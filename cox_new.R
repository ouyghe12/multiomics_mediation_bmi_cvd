setwd("/home/ouyghe/bmi_m_p/")

library(data.table)
library(glmnet)
library(dplyr)
library(DescTools)
library(survival)
library(stats)
library(dplyr)
library(car)
library(parallel)
library(survival)
library(ggplot2)
library(ggpubr)

load("nmr_elastic_tuning_workspace.RData")

final_model_nmr <- final_model_09
scaling_attrs_nmr <- scaling_attrs
X_test = X_test
lambda_09 <- 0.0015848932
predicted <- predict(final_model_nmr, newx = X_test, s = lambda_09)
predicted <- as.numeric(predicted)

coef_09 <- coef(final_model_nmr, s = lambda_09)
coef_09_df <- data.frame(
  variable = rownames(coef_09),
  coefficient = as.numeric(coef_09),
  row.names = NULL
)

nonzero_coef_09_df <- subset(coef_09_df, coefficient != 0)
nonzero_coef_09_no_intercept <- subset(nonzero_coef_09_df, variable != "(Intercept)")
nonzero_count_09 <- nrow(nonzero_coef_09_no_intercept)

mse <- mean((predicted - y_test)^2)
mae <- mean(abs(predicted - y_test))
rsquared <- 1 - sum((predicted - y_test)^2) / sum((y_test - mean(y_test))^2)
rsq <- cor(predicted,y_test)^2
cat("MSE_09: ", mse, "\n")
cat("MAE_09: ", mae, "\n")
cat("R-squared_09: ", rsquared, "\n")
cat("R-squared-cor_09: ", rsq, "\n")
cat("non zero co 09:", nonzero_count_09)

basic_data <- fread("/proj/sens2017538/proj_15152/ukb671784.tab", 
                    sep = "\t", 
                    select = c("f.eid", "f.31.0.0", "f.21003.0.0", "f.48.0.0", "f.49.0.0", 
                               "f.20116.0.0", "f.21000.0.0"))
basic_data$f.31.0.0 <- factor(basic_data$f.31.0.0,levels = c(0,1),labels = c("Female","Male"))
basic_data$f.20116.0.0[basic_data$f.20116.0.0 == -3] <- NA
basic_data$f.20116.0.0 <- factor(basic_data$f.20116.0.0,
                                 levels = c(0, 1, 2),
                                 labels = c("Never", "Previous", "Current"))
basic_data$f.21000.0.0 <- floor(basic_data$f.21000.0.0 / 1000)
basic_data$f.21000.0.0[basic_data$f.21000.0.0 < 0] <- NA
basic_data$f.21000.0.0 <- factor(basic_data$f.21000.0.0,
                                 levels = c(1, 2, 3, 4, 5, 6),
                                 labels = c("White", "Mixed", "Asian", "Black", "Chinese", "Other"))

met_twnsend <- fread("/proj/sens2017538/proj_15152/ukb676391.tab", 
                     sep = "\t", 
                     select = c("f.eid", "f.22040.0.0","f.22189.0.0"))
target_data <- read.csv("target_data.csv")
target_data_nmr <- read.csv("target_data_nmr.csv")
final_model <- readRDS("elastic_net_model_pro.rds")
scaling_attrs <- readRDS("pro_adj_scaling.rds")
pt_medians <- readRDS("pt_train_median.rds")
rm_data <- read.csv("w15152_20250818.csv")
names(rm_data)[names(rm_data) == "X1001943"] <- "f.eid"
basic_data <- basic_data[!basic_data$f.eid %in% rm_data$f.eid, ]
met_twnsend <- met_twnsend[!met_twnsend$f.eid %in% rm_data$f.eid, ]

basic_data <- merge(basic_data, met_twnsend, by = "f.eid", all.x = TRUE, all.y = TRUE)

target_data_nmr <- merge(target_data_nmr, basic_data, by = "f.eid", all.x = TRUE)

target_data <- merge(target_data, basic_data, by = "f.eid", all.x = TRUE)

exclude_cols <- c("f.eid", "f.21001.0.0","f.130708.0.0", 
                  "f.131366.0.0", "CAD","f.40000.0.0",
                  "f.31.0.0", "f.21003.0.0","f.20116.0.0",
                  "f.21000.0.0", "f.22040.0.0", "f.22189.0.0",
                  "f.53.0.0","f.48.0.0","f.49.0.0")
                  
scale_cols <- setdiff(colnames(target_data_nmr), exclude_cols)

fill_na <- function(df, exclude_cols) {
  df[, !(colnames(df) %in% exclude_cols)] <- lapply(
    df[, !(colnames(df) %in% exclude_cols)], 
    function(col) {
      if (all(is.na(col))) {
        return(rep(0, length(col)))
      } else {
        return(ifelse(is.na(col), median(col, na.rm = TRUE), col))
      }
    }
  )
  return(df)
}

target_data_nmr <- fill_na(target_data_nmr, exclude_cols)

min_target_data_nmr <- min(as.matrix(target_data_nmr[, scale_cols, drop = FALSE]), na.rm = TRUE)

cat("Minimum value in target before log1p:", min_target_data_nmr, "\n")

if (min_target_data_nmr <= -1) {
  stop("Some values are <= -1, so log1p is not valid.")
}

target_data_nmr_log <- target_data_nmr

target_data_nmr_log[, scale_cols] <- log1p(target_data_nmr_log[, scale_cols, drop = FALSE])

target_data_nmr <- target_data_nmr_log

target_data_nmr[, scale_cols] <- scale(
  target_data_nmr[, scale_cols],
  center = scaling_attrs_nmr$center[scale_cols],
  scale = scaling_attrs_nmr$scale[scale_cols]
)

check_standardization <- data.frame(
  mean = colMeans(target_data_nmr[, scale_cols], na.rm = TRUE),
  var  = apply(target_data_nmr[, scale_cols], 2, var, na.rm = TRUE)
)

print(check_standardization)

X_target <- as.matrix(target_data_nmr[, scale_cols]) 
target_data_nmr$mBMI <- predict(final_model_nmr, X_target, s = final_model_nmr$lambda)
summary(target_data_nmr$mBMI)
summary(target_data_nmr$f.21001.0.0)

protein_col <- grep("_p$", colnames(target_data), value = TRUE)
pt_col <- gsub("_p$", "_pt", protein_col)

fill_pt_missing <- function(df, pt_medians) {
  for (col in names(pt_medians)) {
    if (col %in% colnames(df)) {
      df[[col]][is.na(df[[col]])] <- pt_medians[[col]]
    }
  }
  return(df)
}

target_data <- fill_pt_missing(target_data, pt_medians)

model_list <- readRDS("protein_lm_models.rds")

for (i in seq_along(protein_col)) {
  prot <- protein_col[i]
  pt <- pt_col[i]
  adj_name <- paste0(prot, "_adj")
  
  model <- model_list[[prot]]
  
  dat_target <- data.frame(x = target_data[[pt]])
  
  pred <- tryCatch({
    predict(model, newdata = dat_target)
  }, error = function(e) rep(NA, nrow(target_data)))
  
  target_data[[adj_name]] <- target_data[[prot]] - pred
}

adj_cols <- grep("_adj$", colnames(target_data), value = TRUE)

target_data[, adj_cols] <- scale(
  target_data[, adj_cols],
  center = scaling_attrs$center,
  scale = scaling_attrs$scale
)

get_protein_stats <- function(df, protein_cols) {
  stats <- data.frame(
    Protein = protein_cols,
    Mean = sapply(df[, protein_cols], mean, na.rm = TRUE),
    SD = sapply(df[, protein_cols], sd, na.rm = TRUE),
    Min = sapply(df[, protein_cols], min, na.rm = TRUE),
    Max = sapply(df[, protein_cols], max, na.rm = TRUE),
    NA_Count = sapply(df[, protein_cols], function(x) sum(is.na(x)))
  )
  return(stats)
}
target_stats <- get_protein_stats(target_data, adj_cols)
summary(target_stats$Mean)
summary(target_stats$SD)

target_data <- target_data[, !grepl("(_p$|_pt$)", colnames(target_data))]

X_target <- as.matrix(target_data[, adj_cols])
target_data$pBMI <- predict(final_model, newx = X_target, s = final_model$lambda)
summary(target_data$pBMI)
summary(target_data$f.21001.0.0)

target_data$baseline_date <- target_data$f.53.0.0
target_data <- target_data %>% select(-f.53.0.0)
target_data$end_date <- as.Date("2021-03-31") 
target_data$f.130708.0.0 <- as.Date(target_data$f.130708.0.0, format = "%Y-%m-%d")
target_data$f.131366.0.0 <- as.Date(target_data$f.131366.0.0, format = "%Y-%m-%d")
target_data$CAD <- as.Date(target_data$CAD, format = "%Y-%m-%d")
target_data$f.40000.0.0 <- as.Date(target_data$f.40000.0.0)
target_data$baseline_date <- as.Date(target_data$baseline_date)

target_data$T2D_followup_end <- pmin(target_data$f.130708.0.0, target_data$f.40000.0.0, target_data$end_date, na.rm = TRUE)
target_data$T2D_followup_time <- as.numeric(difftime(target_data$T2D_followup_end, target_data$baseline_date, units = "days")) / 30.44
target_data$T2D_status <- ifelse(!is.na(target_data$f.130708.0.0) & target_data$f.130708.0.0 <= target_data$end_date, 1, 0)

target_data$stroke_followup_end <- pmin(target_data$f.131366.0.0, target_data$f.40000.0.0, target_data$end_date, na.rm = TRUE)
target_data$stroke_followup_time <- as.numeric(difftime(target_data$stroke_followup_end, target_data$baseline_date, units = "days")) / 30.44
target_data$stroke_status <- ifelse(!is.na(target_data$f.131366.0.0) & target_data$f.131366.0.0 <= target_data$end_date, 1, 0)

target_data$CAD_followup_end <- pmin(target_data$CAD, target_data$f.40000.0.0, target_data$end_date, na.rm = TRUE)
target_data$CAD_followup_time <- as.numeric(difftime(target_data$CAD_followup_end, target_data$baseline_date, units = "days")) / 30.44
target_data$CAD_status <- ifelse(!is.na(target_data$CAD) & target_data$CAD <= target_data$end_date, 1, 0)

summary(target_data$T2D_followup_time)
table(target_data$T2D_status)
summary(target_data$stroke_followup_time)
table(target_data$stroke_status)
summary(target_data$CAD_followup_time)
table(target_data$CAD_status)

target_data_nmr$baseline_date <- target_data_nmr$f.53.0.0
target_data_nmr <- target_data_nmr %>% select(-f.53.0.0)
target_data_nmr$end_date <- as.Date("2021-03-31") 
target_data_nmr$f.130708.0.0 <- as.Date(target_data_nmr$f.130708.0.0, format = "%Y-%m-%d")
target_data_nmr$f.131366.0.0 <- as.Date(target_data_nmr$f.131366.0.0, format = "%Y-%m-%d")
target_data_nmr$CAD <- as.Date(target_data_nmr$CAD, format = "%Y-%m-%d")
target_data_nmr$f.40000.0.0 <- as.Date(target_data_nmr$f.40000.0.0)
target_data_nmr$baseline_date <- as.Date(target_data_nmr$baseline_date)

target_data_nmr$T2D_followup_end <- pmin(target_data_nmr$f.130708.0.0, target_data_nmr$f.40000.0.0, target_data_nmr$end_date, na.rm = TRUE)
target_data_nmr$T2D_followup_time <- as.numeric(difftime(target_data_nmr$T2D_followup_end, target_data_nmr$baseline_date, units = "days")) / 30.44
target_data_nmr$T2D_status <- ifelse(!is.na(target_data_nmr$f.130708.0.0) & target_data_nmr$f.130708.0.0 <= target_data_nmr$end_date, 1, 0)

target_data_nmr$stroke_followup_end <- pmin(target_data_nmr$f.131366.0.0, target_data_nmr$f.40000.0.0, target_data_nmr$end_date, na.rm = TRUE)
target_data_nmr$stroke_followup_time <- as.numeric(difftime(target_data_nmr$stroke_followup_end, target_data_nmr$baseline_date, units = "days")) / 30.44
target_data_nmr$stroke_status <- ifelse(!is.na(target_data_nmr$f.131366.0.0) & target_data_nmr$f.131366.0.0 <= target_data_nmr$end_date, 1, 0)

target_data_nmr$CAD_followup_end <- pmin(target_data_nmr$CAD, target_data_nmr$f.40000.0.0, target_data_nmr$end_date, na.rm = TRUE)
target_data_nmr$CAD_followup_time <- as.numeric(difftime(target_data_nmr$CAD_followup_end, target_data_nmr$baseline_date, units = "days")) / 30.44
target_data_nmr$CAD_status <- ifelse(!is.na(target_data_nmr$CAD) & target_data_nmr$CAD <= target_data_nmr$end_date, 1, 0)

summary(target_data_nmr$T2D_followup_time)
table(target_data_nmr$T2D_status)
summary(target_data_nmr$stroke_followup_time)
table(target_data_nmr$stroke_status)
summary(target_data_nmr$CAD_followup_time)
table(target_data_nmr$CAD_status)

target_data <- target_data %>%
  mutate(
    pBMI        = as.numeric(scale(rank(pBMI, ties.method = "average", na.last = "keep"))),
    `f.21001.0.0` = as.numeric(scale(rank(`f.21001.0.0`, ties.method = "average", na.last = "keep")))
  )

target_data_nmr <- target_data_nmr %>%
  mutate(
    mBMI        = as.numeric(scale(rank(mBMI, ties.method = "average", na.last = "keep"))),
    `f.21001.0.0` = as.numeric(scale(rank(`f.21001.0.0`, ties.method = "average", na.last = "keep")))
  )


new_columns <- setdiff(colnames(target_data), colnames(target_data_nmr))

combined_data <- inner_join(
  target_data_nmr,
  target_data[, c("f.eid", new_columns)],
  by = "f.eid"
)
summary(combined_data$T2D_followup_time)
table(combined_data$T2D_status)
summary(combined_data$stroke_followup_time)
table(combined_data$stroke_status)
summary(combined_data$CAD_followup_time)
table(combined_data$CAD_status)

cox_bmi_T2D <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                       f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                       f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                     data = target_data)
summary(cox_bmi_T2D)

cox_bmi_metabolic_T2D <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                                 pBMI + f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                                 f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                               data = target_data)
summary(cox_bmi_metabolic_T2D)

HR_TE_T2D <- coef(cox_bmi_T2D)["f.21001.0.0"]
HR_DE_T2D <- coef(cox_bmi_metabolic_T2D)["f.21001.0.0"]
IE_T2D <- (HR_TE_T2D - HR_DE_T2D) / (HR_TE_T2D)
print(paste0("metabolic to T2D IE: ", round(IE_T2D * 100, 2), "%"))

cox_bmi_T2D_nmr <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                           f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                           f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                         data = target_data_nmr)
summary(cox_bmi_T2D_nmr)

cox_bmi_metabolic_T2D_nmr <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                                     mBMI + f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                                     f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                   data = target_data_nmr)
summary(cox_bmi_metabolic_T2D_nmr)

HR_TE_T2D_nmr <- coef(cox_bmi_T2D_nmr)["f.21001.0.0"]
HR_DE_T2D_nmr <- coef(cox_bmi_metabolic_T2D_nmr)["f.21001.0.0"]
IE_T2D_nmr <- (HR_TE_T2D_nmr - HR_DE_T2D_nmr) / (HR_TE_T2D_nmr)
print(paste0("metabolic to T2D IE: ", round(IE_T2D_nmr * 100, 2), "%"))


cox_bmi_T2D_co <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                          f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                          f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                        data = combined_data)
summary(cox_bmi_T2D_co)

cox_bmi_metabolic_T2D_co <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                                    mBMI + pBMI + f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                                    f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                  data = combined_data)
summary(cox_bmi_metabolic_T2D_co)

HR_TE_T2D_co <- coef(cox_bmi_T2D_co)["f.21001.0.0"]
HR_DE_T2D_co <- coef(cox_bmi_metabolic_T2D_co)["f.21001.0.0"]
IE_T2D_co <- (HR_TE_T2D_co - HR_DE_T2D_co) / (HR_TE_T2D_co)
print(paste0("metabolic to T2D IE: ", round(IE_T2D_co * 100, 2), "%"))


cox_bmi_stroke <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                          f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                          f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                        data = target_data)
summary(cox_bmi_stroke)

cox_bmi_metabolic_stroke<- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                   pBMI + f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                                   f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                 data = target_data)
summary(cox_bmi_metabolic_stroke)

HR_TE_stroke <- coef(cox_bmi_stroke)["f.21001.0.0"]
HR_DE_stroke <- coef(cox_bmi_metabolic_stroke)["f.21001.0.0"]
IE_stroke <- (HR_TE_stroke - HR_DE_stroke) / (HR_TE_stroke)
print(paste0("metabolic to stroke IE: ", round(IE_stroke * 100, 2), "%"))

cox_bmi_stroke_nmr <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                              f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                              f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                            data = target_data_nmr)
summary(cox_bmi_stroke_nmr)

cox_bmi_metabolic_stroke_nmr <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                        mBMI + f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                                        f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                      data = target_data_nmr)
summary(cox_bmi_metabolic_stroke_nmr)

HR_TE_stroke_nmr <- coef(cox_bmi_stroke_nmr)["f.21001.0.0"]
HR_DE_stroke_nmr <- coef(cox_bmi_metabolic_stroke_nmr)["f.21001.0.0"]
IE_stroke_nmr <- (HR_TE_stroke_nmr - HR_DE_stroke_nmr) / (HR_TE_stroke_nmr)
print(paste0("metabolic to stroke IE: ", round(IE_stroke_nmr * 100, 2), "%"))


cox_bmi_CAD <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                       f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                       f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                     data = target_data)
summary(cox_bmi_CAD)

cox_bmi_metabolic_CAD <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                                 pBMI + f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                                 f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                               data = target_data)
summary(cox_bmi_metabolic_CAD)

HR_TE_CAD <- coef(cox_bmi_CAD)["f.21001.0.0"]
HR_DE_CAD <- coef(cox_bmi_metabolic_CAD)["f.21001.0.0"]
IE_CAD <- (HR_TE_CAD - HR_DE_CAD) / (HR_TE_CAD)
print(paste0("metabolic to CAD IE: ", round(IE_CAD * 100, 2), "%"))

cox_bmi_CAD_nmr <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                           f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                           f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                         data = target_data_nmr)
summary(cox_bmi_CAD_nmr)

cox_bmi_metabolic_CAD_nmr <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                                     mBMI + f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                                     f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                   data = target_data_nmr)
summary(cox_bmi_metabolic_CAD_nmr)

HR_TE_CAD_nmr <- coef(cox_bmi_CAD_nmr)["f.21001.0.0"]
HR_DE_CAD_nmr <- coef(cox_bmi_metabolic_CAD_nmr)["f.21001.0.0"]
IE_CAD_nmr <- (HR_TE_CAD_nmr - HR_DE_CAD_nmr) / (HR_TE_CAD_nmr)
print(paste0("metabolic to CAD IE: ", round(IE_CAD_nmr * 100, 2), "%"))



cox_bmi_stroke_co <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                              f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                              f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                            data = combined_data)
summary(cox_bmi_stroke_co)

cox_bmi_metabolic_stroke_co <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                        mBMI + pBMI + f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                                        f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                      data = combined_data)
summary(cox_bmi_metabolic_stroke_co)

HR_TE_stroke_co <- coef(cox_bmi_stroke_co)["f.21001.0.0"]
HR_DE_stroke_co <- coef(cox_bmi_metabolic_stroke_co)["f.21001.0.0"]
IE_stroke_co <- (HR_TE_stroke_co - HR_DE_stroke_co) / (HR_TE_stroke_co)
print(paste0("metabolic to stroke IE: ", round(IE_stroke_co * 100, 2), "%"))

cox_bmi_CAD_co <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                           f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                           f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                         data = combined_data)
summary(cox_bmi_CAD_co)

cox_bmi_metabolic_CAD_co <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                                     mBMI + pBMI + f.31.0.0 + f.21003.0.0 + f.20116.0.0 + 
                                     f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                   data = combined_data)
summary(cox_bmi_metabolic_CAD_co)

HR_TE_CAD_co <- coef(cox_bmi_CAD_co)["f.21001.0.0"]
HR_DE_CAD_co <- coef(cox_bmi_metabolic_CAD_co)["f.21001.0.0"]
IE_CAD_co <- (HR_TE_CAD_co - HR_DE_CAD_co) / (HR_TE_CAD_co)
print(paste0("metabolic to CAD IE: ", round(IE_CAD_co * 100, 2), "%"))

saveRDS(target_data_nmr, "target_data_nmr.rds")
saveRDS(combined_data, "combined_data.rds")

male <- subset(target_data, f.31.0.0 == "Male")
female <- subset(target_data, f.31.0.0 == "Female")

male_cox_T2D_total <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                              f.21003.0.0 + f.20116.0.0 + 
                              f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                            data = male)
summary(male_cox_T2D_total)
male_cox_T2D_direct <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                               pBMI + f.21003.0.0 + f.20116.0.0 + 
                               f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                             data = male)
summary(male_cox_T2D_direct)
male_T2D_HR_total <- coef(male_cox_T2D_total)["f.21001.0.0"]
male_T2D_HR_direct <- coef(male_cox_T2D_direct)["f.21001.0.0"]
male_T2D_HR_proportion <- (male_T2D_HR_total - male_T2D_HR_direct) / male_T2D_HR_total
print(paste0("metabolic to T2D IE in male: ", round(male_T2D_HR_proportion * 100, 2), "%"))

female_cox_T2D_total <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                                f.21003.0.0 + f.20116.0.0 + 
                                f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                              data = female)
summary(female_cox_T2D_total)
female_cox_T2D_direct <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                                 pBMI + f.21003.0.0 + f.20116.0.0 + 
                                 f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                               data = female)
summary(female_cox_T2D_direct)
female_T2D_HR_total <- coef(female_cox_T2D_total)["f.21001.0.0"]
female_T2D_HR_direct <- coef(female_cox_T2D_direct)["f.21001.0.0"]
female_T2D_HR_proportion <- (female_T2D_HR_total - female_T2D_HR_direct) / female_T2D_HR_total
print(paste0("metabolic to T2D IE in female: ", round(female_T2D_HR_proportion * 100, 2), "%"))

se_male_total <- summary(male_cox_T2D_total)$coefficients["f.21001.0.0", "se(coef)"]
se_male_direct <- summary(male_cox_T2D_direct)$coefficients["f.21001.0.0", "se(coef)"]
se_female_total <- summary(female_cox_T2D_total)$coefficients["f.21001.0.0", "se(coef)"]
se_female_direct <- summary(female_cox_T2D_direct)$coefficients["f.21001.0.0", "se(coef)"]

indirect_male <- male_T2D_HR_total - male_T2D_HR_direct
indirect_female <- female_T2D_HR_total - female_T2D_HR_direct

se_diff <- sqrt(se_male_total^2 + se_male_direct^2 + se_female_total^2 + se_female_direct^2)

z <- (indirect_male - indirect_female) / se_diff
p_value_HR_T2D <- 2 * pnorm(-abs(z))
p_value_HR_T2D

#stroke
male_cox_stroke_total <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                 f.21003.0.0 + f.20116.0.0 + 
                                 f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                               data = male)
summary(male_cox_stroke_total)
male_cox_stroke_direct <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                  pBMI + f.21003.0.0 + f.20116.0.0 + 
                                  f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                data = male)
summary(male_cox_stroke_direct)
male_stroke_HR_total <- coef(male_cox_stroke_total)["f.21001.0.0"]
male_stroke_HR_direct <- coef(male_cox_stroke_direct)["f.21001.0.0"]
male_stroke_HR_proportion <- (male_stroke_HR_total - male_stroke_HR_direct) / male_stroke_HR_total
print(paste0("metabolic to stroke IE in male: ", round(male_stroke_HR_proportion * 100, 2), "%"))

female_cox_stroke_total <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                   f.21003.0.0 + f.20116.0.0 + 
                                   f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                 data = female)
summary(female_cox_stroke_total)
female_cox_stroke_direct <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                    pBMI + f.21003.0.0 + f.20116.0.0 + 
                                    f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                  data = female)
summary(female_cox_stroke_direct)
female_stroke_HR_total <- coef(female_cox_stroke_total)["f.21001.0.0"]
female_stroke_HR_direct <- coef(female_cox_stroke_direct)["f.21001.0.0"]
female_stroke_HR_proportion <- (female_stroke_HR_total - female_stroke_HR_direct) / female_stroke_HR_total
print(paste0("metabolic to stroke IE in female: ", round(female_stroke_HR_proportion * 100, 2), "%"))

se_male_total <- summary(male_cox_stroke_total)$coefficients["f.21001.0.0", "se(coef)"]
se_male_direct <- summary(male_cox_stroke_direct)$coefficients["f.21001.0.0", "se(coef)"]
se_female_total <- summary(female_cox_stroke_total)$coefficients["f.21001.0.0", "se(coef)"]
se_female_direct <- summary(female_cox_stroke_direct)$coefficients["f.21001.0.0", "se(coef)"]

indirect_male <- male_stroke_HR_total - male_stroke_HR_direct
indirect_female <- female_stroke_HR_total - female_stroke_HR_direct

se_diff <- sqrt(se_male_total^2 + se_male_direct^2 + se_female_total^2 + se_female_direct^2)

z <- (indirect_male - indirect_female) / se_diff
p_value_HR_stroke <- 2 * pnorm(-abs(z))
p_value_HR_stroke


#CAD
male_cox_CAD_total <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                              f.21003.0.0 + f.20116.0.0 + 
                              f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                            data = male)
summary(male_cox_CAD_total)
male_cox_CAD_direct <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                               pBMI + f.21003.0.0 + f.20116.0.0 + 
                               f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                             data = male)
summary(male_cox_CAD_direct)
male_CAD_HR_total <- coef(male_cox_CAD_total)["f.21001.0.0"]
male_CAD_HR_direct <- coef(male_cox_CAD_direct)["f.21001.0.0"]
male_CAD_HR_proportion <- (male_CAD_HR_total - male_CAD_HR_direct) / male_CAD_HR_total
print(paste0("metabolic to CAD IE in male: ", round(male_CAD_HR_proportion * 100, 2), "%"))

female_cox_CAD_total <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                                f.21003.0.0 + f.20116.0.0 + 
                                f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                              data = female)
summary(female_cox_CAD_total)
female_cox_CAD_direct <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                                 pBMI + f.21003.0.0 + f.20116.0.0 + 
                                 f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                               data = female)
summary(female_cox_CAD_direct)
female_CAD_HR_total <- coef(female_cox_CAD_total)["f.21001.0.0"]
female_CAD_HR_direct <- coef(female_cox_CAD_direct)["f.21001.0.0"]
female_CAD_HR_proportion <- (female_CAD_HR_total - female_CAD_HR_direct) / female_CAD_HR_total
print(paste0("metabolic to CAD IE in female: ", round(female_CAD_HR_proportion * 100, 2), "%"))

se_male_total <- summary(male_cox_CAD_total)$coefficients["f.21001.0.0", "se(coef)"]
se_male_direct <- summary(male_cox_CAD_direct)$coefficients["f.21001.0.0", "se(coef)"]
se_female_total <- summary(female_cox_CAD_total)$coefficients["f.21001.0.0", "se(coef)"]
se_female_direct <- summary(female_cox_CAD_direct)$coefficients["f.21001.0.0", "se(coef)"]

indirect_male <- male_CAD_HR_total - male_CAD_HR_direct
indirect_female <- female_CAD_HR_total - female_CAD_HR_direct

se_diff <- sqrt(se_male_total^2 + se_male_direct^2 + se_female_total^2 + se_female_direct^2)

z <- (indirect_male - indirect_female) / se_diff
p_value_HR_CAD <- 2 * pnorm(-abs(z))
p_value_HR_CAD

table(male$T2D_status)
table(female$T2D_status)
table(male$stroke_status)
table(female$stroke_status)
table(male$CAD_status)
table(female$CAD_status)

#nmr
male <- subset(target_data_nmr, f.31.0.0 == "Male")
female <- subset(target_data_nmr, f.31.0.0 == "Female")

male_cox_T2D_total <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                              f.21003.0.0 + f.20116.0.0 + 
                              f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                            data = male)
summary(male_cox_T2D_total)
male_cox_T2D_direct <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                               mBMI + f.21003.0.0 + f.20116.0.0 + 
                               f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                             data = male)
summary(male_cox_T2D_direct)
male_T2D_HR_total <- coef(male_cox_T2D_total)["f.21001.0.0"]
male_T2D_HR_direct <- coef(male_cox_T2D_direct)["f.21001.0.0"]
male_T2D_HR_proportion <- (male_T2D_HR_total - male_T2D_HR_direct) / male_T2D_HR_total
print(paste0("metabolic to T2D IE in male: ", round(male_T2D_HR_proportion * 100, 2), "%"))

female_cox_T2D_total <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                                f.21003.0.0 + f.20116.0.0 + 
                                f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                              data = female)
summary(female_cox_T2D_total)
female_cox_T2D_direct <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                                 mBMI + f.21003.0.0 + f.20116.0.0 + 
                                 f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                               data = female)
summary(female_cox_T2D_direct)
female_T2D_HR_total <- coef(female_cox_T2D_total)["f.21001.0.0"]
female_T2D_HR_direct <- coef(female_cox_T2D_direct)["f.21001.0.0"]
female_T2D_HR_proportion <- (female_T2D_HR_total - female_T2D_HR_direct) / female_T2D_HR_total
print(paste0("metabolic to T2D IE in female: ", round(female_T2D_HR_proportion * 100, 2), "%"))

se_male_total <- summary(male_cox_T2D_total)$coefficients["f.21001.0.0", "se(coef)"]
se_male_direct <- summary(male_cox_T2D_direct)$coefficients["f.21001.0.0", "se(coef)"]
se_female_total <- summary(female_cox_T2D_total)$coefficients["f.21001.0.0", "se(coef)"]
se_female_direct <- summary(female_cox_T2D_direct)$coefficients["f.21001.0.0", "se(coef)"]

indirect_male <- male_T2D_HR_total - male_T2D_HR_direct
indirect_female <- female_T2D_HR_total - female_T2D_HR_direct

se_diff <- sqrt(se_male_total^2 + se_male_direct^2 + se_female_total^2 + se_female_direct^2)

z <- (indirect_male - indirect_female) / se_diff
p_value_HR_T2D <- 2 * pnorm(-abs(z))
p_value_HR_T2D

#stroke
male_cox_stroke_total <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                 f.21003.0.0 + f.20116.0.0 + 
                                 f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                               data = male)
summary(male_cox_stroke_total)
male_cox_stroke_direct <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                  mBMI + f.21003.0.0 + f.20116.0.0 + 
                                  f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                data = male)
summary(male_cox_stroke_direct)
male_stroke_HR_total <- coef(male_cox_stroke_total)["f.21001.0.0"]
male_stroke_HR_direct <- coef(male_cox_stroke_direct)["f.21001.0.0"]
male_stroke_HR_proportion <- (male_stroke_HR_total - male_stroke_HR_direct) / male_stroke_HR_total
print(paste0("metabolic to stroke IE in male: ", round(male_stroke_HR_proportion * 100, 2), "%"))

female_cox_stroke_total <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                   f.21003.0.0 + f.20116.0.0 + 
                                   f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                 data = female)
summary(female_cox_stroke_total)
female_cox_stroke_direct <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                    mBMI + f.21003.0.0 + f.20116.0.0 + 
                                    f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                  data = female)
summary(female_cox_stroke_direct)
female_stroke_HR_total <- coef(female_cox_stroke_total)["f.21001.0.0"]
female_stroke_HR_direct <- coef(female_cox_stroke_direct)["f.21001.0.0"]
female_stroke_HR_proportion <- (female_stroke_HR_total - female_stroke_HR_direct) / female_stroke_HR_total
print(paste0("metabolic to stroke IE in female: ", round(female_stroke_HR_proportion * 100, 2), "%"))

se_male_total <- summary(male_cox_stroke_total)$coefficients["f.21001.0.0", "se(coef)"]
se_male_direct <- summary(male_cox_stroke_direct)$coefficients["f.21001.0.0", "se(coef)"]
se_female_total <- summary(female_cox_stroke_total)$coefficients["f.21001.0.0", "se(coef)"]
se_female_direct <- summary(female_cox_stroke_direct)$coefficients["f.21001.0.0", "se(coef)"]

indirect_male <- male_stroke_HR_total - male_stroke_HR_direct
indirect_female <- female_stroke_HR_total - female_stroke_HR_direct

se_diff <- sqrt(se_male_total^2 + se_male_direct^2 + se_female_total^2 + se_female_direct^2)

z <- (indirect_male - indirect_female) / se_diff
p_value_HR_stroke <- 2 * pnorm(-abs(z))
p_value_HR_stroke


#CAD
male_cox_CAD_total <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                              f.21003.0.0 + f.20116.0.0 + 
                              f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                            data = male)
summary(male_cox_CAD_total)
male_cox_CAD_direct <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                               mBMI + f.21003.0.0 + f.20116.0.0 + 
                               f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                             data = male)
summary(male_cox_CAD_direct)
male_CAD_HR_total <- coef(male_cox_CAD_total)["f.21001.0.0"]
male_CAD_HR_direct <- coef(male_cox_CAD_direct)["f.21001.0.0"]
male_CAD_HR_proportion <- (male_CAD_HR_total - male_CAD_HR_direct) / male_CAD_HR_total
print(paste0("metabolic to CAD IE in male: ", round(male_CAD_HR_proportion * 100, 2), "%"))

female_cox_CAD_total <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                                f.21003.0.0 + f.20116.0.0 + 
                                f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                              data = female)
summary(female_cox_CAD_total)
female_cox_CAD_direct <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                                 mBMI + f.21003.0.0 + f.20116.0.0 + 
                                 f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                               data = female)
summary(female_cox_CAD_direct)
female_CAD_HR_total <- coef(female_cox_CAD_total)["f.21001.0.0"]
female_CAD_HR_direct <- coef(female_cox_CAD_direct)["f.21001.0.0"]
female_CAD_HR_proportion <- (female_CAD_HR_total - female_CAD_HR_direct) / female_CAD_HR_total
print(paste0("metabolic to CAD IE in female: ", round(female_CAD_HR_proportion * 100, 2), "%"))

se_male_total <- summary(male_cox_CAD_total)$coefficients["f.21001.0.0", "se(coef)"]
se_male_direct <- summary(male_cox_CAD_direct)$coefficients["f.21001.0.0", "se(coef)"]
se_female_total <- summary(female_cox_CAD_total)$coefficients["f.21001.0.0", "se(coef)"]
se_female_direct <- summary(female_cox_CAD_direct)$coefficients["f.21001.0.0", "se(coef)"]

indirect_male <- male_CAD_HR_total - male_CAD_HR_direct
indirect_female <- female_CAD_HR_total - female_CAD_HR_direct

se_diff <- sqrt(se_male_total^2 + se_male_direct^2 + se_female_total^2 + se_female_direct^2)

z <- (indirect_male - indirect_female) / se_diff
p_value_HR_CAD <- 2 * pnorm(-abs(z))
p_value_HR_CAD

table(male$T2D_status)
table(female$T2D_status)
table(male$stroke_status)
table(female$stroke_status)
table(male$CAD_status)
table(female$CAD_status)

#CO
male <- subset(combined_data, f.31.0.0 == "Male")
female <- subset(combined_data, f.31.0.0 == "Female")

male_cox_T2D_total <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                              f.21003.0.0 + f.20116.0.0 + 
                              f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                            data = male)
summary(male_cox_T2D_total)
male_cox_T2D_direct <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                               mBMI + pBMI + f.21003.0.0 + f.20116.0.0 + 
                               f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                             data = male)
summary(male_cox_T2D_direct)
male_T2D_HR_total <- coef(male_cox_T2D_total)["f.21001.0.0"]
male_T2D_HR_direct <- coef(male_cox_T2D_direct)["f.21001.0.0"]
male_T2D_HR_proportion <- (male_T2D_HR_total - male_T2D_HR_direct) / male_T2D_HR_total
print(paste0("metabolic to T2D IE in male: ", round(male_T2D_HR_proportion * 100, 2), "%"))

female_cox_T2D_total <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                                f.21003.0.0 + f.20116.0.0 + 
                                f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                              data = female)
summary(female_cox_T2D_total)
female_cox_T2D_direct <- coxph(Surv(T2D_followup_time, T2D_status) ~ f.21001.0.0 + 
                                 mBMI + pBMI + f.21003.0.0 + f.20116.0.0 + 
                                 f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                               data = female)
summary(female_cox_T2D_direct)
female_T2D_HR_total <- coef(female_cox_T2D_total)["f.21001.0.0"]
female_T2D_HR_direct <- coef(female_cox_T2D_direct)["f.21001.0.0"]
female_T2D_HR_proportion <- (female_T2D_HR_total - female_T2D_HR_direct) / female_T2D_HR_total
print(paste0("metabolic to T2D IE in female: ", round(female_T2D_HR_proportion * 100, 2), "%"))

se_male_total <- summary(male_cox_T2D_total)$coefficients["f.21001.0.0", "se(coef)"]
se_male_direct <- summary(male_cox_T2D_direct)$coefficients["f.21001.0.0", "se(coef)"]
se_female_total <- summary(female_cox_T2D_total)$coefficients["f.21001.0.0", "se(coef)"]
se_female_direct <- summary(female_cox_T2D_direct)$coefficients["f.21001.0.0", "se(coef)"]

indirect_male <- male_T2D_HR_total - male_T2D_HR_direct
indirect_female <- female_T2D_HR_total - female_T2D_HR_direct

se_diff <- sqrt(se_male_total^2 + se_male_direct^2 + se_female_total^2 + se_female_direct^2)

z <- (indirect_male - indirect_female) / se_diff
p_value_HR_T2D <- 2 * pnorm(-abs(z))
p_value_HR_T2D

#stroke
male_cox_stroke_total <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                 f.21003.0.0 + f.20116.0.0 + 
                                 f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                               data = male)
summary(male_cox_stroke_total)
male_cox_stroke_direct <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                  mBMI + pBMI + f.21003.0.0 + f.20116.0.0 + 
                                  f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                data = male)
summary(male_cox_stroke_direct)
male_stroke_HR_total <- coef(male_cox_stroke_total)["f.21001.0.0"]
male_stroke_HR_direct <- coef(male_cox_stroke_direct)["f.21001.0.0"]
male_stroke_HR_proportion <- (male_stroke_HR_total - male_stroke_HR_direct) / male_stroke_HR_total
print(paste0("metabolic to stroke IE in male: ", round(male_stroke_HR_proportion * 100, 2), "%"))

female_cox_stroke_total <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                   f.21003.0.0 + f.20116.0.0 + 
                                   f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                 data = female)
summary(female_cox_stroke_total)
female_cox_stroke_direct <- coxph(Surv(stroke_followup_time, stroke_status) ~ f.21001.0.0 + 
                                    mBMI + pBMI + f.21003.0.0 + f.20116.0.0 + 
                                    f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                                  data = female)
summary(female_cox_stroke_direct)
female_stroke_HR_total <- coef(female_cox_stroke_total)["f.21001.0.0"]
female_stroke_HR_direct <- coef(female_cox_stroke_direct)["f.21001.0.0"]
female_stroke_HR_proportion <- (female_stroke_HR_total - female_stroke_HR_direct) / female_stroke_HR_total
print(paste0("metabolic to stroke IE in female: ", round(female_stroke_HR_proportion * 100, 2), "%"))

se_male_total <- summary(male_cox_stroke_total)$coefficients["f.21001.0.0", "se(coef)"]
se_male_direct <- summary(male_cox_stroke_direct)$coefficients["f.21001.0.0", "se(coef)"]
se_female_total <- summary(female_cox_stroke_total)$coefficients["f.21001.0.0", "se(coef)"]
se_female_direct <- summary(female_cox_stroke_direct)$coefficients["f.21001.0.0", "se(coef)"]

indirect_male <- male_stroke_HR_total - male_stroke_HR_direct
indirect_female <- female_stroke_HR_total - female_stroke_HR_direct

se_diff <- sqrt(se_male_total^2 + se_male_direct^2 + se_female_total^2 + se_female_direct^2)

z <- (indirect_male - indirect_female) / se_diff
p_value_HR_stroke <- 2 * pnorm(-abs(z))
p_value_HR_stroke


#CAD
male_cox_CAD_total <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                              f.21003.0.0 + f.20116.0.0 + 
                              f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                            data = male)
summary(male_cox_CAD_total)
male_cox_CAD_direct <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                               mBMI + pBMI + f.21003.0.0 + f.20116.0.0 + 
                               f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                             data = male)
summary(male_cox_CAD_direct)
male_CAD_HR_total <- coef(male_cox_CAD_total)["f.21001.0.0"]
male_CAD_HR_direct <- coef(male_cox_CAD_direct)["f.21001.0.0"]
male_CAD_HR_proportion <- (male_CAD_HR_total - male_CAD_HR_direct) / male_CAD_HR_total
print(paste0("metabolic to CAD IE in male: ", round(male_CAD_HR_proportion * 100, 2), "%"))

female_cox_CAD_total <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                                f.21003.0.0 + f.20116.0.0 + 
                                f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                              data = female)
summary(female_cox_CAD_total)
female_cox_CAD_direct <- coxph(Surv(CAD_followup_time, CAD_status) ~ f.21001.0.0 + 
                                 mBMI + pBMI + f.21003.0.0 + f.20116.0.0 + 
                                 f.21000.0.0 + f.22040.0.0 + f.22189.0.0, 
                               data = female)
summary(female_cox_CAD_direct)
female_CAD_HR_total <- coef(female_cox_CAD_total)["f.21001.0.0"]
female_CAD_HR_direct <- coef(female_cox_CAD_direct)["f.21001.0.0"]
female_CAD_HR_proportion <- (female_CAD_HR_total - female_CAD_HR_direct) / female_CAD_HR_total
print(paste0("metabolic to CAD IE in female: ", round(female_CAD_HR_proportion * 100, 2), "%"))

se_male_total <- summary(male_cox_CAD_total)$coefficients["f.21001.0.0", "se(coef)"]
se_male_direct <- summary(male_cox_CAD_direct)$coefficients["f.21001.0.0", "se(coef)"]
se_female_total <- summary(female_cox_CAD_total)$coefficients["f.21001.0.0", "se(coef)"]
se_female_direct <- summary(female_cox_CAD_direct)$coefficients["f.21001.0.0", "se(coef)"]

indirect_male <- male_CAD_HR_total - male_CAD_HR_direct
indirect_female <- female_CAD_HR_total - female_CAD_HR_direct

se_diff <- sqrt(se_male_total^2 + se_male_direct^2 + se_female_total^2 + se_female_direct^2)

z <- (indirect_male - indirect_female) / se_diff
p_value_HR_CAD <- 2 * pnorm(-abs(z))
p_value_HR_CAD

table(male$T2D_status)
table(female$T2D_status)
table(male$stroke_status)
table(female$stroke_status)
table(male$CAD_status)
table(female$CAD_status)

bootstrap_mediation_cox <- function(data, outcome_time, outcome_status, bmi_var, mediator_var, covariates, n_boot = 1000) {
  mediation_proportions <- c()
  
  for (i in 1:n_boot) {
    boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
    
    formula_total <- as.formula(paste0("Surv(", outcome_time, ", ", outcome_status, ") ~ ", bmi_var, " + ", paste(covariates, collapse = " + ")))
    formula_direct <- as.formula(paste0("Surv(", outcome_time, ", ", outcome_status, ") ~ ", bmi_var, " + ", paste(mediator_var, collapse = " + "), " + ", paste(covariates, collapse = " + ")))
    
    cox_total <- coxph(formula_total, data = boot_data)
    cox_direct <- coxph(formula_direct, data = boot_data)
    
    HR_TE <- coef(cox_total)[bmi_var]
    HR_DE <- coef(cox_direct)[bmi_var]
    
    IE <- (HR_TE - HR_DE) / HR_TE
    mediation_proportions <- c(mediation_proportions, IE)
  }
  
  IE_mean <- mean(mediation_proportions)
  CI_low <- quantile(mediation_proportions, 0.025)
  CI_high <- quantile(mediation_proportions, 0.975)
  p_val <- 2 * min(mean(mediation_proportions <= 0), mean(mediation_proportions >= 0))
  
  return(data.frame(IE_mean = IE_mean, CI_low = CI_low, CI_high = CI_high, p_value = p_val))
}



datasets <- list(
  list(name = "mBMI_only", data = target_data_nmr, mediator = "mBMI"),
  list(name = "pBMI_only", data = target_data, mediator = "pBMI"),
  list(name = "combined", data = combined_data, mediator = c("pBMI", "mBMI"))
)

outcomes <- list(
  list(name = "T2D", time = "T2D_followup_time", status = "T2D_status"),
  list(name = "CAD", time = "CAD_followup_time", status = "CAD_status"),
  list(name = "stroke", time = "stroke_followup_time", status = "stroke_status")
)

covariates <- c("f.31.0.0", "f.21003.0.0", "f.20116.0.0", "f.21000.0.0", "f.22040.0.0", "f.22189.0.0")
bmi_var <- "f.21001.0.0"
n_cores <- 4
n_boot <- 1000


combinations <- expand.grid(dataset_index = 1:length(datasets), outcome_index = 1:length(outcomes))

results <- mclapply(1:nrow(combinations), function(i) {
  d <- datasets[[combinations$dataset_index[i]]]
  o <- outcomes[[combinations$outcome_index[i]]]
  
  res <- bootstrap_mediation_cox(
    data = d$data,
    outcome_time = o$time,
    outcome_status = o$status,
    bmi_var = bmi_var,
    mediator_var = d$mediator,
    covariates = covariates,
    n_boot = n_boot
  )
  
  res$dataset <- d$name
  res$mediator <- paste(d$mediator, collapse = " + ")
  res$disease <- o$name
  return(res)
}, mc.cores = n_cores)


final_df <- do.call(rbind, results)
print(final_df)

write.csv(final_df, "mediation_summary_all_datasets_new.csv", row.names = FALSE)


#sex_group

bootstrap_mediation_cox <- function(data, outcome_time, outcome_status, bmi_var, mediator_var, covariates, n_boot = 1000) {
  IE_list <- c()
  HR_TE_list <- c()
  HR_DE_list <- c()
  SE_TE_list <- c()
  SE_DE_list <- c()
  
  for (i in 1:n_boot) {
    boot_data <- data[sample(nrow(data), replace = TRUE), ]
    
    formula_total <- as.formula(paste0("Surv(", outcome_time, ", ", outcome_status, ") ~ ",
                                       bmi_var, " + ", paste(covariates, collapse = " + ")))
    
    formula_direct <- as.formula(paste0("Surv(", outcome_time, ", ", outcome_status, ") ~ ",
                                        bmi_var, " + ", paste(mediator_var, collapse = " + "), " + ",
                                        paste(covariates, collapse = " + ")))
    
    cox_total <- tryCatch(coxph(formula_total, data = boot_data), error = function(e) return(NULL))
    cox_direct <- tryCatch(coxph(formula_direct, data = boot_data), error = function(e) return(NULL))
    
    if (is.null(cox_total) || is.null(cox_direct)) next
    if (!(bmi_var %in% names(coef(cox_total))) || !(bmi_var %in% names(coef(cox_direct)))) next
    
    HR_TE <- coef(cox_total)[bmi_var]
    HR_DE <- coef(cox_direct)[bmi_var]
    IE <- (HR_TE - HR_DE) / HR_TE
    
    SE_TE <- summary(cox_total)$coefficients[bmi_var, "se(coef)"]
    SE_DE <- summary(cox_direct)$coefficients[bmi_var, "se(coef)"]
    
    IE_list <- c(IE_list, IE)
    HR_TE_list <- c(HR_TE_list, HR_TE)
    HR_DE_list <- c(HR_DE_list, HR_DE)
    SE_TE_list <- c(SE_TE_list, SE_TE)
    SE_DE_list <- c(SE_DE_list, SE_DE)
  }
  
  return(data.frame(
    IE_mean = mean(IE_list),
    CI_low = quantile(IE_list, 0.025),
    CI_high = quantile(IE_list, 0.975),
    p_value = 2 * min(mean(IE_list <= 0), mean(IE_list >= 0)),
    HR_TE = mean(HR_TE_list),
    HR_DE = mean(HR_DE_list),
    SE_TE = mean(SE_TE_list),
    SE_DE = mean(SE_DE_list)
  ))
}

z_test_IE_sex_diff <- function(male_res, female_res) {
  IE_male <- male_res$HR_TE - male_res$HR_DE
  IE_female <- female_res$HR_TE - female_res$HR_DE
  se_diff <- sqrt(male_res$SE_TE^2 + male_res$SE_DE^2 + female_res$SE_TE^2 + female_res$SE_DE^2)
  z <- (IE_male - IE_female) / se_diff
  p <- 2 * pnorm(-abs(z))
  return(p)
}

bmi_var <- "f.21001.0.0"
covariates <- c("f.21003.0.0", "f.20116.0.0", "f.21000.0.0", "f.22040.0.0", "f.22189.0.0")
n_boot <- 1000

datasets <- list(
  list(name = "mBMI_only", data = target_data_nmr, mediator = "mBMI"),
  list(name = "pBMI_only", data = target_data, mediator = "pBMI"),
  list(name = "combined", data = combined_data, mediator = c("pBMI", "mBMI"))
)

outcomes <- list(
  list(name = "T2D", time = "T2D_followup_time", status = "T2D_status"),
  list(name = "CAD", time = "CAD_followup_time", status = "CAD_status"),
  list(name = "stroke", time = "stroke_followup_time", status = "stroke_status")
)

results_all <- list()
sex_diff_results <- list()

for (d in datasets) {
  for (o in outcomes) {
    male_data <- subset(d$data, f.31.0.0 == "Male")
    female_data <- subset(d$data, f.31.0.0 == "Female")
    
    male_res <- bootstrap_mediation_cox(male_data, o$time, o$status, bmi_var, d$mediator, covariates, n_boot)
    female_res <- bootstrap_mediation_cox(female_data, o$time, o$status, bmi_var, d$mediator, covariates, n_boot)
    
    male_res$sex <- "Male"
    female_res$sex <- "Female"
    for (res in list(male_res, female_res)) {
      res$disease <- o$name
      res$dataset <- d$name
      res$mediator <- paste(d$mediator, collapse = " + ")
      results_all[[paste0(d$name, "_", o$name, "_", res$sex)]] <- res
    }
    
    p_sex <- z_test_IE_sex_diff(male_res, female_res)
    sex_diff_results[[paste0(d$name, "_", o$name)]] <- data.frame(
      disease = o$name,
      dataset = d$name,
      mediator = paste(d$mediator, collapse = " + "),
      p_value_sex_diff = p_sex
    )
  }
}

df_all_results <- do.call(rbind, results_all)
rownames(df_all_results) <- NULL
df_sex_diff <- do.call(rbind, sex_diff_results)
print(df_all_results)
print(df_sex_diff)

write.csv(df_all_results, "mediation_summary_sex_all_datasets_new.csv", row.names = FALSE)
write.csv(df_sex_diff, "mediation_sex_diff_pvalues_new.csv", row.names = FALSE)
