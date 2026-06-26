setwd("/home/ouyghe/bmi_m_p/")

library(data.table)
library(dplyr)
library(DescTools)
library(survival)
library(stats)
library(car)
library(e1071)
library(glmnet)
library(caret)
library(missRanger)

set.seed(123)

source("/proj/sens2017538/proj_15152/ukb678544.r")  
nmr <- bd
rm(bd)

bmi_date <- fread("/proj/sens2017538/proj_15152/ukb671784.tab", 
                  sep = "\t", 
                  select = c("f.eid", "f.21001.0.0"))
bmi_date$f.21001.0.0 <- as.numeric(bmi_date$f.21001.0.0)

rm_data <- read.csv("w15152_20250818.csv")
names(rm_data)[names(rm_data) == "X1001943"] <- "f.eid"

disease_date <- fread("/proj/sens2017538/proj_15152/ukb671783.tab", 
                      sep = "\t", 
                      select = c("f.eid", "f.130708.0.0", "f.131296.0.0", "f.131298.0.0", "f.131300.0.0", "f.131302.0.0", "f.131304.0.0", "f.131306.0.0", "f.131366.0.0"))
disease_date$f.130708.0.0 <- as.Date(disease_date$f.130708.0.0)
disease_date$f.131296.0.0 <- as.Date(disease_date$f.131296.0.0)
disease_date$f.131298.0.0 <- as.Date(disease_date$f.131298.0.0)
disease_date$f.131300.0.0 <- as.Date(disease_date$f.131300.0.0)
disease_date$f.131302.0.0 <- as.Date(disease_date$f.131302.0.0)
disease_date$f.131304.0.0 <- as.Date(disease_date$f.131304.0.0)
disease_date$f.131306.0.0 <- as.Date(disease_date$f.131306.0.0)
disease_date$f.131366.0.0 <- as.Date(disease_date$f.131366.0.0)
disease_date$f.eid <- as.numeric(disease_date$f.eid)

mr_data <- read.csv("MR_nmr.csv")
mr_data <- mr_data %>% 
  filter(f_eid != "f.23468.0.0")
selected_columns <- mr_data$f_eid
selected_columns <- as.character(selected_columns)
nmr_selected <- nmr[, colnames(nmr) %in% selected_columns, drop = FALSE]
nmr_selected <- cbind(nmr["f.eid"],nmr_selected)

rm(nmr)

nmr_selected <- nmr_selected[!nmr_selected$f.eid %in% rm_data$f.eid, ]
bmi_date <- bmi_date[!bmi_date$f.eid %in% rm_data$f.eid, ]
disease_date <- disease_date[!disease_date$f.eid %in% rm_data$f.eid, ]

nmr_selected_filter <- nmr_selected[
  rowSums(!is.na(nmr_selected[, -which(names(nmr_selected) == "f.eid"), drop = FALSE])) > 0,
]

death <- fread("/proj/sens2017538/proj_15152/ukb671783.tab", 
               sep = "\t", 
               select = c("f.eid", "f.40000.0.0"))
death$f.40000.0.0 <- as.Date(death$f.40000.0.0)

baseline_date <- fread("/proj/sens2017538/proj_15152/ukb671784.tab", 
                       sep = "\t", 
                       select = c("f.eid", "f.53.0.0"))
baseline_date$f.53.0.0 <- as.Date(baseline_date$f.53.0.0)

merged_data <- merge(bmi_date, disease_date, by = "f.eid")
merged_data <- merge(merged_data, nmr_selected_filter, by = "f.eid")
merged_data <- merge(merged_data, death, by = "f.eid")
merged_data <- merge(merged_data, baseline_date, by = "f.eid")

cat("Initial N after merge:", nrow(merged_data), "\n")

n_before <- nrow(merged_data)
TD2_before_2010 <- subset(merged_data,merged_data$f.130708.0.0 <= as.Date(merged_data$f.53.0.0) & !is.na(merged_data$f.130708.0.0))
cat("T2D before baseline excluded:", nrow(TD2_before_2010), "\n")
merged_data <- merged_data %>%
  anti_join(TD2_before_2010, by = "f.eid")
cat("Remaining after excluding prevalent T2D:", nrow(merged_data), "\n")
cat("Removed in this step:", n_before - nrow(merged_data), "\n\n")


n_before <- nrow(merged_data)
stroke_before_2010 <- subset(merged_data,merged_data$f.131366.0.0 <= as.Date(merged_data$f.53.0.0) & !is.na(merged_data$f.131366.0.0))
cat("Stroke before baseline excluded:", nrow(stroke_before_2010), "\n")
merged_data <- merged_data %>%
  anti_join(stroke_before_2010, by = "f.eid")
cat("Remaining after excluding prevalent stroke:", nrow(merged_data), "\n")
cat("Removed in this step:", n_before - nrow(merged_data), "\n\n")

cad_cols <- c("f.131296.0.0","f.131298.0.0", "f.131300.0.0", "f.131302.0.0", "f.131304.0.0", "f.131306.0.0")
merged_data <- as.data.frame(merged_data)
merged_data$CAD <- apply(merged_data[, cad_cols], 1, min, na.rm = TRUE)
merged_data <- merged_data[ , !(names(merged_data) %in% cad_cols)]
n_before <- nrow(merged_data)
CAD_before_2010 <- subset(merged_data,merged_data$CAD <= as.Date(merged_data$f.53.0.0) & !is.na(merged_data$CAD))
cat("CAD before baseline excluded:", nrow(CAD_before_2010), "\n")
merged_data <- merged_data %>%
  anti_join(CAD_before_2010, by = "f.eid")
  
cat("Remaining after excluding prevalent CAD:", nrow(merged_data), "\n")
cat("Removed in this step:", n_before - nrow(merged_data), "\n\n")


cat("Final N after all baseline disease exclusions:", nrow(merged_data), "\n")

columns_to_remove <- c("f.130708.0.0", "f.131366.0.0", "CAD","f.53.0.0","f.40000.0.0")
ml_data <- merged_data[, !(colnames(merged_data) %in% columns_to_remove)]

ml_data <- ml_data %>%
  filter(!is.na(f.21001.0.0))
cat("N with non-missing BMI:", nrow(ml_data), "\n")

ml_data$BMI_group <- factor(
  cut(
    ml_data$f.21001.0.0,
    breaks = c(0, 18.5, 24.9, 29.9, Inf),
    labels = c("0", "1", "2", "3")
  ),
  levels = c("0", "1", "2", "3")
)
bmi_table <- table(ml_data$BMI_group, useNA = "no")
print(bmi_table)

bmi_counts <- as.numeric(bmi_table)
names(bmi_counts) <- names(bmi_table) 
total_samples <- sum(bmi_counts)  
bmi_proportions <- bmi_counts / total_samples  

target_n <- 5000
sample_sizes <- round(target_n * bmi_proportions)  
adjustment <- target_n - sum(sample_sizes)  
sample_sizes[which.max(sample_sizes)] <- sample_sizes[which.max(sample_sizes)] + adjustment  
final_sample_sizes <- setNames(sample_sizes, names(bmi_counts))
print(final_sample_sizes)


ml_test_data <- ml_data %>%
  filter(!is.na(BMI_group)) %>%
  group_by(BMI_group) %>%
  sample_n(size = min(n(), final_sample_sizes[as.character(first(BMI_group))])) %>%
  ungroup() %>%
  select(-BMI_group)

ml_data <- ml_data %>%
  select(-BMI_group)
target_data<- anti_join(merged_data, ml_test_data, by = "f.eid")

write.csv(ml_test_data, "ml_data_nmr.csv", row.names = FALSE) 
write.csv(target_data, "target_data_nmr.csv", row.names = FALSE)

exclude_cols <- c("f.eid", "f.21001.0.0")
cols_to_transform <- setdiff(colnames(ml_test_data), exclude_cols)

set.seed(123)
n <- nrow(ml_test_data)
train_index <- sample(1:n, size = floor(0.7 * n))
test_index <- setdiff(1:n, c(train_index))

train <- ml_test_data[train_index, c(cols_to_transform, exclude_cols), drop = FALSE]
test <- ml_test_data[test_index, c(cols_to_transform, exclude_cols), drop = FALSE]

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

train <- fill_na(train, exclude_cols)
test <- fill_na(test, exclude_cols)

min_train <- min(as.matrix(train[, cols_to_transform, drop = FALSE]), na.rm = TRUE)
min_test  <- min(as.matrix(test[, cols_to_transform, drop = FALSE]), na.rm = TRUE)

cat("Minimum value in train before log1p:", min_train, "\n")
cat("Minimum value in test before log1p:", min_test, "\n")

if (min_train <= -1 || min_test <= -1) {
  stop("Some values are <= -1, so log1p is not valid.")
}

train_log <- train
test_log  <- test

train_log[, cols_to_transform] <- log1p(train_log[, cols_to_transform, drop = FALSE])
test_log[, cols_to_transform]  <- log1p(test_log[, cols_to_transform, drop = FALSE])

write.csv(
  data.frame(variable = cols_to_transform, transform = "log1p_all"),
  "log1p_plan_from_train.csv",
  row.names = FALSE
)

cols_to_scale <- setdiff(colnames(train_log), exclude_cols)
scaled_train <- scale(train_log[, cols_to_scale])
scaling_attrs <- list(
  center = attr(scaled_train, "scaled:center"),
  scale = attr(scaled_train, "scaled:scale")
)
train_log[, cols_to_scale] <- as.data.frame(scaled_train)

test_log[, cols_to_scale] <- scale(
  test_log[, cols_to_scale],
  center = scaling_attrs$center,
  scale = scaling_attrs$scale
)

saveRDS(scaling_attrs, "nmr_scaling_params_missRanger.rds")

train <- train_log
test <- test_log
low_variance_cols <- which(apply(train[, cols_to_transform], 2, var) < 1e-6)
if (length(low_variance_cols) > 0) {
  message("Low variance columns detected: ", paste(colnames(train)[low_variance_cols], collapse = ", "))
} else {
  message("No low variance columns detected.")
}

train <- as.data.frame(train)
test <- as.data.frame(test)

X_train <- as.matrix(train[, !(colnames(train) %in% c("f.eid", "f.21001.0.0"))])
y_train <- train$f.21001.0.0

cor_matrix <- cor(train, use = "pairwise.complete.obs")
model <- lm(y_train ~ ., data = train) 
vif_values <- vif(model)
if(any(vif_values > 10)){
  print("VIF > 10")
} else if(any(vif_values > 5)){
  print("VIF > 5")
} else {
  print("VIF <= 5")
}

X_test <- as.matrix(test[, !(colnames(test) %in% c("f.eid", "f.21001.0.0"))])
y_test <- test$f.21001.0.0

# baseline
baseline_mse <- mean((y_test - mean(y_train))^2)
total_variance <- mean((y_test - mean(y_test))^2)

baseline_r2 <- 1 - (baseline_mse / total_variance)
cat("Baseline MSE:", baseline_mse, "\n")
cat("Baseline R2:", baseline_r2, "\n")

# lambda test
ridge_cv <- cv.glmnet(X_train, y_train, alpha = 0)
lasso_cv <- cv.glmnet(X_train, y_train, alpha = 1)
elastic_cv <- cv.glmnet(X_train, y_train, alpha = 0.8)

best_lambda_ridge <- ridge_cv$lambda.min
best_lambda_lasso <- lasso_cv$lambda.min
best_lambda_elastic <- elastic_cv$lambda.min

cat("Ridge best λ:", best_lambda_ridge, "\n")
cat("LASSO best λ:", best_lambda_lasso, "\n")
cat("Elastic best λ:", best_lambda_elastic, "\n")

# ridge model
ridge_model <- glmnet(X_train, y_train, family = "gaussian", alpha = 0)
print(ridge_model)
plot(ridge_model, label = TRUE)
plot(ridge_model, xvar = "lambda", label = TRUE)
ridge_coef <- predict(ridge_model, s = best_lambda_ridge, type = "coefficients")
print(ridge_coef)
plot(ridge_model, xvar = "dev", label = TRUE)

ridge_predictions <- predict(ridge_model, newx = X_test, s = best_lambda_ridge, type = "response")
plot(ridge_predictions, y_test, 
     xlab = "Predicted BMI", ylab = "Actual BMI", 
     main = "Ridge Regression: Predicted vs Actual", col = "blue", pch = 19)
abline(0, 1, col = "red")
ridge_resid <- ridge_predictions - y_test
ridge_mse <- mean(ridge_resid^2)
cat("Ridge Regression (MSE):", ridge_mse)
ridge_mae <- mean(abs(ridge_resid))
cat("Ridge Regression (MAE):", ridge_mae)
ridge_ss_total <- sum((y_test - mean(y_test))^2)
ridge_ss_res <-sum(ridge_resid^2)
ridge_r_squared <- 1 - (ridge_ss_res / ridge_ss_total)
cat("Ridge Regression (r_squared):", ridge_r_squared)
ridge_n <- 1000
ridge_p <- ncol(X_test)
ridge_adjusted_r_squared <- 1 - (1 - ridge_r_squared) * ((ridge_n -1 ) / (ridge_n - ridge_p - 1))
cat("Ridge Regression (adjusted_r_squared):", ridge_adjusted_r_squared)

#lasso model
lasso_model <- glmnet(X_train, y_train, family = "gaussian", alpha = 1)
print(lasso_model)
plot(lasso_model, xvar = "lambda", label = TRUE)
lasso_coef <- predict(lasso_model, s = best_lambda_lasso, type = "coefficients")
lasso_coef
lasso_predictions <- predict(lasso_model, newx = X_test, s = best_lambda_lasso, type = "response")
plot(lasso_predictions, y_test, 
     xlab = "Predicted BMI", ylab = "Actual BMI", 
     main = "Lasso Regression: Predicted vs Actual", col = "blue", pch = 19)
abline(0, 1, col = "red")
lasso_resid <- lasso_predictions - y_test
lasso_mse <- mean(lasso_resid^2)
cat("Lasso Regression (MSE):", lasso_mse)
lasso_mae <- mean(abs(lasso_resid))
cat("Lasso Regression (MAE):", lasso_mae)
lasso_ss_total <- sum((y_test - mean(y_test))^2)
lasso_ss_res <-sum(lasso_resid^2)
lasso_r_squared <- 1 - (lasso_ss_res / lasso_ss_total)
cat("Lasso Regression (r_squared):", lasso_r_squared)
lasso_n <- 1000
lasso_p <- ncol(X_test)
lasso_adjusted_r_squared <- 1 - (1 - lasso_r_squared) * ((lasso_n -1 ) / (lasso_n - lasso_p - 1))
cat("Lasso Regression (adjusted_r_squared):", lasso_adjusted_r_squared)

nonzero_count <- sum(lasso_coef != 0) - 1 
nonzero_count

# elastic net
set.seed(123)
log_lambda_range <- seq(-4, 0.2, by = 0.1)  
lambda_range <- 10^log_lambda_range        
alpha_values <- seq(0.1, 0.9, by = 0.1) 
grid <- expand.grid(.alpha = alpha_values, .lambda = lambda_range)
control <- trainControl(method = "cv",number = 10) 
train <- train %>% select(-f.eid)
enet_train <- train(
  f.21001.0.0 ~ ., 
  data = train, 
  method = "glmnet", 
  trControl = control, 
  tuneGrid = grid
)
print(enet_train$bestTune)
print(enet_train$results)

write.csv(enet_train$results, file = "enet_train_results_nmr.csv", row.names = FALSE)

saveRDS(enet_train, file = "enet_train_nmr.rds")


best_alpha <- enet_train$bestTune$alpha
best_lambda <- enet_train$bestTune$lambda

final_model_best <- glmnet(
  x = X_train, 
  y = y_train, 
  alpha = best_alpha, 
  lambda = best_lambda, 
  family = "gaussian"
)

predicted <- predict(final_model_best, newx = X_test, s = best_lambda)
coef(final_model_best)
mse <- mean((predicted - y_test)^2)
mae <- mean(abs(predicted - y_test))
rsquared <- 1 - sum((predicted - y_test)^2) / sum((y_test - mean(y_test))^2)
rsq <- cor(predicted,y_test)^2
cat("MSE_best: ", mse, "\n")
cat("MAE_best: ", mae, "\n")
cat("R-squared_best: ", rsquared, "\n")
cat("R-squared-cor_best: ", rsq, "\n")
nonzero_count_best <- sum(predicted != 0) - 1 
cat("non zero co best:", nonzero_count_best)



alpha_09 <- 0.9
lambda_09 <- 0.0015848932

final_model_09 <- glmnet(
  x = X_train, 
  y = y_train, 
  alpha = alpha_09, 
  lambda = lambda_09, 
  family = "gaussian"
)

predicted <- predict(final_model_09, newx = X_test, s = lambda_09)
coef(final_model_09)
mse <- mean((predicted - y_test)^2)
mae <- mean(abs(predicted - y_test))
rsquared <- 1 - sum((predicted - y_test)^2) / sum((y_test - mean(y_test))^2)
rsq <- cor(predicted,y_test)^2
cat("MSE_09: ", mse, "\n")
cat("MAE_09: ", mae, "\n")
cat("R-squared_09: ", rsquared, "\n")
cat("R-squared-cor_09: ", rsq, "\n")
nonzero_count_09 <- sum(predicted != 0) - 1 
cat("non zero co 09:", nonzero_count_09)


alpha_08 <- 0.8
lambda_08 <- 0.0012589254

final_model_08 <- glmnet(
  x = X_train, 
  y = y_train, 
  alpha = alpha_08, 
  lambda = lambda_08, 
  family = "gaussian"
)

predicted <- predict(final_model_08, newx = X_test, s = lambda_08)
coef(final_model_08)
mse <- mean((predicted - y_test)^2)
mae <- mean(abs(predicted - y_test))
rsquared <- 1 - sum((predicted - y_test)^2) / sum((y_test - mean(y_test))^2)
rsq <- cor(predicted,y_test)^2
cat("MSE_08: ", mse, "\n")
cat("MAE_08: ", mae, "\n")
cat("R-squared_08: ", rsquared, "\n")
cat("R-squared-cor_08: ", rsq, "\n")
nonzero_count_08 <- sum(predicted != 0) - 1 
cat("non zero co 08:", nonzero_count_08)

saveRDS(
  list(
    X_train = X_train,
    y_train = y_train,
    X_test = X_test,
    y_test = y_test,
    train_index = train_index,
    test_index = test_index,
    scaling_attrs = scaling_attrs,
    low_variance_cols  = low_variance_cols,
    enet_train = enet_train,
    final_model_best = final_model_best,
    final_model_09 = final_model_09,
    final_model_08 = final_model_08,
    ml_test_data = ml_test_data,
    target_data = target_data
  ),
  file = "nmr_tuning_objects.rds"
)

save.image(file = "nmr_elastic_tuning_workspace.RData")
