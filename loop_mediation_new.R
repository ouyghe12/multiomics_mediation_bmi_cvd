library(survival)
library(dplyr)
library(purrr)
library(tibble)
library(broom)
library(tidyr)
library(parallel)
library(glmnet)

get_nonzero_features <- function(glmnet_obj, s = NULL) {
  if (is.null(s)) {
    if (!is.null(glmnet_obj$lambda.min)) {
      s <- glmnet_obj$lambda.min
    } else if (!is.null(glmnet_obj$lambda.1se)) {
      s <- glmnet_obj$lambda.1se
    } else {
      s <- stats::median(glmnet_obj$lambda)
    }
  }
  co <- as.matrix(coef(glmnet_obj, s = s))
  nz <- rownames(co)[abs(co[, 1]) > 0]
  setdiff(nz, "(Intercept)")
}

get_term_stats <- function(fit, term) {
  sx <- summary(fit)
  co <- as.data.frame(sx$coefficients)
  ci <- as.data.frame(sx$conf.int)
  
  if (!term %in% rownames(co)) {
    stop("Cannot find term in model: ", term)
  }
  
  tibble(
    term = term,
    logHR = co[term, "coef"],
    SE = co[term, "se(coef)"],
    z = co[term, "z"],
    p = co[term, "Pr(>|z|)"],
    HR = ci[term, "exp(coef)"],
    HR_L = ci[term, "lower .95"],
    HR_U = ci[term, "upper .95"]
  )
}

cox_mediation_one <- function(data, outcome_time, outcome_status,
                              exposure, mediator, covariates = NULL,
                              scale_mediator = TRUE) {
  
  cols_need <- c(exposure, mediator, outcome_time, outcome_status, covariates)
  cols_need <- unique(cols_need)
  
  d0 <- data %>%
    dplyr::select(dplyr::all_of(cols_need)) %>%
    tidyr::drop_na()
  
  if (nrow(d0) < 20) return(NULL)
  if (sum(d0[[outcome_status]] == 1) < 5) return(NULL)
  
  d <- d0
  if (scale_mediator) {
    d[[mediator]] <- as.numeric(scale(d[[mediator]]))
  }
  
  covar_part <- if (!is.null(covariates) && length(covariates) > 0) {
    paste("+", paste(covariates, collapse = " + "))
  } else {
    ""
  }
  
  f_total <- as.formula(
    paste0("Surv(", outcome_time, ", ", outcome_status, ") ~ ",
           exposure, " ", covar_part)
  )
  
  f_direct <- as.formula(
    paste0("Surv(", outcome_time, ", ", outcome_status, ") ~ ",
           exposure, " + ", mediator, " ", covar_part)
  )
  
  fit_total <- try(
    coxph(f_total, data = d, ties = "efron", x = FALSE, y = FALSE, model = FALSE),
    silent = TRUE
  )
  fit_direct <- try(
    coxph(f_direct, data = d, ties = "efron", x = FALSE, y = FALSE, model = FALSE),
    silent = TRUE
  )
  
  if (inherits(fit_total, "try-error") || inherits(fit_direct, "try-error")) {
    return(NULL)
  }
  
  st_total  <- get_term_stats(fit_total, exposure)
  st_direct <- get_term_stats(fit_direct, exposure)
  
  d_log  <- st_total$logHR
  SEd    <- st_total$SE
  
  c_log  <- st_direct$logHR
  SEc    <- st_direct$SE
  
  ie_log <- d_log - c_log
  SEie <- sqrt(SEd^2 + SEc^2)
  
  z_ie  <- ie_log / SEie
  p_ie  <- 2 * pnorm(-abs(z_ie))
  
  HR_indirect <- exp(ie_log)
  ie_L_log <- ie_log - 1.96 * SEie
  ie_U_log <- ie_log + 1.96 * SEie
  HR_indirect_L <- exp(ie_L_log)
  HR_indirect_U <- exp(ie_U_log)
  
  if (!is.finite(d_log) || abs(d_log) < .Machine$double.eps) {
    MP <- MP_L <- MP_U <- P_MP <- NA_real_
  } else {
    MP <- ie_log / d_log
    
    var_MP <- (c_log^2 / d_log^4) * (SEd^2) + (1 / d_log^2) * (SEc^2)
    
    if (!is.finite(var_MP) || var_MP <= 0) {
      MP_L <- MP_U <- P_MP <- NA_real_
    } else {
      SE_MP <- sqrt(var_MP)
      MP_L <- MP - 1.96 * SE_MP
      MP_U <- MP + 1.96 * SE_MP
      z_mp <- MP / SE_MP
      P_MP <- 2 * (1 - pnorm(abs(z_mp)))
    }
  }
  
  tibble::tibble(
    mediator = mediator,
    n = nrow(d),
    
    HR_total = st_total$HR,
    HR_total_L = st_total$HR_L,
    HR_total_U = st_total$HR_U,
    P_TOTAL = st_total$p,
    
    HR_direct = st_direct$HR,
    HR_direct_L = st_direct$HR_L,
    HR_direct_U = st_direct$HR_U,
    P_DIRECT = st_direct$p,
    
    HR_indirect = HR_indirect,
    HR_indirect_L = HR_indirect_L,
    HR_indirect_U = HR_indirect_U,
    P_INDIRECT = p_ie,
    
    MP = MP,
    MP_L = MP_L,
    MP_U = MP_U,
    P_MP = P_MP
  )
}

run_dataset_mediation <- function(data, dataset_name,
                                  outcomes, exposure, covariates,
                                  mediators,
                                  scale_mediator = TRUE,
                                  n_cores = 1,
                                  verbose = TRUE) {
  
  mediators <- intersect(mediators, colnames(data))
  if (length(mediators) == 0) {
    stop(dataset_name, ": no mediators found in data")
  }
  
  results <- mclapply(seq_along(outcomes), function(i) {
    o <- outcomes[[i]]
    
    res_list <- lapply(mediators, function(med) {
      if (verbose) {
        message("[", dataset_name, "] outcome=", o$name, " | mediator=", med)
      }
      
      out <- try(
        cox_mediation_one(
          data = data,
          outcome_time = o$time,
          outcome_status = o$status,
          exposure = exposure,
          mediator = med,
          covariates = covariates,
          scale_mediator = scale_mediator
        ),
        silent = TRUE
      )
      
      if (inherits(out, "try-error") || is.null(out)) return(NULL)
      
      out$dataset <- dataset_name
      out$disease <- o$name
      out
    })
    
    dplyr::bind_rows(res_list)
  }, mc.cores = n_cores)
  
  res <- dplyr::bind_rows(results)
  
  if (nrow(res) == 0) return(res)
  
  res %>%
    dplyr::group_by(dataset, disease) %>%
    dplyr::mutate(
      p_value = P_INDIRECT,
      p_bh = p.adjust(P_INDIRECT, method = "BH"),
      `Significant/ns` = dplyr::case_when(
        !is.na(p_bh) & p_bh < 0.05 ~ "Significant",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::relocate(
      dataset, disease, mediator, n,
      HR_total, HR_total_L, HR_total_U, P_TOTAL,
      HR_direct, HR_direct_L, HR_direct_U, P_DIRECT,
      HR_indirect, HR_indirect_L, HR_indirect_U, P_INDIRECT,
      MP, MP_L, MP_U, P_MP,
      p_value, p_bh, `Significant/ns`
    )
}


load("nmr_elastic_tuning_workspace.RData")
final_model_nmr <- final_model_09 

target_data_nmr <- readRDS("target_data_nmr.rds")
target_data <- readRDS("target_data.rds")
final_model <- readRDS("elastic_net_model_pro.rds")


outcomes <- list(
  list(name = "T2D",    time = "T2D_followup_time",    status = "T2D_status"),
  list(name = "CAD",    time = "CAD_followup_time",    status = "CAD_status"),
  list(name = "stroke", time = "stroke_followup_time", status = "stroke_status")
)

covariates <- c("f.31.0.0", "f.21003.0.0", "f.20116.0.0",
                "f.21000.0.0", "f.22040.0.0", "f.22189.0.0")

bmi_var <- "f.21001.0.0"
n_cores <- 4
scale_mediator <- TRUE

mediators_met  <- get_nonzero_features(final_model_nmr)
mediators_prot <- get_nonzero_features(final_model)

res_met <- run_dataset_mediation(
  data = target_data_nmr,
  dataset_name = "metabolomics_only",
  outcomes = outcomes,
  exposure = bmi_var,
  covariates = covariates,
  mediators = mediators_met,
  scale_mediator = scale_mediator,
  n_cores = n_cores,
  verbose = TRUE
)

res_prot <- run_dataset_mediation(
  data = target_data,
  dataset_name = "proteomics_only",
  outcomes = outcomes,
  exposure = bmi_var,
  covariates = covariates,
  mediators = mediators_prot,
  scale_mediator = scale_mediator,
  n_cores = n_cores,
  verbose = TRUE
)

final_df <- dplyr::bind_rows(res_met, res_prot)

print(final_df)
write.csv(final_df, "mediation_individual_full_output.csv", row.names = FALSE)
