################################################################################
# Dependencies: purrr, furrr, tibble, dplyr, WeightIt, MASS
################################################################################

# Load required packages
library(purrr)
library(furrr)
library(tibble)
library(dplyr)

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Calculate Effective Sample Size (ESS)
#'
#' Computes the effective sample size from a vector of weights, defined as
#' ESS = (sum(w))^2 / sum(w^2), which quantifies the information loss due
#' to weighting
#'
#' @param w Numeric vector of weights
#' @return Scalar effective sample size
ESS <- function (w) {
  sum(w)^2 / sum(w^2)
}


################################################################################
# OUTCOME MODELING FUNCTIONS
################################################################################

#' Fit Outcome Regression Models by Treatment Arm
#'
#' Fits separate generalized linear models for treated (A=1) and control (A=0)
#' groups, then predicts outcomes for all observations
#'
#' @param dat Data frame containing outcome, treatment, and covariates
#' @param family GLM family (e.g., "gaussian", "binomial")
#' @param outcome_model Formula specifying the outcome model
#' @return List containing predicted outcomes m1 (treated) and m0 (control)
fit_outcome_model <- function(dat, family, outcome_model) {
  # Fit model on treated subjects and predict for all
  m1 <- glm(outcome_model, family = family, dat %>% filter(A == 1)) %>%
    predict(dat, "response")
  
  # Fit model on control subjects and predict for all
  m0 <- glm(outcome_model, family = family, dat %>% filter(A == 0)) %>%
    predict(dat, "response")
  
  lst(m1, m0)
}


#' AIPW Estimator
#' 
#' @param dat Data frame with outcome Y, treatment A, RCT indicator R
#' @param outcome_model Formula for outcome regression
#' @return Tibble with point estimate (est), standard error (se), and
#'         influence function values (d)
rct_aipw <- function(dat, outcome_model) {
  # Sample sizes
  n_rt <- dat %>% filter(A == 1, R == 1) %>% nrow  # RCT treated
  n_rct <- dat %>% filter(R == 1) %>% nrow          # Total RCT
  dat_rct <- dat %>% filter(R == 1)
  
  # Fit outcome models on RCT data only
  m1 <- glm(outcome_model, family = gaussian, dat_rct %>% filter(A == 1)) %>%
    predict(dat_rct, "response")
  m0 <- glm(outcome_model, family = gaussian, dat_rct %>% filter(A == 0)) %>%
    predict(dat_rct, "response")
  
  # Propensity score (known from randomization)
  pA <- n_rt / n_rct
  
  d <- with(
    dat_rct,
    m1 + A / pA * (Y - m1) - m0 - (1 - A) / (1 - pA) * (Y - m0)
  )
  
  se <- sqrt(sum((d - mean(d))^2) / n_rct^2)
  
  tibble(
    est = mean(d),
    se = se,
    d = list(d)
  )
}


################################################################################
# CONFORMAL SELECTIVE BORROWING FUNCTIONS
################################################################################

#' Conformal Selective Borrow for Single Treatment Arm
#'
#' @param data Full dataset with RCT (R=1) and EC (R=0) observations
#' @param a Treatment arm (0=control, 1=treated)
#' @param gamma Conformity threshold for selection (0-1)
#' @param weight Weighting method: "none" or "residvar" (residual variance)
#' @param qhat Method for computing sampling score ratio: "cw" (calibration),
#'             "pS_est" (estimated), "pS_true" (oracle)
#' @return Tibble with estimate, influence function, number selected, ESS
csb_theata <- function(data, a, gamma, weight, qhat) {
  cv_fold <- 10
  
  # Identify RCT and external control observations
  id_rc <- which(data$R == 1)              # RCT (train & calibration)
  id_test <- which(data$R == 0 & data$A == a)  # EC to test
  y <- data %>% pull(Y)
  x <- data[, 2:3]  # Covariates X1, X2
  
  # Create 10-fold cross-validation splits of RCT data
  split_id_rc <- split(
    id_rc,
    (rep(1:cv_fold, ceiling(length(id_rc) / cv_fold))[1:length(id_rc)]) %>% sample
  )
  
  # Cross-validation: for each fold, train on other folds, test on EC + fold
  compare_all <- map(split_id_rc, function(id_cal) {
    # Split RCT into train and calibration
    id_train <- setdiff(id_rc, id_cal)
    dat_train <- data[id_train, c(2, 3, 6)]
    
    # Predict for EC (test) and calibration set
    id_pred <- c(id_test, id_cal)
    dat_pred <- data[id_pred, c(2, 3, 6)]
    
    # Fit outcome model on training data
    fit <- lm(Y ~ ., data = dat_train)
    yhat_pred <- predict(fit, newdata = dat_pred, type = "response")
    
    # Compute nonconformity scores (absolute residuals)
    score_pred <- abs(yhat_pred - y[id_pred])
    score_test_k <- head(score_pred, length(id_test))   # EC scores
    score_cal_k <- tail(score_pred, length(id_cal))     # Calibration scores
    
    # For each EC, count how many calibration scores are >= EC score
    map_dbl(score_test_k, function(score_test_l) {
      sum(score_cal_k >= score_test_l)
    })
  }) %>% sapply(function(x) x)
  
  # Compute conformal p-values for each EC observation
  # p-value = (# calibration scores >= test score + 1) / (# calibration + 1)
  p_cf <- map_dbl(1:length(id_test), function(l) {
    (sum(compare_all[l, ]) + 1) / (length(id_rc) + 1)
  })
  
  # Select ECs: bias_ec=0 if p_cf >= gamma (exchangeable), 1 otherwise
  bias_ec <- ifelse(p_cf >= gamma, 0, 1)
  
  # Sample size calculations
  n_rt <- data %>% filter(A == 1, R == 1) %>% nrow   # RCT treated
  n_rc <- data %>% filter(A == 0, R == 1) %>% nrow   # RCT control
  n_rct <- data %>% filter(R == 1) %>% nrow          # Total RCT
  n_et <- data %>% filter(A == 1, R == 0) %>% nrow   # EC treated
  n_ec <- data %>% filter(A == 0, R == 0) %>% nrow   # EC control
  
  # Mark selected vs. rejected ECs in full dataset
  data$bias_est <- c(rep(0, n_rct), rep(1, n_et + n_ec))
  data$bias_est[id_test] <- bias_ec
  
  # Combine RCT data with selected ECs (bias_est=0)
  dat <- rbind(data[data$R == 1, ], data[which(data$R == 0 & data$bias_est == 0), ])
  
  n_ra <- dat %>% filter(A == a, R == 1) %>% nrow  # RCT in arm a
  n_e <- dat %>% filter(A == a, R == 0) %>% nrow   # Selected EC in arm a
  
  # Case 1: Too few ECs selected (< 5) - use RCT-only AIPW
  if (sum(bias_ec == 0) < 5) {
    # Remove any ECs from analysis
    if (length(which(dat$R == 0 & dat$A == a)) != 0) {
      dat <- dat[-c(which(dat$R == 0 & dat$A == a)), ]
    }
    
    # Fit outcome model (with or without unmeasured confounder U)
    if (weight == "none") {
      dat$m <- glm(Y ~ X1 + X2, family = gaussian, dat %>% filter(A == a)) %>%
        predict(dat, "response")
    }
    if (weight == "residvar") {
      dat$m <- glm(Y ~ X1 + X2 + U_obs, family = gaussian, dat %>% filter(A == a)) %>%
        predict(dat, "response")
    }
    
    # Known propensity score from randomization
    pA <- n_rt / n_rct
    
    # Standard IPW weights for RCT data
    w <- with(
      dat[which(dat$R == 1), ],
      R * A^a * (1 - A)^(1 - a) / (pA^a * (1 - pA)^(1 - a))
    )
    
    # AIPW influence function
    d <- with(
      dat[which(dat$R == 1), ],
      (R * m + w * (Y - m))
    )
  } else {
    # Case 2: Sufficient ECs selected - implement full borrowing with weighting
    
    # Fit outcome model
    if (weight == "none") {
      # Standard outcome model on observed covariates only
      dat$m <- glm(Y ~ X1 + X2, family = gaussian, dat %>% filter(A == a)) %>%
        predict(dat, "response")
    }
    if (weight == "residvar") {
      # Residual variance-weighted outcome model
      # Fits separate models with/without unmeasured confounder U
      m01 <- glm(Y ~ X1 + X2 + U_obs, family = gaussian, dat %>% filter(A == a, R == 1)) %>%
        predict(dat, "response")
      m0 <- glm(Y ~ X1 + X2, family = gaussian, dat %>% filter(A == a)) %>%
        predict(dat, "response")
      
      # Compute residual variances
      v01 <- glm(Y ~ X1 + X2 + U_obs, family = gaussian, dat %>% filter(A == a, R == 1)) %>%
        resid(type = "response") %>% var
      v0 <- glm(Y ~ X1 + X2, family = gaussian, dat %>% filter(A == a)) %>%
        resid(type = "response") %>% var
      
      # Precision-weighted average of outcome models
      dat$m <- c((v0 * m01[dat$R == 1] + v01 * m0[dat$R == 1]) / (v01 + v0), 
                 m0[dat$R == 0])
    }
    
    # Known propensity score
    pA <- n_rt / n_rct
    
    # Compute residual variance ratio (heteroscedasticity adjustment)
    r1 <- glm(Y ~ X1 + X2 + U_obs, family = gaussian, dat %>% filter(A == a, R == 1)) %>%
      resid(type = "response") %>% var
    r0 <- glm(Y ~ X1 + X2, family = gaussian, dat %>% filter(A == a, R == 0)) %>%
      resid(type = "response") %>% var
    r00 <- r1 / r0
    
    # Compute sampling score ratio qhat (odds of being in RCT vs EC)
    if (qhat == "cw") {
      # Calibration weighting via entropy balancing
      dat$qhat <- compute_cw(dat$R, dat[, 2:3])
    }
    if (qhat == "pS_est") {
      # Estimated via logistic regression
      pS <- glm(R ~ X1 + X2, family = "binomial", dat) %>% 
        predict(dat, "response")
      dat$qhat <- pS / (1 - pS)
    }
    if (qhat == "pS_true") {
      # Oracle (true propensity scores)
      dat$qhat <- dat$pp1 / dat$pp0
    }
    
    # Compute combined weights accounting for RCT/EC sampling and heteroscedasticity
    # Weight formula: qhat * (R*I(A=a) + (1-R)*r00) / (qhat*P(A=a) + r00)
    winit <- with(
      dat,
      qhat * (R * A^a * (1 - A)^(1 - a) + (1 - R) * r00) / 
        (qhat * pA^a * (1 - pA)^(1 - a) + r00)
    )
    
    # Normalize weights to sum to RCT sample size
    w <- winit / sum(winit) * n_rct
    
    # AIPW influence function with EC borrowing
    d <- with(
      dat,
      ((n_rct + n_e) / n_rct) * (R * m + w * (Y - m))
    )
  }
  
  tibble(
    est_theata = mean(d),
    d = list(d),
    n_sel = length(which(dat$R == 0 & dat$A == a)),  # Number of ECs selected
    ess_sel = max(0, ESS(w) - n_ra),                  # Effective sample size gain
    id_sel = list(dat$id)
  )
}


#' Conformal Selective Borrow AIPW for Treatment Effect
#'
#' @param data Full dataset
#' @param gamma0 Conformity threshold for control arm
#' @param gamma1 Conformity threshold for treatment arm
#' @param weight Weighting method
#' @param qhat Sampling score method
#' @return Tibble with ATE estimate, SE, selection counts, and ESS gains
csb_cv_aipw <- function(data, gamma0, gamma1, weight, qhat) {
  
  # Sample sizes
  n_rt <- data %>% filter(A == 1, R == 1) %>% nrow
  n_rc <- data %>% filter(A == 0, R == 1) %>% nrow
  n_rct <- data %>% filter(R == 1) %>% nrow
  n_et <- data %>% filter(A == 1, R == 0) %>% nrow
  n_ec <- data %>% filter(A == 0, R == 0) %>% nrow
  
  # Apply conformal selective borrowing to each arm
  out0 <- csb_theata(data, a = 0, gamma0, weight, qhat)  # Control arm
  out1 <- csb_theata(data, a = 1, gamma1, weight, qhat)  # Treatment arm
  
  # Subset to observations used in either arm
  dat <- data[data$id %in% unique(c(unlist(out0$id_sel), unlist(out1$id_sel))), ]
  
  dat$d1 <- dat$d0 <- 0
  
  # Populate influence functions with appropriate scaling
  dat$d1[which(dat$id %in% unlist(out1$id_sel))] <- 
    nrow(dat) / length(unlist(out1$d)) * unlist(out1$d)
  dat$d0[which(dat$id %in% unlist(out0$id_sel))] <- 
    nrow(dat) / length(unlist(out0$d)) * unlist(out0$d)
  
  # Treatment effect influence function
  d <- dat$d1 - dat$d0
  dof <- nrow(dat)
  
  # Variance estimation accounting for finite RCT sample
  se_i <- d - dat$R / (n_rct / nrow(dat)) * mean(d)
  se <- sqrt(sum(se_i^2) / dof^2)
  
  tibble(
    est = out1$est_theata - out0$est_theata,
    se = se,
    n_sel_t = out1$n_sel,                           # ECs selected in treatment
    n_sel_c = out0$n_sel,                           # ECs selected in control
    n_sel = out1$n_sel + out0$n_sel,                # Total ECs selected
    ess_sel_t = out1$ess_sel,                       # ESS gain in treatment
    ess_sel_c = out0$ess_sel,                       # ESS gain in control
    ess_sel = out1$ess_sel + out0$ess_sel,          # Total ESS gain
    d = list(d),
    id_sel_t = list(which(data$id %in% dat$id[which(dat$R == 0 & dat$A == 1)])),
    id_sel_c = list(which(data$id %in% dat$id[which(dat$R == 0 & dat$A == 0)])),
    id_sel = list(which(data$id %in% dat$id[which(dat$R == 0)]))
  )
}


################################################################################
# CALIBRATION WEIGHTING FUNCTION
################################################################################

#' Compute Calibration Weights via Entropy Balancing
#'
#' @param S Binary indicator (1=RCT, 0=EC)
#' @param X Covariate matrix
#' @return Vector of sampling score ratios (qhat)
compute_cw <- function(S, X) {
  n_rct <- sum(S == 1)
  n_ec <- sum(S == 0)
  
  # Compute entropy balancing weights using WeightIt package
  dat_df <- as.data.frame(X)
  dat_df$S <- S
  W.out <- WeightIt::weightit(
    as.formula(paste("S ~", paste(names(dat_df)[-ncol(dat_df)], collapse = " + "))),
    data = dat_df, estimand = "ATT", method = "ebal", include.obj = TRUE
  )
  w <- W.out$weights / n_ec
  
  # Reproduce calibration weights for all observations
  # (WeightIt constant weights for S=1, need to recover for all)
  Xscale <- WeightIt:::.make_closer_to_1(as.matrix(X))
  Xtarget <- map(1:ncol(Xscale), function(j) {
    Xscale[, j] - colMeans(Xscale[S == 1, , drop = FALSE])[j]
  }) %>% sapply(function(x) x)
  
  # Exponential tilting formula from entropy balancing
  ww_init <- exp(-Xtarget %*% W.out$obj$`0`$par)
  ww <- ww_init / sum(ww_init[S == 0])
  
  # Assign weights to RCT observations
  w[S == 1] <- ww[S == 1]
  qhat <- w * n_rct
  
  # Fallback to logistic regression if entropy balancing fails
  if (any(is.na(qhat) | is.nan(qhat) | is.infinite(qhat))) {
    pS <- glm(S ~ X, family = "binomial", tibble(S, X)) %>%
      predict(tibble(S, X), "response")
    qhat <- pS / (1 - pS)
  }
  
  qhat
}


################################################################################
# ADAPTIVE GAMMA SELECTION
################################################################################

#' Compute Adaptive Gamma via MSE Minimization
#'
#' Selects optimal conformity thresholds (gamma) for treatment and control
#' arms by minimizing estimated mean squared error. Uses bootstrap to estimate
#' variance and bias components across a grid of gamma values
#'
#' @param data Full dataset
#' @param gamma_grid Vector of gamma values to evaluate (default 0 to 1 by 0.1)
#' @param n_rep_gamma Number of bootstrap replications (default 100)
#' @param parallel Use parallel computation
#' @param n_cores Number of cores for parallel processing
#' @param weight Weighting method
#' @param qhat Sampling score method
#' @return List containing selected gamma0 and gamma1
compute_adp_gamma <- function(data,
                              gamma_grid = c(seq(0, 1, by = 0.1)),
                              n_rep_gamma = 100,
                              parallel = FALSE,
                              n_cores = parallel::detectCores(logical = FALSE),
                              weight, qhat,
                              ...) {
  if (!parallel) {n_cores <- 1}
  
  # Ensure gamma=1 (no borrowing) is in grid for bias calculation
  if (any(gamma_grid == 1)) {
    id_nb <- which(gamma_grid == 1)
  } else {
    gamma_grid <- c(gamma_grid, 1)
    id_nb <- which(gamma_grid == 1)
  }
  
  ############################################################################
  # TREATMENT ARM (A=1) GAMMA SELECTION
  ############################################################################
  
  # Prepare data: RCT + EC treatment arm
  dat_rct <- data %>% filter(R == 1)
  dat_ec <- data %>% filter(R == 0, A == 1)
  dat_full <- bind_rows(dat_rct, dat_ec)
  n_rct <- nrow(dat_rct)
  n_full <- nrow(dat_full)
  
  # Generate bootstrap samples for variance estimation
  dat_rep <- map(1:n_rep_gamma, ~ {
    dat_rct_boot <- dat_rct %>%
      group_by(A) %>%
      slice_sample(prop = 1, replace = TRUE)
    bind_rows(dat_rct_boot, dat_ec)
  })
  
  # Evaluate each gamma value
  est_grid_1 <- parallel::mclapply(gamma_grid, function(g) {
    # Estimate on original data
    est_one <- csb_theata(data, 1, g, weight, qhat)$est_theata
    
    # Estimate on bootstrap samples
    est_rep <- map(dat_rep, ~ {
      csb_theata(., 1, g, weight, qhat)$est_theata
    })
    lst(est_one, est_rep)
  }, mc.cores = n_cores)
  
  # Compute MSE = Bias^2 + Variance for each gamma
  res_grid_1 <- map2(est_grid_1, gamma_grid, function(est_g, g) {
    # Variance component from bootstrap
    var_hat <- map_dbl(est_g$est_rep, ~ .) %>% var
    
    if (g == 1) {
      # No borrowing has zero bias by definition
      bias2_hat <- 0
    } else {
      # Bias^2 = (est_gamma - est_no_borrow)^2 - Var(est_gamma - est_no_borrow)
      # Taking max with 0 ensures non-negative bias estimate
      bias2_hat <- max(
        (est_g$est_one - est_grid_1[[id_nb]]$est_one)^2 -
          map2_dbl(est_g$est_rep, est_grid_1[[id_nb]]$est_rep, 
                   function(x, y) {x - y}) %>% var,
        0
      )
    }
    mse_hat <- bias2_hat + var_hat
    lst(mse_hat, bias2_hat, var_hat)
  })
  
  # Select gamma that minimizes MSE
  gamma1_sel <- gamma_grid[which.min(map_dbl(res_grid_1, "mse_hat"))]
  
  ############################################################################
  # CONTROL ARM (A=0) GAMMA SELECTION
  ############################################################################
  
  # Prepare data: RCT + EC control arm
  dat_rct <- data %>% filter(R == 1)
  dat_ec <- data %>% filter(R == 0, A == 0)
  dat_full <- bind_rows(dat_rct, dat_ec)
  n_rct <- nrow(dat_rct)
  n_full <- nrow(dat_full)
  
  # Generate bootstrap samples
  dat_rep <- map(1:n_rep_gamma, ~ {
    dat_rct_boot <- dat_rct %>%
      group_by(A) %>%
      slice_sample(prop = 1, replace = TRUE)
    bind_rows(dat_rct_boot, dat_ec)
  })
  
  # Evaluate each gamma value
  est_grid_0 <- parallel::mclapply(gamma_grid, function(g) {
    est_one <- csb_theata(data, 0, g, weight, qhat)$est_theata
    est_rep <- map(dat_rep, ~ {
      csb_theata(., 0, g, weight, qhat)$est_theata
    })
    lst(est_one, est_rep)
  }, mc.cores = n_cores)
  
  # Compute MSE for each gamma
  res_grid_0 <- map2(est_grid_0, gamma_grid, function(est_g, g) {
    var_hat <- map_dbl(est_g$est_rep, ~ .) %>% var
    if (g == 1) {
      bias2_hat <- 0
    } else {
      bias2_hat <- max(
        (est_g$est_one - est_grid_0[[id_nb]]$est_one)^2 -
          map2_dbl(est_g$est_rep, est_grid_0[[id_nb]]$est_rep, 
                   function(x, y) {x - y}) %>% var,
        0
      )
    }
    mse_hat <- bias2_hat + var_hat
    lst(mse_hat, bias2_hat, var_hat)
  })
  
  # Select gamma that minimizes MSE
  gamma0_sel <- gamma_grid[which.min(map_dbl(res_grid_0, "mse_hat"))]
  
  set <- c(gamma0_sel, gamma1_sel)
  lst(set)
}


################################################################################
# MAIN ESTIMATION FUNCTION
################################################################################

#' External Control Borrowing: Unified Framework
#'
#' Implements multiple methods for incorporating external control data with RCT:
#' - No Borrow: RCT-only analysis (DiM, OM, AIPW)
#' - Borrow Naive: Treat all EC as if from RCT
#' - Borrow: Full borrowing with weighting (IPW, CW, ACW, AIPW, OM)
#' - Conformal Selective Borrow: Selective borrowing via conformal prediction
#'
#' Supports bootstrap and Fisher randomization test inference
#'
#' @param data Data frame with columns: Y (outcome), A (treatment), R (RCT indicator),
#'             X1, X2 (covariates), and optionally U_obs (unmeasured confounder)
#' @param method Character string specifying estimation method
#' @param family GLM family for outcome models
#' @param n_fisher Number of Fisher randomization test permutations (NULL=skip)
#' @param outcome_model Formula for outcome regression
#' @param n_boot Number of bootstrap samples for SE (NULL=skip)
#' @param adp Use adaptive gamma selection
#' @param gamma0 Conformity threshold for control arm
#' @param gamma1 Conformity threshold for treatment arm
#' @param weight Weighting method: "none" or "residvar"
#' @param qhat Sampling score method: "cw", "pS_est", or "pS_true"
#' @param sig_level Significance level for confidence intervals
#' @param parallel Use parallel computation
#' @param n_cores Number of cores for parallel processing
#' @param output_frt Return full Fisher randomization test results
#' @return List containing results table, selected IDs, data info, and estimates
ec_borrow <- function(
    data,
    method = "Conformal Selective Borrow AIPW",
    family = "gaussian",
    n_fisher = NULL,
    outcome_model = paste0("Y~X1+X2"),
    n_boot = NULL,
    adp = FALSE,
    gamma0 = 0.6,
    gamma1 = 0.6,
    weight = "residvar",
    qhat = "pS_true",
    sig_level = 0.05,
    parallel = FALSE,
    n_cores = parallel::detectCores(logical = FALSE),
    output_frt = FALSE
) {
  
  ############################################################################
  # METHOD: NO BORROW - Difference in Means (t-test)
  ############################################################################
  if (identical(method, "No Borrow DiM")) {
    est_fun <- function(dat) {
      fit <- t.test(
        x = dat %>% filter(A == 1, R == 1) %>% pull(Y),
        y = dat %>% filter(A == 0, R == 1) %>% pull(Y)
      )
      tibble(
        est = -(fit$estimate %>% diff %>% unname),
        se = fit$stderr,
        ci_l = fit$conf.int[1],
        ci_u = fit$conf.int[2],
        p_value = fit$p.value,
        ess_sel_t = 0,
        ess_sel_c = 0,
        ess_sel = 0,
        id_sel_t = list(NULL),
        id_sel_c = list(NULL),
        id_sel = list(NULL)
      )
    }
    gamma0 <- 1  # gamma=1 means no EC borrowing
    gamma1 <- 1
  }
  
  ############################################################################
  # METHOD: NO BORROW - Outcome Model only
  ############################################################################
  if (identical(method, "No Borrow OM")) {
    est_fun <- function(dat) {
      dat_naive <- dat[dat$R == 1, ]     
      m10 <- fit_outcome_model(dat_naive, family, outcome_model)
      m1 <- m10$m1
      m0 <- m10$m0
      d <- dat_naive %>%
        mutate(d_i = m1 - m0) %>%
        pull(d_i)
      tibble(
        est = mean(d),
        ess_sel_t = 0,
        ess_sel_c = 0,
        ess_sel = 0,
        id_sel_t = list(NULL),
        id_sel_c = list(NULL),
        id_sel = list(NULL)
      )
    }
    gamma0 <- 1
    gamma1 <- 1
  }
  
  ############################################################################
  # METHOD: NO BORROW - AIPW
  ############################################################################
  if (identical(method, "No Borrow AIPW")) {
    est_fun <- function(dat) {
      rct_aipw(dat, outcome_model) %>%
        mutate(
          ess_sel_t = 0,
          ess_sel_c = 0,
          ess_sel = 0,
          id_sel_t = list(NULL),
          id_sel_c = list(NULL),
          id_sel = list(NULL)
        )
    }
    gamma0 <- 1
    gamma1 <- 1
  }
  
  ############################################################################
  # METHOD: BORROW NAIVE - Treat all EC as RCT
  ############################################################################
  if (identical(method, "Borrow Naive")) {
    est_fun <- function(dat) {
      n_et <- dat %>% filter(A == 1, R == 0) %>% nrow
      n_ec <- dat %>% filter(A == 0, R == 0) %>% nrow
      n_e <- dat %>% filter(R == 0) %>% nrow
      dat_naive <- dat %>% mutate(R = 1)  # Treat all as RCT
      rct_aipw(dat_naive, outcome_model) %>%
        mutate(
          ess_sel_t = n_et,
          ess_sel_c = n_ec,
          ess_sel = n_e,
          id_sel_t = list(which(dat$R == 0 & dat$A == 1)),
          id_sel_c = list(which(dat$R == 0 & dat$A == 0)),
          id_sel = list(which(dat$R == 0))
        )
    }
    gamma0 <- 0  # gamma=0 means borrow all ECs
    gamma1 <- 0
  }
  
  ############################################################################
  # METHOD: BORROW - Difference in Means (t-test)
  ############################################################################
  if (identical(method, "Borrow DiM")) {
    n_et <- dat %>% filter(A == 1, R == 0) %>% nrow
    n_ec <- dat %>% filter(A == 0, R == 0) %>% nrow
    n_e <- dat %>% filter(R == 0) %>% nrow
    dat_naive <- dat %>% mutate(R = 1)
    est_fun <- function(dat) {
      fit <- t.test(
        x = dat_naive %>% filter(A == 1, R == 1) %>% pull(Y),
        y = dat_naive %>% filter(A == 0, R == 1) %>% pull(Y)
      )
      tibble(
        est = -(fit$estimate %>% diff %>% unname),
        se = fit$stderr,
        ci_l = fit$conf.int[1],
        ci_u = fit$conf.int[2],
        p_value = fit$p.value,
        ess_sel_t = n_et,
        ess_sel_c = n_ec,
        ess_sel = n_e,
        id_sel_t = list(which(dat$R == 0 & dat$A == 1)),
        id_sel_c = list(which(dat$R == 0 & dat$A == 0)),
        id_sel = list(which(dat$R == 0))
      )
    }
    gamma0 <- 0
    gamma1 <- 0
  }
  
  ############################################################################
  # METHOD: BORROW - Outcome Model
  ############################################################################
  if (identical(method, "Borrow OM")) {
    est_fun <- function(dat) {
      n_et <- dat %>% filter(A == 1, R == 0) %>% nrow
      n_ec <- dat %>% filter(A == 0, R == 0) %>% nrow
      n_e <- dat %>% filter(R == 0) %>% nrow
      
      # Fit outcome models using all data (RCT + EC)
      m10 <- fit_outcome_model(dat, family, outcome_model)
      m1 <- m10$m1
      m0 <- m10$m0
      
      # Evaluate on RCT data only
      d <- dat %>%
        mutate(d_i = m1 - m0) %>%
        filter(R == 1) %>%
        pull(d_i)
      
      tibble(
        est = mean(d),
        ess_sel_t = n_et,
        ess_sel_c = n_ec,
        ess_sel = n_e,
        id_sel_t = list(which(dat$R == 0 & dat$A == 1)),
        id_sel_c = list(which(dat$R == 0 & dat$A == 0)),
        id_sel = list(which(dat$R == 0))
      )
    }
    gamma0 <- 0
    gamma1 <- 0
  }
  
  ############################################################################
  # METHOD: BORROW - Inverse Probability Weighting (IPW)
  ############################################################################
  if (identical(method, "Borrow IPW")) {
    est_fun <- function(dat) {
      # Sample sizes
      n_rt <- dat %>% filter(A == 1, R == 1) %>% nrow
      n_rc <- dat %>% filter(A == 0, R == 1) %>% nrow
      n_rct <- dat %>% filter(R == 1) %>% nrow
      n_et <- dat %>% filter(A == 1, R == 0) %>% nrow
      n_ec <- dat %>% filter(A == 0, R == 0) %>% nrow
      n_e <- dat %>% filter(R == 0) %>% nrow
      
      # Propensity score for treatment assignment
      pA <- n_rt / n_rct
      
      # Compute sampling score ratio (RCT vs EC)
      if (qhat == "pS_true") {
        dat$qhat <- dat$pp1 / dat$pp0
      } else {
        log_m <- glm(R ~ X1 + X2, dat, family = "binomial")
        dat$qhat <- c(log_m$fitted.values) / ((1 - log_m$fitted.values))
      }
      
      # Compute residual variance ratios for heteroscedasticity adjustment
      # Control arm
      var1 <- lm(outcome_model, dat %>% filter(R == 1 & A == 0)) %>%
        resid(type = "response") %>% var
      var0 <- lm(outcome_model, dat %>% filter(R == 0 & A == 0)) %>%
        resid(type = "response") %>% var
      r0 <- var1 / var0
      
      # Treatment arm
      var1 <- lm(outcome_model, dat %>% filter(R == 1 & A == 1)) %>%
        resid(type = "response") %>% var
      var0 <- lm(outcome_model, dat %>% filter(R == 0 & A == 1)) %>%
        resid(type = "response") %>% var
      r1 <- var1 / var0
      
      # Compute weights for control arm (exclude EC treatment)
      w0init <- with(
        dat[-c(which(dat$R == 0 & dat$A == 1)), ],
        qhat * (R * (1 - A) + (1 - R) * r0) / (qhat * (1 - pA) + r0)
      )
      w0 <- w0init / sum(w0init) * n_rct
      d0 <- with(
        dat[-c(which(dat$R == 0 & dat$A == 1)), ],
        ((n_ec + n_rct) / n_rct) * w0 * Y
      )
      
      # Compute weights for treatment arm (exclude EC control)
      w1init <- with(
        dat[-c(which(dat$R == 0 & dat$A == 0)), ],
        qhat * (R * A + (1 - R) * r1) / (qhat * pA + r1)
      )
      w1 <- w1init / sum(w1init) * n_rct
      d1 <- with(
        dat[-c(which(dat$R == 0 & dat$A == 0)), ],
        ((n_et + n_rct) / n_rct) * w1 * Y
      )
      
      tibble(
        est = mean(d1) - mean(d0),
        ess_sel_t = max(0, ESS(w1) - n_rt),
        ess_sel_c = max(0, ESS(w0) - n_rc),
        ess_sel = max(0, ESS(w1) - n_rt) + max(0, ESS(w0) - n_rc),
        id_sel_t = list(which(dat$R == 0 & dat$A == 1)),
        id_sel_c = list(which(dat$R == 0 & dat$A == 0)),
        id_sel = list(which(dat$R == 0))
      )
    }
    gamma0 <- 0
    gamma1 <- 0
  }
  
  ############################################################################
  # METHOD: BORROW - Calibration Weighting (CW)
  ############################################################################
  if (identical(method, "Borrow CW")) {
    est_fun <- function(dat) {
      # Sample sizes
      n_rt <- dat %>% filter(A == 1, R == 1) %>% nrow
      n_rc <- dat %>% filter(A == 0, R == 1) %>% nrow
      n_rct <- dat %>% filter(R == 1) %>% nrow
      n_et <- dat %>% filter(A == 1, R == 0) %>% nrow
      n_ec <- dat %>% filter(A == 0, R == 0) %>% nrow
      n_e <- dat %>% filter(R == 0) %>% nrow
      
      pA <- n_rt / n_rct
      
      # Compute calibration weights
      dat$qhat <- compute_cw(dat$R, dat[, 2:3])
      
      # Residual variance ratios
      var1 <- lm(outcome_model, dat %>% filter(R == 1 & A == 0)) %>%
        resid(type = "response") %>% var
      var0 <- lm(outcome_model, dat %>% filter(R == 0 & A == 0)) %>%
        resid(type = "response") %>% var
      r0 <- var1 / var0
      
      var1 <- lm(outcome_model, dat %>% filter(R == 1 & A == 1)) %>%
        resid(type = "response") %>% var
      var0 <- lm(outcome_model, dat %>% filter(R == 0 & A == 1)) %>%
        resid(type = "response") %>% var
      r1 <- var1 / var0
      
      # Weights and estimators (same as IPW)
      w0init <- with(
        dat[-c(which(dat$R == 0 & dat$A == 1)), ],
        qhat * (R * (1 - A) + (1 - R) * r0) / (qhat * (1 - pA) + r0)
      )
      w0 <- w0init / sum(w0init) * n_rct
      d0 <- with(
        dat[-c(which(dat$R == 0 & dat$A == 1)), ],
        ((n_ec + n_rct) / n_rct) * w0 * Y
      )
      
      w1init <- with(
        dat[-c(which(dat$R == 0 & dat$A == 0)), ],
        qhat * (R * A + (1 - R) * r1) / (qhat * pA + r1)
      )
      w1 <- w1init / sum(w1init) * n_rct
      d1 <- with(
        dat[-c(which(dat$R == 0 & dat$A == 0)), ],
        ((n_et + n_rct) / n_rct) * w1 * Y
      )
      
      tibble(
        est = mean(d1) - mean(d0),
        ess_sel_t = max(0, ESS(w1) - n_rt),
        ess_sel_c = max(0, ESS(w0) - n_rc),
        ess_sel = max(0, ESS(w1) - n_rt) + max(0, ESS(w0) - n_rc),
        id_sel_t = list(which(dat$R == 0 & dat$A == 1)),
        id_sel_c = list(which(dat$R == 0 & dat$A == 0)),
        id_sel = list(which(dat$R == 0))
      )
    }
    gamma0 <- 0
    gamma1 <- 0
  }
  
  ############################################################################
  # METHOD: BORROW - AIPW (Full borrowing)
  ############################################################################
  if (identical(method, "Borrow AIPW")) {
    est_fun <- function(dat) {
      csb_cv_aipw(dat, gamma0 = 0, gamma1 = 0, weight, qhat) 
    }
    gamma0 <- 0
    gamma1 <- 0
  }
  
  ############################################################################
  # METHOD: BORROW - AIPW with Calibration Weighting
  ############################################################################
  if (identical(method, "Borrow ACW")) {
    est_fun <- function(dat) {
      csb_cv_aipw(dat, gamma0 = 0, gamma1 = 0, weight, qhat = "cw") 
    }
    gamma0 <- 0
    gamma1 <- 0
  }
  
  ############################################################################
  # METHOD: CONFORMAL SELECTIVE BORROW AIPW (Main method)
  ############################################################################
  if (identical(method, "Conformal Selective Borrow AIPW")) {
    est_fun <- function(dat) {
      csb_cv_aipw(dat, gamma0, gamma1, weight, qhat)
    }
  }
  
  
  ############################################################################
  # STEP 1: POINT ESTIMATION
  ############################################################################
  runtime <- system.time(
    out <- est_fun(data)
  )[3] %>% unname
  
  
  ############################################################################
  # STEP 2: STATISTICAL INFERENCE
  ############################################################################
  
  # Methods requiring bootstrap SE
  if (method %in% c("Borrow OM", "No Borrow OM", "Borrow Naive", "Borrow IPW", "Borrow CW")) {
    if (!is.null(n_boot)) {
      # Compute bootstrap standard error
      runtime_boot <- system.time(
        if (parallel) {
          cat(paste0("Parallel computing enabled with ", n_cores,
                     " cores for bootstrap SE of ", method, "\n\n"))
          out_boot <- parallel::mclapply(1:n_boot, function(i) {
            dat_boot <- data %>%
              group_by(A, R) %>%
              slice_sample(prop = 1, replace = TRUE) %>%
              ungroup()
            tryCatch({
              est_fun(dat_boot)$est
            }, error = function(e) {
              NA
            })
          }, mc.cores = n_cores) %>%
            map_dbl(~.)
        } else {
          out_boot <- map_dbl(1:n_boot, ~ {
            dat_boot <- data %>%
              group_by(A, R) %>%
              slice_sample(prop = 1, replace = TRUE) %>%
              ungroup()
            tryCatch({
              est_fun(dat_boot)$est
            }, error = function(e) {
              NA
            })
          })
        }
      )[3] %>% unname
      
      out$se <- sd(out_boot, na.rm = TRUE)
      
      # Warning for NA values in bootstrap
      if (any(is.na(out_boot))) {
        warning(paste0("There are ", sum(is.na(out_boot)),
                       " NA in bootstrap for ", method))
      }
      
      # Construct results with bootstrap SE
      res <- tibble(
        method,
        est = out$est,
        se = out$se,
        ci_l = est - qnorm(1 - sig_level / 2) * se,
        ci_u = est + qnorm(1 - sig_level / 2) * se,
        p_value = (1 - pnorm(abs(est / se))) * 2,
        n_sel_t = map_dbl(out$id_sel_t, length),
        n_sel_c = map_dbl(out$id_sel_c, length),
        n_sel = map_dbl(out$id_sel, length),
        ess_sel_t = out$ess_sel_t,
        ess_sel_c = out$ess_sel_c,
        ess_sel = out$ess_sel,
        runtime = runtime + runtime_boot
      )
    } else {
      # No bootstrap requested
      res <- tibble(
        method,
        est = out$est,
        se = NA,
        ci_l = NA,
        ci_u = NA,
        p_value = NA,
        n_sel_t = map_dbl(out$id_sel_t, length),
        n_sel_c = map_dbl(out$id_sel_c, length),
        n_sel = map_dbl(out$id_sel, length),
        ess_sel_t = out$ess_sel_t,
        ess_sel_c = out$ess_sel_c,
        ess_sel = out$ess_sel,
        runtime = runtime
      )
    }
  } else {
    # Methods with analytic SE (DiM, AIPW, ACW, CSB)
    res <- tibble(
      method,
      est = out$est,
      se = out$se,
      ci_l = est - qnorm(1 - sig_level / 2) * se,
      ci_u = est + qnorm(1 - sig_level / 2) * se,
      p_value = (1 - pnorm(abs(est / se))) * 2,
      n_sel_t = map_dbl(out$id_sel_t, length),
      n_sel_c = map_dbl(out$id_sel_c, length),
      n_sel = map_dbl(out$id_sel, length),
      ess_sel_t = out$ess_sel_t,
      ess_sel_c = out$ess_sel_c,
      ess_sel = out$ess_sel,
      runtime
    )
  }
  
  ############################################################################
  # STEP 3: FISHER RANDOMIZATION TEST
  ############################################################################
  if (!is.null(n_fisher)) {
    runtime_frt <- system.time(
      if (parallel) {
        cat(paste0("Parallel computing enabled with ", n_cores,
                   " cores for ", method, "+FRT\n\n"))
        # Permute treatment assignment within RCT
        out_frt <- parallel::mclapply(1:n_fisher, function(i) {
          dat_rand <- data %>%
            mutate(A = {A[R == 1] <- sample(A[R == 1]); A})
          tryCatch({
            est_fun(dat_rand)
          }, error = function(e) {
            NULL
          })
        }, mc.cores = n_cores) %>%
          map_dfr(~.) %>%
          mutate(
            cond = floor(map_dbl(id_sel, length) / 10) == floor(res$n_sel / 10)
          )
      } else {
        out_frt <- map_dfr(1:n_fisher, ~{
          dat_rand <- data %>%
            mutate(A = {A[R == 1] <- sample(A[R == 1]); A})
          tryCatch({
            est_fun(dat_rand)
          }, error = function(e) {
            NULL
          })
        })
      }
    )[3] %>% unname
    
    # Warnings for missing values
    if (nrow(out_frt) < n_fisher) {
      warning(paste0("nrow(out_frt) is ", nrow(out_frt), " for ",
                     method, "+FRT"))
    }
    if (any(is.na(out_frt$est))) {
      warning(paste0("There are ", sum(is.na(out_frt$est)), " NA in ",
                     method, "+FRT"))
    }
    
    # Fisher randomization test p-value
    res_frt <- tibble(
      method = paste0(method, "+FRT"),
      est = NA,
      se = NA,
      ci_l = NA,
      ci_u = NA,
      p_value = mean(c(abs(out_frt$est) >= abs(out$est), 1), na.rm = TRUE),
      n_sel_t = NA,
      n_sel_c = NA,
      n_sel = NA,
      ess_sel_t = NA,
      ess_sel_c = NA,
      ess_sel = NA,
      runtime = runtime_frt
    )
    
    res <- rbind(res, res_frt)
  } else {
    out_frt <- NULL
  }
  
  # Dataset information
  dat_info <- tibble(
    n_rt = data %>% filter(A == 1, R == 1) %>% nrow,
    n_rc = data %>% filter(A == 0, R == 1) %>% nrow,
    n_rct = n_rt + n_rc,
    n_et = data %>% filter(A == 1, R == 0) %>% nrow,
    n_ec = data %>% filter(A == 0, R == 0) %>% nrow,
    id_ec = list(which(data$R == 0))
  )
  
  # Return results
  if (output_frt) {
    lst(res, id_sel = out$id_sel[[1]], dat_info, gamma0, gamma1, out, out_frt)
  } else {
    lst(res, id_sel = out$id_sel[[1]], dat_info, gamma0, gamma1, out)
  }
}


################################################################################
# DATA GENERATION FUNCTION FOR SIMULATION STUDIES
################################################################################

#' Generate Simulation Data 
#' @param b0 Bias magnitude for control EC subpopulation
#' @param b1 Bias magnitude for treatment EC subpopulation
#' @param coeff Coefficient for unmeasured confounder U
#' @param epsilon Residual SD for EC observations
#' @param corr Use correlated covariates and unmeasured confounder
#' @param rho Proportion of EC with bias
#' @return List containing simulated dataset, true treatment effect, and R-squared values
dat_gen<- function(b0, b1, coeff, epsilon, corr = FALSE, rho = 0.5) {
  
  # Fixed parameters
  pA <- 0.5  # Treatment probability in RCT
  n_r1 <- 600   # RCT sample size
  n_r0 <- 1000  # EC sample size
  N <- 1000000  # Population size for sampling
  
  ############################################################################
  # Generate covariates and unmeasured confounder
  ############################################################################
  if (corr == FALSE) {
    # Independent covariates
    X <- cbind(matrix(rnorm(mean = 0, sd = 1, n = N * 2), N, 2))
    UU <- rnorm(N, 2, 1)
  } else {
    # Correlated covariates and unmeasured confounder
    mu <- c(0, 0, 2)
    Sigma <- matrix(c(1, 0, 0.5,
                      0, 1, 0.5,
                      0.5, 0.5, 1),
                    nrow = 3, byrow = TRUE)
    mvr <- MASS::mvrnorm(n = N, mu = mu, Sigma = Sigma)
    X <- mvr[, 1:2]
    UU <- mvr[, 3]
  }
  
  ############################################################################
  # Generate RCT and EC membership probabilities
  ############################################################################
  eta_1 <- c(0.5, 0.3)
  eta_0 <- c(-0.5, -0.2)
  
  # Selection probabilities (logistic models)
  p1 <- 1 / (1 + exp(X %*% eta_1 + 0.6))
  p0 <- 1 / (1 + exp(X %*% eta_0 - 0.4))
  
  # Sample RCT and EC observations
  R1 <- rbinom(N, 1, p1)
  R0 <- rbinom(N, 1, p0)
  
  # Extract selected observations
  P1 <- p1[which(R1 == 1)[1:n_r1]]
  P0 <- p1[which(R0 == 1)[(length(which(R0 == 1)) - (n_r0 - 1)):length(which(R0 == 1))]]
  pp1 <- c(P1, P0)  # True sampling probabilities (for RCT)
  
  P1 <- p0[which(R1 == 1)[1:n_r1]]
  P0 <- p0[which(R0 == 1)[(length(which(R0 == 1)) - (n_r0 - 1)):length(which(R0 == 1))]]
  pp0 <- c(P1, P0)  # True sampling probabilities (for EC)
  
  ############################################################################
  # Construct dataset
  ############################################################################
  n <- n_r1 + n_r0
  X1 <- X[which(R1 == 1)[1:n_r1], ]
  X0 <- X[which(R0 == 1)[(length(which(R0 == 1)) - (n_r0 - 1)):length(which(R0 == 1))], ]
  X_data <- rbind(X1, X0)
  
  U1 <- UU[which(R1 == 1)[1:n_r1]]
  U0<-UU[which(R0==1)[(length(which(R0==1))-(n_r0-1)):length(which(R0==1))]]
  U<-c(U1,U0)
  
  id<-1:N
  id<-id[c(which(R1==1)[1:n_r1],which(R0==1)[(length(which(R0==1))-(n_r0-1)):length(which(R0==1))])]
  
  A<-c(rep(0:1, each = n_r1 / 2), rep(0:1, each = n_r0 / 2))
  R<-c(rep(1,n_r1),rep(0,n_r0))
  
  beta<-c(2,2,coeff,3,1,1,coeff)
  mat1<-cbind(X_data,U,1,X_data,U)
  mat0<-cbind(X_data,U,0,0,0,0)

  Y1<-rnorm(n,mean=mat0%*%beta,sd=1)
  Y0<-rnorm(n,mean=mat0%*%beta,sd=1)
  
  Y1[R==1]<-Y0[R==1]
  Y1[R==0]<-rnorm(n_r0,mean=mat0[R==0,]%*%beta,sd=epsilon)
  Y0[R==0]<-rnorm(n_r0,mean=mat0[R==0,]%*%beta,sd=epsilon)
  
  samp0<-sample(which(R==0&A==0),round((n_r0/2)*rho),replace = F)
  samp1<-sample(which(R==0&A==1),round((n_r0/2)*rho),replace = F)
  
  Y0[samp0]<-Y0[samp0]-b0
  Y1[samp1]<-Y1[samp1]-b1
  
  samp<-rep(0,n)
  samp[c(samp0,samp1)]<-1
  Y <- A * Y1 + (1 - A) * Y0
  U_obs<-U
  U_obs[R==0]<-NA
  data<-as.data.frame(cbind(id,X_data,Y1,Y0,Y,A,R,samp,U_obs,U,pp1,pp0))
  colnames(data)[2:3]<-c("X1","X2")
  
  true <-0
  
  r2_full<-summary(lm(Y~X1+X2+U+A+A*X1+A*X2+A*U,data))$r.squared
  r2_x<-summary(lm(Y~X1+X2+A+A*X1+A*X2,data))$r.squared
  r2_u<-summary(lm(Y~U+A+A*U,data))$r.squared
  partial_r2_x <- (r2_full - r2_u) / (1 - r2_u)
  partial_r2_u <- (r2_full - r2_x) / (1 - r2_x)
  #partial_r2_x/ partial_r2_u
  
  lst(data, true, b0, b1, coeff, epsilon, corr, rho, partial_r2_x,partial_r2_u)
}

################################################################################
# EXAMPLE RUN
################################################################################
iter=1000
n_cores=40

tie_nb_x<-tie_nb_xu<-tie_fb_x<-tie_fb_xu<-tie_csb_x<-tie_csb_xu<-tie_oracle_x<-tie_oracle_xu<-matrix(rep(NA,2*iter),ncol=2)

for(i in 1:iter){
  print(i)
  data_gen<-dat_gen(b0,b1,coeff,epsilon,corr=F,rho=0.5)
  data<-data_gen$data
  
  tie_nb_x[i,]<-ec_borrow(
    data,
    method = "No Borrow AIPW",
    family = "gaussian",
    n_fisher = 1000,
    # IPW/staIPW/OM/CW
    outcome_model=paste0("Y~X1+X2"),
    n_boot = 1000,
    sig_level = 0.05,
    parallel = T,
    n_cores = n_cores,
    output_frt = T
  )$res$p_value
  
  tie_nb_xu[i,]<-ec_borrow(
    data,
    method = "No Borrow AIPW",
    family = "gaussian",
    n_fisher = 1000,
    # IPW/staIPW/OM/CW
    outcome_model=paste0("Y~X1+X2+U_obs"),
    n_boot = 1000,
    sig_level = 0.05,
    parallel = T,
    n_cores = n_cores,
    output_frt = T
  )$res$p_value
  
  
  tie_fb_x[i,]<-ec_borrow(
    data,
    method = "Borrow AIPW",
    family = "gaussian",
    n_fisher = 1000,
    n_boot = 1000,
    weight="none",
    qhat="pS_true",
    sig_level = 0.05,
    parallel = T,
    n_cores = n_cores,
    output_frt = T
  )$res$p_value
  
  tie_fb_xu[i,]<-ec_borrow(
    data,
    method = "Borrow AIPW",
    family = "gaussian",
    n_fisher = 1000,
    n_boot = 1000,
    weight="residvar",
    qhat="pS_true",
    sig_level = 0.05,
    parallel = T,
    n_cores = n_cores,
    output_frt = T
  )$res$p_value
  
  
  cc<-compute_adp_gamma(data,weight="none",qhat="pS_true")
  gamma_sel<-c(cc$set)
  gamma0_csb=gamma_sel[1]
  gamma1_csb=gamma_sel[2]
  
  tie_csb_x[i,]<-ec_borrow(
    data,
    method = "Conformal Selective Borrow AIPW",
    family = "gaussian",
    n_fisher = 1000,
    n_boot = 1000,
    gamma0 =gamma0_csb,
    gamma1= gamma1_csb,
    weight="none",
    qhat="pS_true",
    sig_level = 0.05,
    parallel = T,
    n_cores = n_cores,
    output_frt = T
  )$res$p_value
  
  cc<-compute_adp_gamma(data,weight="residvar",qhat="pS_true")
  gamma_sel<-c(cc$set)
  gamma0_csb=gamma_sel[1]
  gamma1_csb=gamma_sel[2]
  
  tie_csb_xu[i,]<-ec_borrow(
    data,
    method = "Conformal Selective Borrow AIPW",
    family = "gaussian",
    n_fisher = 1000,
    n_boot = 1000,
    gamma0 =gamma0_csb,
    gamma1= gamma1_csb,
    weight="residvar",
    qhat="pS_true",
    sig_level = 0.05,
    parallel = T,
    n_cores = n_cores,
    output_frt = T
  )$res$p_value
  
  tie_oracle_x<-ec_borrow(
    data[data$samp==0,],
    method = "Borrow AIPW",
    family = "gaussian",
    n_fisher = 1000,
    n_boot = 1000,
    weight="none",
    qhat="pS_true",
    sig_level = 0.05,
    parallel = T,
    n_cores = n_cores,
    output_frt = T
  )$res$p_value
  
  tie_oracle_xu<-ec_borrow(
    data[data$samp==0,],
    method = "Borrow AIPW",
    family = "gaussian",
    n_fisher = 1000,
    n_boot = 1000,
    weight="residvar",
    qhat="pS_true",
    sig_level = 0.05,
    parallel = T,
    n_cores = n_cores,
    output_frt = T
  )$res$p_value
  
  if(i%%50==0){
    res<-cbind(tie_nb_x,tie_nb_xu,tie_fb_x,tie_fb_xu,tie_csb_x,tie_csb_xu,tie_oracle_x,tie_oracle_xu)
  }
  
}

res<-cbind(tie_nb_x,tie_nb_xu,tie_fb_x,tie_fb_xu,tie_csb_x,tie_csb_xu,tie_oracle_x,tie_oracle_xu)




