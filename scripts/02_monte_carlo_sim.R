# SCRIPT 02: MONTE CARLO SIMULATION STUDY FOR GGDSF MLEs
# Reproduces Tables 5 and 6

library(numDeriv)

if (requireNamespace("ggdsf", quietly = TRUE)) {
  library(ggdsf)
} else {
  source("R/ggdsf_distribution.R")
}

# 1. SINGLE-SAMPLE MLE & ASYMPTOTIC COVARIANCE ROUTINE

fit_single_mc_sample <- function(x, init) {
  neg_ll_smooth <- function(par) {
    a <- par[[1L]]; b <- par[[2L]]; c <- par[[3L]]; d <- par[[4L]]
    if (a <= 1e-8 || b <= 1e-8 || c <= 1e-8 || d < -1e-8) return(NaN)
    lr <- log(x / a)
    lz <- (b + d * x) * lr
    if (any(lz > 600)) return(NaN)
    z  <- exp(lz)
    g  <- b + d * x + d * x * lr
    if (any(g <= 0)) return(NaN)
    -sum(c * lz - z - lgamma(c) + log(g) - log(x))
  }
  
  neg_ll_const <- function(par) {
    if (par[[2L]] < (par[[4L]] * par[[1L]]) / exp(2)) return(1e12)
    val <- neg_ll_smooth(par)
    if (is.na(val)) return(1e12)
    val
  }
  
  opt <- tryCatch(
    optim(
      par     = init,
      fn      = neg_ll_const,
      method  = "Nelder-Mead",
      control = list(maxit = 15000, reltol = 1e-10, parscale = pmax(init, 1e-3))
    ),
    error = function(e) NULL
  )
  
  if (is.null(opt) || !is.finite(opt$value)) return(NULL)
  est <- opt$par
  
  # Observed Fisher Information Matrix & Asymptotic Standard Errors
  H <- tryCatch(numDeriv::hessian(neg_ll_smooth, est), error = function(e) NULL)
  if (is.null(H) || any(!is.finite(H))) return(NULL)
  
  cov_mat <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(cov_mat)) return(NULL)
  
  diag_v <- diag(cov_mat)
  if (any(!is.finite(diag_v)) || any(diag_v <= 0)) return(NULL)
  
  se <- sqrt(diag_v)
  list(est = est, se = se)
}

# 2. SCENARIO RUNNER FUNCTION

run_mc_scenario <- function(scenario_name, true_par, sample_sizes = c(50, 100, 250, 500), n_reps = 1000, seed = 2026) {
  set.seed(seed)
  cat(sprintf("\n======================================================================\n"))
  cat(sprintf("RUNNING %s: theta = (%.2f, %.2f, %.2f, %.2f), M = %d\n",
              scenario_name, true_par[[1L]], true_par[[2L]], true_par[[3L]], true_par[[4L]], n_reps))
  cat(sprintf("======================================================================\n"))
  
  mc_results_list <- list()
  
  for (n in sample_sizes) {
    cat(sprintf("  -> Sample size n = %d ...\n", n))
    est_mat <- matrix(NA_real_, nrow = n_reps, ncol = 4)
    se_mat  <- matrix(NA_real_, nrow = n_reps, ncol = 4)
    cov_mat <- matrix(NA_real_, nrow = n_reps, ncol = 4)
    
    valid_count  <- 0L
    attempt      <- 0L
    max_attempts <- n_reps * 3L
    
    while (valid_count < n_reps && attempt < max_attempts) {
      attempt <- attempt + 1L
      sample_x <- rggdsf(n, a = true_par[[1L]], b = true_par[[2L]], c = true_par[[3L]], d = true_par[[4L]])
      fit_res  <- fit_single_mc_sample(sample_x, init = true_par)
      
      if (!is.null(fit_res)) {
        valid_count <- valid_count + 1L
        est_mat[valid_count, ] <- fit_res$est
        se_mat[valid_count, ]  <- fit_res$se
        
        ci_low  <- fit_res$est - qnorm(0.975) * fit_res$se
        ci_high <- fit_res$est + qnorm(0.975) * fit_res$se
        cov_mat[valid_count, ] <- as.numeric(true_par >= ci_low & true_par <= ci_high)
      }
    }
    
    est_mat <- est_mat[seq_len(valid_count), , drop = FALSE]
    se_mat  <- se_mat[seq_len(valid_count), , drop = FALSE]
    cov_mat <- cov_mat[seq_len(valid_count), , drop = FALSE]
    
    mean_est <- colMeans(est_mat)
    bias     <- mean_est - true_par
    rmse     <- sqrt(colMeans((sweep(est_mat, 2, true_par, "-"))^2))
    emp_se   <- apply(est_mat, 2, sd)
    asymp_se <- colMeans(se_mat)
    cp_95    <- colMeans(cov_mat) * 100
    
    mc_results_list[[as.character(n)]] <- data.frame(
      n          = c(n, "", "", ""),
      Parameter  = c("a", "b", "c", "d"),
      True_Value = sprintf("%.2f", true_par),
      Mean_MLE   = sprintf("%.4f", mean_est),
      Bias       = sprintf("%.4f", bias),
      RMSE       = sprintf("%.4f", rmse),
      Emp_SE     = sprintf("%.4f", emp_se),
      Asymp_SE   = sprintf("%.4f", asymp_se),
      CP_95_Pct  = sprintf("%.1f", cp_95),
      row.names  = NULL
    )
  }
  
  do.call(rbind, mc_results_list)
}

# 3. EXECUTE SCENARIO I (TABLE 5) AND SCENARIO II (TABLE 6)

sample_grid <- c(50, 100, 250, 500)
M_reps      <- 1000

# Scenario I: Unimodal shape theta = (3.00, 1.20, 1.50, 0.40)
table5_mc <- run_mc_scenario(
  scenario_name = "SCENARIO I (Unimodal Shape)",
  true_par      = c(a = 3.00, b = 1.20, c = 1.50, d = 0.40),
  sample_sizes  = sample_grid,
  n_reps        = M_reps,
  seed          = 2026
)

# Scenario II: Bathtub shape theta = (5.00, 0.45, 0.90, 0.15)
table6_mc <- run_mc_scenario(
  scenario_name = "SCENARIO II (Bathtub Shape)",
  true_par      = c(a = 5.00, b = 0.45, c = 0.90, d = 0.15),
  sample_sizes  = sample_grid,
  n_reps        = M_reps,
  seed          = 2026
)

# 4. PRINT BOTH EXTENDED TABLES

cat("\n=== TABLE 5: MONTE CARLO RESULTS - SCENARIO I (UNIMODAL, M = 1000) ===\n")
print(table5_mc, row.names = FALSE)

cat("\n=== TABLE 6: MONTE CARLO RESULTS - SCENARIO II (BATHTUB, M = 1000) ===\n")
print(table6_mc, row.names = FALSE)
