# SCRIPT 04: SECTION 5.2 - TYPE-I RIGHT-CENSORED AARSET DATASET (C = 75 WEEKS)
# Reproduces MLE tables, Censored FIMs, KM metrics (KS_cens, MSE_KM), LRT & Figure 7

library(numDeriv)
library(survival)
library(ggplot2)
library(patchwork)

if (requireNamespace("ggdsf", quietly = TRUE)) {
  library(ggdsf)
} else {
  source("R/ggdsf_distribution.R")
}


if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)


# 1. RIGHT-CENSORED MLE, OBSERVED FIM & KAPLAN-MEIER GOODNESS-OF-FIT FUNCTION

fit_ggdsf_censored <- function(time, status, init = c(a = 80.0, b = 0.4, c = 0.8, d = 0.035)) {
  if (length(time) != length(status)) stop("time and status must have equal length.")
  if (any(time <= 0)) stop("All observed times must be strictly positive (y > 0).")
  if (!all(status %in% c(0, 1))) stop("status must be 1 (failure) or 0 (right-censored).")
  
  n <- length(time)
  r <- sum(status)
  
  # Smooth negative censored log-likelihood for 4-parameter GGDSF
  neg_loglik_smooth <- function(par) {
    a <- par[[1L]]; b <- par[[2L]]; c <- par[[3L]]; d <- par[[4L]]
    if (a <= 1e-8 || b <= 1e-8 || c <= 1e-8 || d < -1e-8) return(NaN)
    
    lr <- log(time / a)
    lz <- (b + d * time) * lr
    if (any(lz > 600)) return(NaN)
    z  <- exp(lz)
    g  <- b + d * time + d * time * lr
    
    if (any(g[status == 1] <= 0) || any(!is.finite(z))) return(NaN)
    
    ll_fail <- sum(
      c * lz[status == 1] - z[status == 1] - lgamma(c) +
        log(g[status == 1]) - log(time[status == 1])
    )
    ll_cens <- sum(
      pgamma(z[status == 0], shape = c, rate = 1, lower.tail = FALSE, log.p = TRUE)
    )
    
    ll <- ll_fail + ll_cens
    if (!is.finite(ll)) return(NaN)
    -ll
  }
  
  # Constrained wrapper enforcing b >= d*a/e^2 during optimization
  neg_loglik_constrained <- function(par) {
    a <- par[[1L]]; b <- par[[2L]]; c <- par[[3L]]; d <- par[[4L]]
    if (a <= 1e-6 || b <= 1e-6 || c <= 1e-6 || d < 0) return(1e12)
    if (b < (d * a) / exp(2)) return(1e12)
    val <- neg_loglik_smooth(par)
    if (is.na(val)) return(1e12)
    val
  }
  
  # --- Step 1: Fit 4-parameter GGDSF model ---
  opt1 <- optim(
    par = init, fn = neg_loglik_constrained, method = "Nelder-Mead",
    control = list(maxit = 20000, reltol = 1e-12, parscale = pmax(init, 0.05))
  )
  opt2 <- optim(
    par = opt1$par, fn = neg_loglik_constrained, method = "BFGS",
    control = list(maxit = 10000, reltol = 1e-12, parscale = pmax(opt1$par, 0.05))
  )
  
  est4 <- setNames(opt2$par, c("a", "b", "c", "d"))
  ll4  <- -opt2$value
  
  obs_fisher4 <- numDeriv::hessian(func = neg_loglik_smooth, x = est4)
  colnames(obs_fisher4) <- rownames(obs_fisher4) <- c("a", "b", "c", "d")
  se4 <- sqrt(diag(solve(obs_fisher4)))
  
  df4 <- data.frame(
    Model     = c("GGDSF", "", "", ""),
    Parameter = c("a", "b", "c", "d"),
    MLE       = sprintf("%.4f", est4),
    SE        = sprintf("%.4f", se4),
    CI_Lower  = sprintf("%.4f", est4 - qnorm(0.975) * se4),
    CI_Upper  = sprintf("%.4f", est4 + qnorm(0.975) * se4),
    Wald_z    = sprintf("%.2f", est4 / se4),
    row.names = NULL
  )
  
  # --- Step 2: Fit nested 3-parameter GGD sub-model (d = 0) ---
  neg_loglik_ggd_smooth <- function(par3) {
    neg_loglik_smooth(c(par3[[1L]], par3[[2L]], par3[[3L]], 0))
  }
  neg_loglik_ggd_const <- function(par3) {
    neg_loglik_constrained(c(par3[[1L]], par3[[2L]], par3[[3L]], 0))
  }
  
  opt_ggd1 <- optim(
    par = c(90.0, 2.0, 0.4), fn = neg_loglik_ggd_const, method = "Nelder-Mead",
    control = list(maxit = 20000, reltol = 1e-12, parscale = c(90.0, 2.0, 0.4))
  )
  opt_ggd2 <- optim(
    par = opt_ggd1$par, fn = neg_loglik_ggd_const, method = "BFGS",
    control = list(maxit = 10000, reltol = 1e-12, parscale = pmax(opt_ggd1$par, 0.1))
  )
  
  est3 <- setNames(opt_ggd2$par, c("a", "b", "c"))
  ll3  <- -opt_ggd2$value
  
  obs_fisher3 <- numDeriv::hessian(func = neg_loglik_ggd_smooth, x = est3)
  colnames(obs_fisher3) <- rownames(obs_fisher3) <- c("a", "b", "c")
  se3 <- sqrt(diag(solve(obs_fisher3)))
  
  df3 <- data.frame(
    Model     = c("GGD (d = 0)", "", "", ""),
    Parameter = c("a", "b", "c", "d"),
    MLE       = c(sprintf("%.4f", est3), "0 (fixed)"),
    SE        = c(sprintf("%.4f", se3), "—"),
    CI_Lower  = c(sprintf("%.4f", est3 - qnorm(0.975) * se3), "—"),
    CI_Upper  = c(sprintf("%.4f", est3 + qnorm(0.975) * se3), "—"),
    Wald_z    = c(sprintf("%.2f", est3 / se3), "—"),
    row.names = NULL
  )
  
  # --- Step 3: Non-parametric Kaplan-Meier Discrepancies (KS_cens & MSE_KM) ---
  km_fit <- survfit(Surv(time, status) ~ 1)
  t_km   <- km_fit$time[km_fit$n.event > 0]
  R_km   <- km_fit$surv[km_fit$n.event > 0]
  
  R_ggdsf <- sggdsf(t_km, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]])
  R_ggd   <- sggdsf(t_km, est3[[1L]], est3[[2L]], est3[[3L]], 0)
  
  ks_cens4 <- max(abs(R_km - R_ggdsf))
  ks_cens3 <- max(abs(R_km - R_ggd))
  mse_km4  <- mean((R_km - R_ggdsf)^2)
  mse_km3  <- mean((R_km - R_ggd)^2)
  
  # --- Step 4: Information Criteria & LRT for H0: d = 0 ---
  calc_metrics <- function(ll, k, ks_val, mse_val) {
    aic     <- 2 * k - 2 * ll
    aicc    <- aic + (2 * k * (k + 1)) / (n - k - 1)
    bic     <- k * log(n) - 2 * ll
    bic_eff <- k * log(r) - 2 * ll
    hqic    <- 2 * k * log(log(n)) - 2 * ll
    round(c(
      LogLik  = ll,
      AIC     = aic,
      AICc    = aicc,
      BIC     = bic,
      BIC_eff = bic_eff,
      HQIC    = hqic,
      KS_cens = ks_val,
      MSE_KM  = mse_val
    ), 4)
  }
  
  metrics4 <- calc_metrics(ll4, k = 4, ks_cens4, mse_km4)
  metrics3 <- calc_metrics(ll3, k = 3, ks_cens3, mse_km3)
  
  lrt_stat <- 2 * (ll4 - ll3)
  lrt_pval <- pchisq(max(0, lrt_stat), df = 1, lower.tail = FALSE)
  
  list(
    est4                   = est4,
    est3                   = est3,
    km_fit                 = km_fit,
    censoring_summary      = c(n_total = n, failures = r, censored = n - r, pct_censored = round((n - r) / n * 100, 1)),
    combined_summary_table = rbind(df4, df3),
    fisher_matrix_ggdsf    = round(obs_fisher4, 4),
    fisher_matrix_ggd      = round(obs_fisher3, 4),
    fit_comparison         = data.frame(
      Criterion = names(metrics4),
      GGDSF     = metrics4,
      GGD       = metrics3,
      row.names = NULL
    ),
    LRT_d_zero = c(Statistic = round(lrt_stat, 4), df = 1, p_value = round(lrt_pval, 6))
  )
}

# 2. APPLY TYPE-I RIGHT CENSORING AT C = 75 WEEKS TO AARSET DATASET

aarset_raw <- c(
  0.1, 0.2, 1, 1, 1, 1, 1, 2, 3, 6,
  7, 11, 12, 18, 18, 18, 18, 18, 21, 32,
  36, 40, 45, 46, 47, 50, 55, 60, 63, 63,
  67, 67, 67, 67, 72, 75, 79, 82, 82, 83,
  84, 84, 84, 85, 85, 85, 85, 85, 86, 86
)

C_time   <- 75
y_obs    <- pmin(aarset_raw, C_time)
delta    <- as.numeric(aarset_raw <= C_time)
res_cens <- fit_ggdsf_censored(time = y_obs, status = delta, init = c(80, 0.4, 0.8, 0.035))


# 3. Figure 7

est4 <- res_cens$est4
est3 <- res_cens$est3
km   <- res_cens$km_fit

# Compute Kaplan-Meier weighted histogram heights over [0, 75] (width = 15)
breaks_km <- seq(0, C_time, by = 15)
km_times  <- c(0, km$time)
km_survs  <- c(1, km$surv)

R_at <- function(t_val) {
  idx <- max(which(km_times <= t_val))
  km_survs[idx]
}

km_hist_df <- do.call(rbind, lapply(seq_len(length(breaks_km) - 1), function(i) {
  t0 <- breaks_km[[i]]; t1 <- breaks_km[[i + 1L]]
  prob_mass <- R_at(t0) - R_at(t1)
  data.frame(xmin = t0, xmax = t1, ymin = 0, ymax = prob_mass / (t1 - t0))
}))

t_grid     <- seq(0.2, C_time, length.out = 1000)
all_levels <- c("Kaplan–Meier (95% CI)", "Fitted GGDSF", "Fitted GGD")

df_pdf_c <- rbind(
  data.frame(t = t_grid, dens = dggdsf(t_grid, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]]),
             Model = factor("Fitted GGDSF", levels = all_levels)),
  data.frame(t = t_grid, dens = dggdsf(t_grid, est3[[1L]], est3[[2L]], est3[[3L]], 0),
             Model = factor("Fitted GGD",   levels = all_levels))
)

p_cens_pdf <- ggplot() +
  geom_rect(
    data = km_hist_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey85", color = "grey45", linewidth = 0.5
  ) +
  geom_vline(xintercept = C_time, linetype = "dotted", color = "grey30", linewidth = 0.9) +
  geom_line(data = df_pdf_c, aes(x = t, y = dens, color = Model, linetype = Model), linewidth = 1.15) +
  scale_color_manual(values = c("Kaplan–Meier (95% CI)" = "#222222", "Fitted GGDSF" = "#1f77b4", "Fitted GGD" = "#d62728"), drop = FALSE) +
  scale_linetype_manual(values = c("Kaplan–Meier (95% CI)" = "solid", "Fitted GGDSF" = "solid", "Fitted GGD" = "dashed"), drop = FALSE) +
  coord_cartesian(xlim = c(0, 80), ylim = c(0, 0.026)) +
  labs(title = "(A) KM-Weighted Histogram vs. Fitted PDFs", x = "Time t (weeks)", y = "Density f(t)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"), legend.key.width = unit(1.3, "cm"))

df_km <- data.frame(
  t     = c(0, km$time),
  surv  = c(1, km$surv),
  lower = c(1, ifelse(is.na(km$lower), 0, km$lower)),
  upper = c(1, ifelse(is.na(km$upper), 1, km$upper)),
  Model = factor("Kaplan–Meier (95% CI)", levels = all_levels)
)

df_surv_c <- rbind(
  data.frame(t = t_grid, surv = sggdsf(t_grid, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]]),
             Model = factor("Fitted GGDSF", levels = all_levels)),
  data.frame(t = t_grid, surv = sggdsf(t_grid, est3[[1L]], est3[[2L]], est3[[3L]], 0),
             Model = factor("Fitted GGD",   levels = all_levels))
)

p_cens_surv <- ggplot() +
  geom_ribbon(data = df_km, aes(x = t, ymin = lower, ymax = upper), fill = "grey85", alpha = 0.7) +
  geom_step(data = df_km, aes(x = t, y = surv, color = Model, linetype = Model), linewidth = 0.95) +
  geom_vline(xintercept = C_time, linetype = "dotted", color = "grey30", linewidth = 0.9) +
  geom_line(data = df_surv_c, aes(x = t, y = surv, color = Model, linetype = Model), linewidth = 1.15) +
  scale_color_manual(values = c("Kaplan–Meier (95% CI)" = "#222222", "Fitted GGDSF" = "#1f77b4", "Fitted GGD" = "#d62728"), drop = FALSE) +
  scale_linetype_manual(values = c("Kaplan–Meier (95% CI)" = "solid", "Fitted GGDSF" = "solid", "Fitted GGD" = "dashed"), drop = FALSE) +
  coord_cartesian(xlim = c(0, 80), ylim = c(0, 1.02)) +
  labs(title = "(B) Kaplan–Meier Curve vs. Fitted Reliability R(t)", x = "Time t (weeks)", y = "Reliability R(t)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"), legend.key.width = unit(1.3, "cm"))

fig_aarset_cens <- (p_cens_pdf | p_cens_surv) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
print(fig_aarset_cens)

ggsave("figures/Figure_7.png", plot = fig_aarset_cens,
       width = 13.5, height = 5.6, units = "in", dpi = 300, bg = "white")

#  4. PRINT ALL TABLES IN CONSOLE

cat("\n=== SECTION 5.2: CENSORING SUMMARY (C = 75 WEEKS) ===\n")
print(res_cens$censoring_summary)
cat("\n=== SECTION 5.2: COMBINED CENSORED MLE TABLE ===\n")
print(res_cens$combined_summary_table, row.names = FALSE)
cat("\n=== SECTION 5.2: GOODNESS-OF-FIT & KAPLAN-MEIER METRICS ===\n")
print(res_cens$fit_comparison, row.names = FALSE)
cat("\n=== SECTION 5.2: OBSERVED CENSORED FIM (GGDSF) ===\n")
print(res_cens$fisher_matrix_ggdsf)
cat("\n=== SECTION 5.2: OBSERVED CENSORED FIM (GGD, d = 0) ===\n")
print(res_cens$fisher_matrix_ggd)
cat("\n=== SECTION 5.2: LIKELIHOOD RATIO TEST (H0: d = 0) ===\n")
print(res_cens$LRT_d_zero)