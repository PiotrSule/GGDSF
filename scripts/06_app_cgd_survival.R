# SCRIPT 06: SECTION 5.4 - CLINICAL CGD SURVIVAL DATASET (n = 203)
# Reproduces MLE tables, Observed FIMs, GoF criteria, LRT, and Figure 9

library(survival)
library(numDeriv)
library(goftest)
library(ggplot2)
library(patchwork)

if (requireNamespace("ggdsf", quietly = TRUE)) {
  library(ggdsf)
} else {
  source("R/ggdsf_distribution.R")
}

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)


# 1. CLINICAL SURVIVAL DATASET: CGD FOLLOW-UP TIMES (IN DAYS, n = 203)

cgd_data <- as.numeric(survival::cgd$tstop)
n  <- length(cgd_data)
mx <- max(cgd_data)


# 2. MLE, OBSERVED FISHER MATRICES & GOODNESS-OF-FIT

# GGDSF 
neg_ll4_smooth <- function(par) {
  a <- par[[1L]]; b <- par[[2L]]; c <- par[[3L]]; d <- par[[4L]]
  if (a <= 1e-8 || b <= 1e-8 || c <= 1e-8 || d < -1e-8) return(NaN)
  lr <- log(cgd_data / a)
  lz <- (b + d * cgd_data) * lr
  if (any(lz > 600)) return(NaN)
  z  <- exp(lz)
  g  <- b + d * cgd_data + d * cgd_data * lr
  if (any(g <= 0)) return(NaN)
  -sum(c * lz - z - lgamma(c) + log(g) - log(cgd_data))
}

neg_ll4_const <- function(par) {
  if (par[[2L]] < (par[[4L]] * par[[1L]]) / exp(2)) return(1e12)
  val <- neg_ll4_smooth(par)
  if (is.na(val)) return(1e12)
  val
}

starts4 <- list(
  c(quantile(cgd_data, 0.85), 1.0, 1.0, 0.5 / mx),
  c(quantile(cgd_data, 0.95), 0.6, 0.8, 0.8 / mx),
  c(median(cgd_data),         1.5, 1.2, 0.2 / mx)
)

best_opt4 <- NULL
best_ll4  <- -Inf
for (st in starts4) {
  opt <- tryCatch(
    optim(st, neg_ll4_const, method = "Nelder-Mead",
          control = list(maxit = 40000, reltol = 1e-12, parscale = pmax(st, 1e-4))),
    error = function(e) NULL
  )
  if (!is.null(opt) && is.finite(opt$value) && -opt$value > best_ll4) {
    best_ll4  <- -opt$value
    best_opt4 <- opt
  }
}

est4 <- setNames(best_opt4$par, c("a", "b", "c", "d"))
ll4  <- best_ll4

H4  <- numDeriv::hessian(neg_ll4_smooth, est4)
colnames(H4) <- rownames(H4) <- c("a", "b", "c", "d")
se4 <- sqrt(diag(solve(H4)))

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

# GGD (d = 0)
neg_ll3_smooth <- function(par3) {
  a <- par3[[1L]]; b <- par3[[2L]]; c <- par3[[3L]]
  if (a <= 1e-8 || b <= 1e-8 || c <= 1e-8) return(NaN)
  lr <- log(cgd_data / a)
  lz <- b * lr
  if (any(lz > 600)) return(NaN)
  z  <- exp(lz)
  -sum(c * lz - z - lgamma(c) + log(b) - log(cgd_data))
}

starts3 <- list(
  c(median(cgd_data),         1.5, 1.0),
  c(quantile(cgd_data, 0.75), 2.5, 0.5)
)

best_opt3 <- NULL
best_ll3  <- -Inf
for (st3 in starts3) {
  opt3 <- tryCatch(
    optim(st3, neg_ll3_smooth, method = "Nelder-Mead",
          control = list(maxit = 40000, reltol = 1e-12, parscale = pmax(st3, 1e-4))),
    error = function(e) NULL
  )
  if (!is.null(opt3) && is.finite(opt3$value) && -opt3$value > best_ll3) {
    best_ll3  <- -opt3$value
    best_opt3 <- opt3
  }
}

est3 <- setNames(best_opt3$par, c("a", "b", "c"))
ll3  <- best_ll3

H3  <- numDeriv::hessian(neg_ll3_smooth, est3)
colnames(H3) <- rownames(H3) <- c("a", "b", "c")
se3 <- sqrt(diag(solve(H3)))

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

# --- Information Criteria & Goodness-of-Fit Tests ---
p_ggdsf_fit <- function(t) pggdsf(t, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]])
p_ggd_fit   <- function(t) pggdsf(t, est3[[1L]], est3[[2L]], est3[[3L]], 0)

calc_criteria <- function(ll, k, cdf_fun) {
  aic  <- 2 * k - 2 * ll
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  bic  <- k * log(n) - 2 * ll
  hqic <- 2 * k * log(log(n)) - 2 * ll
  ks   <- suppressWarnings(ks.test(cgd_data, cdf_fun))
  cvm  <- goftest::cvm.test(cgd_data, null = cdf_fun)
  ad   <- goftest::ad.test(cgd_data, null = cdf_fun)
  c(
    sprintf("%.4f", ll), sprintf("%.4f", aic), sprintf("%.4f", aicc),
    sprintf("%.4f", bic), sprintf("%.4f", hqic),
    sprintf("%.4f (%.4f)", ks$statistic,  ks$p.value),
    sprintf("%.4f (%.4f)", cvm$statistic, cvm$p.value),
    sprintf("%.4f (%.4f)", ad$statistic,  ad$p.value)
  )
}

gof_cgd <- data.frame(
  Criterion = c("Log-Likelihood", "AIC", "AICc", "BIC", "HQIC",
                "KS D_n (p-value)", "CvM W^2 (p-value)", "AD A^2 (p-value)"),
  GGDSF = calc_criteria(ll4, k = 4, p_ggdsf_fit),
  GGD   = calc_criteria(ll3, k = 3, p_ggd_fit),
  row.names = NULL
)

lrt_stat <- 2 * (ll4 - ll3)
lrt_pval <- pchisq(max(0, lrt_stat), df = 1, lower.tail = FALSE)

# Figure 9

t_grid     <- seq(1.0, mx + 15, length.out = 1000)
all_levels <- c("Empirical Hazard", "Fitted GGDSF", "Fitted GGD")

df_pdf_c <- rbind(
  data.frame(t = t_grid, dens = dggdsf(t_grid, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]]),
             Model = factor("Fitted GGDSF", levels = all_levels)),
  data.frame(t = t_grid, dens = dggdsf(t_grid, est3[[1L]], est3[[2L]], est3[[3L]], 0),
             Model = factor("Fitted GGD",   levels = all_levels))
)

p_cgd_pdf <- ggplot() +
  geom_histogram(
    data = data.frame(t = cgd_data),
    aes(x = t, y = after_stat(density)),
    bins = 18, fill = "grey85", color = "grey45", linewidth = 0.5
  ) +
  geom_line(data = df_pdf_c, aes(x = t, y = dens, color = Model, linetype = Model), linewidth = 1.15) +
  scale_color_manual(values = c("Empirical Hazard" = "#333333", "Fitted GGDSF" = "#1f77b4", "Fitted GGD" = "#d62728"), drop = FALSE) +
  scale_linetype_manual(values = c("Empirical Hazard" = "dotted", "Fitted GGDSF" = "solid", "Fitted GGD" = "dashed"), drop = FALSE) +
  labs(title = "(A) Empirical Histogram vs. Fitted PDFs", x = "Time t (days)", y = "Density f(t)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"), legend.key.width = unit(1.3, "cm"))

breaks_c <- seq(0, ceiling(mx / 40) * 40, by = 40)
emp_haz_c <- do.call(rbind, lapply(seq_len(length(breaks_c) - 1), function(i) {
  t0 <- breaks_c[[i]]; t1 <- breaks_c[[i + 1L]]
  d_ev <- sum(cgd_data > t0 & cgd_data <= t1)
  n_rk <- sum(cgd_data > t0)
  data.frame(t_mid = (t0 + t1) / 2, haz = if (n_rk > 0) d_ev / (n_rk * (t1 - t0)) else NA)
}))

df_emp_hc <- data.frame(t = emp_haz_c$t_mid, haz = emp_haz_c$haz, Model = factor("Empirical Hazard", levels = all_levels))
df_hrf_c  <- rbind(
  data.frame(t = t_grid, haz = hggdsf(t_grid, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]]),
             Model = factor("Fitted GGDSF", levels = all_levels)),
  data.frame(t = t_grid, haz = hggdsf(t_grid, est3[[1L]], est3[[2L]], est3[[3L]], 0),
             Model = factor("Fitted GGD",   levels = all_levels))
)

p_cgd_hrf <- ggplot() +
  geom_point(data = df_emp_hc, aes(x = t, y = haz, color = Model), size = 2.5) +
  geom_line(data = df_emp_hc, aes(x = t, y = haz, color = Model, linetype = Model), linewidth = 0.8) +
  geom_line(data = df_hrf_c, aes(x = t, y = haz, color = Model, linetype = Model), linewidth = 1.15) +
  scale_color_manual(values = c("Empirical Hazard" = "#333333", "Fitted GGDSF" = "#1f77b4", "Fitted GGD" = "#d62728"), drop = FALSE) +
  scale_linetype_manual(values = c("Empirical Hazard" = "dotted", "Fitted GGDSF" = "solid", "Fitted GGD" = "dashed"), drop = FALSE) +
  coord_cartesian(ylim = c(0, 0.035)) +
  labs(title = "(B) Empirical Hazard vs. Fitted HRFs h(t)", x = "Time t (days)", y = "Hazard Rate h(t)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"), legend.key.width = unit(1.3, "cm"))

fig_cgd <- (p_cgd_pdf | p_cgd_hrf) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
print(fig_cgd)

ggsave("figures/Figure_9.png", plot = fig_cgd,
       width = 13.5, height = 5.6, units = "in", dpi = 300, bg = "white")


# 4. PRINT ALL TABLES IN CONSOLE

cat("\n=== SECTION 5.4: COMBINED MLE TABLE (CGD TSTOP, n = 203) ===\n")
print(rbind(df4, df3), row.names = FALSE)

cat("\n=== SECTION 5.4: GOODNESS-OF-FIT & INFORMATION CRITERIA ===\n")
print(gof_cgd, row.names = FALSE)

cat("\n=== SECTION 5.4: OBSERVED FISHER INFORMATION MATRIX (GGDSF) ===\n")
print(round(H4, 4))

cat("\n=== SECTION 5.4: OBSERVED FISHER INFORMATION MATRIX (GGD, d = 0) ===\n")
print(round(H3, 4))

cat("\n=== SECTION 5.4: LIKELIHOOD RATIO TEST (H0: d = 0) ===\n")
print(c(Statistic = round(lrt_stat, 4), df = 1, p_value = signif(lrt_pval, 4)))
