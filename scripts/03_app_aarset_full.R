# SCRIPT 03: SECTION 5.1 - COMPLETE AARSET (1987) RELIABILITY DATASET (n = 50)
# Reproduces MLE tables, Observed FIMs, GoF criteria, LRT, and Figure 6 (300 DPI)

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

# 1. FULL UN-CENSORED AARSET (1987) DATASET (n = 50, IN WEEKS)

aarset_data <- c(
  0.1, 0.2, 1, 1, 1, 1, 1, 2, 3, 6,
  7, 11, 12, 18, 18, 18, 18, 18, 21, 32,
  36, 40, 45, 46, 47, 50, 55, 60, 63, 63,
  67, 67, 67, 67, 72, 75, 79, 82, 82, 83,
  84, 84, 84, 85, 85, 85, 85, 85, 86, 86
)
n <- length(aarset_data)

# 2. MLE, OBSERVED FISHER INFORMATION MATRICES & GOODNESS-OF-FIT

# GGDSF
neg_ll4_smooth <- function(par) {
  a <- par[[1L]]; b <- par[[2L]]; c <- par[[3L]]; d <- par[[4L]]
  if (a <= 1e-8 || b <= 1e-8 || c <= 1e-8 || d < -1e-8) return(NaN)
  lr <- log(aarset_data / a); lz <- (b + d * aarset_data) * lr; z <- exp(lz)
  g  <- b + d * aarset_data + d * aarset_data * lr
  if (any(g <= 0) || any(!is.finite(z))) return(NaN)
  -sum(c * lz - z - lgamma(c) + log(g) - log(aarset_data))
}

neg_ll4_const <- function(par) {
  if (par[[2L]] < (par[[4L]] * par[[1L]]) / exp(2)) return(1e12)
  val <- neg_ll4_smooth(par)
  if (is.na(val)) return(1e12)
  val
}

opt4 <- optim(
  par     = c(80.0, 0.40, 0.80, 0.034),
  fn      = neg_ll4_const,
  method  = "L-BFGS-B",
  lower   = c(50.0, 0.15, 0.25, 0.005),
  upper   = c(105.0, 1.50, 2.50, 0.150),
  control = list(factr = 1e-12, maxit = 5000)
)

est4 <- setNames(opt4$par, c("a", "b", "c", "d"))
ll4  <- -opt4$value

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
  lr <- log(aarset_data / a); lz <- b * lr; z <- exp(lz)
  if (any(!is.finite(z))) return(NaN)
  -sum(c * lz - z - lgamma(c) + log(b) - log(aarset_data))
}

opt3 <- optim(
  par     = c(83.0, 3.7, 0.22),
  fn      = neg_ll3_smooth,
  method  = "L-BFGS-B",
  lower   = c(50.0, 0.20, 0.05),
  upper   = c(105.0, 10.0, 5.00),
  control = list(factr = 1e-12, maxit = 5000)
)

est3 <- setNames(opt3$par, c("a", "b", "c"))
ll3  <- -opt3$value

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

# Information Criteria & Goodness-of-Fit Tests

p_ggdsf_fit <- function(t) pggdsf(t, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]])
p_ggd_fit   <- function(t) pggdsf(t, est3[[1L]], est3[[2L]], est3[[3L]], 0)

calc_criteria <- function(ll, k, cdf_fun) {
  aic  <- 2 * k - 2 * ll
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  bic  <- k * log(n) - 2 * ll
  hqic <- 2 * k * log(log(n)) - 2 * ll
  ks   <- suppressWarnings(ks.test(aarset_data, cdf_fun))
  cvm  <- goftest::cvm.test(aarset_data, null = cdf_fun)
  ad   <- goftest::ad.test(aarset_data, null = cdf_fun)
  c(
    sprintf("%.4f", ll), sprintf("%.4f", aic), sprintf("%.4f", aicc),
    sprintf("%.4f", bic), sprintf("%.4f", hqic),
    sprintf("%.4f (%.4f)", ks$statistic,  ks$p.value),
    sprintf("%.4f (%.4f)", cvm$statistic, cvm$p.value),
    sprintf("%.4f (%.4f)", ad$statistic,  ad$p.value)
  )
}

gof_aarset <- data.frame(
  Criterion = c("Log-Likelihood", "AIC", "AICc", "BIC", "HQIC",
                "KS D_n (p-value)", "CvM W^2 (p-value)", "AD A^2 (p-value)"),
  GGDSF = calc_criteria(ll4, k = 4, p_ggdsf_fit),
  GGD   = calc_criteria(ll3, k = 3, p_ggd_fit),
  row.names = NULL
)

lrt_stat <- 2 * (ll4 - ll3)
lrt_pval <- pchisq(max(0, lrt_stat), df = 1, lower.tail = FALSE)


# 3. Figure 6

t_grid     <- seq(0.2, 90, length.out = 1000)
all_levels <- c("Empirical Survival", "Fitted GGDSF", "Fitted GGD")

df_pdf <- rbind(
  data.frame(t = t_grid, dens = dggdsf(t_grid, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]]),
             Model = factor("Fitted GGDSF", levels = all_levels)),
  data.frame(t = t_grid, dens = dggdsf(t_grid, est3[[1L]], est3[[2L]], est3[[3L]], 0),
             Model = factor("Fitted GGD",   levels = all_levels))
)

plot_pdf <- ggplot() +
  geom_histogram(
    data = data.frame(t = aarset_data),
    aes(x = t, y = after_stat(density)),
    breaks = seq(0, 90, by = 10),
    fill = "grey85", color = "grey45", linewidth = 0.5
  ) +
  geom_line(data = df_pdf, aes(x = t, y = dens, color = Model, linetype = Model), linewidth = 1.15) +
  scale_color_manual(values = c("Empirical Survival" = "#222222", "Fitted GGDSF" = "#1f77b4", "Fitted GGD" = "#d62728"), drop = FALSE) +
  scale_linetype_manual(values = c("Empirical Survival" = "solid", "Fitted GGDSF" = "solid", "Fitted GGD" = "dashed"), drop = FALSE) +
  coord_cartesian(xlim = c(0, 90), ylim = c(0, 0.045)) +
  labs(title = "(A) Empirical Histogram vs. Fitted PDFs", x = "Time t (weeks)", y = "Density f(t)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"), legend.key.width = unit(1.3, "cm"))

# Empirical survival step function R_n(t) = 1 - F_n(t)
t_sorted <- sort(aarset_data)
df_emp_surv <- data.frame(
  t     = c(0, t_sorted),
  surv  = c(1, 1 - seq_len(n) / n),
  Model = factor("Empirical Survival", levels = all_levels)
)

df_surv <- rbind(
  data.frame(t = t_grid, surv = sggdsf(t_grid, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]]),
             Model = factor("Fitted GGDSF", levels = all_levels)),
  data.frame(t = t_grid, surv = sggdsf(t_grid, est3[[1L]], est3[[2L]], est3[[3L]], 0),
             Model = factor("Fitted GGD",   levels = all_levels))
)

plot_surv <- ggplot() +
  geom_step(data = df_emp_surv, aes(x = t, y = surv, color = Model, linetype = Model), linewidth = 0.95) +
  geom_line(data = df_surv, aes(x = t, y = surv, color = Model, linetype = Model), linewidth = 1.15) +
  scale_color_manual(values = c("Empirical Survival" = "#222222", "Fitted GGDSF" = "#1f77b4", "Fitted GGD" = "#d62728"), drop = FALSE) +
  scale_linetype_manual(values = c("Empirical Survival" = "solid", "Fitted GGDSF" = "solid", "Fitted GGD" = "dashed"), drop = FALSE) +
  coord_cartesian(xlim = c(0, 90), ylim = c(0, 1.02)) +
  labs(title = "(B) Empirical Survival vs. Fitted Reliability R(t)", x = "Time t (weeks)", y = "Reliability R(t)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"), legend.key.width = unit(1.3, "cm"))

fig_aarset_full <- (plot_pdf | plot_surv) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
print(fig_aarset_full)

ggsave("figures/Figure_6.png", plot = fig_aarset_full,
       width = 13.5, height = 5.6, units = "in", dpi = 300, bg = "white")


# 4. PRINT ALL TABLES IN CONSOLE

cat("\n=== SECTION 5.1: COMBINED MLE TABLE (AARSET FULL, n = 50) ===\n")
print(rbind(df4, df3), row.names = FALSE)

cat("\n=== SECTION 5.1: GOODNESS-OF-FIT & INFORMATION CRITERIA ===\n")
print(gof_aarset, row.names = FALSE)

cat("\n=== SECTION 5.1: OBSERVED FISHER INFORMATION MATRIX (GGDSF) ===\n")
print(round(H4, 4))

cat("\n=== SECTION 5.1: OBSERVED FISHER INFORMATION MATRIX (GGD, d = 0) ===\n")
print(round(H3, 4))

cat("\n=== SECTION 5.1: LIKELIHOOD RATIO TEST (H0: d = 0) ===\n")
print(c(Statistic = round(lrt_stat, 4), df = 1, p_value = round(lrt_pval, 6)))
