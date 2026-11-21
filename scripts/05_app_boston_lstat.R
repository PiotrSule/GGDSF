# SCRIPT 05: SECTION 5.3 - SOCIOECONOMIC BOSTON LSTAT DATASET (n = 506)
# Reproduces MLE tables, Observed FIMs, GoF criteria, LRT, and Figure 8

library(MASS)
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


# 1. SOCIOECONOMIC DATASET: BOSTON LOWER-STATUS POPULATION SHARE (n = 506)

econ_data <- MASS::Boston$lstat
n <- length(econ_data)


# 2. INTERIOR MLE, OBSERVED FISHER INFORMATION MATRICES & GOODNESS-OF-FIT

# GGDSF
neg_ll4_smooth <- function(par) {
  a <- par[[1L]]; b <- par[[2L]]; c <- par[[3L]]; d <- par[[4L]]
  if (a <= 1e-8 || b <= 1e-8 || c <= 1e-8 || d < -1e-8) return(NaN)
  lr <- log(econ_data / a); lz <- (b + d * econ_data) * lr
  if (any(lz > 600)) return(NaN)
  z  <- exp(lz)
  g  <- b + d * econ_data + d * econ_data * lr
  if (any(g <= 0) || any(!is.finite(z))) return(NaN)
  -sum(c * lz - z - lgamma(c) + log(g) - log(econ_data))
}

neg_ll4_const <- function(par) {
  if (par[[2L]] < (par[[4L]] * par[[1L]]) / exp(2)) return(1e12)
  val <- neg_ll4_smooth(par)
  if (is.na(val)) return(1e12)
  val
}

init4 <- c(quantile(econ_data, 0.85), 1.0, 1.0, 0.5 / max(econ_data))
opt4  <- optim(
  par     = init4,
  fn      = neg_ll4_const,
  method  = "Nelder-Mead",
  control = list(maxit = 60000, reltol = 1e-13, parscale = pmax(init4, 1e-4))
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
  lr <- log(econ_data / a); lz <- b * lr
  if (any(lz > 600)) return(NaN)
  z  <- exp(lz)
  if (any(!is.finite(z))) return(NaN)
  -sum(c * lz - z - lgamma(c) + log(b) - log(econ_data))
}

init3 <- c(median(econ_data), 1.5, 1.0)
opt3  <- optim(
  par     = init3,
  fn      = neg_ll3_smooth,
  method  = "Nelder-Mead",
  control = list(maxit = 60000, reltol = 1e-13, parscale = pmax(init3, 1e-4))
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

# --- Information Criteria & Goodness-of-Fit Tests ---
p_ggdsf_fit <- function(t) pggdsf(t, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]])
p_ggd_fit   <- function(t) pggdsf(t, est3[[1L]], est3[[2L]], est3[[3L]], 0)

calc_criteria <- function(ll, k, cdf_fun) {
  aic  <- 2 * k - 2 * ll
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  bic  <- k * log(n) - 2 * ll
  hqic <- 2 * k * log(log(n)) - 2 * ll
  ks   <- suppressWarnings(ks.test(econ_data, cdf_fun))
  cvm  <- goftest::cvm.test(econ_data, null = cdf_fun)
  ad   <- goftest::ad.test(econ_data, null = cdf_fun)
  c(
    sprintf("%.4f", ll), sprintf("%.4f", aic), sprintf("%.4f", aicc),
    sprintf("%.4f", bic), sprintf("%.4f", hqic),
    sprintf("%.4f (%.4f)", ks$statistic,  ks$p.value),
    sprintf("%.4f (%.4f)", cvm$statistic, cvm$p.value),
    sprintf("%.4f (%.4f)", ad$statistic,  ad$p.value)
  )
}

gof_econ <- data.frame(
  Criterion = c("Log-Likelihood", "AIC", "AICc", "BIC", "HQIC",
                "KS D_n (p-value)", "CvM W^2 (p-value)", "AD A^2 (p-value)"),
  GGDSF = calc_criteria(ll4, k = 4, p_ggdsf_fit),
  GGD   = calc_criteria(ll3, k = 3, p_ggd_fit),
  row.names = NULL
)

lrt_stat <- 2 * (ll4 - ll3)
lrt_pval <- pchisq(max(0, lrt_stat), df = 1, lower.tail = FALSE)

# Figure 8

t_grid     <- seq(1.5, max(econ_data) + 1.0, length.out = 1000)
all_levels <- c("Empirical CDF", "Fitted GGDSF", "Fitted GGD")

df_pdf_e <- rbind(
  data.frame(t = t_grid, dens = dggdsf(t_grid, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]]),
             Model = factor("Fitted GGDSF", levels = all_levels)),
  data.frame(t = t_grid, dens = dggdsf(t_grid, est3[[1L]], est3[[2L]], est3[[3L]], 0),
             Model = factor("Fitted GGD",   levels = all_levels))
)

p_econ_pdf <- ggplot() +
  geom_histogram(
    data = data.frame(t = econ_data),
    aes(x = t, y = after_stat(density)),
    bins = 22, fill = "grey85", color = "grey45", linewidth = 0.5
  ) +
  geom_line(data = df_pdf_e, aes(x = t, y = dens, color = Model, linetype = Model), linewidth = 1.15) +
  scale_color_manual(values = c("Empirical CDF" = "#333333", "Fitted GGDSF" = "#1f77b4", "Fitted GGD" = "#d62728"), drop = FALSE) +
  scale_linetype_manual(values = c("Empirical CDF" = "dotted", "Fitted GGDSF" = "solid", "Fitted GGD" = "dashed"), drop = FALSE) +
  labs(title = "(A) Empirical Histogram vs. Fitted PDFs", x = "Lower-Status Population Share (%)", y = "Density f(x)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"), legend.key.width = unit(1.3, "cm"))

x_sorted <- sort(econ_data)
df_ecdf  <- data.frame(t = x_sorted, cdf = seq_len(n) / n, Model = factor("Empirical CDF", levels = all_levels))
df_cdf_e <- rbind(
  data.frame(t = t_grid, cdf = pggdsf(t_grid, est4[[1L]], est4[[2L]], est4[[3L]], est4[[4L]]),
             Model = factor("Fitted GGDSF", levels = all_levels)),
  data.frame(t = t_grid, cdf = pggdsf(t_grid, est3[[1L]], est3[[2L]], est3[[3L]], 0),
             Model = factor("Fitted GGD",   levels = all_levels))
)

p_econ_cdf <- ggplot() +
  geom_step(data = df_ecdf, aes(x = t, y = cdf, color = Model, linetype = Model), linewidth = 1.0) +
  geom_line(data = df_cdf_e, aes(x = t, y = cdf, color = Model, linetype = Model), linewidth = 1.15) +
  scale_color_manual(values = c("Empirical CDF" = "#333333", "Fitted GGDSF" = "#1f77b4", "Fitted GGD" = "#d62728"), drop = FALSE) +
  scale_linetype_manual(values = c("Empirical CDF" = "dotted", "Fitted GGDSF" = "solid", "Fitted GGD" = "dashed"), drop = FALSE) +
  labs(title = "(B) Empirical CDF vs. Fitted CDFs", x = "Lower-Status Population Share (%)", y = "Cumulative Probability F(x)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11.5), legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"), legend.key.width = unit(1.3, "cm"))

fig_econ <- (p_econ_pdf | p_econ_cdf) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
print(fig_econ)

ggsave("figures/Figure_8.png", plot = fig_econ,
       width = 13.5, height = 5.6, units = "in", dpi = 300, bg = "white")


# 4. PRINT ALL TABLES IN CONSOLE

cat("\n=== SECTION 5.3: COMBINED MLE TABLE (BOSTON LSTAT, n = 506) ===\n")
print(rbind(df4, df3), row.names = FALSE)

cat("\n=== SECTION 5.3: GOODNESS-OF-FIT & INFORMATION CRITERIA ===\n")
print(gof_econ, row.names = FALSE)

cat("\n=== SECTION 5.3: OBSERVED FISHER INFORMATION MATRIX (GGDSF) ===\n")
print(round(H4, 4))

cat("\n=== SECTION 5.3: OBSERVED FISHER INFORMATION MATRIX (GGD, d = 0) ===\n")
print(round(H3, 4))

cat("\n=== SECTION 5.3: LIKELIHOOD RATIO TEST (H0: d = 0) ===\n")
print(c(Statistic = round(lrt_stat, 4), df = 1, p_value = round(lrt_pval, 6)))
