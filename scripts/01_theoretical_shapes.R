# SCRIPT 01: SECTION 3 - THEORETICAL SHAPES, SIMILARITY, HRF, HRAF & ENTROPIES
# Reproduces Figures 1–5 (300 DPI) with clean, non-overlapping legends

library(ggplot2)
library(patchwork)

if (requireNamespace("ggdsf", quietly = TRUE)) {
  library(ggdsf)
} else {
  source("R/ggdsf_distribution.R")
}

if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)

clean_bottom_legend <- theme(
  legend.position      = "bottom",
  legend.justification = "left",
  legend.title         = element_blank(),
  legend.text          = element_text(size = 8.2),
  legend.key.width     = unit(1.0, "cm"),
  legend.key.height    = unit(0.42, "cm"),
  legend.margin        = margin(t = 2, b = 2)
)

no_legend <- theme(legend.position = "none")

# Figure 1

t_a <- seq(0.02, 4.5, length.out = 800)
specs_pdf_a <- list(
  list(par = c(2.50, 0.70, 0.90, 0.20), lab = "\u03d1 = (2.50, 0.70, 0.90, 0.20) [L-shaped]",       col = "#222222", lty = "solid"),
  list(par = c(2.50, 1.40, 1.50, 0.25), lab = "\u03d1 = (2.50, 1.40, 1.50, 0.25) [Right-skewed]",   col = "#d62728", lty = "dashed"),
  list(par = c(2.50, 2.20, 2.50, 0.30), lab = "\u03d1 = (2.50, 2.20, 2.50, 0.30) [Near-symmetric]", col = "#2ca02c", lty = "dotdash"),
  list(par = c(2.50, 1.20, 4.50, 0.60), lab = "\u03d1 = (2.50, 1.20, 4.50, 0.60) [Sharp peak]",     col = "#1f77b4", lty = "dotted")
)

labs_pdf_a <- vapply(specs_pdf_a, function(x) x$lab, character(1))
cols_pdf_a <- setNames(vapply(specs_pdf_a, function(x) x$col, character(1)), labs_pdf_a)
ltys_pdf_a <- setNames(vapply(specs_pdf_a, function(x) x$lty, character(1)), labs_pdf_a)

df_pdf_a <- do.call(rbind, lapply(specs_pdf_a, function(s) {
  data.frame(
    t     = t_a,
    dens  = dggdsf(t_a, s$par[1], s$par[2], s$par[3], s$par[4]),
    Label = factor(s$lab, levels = labs_pdf_a)
  )
}))

p1_A <- ggplot(df_pdf_a, aes(x = t, y = dens, color = Label, linetype = Label)) +
  geom_line(linewidth = 1.05) +
  scale_color_manual(values = cols_pdf_a) +
  scale_linetype_manual(values = ltys_pdf_a) +
  coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 1.65)) +
  labs(title = "(A) Unimodal & Monotonic PDF Shapes", x = "t", y = "Density f(t)") +
  theme_minimal(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  clean_bottom_legend

specs_pdf_b <- list(
  list(par = c(3.00, 0.85, 1.60, 0.00), lab = "\u03d1 = (3.00, 0.85, 1.60, 0.00) [GGD]",                 col = "#2ca02c", lty = "dotted",  xmax = 16),
  list(par = c(3.00, 0.85, 1.60, 1.20), lab = "\u03d1 = (3.00, 0.85, 1.60, 1.20) [Peak sharpening]",     col = "#9467bd", lty = "dotdash", xmax = 5.0),
  list(par = c(3.00, 0.85, 1.60, 1.85), lab = "\u03d1 = (3.00, 0.85, 1.60, 1.85) [Shoulder / inflection]", col = "#1f77b4", lty = "dashed",  xmax = 4.6),
  list(par = c(3.00, 0.85, 1.60, 2.08), lab = "\u03d1 = (3.00, 0.85, 1.60, 2.08) [Bimodal]",             col = "#d62728", lty = "solid",   xmax = 4.5)
)

labs_pdf_b <- vapply(specs_pdf_b, function(x) x$lab, character(1))
cols_pdf_b <- setNames(vapply(specs_pdf_b, function(x) x$col, character(1)), labs_pdf_b)
ltys_pdf_b <- setNames(vapply(specs_pdf_b, function(x) x$lty, character(1)), labs_pdf_b)

make_sub_pdf_b <- function(idx, show_title = FALSE, keep_legend = FALSE) {
  s  <- specs_pdf_b[[idx]]
  tg <- seq(0.005, s$xmax, length.out = 900)
  df <- data.frame(
    t     = tg,
    dens  = dggdsf(tg, s$par[1], s$par[2], s$par[3], s$par[4]),
    Label = factor(s$lab, levels = labs_pdf_b)
  )
  if (keep_legend) {
    df_dummy <- data.frame(t = NA_real_, dens = NA_real_, Label = factor(labs_pdf_b, levels = labs_pdf_b))
    df <- rbind(df, df_dummy)
  }
  p <- ggplot(df, aes(x = t, y = dens, color = Label, linetype = Label)) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    scale_color_manual(values = cols_pdf_b, drop = FALSE) +
    scale_linetype_manual(values = ltys_pdf_b, drop = FALSE) +
    labs(title = if (show_title) "(B) Unimodal to Bimodal PDF Shapes" else " ", x = "t", y = NULL) +
    theme_minimal(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold", size = 11))
  if (keep_legend) {
    p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) + clean_bottom_legend
  } else {
    p + no_legend
  }
}

p1_B <- (make_sub_pdf_b(1, TRUE, FALSE) + make_sub_pdf_b(2, FALSE, FALSE)) /
  (make_sub_pdf_b(3, FALSE, TRUE) + make_sub_pdf_b(4, FALSE, FALSE))

fig1 <- (p1_A | p1_B) + plot_layout(widths = c(1.0, 1.15), guides = "keep")
ggsave("figures/Figure_1.png", plot = fig1,
       width = 14.0, height = 6.8, units = "in", dpi = 300, bg = "white")

# Figure 2

M_grid <- seq(1.00, 0.50, by = -0.05)
d_vals <- vapply(M_grid, function(m_target) {
  if (m_target == 1) return(0)
  uniroot(function(d) similarity_ggdsf(2.5, 1.5, 2.0, d) - m_target, interval = c(0, 4.4))$root
}, numeric(1))

t_sim <- seq(0.01, 5.2, length.out = 800)
sim_labels <- sprintf("M = %.2f, d = %.3f%s", M_grid, d_vals, ifelse(d_vals == 0, " (GGD)", ""))
df_sim <- do.call(rbind, lapply(seq_along(d_vals), function(i) {
  data.frame(
    t     = t_sim,
    dens  = dggdsf(t_sim, 2.5, 1.5, 2.0, d_vals[i]),
    Curve = factor(sim_labels[i], levels = sim_labels)
  )
}))

fig2 <- ggplot(df_sim, aes(x = t, y = dens, color = Curve, linetype = Curve, linewidth = Curve)) +
  geom_line() +
  scale_color_manual(values = c("#111111", scales::hue_pal()(length(d_vals) - 1))) +
  scale_linetype_manual(values = c("solid", rep(c("dashed", "dotted", "dotdash", "longdash", "twodash"), length.out = length(d_vals) - 1))) +
  scale_linewidth_manual(values = c(1.35, rep(0.85, length(d_vals) - 1))) +
  labs(title = "GGDSF(2.5, 1.5, 2.0, d) deviation from GGD(2.5, 1.5, 2.0)", x = "t", y = "f(t)") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(face = "bold", size = 12),
    legend.position   = c(0.83, 0.62),
    legend.background = element_rect(fill = "white", color = "grey75"),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8.5)
  )

ggsave("figures/Figure_2.png", plot = fig2,
       width = 10.5, height = 6.2, units = "in", dpi = 300, bg = "white")

# Figure 3

t_hrf_a <- seq(0.05, 5.0, length.out = 800)
specs_hrf_a <- list(
  list(par = c(2.50, 1.00, 1.00, 0.00), lab = "\u03d1 = (2.50, 1.00, 1.00, 0.00) [Constant]",     col = "#222222", lty = "dotted"),
  list(par = c(3.00, 0.45, 0.90, 0.00), lab = "\u03d1 = (3.00, 0.45, 0.90, 0.00) [DFR]",          col = "#d62728", lty = "dashed"),
  list(par = c(4.00, 1.35, 1.15, 0.30), lab = "\u03d1 = (4.00, 1.35, 1.15, 0.30) [S-shaped IFR]", col = "#2ca02c", lty = "dotdash"),
  list(par = c(3.00, 2.20, 1.80, 0.35), lab = "\u03d1 = (3.00, 2.20, 1.80, 0.35) [Convex IFR]",   col = "#1f77b4", lty = "solid")
)

labs_hrf_a <- vapply(specs_hrf_a, function(x) x$lab, character(1))
cols_hrf_a <- setNames(vapply(specs_hrf_a, function(x) x$col, character(1)), labs_hrf_a)
ltys_hrf_a <- setNames(vapply(specs_hrf_a, function(x) x$lty, character(1)), labs_hrf_a)

df_hrf_a <- do.call(rbind, lapply(specs_hrf_a, function(s) {
  data.frame(
    t     = t_hrf_a,
    hrf   = hggdsf(t_hrf_a, s$par[1], s$par[2], s$par[3], s$par[4]),
    Label = factor(s$lab, levels = labs_hrf_a)
  )
}))

p3_A <- ggplot(df_hrf_a, aes(x = t, y = hrf, color = Label, linetype = Label)) +
  geom_line(linewidth = 1.05) +
  scale_color_manual(values = cols_hrf_a) +
  scale_linetype_manual(values = ltys_hrf_a) +
  coord_cartesian(xlim = c(0, 5.0), ylim = c(0, 2.3)) +
  labs(title = "(A) Monotonic & Constant HRF Shapes", x = "t", y = "Hazard Rate h(t)") +
  theme_minimal(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  clean_bottom_legend

specs_hrf_b <- list(
  list(par = c(3.00, 0.30, 1.50, 0.40), lab = "\u03d1 = (3.00, 0.30, 1.50, 0.40) [Bathtub]",    col = "#d62728", lty = "dashed",  xmin = 0.03, xmax = 5.4),
  list(par = c(1.00, 0.60, 2.00, 0.00), lab = "\u03d1 = (1.00, 0.60, 2.00, 0.00) [Inverse BT]", col = "#2ca02c", lty = "dotdash", xmin = 0.01, xmax = 6.0),
  list(par = c(6.00, 0.38, 3.60, 0.45), lab = "\u03d1 = (6.00, 0.38, 3.60, 0.45) [W-Shape]",    col = "#1f77b4", lty = "solid",   xmin = 0.005, xmax = 2.3),
  list(par = c(6.00, 2.02, 0.52, 2.45), lab = "\u03d1 = (6.00, 2.02, 0.52, 2.45) [Bimodal]",    col = "#9467bd", lty = "twodash", xmin = 0.005, xmax = 4.2)
)

labs_hrf_b <- vapply(specs_hrf_b, function(x) x$lab, character(1))
cols_hrf_b <- setNames(vapply(specs_hrf_b, function(x) x$col, character(1)), labs_hrf_b)
ltys_hrf_b <- setNames(vapply(specs_hrf_b, function(x) x$lty, character(1)), labs_hrf_b)

make_sub_hrf_b <- function(idx, show_title = FALSE, keep_legend = FALSE) {
  s  <- specs_hrf_b[[idx]]
  tg <- seq(s$xmin, s$xmax, length.out = 900)
  df <- data.frame(
    t     = tg,
    hrf   = hggdsf(tg, s$par[1], s$par[2], s$par[3], s$par[4]),
    Label = factor(s$lab, levels = labs_hrf_b)
  )
  if (keep_legend) {
    df_dummy <- data.frame(t = NA_real_, hrf = NA_real_, Label = factor(labs_hrf_b, levels = labs_hrf_b))
    df <- rbind(df, df_dummy)
  }
  p <- ggplot(df, aes(x = t, y = hrf, color = Label, linetype = Label)) +
    geom_line(linewidth = 1.0, na.rm = TRUE) +
    scale_color_manual(values = cols_hrf_b, drop = FALSE) +
    scale_linetype_manual(values = ltys_hrf_b, drop = FALSE) +
    labs(title = if (show_title) "(B) Non-Monotonic Multi-Phase HRF Shapes" else " ", x = "t", y = NULL) +
    theme_minimal(base_size = 9.5) +
    theme(plot.title = element_text(face = "bold", size = 11))
  if (keep_legend) {
    p + guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) + clean_bottom_legend
  } else {
    p + no_legend
  }
}
p3_B <- (make_sub_hrf_b(1, TRUE, FALSE) + make_sub_hrf_b(2, FALSE, FALSE)) /
  (make_sub_hrf_b(3, FALSE, TRUE) + make_sub_hrf_b(4, FALSE, FALSE))
fig3 <- (p3_A | p3_B) + plot_layout(widths = c(1.0, 1.15), guides = "keep")
ggsave("figures/Figure_3.png", plot = fig3,
       width = 14.0, height = 6.8, units = "in", dpi = 300, bg = "white")

# Figure 4
  
  make_hraf_panel <- function(par, title_str, xmin, xmax) {
    tg   <- seq(xmin, xmax, length.out = 800)
    hrf  <- hggdsf(tg, par[1], par[2], par[3], par[4])
    hraf <- -sggdsf(tg, par[1], par[2], par[3], par[4], log.p = TRUE) / tg
    df   <- rbind(
      data.frame(t = tg, val = hrf,  Type = factor("HRF h(t)",  levels = c("HRF h(t)", "HRAF H(t)"))),
      data.frame(t = tg, val = hraf, Type = factor("HRAF H(t)", levels = c("HRF h(t)", "HRAF H(t)")))
    )
    ggplot(df, aes(x = t, y = val, color = Type, linetype = Type)) +
      geom_line(linewidth = 1.05) +
      scale_color_manual(values = c("HRF h(t)" = "#1f77b4", "HRAF H(t)" = "#d62728")) +
      scale_linetype_manual(values = c("HRF h(t)" = "solid", "HRAF H(t)" = "dashed")) +
      labs(title = title_str, x = "t", y = "Rate Value") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold", size = 10.5), legend.title = element_blank(),
            legend.text = element_text(face = "bold", size = 9.5), legend.key.width = unit(1.3, "cm"))
  }
p4_A <- make_hraf_panel(c(3.00, 0.30, 1.50, 0.40), "(A) \u03d1 = (3.00, 0.30, 1.50, 0.40)", 0.04, 5.4)
p4_B <- make_hraf_panel(c(1.00, 0.60, 2.00, 0.00), "(B) \u03d1 = (1.00, 0.60, 2.00, 0.00)", 0.02, 6.0)
p4_C <- make_hraf_panel(c(6.00, 0.38, 3.60, 0.45), "(C) \u03d1 = (6.00, 0.38, 3.60, 0.45)", 0.005, 2.3)
p4_D <- make_hraf_panel(c(6.00, 2.02, 0.52, 2.45), "(D) \u03d1 = (6.00, 2.02, 0.52, 2.45)", 0.005, 4.2)
fig4 <- ((p4_A + p4_B) / (p4_C + p4_D)) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("figures/Figure_4.png", plot = fig4,
       width = 12.0, height = 7.4, units = "in", dpi = 300, bg = "white")

# Figure 5
  
  specs_ent_a <- list(
    list(par = c(3.00, 0.85, 1.60), lab = "(a, b, c) = (3.00, 0.85, 1.60)", col = "#1f77b4", lty = "solid"),
    list(par = c(2.50, 1.10, 1.50), lab = "(a, b, c) = (2.50, 1.10, 1.50)", col = "#d62728", lty = "dashed"),
    list(par = c(2.50, 1.40, 2.00), lab = "(a, b, c) = (2.50, 1.40, 2.00)", col = "#2ca02c", lty = "dotdash"),
    list(par = c(2.00, 1.20, 2.50), lab = "(a, b, c) = (2.00, 1.20, 2.50)", col = "#9467bd", lty = "dotted")
  )
labs_ent_a <- vapply(specs_ent_a, function(x) x$lab, character(1))
cols_ent_a <- setNames(vapply(specs_ent_a, function(x) x$col, character(1)), labs_ent_a)
ltys_ent_a <- setNames(vapply(specs_ent_a, function(x) x$lty, character(1)), labs_ent_a)
df_ent_a <- do.call(rbind, lapply(specs_ent_a, function(s) {
  dmax <- (s$par[2] * exp(2)) / s$par[1] * 0.99
  dg   <- seq(0, dmax, length.out = 60)
  sh   <- vapply(dg, function(d) entropy_shannon_ggdsf(s$par[1], s$par[2], s$par[3], d), numeric(1))
  data.frame(d = dg, Shannon = sh, Label = factor(s$lab, levels = labs_ent_a))
}))
p5_A <- ggplot(df_ent_a, aes(x = d, y = Shannon, color = Label, linetype = Label)) +
  geom_line(linewidth = 1.05) +
  scale_color_manual(values = cols_ent_a) +
  scale_linetype_manual(values = ltys_ent_a) +
  labs(title = "(A) Shannon Entropy vs. Slope d", x = "Slope Parameter d", y = "Shannon Entropy") +
  theme_minimal(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", size = 11)) +
  guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  clean_bottom_legend
alpha_grid <- setdiff(seq(0.30, 3.00, by = 0.05), 1.0)
d_ent_b    <- c(0.00, 0.60, 1.20, 1.85)
labs_ent_b <- c(
  "\u03d1 = (3.00, 0.85, 1.60, 0.00) [GGD]",
  "\u03d1 = (3.00, 0.85, 1.60, 0.60)",
  "\u03d1 = (3.00, 0.85, 1.60, 1.20)",
  "\u03d1 = (3.00, 0.85, 1.60, 1.85)"
)
cols_ent_b <- setNames(c("#2ca02c", "#9467bd", "#1f77b4", "#d62728"), labs_ent_b)
ltys_ent_b <- setNames(c("dotted", "dotdash", "dashed", "solid"), labs_ent_b)
df_renyi <- do.call(rbind, lapply(seq_along(d_ent_b), function(i) {
  vals <- vapply(alpha_grid, function(al) entropy_renyi_ggdsf(3.00, 0.85, 1.60, d_ent_b[i], alpha = al), numeric(1))
  data.frame(alpha = alpha_grid, val = vals, Label = factor(labs_ent_b[i], levels = labs_ent_b))
}))
df_tsallis <- do.call(rbind, lapply(seq_along(d_ent_b), function(i) {
  vals <- vapply(alpha_grid, function(al) entropy_tsallis_ggdsf(3.00, 0.85, 1.60, d_ent_b[i], q = al), numeric(1))
  data.frame(alpha = alpha_grid, val = vals, Label = factor(labs_ent_b[i], levels = labs_ent_b))
}))
p5_B1 <- ggplot(df_renyi, aes(x = alpha, y = val, color = Label, linetype = Label)) +
  geom_line(linewidth = 1.05) +
  scale_color_manual(values = cols_ent_b) +
  scale_linetype_manual(values = ltys_ent_b) +
  labs(title = "(B) R\u00e9nyi and Tsallis Entropies vs. Order \u03b1", subtitle = "R\u00e9nyi Entropy H_\u03b1(T)",
       x = "Order \u03b1", y = "Entropy Value") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11), plot.subtitle = element_text(face = "bold", hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) +
  clean_bottom_legend
p5_B2 <- ggplot(df_tsallis, aes(x = alpha, y = val, color = Label, linetype = Label)) +
  geom_line(linewidth = 1.05) +
  scale_color_manual(values = cols_ent_b) +
  scale_linetype_manual(values = ltys_ent_b) +
  labs(title = " ", subtitle = "Tsallis Entropy S_\u03b1(T)",
       x = "Order \u03b1", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.subtitle = element_text(face = "bold", hjust = 0.5)) +
  no_legend
p5_B <- p5_B1 | p5_B2
fig5 <- (p5_A | p5_B) + plot_layout(widths = c(0.95, 1.25), guides = "keep")
ggsave("figures/Figure_5.png", plot = fig5,
       width = 14.0, height = 6.2, units = "in", dpi = 300, bg = "white")
