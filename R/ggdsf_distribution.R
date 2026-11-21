# ==============================================================================
# CORE FUNCTIONS FOR THE GENERALIZED GAMMA DISTRIBUTION WITH SHAPE FUNCTION (GGDSF)
# ==============================================================================

#' Check Parametric Admissibility Condition
#'
#' Verifies parameter bounds (\code{a > 0}, \code{b > 0}, \code{c > 0}, \code{d >= 0})
#' and checks the admissibility condition \eqn{b \ge d a e^{-2}}.
#'
#' @param a Positive scale parameter.
#' @param b Positive baseline shape parameter.
#' @param c Positive gamma shape parameter.
#' @param d Non-negative linear shape slope parameter.
#' @return Invisibly returns \code{NULL}; stops on invalid bounds or warns if admissibility is violated.
#' @examples
#' check_ggdsf_params(a = 80, b = 1.5, c = 0.4, d = 0.06)
#' @export
check_ggdsf_params <- function(a, b, c, d) {
  if (any(a <= 0) || any(b <= 0) || any(c <= 0) || any(d < 0)) {
    stop("Parameters must satisfy a > 0, b > 0, c > 0, and d >= 0.")
  }
  if (any(b < (d * a) / exp(2))) {
    warning("Admissibility condition b >= d * a * exp(-2) is violated; PDF may become negative.")
  }
  invisible(NULL)
}

#' Cumulative Distribution Function (CDF) of GGDSF
#'
#' @param q Vector of quantiles.
#' @param a Positive scale parameter.
#' @param b Positive baseline shape parameter.
#' @param c Positive gamma shape parameter.
#' @param d Non-negative linear shape slope parameter.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(T \le q)}, otherwise \eqn{P(T > q)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are given as \code{log(p)}.
#' @return Numeric vector of distribution function values.
#' @importFrom stats pgamma
#' @examples
#' pggdsf(c(5, 20, 50, 80), a = 80, b = 1.5, c = 0.4, d = 0.06)
#' @export
pggdsf <- function(q, a, b, c, d, lower.tail = TRUE, log.p = FALSE) {
  check_ggdsf_params(a, b, c, d)
  res <- rep(if (lower.tail) 0 else 1, length(q))
  valid <- is.finite(q) & (q > 0)
  if (any(valid)) {
    qv <- q[valid]
    lz <- (b + d * qv) * log(qv / a)
    z  <- exp(pmin(lz, 700))
    res[valid] <- pgamma(z, shape = c, rate = 1, lower.tail = lower.tail, log.p = log.p)
  }
  if (any(q == Inf)) {
    val_inf <- if (lower.tail) 1 else 0
    res[q == Inf] <- if (log.p) log(val_inf) else val_inf
  }
  res
}

#' Reliability / Survival Function (RF) of GGDSF
#'
#' @param t Vector of positive time points.
#' @param a Positive scale parameter.
#' @param b Positive baseline shape parameter.
#' @param c Positive gamma shape parameter.
#' @param d Non-negative linear shape slope parameter.
#' @param log.p Logical; if \code{TRUE}, survival probabilities are returned on the log scale.
#' @return Numeric vector of survival probabilities \eqn{R(t) = 1 - F(t)}.
#' @examples
#' sggdsf(c(5, 20, 50, 80), a = 80, b = 1.5, c = 0.4, d = 0.06)
#' @export
sggdsf <- function(t, a, b, c, d, log.p = FALSE) {
  pggdsf(t, a = a, b = b, c = c, d = d, lower.tail = FALSE, log.p = log.p)
}

#' Probability Density Function (PDF) of GGDSF
#'
#' @param x Vector of quantiles.
#' @param a Positive scale parameter.
#' @param b Positive baseline shape parameter.
#' @param c Positive gamma shape parameter.
#' @param d Non-negative linear shape slope parameter.
#' @param log Logical; if \code{TRUE}, densities are returned on the log scale.
#' @return Numeric vector of density values.
#' @examples
#' dggdsf(c(5, 20, 50, 80), a = 80, b = 1.5, c = 0.4, d = 0.06)
#' curve(dggdsf(x, a = 10, b = 0.8, c = 0.9, d = 0.08), from = 0.1, to = 15)
#' @export
dggdsf <- function(x, a, b, c, d, log = FALSE) {
  check_ggdsf_params(a, b, c, d)
  res <- rep(if (log) -Inf else 0, length(x))
  valid <- is.finite(x) & (x > 0)
  if (any(valid)) {
    xv <- x[valid]
    lr <- log(xv / a)
    lz <- (b + d * xv) * lr
    z  <- exp(pmin(lz, 700))
    g  <- b + d * xv + d * xv * lr
    log_dens <- c * lz - z - lgamma(c) + log(pmax(g, .Machine$double.xmin)) - log(xv)
    res[valid] <- if (log) log_dens else exp(log_dens)
  }
  res
}

#' Similarity Measure Relative to Baseline GGD
#'
#' Computes the overlapping coefficient \eqn{M \in (0, 1]} between the GGDSF
#' density and a reference Generalized Gamma density (\code{d = 0}).
#'
#' @param a Positive scale parameter of GGDSF.
#' @param b Positive baseline shape parameter of GGDSF.
#' @param c Positive gamma shape parameter of GGDSF.
#' @param d Non-negative linear shape slope parameter of GGDSF.
#' @param a0 Scale parameter of the reference GGD (defaults to \code{a}).
#' @param b0 Shape parameter of the reference GGD (defaults to \code{b}).
#' @param c0 Gamma shape parameter of the reference GGD (defaults to \code{c}).
#' @return Scalar similarity index \eqn{M}.
#' @importFrom stats integrate
#' @examples
#' # Returns 1 when d = 0 (identical to baseline GGD)
#' similarity_ggdsf(a = 2.5, b = 1.5, c = 2.0, d = 0)
#'
#' # Decreases as slope d increases
#' similarity_ggdsf(a = 2.5, b = 1.5, c = 2.0, d = 0.25)
#' @export
similarity_ggdsf <- function(a, b, c, d, a0 = a, b0 = b, c0 = c) {
  check_ggdsf_params(a, b, c, d)
  check_ggdsf_params(a0, b0, c0, 0)
  integrand <- function(t) {
    f1 <- dggdsf(t, a = a,  b = b,  c = c,  d = d, log = FALSE)
    f0 <- dggdsf(t, a = a0, b = b0, c = c0, d = 0, log = FALSE)
    pmin(f1, f0)
  }
  integrate(integrand, lower = 0, upper = Inf, rel.tol = 1e-10, subdivisions = 500)$value
}

#' Hazard Rate Function (HRF) of GGDSF
#'
#' @param t Vector of positive time points.
#' @param a Positive scale parameter.
#' @param b Positive baseline shape parameter.
#' @param c Positive gamma shape parameter.
#' @param d Non-negative linear shape slope parameter.
#' @return Numeric vector of hazard rate values \eqn{h(t) = f(t) / R(t)}.
#' @examples
#' # Bathtub-shaped hazard trajectory
#' hggdsf(c(1, 10, 40, 75, 85), a = 83.3, b = 1.5, c = 0.37, d = 0.06)
#' @export
hggdsf <- function(t, a, b, c, d) {
  f_val <- dggdsf(t, a = a, b = b, c = c, d = d, log = FALSE)
  r_val <- sggdsf(t, a = a, b = b, c = c, d = d, log.p = FALSE)
  f_val / pmax(r_val, .Machine$double.xmin)
}

#' Reversed Hazard Rate Function (RHRF) of GGDSF
#'
#' @param t Vector of positive time points.
#' @param a Positive scale parameter.
#' @param b Positive baseline shape parameter.
#' @param c Positive gamma shape parameter.
#' @param d Non-negative linear shape slope parameter.
#' @return Numeric vector of reversed hazard rate values \eqn{r(t) = f(t) / F(t)}.
#' @examples
#' rhggdsf(c(5, 20, 50, 80), a = 80, b = 1.5, c = 0.4, d = 0.06)
#' @export
rhggdsf <- function(t, a, b, c, d) {
  f_val <- dggdsf(t, a = a, b = b, c = c, d = d, log = FALSE)
  F_val <- pggdsf(t, a = a, b = b, c = c, d = d, lower.tail = TRUE, log.p = FALSE)
  f_val / pmax(F_val, .Machine$double.xmin)
}

#' Quantile Function of GGDSF
#'
#' @param p Vector of probabilities in \eqn{[0, 1]}.
#' @param a Positive scale parameter.
#' @param b Positive baseline shape parameter.
#' @param c Positive gamma shape parameter.
#' @param d Non-negative linear shape slope parameter.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities correspond to \eqn{P(T \le t)}.
#' @return Numeric vector of quantiles.
#' @importFrom stats qgamma uniroot
#' @examples
#' q_vals <- qggdsf(c(0.1, 0.5, 0.9), a = 10, b = 1.2, c = 1.5, d = 0.08)
#' pggdsf(q_vals, a = 10, b = 1.2, c = 1.5, d = 0.08)
#' @export
qggdsf <- function(p, a, b, c, d, lower.tail = TRUE) {
  check_ggdsf_params(a, b, c, d)
  if (!lower.tail) p <- 1 - p
  vapply(p, function(u) {
    if (is.na(u) || u < 0 || u > 1) return(NaN)
    if (u == 0) return(0)
    if (u == 1) return(Inf)
    upper <- a * max(1.5, qgamma(u, shape = c, rate = 1)^(1 / b))
    while (pggdsf(upper, a, b, c, d) < u) upper <- upper * 2
    uniroot(function(t) pggdsf(t, a, b, c, d) - u, interval = c(0, upper), tol = 1e-11)$root
  }, numeric(1))
}

#' Random Variate Generation for GGDSF
#'
#' @param n Number of observations to generate.
#' @param a Positive scale parameter.
#' @param b Positive baseline shape parameter.
#' @param c Positive gamma shape parameter.
#' @param d Non-negative linear shape slope parameter.
#' @return Numeric vector of \code{n} simulated observations.
#' @importFrom stats runif
#' @examples
#' set.seed(123)
#' x <- rggdsf(20, a = 10, b = 1.2, c = 1.5, d = 0.08)
#' summary(x)
#' @export
rggdsf <- function(n, a, b, c, d) {
  qggdsf(runif(n), a = a, b = b, c = c, d = d)
}

#' Shannon Differential Entropy of GGDSF
#'
#' @param a Positive scale parameter.
#' @param b Positive baseline shape parameter.
#' @param c Positive gamma shape parameter.
#' @param d Non-negative linear shape slope parameter.
#' @return Scalar Shannon differential entropy \eqn{H(T) = -\mathbb{E}[\ln f(T)]}.
#' @importFrom stats integrate
#' @examples
#' entropy_shannon_ggdsf(a = 3.0, b = 0.85, c = 1.6, d = 0.4)
#' @export
entropy_shannon_ggdsf <- function(a, b, c, d) {
  check_ggdsf_params(a, b, c, d)
  integrand <- function(t) {
    f  <- dggdsf(t, a, b, c, d, log = FALSE)
    lf <- dggdsf(t, a, b, c, d, log = TRUE)
    ifelse(f > 0 & is.finite(lf), -f * lf, 0)
  }
  integrate(integrand, lower = 0, upper = Inf, rel.tol = 1e-10, subdivisions = 500)$value
}

#' Renyi Entropy of GGDSF
#'
#' @param a Positive scale parameter.
#' @param b Positive baseline shape parameter.
#' @param c Positive gamma shape parameter.
#' @param d Non-negative linear shape slope parameter.
#' @param alpha Positive entropy order (\code{alpha != 1}).
#' @return Scalar Renyi entropy of order \code{alpha}.
#' @importFrom stats integrate
#' @examples
#' entropy_renyi_ggdsf(a = 3.0, b = 0.85, c = 1.6, d = 0.4, alpha = 2)
#' @export
entropy_renyi_ggdsf <- function(a, b, c, d, alpha = 2) {
  check_ggdsf_params(a, b, c, d)
  if (alpha <= 0 || alpha == 1) stop("Order alpha must satisfy alpha > 0 and alpha != 1.")
  integrand <- function(t) {
    lf <- dggdsf(t, a, b, c, d, log = TRUE)
    ifelse(is.finite(lf), exp(alpha * lf), 0)
  }
  val <- integrate(integrand, lower = 0, upper = Inf, rel.tol = 1e-10, subdivisions = 500)$value
  log(val) / (1 - alpha)
}
#' Tsallis Entropy of GGDSF
#'
#' @param a Positive scale parameter.
#' @param b Positive baseline shape parameter.
#' @param c Positive gamma shape parameter.
#' @param d Non-negative linear shape slope parameter.
#' @param q Positive entropy order (\code{q != 1}).
#' @return Scalar Tsallis entropy of order \code{q}.
#' @importFrom stats integrate
#' @examples
#' entropy_tsallis_ggdsf(a = 3.0, b = 0.85, c = 1.6, d = 0.4, q = 2)
#' @export
entropy_tsallis_ggdsf <- function(a, b, c, d, q = 2) {
  check_ggdsf_params(a, b, c, d)
  if (q <= 0 || q == 1) stop("Order q must satisfy q > 0 and q != 1.")
  integrand <- function(t) {
    lf <- dggdsf(t, a, b, c, d, log = TRUE)
    ifelse(is.finite(lf), exp(q * lf), 0)
  }
  val <- integrate(integrand, lower = 0, upper = Inf, rel.tol = 1e-10, subdivisions = 500)$value
  (1 - val) / (q - 1)
}