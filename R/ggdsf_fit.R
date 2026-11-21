# ==============================================================================
# MAXIMUM LIKELIHOOD ESTIMATION FOR COMPLETE AND RIGHT-CENSORED GGDSF SAMPLES
# ==============================================================================

#' Maximum Likelihood Estimation for Complete GGDSF Samples
#'
#' Fits the four-parameter GGDSF distribution and the nested three-parameter
#' Generalized Gamma sub-model (\code{d = 0}) via maximum likelihood, returning
#' parameter estimates, asymptotic standard errors from the observed Fisher
#' information matrix, information criteria, and a likelihood-ratio test.
#'
#' @param x Numeric vector of strictly positive observations.
#' @param init Optional numeric vector of initial parameter values \code{c(a, b, c, d)}.
#' @return A list containing:
#' \describe{
#'   \item{\code{estimates}}{Named vector of GGDSF MLEs \code{(a, b, c, d)}.}
#'   \item{\code{ggd_estimates}}{Named vector of nested GGD MLEs \code{(a, b, c)}.}
#'   \item{\code{summary_table}}{Data frame with MLEs, SEs, 95\% Wald CIs, and Wald z-statistics.}
#'   \item{\code{fisher_matrix}}{Observed Fisher information matrix for GGDSF.}
#'   \item{\code{criteria}}{Data frame comparing Log-Likelihood, AIC, AICc, BIC, and HQIC.}
#'   \item{\code{lrt}}{Named vector with LRT statistic, degrees of freedom, and p-value for \eqn{H_0: d = 0}.}
#' }
#' @importFrom stats optim qnorm pchisq median quantile
#' @importFrom numDeriv hessian
#' @examples
#' data(aarset)
#' res <- fit_ggdsf(aarset, init = c(80, 1.5, 0.4, 0.06))
#' print(res$summary_table)
#' print(res$lrt)
#' @export
fit_ggdsf <- function(x, init = NULL) {
  x <- as.numeric(x[is.finite(x)])
  if (length(x) < 5 || any(x <= 0)) {
    stop("Input 'x' must contain at least 5 strictly positive finite observations.")
  }
  n  <- length(x)
  mx <- max(x)
  
  neg_ll4_smooth <- function(par) {
    a <- par[[1L]]; b <- par[[2L]]; c <- par[[3L]]; d <- par[[4L]]
    if (a <= 1e-8 || b <= 1e-8 || c <= 1e-8 || d < -1e-8) return(NaN)
    lr <- log(x / a); lz <- (b + d * x) * lr
    if (any(lz > 600)) return(NaN)
    z  <- exp(lz)
    g  <- b + d * x + d * x * lr
    if (any(g <= 0) || any(!is.finite(z))) return(NaN)
    -sum(c * lz - z - lgamma(c) + log(g) - log(x))
  }
  
  neg_ll4_const <- function(par) {
    if (par[[2L]] < (par[[4L]] * par[[1L]]) / exp(2)) return(1e12)
    val <- neg_ll4_smooth(par)
    if (is.na(val)) 1e12 else val
  }
  
  starts4 <- if (!is.null(init)) {
    list(as.numeric(init))
  } else {
    list(
      c(quantile(x, 0.85), 1.0, 1.0, 0.5 / mx),
      c(quantile(x, 0.90), 0.5, 0.8, 0.8 / mx),
      c(median(x),         1.5, 1.2, 0.2 / mx)
    )
  }
  
  best_opt4 <- NULL; best_ll4 <- -Inf
  for (st in starts4) {
    opt <- tryCatch(
      optim(st, neg_ll4_const, method = "Nelder-Mead",
            control = list(maxit = 40000, reltol = 1e-12, parscale = pmax(st, 1e-4))),
      error = function(e) NULL
    )
    if (!is.null(opt) && is.finite(opt$value) && -opt$value > best_ll4) {
      best_ll4 <- -opt$value; best_opt4 <- opt
    }
  }
  if (is.null(best_opt4)) stop("Optimization failed to converge for GGDSF.")
  
  est4 <- setNames(best_opt4$par, c("a", "b", "c", "d"))
  ll4  <- best_ll4
  
  H4  <- numDeriv::hessian(neg_ll4_smooth, est4)
  colnames(H4) <- rownames(H4) <- c("a", "b", "c", "d")
  se4 <- tryCatch(sqrt(diag(solve(H4))), error = function(e) rep(NA_real_, 4))
  
  df4 <- data.frame(
    Model     = c("GGDSF", "", "", ""),
    Parameter = c("a", "b", "c", "d"),
    MLE       = round(est4, 4),
    SE        = round(se4, 4),
    CI_Lower  = round(est4 - qnorm(0.975) * se4, 4),
    CI_Upper  = round(est4 + qnorm(0.975) * se4, 4),
    Wald_z    = round(est4 / se4, 2),
    row.names = NULL
  )
  
  # Nested 3-parameter GGD (d = 0)
  neg_ll3_smooth <- function(par3) {
    a <- par3[[1L]]; b <- par3[[2L]]; c <- par3[[3L]]
    if (a <= 1e-8 || b <= 1e-8 || c <= 1e-8) return(NaN)
    lr <- log(x / a); lz <- b * lr
    if (any(lz > 600)) return(NaN)
    z  <- exp(lz)
    -sum(c * lz - z - lgamma(c) + log(b) - log(x))
  }
  
  starts3 <- list(c(median(x), 1.5, 1.0), c(quantile(x, 0.80), 2.5, 0.5))
  best_opt3 <- NULL; best_ll3 <- -Inf
  for (st3 in starts3) {
    opt3 <- tryCatch(
      optim(st3, neg_ll3_smooth, method = "Nelder-Mead",
            control = list(maxit = 40000, reltol = 1e-12, parscale = pmax(st3, 1e-4))),
      error = function(e) NULL
    )
    if (!is.null(opt3) && is.finite(opt3$value) && -opt3$value > best_ll3) {
      best_ll3 <- -opt3$value; best_opt3 <- opt3
    }
  }
  est3 <- setNames(best_opt3$par, c("a", "b", "c"))
  ll3  <- best_ll3
  
  calc_ic <- function(ll, k) {
    aic  <- 2 * k - 2 * ll
    aicc <- aic + (2 * k * (k + 1)) / max(1, n - k - 1)
    bic  <- k * log(n) - 2 * ll
    hqic <- 2 * k * log(log(n)) - 2 * ll
    round(c(LogLik = ll, AIC = aic, AICc = aicc, BIC = bic, HQIC = hqic), 4)
  }
  
  lrt_stat <- max(0, 2 * (ll4 - ll3))
  lrt_pval <- pchisq(lrt_stat, df = 1, lower.tail = FALSE)
  
  list(
    estimates     = est4,
    ggd_estimates = est3,
    summary_table = df4,
    fisher_matrix = round(H4, 4),
    criteria      = data.frame(
      Criterion = c("LogLik", "AIC", "AICc", "BIC", "HQIC"),
      GGDSF     = calc_ic(ll4, 4),
      GGD       = calc_ic(ll3, 3),
      row.names = NULL
    ),
    lrt           = c(Statistic = round(lrt_stat, 4), df = 1, p_value = signif(lrt_pval, 5))
  )
}

#' Maximum Likelihood Estimation for Right-Censored GGDSF Samples
#'
#' Fits the four-parameter GGDSF distribution and the nested three-parameter
#' GGD sub-model (\code{d = 0}) under right censoring.
#'
#' @param time Numeric vector of strictly positive observed times (failures or censoring times).
#' @param status Numeric event indicator vector (\code{1} = failure, \code{0} = right-censored).
#' @param init Numeric vector of initial parameter values \code{c(a, b, c, d)}.
#' @return A list containing \code{estimates}, \code{ggd_estimates}, \code{summary_table},
#'   \code{fisher_matrix}, \code{criteria}, and \code{lrt}.
#' @importFrom stats optim qnorm pchisq pgamma
#' @importFrom numDeriv hessian
#' @examples
#' data(aarset)
#' y_obs  <- pmin(aarset, 75)
#' delta  <- as.numeric(aarset <= 75)
#' res_c  <- fit_ggdsf_censored(y_obs, delta)
#' print(res_c$summary_table)
#' print(res_c$lrt)
#' @export
fit_ggdsf_censored <- function(time, status, init = c(80.0, 0.4, 0.8, 0.035)) {
  if (length(time) != length(status)) stop("'time' and 'status' must have equal length.")
  if (any(time <= 0)) stop("All observed times must be strictly positive.")
  if (!all(status %in% c(0, 1))) stop("'status' must be 1 (failure) or 0 (right-censored).")
  
  n <- length(time)
  r <- sum(status)
  
  neg_ll_smooth <- function(par) {
    a <- par[[1L]]; b <- par[[2L]]; c <- par[[3L]]; d <- par[[4L]]
    if (a <= 1e-8 || b <= 1e-8 || c <= 1e-8 || d < -1e-8) return(NaN)
    lr <- log(time / a); lz <- (b + d * time) * lr
    if (any(lz > 600)) return(NaN)
    z  <- exp(lz)
    g  <- b + d * time + d * time * lr
    if (any(g[status == 1] <= 0) || any(!is.finite(z))) return(NaN)
    
    ll_fail <- sum(c * lz[status == 1] - z[status == 1] - lgamma(c) +
                     log(g[status == 1]) - log(time[status == 1]))
    ll_cens <- sum(pgamma(z[status == 0], shape = c, rate = 1, lower.tail = FALSE, log.p = TRUE))
    ll <- ll_fail + ll_cens
    if (!is.finite(ll)) return(NaN)
    -ll
  }
  
  neg_ll_const <- function(par) {
    if (par[[2L]] < (par[[4L]] * par[[1L]]) / exp(2)) return(1e12)
    val <- neg_ll_smooth(par)
    if (is.na(val)) 1e12 else val
  }
  
  opt1 <- optim(init, neg_ll_const, method = "Nelder-Mead",
                control = list(maxit = 20000, reltol = 1e-12, parscale = pmax(init, 0.05)))
  opt2 <- optim(opt1$par, neg_ll_const, method = "BFGS",
                control = list(maxit = 10000, reltol = 1e-12, parscale = pmax(opt1$par, 0.05)))
  
  est4 <- setNames(opt2$par, c("a", "b", "c", "d"))
  ll4  <- -opt2$value
  
  H4  <- numDeriv::hessian(neg_ll_smooth, est4)
  colnames(H4) <- rownames(H4) <- c("a", "b", "c", "d")
  se4 <- tryCatch(sqrt(diag(solve(H4))), error = function(e) rep(NA_real_, 4))
  
  df4 <- data.frame(
    Model     = c("GGDSF", "", "", ""),
    Parameter = c("a", "b", "c", "d"),
    MLE       = round(est4, 4),
    SE        = round(se4, 4),
    CI_Lower  = round(est4 - qnorm(0.975) * se4, 4),
    CI_Upper  = round(est4 + qnorm(0.975) * se4, 4),
    Wald_z    = round(est4 / se4, 2),
    row.names = NULL
  )
  
  # Nested GGD (d = 0)
  opt_g1 <- optim(c(90.0, 2.0, 0.4), function(p) neg_ll_const(c(p[[1L]], p[[2L]], p[[3L]], 0)),
                  method = "Nelder-Mead", control = list(maxit = 20000, reltol = 1e-12))
  opt_g2 <- optim(opt_g1$par, function(p) neg_ll_const(c(p[[1L]], p[[2L]], p[[3L]], 0)),
                  method = "BFGS", control = list(maxit = 10000, reltol = 1e-12))
  est3 <- setNames(opt_g2$par, c("a", "b", "c"))
  ll3  <- -opt_g2$value
  
  calc_ic <- function(ll, k) {
    aic     <- 2 * k - 2 * ll
    aicc    <- aic + (2 * k * (k + 1)) / max(1, n - k - 1)
    bic     <- k * log(n) - 2 * ll
    bic_eff <- k * log(max(1, r)) - 2 * ll
    hqic    <- 2 * k * log(log(n)) - 2 * ll
    round(c(LogLik = ll, AIC = aic, AICc = aicc, BIC = bic, BIC_eff = bic_eff, HQIC = hqic), 4)
  }
  
  lrt_stat <- max(0, 2 * (ll4 - ll3))
  lrt_pval <- pchisq(lrt_stat, df = 1, lower.tail = FALSE)
  
  list(
    estimates     = est4,
    ggd_estimates = est3,
    summary_table = df4,
    fisher_matrix = round(H4, 4),
    criteria      = data.frame(
      Criterion = c("LogLik", "AIC", "AICc", "BIC", "BIC_eff", "HQIC"),
      GGDSF     = calc_ic(ll4, 4),
      GGD       = calc_ic(ll3, 3),
      row.names = NULL
    ),
    lrt           = c(Statistic = round(lrt_stat, 4), df = 1, p_value = signif(lrt_pval, 5))
  )
}
