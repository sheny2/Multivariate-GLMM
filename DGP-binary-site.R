expit <- function(z) 1 / (1 + exp(-z))

generate_data_mv_bin_site <- function(
    seed    = 123,
    H       = 5,                                   # sites
    m_hosp  = sample(10:30, H, replace = TRUE),    # patients per site
    px      = 6,                                   # patient-level covariates
    p_bin   = 3,                                   # of which are binary
    py      = 2,                                   # outcomes
    
    # NEW: site-level covariates (qz total; qz_bin binary, qz_cont continuous)
    qz      = 2,
    qz_bin  = 1,
    
    beta    = matrix(runif(px*py, -1.0, 1.0), nrow = px, ncol = py),
    gamma   = matrix(runif(max(qz,1)*py, -0.5, 0.5), nrow = max(qz,1), ncol = py), # NEW: site-level coeffs
    
    sigma_u = c(0.3, 0.6),                         # site RI SDs (length py or scalar)
    sigma   = 0.9,                                 # residual SD (logit scale)
    rho     = 0.4                                  # residual corr across outcomes
) {
  set.seed(seed)
  stopifnot(is.matrix(beta), ncol(beta) == py)
  stopifnot(qz_bin <= qz)
  qz_cont <- qz - qz_bin
  
  # Sizes and IDs
  nn <- m_hosp
  N  <- sum(nn)
  site_id <- rep(seq_len(H), times = nn)
  pat_id  <- sequence(nn)
  
  # Patient-level X
  p_cont <- px - p_bin
  X_bin  <- if (p_bin > 0) matrix(rbinom(N * p_bin, 1, 0.3), nrow = N) else NULL
  X_cont <- if (p_cont > 0) matrix(rnorm(N * p_cont), nrow = N) else NULL
  X      <- cbind(X_bin, X_cont)
  if (is.null(X)) X <- matrix(numeric(N*0), nrow = N, ncol = 0)
  
  # NEW: Site-level Z (defined once per site, then mapped to patients)
  Zb <- if (qz_bin > 0) matrix(rbinom(H * qz_bin, 1, 0.4), nrow = H) else NULL
  Zc <- if (qz_cont > 0) matrix(rnorm(H * qz_cont), nrow = H) else NULL
  ZH <- cbind(Zb, Zc)                                  # H x qz (may be 0 cols)
  if (is.null(ZH)) ZH <- matrix(numeric(H*0), nrow = H, ncol = 0)
  Z  <- ZH[site_id, , drop = FALSE]                    # N x qz
  if (qz == 0L) {
    # Keep shapes consistent when qz=0
    gamma <- matrix(numeric(0), nrow = 0, ncol = py)
  } else {
    stopifnot(nrow(gamma) == qz, ncol(gamma) == py)
  }
  
  # Site random intercepts u_{hr}
  if (length(sigma_u) == 1L) sigma_u <- rep(sigma_u, py)
  stopifnot(length(sigma_u) == py, all(sigma_u >= 0))
  u_h <- sapply(seq_len(py), function(r) rnorm(H, 0, sigma_u[r]))
  if (is.vector(u_h)) u_h <- matrix(u_h, nrow = H, ncol = py)
  u <- u_h[site_id, , drop = FALSE]                   # N x py
  
  # Residuals across outcomes: exchangeable covariance
  D <- sigma^2 * ((1 - rho) * diag(py) + rho * matrix(1, py, py))
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install MASS.")
  res_mat <- MASS::mvrnorm(n = N, mu = rep(0, py), Sigma = D)  # N x py
  
  # Linear predictors & probabilities
  eta <- matrix(NA_real_, nrow = N, ncol = py)
  for (r in seq_len(py)) {
    xb <- if (ncol(X) > 0) X %*% beta[, r, drop = FALSE] else 0
    zb <- if (ncol(Z) > 0) Z %*% gamma[, r, drop = FALSE] else 0   # NEW: site covariate effect
    eta[, r] <- xb + zb + u[, r] + res_mat[, r]
  }
  P <- expit(eta)
  
  # Binary outcomes
  Y <- matrix(rbinom(N * py, size = 1, prob = as.vector(P)), nrow = N, ncol = py)
  
  # Assemble data.table (replicate site-level Z to patient rows for convenience)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table.")
  dat <- data.table::data.table(site = site_id, patient = pat_id, X, Z, Y)
  data.table::setnames(dat,
                       c("site","patient",
                         if (ncol(X) > 0) paste0("X", seq_len(ncol(X))) else character(),
                         if (ncol(Z) > 0) paste0("Z", seq_len(ncol(Z))) else character(),
                         paste0("Y", seq_len(py))))
  
  # Attributes
  attr(dat, "P")        <- P
  attr(dat, "eta")      <- eta
  attr(dat, "u_site")   <- u
  attr(dat, "u_site_H") <- u_h
  attr(dat, "res_pat")  <- res_mat
  attr(dat, "D")        <- D
  attr(dat, "sigma_u")  <- sigma_u
  attr(dat, "Z_site")   <- ZH          # H x qz (canonical site-level values)
  attr(dat, "gamma")    <- gamma
  return(dat[])
}

# Example:
dat <- generate_data_mv_bin_site(H=8, m_hosp=sample(80:120, 8, TRUE), px=5, p_bin=2, py=3, sigma_u = c(0.4, 0.7, 1.0),
                            qz=2, qz_bin=1) 
head(dat)
