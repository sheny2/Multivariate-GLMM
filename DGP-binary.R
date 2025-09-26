# ---- Two-level multivariate DGP (site -> patient) with LOGIT outcomes ----
# Allows outcome-specific site RI variances: sigma_u can be scalar or length-py vector

expit <- function(z) 1 / (1 + exp(-z))

generate_data_mv_bin <- function(
    seed   = 123,
    H      = 5,                                   # number of sites
    m_hosp = sample(10:30, H, replace = TRUE),    # patients per site (length H)
    px     = 6,                                   # number of covariates
    p_bin  = 3,                                   # number of binary covariates
    py     = 2,                                   # number of outcomes
    beta   = matrix(runif(px*py, -1.0, 1.0), nrow = px, ncol = py),  # fixed effects
    sigma_u = c(0.3,0.6),                                # site RI SD: scalar OR length-py vector
    sigma   = 0.9,                                # residual SD (logit scale)
    rho     = 0.4                                 # exchangeable residual corr across outcomes
) {
  set.seed(seed)
  stopifnot(is.matrix(beta), ncol(beta) == py)

  # Sizes and IDs
  nn <- m_hosp
  N  <- sum(nn)
  site_id <- rep(seq_len(H), times = nn)
  pat_id  <- sequence(nn)

  # Design matrix X (patient-level)
  p_cont <- px - p_bin
  X_bin  <- matrix(rbinom(N * p_bin, 1, 0.3), nrow = N, ncol = p_bin)
  X_cont <- matrix(rnorm(N * p_cont, 0, 1),   nrow = N, ncol = p_cont)
  X      <- cbind(X_bin, X_cont)

  # -----------------------------
  # Site random intercepts u_{hr}
  # -----------------------------
  # Allow sigma_u to be scalar or length-py vector
  if (length(sigma_u) == 1L) sigma_u <- rep(sigma_u, py)
  stopifnot(length(sigma_u) == py, all(sigma_u >= 0))

  # For simplicity: independent across outcomes (diag covariance per site)
  # u_h: H x py, with column r ~ N(0, sigma_u[r]^2)
  u_h <- matrix(0, nrow = H, ncol = py)
  for (r in seq_len(py)) u_h[, r] <- rnorm(H, 0, sigma_u[r])
  # Map to patients: u is N x py where each row i uses the site of that patient
  u <- u_h[site_id, , drop = FALSE]  # N x py

  # ------------------------------------------------
  # Residuals across outcomes: exchangeable covariance
  # D = sigma^2 * [(1 - rho) I + rho * 11^T]
  # ------------------------------------------------
  D <- sigma^2 * ((1 - rho) * diag(py) + rho * matrix(1, py, py))
  # mvrnorm requires MASS
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install MASS.")
  res_mat <- MASS::mvrnorm(n = N, mu = rep(0, py), Sigma = D)  # N x py

  # Linear predictors & probabilities (logit link)
  eta <- matrix(NA_real_, nrow = N, ncol = py)
  for (r in seq_len(py)) {
    eta[, r] <- X %*% beta[, r] + u[, r] + res_mat[, r]
  }
  P <- expit(eta)

  # Binary outcomes
  Y <- matrix(rbinom(N * py, size = 1, prob = as.vector(P)), nrow = N, ncol = py)

  # Assemble data.table
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table.")
  dat <- data.table::data.table(site = site_id, patient = pat_id, X, Y)
  data.table::setnames(dat, c("site","patient", paste0("X", 1:px), paste0("Y", 1:py)))

  # Attach useful attributes
  attr(dat, "P")        <- P
  attr(dat, "eta")      <- eta
  attr(dat, "u_site")   <- u            # N x py
  attr(dat, "u_site_H") <- u_h          # H x py
  attr(dat, "res_pat")  <- res_mat
  attr(dat, "D")        <- D
  attr(dat, "sigma_u")  <- sigma_u
  return(dat[])
}

# Example
# dat_bin <- generate_data_mv_bin(seed = 123, py = 3, sigma_u = c(0.4, 0.7, 1.0))
# head(dat_bin)
