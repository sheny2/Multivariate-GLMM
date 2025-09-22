# ---- Two-level multivariate DGP (site -> patient),
library(MASS)
library(data.table)

generate_data_mv <- function(
    seed = 123,
    H = 5,                                   # number of sites
    m_hosp = sample(10:30, H, replace = TRUE),  # patients per site (length H)
    px = 6,                                  # number of covariates
    p_bin = 3,                               # number of binary covariates
    py = 3,                                  # number of outcomes
    beta = matrix(runif(px*py, 5, 10), nrow = px, ncol = py), # fixed effects (px x py)
    sigma_u = 0.2,                           # site RE sd (intercept)
    # ---- residual structure (across outcomes) ----
    sigma_e = 0.3,                           # residual sd (per outcome)
    rho = 0.3,                               # exchangeable residual correlation across outcomes
    # ---- random slopes at site level ----
    rs_idx = 1:px,                            # indices of X with site-level random slopes
    tau = 0.15,                              # sd(s) for random slopes (scalar or length(rs_idx))
    slope_by_outcome = FALSE                 # if TRUE, slopes differ by outcome; else shared across outcomes
){
  set.seed(seed)
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install 'MASS'.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install 'data.table'.")

  # Sizes and IDs
  nn <- m_hosp
  if (length(nn) != H) stop("m_hosp must have length H.")
  N  <- sum(nn)
  site_id <- rep(seq_len(H), times = nn)
  pat_id  <- sequence(nn)

  # Covariates (patient-level)
  p_cont <- px - p_bin
  X_bin  <- matrix(rbinom(N * p_bin, 1, 0.3), nrow = N, ncol = p_bin)
  X_cont <- matrix(rnorm(N * p_cont, 0, 1),   nrow = N, ncol = p_cont)
  X      <- cbind(X_bin, X_cont)  # N x px
  colnames(X) <- paste0("X", seq_len(px))

  # Site-level random intercepts
  u_h <- rnorm(H, mean = 0, sd = sigma_u)    # length H
  u_per_row <- u_h[site_id]                  # length N

  # Handle random slope settings
  rs_idx <- as.integer(rs_idx)
  if (length(rs_idx) > 0 && (any(rs_idx < 1) || any(rs_idx > px))) {
    stop("rs_idx must be in 1:px.")
  }
  if (length(rs_idx) == 0) tau <- numeric(0)
  if (length(tau) == 1 && length(rs_idx) > 1) tau <- rep(tau, length(rs_idx))
  if (length(tau) != length(rs_idx)) stop("tau must be scalar or length(rs_idx).")

  # Build residual covariance (exchangeable across outcomes)
  rho_mat <- matrix(rho, nrow = py, ncol = py); diag(rho_mat) <- 1
  Sigma_eps <- (sigma_e^2) * rho_mat

  # Multivariate residuals per patient
  eps <- MASS::mvrnorm(n = N, mu = rep(0, py), Sigma = Sigma_eps)  # N x py

  # Site-level random slopes
  # If slope_by_outcome = FALSE:
  #   for each site h and each j in rs_idx, draw b_{h,j} ~ N(0, tau_j^2),
  #   and apply the SAME slope across all outcomes.
  # If slope_by_outcome = TRUE:
  #   for each site h, each outcome k, and each j in rs_idx,
  #   draw b_{h,j,k} ~ N(0, tau_j^2) independently.
  #
  # (You can later extend to correlated slopes by sampling from a multivariate normal.)
  if (length(rs_idx) > 0) {
    if (!slope_by_outcome) {
      # H x |rs_idx|
      B_h <- matrix(rnorm(H * length(rs_idx), 0, rep(tau, each = H)),
                    nrow = H, ncol = length(rs_idx), byrow = FALSE)
      # For each row, compute sum_j X_ij * b_{site(i), j}
      X_rs <- X[, rs_idx, drop = FALSE]                # N x R
      b_per_row <- B_h[site_id, , drop = FALSE]        # N x R
      RS_effect <- X_rs * b_per_row                    # elementwise
      RS_sum <- rowSums(RS_effect)                     # N
      # Shared across outcomes -> add RS_sum to each outcome column later
    } else {
      # outcome-specific slopes: array H x |rs_idx| x py
      B_hk <- array(
        rnorm(H * length(rs_idx) * py, 0,
              rep(tau, each = H * py)),
        dim = c(H, length(rs_idx), py)
      )
      X_rs <- X[, rs_idx, drop = FALSE]               # N x R
      # For each outcome k, compute RS contribution
      RS_sum_k <- matrix(0, nrow = N, ncol = py)
      for (k in seq_len(py)) {
        b_h_k <- B_hk[, , k, drop = FALSE]            # H x R
        b_per_row_k <- b_h_k[site_id, , drop = FALSE] # N x R
        RS_effect_k <- X_rs * b_per_row_k             # N x R
        RS_sum_k[, k] <- rowSums(RS_effect_k)         # N
      }
    }
  }

  # Outcomes
  Y <- matrix(0, nrow = N, ncol = py)
  for (k in seq_len(py)) {
    mu_fixed <- as.numeric(X %*% beta[, k]) + u_per_row
    if (length(rs_idx) == 0) {
      mu <- mu_fixed
    } else if (!slope_by_outcome) {
      mu <- mu_fixed + RS_sum
    } else {
      mu <- mu_fixed + RS_sum_k[, k]
    }
    Y[, k] <- mu + eps[, k]
  }
  colnames(Y) <- paste0("Y", seq_len(py))

  # Final dataset: one row per patient
  dat <- data.table::data.table(
    site = site_id,
    patient = pat_id
  )
  dat <- cbind(dat, X, Y)
  data.table::setnames(dat,
                       old = c("site", "patient", colnames(X), colnames(Y)),
                       new = c("site", "patient", paste0("X", 1:px), paste0("Y", 1:py)))
  return(dat[])
}
