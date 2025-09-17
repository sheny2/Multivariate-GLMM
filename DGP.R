# ---- Two-level multivariate DGP (site -> patient),
library(MASS)
library(data.table)

generate_data_mv <- function(
    seed = 123,
    H = 5,                                # number of sites
    m_hosp = sample(10:30, H, replace=TRUE), # patients per site (vector length H)
    px = 6,                               # number of covariates
    p_bin = 3,                            # number of binary covariates
    py = 3,                               # number of outcomes
    beta = matrix(runif(px*py, 5, 10), nrow = px, ncol = py), # fixed effects
    sigma_u = 0.2,                        # site RE sd
    sigma_e = 0.3,                        # residual sd (per outcome)
    rho = 0.3                             # exchangeable residual correlation across outcomes
){
  set.seed(seed)

  # Sizes and IDs
  nn <- m_hosp                    # patients per site
  N  <- sum(nn)                   # total patients
  site_id <- rep(1:H, times = nn)
  pat_id  <- sequence(nn)

  # Random effects
  u_h <- rnorm(H, mean = 0, sd = sigma_u)                # site RE
  uh <- u_h[site_id]

  # Covariates (patient-level)
  p_cont <- px - p_bin
  X_bin  <- matrix(rbinom(N * p_bin, 1, 0.3), nrow = N, ncol = p_bin)
  X_cont <- matrix(rnorm(N * p_cont, 0, 1), nrow = N, ncol = p_cont)
  X      <- cbind(X_bin, X_cont)

  # Residual covariance (exchangeable)
  rho_mat <- matrix(rho, nrow = py, ncol = py); diag(rho_mat) <- 1
  Sigma_eps <- (sigma_e^2) * rho_mat

  # Multivariate residuals
  eps <- MASS::mvrnorm(n = N, mu = rep(0, py), Sigma = Sigma_eps)

  # Outcomes
  Y <- matrix(0, nrow = N, ncol = py)
  for(k in 1:py){
    Y[, k] <- X %*% beta[, k] + uh + eps[, k]
  }

  # Final dataset: one row per patient
  dat <- data.table(site = site_id, patient = pat_id, X, Y)
  setnames(dat, c("site", "patient", paste0("X", 1:px), paste0("Y", 1:py)))
  return(dat[])
}

# Example:
dat <- generate_data_mv(seed = 123)
head(dat)

