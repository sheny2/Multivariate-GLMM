# ---- Two-level multivariate DGP (site -> patient) for two outcomes

generate_log_data_mv <- function(
    seed   = 123,
    H      = 5,
    m_hosp = sample(10:30, H, replace = TRUE),
    px     = 6,
    p_bin  = 3,
    beta   = NULL,           
    sigma_a = 1,
    sigma_b = 1,
    rho     = 0.6
){
  stopifnot(p_bin <= px)
  set.seed(seed)
  
  nn      <- m_hosp
  N       <- sum(nn)
  site_id <- rep(seq_len(H), times = nn)
  pat_id  <- sequence(nn)
  p_cont  <- px - p_bin
  
  # site-level binary covariates
  X_site <- if (p_bin > 0) {
    site_cov <- matrix(rbinom(H * p_bin, 1, 0.3), nrow = H, ncol = p_bin)
    site_cov[site_id, , drop = FALSE]
  } else NULL
  
  # patient-level continuous covariates
  X_pat <- if (p_cont > 0) matrix(rnorm(N * p_cont, 0, 1), nrow = N, ncol = p_cont) else NULL
  
  X <- cbind(X_pat, X_site)
  if (is.null(X)) X <- matrix(nrow = N, ncol = 0)
  colnames(X) <- if (ncol(X) > 0) paste0("X", seq_len(ncol(X))) else character(0)
  
  if (is.null(beta)) {
    beta <- matrix(runif(px * 2, -0.7, 0.7), nrow = px, ncol = 2)
  }
  if (ncol(X) != nrow(beta)) stop("ncol(X) must equal nrow(beta) = px.")
  
  # site random effects: U ~ MVN(0, Sigma)
  Sigma <- matrix(c(sigma_a^2, rho * sigma_a * sigma_b,
                    rho * sigma_a * sigma_b, sigma_b^2), 2, 2)
  u_h <- MASS::mvrnorm(H, mu = c(0, 0), Sigma = Sigma)
  
  # random effect per outcome
  u1 <- u_h[site_id, 1]
  u2 <- u_h[site_id, 2]
  
  eta1 <- u1 + as.vector(X %*% beta[, 1])
  eta2 <- u2 + as.vector(X %*% beta[, 2])
  
  y1 <- rbinom(N, 1, plogis(eta1))
  y2 <- rbinom(N, 1, plogis(eta2))
  
  dat <- data.frame(site = site_id, patient = pat_id, X, Y1 = y1, Y2 = y2, check.names = FALSE)
  return(dat)
}

# Example:
dat <- generate_log_data_mv(seed = 123, H = 10, m_hosp = rep(100,10))
head(dat)