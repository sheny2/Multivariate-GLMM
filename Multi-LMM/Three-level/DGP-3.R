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


generate_data_mv3 <- function(
    seed = 123,
    H = 5,                                   # number of sites
    m_hosp = sample(10:30, H, replace = TRUE),  # patients per site (length H)
    px = 6,                                  # number of covariates
    p_bin = 3,                               # number of binary covariates
    py = 3,                                  # number of outcomes
    beta = matrix(stats::runif(px*py, 5, 10), nrow = px, ncol = py), # fixed effects (px x py)
    sigma_u = 0.3,                           # site RE sd (intercept)
    # ---- patient random intercept (site-specific SDs) ----
    sigma_v_hosp = NULL,                     # length-H vector; if NULL, runif(H, 0.5, 0.7)
    # ---- residual structure (across outcomes) ----
    sigma_e = 0.5,                           # residual sd (same for all outcomes)
    rho = 0.5,                               # exchangeable residual correlation across outcomes
    # ---- visits ----
    visits_range = c(1, 20),                 # inclusive range for number of visits per patient
    # ---- optional random slopes at site level (like your previous template) ----
    rs_idx = integer(0),                     # indices of X with site-level random slopes (subset of 1:px)
    tau = 0.15,                              # sd(s) for random slopes (scalar or length(rs_idx))
    slope_by_outcome = FALSE,                # if TRUE, slopes differ by outcome; else shared across outcomes
    # ---- outcome-specific X? (your original DGP shared X across outcomes) ----
    outcome_specific_X = FALSE               # if TRUE, draw a different X for each outcome (still stored once)
) {
  set.seed(seed)
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install 'MASS'.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install 'data.table'.")

  # --- sizes & ids ---
  nn <- m_hosp
  if (length(nn) != H) stop("m_hosp must have length H.")
  N_pat  <- sum(nn)                               # number of patients across all sites
  site_id_pat <- rep(seq_len(H), times = nn)      # length N_pat
  pat_id_within_site <- sequence(nn)              # length N_pat

  # visits per patient
  if (length(visits_range) != 2 || visits_range[1] < 1) stop("visits_range must be c(min,max), min>=1.")
  n_visits <- sample(visits_range[1]:visits_range[2], N_pat, replace = TRUE)
  N <- sum(n_visits)                              # total rows (visits)

  # expand to visit level
  id.visit <- sequence(n_visits)
  id.hosp.expanded <- rep(site_id_pat, times = n_visits)         # length N
  id.pat.expanded  <- rep(pat_id_within_site, times = n_visits)  # length N

  # --- random effects ---
  # site intercepts
  u_h <- stats::rnorm(H, mean = 0, sd = sigma_u)                     # length H
  u_h_patient <- u_h[site_id_pat]                                    # per patient
  u_h_expanded <- rep(u_h_patient, times = n_visits)                 # per visit

  # patient intercepts with site-specific SDs
  if (is.null(sigma_v_hosp)) sigma_v_hosp <- stats::runif(H, 0.5, 0.7)
  if (length(sigma_v_hosp) != H) stop("sigma_v_hosp must be length H.")
  sd_per_patient <- sigma_v_hosp[site_id_pat]                        # length N_pat
  v_hi <- stats::rnorm(N_pat, mean = 0, sd = sd_per_patient)         # per patient
  v_hi_expanded <- rep(v_hi, times = n_visits)                       # per visit

  # --- covariates at visit level ---
  p_cont <- px - p_bin
  if (p_cont < 0) stop("p_bin must be <= px.")
  # shared X across outcomes (your original DGP)
  X_bin  <- matrix(stats::rbinom(N * p_bin, size = 1, prob = 0.3), nrow = N, ncol = p_bin)
  X_cont <- matrix(stats::rnorm(N * p_cont, mean = 0, sd = 1),      nrow = N, ncol = p_cont)
  X_hij  <- cbind(X_bin, X_cont)                                    # N x px
  colnames(X_hij) <- paste0("X", seq_len(px))

  # (optional) outcome-specific X (values differ by outcome)
  if (isTRUE(outcome_specific_X)) {
    X_list <- vector("list", py)
    for (k in seq_len(py)) {
      Xk_bin  <- matrix(stats::rbinom(N * p_bin, size = 1, prob = 0.3), nrow = N, ncol = p_bin)
      Xk_cont <- matrix(stats::rnorm(N * p_cont, mean = 0, sd = 1),      nrow = N, ncol = p_cont)
      X_list[[k]] <- cbind(Xk_bin, Xk_cont)  # N x px
      colnames(X_list[[k]]) <- paste0("X", seq_len(px))
    }
  }

  # --- residual errors across outcomes (exchangeable) ---
  rho_mat <- matrix(rho, nrow = py, ncol = py); diag(rho_mat) <- 1
  Sigma_eps <- (sigma_e^2) * rho_mat
  epsilon_hij <- MASS::mvrnorm(N, mu = rep(0, py), Sigma = Sigma_eps) # N x py

  # --- optional site-level random slopes (on covariates in rs_idx) ---
  RS_sum_shared <- rep(0, N)              # used if slope_by_outcome = FALSE
  RS_sum_k      <- matrix(0, N, py)       # used if slope_by_outcome = TRUE
  if (length(rs_idx) > 0) {
    if (any(rs_idx < 1 | rs_idx > px)) stop("rs_idx must be subset of 1:px.")
    if (length(tau) == 1) tau <- rep(tau, length(rs_idx))
    if (length(tau) != length(rs_idx)) stop("tau must be scalar or length(rs_idx).")

    if (!slope_by_outcome) {
      # draw site-by-slope effects, shared across outcomes
      B_h <- matrix(stats::rnorm(H * length(rs_idx), 0, rep(tau, each = H)),
                    nrow = H, ncol = length(rs_idx), byrow = FALSE)      # H x |rs_idx|
      b_row <- B_h[id.hosp.expanded, , drop = FALSE]                      # N x |rs_idx|
      RS_sum_shared <- rowSums(X_hij[, rs_idx, drop = FALSE] * b_row)
    } else {
      # outcome-specific slopes
      B_hk <- array(
        stats::rnorm(H * length(rs_idx) * py, 0, rep(tau, each = H * py)),
        dim = c(H, length(rs_idx), py)
      ) # H x |rs_idx| x py
      for (k in seq_len(py)) {
        b_row_k <- B_hk[id.hosp.expanded, , k, drop = FALSE]              # N x |rs_idx|
        RS_sum_k[, k] <- rowSums(X_hij[, rs_idx, drop = FALSE] * b_row_k)
      }
    }
  }

  # --- outcomes per visit ---
  Y_hij <- matrix(0, nrow = N, ncol = py)
  for (k in seq_len(py)) {
    Xk <- if (isTRUE(outcome_specific_X)) X_list[[k]] else X_hij
    mu_fixed <- as.numeric(Xk %*% beta[, k]) + u_h_expanded + v_hi_expanded
    mu <- if (length(rs_idx) == 0) {
      mu_fixed
    } else if (!slope_by_outcome) {
      mu_fixed + RS_sum_shared
    } else {
      mu_fixed + RS_sum_k[, k]
    }
    Y_hij[, k] <- mu + epsilon_hij[, k]
  }
  colnames(Y_hij) <- paste0("Y", seq_len(py))

  # --- build data.tables (wide per-visit + rearranged like yours) ---
  DT <- data.table::data.table(
    site    = id.hosp.expanded,
    patient = id.pat.expanded,
    visit   = id.visit
  )
  DT <- cbind(DT, X_hij, Y_hij)

  # visit counts per patient (for your downstream ordering)
  visit_count <- DT[, .(total_visits = .N), by = .(site, patient)]

  # merge & arrange
  rearranged_data <- merge(DT, visit_count, by = c("site","patient"))
  setorder(rearranged_data, site, total_visits, patient, visit)
  rearranged_data[, site := factor(site)]

  # --- also provide a long version (one y column) if you want it ready-to-go ---
  # (uses the X shared across outcomes if outcome_specific_X = FALSE;
  #  if TRUE, we store X columns as X1..Xpx from the shared baseline X_hij to keep a simple schema)
  Y_long <- as.vector(Y_hij)
  outcome <- rep(paste0("Y", seq_len(py)), each = N)
  long_df <- data.frame(
    site    = rearranged_data$site,
    patient = rearranged_data$patient,
    visit   = rearranged_data$visit,
    outcome = factor(outcome, levels = paste0("Y", seq_len(py))),
    y       = Y_long
  )
  long_df <- cbind(long_df, as.data.frame(X_hij))   # attach X columns once
  names(long_df)[names(long_df) %in% colnames(X_hij)] <- paste0("x", seq_len(px)) # lowercase x*

  # return
  list(
    data              = rearranged_data,   # wide per visit (X1..Xpx, Y1..Ypy)
    long_df           = long_df,           # stacked outcomes (one y column)
    beta              = beta,
    sigma_u           = sigma_u,
    sigma_v_hosp      = sigma_v_hosp,
    Sigma_eps         = Sigma_eps,
    rs_idx            = rs_idx,
    tau               = tau,
    slope_by_outcome  = slope_by_outcome,
    outcome_specific_X = outcome_specific_X,
    meta = list(
      H = H, px = px, p_bin = p_bin, py = py,
      visits_range = visits_range,
      N_patients = N_pat, N_rows = N
    )
  )
}



generate_data_mv3 <- function(
    seed = 123,
    H = 5,                                   # number of sites
    m_hosp = sample(10:30, H, replace = TRUE),  # patients per site (length H)
    px = 6,                                  # number of covariates
    p_bin = 3,                               # number of binary covariates
    py = 3,                                  # number of outcomes
    beta = matrix(stats::runif(px*py, 5, 10), nrow = px, ncol = py), # fixed effects (px x py)
    sigma_u = 0.3,                           # site RE sd (intercept)
    # ---- patient random intercept (site-specific SDs) ----
    sigma_v_hosp = NULL,                     # length-H vector; if NULL, runif(H, 0.5, 0.7)
    # ---- residual structure (across outcomes) ----
    sigma_e = 0.5,                           # residual sd (same for all outcomes)
    rho = 0.5,                               # exchangeable residual correlation across outcomes
    # ---- visits ----
    visits_range = c(1, 20),                 # inclusive range for number of visits per patient
    # ---- optional random slopes at site level (like your previous template) ----
    rs_idx = integer(0),                     # indices of X with site-level random slopes (subset of 1:px)
    tau = 0.15,                              # sd(s) for random slopes (scalar or length(rs_idx))
    slope_by_outcome = FALSE,                # if TRUE, slopes differ by outcome; else shared across outcomes
    # ---- outcome-specific X? (your original DGP shared X across outcomes) ----
    outcome_specific_X = FALSE               # if TRUE, draw a different X for each outcome (still stored once)
) {
  set.seed(seed)
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Please install 'MASS'.")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install 'data.table'.")

  # --- sizes & ids ---
  nn <- m_hosp
  if (length(nn) != H) stop("m_hosp must have length H.")
  N_pat  <- sum(nn)                               # number of patients across all sites
  site_id_pat <- rep(seq_len(H), times = nn)      # length N_pat
  pat_id_within_site <- sequence(nn)              # length N_pat

  # visits per patient
  if (length(visits_range) != 2 || visits_range[1] < 1) stop("visits_range must be c(min,max), min>=1.")
  n_visits <- sample(visits_range[1]:visits_range[2], N_pat, replace = TRUE)
  N <- sum(n_visits)                              # total rows (visits)

  # expand to visit level
  id.visit <- sequence(n_visits)
  id.hosp.expanded <- rep(site_id_pat, times = n_visits)         # length N
  id.pat.expanded  <- rep(pat_id_within_site, times = n_visits)  # length N

  # --- random effects ---
  # site intercepts
  u_h <- stats::rnorm(H, mean = 0, sd = sigma_u)                     # length H
  u_h_patient <- u_h[site_id_pat]                                    # per patient
  u_h_expanded <- rep(u_h_patient, times = n_visits)                 # per visit

  # patient intercepts with site-specific SDs
  if (is.null(sigma_v_hosp)) sigma_v_hosp <- stats::runif(H, 0.5, 0.7)
  if (length(sigma_v_hosp) != H) stop("sigma_v_hosp must be length H.")
  sd_per_patient <- sigma_v_hosp[site_id_pat]                        # length N_pat
  v_hi <- stats::rnorm(N_pat, mean = 0, sd = sd_per_patient)         # per patient
  v_hi_expanded <- rep(v_hi, times = n_visits)                       # per visit

  # --- covariates at visit level ---
  p_cont <- px - p_bin
  if (p_cont < 0) stop("p_bin must be <= px.")
  # shared X across outcomes (your original DGP)
  X_bin  <- matrix(stats::rbinom(N * p_bin, size = 1, prob = 0.3), nrow = N, ncol = p_bin)
  X_cont <- matrix(stats::rnorm(N * p_cont, mean = 0, sd = 1),      nrow = N, ncol = p_cont)
  X_hij  <- cbind(X_bin, X_cont)                                    # N x px
  colnames(X_hij) <- paste0("X", seq_len(px))

  # (optional) outcome-specific X (values differ by outcome)
  if (isTRUE(outcome_specific_X)) {
    X_list <- vector("list", py)
    for (k in seq_len(py)) {
      Xk_bin  <- matrix(stats::rbinom(N * p_bin, size = 1, prob = 0.3), nrow = N, ncol = p_bin)
      Xk_cont <- matrix(stats::rnorm(N * p_cont, mean = 0, sd = 1),      nrow = N, ncol = p_cont)
      X_list[[k]] <- cbind(Xk_bin, Xk_cont)  # N x px
      colnames(X_list[[k]]) <- paste0("X", seq_len(px))
    }
  }

  # --- residual errors across outcomes (exchangeable) ---
  rho_mat <- matrix(rho, nrow = py, ncol = py); diag(rho_mat) <- 1
  Sigma_eps <- (sigma_e^2) * rho_mat
  epsilon_hij <- MASS::mvrnorm(N, mu = rep(0, py), Sigma = Sigma_eps) # N x py

  # --- optional site-level random slopes (on covariates in rs_idx) ---
  RS_sum_shared <- rep(0, N)              # used if slope_by_outcome = FALSE
  RS_sum_k      <- matrix(0, N, py)       # used if slope_by_outcome = TRUE
  if (length(rs_idx) > 0) {
    if (any(rs_idx < 1 | rs_idx > px)) stop("rs_idx must be subset of 1:px.")
    if (length(tau) == 1) tau <- rep(tau, length(rs_idx))
    if (length(tau) != length(rs_idx)) stop("tau must be scalar or length(rs_idx).")

    if (!slope_by_outcome) {
      # draw site-by-slope effects, shared across outcomes
      B_h <- matrix(stats::rnorm(H * length(rs_idx), 0, rep(tau, each = H)),
                    nrow = H, ncol = length(rs_idx), byrow = FALSE)      # H x |rs_idx|
      b_row <- B_h[id.hosp.expanded, , drop = FALSE]                      # N x |rs_idx|
      RS_sum_shared <- rowSums(X_hij[, rs_idx, drop = FALSE] * b_row)
    } else {
      # outcome-specific slopes
      B_hk <- array(
        stats::rnorm(H * length(rs_idx) * py, 0, rep(tau, each = H * py)),
        dim = c(H, length(rs_idx), py)
      ) # H x |rs_idx| x py
      for (k in seq_len(py)) {
        b_row_k <- B_hk[id.hosp.expanded, , k, drop = FALSE]              # N x |rs_idx|
        RS_sum_k[, k] <- rowSums(X_hij[, rs_idx, drop = FALSE] * b_row_k)
      }
    }
  }

  # --- outcomes per visit ---
  Y_hij <- matrix(0, nrow = N, ncol = py)
  for (k in seq_len(py)) {
    Xk <- if (isTRUE(outcome_specific_X)) X_list[[k]] else X_hij
    mu_fixed <- as.numeric(Xk %*% beta[, k]) + u_h_expanded + v_hi_expanded
    mu <- if (length(rs_idx) == 0) {
      mu_fixed
    } else if (!slope_by_outcome) {
      mu_fixed + RS_sum_shared
    } else {
      mu_fixed + RS_sum_k[, k]
    }
    Y_hij[, k] <- mu + epsilon_hij[, k]
  }
  colnames(Y_hij) <- paste0("Y", seq_len(py))

  # --- build data.tables (wide per-visit + rearranged like yours) ---
  DT <- data.table::data.table(
    site    = id.hosp.expanded,
    patient = id.pat.expanded,
    visit   = id.visit
  )
  DT <- cbind(DT, X_hij, Y_hij)

  # visit counts per patient (for your downstream ordering)
  visit_count <- DT[, .(total_visits = .N), by = .(site, patient)]

  # merge & arrange
  rearranged_data <- merge(DT, visit_count, by = c("site","patient"))
  setorder(rearranged_data, site, total_visits, patient, visit)
  rearranged_data[, site := factor(site)]

  # --- also provide a long version (one y column) if you want it ready-to-go ---
  # (uses the X shared across outcomes if outcome_specific_X = FALSE;
  #  if TRUE, we store X columns as X1..Xpx from the shared baseline X_hij to keep a simple schema)
  Y_long <- as.vector(Y_hij)
  outcome <- rep(paste0("Y", seq_len(py)), each = N)
  long_df <- data.frame(
    site    = rearranged_data$site,
    patient = rearranged_data$patient,
    visit   = rearranged_data$visit,
    outcome = factor(outcome, levels = paste0("Y", seq_len(py))),
    y       = Y_long
  )
  long_df <- cbind(long_df, as.data.frame(X_hij))   # attach X columns once
  names(long_df)[names(long_df) %in% colnames(X_hij)] <- paste0("x", seq_len(px)) # lowercase x*

  # return
  list(
    data              = rearranged_data,   # wide per visit (X1..Xpx, Y1..Ypy)
    long_df           = long_df,           # stacked outcomes (one y column)
    beta              = beta,
    sigma_u           = sigma_u,
    sigma_v_hosp      = sigma_v_hosp,
    Sigma_eps         = Sigma_eps,
    rs_idx            = rs_idx,
    tau               = tau,
    slope_by_outcome  = slope_by_outcome,
    outcome_specific_X = outcome_specific_X,
    meta = list(
      H = H, px = px, p_bin = p_bin, py = py,
      visits_range = visits_range,
      N_patients = N_pat, N_rows = N
    )
  )
}


