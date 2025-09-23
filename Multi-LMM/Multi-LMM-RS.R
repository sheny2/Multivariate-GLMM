# Profile-likelihood fit (two-level multivariate LMM)
# Site random INTERCEPT + shared random SLOPES across outcomes (direct Sigma)
# Residuals: exchangeable across outcomes, iid across patients within site
peal_profile_exch_RI_RS <- function(
    Y, X, id.site,
    X_rs = NULL,            # optional: N x p matrix used for random slopes (patient-space)
    rs_idx = NULL,          # optional: which columns of X_rs have random slopes; default = all
    init = list(beta = NULL,
                log_sigma2 = NULL,
                log_tau0_2 = NULL,
                log_tau2 = NULL,
                eta_rho = 0),
    control = list()
) {

  if (!is.matrix(Y)) stop("Y must be an N x R numeric matrix.")
  storage.mode(Y) <- "double"
  N <- nrow(Y); R <- ncol(Y)
  if (length(id.site) != N) stop("id.site must have length nrow(Y).")

  # ---- Fixed-effect design (block-diagonal across outcomes) ----
  if (is.list(X)) {
    if (length(X) != R) stop("If X is a list, length must equal ncol(Y).")
    X_list <- X
    # If X_rs not provided, require identical X across outcomes to reuse it
    if (is.null(X_rs)) {
      all_same <- all(vapply(X_list, function(A) identical(A, X_list[[1]]), TRUE))
      if (!all_same)
        stop("Provide X_rs or pass identical X across outcomes to reuse as random-slope design.")
      X_rs_common <- X_list[[1]]
    } else {
      if (!is.matrix(X_rs) || nrow(X_rs) != N) stop("X_rs must be an N x p matrix.")
      X_rs_common <- X_rs
    }
  } else if (is.matrix(X)) {
    X_list <- replicate(R, X, simplify = FALSE)
    X_rs_common <- if (is.null(X_rs)) X else X_rs
    if (nrow(X_rs_common) != N) stop("X_rs must be an N x p matrix.")
  } else stop("X must be a matrix or a list of length R.")

  p_r <- vapply(X_list, ncol, 0L)
  for (r in seq_len(R)) if (nrow(X_list[[r]]) != N) stop("Each X[[r]] must have nrow == nrow(Y).")
  p_tot <- sum(p_r)

  # Determine which columns have random slopes
  p_rs <- ncol(X_rs_common)
  if (is.null(rs_idx)) {
    rs_idx <- if (p_rs > 0) seq_len(p_rs) else integer(0)
  } else {
    if (length(rs_idx) && (max(rs_idx) > p_rs || min(rs_idx) < 1))
      stop("rs_idx contains out-of-range column indices for X_rs.")
  }
  k_rs <- length(rs_idx)  # number of slope covariates with random slopes

  # ---- Site indexing ----
  site_fac <- as.factor(id.site)
  sites    <- levels(site_fac)
  H        <- length(sites)
  site_idx_list <- split(seq_len(N), site_fac)

  # ---- Build per-site blocks ----
  build_site_blocks <- function(idx_rows) {
    m  <- length(idx_rows)
    Yh <- Y[idx_rows, , drop = FALSE]         # m x R
    yv <- as.numeric(Yh)                      # outcome-major vec (length m*R)

    Xh_list <- lapply(X_list, function(Xr) Xr[idx_rows, , drop = FALSE])
    nc <- vapply(Xh_list, ncol, 0L)
    Xb <- matrix(0.0, nrow = m * R, ncol = sum(nc))
    col_off <- 0L
    for (r in seq_len(R)) {
      rows <- ((r - 1L) * m + 1L):(r * m)
      cols <- (col_off + 1L):(col_off + nc[r])
      Xb[rows, cols] <- Xh_list[[r]]
      col_off <- col_off + nc[r]
    }

    # Random-effect (patient-space) design: Z = [1, X_rs[, rs_idx]]
    Z_h <- if (k_rs > 0) {
      cbind(1.0, X_rs_common[idx_rows, rs_idx, drop = FALSE])
    } else {
      matrix(1.0, nrow = m, ncol = 1L)  # random intercept only
    }

    list(y = yv, X = Xb, Z = Z_h, m = m)
  }
  site_blocks <- lapply(site_idx_list, build_site_blocks)

  # ---- Initial values ----
  if (is.null(init$beta)) {
    beta0 <- numeric(p_tot)
    off <- 0L
    for (r in seq_len(R)) {
      br <- tryCatch(qr.coef(qr(X_list[[r]]), Y[, r]),
                     error = function(e) rep(0, ncol(X_list[[r]])))
      beta0[(off + 1L):(off + length(br))] <- br
      off <- off + length(br)
    }
  } else {
    beta0 <- as.numeric(init$beta)
    if (length(beta0) != p_tot) stop("init$beta has wrong length.")
  }

  y_all <- as.numeric(Y)
  vy <- stats::var(y_all)
  log_sigma2_0 <- if (is.null(init$log_sigma2)) log(vy * 0.6 + 1e-6) else init$log_sigma2
  log_tau0_2_0 <- if (is.null(init$log_tau0_2)) log(vy * 0.2 + 1e-6) else init$log_tau0_2
  log_tau2_0   <- if (is.null(init$log_tau2))   log(vy * 0.2 + 1e-6) else init$log_tau2

  # Map eta_rho -> rho in (L, U) to ensure PD of exchangeable residuals
  Lrho <- -1 / (R - 1) + 1e-8
  Urho <-  1 - 1e-8
  eta_rho_0 <- if (is.null(init$eta_rho)) 0 else init$eta_rho

  unpack_par <- function(phi) {
    log_sigma2 <- phi[1]; log_tau0_2 <- phi[2]; log_tau2 <- phi[3]; eta_rho <- phi[4]
    sigma2 <- exp(log_sigma2)
    tau0_2 <- exp(log_tau0_2)      # random intercept variance
    tau2   <- exp(log_tau2)        # shared random-slope variance
    rho    <- Lrho + (Urho - Lrho) * stats::plogis(eta_rho)
    list(sigma2 = sigma2, tau0_2 = tau0_2, tau2 = tau2, rho = rho)
  }

  # ---- Helpers: build Sigma_h (direct), safe Cholesky, solve, logdet ----
  Sigma_site_RI_RS <- function(m, R, sigma2, rho, tau0_2, tau2, Z) {
    # Z: m x q, with first column for intercept, remaining for selected slopes
    q <- ncol(Z)
    J_R <- matrix(1.0, R, R)
    I_R <- diag(1.0, R)
    I_m <- diag(1.0, m)

    # residual part
    B_R <- sigma2 * ((1 - rho) * I_R + rho * J_R)            # R x R

    # random effects part: Z D Z^T
    if (q == 1L) {
      # intercept-only
      Srand <- tau0_2 * (Z %*% t(Z))                         # m x m
    } else {
      D     <- diag(c(tau0_2, rep(tau2, q - 1L)), nrow = q)  # q x q
      Srand <- Z %*% D %*% t(Z)                              # m x m
    }

    kronecker(B_R, I_m) + kronecker(J_R, Srand)              # (mR) x (mR)
  }

  safe_chol <- function(S) {
    jitter <- 0.0
    for (k in 0:6) {
      out <- try(chol(S + jitter * diag(nrow(S))), silent = TRUE)
      if (!inherits(out, "try-error")) return(list(ch = out, jitter = jitter))
      jitter <- if (k == 0) 1e-8 else jitter * 10
    }
    stop("Cholesky failed for Sigma_h; matrix may be ill-conditioned.")
  }

  chol_solve <- function(ch, B) {  # solve (t(ch) %*% ch) x = B
    y <- forwardsolve(t(ch), B, upper.tri = FALSE, transpose = FALSE)
    backsolve(ch, y, upper.tri = TRUE, transpose = FALSE)
  }

  chol_logdet <- function(ch) 2 * sum(log(diag(ch)))

  # ---- Negative profiled log-likelihood ----
  npll <- function(phi, return_beta = FALSE) {
    par <- unpack_par(phi)
    sigma2 <- par$sigma2; tau0_2 <- par$tau0_2; tau2 <- par$tau2; rho <- par$rho

    XtVinvX <- matrix(0.0, nrow = p_tot, ncol = p_tot)
    XtVinvy <- numeric(p_tot)
    logdet  <- 0.0
    const   <- 0.0

    for (h in seq_len(H)) {
      bh <- site_blocks[[h]]
      yv <- bh$y
      Xb <- bh$X
      m  <- bh$m
      Z  <- bh$Z

      Sh   <- Sigma_site_RI_RS(m, R, sigma2, rho, tau0_2, tau2, Z)
      ch   <- safe_chol(Sh)$ch
      ldSh <- chol_logdet(ch)

      Wy <- chol_solve(ch, yv)
      WX <- chol_solve(ch, Xb)

      XtVinvX <- XtVinvX + crossprod(Xb, WX)
      XtVinvy <- XtVinvy + crossprod(Xb, Wy)
      logdet  <- logdet + ldSh
      const   <- const + m * R * log(2 * pi)
    }

    # GLS beta-hat
    beta_hat <- tryCatch(solve(XtVinvX, XtVinvy),
                         error = function(e) qr.solve(XtVinvX, XtVinvy))

    # Quadratic term
    quad <- 0.0
    for (h in seq_len(H)) {
      bh <- site_blocks[[h]]
      yv <- bh$y
      Xb <- bh$X
      m  <- bh$m
      Z  <- bh$Z

      Sh <- Sigma_site_RI_RS(m, R, sigma2, rho, tau0_2, tau2, Z)
      ch <- safe_chol(Sh)$ch

      rv  <- yv - as.numeric(Xb %*% beta_hat)
      Wrv <- chol_solve(ch, rv)
      quad <- quad + sum(rv * Wrv)
    }

    ll <- -0.5 * (const + logdet + quad)
    if (return_beta) return(list(npll = -ll, beta_hat = beta_hat))
    -ll
  }

  # ---- Optimize variance parameters (profile over beta) ----
  phi0 <- c(log_sigma2_0, log_tau0_2_0, log_tau2_0, eta_rho_0)
  opt  <- optim(phi0, npll, method = "BFGS", control = control)

  par_hat    <- unpack_par(opt$par)
  beta_hat   <- npll(opt$par, return_beta = TRUE)$beta_hat
  sigma2_hat <- par_hat$sigma2
  tau0_2_hat <- par_hat$tau0_2
  tau2_hat   <- par_hat$tau2
  rho_hat    <- par_hat$rho
  logLik_hat <- -opt$value

  structure(list(
    beta        = beta_hat,
    sigma2      = sigma2_hat,
    tau0_2_int  = tau0_2_hat,   # site random-intercept variance
    tau2_slope  = tau2_hat,     # shared random-slope variance (for columns in rs_idx)
    rho         = rho_hat,
    logLik      = as.numeric(logLik_hat),
    convergence = opt$convergence,
    counts      = opt$counts,
    message     = opt$message,
    R = R,
    p_per_outcome = p_r,
    rs_idx = rs_idx
  ), class = "peal_profile_exch_RI_RS_direct")
}

# Pretty printer
print.peal_profile_exch_RI_RS_direct <- function(x, ...) {
  cat("peal_profile_exch_RI_RS_direct fit (site RI + shared RS; direct Sigma)\n")
  cat(sprintf("  Outcomes (R): %d  |  Random slopes on %d columns\n",
              x$R, length(x$rs_idx)))
  cat(sprintf("  logLik: %.3f  |  convergence: %s\n",
              x$logLik, as.character(x$convergence)))
  cat(sprintf("  sigma^2: %.6f,  tau0^2(RI): %.6f,  tau^2(RS): %.6f,  rho: %.6f\n",
              x$sigma2, x$tau0_2_int, x$tau2_slope, x$rho))
  cat("  beta (stacked by outcomes):\n")
  print(drop(x$beta))
  invisible(x)
}
