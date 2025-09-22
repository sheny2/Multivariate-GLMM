peal_profile_exch_RI_RS <- function(
    Y, X, id.site,
    X_rs = NULL,                     # optional: explicit random-slope design; if NULL, deduce from X
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

  # ----- Fixed-effect design parsing -----
  if (is.list(X)) {
    if (length(X) != R) stop("If X is a list, its length must equal ncol(Y).")
    X_list <- X
    # If random-slope design not provided explicitly, require identical X across outcomes
    if (is.null(X_rs)) {
      all_same <- all(vapply(X_list, function(A) identical(A, X_list[[1]]), TRUE))
      if (!all_same) stop("Provide X_rs or pass identical X for all outcomes for random slopes.")
      X_rs_common <- X_list[[1]]
    } else {
      if (!is.matrix(X_rs) || nrow(X_rs) != N) stop("X_rs must be an N x p matrix.")
      X_rs_common <- X_rs
    }
  } else if (is.matrix(X)) {
    X_list <- replicate(R, X, simplify = FALSE)
    X_rs_common <- if (is.null(X_rs)) X else X_rs
  } else {
    stop("X must be a matrix or a list of length R.")
  }

  p_r  <- vapply(X_list, ncol, 0L)
  for (r in seq_len(R)) if (nrow(X_list[[r]]) != N) stop("Each X[[r]] must have nrow == nrow(Y).")
  p_tot <- sum(p_r)

  # Random-slope design Z = [1, X_rs]
  if (!is.matrix(X_rs_common)) stop("Random-slope design (X_rs or common X) must be a matrix.")
  p  <- ncol(X_rs_common)              # number of slope covariates
  # D = diag(tau0^2, tau^2, ..., tau^2) has size (1+p)
  # We'll pass D via its diagonal vector ddiag = c(tau0^2, rep(tau^2, p))

  # ----- Site indexing -----
  site_fac <- as.factor(id.site)
  sites <- levels(site_fac)
  H <- length(sites)
  site_idx_list <- split(seq_len(N), site_fac)

  # ----- Per-site blocks (y_h, X_h^big, Z_h) -----
  build_site_blocks <- function(idx_rows) {
    m <- length(idx_rows)
    Yh <- Y[idx_rows, , drop = FALSE]
    y_h_vec <- as.numeric(Yh) # outcome-major (column-major) vec of m x R

    # Fixed-effect block-diagonal across outcomes
    Xh_list <- lapply(X_list, function(Xr) Xr[idx_rows, , drop = FALSE])
    nr <- vapply(Xh_list, nrow, 0L); stopifnot(all(nr == m))
    nc <- vapply(Xh_list, ncol, 0L)
    X_h_big <- matrix(0.0, nrow = m * R, ncol = sum(nc))
    col_off <- 0L
    for (r in seq_len(R)) {
      rows <- ((r - 1L) * m + 1L):(r * m)
      cols <- (col_off + 1L):(col_off + nc[r])
      X_h_big[rows, cols] <- Xh_list[[r]]
      col_off <- col_off + nc[r]
    }

    # Random-effect (site) design in patient space: Z = [1, X_rs]
    Xh_rs <- X_rs_common[idx_rows, , drop = FALSE]  # m x p
    Z_h   <- cbind(1.0, Xh_rs)                      # m x (1+p)

    list(y = y_h_vec, X = X_h_big, Z = Z_h, m = m)
  }
  site_blocks <- lapply(site_idx_list, build_site_blocks)

  # ----- Initial values -----
  if (is.null(init$beta)) {
    beta0 <- numeric(p_tot)
    off <- 0L
    for (r in seq_len(R)) {
      Xr <- X_list[[r]]
      yr <- Y[, r]
      br <- tryCatch(qr.coef(qr(Xr), yr), error = function(e) rep(0, ncol(Xr)))
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

  # Map eta_rho -> rho in (L,U) to ensure PD of exchangeable residuals
  Lrho <- -1 / (R - 1) + 1e-8
  Urho <-  1 - 1e-8
  eta_rho_0 <- if (is.null(init$eta_rho)) 0 else init$eta_rho

  unpack_par <- function(phi) {
    log_sigma2 <- phi[1]; log_tau0_2 <- phi[2]; log_tau2 <- phi[3]; eta_rho <- phi[4]
    sigma2 <- exp(log_sigma2)
    tau0_2 <- exp(log_tau0_2)
    tau2   <- exp(log_tau2)
    rho    <- Lrho + (Urho - Lrho) * stats::plogis(eta_rho)
    list(sigma2 = sigma2, tau0_2 = tau0_2, tau2 = tau2, rho = rho)
  }

  # ----- Core linear ops on a single site's (m x R) matrix U -----
  # Apply (c I_m + Z D Z^T)^{-1} to U (column-wise)
  apply_block_inv_mat <- function(U, Z, c, ddiag) {
    # ddiag: length (1+p) diagonal of D
    m <- nrow(U); q <- ncol(Z)
    if (q != length(ddiag)) stop("Z and D diag length mismatch.")
    # M = D^{-1} + (1/c) Z^T Z  (size (1+p) x (1+p))
    # do everything in double; guard tiny variances
    eps <- 1e-10
    Dinv <- diag(1 / pmax(ddiag, eps), nrow = q, ncol = q)
    ZtZ  <- crossprod(Z)                   # (1+p) x (1+p)
    M    <- Dinv + (1 / c) * ZtZ
    Minv <- tryCatch(solve(M), error = function(e) qr.solve(M))
    XtU  <- crossprod(Z, U)                # (1+p) x R
    corr <- Z %*% (Minv %*% XtU)           # m x R
    (U / c) - (corr / (c * c))
  }

  Vinv_apply_site_vec <- function(v, m, R, sigma2, rho, tau0_2, tau2, Z) {
    U <- matrix(v, nrow = m, ncol = R)
    alpha <- sigma2 * (1 - rho)
    beta  <- sigma2 * rho
    ddiag <- c(tau0_2, rep(tau2, ncol(Z) - 1))

    row_mean <- rowMeans(U)
    Ubar <- matrix(row_mean, nrow = m, ncol = R)   # mean across outcomes (P space)
    Ures <- U - Ubar                               # residual across outcomes (Q space)

    Ures_out <- apply_block_inv_mat(Ures, Z, c = alpha, ddiag = ddiag)
    Ubar_out <- apply_block_inv_mat(Ubar, Z, c = alpha + beta * R, ddiag = ddiag)

    as.numeric(Ures_out + Ubar_out)
  }

  Vinv_apply_site_mat <- function(M, m, R, sigma2, rho, tau0_2, tau2, Z) {
    k <- ncol(M)
    out <- matrix(0.0, nrow = m * R, ncol = k)
    for (j in seq_len(k)) out[, j] <- Vinv_apply_site_vec(M[, j], m, R, sigma2, rho, tau0_2, tau2, Z)
    out
  }

  # log|V_h| using small (1+p)-dim determinants
  logdet_site <- function(m, R, sigma2, rho, tau0_2, tau2, Z) {
    alpha <- sigma2 * (1 - rho)
    beta  <- sigma2 * rho
    ddiag <- c(tau0_2, rep(tau2, ncol(Z) - 1))
    ZtZ   <- crossprod(Z)
    # helper: log|c I + Z D Z^T| = m log c + log| I + (1/c) D Z^T Z |
    logdet_c <- function(cval) {
      A <- diag(1 + 0 * ddiag, nrow = length(ddiag)) + (1 / cval) * (diag(ddiag) %*% ZtZ)
      m * log(cval) + as.numeric(determinant(A, logarithm = TRUE)$modulus)
    }
    (R - 1) * logdet_c(alpha) + logdet_c(alpha + beta * R)
  }

  # ----- Negative profiled log-likelihood -----
  npll <- function(phi, return_beta = FALSE) {
    par <- unpack_par(phi)
    sigma2 <- par$sigma2; tau0_2 <- par$tau0_2; tau2 <- par$tau2; rho <- par$rho

    XtVinvX <- matrix(0.0, nrow = p_tot, ncol = p_tot)
    XtVinvy <- numeric(p_tot)
    logdet  <- 0.0

    for (h in seq_len(H)) {
      bh <- site_blocks[[h]]
      yv <- bh$y
      Xb <- bh$X
      m  <- bh$m
      Z  <- bh$Z

      Wy <- Vinv_apply_site_vec(yv, m, R, sigma2, rho, tau0_2, tau2, Z)
      WX <- Vinv_apply_site_mat(Xb, m, R, sigma2, rho, tau0_2, tau2, Z)

      XtVinvX <- XtVinvX + crossprod(Xb, WX)
      XtVinvy <- XtVinvy + crossprod(Xb, Wy)
      logdet  <- logdet + logdet_site(m, R, sigma2, rho, tau0_2, tau2, Z)
    }

    beta_hat <- tryCatch(solve(XtVinvX, XtVinvy),
                         error = function(e) qr.solve(XtVinvX, XtVinvy))

    quad <- 0.0
    const <- 0.0
    for (h in seq_len(H)) {
      bh <- site_blocks[[h]]
      yv <- bh$y
      Xb <- bh$X
      m  <- bh$m
      Z  <- bh$Z
      rv <- yv - as.numeric(Xb %*% beta_hat)
      Wrv <- Vinv_apply_site_vec(rv, m, R, sigma2, rho, tau0_2, tau2, Z)
      quad <- quad + sum(rv * Wrv)
      const <- const + m * R * log(2 * pi)
    }

    ll <- -0.5 * (const + logdet + quad)
    if (return_beta) return(list(npll = -ll, beta_hat = beta_hat))
    -ll
  }

  # ----- Optimize variance parameters (profile over beta) -----
  phi0 <- c(log_sigma2_0, log_tau0_2_0, log_tau2_0, eta_rho_0)
  opt <- optim(phi0, npll, method = "BFGS", control = control)

  par_hat    <- unpack_par(opt$par)
  beta_hat   <- npll(opt$par, return_beta = TRUE)$beta_hat
  sigma2_hat <- par_hat$sigma2
  tau0_2_hat <- par_hat$tau0_2
  tau2_hat   <- par_hat$tau2
  rho_hat    <- par_hat$rho
  logLik_hat <- -opt$value

  structure(list(
    beta       = beta_hat,
    sigma2     = sigma2_hat,
    tau0_2_int = tau0_2_hat,    # site random intercept variance
    tau2_slope = tau2_hat,      # site random slope variance (common across p slopes)
    rho        = rho_hat,
    logLik     = as.numeric(logLik_hat),
    convergence = opt$convergence,
    counts      = opt$counts,
    message     = opt$message,
    R = R,
    p_per_outcome = p_r,
    p_random = 1 + p
  ), class = "peal_profile_exch_RI_RS")
}

print.peal_profile_exch_RI_RS <- function(x, ...) {
  cat("peal_profile_exch_RI_RS fit (site RI + site RS shared across outcomes)\n")
  cat(sprintf("  Outcomes (R): %d | Random dim (1+p): %d\n", x$R, x$p_random))
  cat(sprintf("  logLik: %.3f  |  convergence: %s\n",
              x$logLik, as.character(x$convergence)))
  cat(sprintf("  sigma^2: %.6f,  tau0^2(RI): %.6f,  tau^2(RS): %.6f,  rho: %.6f\n",
              x$sigma2, x$tau0_2_int, x$tau2_slope, x$rho))
  cat("  beta (stacked by outcomes):\n")
  print(drop(x$beta))
  invisible(x)
}
