#' Profile-likelihood fit of a two-level multivariate LMM
#' - Random intercept at site only: y_hir = x_hir^T beta_r + b_h + e_hir
#' - Exchangeable residual covariance across outcomes:
#'     Cov(e_hi) = sigma^2 * [(1 - rho) I_R + rho * 1_R 1_R^T]
#'
#' INPUTS
#'   Y        : N x R numeric matrix (rows = patients, cols = outcomes)
#'   X        : Either
#'                (a) a list of length R, each X[[r]] is N x p_r (outcome-specific design), or
#'                (b) a single N x p matrix (same covariates for all outcomes)
#'   id.site  : length-N vector/factor of site IDs (patients -> sites)
#'   Z        : (ignored; kept for API compatibility; random intercept is fixed by id.site)
#'   init     : optional list with init$beta, init$log_sigma2, init$log_tau2, init$eta_rho
#'   control  : list passed to optim (e.g., list(maxit=200))
#'
#' RETURNS
#'   list with: beta, sigma2, tau2, rho, logLik, convergence, counts, message
peal_profile_exch_RI <- function(Y, X, id.site, Z = NULL,
                                 init = list(beta = NULL,
                                             log_sigma2 = NULL,
                                             log_tau2 = NULL,
                                             eta_rho = 0),
                                 control = list()) {

  if (!is.matrix(Y)) stop("Y must be an N x R numeric matrix.")
  storage.mode(Y) <- "double"
  N <- nrow(Y); R <- ncol(Y)
  if (length(id.site) != N) stop("id.site must have length nrow(Y).")

  # Parse X into a list over outcomes
  X_list <- NULL
  if (is.list(X)) {
    if (length(X) != R) stop("If X is a list, it must have length equal to ncol(Y).")
    X_list <- X
  } else if (is.matrix(X)) {
    # same design for all outcomes
    X_list <- replicate(R, X, simplify = FALSE)
  } else {
    stop("X must be either a matrix (common design) or a list of length R (per-outcome design).")
  }
  # Basic checks
  p_r <- vapply(X_list, ncol, 0L)
  for (r in seq_len(R)) {
    if (nrow(X_list[[r]]) != N) stop("Each X[[r]] must have nrow == nrow(Y).")
  }
  p_tot <- sum(p_r)

  # Make site index
  site_fac <- as.factor(id.site)
  sites <- levels(site_fac)
  H <- length(sites)

  # ---- Helpers: outcome-major stacking and P/Q projections ----
  # We treat vectors in outcome-major order: [y^(1); ...; y^(R)] where each y^(r) is m_h-long.
  # For a site with m rows (patients) and R outcomes, represent v as m x R matrix U (columns = outcomes).
  Vinv_apply_site_vec <- function(v, m, R, sigma2, rho, tau2) {
    # v is length m*R in outcome-major order (column-wise vec of m x R matrix)
    U <- matrix(v, nrow = m, ncol = R)

    # Projections:
    # Pr(U): row-wise mean across outcomes, replicated across columns
    PrU <- matrix(rowMeans(U), nrow = m, ncol = R, byrow = FALSE)
    # Pm(U): column-wise mean across patients, replicated across rows
    PmU <- matrix(colMeans(U), nrow = m, ncol = R, byrow = TRUE)
    # grand mean
    g <- mean(U)

    a <- matrix(g, nrow = m, ncol = R)            # Pm(Pr(U))
    b <- PrU - a                                   # Qm(Pr(U))
    c <- PmU - a                                   # Pm(Qr(U))
    d <- U - PrU - PmU + a                         # Qm(Qr(U))

    alpha  <- sigma2 * (1 - rho)
    beta   <- sigma2 * rho
    lambdaA <- alpha + beta * R + tau2 * R * m     # (P_R ⊗ P_m), mult 1
    lambdaB <- alpha + beta * R                    # (P_R ⊗ Q_m), mult m-1
    lambdaC <- alpha                               # (Q_R ⊗ P_m), mult R-1
    lambdaD <- alpha                               # (Q_R ⊗ Q_m), mult (R-1)(m-1)

    Uout <- a / lambdaA + b / lambdaB + c / lambdaC + d / lambdaD
    as.numeric(Uout)
  }

  Vinv_apply_site_mat <- function(M, m, R, sigma2, rho, tau2) {
    # Apply Vinv to each column of M (dimension m*R x k) using the vector function
    k <- ncol(M)
    out <- matrix(0.0, nrow = m * R, ncol = k)
    for (j in seq_len(k)) out[, j] <- Vinv_apply_site_vec(M[, j], m, R, sigma2, rho, tau2)
    out
  }

  # Log |Sigma_h| via eigenvalues (see derivation)
  logdet_site <- function(m, R, sigma2, rho, tau2) {
    alpha  <- sigma2 * (1 - rho)
    beta   <- sigma2 * rho
    if (alpha <= 0) return(-Inf) # invalid
    # multiplicities: (R-1)*m, (m-1), 1
    (R - 1) * m * log(alpha) +
      (m - 1) * log(alpha + beta * R) +
      log(alpha + beta * R + tau2 * R * m)
  }

  # Build per-site blocks of (y_h, X_h) in outcome-major order
  # y_h_vec: length m_h*R
  # X_h_big: (m_h*R) x (sum p_r) block-diagonal across outcomes
  build_site_blocks <- function(idx_rows) {
    m <- length(idx_rows)
    # y (m x R) => outcome-major vector
    Yh <- Y[idx_rows, , drop = FALSE]
    y_h_vec <- as.numeric(Yh) # column-major vec = outcome-major stacking
    # X block-diag
    Xh_list <- lapply(X_list, function(Xr) Xr[idx_rows, , drop = FALSE])
    # Build block-diagonal dense matrix
    # Use base R to avoid extra deps
    blocks <- Xh_list
    nr <- vapply(blocks, nrow, 0L)
    nc <- vapply(blocks, ncol, 0L)
    stopifnot(all(nr == m))
    X_h_big <- matrix(0.0, nrow = m * R, ncol = sum(nc))
    col_off <- 0L
    for (r in seq_len(R)) {
      cols <- seq_len(nc[r]) + col_off
      rows <- ((r - 1L) * m + 1L):(r * m)
      X_h_big[rows, cols] <- blocks[[r]]
      col_off <- col_off + nc[r]
    }
    list(y = y_h_vec, X = X_h_big, m = m)
  }

  site_idx_list <- split(seq_len(N), site_fac)
  site_blocks <- lapply(site_idx_list, build_site_blocks)  # list of length H

  # ---- Initial values ----
  if (is.null(init$beta)) {
    # OLS per outcome
    beta0 <- numeric(p_tot)
    off <- 0L
    for (r in seq_len(R)) {
      Xr <- X_list[[r]]
      yr <- Y[, r]
      br <- tryCatch({
        qr.coef(qr(Xr), yr)
      }, error = function(e) rep(0, ncol(Xr)))
      beta0[(off + 1L):(off + length(br))] <- br
      off <- off + length(br)
    }
  } else {
    beta0 <- as.numeric(init$beta)
    if (length(beta0) != p_tot) stop("init$beta has wrong length.")
  }

  # Variance params:
  y_all <- as.numeric(Y)
  vy <- stats::var(y_all)
  log_sigma2_0 <- if (is.null(init$log_sigma2)) log(vy * 0.8 + 1e-6) else init$log_sigma2
  log_tau2_0   <- if (is.null(init$log_tau2))   log(vy * 0.2 + 1e-6) else init$log_tau2
  # Map eta_rho -> rho in (L,U) to respect PD of exchangeable residuals
  Lrho <- -1 / (R - 1) + 1e-8
  Urho <-  1 - 1e-8
  eta_rho_0 <- if (is.null(init$eta_rho)) 0 else init$eta_rho

  # Transformations
  unpack_par <- function(phi) {
    log_sigma2 <- phi[1]; log_tau2 <- phi[2]; eta_rho <- phi[3]
    sigma2 <- exp(log_sigma2)
    tau2   <- exp(log_tau2)
    rho    <- Lrho + (Urho - Lrho) * stats::plogis(eta_rho)
    list(sigma2 = sigma2, tau2 = tau2, rho = rho)
  }

  # Evaluate negative profiled log-likelihood at phi
  # Also optionally return beta_hat if needed
  npll <- function(phi, return_beta = FALSE) {
    par <- unpack_par(phi)
    sigma2 <- par$sigma2; tau2 <- par$tau2; rho <- par$rho

    # Build X^T V^{-1} X and X^T V^{-1} y by summing site-wise
    XtVinvX <- matrix(0.0, nrow = p_tot, ncol = p_tot)
    XtVinvy <- numeric(p_tot)
    logdet  <- 0.0

    for (h in seq_len(H)) {
      bh <- site_blocks[[h]]
      yv <- bh$y
      Xb <- bh$X
      m  <- bh$m

      Wy <- Vinv_apply_site_vec(yv, m, R, sigma2, rho, tau2)
      WX <- Vinv_apply_site_mat(Xb, m, R, sigma2, rho, tau2)

      XtVinvX <- XtVinvX + crossprod(Xb, WX)
      XtVinvy <- XtVinvy + crossprod(Xb, Wy)
      logdet  <- logdet + logdet_site(m, R, sigma2, rho, tau2)
    }

    # GLS beta-hat (use qr for stability)
    beta_hat <- tryCatch({
      solve(XtVinvX, XtVinvy)
    }, error = function(e) {
      qr.solve(XtVinvX, XtVinvy)
    })

    # Compute quadratic form (y - X beta)^T V^{-1} (y - X beta)
    quad <- 0.0
    const <- 0.0
    off <- 0L
    for (h in seq_len(H)) {
      bh <- site_blocks[[h]]
      yv <- bh$y
      Xb <- bh$X
      m  <- bh$m
      rv <- yv - as.numeric(Xb %*% beta_hat)
      Wrv <- Vinv_apply_site_vec(rv, m, R, sigma2, rho, tau2)
      quad <- quad + sum(rv * Wrv)
      const <- const + m * R * log(2 * pi)
    }

    # log-likelihood
    ll <- -0.5 * (const + logdet + quad)
    if (return_beta) return(list(npll = -ll, beta_hat = beta_hat))
    -ll
  }

  # Optimize variance parameters (profile over beta)
  phi0 <- c(log_sigma2_0, log_tau2_0, eta_rho_0)
  opt <- optim(phi0, npll, method = "BFGS", control = control)

  # Final estimates
  par_hat <- unpack_par(opt$par)
  beta_hat <- npll(opt$par, return_beta = TRUE)$beta_hat
  sigma2_hat <- par_hat$sigma2
  tau2_hat   <- par_hat$tau2
  rho_hat    <- par_hat$rho

  # Final logLik
  logLik_hat <- -opt$value

  structure(list(
    beta   = beta_hat,
    sigma2 = sigma2_hat,
    sigma_u2 = tau2_hat,
    rho    = rho_hat,
    logLik = as.numeric(logLik_hat),
    convergence = opt$convergence,
    counts = opt$counts,
    message = opt$message,
    R = R,
    p_per_outcome = p_r
  ), class = "peal_profile_exch_RI")
}

# ---- (Optional) print method ----
print.peal_profile_exch_RI <- function(x, ...) {
  cat("peal_profile_exch_RI fit\n")
  cat(sprintf("  Outcomes (R): %d\n", x$R))
  cat(sprintf("  logLik: %.3f  |  convergence: %s\n",
              x$logLik, as.character(x$convergence)))
  cat(sprintf("  sigma^2: %.6f,  sigma_u^2: %.6f,  rho: %.6f\n",
              x$sigma2, x$sigma_u2, x$rho))
  cat("  beta (stacked by outcomes):\n")
  print(drop(x$beta))
  invisible(x)
}
