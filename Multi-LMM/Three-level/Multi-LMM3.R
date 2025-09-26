# -------------------------
# helpers for Z_hv
# -------------------------
generate_record_count <- function(data) {
  counts <- table(data[, "n_hi"])
  result_matrix <- cbind(as.numeric(names(counts)), as.numeric(counts))
  colnames(result_matrix) <- c("n_hi", "frequency")
  return(result_matrix)
}

generate_Zhv_matrix <- function(data) {
  record_count_matrix <- generate_record_count(data)
  diagonal_blocks <- list()
  for (i in 1:nrow(record_count_matrix)) {
    n_hi <- record_count_matrix[i, "n_hi"]
    frequency <- record_count_matrix[i, "frequency"]
    identity_block <- diag(frequency)
    ones_vector <- matrix(1, nrow = n_hi, ncol = 1)
    diagonal_blocks[[i]] <- kronecker(identity_block, ones_vector)
  }
  big_matrix <- do.call(Matrix::bdiag, diagonal_blocks)
  big_matrix <- cbind(1, big_matrix)
  return(as.matrix(big_matrix))
}

# -------------------------
# Exchangeable residual correlation utilities
# Rcorr = (1 - rho) I_R + rho * 11^T (scale-free)
# We profile sigma^2 separately; all whitening uses Rcorr^{-1}
# -------------------------
.Rcorr_inv_and_logdet <- function(R, rho, eps = 1e-8) {
  # bounds: rho in (-1/(R-1), 1)
  lower <- -1/(R-1) + eps
  upper <- 1 - eps
  if (rho <= lower || rho >= upper) stop("rho out of admissible range for exchangeable R-corr.")
  # eigenvalues of Rcorr: (1-ρ) with multiplicity (R-1); (1-ρ+ρR) once
  lam1 <- (1 - rho)               # multiplicity R-1
  lam2 <- (1 - rho + rho * R)     # multiplicity 1
  logdet_Rcorr <- (R-1)*log(lam1) + log(lam2)

  # closed form inverse: Rcorr^{-1} = (1/(1-ρ)) [ I - (ρ/(1-ρ+ρR)) 11^T ]
  a <- 1/(1 - rho)
  b <- rho / (1 - rho + rho * R)
  Rinve <- a * (diag(R) - b * matrix(1, R, R))
  list(Rinve = Rinve, logdet_Rcorr = logdet_Rcorr, lower = lower, upper = upper)
}


# -------------------------
# Direct-profile REML using X, Y, Z, id.site (no ShXYZ) but similar profile llk as peal
# Same parameterization:
#   par[1]     = theta_u          (site RI variance / sigma^2)
#   par[1 + h] = theta_vh (h=1..K) (patient RI variance by site / sigma^2)
#   [optional last] rho (if estimate_rho = TRUE)
# Returns list with lp, b, s2, rho, allterms (bterm1/bterm2/qterm/etc.)
# -------------------------
lmm.profile3_mv_DIRECT <- function(par,
                                   Y, X, Z, id.site,
                                   reml = TRUE,
                                   estimate_rho = TRUE, rho_fixed = 0,
                                   verbose = FALSE) {

  if (!is.matrix(Y)) stop("Y must be an N x R matrix.")
  N <- nrow(Y); R <- ncol(Y)

  # Accept X as N x p (shared across outcomes) or list length R of N x p_r
  if (is.list(X)) {
    if (length(X) != R) stop("If X is a list, length(X) must equal ncol(Y).")
    X_list <- X
  } else {
    X_list <- replicate(R, as.matrix(X), simplify = FALSE)
  }
  for (r in seq_len(R)) if (nrow(X_list[[r]]) != N) stop("Each X[[r]] must have nrow == nrow(Y).")

  # Sites and indices
  site_fac <- as.factor(id.site)
  sites    <- levels(site_fac)
  K        <- length(sites)
  idx_list <- split(seq_len(N), site_fac)

  # Parse parameter vector
  pz <- K + 1
  if (estimate_rho) {
    if (length(par) != pz + 1) stop("par must be K+2 when estimate_rho=TRUE.")
    rho <- par[pz + 1]
  } else {
    rho <- rho_fixed
    if (length(par) != pz) stop("par must be K+1 when rho is fixed.")
  }
  theta_u  <- par[1]
  theta_vh <- par[1 + seq_len(K)]

  # Exchangeable outcome correlation (scale-free)
  rc <- .Rcorr_inv_and_logdet(R, rho)
  Rinve <- rc$Rinve
  logdet_Rcorr <- rc$logdet_Rcorr

  # Dimensions for fixed effects (assume same px across outcomes for stacking)
  px_vec <- vapply(X_list, ncol, 0L)
  if (length(unique(px_vec)) != 1L) stop("All X[[r]] must have the same #columns for stacking.")
  px <- px_vec[1]
  pxR <- px * R
  n_tot <- N * R

  # Accumulators
  lpterm1 <- 0.0         # log|Gamma| pieces
  lpterm2 <- 0.0         # y~'Gamma^{-1} y~
  bterm1  <- matrix(0.0, pxR, pxR)  # X~'Gamma^{-1} X~
  bterm2  <- numeric(pxR)           # X~'Gamma^{-1} y~
  Nsum    <- 0L

  # Loop over sites; compute per-site crossproducts and apply Woodbury
  for (h in seq_len(K)) {
    idx <- idx_list[[h]]
    N_h <- length(idx); Nsum <- Nsum + N_h

    Yh <- Y[idx, , drop = FALSE]        # N_h x R
    # Build per-outcome X for site
    Xh_list <- lapply(X_list, function(Xr) Xr[idx, , drop = FALSE])

    # Z_h: N_h x (1 + m_h), first column = site RI, others = patient RIs
    Zh <- as.matrix(Z[[h]])
    if (nrow(Zh) != N_h) stop(sprintf("Z[[%d]] rows must equal site-%d record count.", h, h))
    pzh <- ncol(Zh)

    # Per-site crossproducts
    ShX  <- crossprod(Xh_list[[1]], Xh_list[[1]])       # px x px (same across outcomes)
    ShXZ <- crossprod(Xh_list[[1]], Zh)                 # px x pzh
    ShZ  <- crossprod(Zh, Zh)                           # pzh x pzh
    ShYY <- crossprod(Yh, Yh)                           # R x R
    # X'Y and Z'Y by outcome
    ShXY <- matrix(0.0, nrow = px, ncol = R)
    ShZY <- matrix(0.0, nrow = pzh, ncol = R)
    for (r in seq_len(R)) {
      ShXY[, r] <- crossprod(Xh_list[[r]], Yh[, r])
      ShZY[, r] <- crossprod(Zh,            Yh[, r])
    }

    # Outcome whitening (scale-free)
    ShX_tilde   <- kronecker(Rinve, ShX)               # (pxR) x (pxR)
    ShXZ_tilde  <- kronecker(Rinve, ShXZ)              # (pxR) x (pzhR)
    ShZ_tilde   <- kronecker(Rinve, ShZ)               # (pzhR) x (pzhR)
    ShXY_t_vec  <- as.vector(ShXY %*% Rinve)           # length pxR
    ShZY_t_vec  <- as.vector(ShZY %*% Rinve)           # length pzhR
    YtildeY     <- sum(Rinve * ShYY)                   # tr(Rinve * Y'Y)

    # Theta_h inverse (block-diag across outcomes)
    V0_inv_diag <- c(1/theta_u, rep(1/theta_vh[h], pzh - 1L))
    Theta_h_inv <- as.matrix(Matrix::bdiag(replicate(R, diag(V0_inv_diag, pzh, pzh), simplify = FALSE)))
    logdet_Theta_h <- R * sum(log(c(theta_u, rep(theta_vh[h], pzh - 1L))))

    # A_h = Theta_h^{-1} + Z~'Z~
    A_h <- Theta_h_inv + ShZ_tilde

    # log|Gamma_h| = N_h*log|Rcorr| + log|Theta_h| + log|A_h|
    lpterm1 <- lpterm1 + N_h * logdet_Rcorr +
      logdet_Theta_h + as.numeric(determinant(A_h, logarithm = TRUE)$modulus)

    # Woodbury reductions
    Ainv_tXZ <- solve(A_h, t(ShXZ_tilde))  # (pzhR) x (pxR)
    bterm1   <- bterm1 + (ShX_tilde - ShXZ_tilde %*% Ainv_tXZ)

    Ainv_ZY  <- solve(A_h, ShZY_t_vec)     # (pzhR) x 1
    bterm2   <- bterm2 + (ShXY_t_vec - as.vector(ShXZ_tilde %*% Ainv_ZY))

    lpterm2  <- lpterm2 + (YtildeY - drop(t(ShZY_t_vec) %*% Ainv_ZY))
  }

  # GLS beta (profiled over sigma^2)
  L <- chol(bterm1)
  beta_hat <- backsolve(L, forwardsolve(t(L), bterm2))

  # Quadratic for sigma^2 profiling
  qterm <- as.numeric(lpterm2 - 2 * sum(bterm2 * beta_hat) + t(beta_hat) %*% bterm1 %*% beta_hat)

  # Degrees of freedom
  NtotR <- Nsum * R
  if (NtotR != n_tot) stop("Book-keeping error: NtotR mismatch.")
  pxR   <- px * R

  if (reml) {
    s2 <- qterm / (NtotR - pxR)
    # REML penalty: log|X'Gamma^{-1}X|
    # note: bterm1 is X~'Gamma^{-1}X~ (scale-free); |X'G^{-1}X| = |bterm1| / s2^{pxR}
    logdet_XtGinvX <- as.numeric(determinant(bterm1, logarithm = TRUE)$modulus) - pxR * log(s2)
    lp <- -0.5 * ( (NtotR - pxR)*log(2*pi) + lpterm1 + (NtotR - pxR)*log(s2) + qterm/s2 + logdet_XtGinvX )
  } else {
    s2 <- qterm / NtotR
    lp <- -0.5 * ( NtotR*log(2*pi) + lpterm1 + NtotR*log(s2) + qterm/s2 )
  }

  list(
    lp = as.numeric(lp),
    b  = beta_hat,
    s2 = as.numeric(s2),
    rho = rho,
    allterms = list(
      bterm1 = bterm1, bterm2 = bterm2,
      qterm = qterm, lpterm1 = lpterm1, lpterm2 = lpterm2,
      NtotR = NtotR, pxR = pxR
    )
  )
}



# -------------------------
# Top-level fit: now supports direct X/Y/Z path with REML profile
# Set use_direct = TRUE to bypass ShXYZ and call lmm.profile3_mv_DIRECT
# -------------------------
MLMM.fit.RI <- function(Y, X, Z, id.site, weights = NULL,
                        pooled = FALSE, reml = TRUE,
                        common.s2 = TRUE,
                        ShXYZ = NULL,
                        corstr = c('exchangeable','independence'),
                        mypar.init = NULL,
                        estimate_rho = TRUE, rho_init = 0.1,
                        hessian = TRUE, verbose = TRUE,
                        use_direct = TRUE) {

  corstr <- match.arg(corstr)
  if (!is.matrix(Y)) stop("For multivariate fit, supply Y as N x R matrix.")
  R <- ncol(Y)

  # Sites / K / px
  site_fac <- as.factor(id.site)
  K  <- length(levels(site_fac))
  px <- if (is.list(X)) ncol(X[[1]]) else ncol(X)

  # Parameter init: theta_u + theta_vh (K) [+ rho]
  if (is.null(mypar.init)) {
    thetas <- rep(1, K + 1)  # theta_u=1, theta_vh[h]=1
    mypar.init <- if (estimate_rho && corstr == "exchangeable") c(thetas, rho_init) else thetas
    if (verbose) cat('Default mypar.init (theta_u + theta_v[h], [rho]) =', mypar.init, '\n')
  } else if (estimate_rho && corstr == "exchangeable" && length(mypar.init) == (K + 1)) {
    mypar.init <- c(mypar.init, rho_init)
  }

  # Bounds
  lower <- rep(1e-6, length(mypar.init))
  upper <- rep(Inf,  length(mypar.init))
  if (corstr == "exchangeable" && estimate_rho) {
    rc <- .Rcorr_inv_and_logdet(R, rho = 0)  # for bounds only
    lower[length(lower)] <- rc$lower
    upper[length(upper)] <- rc$upper - 1e-8
  }

  # Objective (switches to direct profile if requested)
  fn <- function(parameter) {
    if (corstr == "independence") {
      if (use_direct) {
        -lmm.profile3_mv_DIRECT(par = parameter[1:(K+1)],
                                Y = Y, X = X, Z = Z, id.site = id.site,
                                reml = reml, estimate_rho = FALSE, rho_fixed = 0)$lp
      } else {
        stop("independence + summaries path not wired here; use use_direct=TRUE.")
      }
    } else {
      if (use_direct) {
        -lmm.profile3_mv_DIRECT(par = parameter,
                                Y = Y, X = X, Z = Z, id.site = id.site,
                                reml = reml, estimate_rho = TRUE)$lp
      } else {
        # (legacy summaries path) – left intact if you still want it
        -lmm.profile3_mv(par = parameter,
                         Y = Y, X = X, Z = Z, id.site = id.site, ShXYZ = ShXYZ %||% peal.get.summary_mv(Y,X,Z,id.site,weights),
                         reml = reml, pooled = FALSE,
                         estimate_rho = TRUE, verbose = FALSE)$lp
      }
    }
  }

  opt <- optim(mypar.init, fn, method = "L-BFGS-B",
               hessian = hessian, lower = lower, upper = upper)

  # Recompute at optimum for estimates
  if (corstr == "independence") {
    prof <- lmm.profile3_mv_DIRECT(par = opt$par[1:(K+1)],
                                   Y = Y, X = X, Z = Z, id.site = id.site,
                                   reml = reml, estimate_rho = FALSE, rho_fixed = 0)
    rho_hat <- 0
  } else {
    prof <- lmm.profile3_mv_DIRECT(par = opt$par,
                                   Y = Y, X = X, Z = Z, id.site = id.site,
                                   reml = reml, estimate_rho = TRUE)
    rho_hat <- prof$rho
  }

  # Inference for beta: Var(vec(beta)) = s2 * (X~'Gamma^{-1}X~)^{-1}
  XtGinvX <- prof$allterms$bterm1
  Vbeta   <- tryCatch(solve(XtGinvX), error = function(e) MASS::ginv(XtGinvX)) * prof$s2
  se      <- sqrt(pmax(diag(Vbeta), 0))
  wald    <- prof$b / se
  lb      <- prof$b - 1.96 * se
  ub      <- prof$b + 1.96 * se

  # Residual covariance: Sigma_e = s2 * Rcorr(rho)
  Rcorr   <- (1 - rho_hat) * diag(R) + rho_hat * matrix(1, R, R)
  Sigma_e <- prof$s2 * Rcorr

  if (verbose) cat(ifelse(opt$convergence == 0, "Convergence reached.", "Non-convergence!"),
                   " Function evals:", opt$counts[1], "\n")

  list(
    b = prof$b, b.sd = se, wald = wald, lb = lb, ub = ub,
    theta = opt$par[1:(K+1)],     # scaled (theta_u, theta_vh[h])
    rho = rho_hat,
    Sigma_e = Sigma_e,
    s2 = prof$s2,
    opt = opt,
    res.profile = prof
  )
}
