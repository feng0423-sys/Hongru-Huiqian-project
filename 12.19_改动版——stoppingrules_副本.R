## =========================================================
##  Shared utilities
## =========================================================

symmetrize <- function(A) 0.5 * (A + t(A))

clip_eig <- function(A, eps = 0) {
  ea <- eigen(symmetrize(A), symmetric = TRUE)
  V  <- ea$vectors
  d  <- pmax(ea$values, eps)
  V %*% diag(d, length(d)) %*% t(V)
}

proj_psd <- function(A) clip_eig(A, eps = 0)

is_pd_by_chol_rank <- function(A, eps = 1e-8) {
  A_sym <- symmetrize(A)
  n <- nrow(A_sym)
  R <- tryCatch(
    suppressWarnings(chol(A_sym + diag(eps, n), pivot = TRUE)),
    error = function(e) NULL
  )
  if (is.null(R)) return(FALSE)
  attr(R, "rank") == n
}

soft_thresh <- function(A, tau) sign(A) * pmax(abs(A) - tau, 0)

soft_thresh_offdiag <- function(A, tau) {
  B <- symmetrize(A)
  idx <- matrix(TRUE, nrow(A), ncol(A)); diag(idx) <- FALSE
  B[idx] <- soft_thresh(B[idx], tau)
  B
}

safe_inv_psd <- function(A, eps = 1e-8, jitter_init = 0, jitter_max = 1e-5) {
  A_sym <- symmetrize(A)
  n <- nrow(A_sym)
  jitter <- jitter_init
  while (jitter <= jitter_max) {
    chol_try <- tryCatch(
      suppressWarnings(chol(A_sym + diag(jitter, n))),
      error = function(e) NULL
    )
    if (!is.null(chol_try)) {
      invA <- chol2inv(chol_try)
      return(symmetrize(invA))
    }
    jitter <- if (jitter == 0) eps else jitter * 10
  }
  ea <- eigen(A_sym, symmetric = TRUE)
  V  <- ea$vectors
  d  <- pmax(ea$values, eps)
  symmetrize(V %*% (t(V) / d))
}

safe_logdet_spd <- function(A, jitter_init = 0, jitter_max = 1e-5, eps = 1e-10) {
  A_sym <- symmetrize(A)
  n <- nrow(A_sym)
  jitter <- jitter_init
  while (jitter <= jitter_max) {
    R <- tryCatch(
      suppressWarnings(chol(A_sym + diag(jitter, n))),
      error = function(e) NULL
    )
    if (!is.null(R)) {
      return(2 * sum(log(diag(R))))
    }
    jitter <- if (jitter == 0) eps else jitter * 10
  }
  ea <- eigen(A_sym, symmetric = TRUE, only.values = TRUE)$values
  sum(log(pmax(ea, eps)))
}

eig_max_symmetric <- function(A) {
  eigvals <- suppressWarnings(eigen(symmetrize(A), symmetric = TRUE, only.values = TRUE)$values)
  max(eigvals)
}

## =========================================================
##  Shared simulation setup
## =========================================================

set.seed(20251025)
shared_params <- list(
  n = 1000,
  p = 100,
  q = 100,
  ratio_lambda = 1e-4,
  gamma_sparsity = 0.05,
  gamma_mag = c(1, 2),
  omega_kappa = 60,
  omega_sparsity = 0.05
)

generate_sparse_gamma <- function(p_aug, q,
                                  sparsity = 0.05,
                                  nnz_per_col = NULL,
                                  mag_range = c(1, 2),
                                  intercept_zero = TRUE,
                                  signed = TRUE) {
  if (length(mag_range) != 2) stop("mag_range must have length 2")
  lo <- min(mag_range); hi <- max(mag_range)
  gamma <- matrix(0, p_aug, q)
  idx_pool <- if (intercept_zero && p_aug > 1) 2:p_aug else seq_len(p_aug)
  if (length(idx_pool) == 0) return(gamma)
  if (is.null(nnz_per_col)) {
    nnz_per_col <- max(1, round(sparsity * length(idx_pool)))
  } else {
    nnz_per_col <- max(1, min(nnz_per_col, length(idx_pool)))
  }
  for (j in seq_len(q)) {
    nnz <- min(nnz_per_col, length(idx_pool))
    active <- sample(idx_pool, nnz)
    vals <- runif(nnz, lo, hi)
    if (signed) vals <- vals * sample(c(-1, 1), nnz, replace = TRUE)
    gamma[active, j] <- vals
  }
  if (intercept_zero && p_aug >= 1) gamma[1, ] <- 0
  gamma
}

generate_sparse_omega <- function(q,
                                  kappa_target = 50,
                                  offdiag_prob = 0.05,
                                  base_val = 0.5,
                                  ensure_pd = TRUE) {
  if (q == 1) return(matrix(1, 1, 1))
  B <- matrix(0, q, q)
  upper_mask <- matrix(runif(q^2) < offdiag_prob, q, q)
  upper_mask[lower.tri(upper_mask, diag = TRUE)] <- FALSE
  B[upper_mask] <- base_val
  B <- B + t(B)
  eig_B <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
  lambda_min <- min(eig_B)
  lambda_max <- max(eig_B)
  if (abs(lambda_min - lambda_max) < 1e-12) {
    delta <- 1
  } else {
    delta <- (lambda_max - kappa_target * lambda_min) / (kappa_target - 1)
  }
  if (!is.finite(delta) || delta <= -lambda_min + 1e-8) {
    delta <- -lambda_min + 1e-4
  }
  Omega <- B + delta * diag(q)
  if (ensure_pd) {
    eig_vals <- suppressWarnings(eigen(symmetrize(Omega), symmetric = TRUE, only.values = TRUE)$values)
    eig_vals <- pmax(eig_vals, 1e-8)
    attr(Omega, "cond_est") <- max(eig_vals) / min(eig_vals)
  }
  Omega
}

p_aug <- shared_params$p + 1
gamma_true <- generate_sparse_gamma(
  p_aug,
  shared_params$q,
  sparsity = shared_params$gamma_sparsity,
  mag_range = shared_params$gamma_mag
)
shared_Beta_true <- gamma_true[-1, , drop = FALSE]
shared_Omega_true <- generate_sparse_omega(
  shared_params$q,
  kappa_target = shared_params$omega_kappa,
  offdiag_prob = shared_params$omega_sparsity
)
shared_Sigma_true <- safe_inv_psd(shared_Omega_true)
shared_truth <- list(
  gamma = gamma_true,
  Omega = shared_Omega_true,
  cond_number = attr(shared_Omega_true, "cond_est")
)

shared_X <- matrix(rnorm(shared_params$n * shared_params$p), shared_params$n, shared_params$p)
suppressPackageStartupMessages(library(MASS))
shared_Y <- shared_X %*% shared_Beta_true + MASS::mvrnorm(
  shared_params$n,
  mu = rep(0, shared_params$q),
  Sigma = shared_Sigma_true
)
shared_Xc <- cbind(1, scale(shared_X, center = TRUE, scale = FALSE))
shared_p_aug <- ncol(shared_Xc)

Sxy_shared <- crossprod(shared_Xc, shared_Y) / shared_params$n
Syy_shared <- crossprod(shared_Y) / shared_params$n
lambda_gamma_0 <- max(abs(Sxy_shared[-1, , drop = FALSE]))
tmp_off_shared <- Syy_shared - diag(diag(Syy_shared))
lambda_Omega_0 <- max(abs(tmp_off_shared))
path_len <- 40
lambda_gamma_seq <- exp(seq(log(lambda_gamma_0),
                            log(lambda_gamma_0 * shared_params$ratio_lambda),
                            length.out = path_len))
lambda_Omega_seq <- exp(seq(log(lambda_Omega_0),
                            log(lambda_Omega_0 * shared_params$ratio_lambda),
                            length.out = path_len))
lambda_path <- data.frame(
  step = seq_len(path_len),
  lambda_gamma = lambda_gamma_seq,
  lambda_Omega = lambda_Omega_seq
)

compute_null_corner_init <- function(Sxx, Sxy, Syy, unpenalized = 1L, eps = 1e-8) {
  p_aug <- nrow(Sxx); q <- ncol(Sxy)
  gamma_init <- matrix(0, p_aug, q)
  U <- unpenalized
  Sxx_UU <- Sxx[U, U, drop = FALSE]
  Sxy_U <- Sxy[U, , drop = FALSE]
  Sxx_UU_inv <- solve(Sxx_UU)
  B <- t(Sxy_U) %*% Sxx_UU_inv %*% Sxy_U
  diag_term <- diag(Syy - B)
  diag_term <- pmax(diag_term, eps)
  Omega_diag <- 1 / diag_term
  Omega_init <- diag(Omega_diag, nrow = length(Omega_diag))
  gamma_U <- Sxx_UU_inv %*% Sxy_U %*% Omega_init
  gamma_init[U, ] <- gamma_U
  list(gamma = gamma_init, Omega = Omega_init)
}

tol_env <- suppressWarnings(as.numeric(Sys.getenv("SOLVER_TOL", unset = NA_character_)))
if (is.na(tol_env)) tol_env <- 1e-5
maxit_env <- suppressWarnings(as.integer(Sys.getenv("SOLVER_MAXIT", unset = NA_character_)))
if (is.na(maxit_env) || maxit_env <= 0) maxit_env <- 1000

ctrl <- list(
  lambda = list(gamma = lambda_gamma_seq[1], Omega = lambda_Omega_seq[1]),
  maxit = maxit_env,
  tol = tol_env,
  tol_obj = tol_env,
  tol_gm = tol_env,
  use_relative_obj = TRUE
)

prox_tol_env <- suppressWarnings(as.numeric(Sys.getenv("PROX_TOL", unset = NA_character_)))
if (is.na(prox_tol_env)) prox_tol_env <- 0.1 * tol_env
prox_max_admm_env <- suppressWarnings(as.integer(Sys.getenv("PROX_MAX_ADMM", unset = NA_character_)))
if (is.na(prox_max_admm_env) || prox_max_admm_env <= 0) prox_max_admm_env <- 1000
prox_eta_env <- suppressWarnings(as.numeric(Sys.getenv("PROX_ETA", unset = NA_character_)))
if (is.na(prox_eta_env) || prox_eta_env < 0) prox_eta_env <- 1e-6
prox_ctrl <- list(max_admm = prox_max_admm_env, tol = prox_tol_env, eta = prox_eta_env)

pn_eta_env <- suppressWarnings(as.numeric(Sys.getenv("PN_PROX_ETA", unset = NA_character_)))
if (is.na(pn_eta_env) || pn_eta_env < 0) pn_eta_env <- prox_eta_env
pn_prox_ctrl <- list(max_admm = prox_max_admm_env, tol = prox_tol_env, eta = pn_eta_env)

armijo_c <- 0.1

null_corner_init <- compute_null_corner_init(
  crossprod(shared_Xc) / shared_params$n,
  Sxy_shared,
  Syy_shared,
  unpenalized = 1L
)

shared_data <- list(
  X = shared_Xc,
  Y = shared_Y,
  lambda = ctrl$lambda,
  lambda_path = lambda_path,
  gamma_init = null_corner_init$gamma,
  Omega_init = null_corner_init$Omega,
  params = shared_params,
  truth = shared_truth
)

## =========================================================
##  Objective and derivatives
## =========================================================

.ensure_Oinv <- function(Omega, Oinv = NULL) {
  if (is.null(Oinv)) safe_inv_psd(Omega) else symmetrize(Oinv)
}
.ensure_beta <- function(gamma, Omega, Oinv = NULL, beta = NULL) {
  if (!is.null(beta)) return(beta)
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  gamma %*% Oinv_eff
}

smooth_obj <- function(gamma, Omega, Sxx, Sxy, Syy, Oinv = NULL, beta = NULL) {
  logdet <- safe_logdet_spd(Omega)
  term1 <- -logdet
  term3 <- -2 * sum(Sxy * gamma)
  term4 <- sum(Omega * Syy)
  
  if (!is.null(beta)) {
    BSB  <- t(beta) %*% Sxx %*% beta
    term2 <- sum(BSB * Omega)
  } else {
    Oinv_eff <- .ensure_Oinv(Omega, Oinv)
    M <- t(gamma) %*% Sxx %*% gamma
    term2 <- sum(Oinv_eff * M)
  }
  0.5 * (term1 + term2 + term3 + term4)
}

loss_fn <- function(gamma, Omega, Sxx, Sxy, Syy, n,
                    lambda_gamma, lambda_Omega, lambda_Omega_diag = 0,
                    Oinv = NULL, beta = NULL) {
  pen_g <- lambda_gamma * sum(abs(gamma[-1, , drop = FALSE]))
  pen_o <- lambda_Omega * (sum(abs(Omega)) - sum(abs(diag(Omega)))) +
    lambda_Omega_diag * sum(abs(diag(Omega)))
  smooth_obj(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv, beta = beta) + pen_g + pen_o
}

nll_fn <- function(gamma, Omega, Sxx, Sxy, Syy, n, Oinv = NULL, beta = NULL) {
  smooth_obj(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv, beta = beta)
}

grad_gamma <- function(gamma, Omega, Sxx, Sxy, Oinv = NULL, beta = NULL) {
  if (is.null(beta)) {
    Oinv_eff <- .ensure_Oinv(Omega, Oinv)
    Sxx %*% gamma %*% Oinv_eff - Sxy
  } else {
    Sxx %*% beta - Sxy
  }
}

grad_Omega <- function(gamma, Omega, Sxx, Sxy, Syy, Oinv = NULL, beta = NULL) {
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  if (is.null(beta)) {
    M <- t(gamma) %*% Sxx %*% gamma
    0.5 * ( Syy - Oinv_eff %*% M %*% Oinv_eff - Oinv_eff )
  } else {
    BSB <- t(beta) %*% Sxx %*% beta
    0.5 * ( Syy - BSB - Oinv_eff )
  }
}

H_blocks <- function(gamma, Omega, Sxx, Oinv = NULL, beta = NULL) {
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  beta_eff <- .ensure_beta(gamma, Omega, Oinv_eff, beta)
  
  Hgg <- Sxx
  Hgo <- - Sxx %*% beta_eff
  Hog <- - t(Hgo)
  Hoo <- 0.5 * Oinv_eff + t(beta_eff) %*% Sxx %*% beta_eff
  list(Hgg = Hgg, Hgo = Hgo, Hog = Hog, Hoo = Hoo)
}

H_matrix <- function(Hb) {
  rbind(cbind(Hb$Hgg, Hb$Hgo),
        cbind(Hb$Hog, Hb$Hoo))
}

pack_xi <- function(gamma, Omega) rbind(gamma, Omega)
unpack_xi <- function(Xi, p) list(
  gamma = Xi[1:p, , drop = FALSE],
  Omega = symmetrize(Xi[(p + 1):nrow(Xi), , drop = FALSE])
)

local_norm <- function(d_gamma, d_Omega, gamma, Omega, Sxx, Oinv = NULL, beta = NULL) {
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  Hb <- H_blocks(gamma, Omega, Sxx, Oinv = Oinv_eff, beta = beta)
  H  <- H_matrix(Hb)
  Xi_d <- pack_xi(d_gamma, d_Omega)
  val <- sum(Xi_d * (H %*% Xi_d %*% Oinv_eff))
  sqrt(max(val, 0))
}

## =========================================================
##  Prox operator for Omega
## =========================================================

prox_psd_offdiag_l1 <- function(V,
                                tau,
                                mu       = 1,
                                rho      = 1,
                                max_admm = 1000,
                                tol      = 1e-5,
                                eig_floor = 0,
                                Omega_init = NULL,
                                A_init = NULL) {
  
  Z <- symmetrize(V)
  omega_hat <- symmetrize( soft_thresh_offdiag(Z, tau / rho) )
  omega_hat <- clip_eig(omega_hat, eps = eig_floor)
  attr(omega_hat, "state") <- list(Omega = omega_hat, A = matrix(0, nrow(omega_hat), ncol(omega_hat)))
  if (is_pd_by_chol_rank(omega_hat)) {
    return(omega_hat)
  }
  
  Omega <- if (is.null(Omega_init)) omega_hat else symmetrize(Omega_init)
  A     <- if (is.null(A_init)) matrix(0, nrow(Omega), ncol(Omega)) else symmetrize(A_init)
  
  for (l in seq_len(max_admm)) {
    K <- clip_eig(Omega + mu * A, eps = eig_floor)
    if (is_pd_by_chol_rank(K)) {
      attr(K, "state") <- list(Omega = K, A = A)
      return(K)
    }
    
    T <- (K + mu * (Z - A)) / (1 + mu)
    Omega_new <- symmetrize(
      soft_thresh_offdiag(T, (mu * tau) / ((1 + mu) * rho))
    )
    
    A <- A - (K - Omega_new) / mu
    
    if (max(abs(Omega_new - Omega)) < tol &&
        max(abs(K - Omega_new))     < tol) {
      attr(Omega_new, "state") <- list(Omega = Omega_new, A = A)
      return(Omega_new)
    }
    Omega <- Omega_new
  }
  
  Omega <- symmetrize(Omega)
  attr(Omega, "state") <- list(Omega = Omega, A = A)
  Omega
}

## =========================================================
##  Stopping helpers
## =========================================================

rel_obj_change <- function(prev, curr) {
  abs(curr - prev) / (1 + abs(prev))
}

oracle_gap_ok <- function(curr, oracle, gap_tol, use_relative_gap = FALSE) {
  if (is.null(oracle) || is.null(gap_tol) || !is.finite(oracle) || !is.finite(curr)) return(FALSE)
  gap <- curr - oracle
  if (!use_relative_gap) return(gap <= gap_tol)
  gap / max(1, abs(oracle)) <= gap_tol
}

## =========================================================
##  Prox-Gradient with improved stopping
## =========================================================

pg_sparse_mvreg <- function(Xc, Y, lambda_gamma, lambda_Omega,
                            gamma_init = NULL, Omega_init = NULL,
                            max_iter = 200, tol = 1e-4, L0 = 1,
                            track_loss = TRUE,
                            ctrl = NULL, warm_state = NULL,
                            oracle_value = NULL, gap_tol = NULL, use_relative_gap = FALSE,
                            force_oracle_only = FALSE) {
  
  if (!is.null(ctrl)) {
    if (!is.null(ctrl$lambda$gamma)) lambda_gamma <- ctrl$lambda$gamma
    if (!is.null(ctrl$lambda$Omega)) lambda_Omega <- ctrl$lambda$Omega
    if (!is.null(ctrl$maxit)) max_iter <- ctrl$maxit
    if (!is.null(ctrl$tol)) tol <- ctrl$tol
  }
  
  tol_obj <- if (!is.null(ctrl$tol_obj)) ctrl$tol_obj else tol
  tol_gm  <- if (!is.null(ctrl$tol_gm))  ctrl$tol_gm  else tol
  use_rel_obj <- isTRUE(ctrl$use_relative_obj)
  
  n <- nrow(Xc); p <- ncol(Xc); q <- ncol(Y)
  Sxx <- crossprod(Xc) / n
  Sxy <- crossprod(Xc, Y) / n
  Syy <- crossprod(Y) / n
  
  gamma <- if (is.null(gamma_init)) matrix(0, p, q) else gamma_init
  Omega <- if (is.null(Omega_init)) diag(q)          else symmetrize(Omega_init)
  Oinv  <- safe_inv_psd(Omega)
  
  t0 <- Sys.time()
  loss_curr <- loss_fn(gamma, Omega, Sxx, Sxy, Syy, n,
                       lambda_gamma, lambda_Omega, Oinv = Oinv)
  
  if (track_loss) {
    loss_hist <- numeric(max_iter + 1)
    time_hist <- numeric(max_iter + 1)
    loss_hist[1] <- loss_curr
    time_hist[1] <- 0
  }
  
  L_prev <- L0
  converged <- FALSE
  stop_reason <- NA_character_
  backtrack_stalls <- 0L
  pg_psd_state <- if (!is.null(warm_state) && !is.null(warm_state$psd_state)) warm_state$psd_state else NULL
  
  for (m in seq_len(max_iter)) {
    Lm <- L_prev
    loss_prev <- loss_curr
    accepted_step <- FALSE
    
    repeat {
      Gg_m <- grad_gamma(gamma, Omega, Sxx, Sxy, Oinv = Oinv)
      Go_m <- grad_Omega(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv)
      
      Vg_m <- gamma - Gg_m / Lm
      Vo_m <- symmetrize(Omega - Go_m / Lm)
      
      s_gamma_m <- soft_thresh(Vg_m, lambda_gamma / Lm)
      s_gamma_m[1, ] <- Vg_m[1, ]
      
      s_Omega_m <- prox_psd_offdiag_l1(
        Vo_m,
        tau = lambda_Omega / Lm,
        Omega_init = if (is.null(pg_psd_state)) Omega else pg_psd_state$Omega,
        A_init = if (is.null(pg_psd_state)) NULL else pg_psd_state$A,
        max_admm = prox_ctrl$max_admm,
        tol = prox_ctrl$tol,
        eig_floor = prox_ctrl$eta
      )
      pg_psd_state <- attr(s_Omega_m, "state")
      
      d_gamma_m <- s_gamma_m - gamma
      d_Omega_m <- s_Omega_m - Omega
      
      lambda_m <- local_norm(d_gamma_m, d_Omega_m, gamma, Omega, Sxx, Oinv = Oinv)
      beta_m   <- sqrt(Lm) * sqrt(sum(d_gamma_m^2) + sum(d_Omega_m^2))
      
      if (lambda_m^2 / max(beta_m^2, .Machine$double.eps) + lambda_m > 1) {
        Lm <- 2 * Lm
        next
      }
      
      alpha_m <- (beta_m^2) / (lambda_m * (lambda_m + beta_m^2))
      gamma_old <- gamma
      Omega_old <- Omega
      
      gamma_candidate <- gamma_old + alpha_m * d_gamma_m
      Omega_candidate <- symmetrize(Omega_old + alpha_m * d_Omega_m)
      
      Oinv_try <- tryCatch(safe_inv_psd(Omega_candidate), error = function(e) NULL)
      if (is.null(Oinv_try)) {
        backtrack_stalls <- backtrack_stalls + 1L
        L_prev <- Lm
        break
      }
      
      loss_try <- loss_fn(gamma_candidate, Omega_candidate, Sxx, Sxy, Syy, n,
                          lambda_gamma, lambda_Omega, Oinv = Oinv_try)
      if (!is.finite(loss_try)) {
        backtrack_stalls <- backtrack_stalls + 1L
        L_prev <- Lm
        break
      }
      
      gamma <- gamma_candidate
      Omega <- Omega_candidate
      Oinv  <- Oinv_try
      loss_curr <- loss_try
      L_prev <- Lm
      accepted_step <- TRUE
      
      ## Oracle gap stopping if requested
      if (oracle_gap_ok(loss_curr, oracle_value, gap_tol, use_relative_gap)) {
        converged <- TRUE
        stop_reason <- "oracle gap reached"
        break
      }
      if (force_oracle_only) {
        # keep iterating until oracle gap or max_iter
        break
      }
      
      ## Objective stabilization stopping
      if (use_rel_obj) {
        if (rel_obj_change(loss_prev, loss_curr) <= tol_obj) {
          converged <- TRUE
          stop_reason <- "objective stabilized"
          break
        }
      } else {
        if (abs(loss_curr - loss_prev) <= tol_obj) {
          converged <- TRUE
          stop_reason <- "objective stabilized"
          break
        }
      }
      
      break
    }
    
    if (track_loss) {
      loss_hist[m + 1] <- loss_curr
      time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    
    if (converged) break
    if (!accepted_step) next
  }
  
  elapsed_time_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out <- list(gamma = gamma, Omega = Omega, iters = m,
              converged = converged,
              stop_reason = ifelse(is.na(stop_reason), "max_iter reached", stop_reason),
              elapsed_time_sec = elapsed_time_sec,
              Sxx = Sxx, Sxy = Sxy, Syy = Syy, n = n,
              backtrack_stalls = backtrack_stalls,
              psd_state = pg_psd_state)
  if (track_loss) {
    out$loss_hist <- loss_hist[1:(m + 1)]
    out$time_hist <- time_hist[1:(m + 1)]
  }
  out
}

## =========================================================
##  Prox-Newton subproblem ADMM
## =========================================================

pn_subproblem_admm <- function(gamma, Omega, Sxx, Sxy, Syy,
                               lambda_gamma, lambda_Omega,
                               Oinv = NULL, beta = NULL,
                               rho = 1, mu = 1,
                               max_admm = 150, admm_tol = 1e-4,
                               eig_floor = 1e-12,
                               state = NULL) {
  p <- nrow(gamma); q <- ncol(gamma)
  
  beta_eff <- .ensure_beta(gamma, Omega, Oinv, beta)
  Oinv_eff <- .ensure_Oinv(Omega, Oinv)
  
  Gg <- grad_gamma(gamma, Omega, Sxx, Sxy, Oinv = Oinv_eff, beta = beta_eff)
  Go <- grad_Omega(gamma, Omega, Sxx, Sxy, Syy, Oinv = Oinv_eff, beta = beta_eff)
  P  <- pack_xi(Gg, Go)
  
  Hb <- H_blocks(gamma, Omega, Sxx, Oinv = Oinv_eff, beta = beta_eff)
  H  <- symmetrize(H_matrix(Hb))
  
  Xi <- pack_xi(gamma, Omega)
  Gm <- P - H %*% Xi %*% Oinv_eff
  
  eigH <- eigen(symmetrize(H), symmetric = TRUE)
  UH   <- eigH$vectors
  lamH <- pmax(eigH$values, eig_floor)
  
  eigS <- eigen(symmetrize(Oinv_eff), symmetric = TRUE)
  US   <- eigS$vectors
  lamS <- pmax(eigS$values, eig_floor)
  
  Xi_var <- Xi
  Z_var  <- Xi
  Gam    <- matrix(0, p + q, q)
  if (!is.null(state)) {
    if (!is.null(state$Xi)) Xi_var <- state$Xi
    if (!is.null(state$Z))  Z_var  <- state$Z
    if (!is.null(state$Gam)) Gam   <- state$Gam
  }
  prox_state <- if (!is.null(state)) state$psd_state else NULL
  
  Den <- outer(lamH, lamS, "*") + rho
  
  for (k in 1:max_admm) {
    R  <- rho * Z_var - Gam - Gm
    Rt <- crossprod(UH, R) %*% US
    Y  <- Rt / Den
    Xi_var <- UH %*% Y %*% t(US)
    
    Xi_parts <- unpack_xi(Xi_var + Gam / rho, p)
    Zg <- soft_thresh(Xi_parts$gamma, lambda_gamma / rho)
    Zg[1, ] <- Xi_parts$gamma[1, ]
    
    Zo <- prox_psd_offdiag_l1(
      Xi_parts$Omega,
      tau = lambda_Omega / rho,
      mu = mu,
      rho = rho,
      Omega_init = if (is.null(prox_state)) Omega else prox_state$Omega,
      A_init = if (is.null(prox_state)) NULL else prox_state$A,
      max_admm = pn_prox_ctrl$max_admm,
      tol = pn_prox_ctrl$tol,
      eig_floor = pn_prox_ctrl$eta
    )
    prox_state <- attr(Zo, "state")
    Z_new <- pack_xi(Zg, Zo)
    
    Gam <- Gam + rho * (Xi_var - Z_new)
    
    r_norm <- sqrt(sum((Xi_var - Z_new)^2))
    s_norm <- rho * sqrt(sum((Z_new - Z_var)^2))
    Z_var  <- Z_new
    
    scale_pri <- max(sqrt(length(Xi_var)), sqrt(sum(Xi_var^2)), sqrt(sum(Z_var^2)))
    scale_dual <- max(sqrt(length(Gam)), sqrt(sum(Gam^2)))
    eps_pri <- admm_tol * scale_pri
    eps_dual <- admm_tol * scale_dual
    
    if (r_norm <= eps_pri && s_norm <= eps_dual) break
  }
  
  out_sub <- unpack_xi(Z_var, p)
  attr(out_sub$Omega, "state") <- list(
    psd_state = prox_state,
    Xi = Xi_var,
    Z = Z_var,
    Gam = Gam
  )
  out_sub
}

## =========================================================
##  Prox-Newton outer loop with optional oracle stopping
## =========================================================

pn_sparse_mvreg <- function(Xc, Y, lambda_gamma, lambda_Omega,
                            gamma_init = NULL, Omega_init = NULL,
                            max_iter = 100, tol = 1e-4, track_loss = TRUE,
                            ctrl = NULL, warm_state = NULL,
                            oracle_value = NULL, gap_tol = NULL, use_relative_gap = FALSE,
                            force_oracle_only = FALSE) {
  
  if (!is.null(ctrl)) {
    if (!is.null(ctrl$lambda$gamma)) lambda_gamma <- ctrl$lambda$gamma
    if (!is.null(ctrl$lambda$Omega)) lambda_Omega <- ctrl$lambda$Omega
    if (!is.null(ctrl$maxit)) max_iter <- ctrl$maxit
    if (!is.null(ctrl$tol)) tol <- ctrl$tol
  }
  
  tol_obj <- if (!is.null(ctrl$tol_obj)) ctrl$tol_obj else tol
  use_rel_obj <- isTRUE(ctrl$use_relative_obj)
  
  n <- nrow(Xc); p <- ncol(Xc); q <- ncol(Y)
  Sxx <- crossprod(Xc) / n
  Sxy <- crossprod(Xc, Y) / n
  Syy <- crossprod(Y) / n
  
  gamma <- if (is.null(gamma_init)) matrix(0, p, q) else gamma_init
  Omega <- if (is.null(Omega_init)) diag(q)          else symmetrize(Omega_init)
  Oinv  <- safe_inv_psd(Omega)
  
  t0 <- Sys.time()
  loss_curr <- loss_fn(gamma, Omega, Sxx, Sxy, Syy, n,
                       lambda_gamma, lambda_Omega, Oinv = Oinv)
  
  if (track_loss) {
    loss_hist <- numeric(max_iter + 1)
    time_hist <- numeric(max_iter + 1)
    loss_hist[1] <- loss_curr
    time_hist[1] <- 0
  }
  
  converged <- FALSE
  stop_reason <- NA_character_
  lambda_hist <- numeric(max_iter)
  pn_psd_state <- if (!is.null(warm_state) && !is.null(warm_state$psd_state)) warm_state$psd_state else NULL
  
  for (m in seq_len(max_iter)) {
    loss_prev <- loss_curr
    beta <- gamma %*% Oinv
    
    sub <- pn_subproblem_admm(gamma, Omega, Sxx, Sxy, Syy,
                              lambda_gamma, lambda_Omega,
                              Oinv = Oinv, beta = beta,
                              rho = 1, mu = 1,
                              max_admm = 150, admm_tol = 1e-4,
                              state = pn_psd_state)
    pn_psd_state <- attr(sub$Omega, "state")
    
    d_gamma <- sub$gamma - gamma
    d_Omega <- sub$Omega - Omega
    
    lam <- local_norm(d_gamma, d_Omega, gamma, Omega, Sxx, Oinv = Oinv)
    lambda_hist[m] <- lam
    
    if (force_oracle_only) {
      # skip step-based stop; rely on oracle gap or max_iter
      alpha <- (1 + lam)^(-1)
      gamma_new <- gamma + alpha * d_gamma
      Omega_new <- symmetrize(Omega + alpha * d_Omega)
      Oinv_new <- safe_inv_psd(Omega_new)
      loss_new <- loss_fn(gamma_new, Omega_new, Sxx, Sxy, Syy, n,
                          lambda_gamma, lambda_Omega, Oinv = Oinv_new)
      gamma <- gamma_new; Omega <- Omega_new; Oinv <- Oinv_new; loss_curr <- loss_new
      if (oracle_gap_ok(loss_curr, oracle_value, gap_tol, use_relative_gap)) {
        converged <- TRUE; stop_reason <- "oracle gap reached"
      }
    } else {
      alpha <- (1 + lam)^(-1)
      gamma_new <- gamma + alpha * d_gamma
      Omega_new <- symmetrize(Omega + alpha * d_Omega)
      Oinv_new <- safe_inv_psd(Omega_new)
      
      loss_new <- loss_fn(gamma_new, Omega_new, Sxx, Sxy, Syy, n,
                          lambda_gamma, lambda_Omega, Oinv = Oinv_new)
      
      gamma <- gamma_new
      Omega <- Omega_new
      Oinv  <- Oinv_new
      loss_curr <- loss_new
      
      if (oracle_gap_ok(loss_curr, oracle_value, gap_tol, use_relative_gap)) {
        converged <- TRUE
        stop_reason <- "oracle gap reached"
      } else if (use_rel_obj) {
        if (rel_obj_change(loss_prev, loss_curr) <= tol_obj) {
          converged <- TRUE
          stop_reason <- "objective stabilized"
        }
      } else {
        if (abs(loss_curr - loss_prev) <= tol_obj) {
          converged <- TRUE
          stop_reason <- "objective stabilized"
        }
      }
    }
    
    if (track_loss) {
      loss_hist[m + 1] <- loss_curr
      time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    
    if (converged) break
  }
  
  elapsed_time_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out <- list(gamma = gamma, Omega = Omega, iters = m,
              converged = converged,
              stop_reason = ifelse(is.na(stop_reason), "max_iter reached", stop_reason),
              elapsed_time_sec = elapsed_time_sec,
              Sxx = Sxx, Sxy = Sxy, Syy = Syy, n = n,
              lambda_hist = lambda_hist[seq_len(m)],
              psd_state = pn_psd_state)
  if (track_loss) {
    out$loss_hist <- loss_hist[1:(m + 1)]
    out$time_hist <- time_hist[1:(m + 1)]
  }
  out
}

## =========================================================
##  Lambda-path runners with optional oracle stopping
## =========================================================

run_pg_path <- function(data_list, ctrl, ref_path = NULL,
                        oracle_vec = NULL, gap_tol = NULL, use_relative_gap = FALSE,
                        force_oracle_only = FALSE) {
  lambda_path <- data_list$lambda_path
  gamma_init <- data_list$gamma_init
  Omega_init <- data_list$Omega_init
  warm_state <- NULL
  results <- vector("list", nrow(lambda_path))
  
  for (i in seq_len(nrow(lambda_path))) {
    ctrl_i <- ctrl
    ctrl_i$lambda$gamma <- lambda_path$lambda_gamma[i]
    ctrl_i$lambda$Omega <- lambda_path$lambda_Omega[i]
    
    gamma_init_step <- gamma_init
    Omega_init_step <- Omega_init
    warm_state_step <- warm_state
    
    if (i == 1 && !is.null(ref_path) && length(ref_path) >= 1) {
      seed_fit <- ref_path[[1]]
      if (!is.null(seed_fit$gamma)) gamma_init_step <- seed_fit$gamma
      if (!is.null(seed_fit$Omega)) Omega_init_step <- seed_fit$Omega
      if (!is.null(seed_fit$psd_state)) warm_state_step <- list(psd_state = seed_fit$psd_state)
    }
    
    oracle_value <- if (!is.null(oracle_vec)) oracle_vec[i] else NULL
    
    fit <- pg_sparse_mvreg(data_list$X, data_list$Y,
                           lambda_path$lambda_gamma[i],
                           lambda_path$lambda_Omega[i],
                           gamma_init = gamma_init_step,
                           Omega_init = Omega_init_step,
                           max_iter = ctrl$maxit, tol = ctrl$tol, L0 = 1,
                           track_loss = TRUE,
                           ctrl = ctrl_i, warm_state = warm_state_step,
                           oracle_value = oracle_value, gap_tol = gap_tol,
                           use_relative_gap = use_relative_gap,
                           force_oracle_only = force_oracle_only)
    
    fit$lambda_gamma <- lambda_path$lambda_gamma[i]
    fit$lambda_Omega <- lambda_path$lambda_Omega[i]
    fit$path_index <- lambda_path$step[i]
    results[[i]] <- fit
    
    gamma_init <- fit$gamma
    Omega_init <- fit$Omega
    warm_state <- list(psd_state = fit$psd_state)
  }
  results
}

run_pn_path <- function(data_list, ctrl,
                        oracle_vec = NULL, gap_tol = NULL, use_relative_gap = FALSE,
                        force_oracle_only = FALSE) {
  lambda_path <- data_list$lambda_path
  gamma_init <- data_list$gamma_init
  Omega_init <- data_list$Omega_init
  warm_state <- NULL
  results <- vector("list", nrow(lambda_path))
  
  for (i in seq_len(nrow(lambda_path))) {
    ctrl_i <- ctrl
    ctrl_i$lambda$gamma <- lambda_path$lambda_gamma[i]
    ctrl_i$lambda$Omega <- lambda_path$lambda_Omega[i]
    
    oracle_value <- if (!is.null(oracle_vec)) oracle_vec[i] else NULL
    
    fit <- pn_sparse_mvreg(data_list$X, data_list$Y,
                           lambda_path$lambda_gamma[i],
                           lambda_path$lambda_Omega[i],
                           gamma_init = gamma_init,
                           Omega_init = Omega_init,
                           max_iter = ctrl$maxit, tol = ctrl$tol, track_loss = TRUE,
                           ctrl = ctrl_i, warm_state = warm_state,
                           oracle_value = oracle_value, gap_tol = gap_tol,
                           use_relative_gap = use_relative_gap,
                           force_oracle_only = force_oracle_only)
    
    fit$lambda_gamma <- lambda_path$lambda_gamma[i]
    fit$lambda_Omega <- lambda_path$lambda_Omega[i]
    fit$path_index <- lambda_path$step[i]
    results[[i]] <- fit
    
    gamma_init <- fit$gamma
    Omega_init <- fit$Omega
    warm_state <- list(psd_state = fit$psd_state)
  }
  results
}

## =========================================================
##  Oracle path utility for fair comparisons
## =========================================================

compute_oracle_path <- function(data_list, ctrl_oracle) {
  pn_oracle <- run_pn_path(data_list, ctrl_oracle)
  oracle_vals <- vapply(pn_oracle, function(fit) tail(fit$loss_hist, 1), numeric(1))
  list(pn_path = pn_oracle, oracle_vals = oracle_vals)
}

## =========================================================
##  Benchmark wrapper
## =========================================================

benchmark_pg_vs_pn <- function(data_list, ctrl,
                               use_oracle = FALSE,
                               gap_tol = NULL,
                               use_relative_gap = FALSE,
                               ctrl_oracle = NULL,
                               force_oracle_only = FALSE) {
  
  oracle_vals <- NULL
  pn_oracle_path <- NULL
  if (use_oracle) {
    if (is.null(ctrl_oracle)) stop("ctrl_oracle must be provided when use_oracle is TRUE")
    oracle_out <- compute_oracle_path(data_list, ctrl_oracle)
    pn_oracle_path <- oracle_out$pn_path
    oracle_vals <- oracle_out$oracle_vals
  }
  
  pn_path <- run_pn_path(data_list, ctrl, oracle_vec = oracle_vals, gap_tol = gap_tol,
                         use_relative_gap = use_relative_gap,
                         force_oracle_only = force_oracle_only)
  pg_path <- run_pg_path(data_list, ctrl, oracle_vec = oracle_vals, gap_tol = gap_tol,
                         use_relative_gap = use_relative_gap,
                         force_oracle_only = force_oracle_only)
  
  pn_fit <- pn_path[[length(pn_path)]]
  pg_fit <- pg_path[[length(pg_path)]]
  
  summarize_path <- function(res_list, method_name) {
    do.call(rbind, lapply(res_list, function(fit) {
      data.frame(method = method_name,
                 path_index = fit$path_index,
                 lambda_gamma = fit$lambda_gamma,
                 lambda_Omega = fit$lambda_Omega,
                 elapsed_sec = fit$elapsed_time_sec,
                 iters = fit$iters,
                 converged = fit$converged,
                 stop_reason = fit$stop_reason,
                 final_loss = tail(fit$loss_hist, 1),
                 stringsAsFactors = FALSE)
    }))
  }
  
  path_summary <- rbind(
    summarize_path(pn_path, "Prox-Newton"),
    summarize_path(pg_path, "Prox-Gradient")
  )
  
  list(pg = pg_fit,
       pn = pn_fit,
       pg_path = pg_path,
       pn_path = pn_path,
       oracle_path = pn_oracle_path,
       oracle_vals = oracle_vals,
       path_summary = path_summary)
}

## =========================================================
##  Example run and plots
## =========================================================

suppressPackageStartupMessages(library(ggplot2))

## Example 1: run with stabilization stopping only
res <- benchmark_pg_vs_pn(shared_data, ctrl, use_oracle = FALSE)

## Example 2: fair comparison via oracle gap
## Choose a tight oracle control and then stop both methods by objective gap
ctrl_oracle <- ctrl
ctrl_oracle$tol <- 1e-10
ctrl_oracle$tol_obj <- 1e-12
ctrl_oracle$tol_gm <- 1e-10
ctrl_oracle$maxit <- 5000

## For fair stopping, pick gap_tol values like 1e-2, 1e-4, 1e-6 in a loop
## Here is one run at gap_tol = 1e-4
# res_gap <- benchmark_pg_vs_pn(shared_data, ctrl, use_oracle = TRUE,
#                              gap_tol = 1e-4, use_relative_gap = TRUE,
#                              ctrl_oracle = ctrl_oracle)

## Plot final loss vs lambda index for the stabilization stopping run
loss_plot_df <- res$path_summary

tol_tag <- Sys.getenv("SOLVER_TOL", unset = "")
suffix <- if (nzchar(tol_tag)) paste0("_tol", tol_tag) else ""

final_loss_plot <- ggplot(loss_plot_df,
                          aes(x = log(lambda_gamma), y = final_loss,
                              colour = method, shape = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_shape_manual(values = c("Prox-Newton" = 16, "Prox-Gradient" = 4)) +
  labs(x = expression(log(lambda)),
       y = "Training Loss",
       colour = NULL,
       shape = NULL) +
  theme_minimal(base_size = 13)

iter_plot <- ggplot(loss_plot_df,
                    aes(x = log(lambda_gamma), y = iters,
                        colour = method, shape = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  scale_shape_manual(values = c("Prox-Newton" = 16, "Prox-Gradient" = 4)) +
  labs(x = expression(log(lambda[gamma])),
       y = "number of iterations",
       colour = NULL,
       shape = NULL) +
  theme_minimal(base_size = 13)

ggsave(paste0("pn_pg_final_loss_vs_lambda", suffix, ".png"), final_loss_plot, width = 9, height = 4.5, dpi = 300)
ggsave(paste0("pn_pg_iters_vs_lambda", suffix, ".png"), iter_plot, width = 9, height = 4.5, dpi = 300)

print(final_loss_plot)
print(iter_plot)

total_time_by_method <- aggregate(elapsed_sec ~ method, data = res$path_summary, sum)
print(total_time_by_method)

saveRDS(res, file = "res_latest.rds")
saveRDS(list(
  pn_path = res$pn_path,
  pg_path = res$pg_path,
  lambda_path = shared_data$lambda_path
), file = "lambda_path_estimators.rds")
