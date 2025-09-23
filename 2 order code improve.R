# ================================================================
# Sparse multivariate regression (gamma, Omega)
# First-order Prox-Gradient vs Second-order Proximal-Newton (inexact)
#
# Goal: implement a proximal-Newton outer loop (per the notes) with
#  - quadratic model using block-diagonal Hessian
#  - gamma-block solved by a short inner prox-gradient on the Newton quadratic
#  - Omega-block solved via an ADMM-based PSD + L1 proximal around the
#    Newton center (identity-metric prox, inexact but effective)
#  - self-concordant step size alpha = 1/(1 + lambda_local)
# and compare convergence speed vs our first-order prox-gradient under a
# fixed wall-clock time budget.
#
# Reference: "Notes on a Proximal–Newton Framework for Matrix-Valued
# Gaussian Likelihoods" (July 2025): Algorithm 1 and the ADMM inner loop.
# ================================================================

############################  Helper functions  #############################

symmetrize <- function(A) 0.5 * (A + t(A))                   # Symmetrize a square matrix

is_psd_chol <- function(A, tol = 1e-8) {                     # PSD check via pivoted Cholesky
  p <- tryCatch(chol(A, pivot = TRUE, tol = tol), error = function(e) NULL)
  !is.null(p) && attr(p, "rank") == nrow(A)
}

clip_eig <- function(A, eps = 0) {                           # Eigenvalue clipping to enforce PSD
  ea <- eigen(symmetrize(A), symmetric = TRUE)
  V  <- ea$vectors
  d  <- pmax(ea$values, eps)
  V %*% diag(d, length(d)) %*% t(V)
}

proj_psd <- function(A) clip_eig(A, eps = 0)                 # PSD projection

safe_logdet <- function(A, eps0 = 1e-6) {                    # Stable log|A| (adds eps0*I)
  (determinant(A + diag(eps0, nrow(A)), logarithm = TRUE)$modulus)
}

safe_inv_psd <- function(A, eps = 1e-6) {                    # Stable inverse for PSD via eig floor
  ea <- eigen(symmetrize(A), symmetric = TRUE)
  V  <- ea$vectors
  d  <- pmax(ea$values, eps)
  V %*% diag(1 / d, length(d)) %*% t(V)
}

soft_thresh <- function(A, tau) sign(A) * pmax(abs(A) - tau, 0)   # L1 proximal (elementwise)

eig_max_symmetric <- function(A) {                                # Largest eigenvalue of symmetric A
  eigvals <- suppressWarnings(eigen(symmetrize(A), symmetric = TRUE, only.values = TRUE)$values)
  max(eigvals)
}

####################  ADMM proximal operator for Omega  #####################
# Solves:  min_Omega  0.5 * ||Omega - V||_F^2 + lambda * ||offdiag(Omega)||_1
#          s.t. Omega is PSD. (Identity-metric prox around V.)
# This is used as the inexact prox for the Omega block around the Newton center.

admm_psd <- function(V, lambda, mu = 1, max_admm = 1000, tol = 1e-4, eps_psd = 1e-8) {
  p      <- nrow(V)
  Omega  <- proj_psd(V)                         # start at PSD projection of V
  A      <- matrix(0, p, p)                     # scaled dual variable
  
  for (l in 1:max_admm) {
    K          <- proj_psd(Omega + mu * A)      # PSD projection of auxiliary variable
    tmp        <- mu * symmetrize(V) - mu * A + K
    Omega_new  <- symmetrize(soft_thresh(tmp, mu * lambda) / (mu + 1))  # L1 prox + averaging
    A          <- A - (K - Omega_new) / mu       # dual ascent (scaled form)
    
    if (max(abs(Omega_new - Omega)) < tol &&
        max(abs(K - Omega_new))      < tol &&
        is_psd_chol(Omega_new, tol = eps_psd))
      return(Omega_new)
    
    Omega <- Omega_new
  }
  K                                           # return last PSD projection if not converged
}

#########################  Loss, gradient, NLL  #############################
# Moments use averages: Sxx = X^T X / n, Sxy = X^T Y / n, Syy = Y^T Y / n.

loss_fn <- function(gamma, Omega, Sxx, Sxy, Syy, n, lambda_gamma, lambda_Omega, eps = 1e-6) {
  logdet <- safe_logdet(Omega, eps)
  resid  <- Syy - t(gamma) %*% Sxy - t(Sxy) %*% gamma + t(gamma) %*% Sxx %*% gamma
  pen_g  <- lambda_gamma * sum(abs(gamma[-1, , drop = FALSE]))       # no penalty on intercept row
  off_d  <- abs(Omega) - abs(diag(diag(Omega)))
  pen_o  <- lambda_Omega * sum(off_d)                                # off-diagonal L1 on Omega
  0.5 * (-logdet + sum(Omega * resid)) + pen_g + pen_o
}

nll_fn <- function(gamma, Omega, Sxx, Sxy, Syy, n, eps = 1e-6) {
  logdet <- safe_logdet(Omega, eps)
  resid  <- Syy - t(gamma) %*% Sxy - t(Sxy) %*% gamma + t(gamma) %*% Sxx %*% gamma
  0.5 * (-logdet + sum(Omega * resid))
}

grad_gamma <- function(gamma, Omega, Sxx, Sxy) {
  Omega_inv <- safe_inv_psd(Omega)
  (Sxx %*% gamma - Sxy) %*% Omega_inv
}

grad_Omega <- function(gamma, Omega, Sxx, Sxy, Syy, eps = 1e-6) {
  Omega_inv <- safe_inv_psd(Omega, eps)
  Sigma_res <- Syy - t(gamma) %*% Sxy - t(Sxy) %*% gamma + t(gamma) %*% Sxx %*% gamma
  0.5 * (Sigma_res - Omega_inv)
}

########################  First-order: Prox-Gradient  #######################
# Includes optional time_budget (seconds). Stops when time exceeded.

pg_sparse_mvreg <- function(Xc, Y, lambda_gamma, lambda_Omega,
                            gamma_init = NULL, Omega_init = NULL,
                            max_iter = 200, tol = 1e-4, L0 = 1,
                            track_loss = TRUE, time_budget = Inf,
                            X_test = NULL, Y_test = NULL) {
  n   <- nrow(Xc); p <- ncol(Xc); q <- ncol(Y)
  Sxx <- crossprod(Xc) / n
  Sxy <- crossprod(Xc, Y) / n
  Syy <- crossprod(Y) / n
  
  gamma <- if (is.null(gamma_init)) matrix(0, p, q) else gamma_init
  Omega <- if (is.null(Omega_init)) diag(q)          else Omega_init
  L     <- L0
  
  if (track_loss) {
    loss_hist <- numeric(max_iter + 1)
    time_hist <- numeric(max_iter + 1)
    t0        <- Sys.time()
    loss_hist[1] <- loss_fn(gamma, Omega, Sxx, Sxy, Syy, n, lambda_gamma, lambda_Omega)
    time_hist[1] <- 0
  }
  
  for (m in 1:max_iter) {
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > time_budget) {
      m <- m - 1; break
    }
    
    Gg <- grad_gamma(gamma, Omega, Sxx, Sxy)
    Go <- grad_Omega(gamma, Omega, Sxx, Sxy, Syy)
    
    repeat { # backtracking
      Vg <- gamma - Gg / L
      Vo <- Omega - Go  / L
      
      gamma_new <- soft_thresh(Vg, lambda_gamma / L)
      gamma_new[1, ] <- Vg[1, ]        # keep intercept unpenalised
      
      Vo_sym    <- symmetrize(Vo)
      Omega_tmp <- soft_thresh(Vo_sym, lambda_Omega / L)
      Omega_new <- if (is_psd_chol(Omega_tmp)) Omega_tmp else admm_psd(Vo_sym, lambda_Omega / L)
      
      xi_old <- c(gamma, Omega); xi_new <- c(gamma_new, Omega_new)
      diff   <- xi_new - xi_old
      
      lhs <- loss_fn(gamma_new, Omega_new, Sxx, Sxy, Syy, n, lambda_gamma, lambda_Omega)
      rhs <- loss_fn(gamma, Omega, Sxx, Sxy, Syy, n, lambda_gamma, lambda_Omega) +
        sum(c(Gg, Go) * diff) + (L / 2) * sum(diff^2)
      if (lhs <= rhs) break
      L <- L * 2
    }
    
    if (track_loss) {
      loss_hist[m + 1] <- lhs
      time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    
    gamma <- gamma_new; Omega <- symmetrize(Omega_new)
    
    if (sum(abs(diff)) < tol) break
  }
  
  out <- list(gamma = gamma, Omega = Omega, iters = m, Sxx = Sxx, Sxy = Sxy, Syy = Syy, n = n)
  if (track_loss) { out$loss_hist <- loss_hist[1:(m + 1)]; out$time_hist <- time_hist[1:(m + 1)] }
  out
}

#######################  Second-order: Proximal-Newton  #####################
# Block-diagonal Hessian approximation H = diag(H_gamma, H_Omega)
# - gamma-block Hessian: H_gamma(Delta) = Sxx * Delta * Omega^{-1}
# - Omega-block Hessian:  H_Omega(Delta) = 0.5 * Omega^{-1} * Delta * Omega^{-1}
# Subproblem is solved inexactly:
#   * gamma: a few inner prox-gradient steps on the quadratic model
#   * Omega: ADMM-based identity-metric proximal around Newton center
# Step size: self-concordant alpha = 1 / (1 + lambda_local)

solve_gamma_subproblem <- function(gamma, Omega_inv, Sxx, Gg, lambda_gamma,
                                   inner_max = 10, inner_tol = 1e-6) {
  # Quadratic model: 0.5 * tr( Omega_inv * (Gamma - gamma)^T * Sxx * (Gamma - gamma) ) + <Gg, (Gamma - gamma)>
  # Prox step: soft-threshold on Gamma with step 1/Lg, where Lg = lambda_max(Sxx) * lambda_max(Omega_inv)
  Lg <- max(1e-8, eig_max_symmetric(Sxx)) * max(1e-8, eig_max_symmetric(Omega_inv))
  Gamma <- gamma
  for (t in 1:inner_max) {
    Grad_q <- Sxx %*% (Gamma - gamma) %*% Omega_inv + Gg      # gradient of the quadratic model
    V      <- Gamma - Grad_q / Lg                              # gradient step in quadratic metric
    Gamma_new <- soft_thresh(V, lambda_gamma / Lg)             # L1 proximal (elementwise)
    Gamma_new[1, ] <- V[1, ]                                   # do not penalise intercept row
    if (sum(abs(Gamma_new - Gamma)) < inner_tol) { Gamma <- Gamma_new; break }
    Gamma <- Gamma_new
  }
  Gamma
}

pn_sparse_mvreg <- function(Xc, Y, lambda_gamma, lambda_Omega,
                            gamma_init = NULL, Omega_init = NULL,
                            max_iter = 100, tol = 1e-4, track_loss = TRUE,
                            time_budget = Inf, backtrack = TRUE) {
  n   <- nrow(Xc); p <- ncol(Xc); q <- ncol(Y)
  Sxx <- crossprod(Xc) / n
  Sxy <- crossprod(Xc, Y) / n
  Syy <- crossprod(Y) / n
  
  gamma <- if (is.null(gamma_init)) matrix(0, p, q) else gamma_init
  Omega <- if (is.null(Omega_init)) diag(q)          else Omega_init
  
  if (track_loss) {
    loss_hist <- numeric(max_iter + 1)
    time_hist <- numeric(max_iter + 1)
    t0        <- Sys.time()
    loss_hist[1] <- loss_fn(gamma, Omega, Sxx, Sxy, Syy, n, lambda_gamma, lambda_Omega)
    time_hist[1] <- 0
  }
  
  for (m in 1:max_iter) {
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > time_budget) {
      m <- m - 1; break
    }
    
    # Gradients at current iterate
    Gg <- grad_gamma(gamma, Omega, Sxx, Sxy)
    Go <- grad_Omega(gamma, Omega, Sxx, Sxy, Syy)
    Omega_inv <- safe_inv_psd(Omega)
    
    # === Solve gamma-subproblem (quadratic model + L1) ===
    gamma_star <- solve_gamma_subproblem(gamma, Omega_inv, Sxx, Gg, lambda_gamma,
                                         inner_max = 10, inner_tol = 1e-6)
    d_gamma <- gamma_star - gamma
    
    # === Solve Omega-subproblem (around Newton center) ===
    # Newton (unpenalised) direction: d_Omega_newton = - H_Omega^{-1}(Go) = -2 * Omega * Go * Omega
    d_Omega_newton <- -2 * Omega %*% Go %*% Omega
    Omega_center   <- symmetrize(Omega + d_Omega_newton)     # Newton center
    
    # Inexact prox around the center using identity-metric ADMM
    Omega_star <- admm_psd(Omega_center, lambda_Omega)
    d_Omega    <- symmetrize(Omega_star - Omega)
    
    # === Self-concordant step length alpha = 1/(1 + lambda_local) ===
    # lambda_local^2 = <d_gamma, H_gamma d_gamma> + <d_Omega, H_Omega d_Omega>
    # where <d_gamma, H_gamma d_gamma> = tr( Omega^{-1} d_gamma^T Sxx d_gamma )
    # and   <d_Omega, H_Omega d_Omega> = 0.5 * || Omega^{-1/2} d_Omega Omega^{-1/2} ||_F^2
    # Compute via Cholesky factors for stability
    # Ensure Sxx and Omega_inv are PD for chol; add tiny jitter if needed
    jitter <- 1e-10
    Sxx_j   <- Sxx + diag(jitter, nrow(Sxx))
    Oinv_j  <- Omega_inv + diag(jitter, nrow(Omega_inv))
    Ax <- tryCatch(chol(Sxx_j), error = function(e) NULL)
    Ao <- tryCatch(chol(Oinv_j), error = function(e) NULL)
    if (is.null(Ax) || is.null(Ao)) { Ax <- diag(nrow(Sxx)); Ao <- diag(nrow(Omega_inv)) }
    
    W1 <- Ax %*% d_gamma %*% t(Ao)
    term1 <- sum(W1 * W1)
    M  <- Ao %*% d_Omega %*% t(Ao)
    term2 <- 0.5 * sum(M * M)
    lambda_local <- sqrt(max(0, term1 + term2))
    alpha <- 1 / (1 + lambda_local)
    
    # Candidate update with possible backtracking safeguard
    gamma_cand <- gamma + alpha * d_gamma
    Omega_cand <- symmetrize(Omega + alpha * d_Omega)
    if (!is_psd_chol(Omega_cand)) Omega_cand <- proj_psd(Omega_cand)
    
    L_old <- loss_fn(gamma, Omega, Sxx, Sxy, Syy, n, lambda_gamma, lambda_Omega)
    L_new <- loss_fn(gamma_cand, Omega_cand, Sxx, Sxy, Syy, n, lambda_gamma, lambda_Omega)
    
    if (backtrack && (is.na(L_new) || L_new > L_old)) {
      bt <- 0
      while (bt < 10 && (is.na(L_new) || L_new > L_old)) {
        alpha <- 0.5 * alpha
        gamma_cand <- gamma + alpha * d_gamma
        Omega_cand <- symmetrize(Omega + alpha * d_Omega)
        if (!is_psd_chol(Omega_cand)) Omega_cand <- proj_psd(Omega_cand)
        L_new <- loss_fn(gamma_cand, Omega_cand, Sxx, Sxy, Syy, n, lambda_gamma, lambda_Omega)
        bt <- bt + 1
      }
    }
    
    if (track_loss) {
      loss_hist[m + 1] <- L_new
      time_hist[m + 1] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    
    gamma <- gamma_cand; Omega <- Omega_cand
    
    if (abs(L_old - L_new) < tol && lambda_local < sqrt(tol)) break
  }
  
  out <- list(gamma = gamma, Omega = Omega, iters = m, Sxx = Sxx, Sxy = Sxy, Syy = Syy, n = n)
  if (track_loss) { out$loss_hist <- loss_hist[1:(m + 1)]; out$time_hist <- time_hist[1:(m + 1)] }
  out
}

############################  Benchmark routine  ############################
# Compare Prox-Gradient (PG) vs Prox-Newton (PN) under the SAME time budget.

benchmark_pg_vs_pn <- function(n = 200, p = 30, q = 6, seed = 1,
                               time_budget = 2.0,   # seconds per method
                               ratio_lambda = 0.2   # choose lambdas as this fraction of max
) {
  set.seed(seed)
  # Simulate sparse Beta and SPD Omega
  Beta_true <- matrix(0, p, q); Beta_true[1:5, c(1,3,5)] <- rnorm(15)
  W <- symmetrize(matrix(runif(q^2, -0.2, 0.2), q, q)); diag(W) <- runif(q, 1.0, 2.0)
  Omega_true <- proj_psd(W + q * diag(q) * 0)  # ensure PSD
  Sigma_true <- safe_inv_psd(Omega_true)
  
  X <- matrix(rnorm(n * p), n, p)
  Y <- X %*% Beta_true + MASS::mvrnorm(n, mu = rep(0, q), Sigma = Sigma_true)
  
  # Center X and add intercept for training design
  Xc <- cbind(1, scale(X, center = TRUE, scale = FALSE))
  n_ <- nrow(Xc)
  Sxy_g <- crossprod(Xc, Y) / n_
  Syy_g <- crossprod(Y) / n_
  
  # Max lambdas that zero out (except intercept) per notes: Lambda_gamma = ||Sxy||_inf (excluding intercept row);
  # Lambda_Omega = 0.5 * ||Syy - diag(diag(Syy))||_inf
  lambda_gamma_0 <- max(abs(Sxy_g[-1, , drop = FALSE]))
  lambda_Omega_0 <- 0.5 * max(abs(Syy_g - diag(diag(Syy_g))))
  lambda_gamma <- ratio_lambda * lambda_gamma_0
  lambda_Omega <- ratio_lambda * lambda_Omega_0
  
  # Common init
  p_ <- ncol(Xc); q_ <- ncol(Y)
  gamma0 <- matrix(0, p_, q_)
  Omega0 <- diag(q_)
  
  # Run Prox-Gradient with time budget
  pg_fit <- pg_sparse_mvreg(Xc, Y, lambda_gamma, lambda_Omega,
                            gamma_init = gamma0, Omega_init = Omega0,
                            max_iter = 5000, tol = 1e-6, L0 = 1,
                            track_loss = TRUE, time_budget = time_budget)
  
  # Run Prox-Newton with same time budget
  pn_fit <- pn_sparse_mvreg(Xc, Y, lambda_gamma, lambda_Omega,
                            gamma_init = gamma0, Omega_init = Omega0,
                            max_iter = 1000, tol = 1e-6, track_loss = TRUE,
                            time_budget = time_budget, backtrack = TRUE)
  
  # Prepare tidy data for plotting
  df <- rbind(
    data.frame(method = "Prox-Gradient (1st)", time = pg_fit$time_hist, NLL = pg_fit$loss_hist),
    data.frame(method = "Prox-Newton (2nd)",   time = pn_fit$time_hist, NLL = pn_fit$loss_hist)
  )
  
  list(pg = pg_fit, pn = pn_fit, df = df,
       lambdas = c(lambda_gamma = lambda_gamma, lambda_Omega = lambda_Omega))
}

########################  Diagnostic contour plotters  ######################
# Visual check near a fitted solution: contour the penalised objective along
# two random directions (gamma and Omega blocks jointly) as in the 1st-order code.

plot_objective_contour_at_fit <- function(fit, lambda_gamma, lambda_Omega,
                                          t_range = 0.3, ngrid = 61, nlevels = 15,
                                          seed = 2025, main = "Objective contours near solution") {
  set.seed(seed)
  gamma_sol <- fit$gamma; Omega_sol <- fit$Omega
  Sxx <- fit$Sxx; Sxy <- fit$Sxy; Syy <- fit$Syy; n <- fit$n
  p <- nrow(gamma_sol); q <- ncol(gamma_sol)
  
  d1_g <- matrix(rnorm(p * q), p, q)
  d2_g <- matrix(rnorm(p * q), p, q)
  d1_o <- symmetrize(matrix(rnorm(q * q), q, q))
  d2_o <- symmetrize(matrix(rnorm(q * q), q, q))
  
  norm1 <- sqrt(sum(d1_g^2) + sum(d1_o^2))
  norm2 <- sqrt(sum(d2_g^2) + sum(d2_o^2))
  d1_g <- d1_g / norm1; d1_o <- d1_o / norm1
  d2_g <- d2_g / norm2; d2_o <- d2_o / norm2
  
  tgrid <- seq(-t_range, t_range, length.out = ngrid)
  Loss_grid <- outer(tgrid, tgrid, Vectorize(function(t1, t2) {
    gamma_p <- gamma_sol + t1 * d1_g + t2 * d2_g
    Omega_p <- symmetrize(Omega_sol + t1 * d1_o + t2 * d2_o)
    Omega_p <- clip_eig(Omega_p)
    loss_fn(gamma_p, Omega_p, Sxx, Sxy, Syy, n,
            lambda_gamma = lambda_gamma, lambda_Omega = lambda_Omega)
  }))
  Loss_grid[!is.finite(Loss_grid)] <- NA
  contour(tgrid, tgrid, Loss_grid, nlevels = nlevels,
          xlab = expression(t[1]), ylab = expression(t[2]),
          main = main)
  points(0, 0, pch = 19, col = "red")
  abline(h = 0, v = 0, lty = 3, col = "grey70")
  grid()
}

###############################  Demo script  ###############################
# Example usage: run the benchmark, draw NLL-vs-time, and make contour plots
# at the fitted PG and PN solutions (to visually assess local convex bowl and
# convergence vicinity).
#
library(MASS); library(ggplot2)
res <- benchmark_pg_vs_pn(n = 200, p = 30, q = 6, time_budget = 2.0, ratio_lambda = 0.2)
#
# ## 1) NLL vs time under SAME budget
ggplot(res$df, aes(time, NLL, colour = method)) +
  geom_line(size = 1.2) + theme_minimal(base_size = 14) +
  labs(title = "Train objective vs time (same budget)",
       x = "Time (sec)", y = "Penalized NLL", colour = NULL)

# ## 2) Contours near PN and PG solutions (visual convergence check)
par(mfrow = c(1, 2))
plot_objective_contour_at_fit(res$pn, res$lambdas["lambda_gamma"], res$lambdas["lambda_Omega"],
                               main = "PN: objective contours near solution")
plot_objective_contour_at_fit(res$pg, res$lambdas["lambda_gamma"], res$lambdas["lambda_Omega"],
                               main = "PG: objective contours near solution")
par(mfrow = c(1, 1))

# Example: run the benchmark and draw the NLL-vs-time comparison.
# You can source() this file and then run:
res <- benchmark_pg_vs_pn(n = 200, p = 30, q = 6, time_budget = 2.0, ratio_lambda = 0.2)
library(ggplot2)
ggplot(res$df, aes(time, NLL, colour = method)) + geom_line(size = 1.2) +
  theme_minimal(base_size = 14) + 
  labs(title = "Train objective vs time (same budget)",
                                       x = "Time (sec)", y = "Penalized NLL", colour = NULL)

# Notes:
# - The PN implementation uses an inexact Omega-prox (identity-metric ADMM) around the Newton center,
#   but with self-concordant damping it typically converges much faster than PG in practice.
# - If you want stricter inner solves, increase inner_max in solve_gamma_subproblem() and max_admm in admm_psd().
# - For very large q, consider using faster eigen/Cholesky or low-rank approximations.
