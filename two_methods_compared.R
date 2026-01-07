setwd("/Users/nominora/Desktop/MRCE")
pkgload::load_all(".")

## MRCE vs PN comparison with per-mouse 10/5 splits, 5-fold CV (8/2 within the 10), 100 reps
suppressPackageStartupMessages({
  library(readxl)
  library(MRCE)
  library(ggplot2)
})


t_start <- Sys.time()

message("Reading data...")
df <- read_excel("Data_Cortex_Nuclear.xls")
df <- as.data.frame(df)

# Factorize
for (nm in c("Genotype", "Treatment", "Behavior")) df[[nm]] <- factor(df[[nm]])

# Protein cols 2:78
stopifnot(ncol(df) >= 78)
y_cols <- 2:78

# Median impute Y
for (j in y_cols) {
  v <- as.numeric(df[[j]])
  v[!is.finite(v)] <- NA
  med <- median(v, na.rm = TRUE)
  if (!is.finite(med)) med <- 0
  v[is.na(v)] <- med
  df[[j]] <- v
}

Y_full <- as.matrix(df[, y_cols, drop = FALSE])
X_full <- model.matrix(~ Genotype + Treatment + Behavior - 1, data = df)

stopifnot(all(is.finite(X_full)), all(is.finite(Y_full)))

# Mouse grouping: strip suffix after underscore
df$MouseID_base <- sub("_.*$", "", as.character(df$MouseID))

# Train/test split helper: for every mouse with >=15 rows, sample 10 for train (for CV) and 5 for test
draw_train_test <- function() {
  mouse_tab <- table(df$MouseID_base)
  mouse_ids <- names(mouse_tab)[mouse_tab >= 15]
  if (!length(mouse_ids)) stop("No MouseID_base with >=15 measurements; check data.")
  train_idx <- integer(0)
  test_idx <- integer(0)
  train_mouse <- character(0)
  for (id in mouse_ids) {
    rows <- which(df$MouseID_base == id)
    rows <- rows[1:15]
    tr <- sample(rows, 10, replace = FALSE) # 10 used for CV training
    te <- setdiff(rows, tr) # remaining 5 for held-out test
    train_idx <- c(train_idx, tr)
    test_idx <- c(test_idx, te)
    train_mouse <- c(train_mouse, rep(id, length(tr)))
  }
  list(train_idx = train_idx, test_idx = test_idx, train_mouse = train_mouse, mouse_ids = mouse_ids)
}

# Helper: scale using training stats
scale_with <- function(X, Y, rows) {
  centerY <- colMeans(Y[rows, , drop = FALSE])
  scaleY <- apply(Y[rows, , drop = FALSE], 2, sd)
  centerX <- colMeans(X[rows, , drop = FALSE])
  scaleX <- apply(X[rows, , drop = FALSE], 2, sd)
  scaleY[scaleY == 0] <- 1
  scaleX[scaleX == 0] <- 1
  list(
    X_tr = scale(X[rows, , drop = FALSE], center = centerX, scale = scaleX),
    Y_tr = scale(Y[rows, , drop = FALSE], center = centerY, scale = scaleY),
    centerX = centerX, scaleX = scaleX,
    centerY = centerY, scaleY = scaleY
  )
}

# Fold assignment on training mice: each mouse’s 10 rows -> 5 folds (2 rows per fold)
make_folds <- function(mouse_ids, k = 5) {
  folds <- vector("list", k)
  uniq_ids <- unique(mouse_ids)
  for (id in uniq_ids) {
    rows <- which(mouse_ids == id)
    if (length(rows) %% k != 0) stop("Cannot split mouse ", id, " rows evenly into ", k, " folds.")
    per_fold <- length(rows) / k # here 10 / 5 = 2
    assign_seq <- sample(rep(seq_len(k), each = per_fold))
    for (f in seq_len(k)) folds[[f]] <- c(folds[[f]], rows[assign_seq == f])
  }
  folds
}

# Lambda grids (MRCE vs PN)
# MRCE: lam1.vec, lam2.vec (decreasing, independent)
lam_mrce_vec <- rev(exp(seq(log(1e-2), log(1), length.out = 8)))  # decreasing as required by mrce(cv)
# PN: data-driven path per split (computed inside loop) with length 60, ratio 1e-3
pn_path_len <- 60
pn_ratio <- 1e-3

# Load PN functions only
load_pn_funs <- function(path) {
  if (!file.exists(path)) stop("missing ", path)
  env <- new.env(parent = baseenv())
  exprs <- parse(path)
  for (e in exprs) {
    if (is.call(e) && length(e) >= 3 && as.character(e[[1]]) %in% c("<-", "=")) {
      rhs <- e[[3]]
      if (is.call(rhs) && identical(rhs[[1]], as.symbol("function"))) {
        eval(e, env)
      }
    }
  }
  env
}
pn_env <- load_pn_funs("11.7_need_to_change.R")
assign("pn_prox_ctrl", list(max_admm = 1000, tol = 1e-5, eta = 1e-8), envir = pn_env)
assign("armijo_c", 0.1, envir = pn_env)
pn_sparse_mvreg <- get("pn_sparse_mvreg", envir = pn_env)

# Configuration
target_reps <- 100L
results_df <- data.frame(
  rep = integer(target_reps),
  mse_mrce = numeric(target_reps),
  mse_pn_mse = numeric(target_reps),
  time_mrce = numeric(target_reps),
  time_pn_mse = numeric(target_reps)
)
better_mse_count <- 0L

message(sprintf("Starting %d repetitions...", target_reps))

for (rep_i in 1:target_reps) {
  set.seed(12345 + rep_i)
  
  split <- draw_train_test()
  train_idx <- split$train_idx
  test_idx <- split$test_idx
  train_mouse <- split$train_mouse
  
  X_tr_raw <- X_full[train_idx, , drop = FALSE]
  Y_tr_raw <- Y_full[train_idx, , drop = FALSE]
  X_te_raw <- X_full[test_idx, , drop = FALSE]
  Y_te_raw <- Y_full[test_idx, , drop = FALSE]
  mouse_tr <- train_mouse
  
  sdY <- apply(Y_tr_raw, 2, sd)
  keep_y <- which(is.finite(sdY) & sdY > 0)
  Y_tr_raw <- Y_tr_raw[, keep_y, drop = FALSE]
  Y_te_raw <- Y_te_raw[, keep_y, drop = FALSE]
  
  sdX <- apply(X_tr_raw, 2, sd)
  keep_x <- which(is.finite(sdX) & sdX > 0)
  X_tr_raw <- X_tr_raw[, keep_x, drop = FALSE]
  X_te_raw <- X_te_raw[, keep_x, drop = FALSE]
  
  folds <- make_folds(mouse_tr, k = 5)
  
  # MRCE CV
  mrce_cv_fit <- mrce(
    X = X_tr_raw, Y = Y_tr_raw,
    lam1.vec = lam_mrce_vec, lam2.vec = lam_mrce_vec,
    method = "cv", kfold = 5,folds = folds,
    standardize = FALSE, silent = TRUE
  )
  mrce_cv_err <- mrce_cv_fit$cv.err
  best_idx <- which(mrce_cv_err == min(mrce_cv_err), arr.ind = TRUE)[1, ]
  lam_best_mrce1 <- lam_mrce_vec[best_idx[1]]
  lam_best_mrce2 <- lam_mrce_vec[best_idx[2]]
  
  # PN CV
  # PN lambda sequences derived from training stats (following 11.7_need_to_change.R style)
  n_tr <- nrow(X_tr_raw)
  X_tr_pn_center <- scale(X_tr_raw, center = TRUE, scale = FALSE)
  X_tr_pn_full <- cbind(Intercept = 1, X_tr_pn_center)
  Sxy_tr <- crossprod(X_tr_pn_full, Y_tr_raw) / n_tr
  Syy_tr <- crossprod(Y_tr_raw) / n_tr
  lambda_gamma_0 <- max(abs(Sxy_tr[-1, , drop = FALSE]))
  tmp_off <- Syy_tr - diag(diag(Syy_tr))
  lambda_Omega_0 <- max(abs(tmp_off))
  lam_gamma_vec <- exp(seq(log(lambda_gamma_0), log(lambda_gamma_0 * pn_ratio), length.out = pn_path_len))
  lam_omega_vec <- exp(seq(log(lambda_Omega_0), log(lambda_Omega_0 * pn_ratio), length.out = pn_path_len))

  pn_cv_mse <- numeric(pn_path_len)
  pn_warm_gamma <- vector("list", length(folds))
  pn_warm_Omega <- vector("list", length(folds))
  for (li in seq_len(pn_path_len)) {
    lam_g <- lam_gamma_vec[li]
    lam_o <- lam_omega_vec[li]
    fold_mse <- numeric(length(folds))
    for (i in seq_along(folds)) {
      val_rows <- folds[[i]]
      tr_rows <- setdiff(seq_len(nrow(X_tr_raw)), val_rows)
      sc <- scale_with(X_tr_raw, Y_tr_raw, tr_rows)
      X_tr_f <- sc$X_tr
      Y_tr_f <- sc$Y_tr
      X_val_f <- scale(X_tr_raw[val_rows, , drop = FALSE], center = sc$centerX, scale = sc$scaleX)
      Y_val_f <- scale(Y_tr_raw[val_rows, , drop = FALSE], center = sc$centerY, scale = sc$scaleY)
      X_tr_pn <- cbind(Intercept = rep(1, nrow(X_tr_f)), X_tr_f)
      X_val_pn <- cbind(Intercept = rep(1, nrow(X_val_f)), X_val_f)
      fit <- pn_sparse_mvreg(
        Xc = X_tr_pn,
        Y = Y_tr_f,
        lambda_gamma = lam_g,
        lambda_Omega = lam_o,
        gamma_init = if (li > 1) pn_warm_gamma[[i]] else NULL,
        Omega_init = if (li > 1) pn_warm_Omega[[i]] else diag(ncol(Y_tr_f)),
        max_iter = 30,
        tol = 1e-4
      )
      pn_warm_gamma[[i]] <- fit$gamma
      pn_warm_Omega[[i]] <- fit$Omega
      beta <- fit$gamma %*% solve(fit$Omega)
      mu <- as.numeric(colMeans(Y_tr_f) - t(beta) %*% colMeans(X_tr_pn))
      Yhat <- tcrossprod(matrix(1, nrow = nrow(X_val_pn), ncol = 1), mu) + X_val_pn %*% beta
      fold_mse[i] <- mean((Y_val_f - Yhat)^2)
    }
    pn_cv_mse[li] <- mean(fold_mse)
  }
  best_idx_pn <- which.min(pn_cv_mse)
  lam_best_pn_mse <- lam_gamma_vec[best_idx_pn]
  lam_best_pn_mse_Omega <- lam_omega_vec[best_idx_pn]
  
  # Refit & Test
  sc_full <- scale_with(X_tr_raw, Y_tr_raw, rows = seq_len(nrow(X_tr_raw)))
  X_tr <- sc_full$X_tr
  Y_tr <- sc_full$Y_tr
  X_te <- scale(X_te_raw, center = sc_full$centerX, scale = sc_full$scaleX)
  Y_te <- scale(Y_te_raw, center = sc_full$centerY, scale = sc_full$scaleY)
  
  # MRCE Refit
  t0_mrce <- Sys.time()
  mrce_fit <- mrce(
    X = X_tr, Y = Y_tr, lam1 = lam_best_mrce1, lam2 = lam_best_mrce2,
    method = "single", standardize = FALSE, silent = TRUE
  )
  time_mrce <- as.numeric(difftime(Sys.time(), t0_mrce, units = "secs"))
  
  beta_mrce <- mrce_fit$Bhat
  mu_mrce <- as.numeric(mrce_fit$muhat)
  Yhat_mrce <- tcrossprod(matrix(1, nrow = nrow(X_te), ncol = 1), mu_mrce) + X_te %*% beta_mrce
  mse_mrce_test <- mean((Y_te - Yhat_mrce)^2)
  
  # PN Refit
  X_tr_pn <- cbind(Intercept = rep(1, nrow(X_tr)), X_tr)
  X_te_pn <- cbind(Intercept = rep(1, nrow(X_te)), X_te)
  
  t0_pn_mse <- Sys.time()
  pn_fit_mse <- pn_sparse_mvreg(
    Xc = X_tr_pn,
    Y = Y_tr,
    lambda_gamma = lam_best_pn_mse,
    lambda_Omega = lam_best_pn_mse_Omega,
    gamma_init = NULL,
    Omega_init = diag(ncol(Y_tr)),
    max_iter = 30,
    tol = 1e-4
  )
  time_pn_mse <- as.numeric(difftime(Sys.time(), t0_pn_mse, units = "secs"))
  
  beta_pn_mse <- pn_fit_mse$gamma %*% solve(pn_fit_mse$Omega)
  mu_pn_mse <- as.numeric(colMeans(Y_tr) - t(beta_pn_mse) %*% colMeans(X_tr_pn))
  Yhat_pn_mse <- tcrossprod(matrix(1, nrow = nrow(X_te_pn), ncol = 1), mu_pn_mse) + X_te_pn %*% beta_pn_mse
  mse_pn_mse_test <- mean((Y_te - Yhat_pn_mse)^2)
  
  if (mse_pn_mse_test < mse_mrce_test) better_mse_count <- better_mse_count + 1L
  
  results_df[rep_i, ] <- c(rep_i, mse_mrce_test, mse_pn_mse_test,
                           time_mrce, time_pn_mse)
  
  if (rep_i %% 10 == 0) {
    curr_mean_mrce <- mean(results_df$mse_mrce[1:rep_i])
    curr_mean_pn <- mean(results_df$mse_pn_mse[1:rep_i])
    message(sprintf(
      "Rep %d/%d | Mean MSE: MRCE=%.4f, PN=%.4f | PN Better: %d/%d",
      rep_i, target_reps, curr_mean_mrce, curr_mean_pn,
      better_mse_count, rep_i
    ))
  }
}

t_end <- Sys.time()
total_time <- as.numeric(difftime(t_end, t_start, units = "secs"))

# Final Statistics
calc_stats <- function(x) {
  m <- mean(x)
  s <- sd(x)
  se <- s / sqrt(length(x))
  c(mean = m, sd = s, se = se)
}

stats_mrce_mse <- calc_stats(results_df$mse_mrce)
stats_pn_mse <- calc_stats(results_df$mse_pn_mse)

cat("\n=== Final Results (100 Repetitions) ===\n")
cat(sprintf("Total Runtime: %.1f seconds\n\n", total_time))

cat("--- Test MSE ---\n")
cat(sprintf("MRCE: Mean=%.4f, SD=%.4f, SE=%.4f\n", stats_mrce_mse["mean"], stats_mrce_mse["sd"], stats_mrce_mse["se"]))
cat(sprintf("PN  : Mean=%.4f, SD=%.4f, SE=%.4f\n", stats_pn_mse["mean"], stats_pn_mse["sd"], stats_pn_mse["se"]))
cat(sprintf(
  "PN better than MRCE in %d / %d runs (%.1f%%)\n",
  better_mse_count, target_reps, 100 * better_mse_count / target_reps
))

saveRDS(results_df, "compare_results_modified.rds")
