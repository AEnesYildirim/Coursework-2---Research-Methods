# =========================================================
# Coursework 2: submission-ready main estimation script
# =========================================================
# Model:
#   - one-dimensional latent bargaining factor mu_i
#   - anchor on observed consumption-share proxy s_m
#   - logit-transform s_m
#   - joint MSL for measurement block, bargaining equation,
#     female participation, and female wage/income with selection
#
# Estimation design:
#   - Halton draws only
#   - one estimation path only: Nelder-Mead then BFGS
#   - optional warm start from a previously saved theta vector
#   - repaired positive-definite Hessian-based VCOV
#   - repaired robust sandwich VCOV for main reported SEs
#
# Notes:
#   - This is a standalone script intended for the main submission run.
#   - Robustness checks should be run separately later.
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

# -----------------------------
# 0. Utility helpers
# -----------------------------

clip01 <- function(x, eps = 1e-6) {
  pmin(pmax(x, eps), 1 - eps)
}

row_log_mean_exp <- function(log_mat) {
  m <- apply(log_mat, 1L, max)
  m + log(rowMeans(exp(log_mat - m)))
}

halton_1d <- function(n, base = 2L, start_index = 11L) {
  radical_inverse <- function(index, base) {
    f <- 1 / base
    r <- 0
    i <- index
    while (i > 0) {
      r <- r + f * (i %% base)
      i <- floor(i / base)
      f <- f / base
    }
    r
  }

  vapply(start_index + seq_len(n) - 1L, radical_inverse, numeric(1), base = base)
}

make_halton_draws <- function(R, seed = 125, eps = 1e-10) {
  if (R < 2L) stop("R must be at least 2.")
  set.seed(seed)
  shift <- runif(1)
  u <- (halton_1d(R, base = 2L, start_index = 11L) + shift) %% 1
  qnorm(clip01(u, eps = eps))
}

repeat_draws_by_row <- function(draws, N) {
  matrix(draws, nrow = N, ncol = length(draws), byrow = TRUE)
}

is_binary_like <- function(x) {
  x_no_na <- unique(stats::na.omit(x))
  length(x_no_na) <= 2L && all(x_no_na %in% c(0, 1))
}

standardise_measurement_block <- function(Y_raw, anchor_name) {
  Y <- Y_raw
  centers <- rep(0, ncol(Y_raw))
  scales <- rep(1, ncol(Y_raw))
  names(centers) <- colnames(Y_raw)
  names(scales) <- colnames(Y_raw)

  for (j in seq_len(ncol(Y_raw))) {
    nm <- colnames(Y_raw)[j]
    if (nm == anchor_name) next

    mu <- mean(Y_raw[, j], na.rm = TRUE)
    sdj <- sd(Y_raw[, j], na.rm = TRUE)
    if (!is.finite(sdj) || sdj <= 0) sdj <- 1

    Y[, j] <- (Y_raw[, j] - mu) / sdj
    centers[nm] <- mu
    scales[nm] <- sdj
  }

  list(Y = Y, centers = centers, scales = scales)
}

standardise_continuous_df <- function(df, binary_names = character()) {
  out <- as.data.frame(df)
  centers <- rep(0, ncol(out))
  scales <- rep(1, ncol(out))
  standardized <- rep(FALSE, ncol(out))
  names(centers) <- colnames(out)
  names(scales) <- colnames(out)
  names(standardized) <- colnames(out)

  for (nm in colnames(out)) {
    x <- out[[nm]]
    if (!is.numeric(x)) next
    if (nm %in% binary_names) next
    if (is_binary_like(x)) next

    mu <- mean(x, na.rm = TRUE)
    sdj <- sd(x, na.rm = TRUE)
    if (!is.finite(sdj) || sdj <= 0) sdj <- 1

    out[[nm]] <- (x - mu) / sdj
    centers[nm] <- mu
    scales[nm] <- sdj
    standardized[nm] <- TRUE
  }

  list(data = out, centers = centers, scales = scales, standardized = standardized)
}

make_pd_inverse <- function(H,
                            eig_floor_abs = 1e-6,
                            eig_floor_rel = 1e-8) {
  Hs <- 0.5 * (H + t(H))
  ee <- eigen(Hs, symmetric = TRUE)
  vals <- ee$values
  vecs <- ee$vectors

  floor_val <- max(eig_floor_abs, eig_floor_rel * max(abs(vals), na.rm = TRUE))
  vals_fixed <- pmax(vals, floor_val)

  H_pd <- vecs %*% diag(vals_fixed, nrow = length(vals_fixed)) %*% t(vecs)
  V_pd <- vecs %*% diag(1 / vals_fixed, nrow = length(vals_fixed)) %*% t(vecs)

  list(
    hessian_pd = H_pd,
    vcov_pd = V_pd,
    eigen_original = vals,
    eigen_fixed = vals_fixed,
    floor_val = floor_val,
    n_clipped = sum(vals < floor_val)
  )
}

make_psd_matrix <- function(M,
                            eig_floor_abs = 1e-10,
                            eig_floor_rel = 1e-10) {
  Ms <- 0.5 * (M + t(M))
  ee <- eigen(Ms, symmetric = TRUE)
  vals <- ee$values
  vecs <- ee$vectors

  floor_val <- max(eig_floor_abs, eig_floor_rel * max(abs(vals), na.rm = TRUE))
  vals_fixed <- pmax(vals, floor_val)

  list(
    matrix_psd = vecs %*% diag(vals_fixed, nrow = length(vals_fixed)) %*% t(vecs),
    eigen_original = vals,
    eigen_fixed = vals_fixed,
    floor_val = floor_val,
    n_clipped = sum(vals < floor_val)
  )
}

# -----------------------------
# 1. Variable lists and sample construction
# -----------------------------

get_variable_lists <- function(include_lr_inc = FALSE) {
  measurement_vars <- c(
    "tasksunfair.proj.m", "tasksunfair.proj.f",
    "subjective.rho.proj.m", "subjective.rho.proj.f",
    "yieldfor.proj.m", "yieldfor.proj.f",
    "helps.proj.m", "helps.proj.f",
    "harderfor.proj.m", "harderfor.proj.f",
    "soledecisions.proj.m", "soledecisions.proj.f"
  )

  trait_vars <- c(
    "friendly.proj.m", "friendly.proj.f",
    "open.proj.m", "open.proj.f",
    "patient.proj.m", "patient.proj.f",
    "understanding.proj.m", "understanding.proj.f",
    "empathy.proj.m", "empathy.proj.f",
    "tolerant.proj.m", "tolerant.proj.f",
    "critical.proj.m", "critical.proj.f",
    "lazy.proj.m", "lazy.proj.f",
    "dominant.proj.m", "dominant.proj.f",
    "emotional.proj.m", "emotional.proj.f",
    "moody.proj.m", "moody.proj.f",
    "thoughtless.proj.m", "thoughtless.proj.f",
    "unreasonable.proj.m", "unreasonable.proj.f",
    "distant.proj.m", "distant.proj.f",
    "complaining.proj.m", "complaining.proj.f"
  )

  z_vars_all <- c("duration", "abouttobeparents", "lr_inc", "d_educ", "d_age", "both_m")
  z_vars_mu <- if (include_lr_inc) z_vars_all else setdiff(z_vars_all, "lr_inc")

  list(
    measurement_vars = measurement_vars,
    trait_vars = trait_vars,
    z_vars_mu = z_vars_mu,
    z_vars_all = z_vars_all
  )
}

prepare_base_df <- function(df, vars, use_logit_sm = TRUE) {
  df %>%
    mutate(
      both_m = as.integer(both_m),
      D_f = case_when(
        is.na(income.f) ~ NA_integer_,
        income.f > 0    ~ 1L,
        TRUE            ~ 0L
      ),
      W_f = case_when(
        !is.na(income.f) & income.f > 0 ~ log(pmax(income.f, 1e-8)),
        TRUE                            ~ NA_real_
      ),
      s_m_input = if (use_logit_sm) qlogis(clip01(s_m)) else s_m
    )
}

make_sample_flags <- function(df, vars) {
  base_y_names <- c("s_m_input", vars$measurement_vars)

  tibble(
    keep_measurement = complete.cases(df[, base_y_names, drop = FALSE]),
    keep_bargaining = complete.cases(df[, c(vars$trait_vars, vars$z_vars_mu), drop = FALSE]),
    keep_labour = complete.cases(df[, c("age.f", "educ.f", "D_f"), drop = FALSE]),
    keep_wage = df$D_f == 0L | (df$D_f == 1L & !is.na(df$W_f))
  )
}

make_accounting_table <- function(flags) {
  tibble(
    stage = c(
      "Raw rows",
      "After measurement block completeness",
      "After bargaining-covariate completeness",
      "After labour-covariate completeness",
      "Final joint-MSL sample"
    ),
    n = c(
      nrow(flags),
      sum(flags$keep_measurement),
      sum(flags$keep_measurement & flags$keep_bargaining),
      sum(flags$keep_measurement & flags$keep_bargaining & flags$keep_labour),
      sum(flags$keep_measurement & flags$keep_bargaining & flags$keep_labour & flags$keep_wage)
    )
  )
}

build_joint_sample <- function(data_path = "prepdata.csv",
                               include_lr_inc = FALSE,
                               use_logit_sm = TRUE,
                               anchor_indicator = "s_m_input") {
  vars <- get_variable_lists(include_lr_inc = include_lr_inc)
  df <- readr::read_csv(data_path, show_col_types = FALSE)

  if (!("hid" %in% names(df))) {
    df <- df %>% mutate(hid = row_number())
  }

  required_vars <- c(
    "hid", "s_m", "age.f", "educ.f", "income.f",
    vars$measurement_vars, vars$trait_vars, vars$z_vars_all
  )
  missing_cols <- setdiff(required_vars, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- prepare_base_df(df, vars, use_logit_sm = use_logit_sm)
  flags <- make_sample_flags(df, vars)
  accounting_tbl <- make_accounting_table(flags)

  keep_final <- with(flags, keep_measurement & keep_bargaining & keep_labour & keep_wage)
  base_y_names <- c("s_m_input", vars$measurement_vars)
  if (!(anchor_indicator %in% base_y_names)) {
    stop("anchor_indicator must be one of: ", paste(base_y_names, collapse = ", "))
  }
  y_names <- c(anchor_indicator, setdiff(base_y_names, anchor_indicator))

  joint_df <- df[keep_final, c(
    "hid", "s_m", "s_m_input",
    vars$measurement_vars,
    vars$trait_vars,
    vars$z_vars_all,
    "age.f", "educ.f", "D_f", "W_f"
  )]

  Y_raw <- as.matrix(joint_df[, y_names, drop = FALSE])
  std_obj <- standardise_measurement_block(Y_raw, anchor_name = anchor_indicator)
  Y <- std_obj$Y
  colnames(Y) <- y_names

  X_mu_df_raw <- joint_df %>% select(all_of(c(vars$trait_vars, vars$z_vars_mu)))
  X_mu_std <- standardise_continuous_df(
    X_mu_df_raw,
    binary_names = c("abouttobeparents", "both_m")
  )
  X_mu <- model.matrix(~ ., data = X_mu_std$data)

  Z_sel_df_raw <- joint_df %>% select(all_of(c("age.f", "educ.f", "duration", "abouttobeparents", "d_educ", "d_age", "both_m")))
  Z_sel_std <- standardise_continuous_df(
    Z_sel_df_raw,
    binary_names = c("abouttobeparents", "both_m")
  )
  Z_sel <- model.matrix(~ ., data = Z_sel_std$data)

  X_w_df_raw <- joint_df %>% select(all_of(c("age.f", "educ.f")))
  X_w_std <- standardise_continuous_df(X_w_df_raw)
  X_w <- model.matrix(~ ., data = X_w_std$data)

  list(
    raw_df = df,
    accounting = accounting_tbl,
    joint_df = joint_df,
    Y = Y,
    Y_raw = Y_raw,
    y_names = y_names,
    y_centers = std_obj$centers,
    y_scales = std_obj$scales,
    X_mu = X_mu,
    X_mu_scaling = list(
      centers = X_mu_std$centers,
      scales = X_mu_std$scales,
      standardized = X_mu_std$standardized
    ),
    Z_sel = Z_sel,
    Z_sel_scaling = list(
      centers = Z_sel_std$centers,
      scales = Z_sel_std$scales,
      standardized = Z_sel_std$standardized
    ),
    X_w = X_w,
    X_w_scaling = list(
      centers = X_w_std$centers,
      scales = X_w_std$scales,
      standardized = X_w_std$standardized
    ),
    D = as.integer(joint_df$D_f),
    W = joint_df$W_f,
    workers = which(joint_df$D_f == 1L),
    vars = vars,
    settings = list(
      include_lr_inc = include_lr_inc,
      use_logit_sm = use_logit_sm,
      anchor_indicator = anchor_indicator
    )
  )
}

# -----------------------------
# 2. Parameter mapping
# -----------------------------

make_dim_list <- function(prepped) {
  list(
    J = ncol(prepped$Y),
    K_mu = ncol(prepped$X_mu),
    K_sel = ncol(prepped$Z_sel),
    K_w = ncol(prepped$X_w)
  )
}

pack_theta_joint <- function(beta_free,
                             lambda_free,
                             log_sigma_y,
                             gamma,
                             log_sigma_u,
                             alpha_mu,
                             alpha_sel,
                             delta,
                             log_sigma_eps,
                             atanh_rho) {
  c(
    beta_free,
    lambda_free,
    log_sigma_y,
    gamma,
    log_sigma_u,
    alpha_mu,
    alpha_sel,
    delta,
    log_sigma_eps,
    atanh_rho
  )
}

unpack_theta_joint <- function(theta, dims) {
  J <- dims$J
  K_mu <- dims$K_mu
  K_sel <- dims$K_sel
  K_w <- dims$K_w

  idx <- 1L
  beta_free <- theta[idx:(idx + J - 2L)]; idx <- idx + J - 1L
  lambda_free <- theta[idx:(idx + J - 2L)]; idx <- idx + J - 1L
  log_sigma_y <- theta[idx:(idx + J - 1L)]; idx <- idx + J
  gamma <- theta[idx:(idx + K_mu - 1L)]; idx <- idx + K_mu
  log_sigma_u <- theta[idx]; idx <- idx + 1L
  alpha_mu <- theta[idx]; idx <- idx + 1L
  alpha_sel <- theta[idx:(idx + K_sel - 1L)]; idx <- idx + K_sel
  delta <- theta[idx:(idx + K_w - 1L)]; idx <- idx + K_w
  log_sigma_eps <- theta[idx]; idx <- idx + 1L
  atanh_rho <- theta[idx]

  list(
    beta = c(0, beta_free),
    lambda = c(1, lambda_free),
    sigma_y = pmax(exp(log_sigma_y), 1e-6),
    gamma = gamma,
    sigma_u = pmax(exp(log_sigma_u), 1e-6),
    alpha_mu = alpha_mu,
    alpha_sel = alpha_sel,
    delta = delta,
    sigma_eps = pmax(exp(log_sigma_eps), 1e-6),
    rho = 0.999 * tanh(atanh_rho)
  )
}

expected_theta_length <- function(prepped) {
  dims <- make_dim_list(prepped)
  (dims$J - 1L) + (dims$J - 1L) + dims$J + dims$K_mu + 1L + 1L + dims$K_sel + dims$K_w + 1L + 1L
}

# -----------------------------
# 3. Initial values
# -----------------------------

get_initial_values_joint <- function(prepped) {
  fa_start <- factanal(x = prepped$Y, factors = 1, scores = "regression")

  load0 <- as.numeric(fa_start$loadings[, 1])
  uniq0 <- pmax(as.numeric(fa_start$uniquenesses), 0.05)
  beta0 <- colMeans(prepped$Y)
  mu0 <- as.numeric(fa_start$scores[, 1])

  if (load0[1] < 0) {
    load0 <- -load0
    mu0 <- -mu0
  }
  if (abs(load0[1]) < 1e-4) load0[1] <- 1

  mu0 <- mu0 * load0[1]
  lambda0 <- load0 / load0[1]
  lambda0[1] <- 1

  cor_check <- suppressWarnings(cor(mu0, prepped$joint_df$s_m, use = "complete.obs"))
  if (is.finite(cor_check) && cor_check < 0) {
    mu0 <- -mu0
    lambda0 <- -lambda0
    lambda0[1] <- 1
  }

  beta_free0 <- beta0[-1]
  lambda_free0 <- lambda0[-1]
  log_sigma_y0 <- log(sqrt(uniq0))

  gamma_fit <- lm(mu0 ~ prepped$X_mu - 1)
  gamma0 <- coef(gamma_fit)
  gamma0[is.na(gamma0)] <- 0

  mu_mean0 <- as.numeric(prepped$X_mu %*% gamma0)
  sigma_u0 <- sd(mu0 - mu_mean0)
  if (!is.finite(sigma_u0) || sigma_u0 <= 0) sigma_u0 <- 1
  log_sigma_u0 <- log(sigma_u0)

  part_df <- bind_cols(
    tibble(D_f = prepped$D, mu0 = mu0),
    as.data.frame(prepped$Z_sel[, -1, drop = FALSE])
  )
  part_fit <- glm(D_f ~ mu0 + ., family = binomial(link = "probit"), data = part_df)
  part_coef <- coef(part_fit)
  part_coef[is.na(part_coef)] <- 0

  alpha_mu0 <- unname(part_coef["mu0"])
  alpha_sel0 <- rep(0, ncol(prepped$Z_sel))
  names(alpha_sel0) <- colnames(prepped$Z_sel)
  overlap_sel <- intersect(names(part_coef), names(alpha_sel0))
  alpha_sel0[overlap_sel] <- part_coef[overlap_sel]

  wage_df <- bind_cols(
    tibble(W_f = prepped$W, D_f = prepped$D),
    as.data.frame(prepped$X_w[, -1, drop = FALSE])
  )
  wage_fit <- lm(W_f ~ ., data = wage_df %>% filter(D_f == 1L))
  wage_coef <- coef(wage_fit)
  wage_coef[is.na(wage_coef)] <- 0

  delta0 <- rep(0, ncol(prepped$X_w))
  names(delta0) <- colnames(prepped$X_w)
  overlap_w <- intersect(names(wage_coef), names(delta0))
  delta0[overlap_w] <- wage_coef[overlap_w]

  sigma_eps0 <- summary(wage_fit)$sigma
  if (!is.finite(sigma_eps0) || sigma_eps0 <= 0) sigma_eps0 <- 1
  log_sigma_eps0 <- log(sigma_eps0)

  rho0 <- 0
  if (requireNamespace("sampleSelection", quietly = TRUE)) {
    selection_formula <- as.formula(
      paste("D_f ~ mu0 +", paste(colnames(prepped$Z_sel)[-1], collapse = " + "))
    )
    outcome_formula <- as.formula(
      paste("W_f ~", paste(colnames(prepped$X_w)[-1], collapse = " + "))
    )

    heckman_try <- try(
      sampleSelection::selection(
        selection = selection_formula,
        outcome = outcome_formula,
        data = bind_cols(
          tibble(W_f = prepped$W, D_f = prepped$D, mu0 = mu0),
          as.data.frame(prepped$Z_sel[, -1, drop = FALSE]),
          as.data.frame(prepped$X_w[, -1, drop = FALSE])
        ),
        method = "ml"
      ),
      silent = TRUE
    )

    if (!inherits(heckman_try, "try-error")) {
      heck_est <- summary(heckman_try)$estimate
      if ("rho" %in% rownames(heck_est)) {
        rho_guess <- heck_est["rho", "Estimate"]
        if (is.finite(rho_guess)) rho0 <- max(min(rho_guess, 0.95), -0.95)
      }
    }
  }

  list(
    theta = pack_theta_joint(
      beta_free = beta_free0,
      lambda_free = lambda_free0,
      log_sigma_y = log_sigma_y0,
      gamma = gamma0,
      log_sigma_u = log_sigma_u0,
      alpha_mu = alpha_mu0,
      alpha_sel = alpha_sel0,
      delta = delta0,
      log_sigma_eps = log_sigma_eps0,
      atanh_rho = atanh(rho0)
    ),
    mu0 = mu0
  )
}

# -----------------------------
# 4. Joint simulated likelihood
# -----------------------------

joint_loglik_obs <- function(theta, prepped, draws) {
  dims <- make_dim_list(prepped)
  p <- unpack_theta_joint(theta, dims)

  N <- nrow(prepped$Y)
  R <- length(draws)
  J <- ncol(prepped$Y)

  draw_mat <- repeat_draws_by_row(draws, N)
  mu_mean <- as.numeric(prepped$X_mu %*% p$gamma)
  MU <- matrix(mu_mean, nrow = N, ncol = R) + p$sigma_u * draw_mat

  log_meas <- matrix(0, nrow = N, ncol = R)
  for (j in seq_len(J)) {
    yj <- matrix(prepped$Y[, j], nrow = N, ncol = R)
    mean_j <- p$beta[j] + p$lambda[j] * MU
    log_meas <- log_meas + dnorm(yj, mean = mean_j, sd = p$sigma_y[j], log = TRUE)
  }

  sel_base <- as.numeric(prepped$Z_sel %*% p$alpha_sel)
  sel_index <- matrix(sel_base, nrow = N, ncol = R) + p$alpha_mu * MU

  log_dw <- matrix(0, nrow = N, ncol = R)
  nonworkers <- which(prepped$D == 0L)
  workers <- prepped$workers

  if (length(nonworkers) > 0) {
    log_dw[nonworkers, ] <- pnorm(
      sel_index[nonworkers, , drop = FALSE],
      lower.tail = FALSE,
      log.p = TRUE
    )
  }

  if (length(workers) > 0) {
    w_mean <- as.numeric(prepped$X_w[workers, , drop = FALSE] %*% p$delta)
    v <- (prepped$W[workers] - w_mean) / p$sigma_eps
    v_mat <- matrix(v, nrow = length(workers), ncol = R)
    denom <- sqrt(pmax(1 - p$rho^2, 1e-8))
    cond_arg <- (sel_index[workers, , drop = FALSE] + p$rho * v_mat) / denom

    log_dw[workers, ] <-
      dnorm(v_mat, mean = 0, sd = 1, log = TRUE) -
      log(p$sigma_eps) +
      pnorm(cond_arg, log.p = TRUE)
  }

  row_log_mean_exp(log_meas + log_dw)
}

objective_joint_msl <- function(theta, prepped, draws) {
  val <- -sum(joint_loglik_obs(theta, prepped = prepped, draws = draws))
  if (!is.finite(val)) return(1e12)
  val
}

# -----------------------------
# 5. One-path estimation
# -----------------------------

load_theta_start <- function(theta_start = NULL,
                             theta_start_path = NULL,
                             prepped = NULL) {
  if (!is.null(theta_start_path)) {
    saved <- readRDS(theta_start_path)
    theta_start <- if (is.list(saved) && !is.null(saved$theta)) saved$theta else saved
  }

  if (is.null(theta_start)) return(NULL)
  if (length(theta_start) != expected_theta_length(prepped)) {
    stop("theta_start length does not match the current specification.")
  }
  theta_start
}

fit_joint_msl_submission <- function(data_path = "prepdata.csv",
                                     include_lr_inc = FALSE,
                                     use_logit_sm = TRUE,
                                     anchor_indicator = "s_m_input",
                                     R = 500,
                                     seed = 125,
                                     theta_start = NULL,
                                     theta_start_path = NULL,
                                     nm_maxit = 300,
                                     bfgs_maxit = 600,
                                     trace = 1) {
  prepped <- build_joint_sample(
    data_path = data_path,
    include_lr_inc = include_lr_inc,
    use_logit_sm = use_logit_sm,
    anchor_indicator = anchor_indicator
  )

  theta_loaded <- load_theta_start(theta_start, theta_start_path, prepped)
  init <- get_initial_values_joint(prepped)
  theta0 <- if (is.null(theta_loaded)) init$theta else theta_loaded
  draws <- make_halton_draws(R = R, seed = seed)

  cat("\n============================\n")
  cat("Joint MSL submission run\n")
  cat("============================\n")
  cat("N =", nrow(prepped$joint_df), "\n")
  cat("Workers =", sum(prepped$D == 1L), "\n")
  cat("Indicators =", ncol(prepped$Y), "\n")
  cat("Anchor =", anchor_indicator, "\n")
  cat("R draws =", R, "\n")
  cat("Draw type = Halton\n")
  cat("Warm start from saved theta =", !is.null(theta_loaded), "\n")

  nm_fit <- optim(
    par = theta0,
    fn = objective_joint_msl,
    prepped = prepped,
    draws = draws,
    method = "Nelder-Mead",
    control = list(maxit = nm_maxit, trace = trace, REPORT = 25)
  )

  bfgs_fit <- optim(
    par = nm_fit$par,
    fn = objective_joint_msl,
    prepped = prepped,
    draws = draws,
    method = "BFGS",
    control = list(maxit = bfgs_maxit, trace = trace, REPORT = 25)
  )

  out <- list(
    prepped = prepped,
    init = init,
    draws = draws,
    draw_type = "halton",
    R = R,
    seed = seed,
    nm_fit = nm_fit,
    fit = bfgs_fit,
    hessian = NULL,
    hessian_pd = NULL,
    vcov = NULL,
    vcov_robust = NULL,
    hessian_diagnostics = NULL,
    robust_diagnostics = NULL
  )
  class(out) <- "cw2_joint_msl"
  out
}

# -----------------------------
# 6. Inference: repaired Hessian + robust VCOV
# -----------------------------

compute_hessian_joint_msl <- function(object,
                                      hessian_method = c("optimHess", "numDeriv")) {
  if (!inherits(object, "cw2_joint_msl")) stop("object must inherit from 'cw2_joint_msl'.")
  hessian_method <- match.arg(hessian_method)

  if (hessian_method == "numDeriv") {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("numDeriv is required for hessian_method = 'numDeriv'.")
    }
    return(numDeriv::hessian(
      func = objective_joint_msl,
      x = object$fit$par,
      prepped = object$prepped,
      draws = object$draws
    ))
  }

  optimHess(
    par = object$fit$par,
    fn = objective_joint_msl,
    prepped = object$prepped,
    draws = object$draws
  )
}

compute_inference_joint_msl <- function(object,
                                        hessian_method = c("optimHess", "numDeriv"),
                                        use_robust = TRUE,
                                        eig_floor_abs = 1e-6,
                                        eig_floor_rel = 1e-8,
                                        robust_floor_abs = 1e-10,
                                        robust_floor_rel = 1e-10) {
  if (!inherits(object, "cw2_joint_msl")) stop("object must inherit from 'cw2_joint_msl'.")
  hessian_method <- match.arg(hessian_method)

  H <- compute_hessian_joint_msl(object, hessian_method = hessian_method)
  fix_obj <- make_pd_inverse(H, eig_floor_abs = eig_floor_abs, eig_floor_rel = eig_floor_rel)

  out <- object
  out$hessian <- H
  out$hessian_pd <- fix_obj$hessian_pd
  out$vcov <- fix_obj$vcov_pd
  out$hessian_diagnostics <- tibble(
    min_eigenvalue_original = min(fix_obj$eigen_original, na.rm = TRUE),
    min_eigenvalue_fixed = min(fix_obj$eigen_fixed, na.rm = TRUE),
    max_eigenvalue_original = max(fix_obj$eigen_original, na.rm = TRUE),
    clipping_floor = fix_obj$floor_val,
    n_eigenvalues_clipped = fix_obj$n_clipped,
    has_non_finite = any(!is.finite(H))
  )

  if (use_robust) {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("numDeriv is required for robust standard errors.")
    }

    score_mat <- numDeriv::jacobian(
      func = function(th) joint_loglik_obs(th, prepped = out$prepped, draws = out$draws),
      x = out$fit$par
    )
    meat <- crossprod(score_mat)
    V_robust_raw <- out$vcov %*% meat %*% out$vcov
    robust_fix <- make_psd_matrix(V_robust_raw,
                                  eig_floor_abs = robust_floor_abs,
                                  eig_floor_rel = robust_floor_rel)

    out$score_matrix <- score_mat
    out$vcov_robust_raw <- V_robust_raw
    out$vcov_robust <- robust_fix$matrix_psd
    out$robust_diagnostics <- tibble(
      min_eigenvalue_original = min(robust_fix$eigen_original, na.rm = TRUE),
      min_eigenvalue_fixed = min(robust_fix$eigen_fixed, na.rm = TRUE),
      max_eigenvalue_original = max(robust_fix$eigen_original, na.rm = TRUE),
      clipping_floor = robust_fix$floor_val,
      n_eigenvalues_clipped = robust_fix$n_clipped,
      has_non_finite = any(!is.finite(V_robust_raw))
    )
  }

  out
}

# -----------------------------
# 7. Output tables
# -----------------------------

tidy_joint_msl <- function(object, transformed = FALSE, use_robust = TRUE) {
  if (!inherits(object, "cw2_joint_msl")) stop("object must inherit from 'cw2_joint_msl'.")

  raw_names <- c(
    paste0("beta_", object$prepped$y_names[-1]),
    paste0("lambda_", object$prepped$y_names[-1]),
    paste0("log_sigma_y_", object$prepped$y_names),
    paste0("gamma_", colnames(object$prepped$X_mu)),
    "log_sigma_u",
    "alpha_mu",
    paste0("alpha_", colnames(object$prepped$Z_sel)),
    paste0("delta_", colnames(object$prepped$X_w)),
    "log_sigma_eps",
    "atanh_rho"
  )

  vcov_to_use <- if (use_robust && !is.null(object$vcov_robust)) object$vcov_robust else object$vcov
  raw_tbl <- tibble(parameter = raw_names, estimate = object$fit$par)

  if (!is.null(vcov_to_use)) {
    se <- sqrt(pmax(diag(vcov_to_use), 0))
    raw_tbl <- raw_tbl %>%
      mutate(
        std_error = se,
        z_value = if_else(std_error > 0, estimate / std_error, NA_real_),
        p_value = if_else(!is.na(z_value), 2 * (1 - pnorm(abs(z_value))), NA_real_)
      )
  } else {
    raw_tbl <- raw_tbl %>% mutate(std_error = NA_real_, z_value = NA_real_, p_value = NA_real_)
  }

  if (!transformed) return(raw_tbl)

  trans_tbl <- raw_tbl
  is_log_sigma <- grepl("^log_sigma_", trans_tbl$parameter)
  is_rho <- trans_tbl$parameter == "atanh_rho"

  raw_est <- trans_tbl$estimate
  raw_se <- trans_tbl$std_error

  trans_tbl$estimate[is_log_sigma] <- exp(raw_est[is_log_sigma])
  trans_tbl$estimate[is_rho] <- 0.999 * tanh(raw_est[is_rho])

  trans_tbl$std_error[is_log_sigma] <- exp(raw_est[is_log_sigma]) * raw_se[is_log_sigma]
  trans_tbl$std_error[is_rho] <- 0.999 * (1 - tanh(raw_est[is_rho])^2) * raw_se[is_rho]

  trans_tbl$parameter[is_log_sigma] <- sub("^log_", "", trans_tbl$parameter[is_log_sigma])
  trans_tbl$parameter[is_rho] <- "rho"

  trans_tbl <- trans_tbl %>%
    mutate(
      z_value = if_else(std_error > 0, estimate / std_error, NA_real_),
      p_value = if_else(!is.na(z_value), 2 * (1 - pnorm(abs(z_value))), NA_real_)
    )

  anchor_tbl <- tibble(
    parameter = c(
      paste0("beta_", object$prepped$y_names[1]),
      paste0("lambda_", object$prepped$y_names[1])
    ),
    estimate = c(0, 1),
    std_error = c(NA_real_, NA_real_),
    z_value = c(NA_real_, NA_real_),
    p_value = c(NA_real_, NA_real_)
  )

  bind_rows(anchor_tbl, trans_tbl)
}

summarize_joint_msl <- function(object) {
  if (!inherits(object, "cw2_joint_msl")) stop("object must inherit from 'cw2_joint_msl'.")
  dims <- make_dim_list(object$prepped)
  p <- unpack_theta_joint(object$fit$par, dims)

  list(
    optimization = tibble(
      convergence = object$fit$convergence,
      neg_loglik = object$fit$value,
      N = nrow(object$prepped$joint_df),
      workers = sum(object$prepped$D == 1L),
      draws = object$R,
      draw_type = object$draw_type,
      seed = object$seed
    ),
    accounting = object$prepped$accounting,
    measurement = tibble(
      indicator = object$prepped$y_names,
      beta = p$beta,
      lambda = p$lambda,
      sigma_y = p$sigma_y
    ),
    bargaining = tibble(
      term = colnames(object$prepped$X_mu),
      gamma = p$gamma
    ),
    participation = tibble(
      term = c("mu", colnames(object$prepped$Z_sel)),
      estimate = c(p$alpha_mu, p$alpha_sel)
    ),
    wage = tibble(
      term = colnames(object$prepped$X_w),
      estimate = p$delta
    ),
    selection_scales = tibble(
      parameter = c("sigma_u", "sigma_eps", "rho"),
      estimate = c(p$sigma_u, p$sigma_eps, p$rho)
    ),
    hessian_diagnostics = object$hessian_diagnostics,
    robust_diagnostics = object$robust_diagnostics
  )
}

posterior_mu_mean <- function(object) {
  if (!inherits(object, "cw2_joint_msl")) stop("object must inherit from 'cw2_joint_msl'.")
  dims <- make_dim_list(object$prepped)
  p <- unpack_theta_joint(object$fit$par, dims)

  N <- nrow(object$prepped$Y)
  R <- length(object$draws)
  J <- ncol(object$prepped$Y)

  draw_mat <- repeat_draws_by_row(object$draws, N)
  mu_mean <- as.numeric(object$prepped$X_mu %*% p$gamma)
  MU <- matrix(mu_mean, nrow = N, ncol = R) + p$sigma_u * draw_mat

  log_weights <- matrix(0, nrow = N, ncol = R)
  for (j in seq_len(J)) {
    yj <- matrix(object$prepped$Y[, j], nrow = N, ncol = R)
    mean_j <- p$beta[j] + p$lambda[j] * MU
    log_weights <- log_weights + dnorm(yj, mean = mean_j, sd = p$sigma_y[j], log = TRUE)
  }

  sel_base <- as.numeric(object$prepped$Z_sel %*% p$alpha_sel)
  sel_index <- matrix(sel_base, nrow = N, ncol = R) + p$alpha_mu * MU
  nonworkers <- which(object$prepped$D == 0L)
  workers <- object$prepped$workers

  if (length(nonworkers) > 0) {
    log_weights[nonworkers, ] <- log_weights[nonworkers, ] + pnorm(
      sel_index[nonworkers, , drop = FALSE],
      lower.tail = FALSE,
      log.p = TRUE
    )
  }

  if (length(workers) > 0) {
    w_mean <- as.numeric(object$prepped$X_w[workers, , drop = FALSE] %*% p$delta)
    v <- (object$prepped$W[workers] - w_mean) / p$sigma_eps
    v_mat <- matrix(v, nrow = length(workers), ncol = R)
    denom <- sqrt(pmax(1 - p$rho^2, 1e-8))
    cond_arg <- (sel_index[workers, , drop = FALSE] + p$rho * v_mat) / denom

    log_weights[workers, ] <- log_weights[workers, ] +
      dnorm(v_mat, mean = 0, sd = 1, log = TRUE) -
      log(p$sigma_eps) +
      pnorm(cond_arg, log.p = TRUE)
  }

  log_den <- matrix(row_log_mean_exp(log_weights), nrow = N, ncol = R)
  weights <- exp(log_weights - log_den)
  rowMeans(weights * MU)
}

append_posterior_mu <- function(object) {
  object$prepped$joint_df %>% mutate(mu_post = posterior_mu_mean(object))
}

extract_best_theta <- function(result_object) {
  if (!is.list(result_object) || is.null(result_object$fit) || is.null(result_object$fit$fit$par)) {
    stop("result_object must be the list returned by run_cw2_submission_main().")
  }
  result_object$fit$fit$par
}

save_best_theta <- function(result_object, path = "cw2_best_theta_submission.rds") {
  out <- list(
    theta = extract_best_theta(result_object),
    settings = result_object$fit$prepped$settings
  )
  saveRDS(out, path)
  invisible(path)
}

build_submission_output <- function(fit_object, use_robust = TRUE) {
  list(
    fit = fit_object,
    summary = summarize_joint_msl(fit_object),
    coefficients_raw = tidy_joint_msl(fit_object, transformed = FALSE, use_robust = FALSE),
    coefficients_structural = tidy_joint_msl(fit_object, transformed = TRUE, use_robust = use_robust),
    bargaining_index = append_posterior_mu(fit_object)
  )
}

# -----------------------------
# 8. Main submission wrapper
# -----------------------------

run_cw2_submission_main <- function(data_path = "prepdata.csv",
                                    theta_start_path = NULL,
                                    include_lr_inc = FALSE,
                                    use_logit_sm = TRUE,
                                    anchor_indicator = "s_m_input",
                                    R = 500,
                                    seed = 125,
                                    trace = 1,
                                    hessian_method = "optimHess",
                                    use_robust = TRUE,
                                    eig_floor_abs = 1e-6,
                                    eig_floor_rel = 1e-8,
                                    robust_floor_abs = 1e-10,
                                    robust_floor_rel = 1e-10) {
  fit0 <- fit_joint_msl_submission(
    data_path = data_path,
    include_lr_inc = include_lr_inc,
    use_logit_sm = use_logit_sm,
    anchor_indicator = anchor_indicator,
    R = R,
    seed = seed,
    theta_start_path = theta_start_path,
    trace = trace
  )

  fit1 <- compute_inference_joint_msl(
    fit0,
    hessian_method = hessian_method,
    use_robust = use_robust,
    eig_floor_abs = eig_floor_abs,
    eig_floor_rel = eig_floor_rel,
    robust_floor_abs = robust_floor_abs,
    robust_floor_rel = robust_floor_rel
  )

  build_submission_output(fit1, use_robust = use_robust)
}

# -----------------------------
# 9. Recommended run for submission
# -----------------------------
# Install once if needed:
# install.packages(c("dplyr", "readr", "tibble", "numDeriv"))
# Optional for better rho initialisation:
# install.packages("sampleSelection")
#
# Recommended main run from scratch:
# main_submission <- run_cw2_submission_main(
#   data_path = "prepdata.csv",
#   theta_start_path = NULL,
#   include_lr_inc = FALSE,
#   use_logit_sm = TRUE,
#   anchor_indicator = "s_m_input",
#   R = 500,
#   seed = 125,
#   trace = 1,
#   hessian_method = "optimHess",
#   use_robust = TRUE,
#   eig_floor_abs = 1e-6,
#   eig_floor_rel = 1e-8,
#   robust_floor_abs = 1e-10,
#   robust_floor_rel = 1e-10
# )
#
# If you already have a saved best theta from the same specification:
# main_submission <- run_cw2_submission_main(
#   data_path = "prepdata.csv",
#   theta_start_path = "cw2_best_theta_submission.rds",
#   include_lr_inc = FALSE,
#   use_logit_sm = TRUE,
#   anchor_indicator = "s_m_input",
#   R = 500,
#   seed = 125,
#   trace = 1,
#   hessian_method = "optimHess",
#   use_robust = TRUE
# )
#
# Print core outputs:
# print(main_submission$summary$optimization)
# print(main_submission$summary$accounting, n = Inf)
# print(main_submission$summary$hessian_diagnostics)
# print(main_submission$summary$robust_diagnostics)
# print(main_submission$summary$measurement, n = Inf)
# print(main_submission$summary$participation, n = Inf)
# print(main_submission$summary$wage, n = Inf)
# print(main_submission$summary$selection_scales, n = Inf)
# print(main_submission$coefficients_structural, n = Inf)
# print(cor(main_submission$bargaining_index$mu_post,
#           main_submission$bargaining_index$s_m,
#           use = "complete.obs"))
# print(summary(main_submission$bargaining_index$mu_post))
#
# Save outputs:
# saveRDS(main_submission, "cw2_submission_main_results.rds")
# write.csv(main_submission$coefficients_structural,
#           "cw2_submission_main_coefficients.csv",
#           row.names = FALSE)
# write.csv(main_submission$bargaining_index,
#           "cw2_submission_main_bargaining_index.csv",
#           row.names = FALSE)
#
# Save theta for fast same-spec reruns:
# save_best_theta(main_submission, "cw2_best_theta_submission.rds")
