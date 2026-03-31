# =========================================================
# Coursework 2: Full joint MSL for bargaining, participation, and wages
# Rewrites the two-step / reduced-MSL approach into a single joint model.
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

# -----------------------------
# 0. Utilities
# -----------------------------

clip01 <- function(x, eps = 1e-6) {
  pmin(pmax(x, eps), 1 - eps)
}

row_log_mean_exp <- function(log_mat) {
  m <- apply(log_mat, 1L, max)
  m + log(rowMeans(exp(log_mat - m)))
}

safe_inverse <- function(M, ridge = 1e-8) {
  out <- tryCatch(solve(M), error = function(e) NULL)
  if (!is.null(out)) return(out)
  solve(M + diag(ridge, nrow(M)))
}

make_antithetic_draws <- function(N, R, seed = 123) {
  if (R < 2) stop("R must be at least 2.")
  half_R <- ceiling(R / 2)
  set.seed(seed)
  U <- matrix(rnorm(N * half_R), nrow = N, ncol = half_R)
  draws <- cbind(U, -U)
  draws[, seq_len(R), drop = FALSE]
}

# -----------------------------
# 1. Variable lists
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
  z_vars_mu  <- if (include_lr_inc) z_vars_all else setdiff(z_vars_all, "lr_inc")

  list(
    measurement_vars = measurement_vars,
    trait_vars = trait_vars,
    z_vars_mu = z_vars_mu,
    z_vars_all = z_vars_all
  )
}

# -----------------------------
# 2. Data preparation
# -----------------------------

build_joint_sample <- function(data_path = "prepdata.csv",
                               include_lr_inc = FALSE,
                               use_logit_sm = FALSE) {

  vars <- get_variable_lists(include_lr_inc = include_lr_inc)

  df <- readr::read_csv(data_path, show_col_types = FALSE)

  if (!("hid" %in% names(df))) {
    df <- dplyr::mutate(df, hid = dplyr::row_number())
  }

  required_vars <- c(
    "hid", "s_m", "age.f", "educ.f", "income.f",
    vars$measurement_vars, vars$trait_vars, vars$z_vars_all
  )

  missing_cols <- setdiff(required_vars, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- df %>%
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

  # Joint MSL sample:
  # - full measurement block observed
  # - full bargaining covariates observed
  # - full participation covariates observed
  # - if D_f = 1, W_f must be observed
  joint_df <- df %>%
    select(
      hid,
      s_m, s_m_input,
      all_of(vars$measurement_vars),
      all_of(vars$trait_vars),
      all_of(vars$z_vars_all),
      age.f, educ.f,
      D_f, W_f
    ) %>%
    filter(
      if_all(
        c("s_m_input", vars$measurement_vars,
          vars$trait_vars, vars$z_vars_mu,
          "age.f", "educ.f", "D_f"),
        ~ !is.na(.)
      )
    ) %>%
    filter(D_f == 0L | (D_f == 1L & !is.na(W_f)))

  y_names <- c("s_m_input", vars$measurement_vars)

  Y_raw <- as.matrix(joint_df[, y_names])
  Y <- scale(Y_raw)
  colnames(Y) <- y_names

  X_mu_df <- joint_df %>% select(all_of(c(vars$trait_vars, vars$z_vars_mu)))
  X_mu <- model.matrix(~ ., data = X_mu_df)

  Z_sel_df <- joint_df %>%
    select(age.f, educ.f, duration, abouttobeparents, d_educ, d_age, both_m)
  Z_sel <- model.matrix(~ ., data = Z_sel_df)

  X_w_df <- joint_df %>% select(age.f, educ.f)
  X_w <- model.matrix(~ ., data = X_w_df)

  workers <- which(joint_df$D_f == 1L)

  list(
    raw_df = df,
    joint_df = joint_df,
    Y = Y,
    Y_raw = Y_raw,
    y_names = y_names,
    X_mu = X_mu,
    Z_sel = Z_sel,
    X_w = X_w,
    D = as.integer(joint_df$D_f),
    W = joint_df$W_f,
    workers = workers,
    vars = vars,
    settings = list(
      include_lr_inc = include_lr_inc,
      use_logit_sm = use_logit_sm
    )
  )
}

# -----------------------------
# 3. Parameter mapping
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

  beta_free <- theta[idx:(idx + J - 2L)]
  idx <- idx + J - 1L

  lambda_free <- theta[idx:(idx + J - 2L)]
  idx <- idx + J - 1L

  log_sigma_y <- theta[idx:(idx + J - 1L)]
  idx <- idx + J

  gamma <- theta[idx:(idx + K_mu - 1L)]
  idx <- idx + K_mu

  log_sigma_u <- theta[idx]
  idx <- idx + 1L

  alpha_mu <- theta[idx]
  idx <- idx + 1L

  alpha_sel <- theta[idx:(idx + K_sel - 1L)]
  idx <- idx + K_sel

  delta <- theta[idx:(idx + K_w - 1L)]
  idx <- idx + K_w

  log_sigma_eps <- theta[idx]
  idx <- idx + 1L

  atanh_rho <- theta[idx]

  beta <- c(0, beta_free)         # anchor intercept fixed at 0
  lambda <- c(1, lambda_free)     # anchor loading fixed at 1
  sigma_y <- exp(log_sigma_y)
  sigma_u <- exp(log_sigma_u)
  sigma_eps <- exp(log_sigma_eps)
  rho <- tanh(atanh_rho)

  list(
    beta = beta,
    lambda = lambda,
    sigma_y = sigma_y,
    gamma = gamma,
    sigma_u = sigma_u,
    alpha_mu = alpha_mu,
    alpha_sel = alpha_sel,
    delta = delta,
    sigma_eps = sigma_eps,
    rho = rho
  )
}

# -----------------------------
# 4. Initial values
# -----------------------------

get_initial_values_joint <- function(prepped) {

  dims <- make_dim_list(prepped)
  J <- dims$J

  # Measurement starts from the user's current factor-analysis workflow
  fa_start <- factanal(
    x = prepped$Y,
    factors = 1,
    scores = "regression"
  )

  load0 <- as.numeric(fa_start$loadings[, 1])
  uniq0 <- pmax(as.numeric(fa_start$uniquenesses), 0.05)
  beta0 <- colMeans(prepped$Y)
  mu0 <- as.numeric(fa_start$scores[, 1])

  # orient factor so the anchor loads positively
  if (load0[1] < 0) {
    load0 <- -load0
    mu0 <- -mu0
  }

  # rescale so anchor loading is exactly 1
  if (abs(load0[1]) < 1e-4) load0[1] <- 1
  mu0 <- mu0 * load0[1]
  lambda0 <- load0 / load0[1]
  lambda0[1] <- 1

  beta_free0 <- beta0[-1]
  lambda_free0 <- lambda0[-1]
  log_sigma_y0 <- log(sqrt(uniq0))

  # Bargaining equation starts from OLS on factor scores
  gamma_fit <- lm(mu0 ~ prepped$X_mu - 1)
  gamma0 <- coef(gamma_fit)
  gamma0[is.na(gamma0)] <- 0

  mu_mean0 <- as.numeric(prepped$X_mu %*% gamma0)
  sigma_u0 <- sd(mu0 - mu_mean0)
  if (!is.finite(sigma_u0) || sigma_u0 <= 0) sigma_u0 <- 1
  log_sigma_u0 <- log(sigma_u0)

  # Participation starts from probit using generated mu only as an initializer
  init_df <- prepped$joint_df %>%
    mutate(mu0 = mu0)

  part_fit <- glm(
    D_f ~ mu0 + age.f + educ.f + duration + abouttobeparents + d_educ + d_age + both_m,
    family = binomial(link = "probit"),
    data = init_df
  )

  part_coef <- coef(part_fit)
  part_coef[is.na(part_coef)] <- 0

  alpha_mu0 <- unname(part_coef["mu0"])
  alpha_sel0 <- rep(0, ncol(prepped$Z_sel))
  names(alpha_sel0) <- colnames(prepped$Z_sel)
  overlap_sel <- intersect(names(part_coef), names(alpha_sel0))
  alpha_sel0[overlap_sel] <- part_coef[overlap_sel]

  # Wage starts from worker-only OLS
  wage_fit <- lm(
    W_f ~ age.f + educ.f,
    data = init_df %>% filter(D_f == 1L)
  )

  wage_coef <- coef(wage_fit)
  wage_coef[is.na(wage_coef)] <- 0

  delta0 <- rep(0, ncol(prepped$X_w))
  names(delta0) <- colnames(prepped$X_w)
  overlap_w <- intersect(names(wage_coef), names(delta0))
  delta0[overlap_w] <- wage_coef[overlap_w]

  sigma_eps0 <- summary(wage_fit)$sigma
  if (!is.finite(sigma_eps0) || sigma_eps0 <= 0) sigma_eps0 <- 1
  log_sigma_eps0 <- log(sigma_eps0)

  # Optional Heckman-based rho initializer
  rho0 <- 0
  if (requireNamespace("sampleSelection", quietly = TRUE)) {
    heckman_try <- try(
      sampleSelection::selection(
        selection = D_f ~ mu0 + age.f + educ.f + duration + abouttobeparents + d_educ + d_age + both_m,
        outcome   = W_f ~ age.f + educ.f,
        data = init_df,
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
  atanh_rho0 <- atanh(rho0)

  theta0 <- pack_theta_joint(
    beta_free = beta_free0,
    lambda_free = lambda_free0,
    log_sigma_y = log_sigma_y0,
    gamma = gamma0,
    log_sigma_u = log_sigma_u0,
    alpha_mu = alpha_mu0,
    alpha_sel = alpha_sel0,
    delta = delta0,
    log_sigma_eps = log_sigma_eps0,
    atanh_rho = atanh_rho0
  )

  list(
    theta = theta0,
    mu0 = mu0,
    dims = dims
  )
}

# -----------------------------
# 5. Joint simulated likelihood
# -----------------------------

joint_loglik_obs <- function(theta, prepped, draws) {

  dims <- make_dim_list(prepped)
  p <- unpack_theta_joint(theta, dims)

  N <- nrow(prepped$Y)
  R <- ncol(draws)
  J <- ncol(prepped$Y)

  mu_mean <- as.numeric(prepped$X_mu %*% p$gamma)
  MU <- matrix(mu_mean, nrow = N, ncol = R) + p$sigma_u * draws

  # Measurement block
  log_meas <- matrix(0, nrow = N, ncol = R)

  for (j in seq_len(J)) {
    yj <- matrix(prepped$Y[, j], nrow = N, ncol = R)
    mean_j <- p$beta[j] + p$lambda[j] * MU
    log_meas <- log_meas + dnorm(yj, mean = mean_j, sd = p$sigma_y[j], log = TRUE)
  }

  # Labour outcome block
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
    cond_arg <- (sel_index[workers, , drop = FALSE] + p$rho * v_mat) /
      sqrt(1 - p$rho^2)

    log_dw[workers, ] <-
      dnorm(v_mat, mean = 0, sd = 1, log = TRUE) -
      log(p$sigma_eps) +
      pnorm(cond_arg, log.p = TRUE)
  }

  row_log_mean_exp(log_meas + log_dw)
}

objective_joint_msl <- function(theta, prepped, draws) {
  -sum(joint_loglik_obs(theta, prepped = prepped, draws = draws))
}

# -----------------------------
# 6. Estimation wrapper
# -----------------------------

fit_joint_msl <- function(data_path = "prepdata.csv",
                          include_lr_inc = FALSE,
                          use_logit_sm = FALSE,
                          R = 200,
                          seed = 123,
                          nm_maxit = 300,
                          bfgs_maxit = 600,
                          trace = 1,
                          compute_hessian = FALSE) {

  prepped <- build_joint_sample(
    data_path = data_path,
    include_lr_inc = include_lr_inc,
    use_logit_sm = use_logit_sm
  )

  init <- get_initial_values_joint(prepped)
  draws <- make_antithetic_draws(nrow(prepped$Y), R = R, seed = seed)

  cat("\n============================\n")
  cat("Joint MSL sample summary\n")
  cat("============================\n")
  cat("N =", nrow(prepped$joint_df), "\n")
  cat("Workers =", sum(prepped$D == 1L), "\n")
  cat("Indicators =", ncol(prepped$Y), "\n")
  cat("R draws =", R, "\n")
  cat("include_lr_inc =", include_lr_inc, "\n")
  cat("use_logit_sm =", use_logit_sm, "\n")

  cat("\nStarting Nelder-Mead...\n")
  nm_fit <- optim(
    par = init$theta,
    fn = objective_joint_msl,
    prepped = prepped,
    draws = draws,
    method = "Nelder-Mead",
    control = list(maxit = nm_maxit, trace = trace, REPORT = 10)
  )

  cat("\nStarting BFGS...\n")
  bfgs_fit <- optim(
    par = nm_fit$par,
    fn = objective_joint_msl,
    prepped = prepped,
    draws = draws,
    method = "BFGS",
    hessian = compute_hessian,
    control = list(maxit = bfgs_maxit, trace = trace, REPORT = 10)
  )

  out <- list(
    prepped = prepped,
    init = init,
    draws = draws,
    nm_fit = nm_fit,
    fit = bfgs_fit
  )
  class(out) <- "cw2_joint_msl"
  out
}

# -----------------------------
# 7. Post-estimation helpers
# -----------------------------

compute_vcov_joint_msl <- function(object) {
  if (!inherits(object, "cw2_joint_msl")) stop("object must inherit from 'cw2_joint_msl'.")

  H <- object$fit$hessian
  if (is.null(H) || any(!is.finite(H))) {
    H <- optimHess(
      par = object$fit$par,
      fn = objective_joint_msl,
      prepped = object$prepped,
      draws = object$draws
    )
  }

  object$hessian <- H
  object$vcov <- safe_inverse(H)
  object
}

tidy_joint_msl <- function(object, transformed = FALSE) {
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

  raw_tbl <- tibble(
    parameter = raw_names,
    estimate = object$fit$par
  )

  if (!is.null(object$vcov)) {
    se <- sqrt(diag(object$vcov))
    raw_tbl <- raw_tbl %>%
      mutate(
        std_error = se,
        z_value = estimate / std_error,
        p_value = 2 * (1 - pnorm(abs(z_value)))
      )
  } else {
    raw_tbl <- raw_tbl %>%
      mutate(
        std_error = NA_real_,
        z_value = NA_real_,
        p_value = NA_real_
      )
  }

  if (!transformed) return(raw_tbl)

  trans_tbl <- raw_tbl %>%
    mutate(
      estimate = case_when(
        grepl("^log_sigma_", parameter) ~ exp(estimate),
        parameter == "atanh_rho" ~ tanh(estimate),
        TRUE ~ estimate
      ),
      std_error = case_when(
        grepl("^log_sigma_", parameter) & !is.na(std_error) ~ estimate * std_error,
        parameter == "atanh_rho" & !is.na(std_error) ~ (1 - estimate^2) * std_error,
        TRUE ~ std_error
      ),
      parameter = case_when(
        grepl("^log_sigma_", parameter) ~ sub("^log_", "", parameter),
        parameter == "atanh_rho" ~ "rho",
        TRUE ~ parameter
      )
    ) %>%
    mutate(
      z_value = if_else(!is.na(std_error), estimate / std_error, NA_real_),
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
      draws = ncol(object$draws)
    ),
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
    )
  )
}

posterior_mu_mean <- function(object) {
  if (!inherits(object, "cw2_joint_msl")) stop("object must inherit from 'cw2_joint_msl'.")

  dims <- make_dim_list(object$prepped)
  p <- unpack_theta_joint(object$fit$par, dims)

  N <- nrow(object$prepped$Y)
  R <- ncol(object$draws)

  mu_mean <- as.numeric(object$prepped$X_mu %*% p$gamma)
  MU <- matrix(mu_mean, nrow = N, ncol = R) + p$sigma_u * object$draws

  # Un-normalised posterior weights proportional to measurement * labour block
  log_weights <- matrix(0, nrow = N, ncol = R)

  J <- ncol(object$prepped$Y)
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
    cond_arg <- (sel_index[workers, , drop = FALSE] + p$rho * v_mat) /
      sqrt(1 - p$rho^2)

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
  object$prepped$joint_df %>%
    mutate(mu_post = posterior_mu_mean(object))
}

# -----------------------------
# 8. Example usage
# -----------------------------
# fit <- fit_joint_msl(
#   data_path = "prepdata.csv",
#   include_lr_inc = FALSE,   # critique-consistent default
#   use_logit_sm = FALSE,     # set TRUE if you want a logit transform of s_m
#   R = 200,
#   seed = 123,
#   nm_maxit = 300,
#   bfgs_maxit = 600,
#   trace = 1,
#   compute_hessian = FALSE
# )
#
# fit <- compute_vcov_joint_msl(fit)
# summary_list <- summarize_joint_msl(fit)
# raw_table <- tidy_joint_msl(fit, transformed = FALSE)
# structural_table <- tidy_joint_msl(fit, transformed = TRUE)
# bargaining_index_df <- append_posterior_mu(fit)
#
# print(summary_list$optimization)
# print(summary_list$measurement, n = Inf)
# print(summary_list$participation, n = Inf)
# print(summary_list$wage, n = Inf)
