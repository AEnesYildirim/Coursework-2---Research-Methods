# =========================================================
# Coursework 2: Main joint-MSL estimator for bargaining,
# participation, and wages with sample selection
#
# Main-specification defaults:
#   - anchor latent factor on the observed consumption-share proxy
#   - logit-transform s_m to respect its (0,1) support
#   - keep the anchor in transformed units; standardise only the
#     non-anchor measurement indicators
#   - standardise continuous structural covariates to stabilise MSL
#   - estimate the measurement block, bargaining equation,
#     participation equation, and wage-selection block jointly
#   - use fixed common simulation draws (randomised Halton by default)
#   - use staged continuation in R and multi-start optimisation
#
# This file is intended to be the single baseline script for the
# coursework main result. Robustness variants can be layered on top.
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
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

  warning(
    "Hessian inversion failed; applying ridge regularisation. ",
    "Standard errors based on the regularised inverse should be interpreted cautiously."
  )
  solve(M + diag(ridge, nrow(M)))
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

make_common_draws <- function(R,
                              draw_type = c("halton", "antithetic"),
                              seed = 123,
                              eps = 1e-10) {
  draw_type <- match.arg(draw_type)

  if (R < 2L) stop("R must be at least 2.")

  if (draw_type == "halton") {
    set.seed(seed)
    shift <- runif(1)
    u <- (halton_1d(R, base = 2L, start_index = 11L) + shift) %% 1
    u <- clip01(u, eps = eps)
    return(qnorm(u))
  }

  half_R <- ceiling(R / 2)
  set.seed(seed)
  z <- rnorm(half_R)
  draws <- c(z, -z)
  draws[seq_len(R)]
}

repeat_draws_by_row <- function(draws, N) {
  matrix(draws, nrow = N, ncol = length(draws), byrow = TRUE)
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

is_binary_like <- function(x) {
  x_no_na <- unique(stats::na.omit(x))
  length(x_no_na) <= 2L && all(x_no_na %in% c(0, 1))
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

normalise_R_schedule <- function(R = 500, R_schedule = NULL) {
  if (!is.null(R_schedule)) {
    out <- sort(unique(as.integer(R_schedule)))
    out <- out[out >= 2L]
    if (length(out) == 0L) stop("R_schedule must contain at least one integer >= 2.")
    return(out)
  }

  if (length(R) > 1L) {
    out <- sort(unique(as.integer(R)))
    out <- out[out >= 2L]
    if (length(out) == 0L) stop("R must contain at least one integer >= 2.")
    return(out)
  }

  R <- as.integer(R[1])
  if (R < 2L) stop("R must be at least 2.")
  if (R <= 100L) return(R)
  if (R <= 300L) return(c(100L, R))
  c(100L, 300L, R)
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
  z_vars_mu <- if (include_lr_inc) z_vars_all else setdiff(z_vars_all, "lr_inc")

  list(
    measurement_vars = measurement_vars,
    trait_vars = trait_vars,
    z_vars_mu = z_vars_mu,
    z_vars_all = z_vars_all
  )
}

sample_accounting <- function(df, vars, use_logit_sm = TRUE) {
  tmp <- df %>%
    mutate(
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

  base_y_names <- c("s_m_input", vars$measurement_vars)

  keep0 <- rep(TRUE, nrow(tmp))
  keep1 <- keep0 & complete.cases(tmp[, base_y_names, drop = FALSE])
  keep2 <- keep1 & complete.cases(tmp[, c(vars$trait_vars, vars$z_vars_mu), drop = FALSE])
  keep3 <- keep2 & complete.cases(tmp[, c("age.f", "educ.f", "D_f"), drop = FALSE])
  keep4 <- keep3 & (tmp$D_f == 0L | (tmp$D_f == 1L & !is.na(tmp$W_f)))

  tibble(
    stage = c(
      "Raw rows",
      "After measurement block completeness",
      "After bargaining-covariate completeness",
      "After labour-covariate completeness",
      "Final joint-MSL sample"
    ),
    n = c(sum(keep0), sum(keep1), sum(keep2), sum(keep3), sum(keep4))
  )
}

# -----------------------------
# 2. Data preparation
# -----------------------------

build_joint_sample <- function(data_path = "prepdata.csv",
                               include_lr_inc = FALSE,
                               use_logit_sm = TRUE,
                               anchor_indicator = "s_m_input") {
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

  accounting_tbl <- sample_accounting(df, vars, use_logit_sm = use_logit_sm)

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

  base_y_names <- c("s_m_input", vars$measurement_vars)
  if (!(anchor_indicator %in% base_y_names)) {
    stop("anchor_indicator must be one of: ", paste(base_y_names, collapse = ", "))
  }
  y_names <- c(anchor_indicator, setdiff(base_y_names, anchor_indicator))

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
      if_all(c(base_y_names, vars$trait_vars, vars$z_vars_mu, "age.f", "educ.f", "D_f"), ~ !is.na(.))
    ) %>%
    filter(D_f == 0L | (D_f == 1L & !is.na(W_f)))

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

  # Exclusion restrictions: duration, abouttobeparents, d_educ, d_age,
  # and both_m enter selection but not the wage equation. Defend this in the report.
  Z_sel_df_raw <- joint_df %>%
    select(age.f, educ.f, duration, abouttobeparents, d_educ, d_age, both_m)
  Z_sel_std <- standardise_continuous_df(
    Z_sel_df_raw,
    binary_names = c("abouttobeparents", "both_m")
  )
  Z_sel <- model.matrix(~ ., data = Z_sel_std$data)

  X_w_df_raw <- joint_df %>% select(age.f, educ.f)
  X_w_std <- standardise_continuous_df(X_w_df_raw)
  X_w <- model.matrix(~ ., data = X_w_std$data)

  workers <- which(joint_df$D_f == 1L)

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
    X_mu_data = X_mu_std$data,
    X_mu_scaling = list(
      centers = X_mu_std$centers,
      scales = X_mu_std$scales,
      standardized = X_mu_std$standardized
    ),
    Z_sel = Z_sel,
    Z_sel_data = Z_sel_std$data,
    Z_sel_scaling = list(
      centers = Z_sel_std$centers,
      scales = Z_sel_std$scales,
      standardized = Z_sel_std$standardized
    ),
    X_w = X_w,
    X_w_data = X_w_std$data,
    X_w_scaling = list(
      centers = X_w_std$centers,
      scales = X_w_std$scales,
      standardized = X_w_std$standardized
    ),
    D = as.integer(joint_df$D_f),
    W = joint_df$W_f,
    workers = workers,
    vars = vars,
    settings = list(
      include_lr_inc = include_lr_inc,
      use_logit_sm = use_logit_sm,
      anchor_indicator = anchor_indicator
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

  sigma_y <- pmax(exp(log_sigma_y), 1e-6)
  sigma_u <- pmax(exp(log_sigma_u), 1e-6)
  sigma_eps <- pmax(exp(log_sigma_eps), 1e-6)
  rho <- 0.999 * tanh(atanh_rho)

  list(
    beta = c(0, beta_free),
    lambda = c(1, lambda_free),
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

  part_df <- dplyr::bind_cols(
    tibble(D_f = prepped$D, mu0 = mu0),
    as.data.frame(prepped$Z_sel[, -1, drop = FALSE])
  )
  part_fit <- glm(
    D_f ~ mu0 + .,
    family = binomial(link = "probit"),
    data = part_df
  )
  part_coef <- coef(part_fit)
  part_coef[is.na(part_coef)] <- 0

  alpha_mu0 <- unname(part_coef["mu0"])
  alpha_sel0 <- rep(0, ncol(prepped$Z_sel))
  names(alpha_sel0) <- colnames(prepped$Z_sel)
  overlap_sel <- intersect(names(part_coef), names(alpha_sel0))
  alpha_sel0[overlap_sel] <- part_coef[overlap_sel]

  wage_df <- dplyr::bind_cols(
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
        data = dplyr::bind_cols(
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
    mu0 = mu0,
    dims = make_dim_list(prepped)
  )
}

# -----------------------------
# 5. Joint simulated likelihood
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
# 6. Estimation wrapper
# -----------------------------

make_starting_values <- function(theta0,
                                 n_starts = 5,
                                 jitter_scale = 0.05,
                                 seed = 456) {
  starts <- vector("list", n_starts)
  starts[[1]] <- theta0

  if (n_starts <= 1L) return(starts)

  set.seed(seed)
  for (s in 2:n_starts) {
    shock <- rnorm(length(theta0), mean = 0, sd = jitter_scale)
    starts[[s]] <- theta0 + shock
  }
  starts
}

run_one_start <- function(start_par,
                          prepped,
                          draws,
                          nm_maxit,
                          bfgs_maxit,
                          trace,
                          compute_hessian) {
  nm_fit <- optim(
    par = start_par,
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
    hessian = compute_hessian,
    control = list(maxit = bfgs_maxit, trace = trace, REPORT = 25)
  )

  list(
    nm_fit = nm_fit,
    fit = bfgs_fit
  )
}

run_multistart_block <- function(start_list,
                                 prepped,
                                 draws,
                                 nm_maxit,
                                 bfgs_maxit,
                                 trace,
                                 compute_hessian,
                                 parallel_starts = FALSE) {
  runner <- function(s) {
    fit_obj <- run_one_start(
      start_par = start_list[[s]],
      prepped = prepped,
      draws = draws,
      nm_maxit = nm_maxit,
      bfgs_maxit = bfgs_maxit,
      trace = trace,
      compute_hessian = compute_hessian
    )

    list(
      fit_obj = fit_obj,
      summary = tibble(
        start_id = s,
        nm_value = fit_obj$nm_fit$value,
        nm_conv = fit_obj$nm_fit$convergence,
        bfgs_value = fit_obj$fit$value,
        bfgs_conv = fit_obj$fit$convergence
      )
    )
  }

  if (parallel_starts && requireNamespace("future.apply", quietly = TRUE)) {
    fits <- future.apply::future_lapply(seq_along(start_list), runner, future.seed = TRUE)
  } else {
    fits <- lapply(seq_along(start_list), runner)
  }

  fit_summaries <- bind_rows(lapply(fits, `[[`, "summary"))
  fit_objects <- lapply(fits, `[[`, "fit_obj")
  best_idx <- which.min(fit_summaries$bfgs_value)

  list(
    fits = fit_objects,
    fit_summaries = fit_summaries,
    best_idx = best_idx,
    best_fit = fit_objects[[best_idx]]
  )
}

fit_joint_msl <- function(data_path = "prepdata.csv",
                          include_lr_inc = FALSE,
                          use_logit_sm = TRUE,
                          anchor_indicator = "s_m_input",
                          R = 500,
                          R_schedule = NULL,
                          draw_type = c("halton", "antithetic"),
                          seed = 123,
                          n_starts = 5,
                          start_jitter = 0.05,
                          start_seed = 456,
                          nm_maxit = 300,
                          bfgs_maxit = 600,
                          trace = 1,
                          compute_hessian = FALSE,
                          parallel_starts = FALSE) {
  draw_type <- match.arg(draw_type)
  draw_schedule <- normalise_R_schedule(R = R, R_schedule = R_schedule)

  prepped <- build_joint_sample(
    data_path = data_path,
    include_lr_inc = include_lr_inc,
    use_logit_sm = use_logit_sm,
    anchor_indicator = anchor_indicator
  )

  init <- get_initial_values_joint(prepped)
  current_start_list <- make_starting_values(
    theta0 = init$theta,
    n_starts = n_starts,
    jitter_scale = start_jitter,
    seed = start_seed
  )

  stage_results <- vector("list", length(draw_schedule))
  final_draws <- NULL
  final_block <- NULL

  cat("\n============================\n")
  cat("Joint MSL sample summary\n")
  cat("============================\n")
  cat("N =", nrow(prepped$joint_df), "\n")
  cat("Workers =", sum(prepped$D == 1L), "\n")
  cat("Indicators =", ncol(prepped$Y), "\n")
  cat("Anchor =", anchor_indicator, "\n")
  cat("draw_schedule =", paste(draw_schedule, collapse = " -> "), "\n")
  cat("draw_type =", draw_type, "\n")
  cat("include_lr_inc =", include_lr_inc, "\n")
  cat("use_logit_sm =", use_logit_sm, "\n")
  cat("n_starts =", n_starts, "\n")
  cat("parallel_starts =", parallel_starts, "\n")

  for (stage_idx in seq_along(draw_schedule)) {
    R_stage <- draw_schedule[stage_idx]
    stage_seed <- seed + stage_idx - 1L
    final_draws <- make_common_draws(R = R_stage, draw_type = draw_type, seed = stage_seed)

    cat("\nStage", stage_idx, "of", length(draw_schedule), ": R =", R_stage, "\n")

    block <- run_multistart_block(
      start_list = current_start_list,
      prepped = prepped,
      draws = final_draws,
      nm_maxit = nm_maxit,
      bfgs_maxit = bfgs_maxit,
      trace = trace,
      compute_hessian = compute_hessian && stage_idx == length(draw_schedule),
      parallel_starts = parallel_starts
    )

    stage_tbl <- block$fit_summaries %>%
      mutate(stage_id = stage_idx, R = R_stage, stage_seed = stage_seed)

    stage_results[[stage_idx]] <- list(
      R = R_stage,
      draws = final_draws,
      results = block,
      summary = stage_tbl
    )

    final_block <- block

    if (stage_idx < length(draw_schedule)) {
      current_start_list <- make_starting_values(
        theta0 = block$best_fit$fit$par,
        n_starts = n_starts,
        jitter_scale = start_jitter / (stage_idx + 1),
        seed = start_seed + stage_idx
      )
    }
  }

  multi_start_tbl <- bind_rows(lapply(stage_results, `[[`, "summary"))
  final_stage_summary <- stage_results[[length(stage_results)]]$summary
  final_best_idx <- which.min(final_stage_summary$bfgs_value)

  out <- list(
    prepped = prepped,
    init = init,
    draws = final_draws,
    draw_schedule = draw_schedule,
    stage_results = stage_results,
    nm_fit = final_block$best_fit$nm_fit,
    fit = final_block$best_fit$fit,
    multi_start_results = multi_start_tbl,
    all_fits = lapply(stage_results, function(x) x$results$fits),
    best_start_id = final_best_idx,
    best_stage_id = length(stage_results),
    hessian = NULL,
    vcov = NULL,
    vcov_robust = NULL,
    hessian_diagnostics = NULL
  )
  class(out) <- "cw2_joint_msl"
  out
}

# -----------------------------
# 7. Hessian and VCOV
# -----------------------------

compute_hessian_joint_msl <- function(object,
                                      hessian_method = c("optimHess", "numDeriv")) {
  if (!inherits(object, "cw2_joint_msl")) {
    stop("object must inherit from 'cw2_joint_msl'.")
  }

  hessian_method <- match.arg(hessian_method)

  if (hessian_method == "numDeriv") {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      warning("numDeriv not installed; falling back to optimHess.")
      hessian_method <- "optimHess"
    }
  }

  if (hessian_method == "numDeriv") {
    H <- numDeriv::hessian(
      func = objective_joint_msl,
      x = object$fit$par,
      prepped = object$prepped,
      draws = object$draws
    )
  } else {
    H <- optimHess(
      par = object$fit$par,
      fn = objective_joint_msl,
      prepped = object$prepped,
      draws = object$draws
    )
  }

  H
}

compute_vcov_joint_msl <- function(object,
                                   hessian_method = c("optimHess", "numDeriv"),
                                   robust = FALSE) {
  if (!inherits(object, "cw2_joint_msl")) {
    stop("object must inherit from 'cw2_joint_msl'.")
  }

  hessian_method <- match.arg(hessian_method)
  H <- compute_hessian_joint_msl(object, hessian_method = hessian_method)
  vcov_hessian <- safe_inverse(H)

  object$hessian <- H
  object$vcov <- vcov_hessian

  eigvals <- tryCatch(
    eigen(H, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) rep(NA_real_, nrow(H))
  )
  object$hessian_diagnostics <- tibble(
    min_eigenvalue = min(eigvals, na.rm = TRUE),
    max_eigenvalue = max(eigvals, na.rm = TRUE),
    has_non_finite = any(!is.finite(H)),
    ridge_warning = is.null(tryCatch(solve(H), error = function(e) NULL))
  )

  if (robust) {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      warning("numDeriv not installed; robust vcov not computed.")
      return(object)
    }

    score_mat <- numDeriv::jacobian(
      func = function(th) joint_loglik_obs(th, prepped = object$prepped, draws = object$draws),
      x = object$fit$par
    )
    meat <- crossprod(score_mat)
    object$score_matrix <- score_mat
    object$vcov_robust <- vcov_hessian %*% meat %*% vcov_hessian
  }

  object
}

# -----------------------------
# 8. Output helpers
# -----------------------------

tidy_joint_msl <- function(object, transformed = FALSE, use_robust = FALSE) {
  if (!inherits(object, "cw2_joint_msl")) {
    stop("object must inherit from 'cw2_joint_msl'.")
  }

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

  raw_tbl <- tibble(parameter = raw_names, estimate = object$fit$par)

  vcov_to_use <- NULL
  if (use_robust && !is.null(object$vcov_robust)) {
    vcov_to_use <- object$vcov_robust
  } else if (!is.null(object$vcov)) {
    vcov_to_use <- object$vcov
  }

  if (!is.null(vcov_to_use)) {
    se <- sqrt(diag(vcov_to_use))
    raw_tbl <- raw_tbl %>%
      mutate(
        std_error = se,
        z_value = estimate / std_error,
        p_value = 2 * (1 - pnorm(abs(z_value)))
      )
  } else {
    raw_tbl <- raw_tbl %>%
      mutate(std_error = NA_real_, z_value = NA_real_, p_value = NA_real_)
  }

  if (!transformed) return(raw_tbl)

  trans_tbl <- raw_tbl
  is_log_sigma <- grepl("^log_sigma_", trans_tbl$parameter)
  is_rho <- trans_tbl$parameter == "atanh_rho"

  raw_est <- trans_tbl$estimate
  raw_se <- trans_tbl$std_error

  trans_tbl$estimate[is_log_sigma] <- exp(raw_est[is_log_sigma])
  trans_tbl$estimate[is_rho] <- 0.999 * tanh(raw_est[is_rho])

  if (!all(is.na(raw_se))) {
    trans_tbl$std_error[is_log_sigma] <- exp(raw_est[is_log_sigma]) * raw_se[is_log_sigma]
    trans_tbl$std_error[is_rho] <- 0.999 * (1 - tanh(raw_est[is_rho])^2) * raw_se[is_rho]
  }

  trans_tbl$parameter[is_log_sigma] <- sub("^log_", "", trans_tbl$parameter[is_log_sigma])
  trans_tbl$parameter[is_rho] <- "rho"

  trans_tbl <- trans_tbl %>%
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
  if (!inherits(object, "cw2_joint_msl")) {
    stop("object must inherit from 'cw2_joint_msl'.")
  }

  dims <- make_dim_list(object$prepped)
  p <- unpack_theta_joint(object$fit$par, dims)

  list(
    optimization = tibble(
      convergence = object$fit$convergence,
      neg_loglik = object$fit$value,
      N = nrow(object$prepped$joint_df),
      workers = sum(object$prepped$D == 1L),
      draws = length(object$draws),
      draw_schedule = paste(object$draw_schedule, collapse = " -> "),
      best_stage_id = object$best_stage_id,
      best_start_id = object$best_start_id
    ),
    continuation = bind_rows(lapply(object$stage_results, `[[`, "summary")),
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
    X_mu_scaling = tibble(
      variable = names(object$prepped$X_mu_scaling$centers),
      center = unname(object$prepped$X_mu_scaling$centers),
      scale = unname(object$prepped$X_mu_scaling$scales),
      standardized = unname(object$prepped$X_mu_scaling$standardized)
    ),
    Z_sel_scaling = tibble(
      variable = names(object$prepped$Z_sel_scaling$centers),
      center = unname(object$prepped$Z_sel_scaling$centers),
      scale = unname(object$prepped$Z_sel_scaling$scales),
      standardized = unname(object$prepped$Z_sel_scaling$standardized)
    ),
    X_w_scaling = tibble(
      variable = names(object$prepped$X_w_scaling$centers),
      center = unname(object$prepped$X_w_scaling$centers),
      scale = unname(object$prepped$X_w_scaling$scales),
      standardized = unname(object$prepped$X_w_scaling$standardized)
    )
  )
}

posterior_mu_mean <- function(object) {
  if (!inherits(object, "cw2_joint_msl")) {
    stop("object must inherit from 'cw2_joint_msl'.")
  }

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

print.cw2_joint_msl <- function(x, ...) {
  cat("cw2_joint_msl object\n")
  cat("  N =", nrow(x$prepped$joint_df), "\n")
  cat("  Workers =", sum(x$prepped$D == 1L), "\n")
  cat("  neg_loglik =", x$fit$value, "\n")
  cat("  convergence =", x$fit$convergence, "\n")
  cat("  draw_schedule =", paste(x$draw_schedule, collapse = " -> "), "\n")
  cat("  best_stage_id =", x$best_stage_id, "\n")
  cat("  best_start_id =", x$best_start_id, "\n")
  invisible(x)
}

# -----------------------------
# 9. Main-result convenience wrapper
# -----------------------------

run_cw2_main_result <- function(data_path = "prepdata.csv",
                                R = 500,
                                R_schedule = NULL,
                                draw_type = "halton",
                                n_starts = 5,
                                seed = 123,
                                start_seed = 456,
                                compute_hessian = TRUE,
                                hessian_method = "numDeriv",
                                robust_vcov = FALSE,
                                trace = 1,
                                parallel_starts = FALSE) {
  fit <- fit_joint_msl(
    data_path = data_path,
    include_lr_inc = FALSE,
    use_logit_sm = TRUE,
    anchor_indicator = "s_m_input",
    R = R,
    R_schedule = R_schedule,
    draw_type = draw_type,
    seed = seed,
    n_starts = n_starts,
    start_jitter = 0.05,
    start_seed = start_seed,
    nm_maxit = 300,
    bfgs_maxit = 600,
    trace = trace,
    compute_hessian = FALSE,
    parallel_starts = parallel_starts
  )

  if (compute_hessian || robust_vcov) {
    fit <- compute_vcov_joint_msl(
      fit,
      hessian_method = hessian_method,
      robust = robust_vcov
    )
  }

  list(
    fit = fit,
    summary = summarize_joint_msl(fit),
    coefficients_raw = tidy_joint_msl(fit, transformed = FALSE, use_robust = FALSE),
    coefficients_structural = tidy_joint_msl(fit, transformed = TRUE, use_robust = robust_vcov),
    bargaining_index = append_posterior_mu(fit)
  )
}

# -----------------------------
# 10. Example usage
# -----------------------------
# main_res <- run_cw2_main_result(
#   data_path = "prepdata.csv",
#   R = 500,
#   R_schedule = c(100, 300, 500),
#   draw_type = "halton",
#   n_starts = 5,
#   seed = 123,
#   start_seed = 456,
#   compute_hessian = TRUE,
#   hessian_method = "numDeriv",
#   robust_vcov = FALSE,
#   trace = 1,
#   parallel_starts = FALSE
# )
#
# fit <- main_res$fit
# print(fit)
# print(main_res$summary$optimization)
# print(main_res$summary$continuation, n = Inf)
# print(main_res$summary$accounting, n = Inf)
# print(main_res$summary$measurement, n = Inf)
# print(main_res$summary$participation, n = Inf)
# print(main_res$summary$wage, n = Inf)
# print(main_res$coefficients_structural, n = Inf)
#
# saveRDS(main_res, "cw2_main_joint_msl_results.rds")
