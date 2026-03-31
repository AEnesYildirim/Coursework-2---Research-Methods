# =========================================================
# Coursework 2: Re-anchor patch for the joint MSL script
# Source this AFTER Coursework2_joint_msl_deltafix.R
# =========================================================

# This patch changes the latent-factor normalization so that the anchor
# indicator can be chosen explicitly. The default is a strong survey item,
# yieldfor.proj.m, rather than s_m_input.

build_joint_sample <- function(data_path = "prepdata.csv",
                               include_lr_inc = FALSE,
                               use_logit_sm = FALSE,
                               anchor_indicator = "yieldfor.proj.m") {

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

  base_y_names <- c("s_m_input", vars$measurement_vars)

  if (!(anchor_indicator %in% base_y_names)) {
    stop(
      "anchor_indicator must be one of: ",
      paste(base_y_names, collapse = ", ")
    )
  }

  # Reorder the measurement system so the chosen anchor is first.
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
      if_all(
        c(base_y_names,
          vars$trait_vars, vars$z_vars_mu,
          "age.f", "educ.f", "D_f"),
        ~ !is.na(.)
      )
    ) %>%
    filter(D_f == 0L | (D_f == 1L & !is.na(W_f)))

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
      use_logit_sm = use_logit_sm,
      anchor_indicator = anchor_indicator
    )
  )
}

fit_joint_msl <- function(data_path = "prepdata.csv",
                          include_lr_inc = FALSE,
                          use_logit_sm = FALSE,
                          anchor_indicator = "yieldfor.proj.m",
                          R = 200,
                          seed = 123,
                          nm_maxit = 300,
                          bfgs_maxit = 600,
                          trace = 1,
                          compute_hessian = FALSE) {

  prepped <- build_joint_sample(
    data_path = data_path,
    include_lr_inc = include_lr_inc,
    use_logit_sm = use_logit_sm,
    anchor_indicator = anchor_indicator
  )

  init <- get_initial_values_joint(prepped)
  draws <- make_antithetic_draws(nrow(prepped$Y), R = R, seed = seed)

  cat("\n============================\n")
  cat("Joint MSL sample summary\n")
  cat("============================\n")
  cat("N =", nrow(prepped$joint_df), "\n")
  cat("Workers =", sum(prepped$D == 1L), "\n")
  cat("Indicators =", ncol(prepped$Y), "\n")
  cat("Anchor =", anchor_indicator, "\n")
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
# Example usage
# -----------------------------
# source("Coursework2_joint_msl_deltafix.R")
# source("Coursework2_joint_msl_reanchor_patch.R")
#
# fit <- fit_joint_msl(
#   data_path = "prepdata.csv",
#   include_lr_inc = FALSE,
#   use_logit_sm = FALSE,
#   anchor_indicator = "yieldfor.proj.m",
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
# structural_table <- tidy_joint_msl(fit, transformed = TRUE)
# bargaining_index_df <- append_posterior_mu(fit)
