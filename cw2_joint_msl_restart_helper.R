# =========================================================
# Coursework 2: Restart helper for the joint-MSL estimator
#
# Purpose:
#   - save the best parameter vector from a previous successful run
#   - restart the model next time from that vector only
#   - optionally skip multi-start and staged continuation for speed
#
# IMPORTANT:
#   This should be used only when the model specification, sample,
#   and variable construction are unchanged. If you change the spec,
#   go back to the full multi-start continuation routine.
# =========================================================

source("cw2_joint_msl_main_msl_upgrade.R")

extract_best_theta <- function(main_res) {
  if (!is.list(main_res) || is.null(main_res$fit) || is.null(main_res$fit$fit$par)) {
    stop("main_res must be the object returned by run_cw2_main_result().")
  }
  main_res$fit$fit$par
}

save_best_theta <- function(main_res,
                            path = "cw2_best_theta.rds",
                            note = NULL) {
  theta <- extract_best_theta(main_res)

  out <- list(
    theta = theta,
    settings = main_res$fit$prepped$settings,
    draw_schedule = main_res$fit$draw_schedule,
    note = note,
    timestamp = Sys.time()
  )

  saveRDS(out, path)
  invisible(path)
}

fit_joint_msl_from_theta <- function(data_path = "prepdata.csv",
                                     theta_start,
                                     include_lr_inc = FALSE,
                                     use_logit_sm = TRUE,
                                     anchor_indicator = "s_m_input",
                                     R = 500,
                                     draw_type = c("halton", "antithetic"),
                                     seed = 123,
                                     method = c("BFGS", "NM_BFGS"),
                                     nm_maxit = 200,
                                     bfgs_maxit = 600,
                                     trace = 1,
                                     compute_hessian = FALSE) {
  draw_type <- match.arg(draw_type)
  method <- match.arg(method)

  prepped <- build_joint_sample(
    data_path = data_path,
    include_lr_inc = include_lr_inc,
    use_logit_sm = use_logit_sm,
    anchor_indicator = anchor_indicator
  )

  dims <- make_dim_list(prepped)
  expected_len <- (dims$J - 1L) + (dims$J - 1L) + dims$J + dims$K_mu + 1L + 1L + dims$K_sel + dims$K_w + 1L + 1L
  if (length(theta_start) != expected_len) {
    stop("theta_start has length ", length(theta_start),
         " but this specification expects length ", expected_len, ".")
  }

  draws <- make_common_draws(R = R, draw_type = draw_type, seed = seed)

  cat("\n============================\n")
  cat("Restart-from-best MSL summary\n")
  cat("============================\n")
  cat("N =", nrow(prepped$joint_df), "\n")
  cat("Workers =", sum(prepped$D == 1L), "\n")
  cat("Indicators =", ncol(prepped$Y), "\n")
  cat("Anchor =", anchor_indicator, "\n")
  cat("R draws =", R, "\n")
  cat("draw_type =", draw_type, "\n")
  cat("method =", method, "\n")

  if (method == "NM_BFGS") {
    nm_fit <- optim(
      par = theta_start,
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
  } else {
    nm_fit <- NULL
    bfgs_fit <- optim(
      par = theta_start,
      fn = objective_joint_msl,
      prepped = prepped,
      draws = draws,
      method = "BFGS",
      hessian = compute_hessian,
      control = list(maxit = bfgs_maxit, trace = trace, REPORT = 25)
    )
  }

  out <- list(
    prepped = prepped,
    init = list(theta = theta_start),
    draws = draws,
    draw_schedule = R,
    stage_results = NULL,
    nm_fit = nm_fit,
    fit = bfgs_fit,
    multi_start_results = tibble(
      start_id = 1L,
      nm_value = if (is.null(nm_fit)) NA_real_ else nm_fit$value,
      nm_conv = if (is.null(nm_fit)) NA_integer_ else nm_fit$convergence,
      bfgs_value = bfgs_fit$value,
      bfgs_conv = bfgs_fit$convergence,
      stage_id = 1L,
      R = R,
      stage_seed = seed
    ),
    all_fits = list(list(nm_fit = nm_fit, fit = bfgs_fit)),
    best_start_id = 1L,
    best_stage_id = 1L,
    hessian = NULL,
    vcov = NULL,
    vcov_robust = NULL,
    hessian_diagnostics = NULL
  )
  class(out) <- "cw2_joint_msl"
  out
}

run_cw2_from_saved_best <- function(theta_path = "cw2_best_theta.rds",
                                    data_path = "prepdata.csv",
                                    R = 500,
                                    draw_type = "halton",
                                    seed = 123,
                                    method = "BFGS",
                                    compute_hessian = TRUE,
                                    hessian_method = "numDeriv",
                                    robust_vcov = FALSE,
                                    trace = 1) {
  saved <- readRDS(theta_path)
  if (is.null(saved$theta)) stop("Saved file does not contain a theta vector.")

  fit <- fit_joint_msl_from_theta(
    data_path = data_path,
    theta_start = saved$theta,
    include_lr_inc = saved$settings$include_lr_inc,
    use_logit_sm = saved$settings$use_logit_sm,
    anchor_indicator = saved$settings$anchor_indicator,
    R = R,
    draw_type = draw_type,
    seed = seed,
    method = method,
    trace = trace,
    compute_hessian = FALSE
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
# Example usage
# -----------------------------
# source("cw2_joint_msl_main_msl_upgrade.R")
# source("cw2_joint_msl_restart_helper.R")
#
# # First full run
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
# # Save the best theta for later reuse
# save_best_theta(main_res, path = "cw2_best_theta.rds")
#
# # Faster rerun next time from that saved theta only
# fast_res <- run_cw2_from_saved_best(
#   theta_path = "cw2_best_theta.rds",
#   data_path = "prepdata.csv",
#   R = 500,
#   draw_type = "halton",
#   seed = 123,
#   method = "BFGS",
#   compute_hessian = TRUE,
#   hessian_method = "numDeriv",
#   robust_vcov = FALSE,
#   trace = 1
# )
