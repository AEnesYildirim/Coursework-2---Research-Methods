# =========================================================
# Fast rerun helper for the SAME joint-MSL specification
#
# Purpose:
#   - reuse the previously best theta vector
#   - run only one BFGS refinement at the final draw count
#   - keep the same draw structure and seed for coefficients that
#     are as close as possible to the previous full run
#
# IMPORTANT:
#   Use this only when the specification is unchanged:
#   same data, same sample construction, same anchor, same variables,
#   same draw_type, and same seed.
# =========================================================

source("cw2_joint_msl_main_msl_upgrade.R")
source("cw2_joint_msl_restart_helper.R")

run_cw2_fast_same_spec <- function(theta_path = "cw2_best_theta.rds",
                                   data_path = "prepdata.csv",
                                   R = 500,
                                   draw_type = "halton",
                                   seed = 125,
                                   compute_hessian = FALSE,
                                   hessian_method = "optimHess",
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
    method = "BFGS",
    bfgs_maxit = 400,
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
# Step 1: after one full run, save the best theta
# save_best_theta(main_res, path = "cw2_best_theta.rds")
#
# Step 2: fast rerun next time from the same solution basin
# fast_res <- run_cw2_fast_same_spec(
#   theta_path = "cw2_best_theta.rds",
#   data_path = "prepdata.csv",
#   R = 500,
#   draw_type = "halton",
#   seed = 125,
#   compute_hessian = FALSE,
#   robust_vcov = FALSE,
#   trace = 1
# )
#
# print(fast_res$summary$optimization)
# print(fast_res$summary$selection_scales, n = Inf)
# print(fast_res$summary$participation, n = Inf)
# print(fast_res$summary$wage, n = Inf)
