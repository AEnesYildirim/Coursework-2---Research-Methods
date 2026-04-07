# =========================================================
# Coursework 2: FINAL MAIN ESTIMATION SCRIPT
# One estimation only, Halton draws only, robust SEs repaired
#
# Design:
#   - single final-stage estimation only: R = 500 by default
#   - Halton draws only
#   - prefers restarting from a previously saved best theta
#     to stay as close as possible to the same coefficients
#   - if no saved theta exists, runs one direct estimation from scratch
#   - computes repaired positive-definite Hessian-based VCOV
#   - computes robust sandwich VCOV on top of the repaired Hessian
#
# Notes:
#   - this is intended to be the main estimation code for now
#   - robustness checks can be run later in separate scripts
# =========================================================

source("cw2_joint_msl_main_msl_upgrade.R")
source("cw2_joint_msl_restart_helper.R")
source("cw2_joint_msl_se_fix.R")

run_cw2_main_final_one_run <- function(data_path = "prepdata.csv",
                                       theta_path = "cw2_best_theta.rds",
                                       use_saved_best = TRUE,
                                       R = 500,
                                       draw_type = "halton",
                                       seed = 125,
                                       trace = 1,
                                       hessian_method = "optimHess",
                                       eig_floor_abs = 1e-6,
                                       eig_floor_rel = 1e-8) {
  if (draw_type != "halton") {
    stop("This final script is designed to use Halton draws only.")
  }

  has_saved_theta <- isTRUE(use_saved_best) && file.exists(theta_path)

  if (has_saved_theta) {
    message("Using saved best theta for fastest same-spec rerun.")

    saved <- readRDS(theta_path)
    if (is.null(saved$theta)) stop("Saved theta file does not contain a theta vector.")

    fit0 <- fit_joint_msl_from_theta(
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
  } else {
    message("No saved theta found. Running one direct estimation from scratch.")

    fit0 <- fit_joint_msl(
      data_path = data_path,
      include_lr_inc = FALSE,
      use_logit_sm = TRUE,
      anchor_indicator = "s_m_input",
      R = R,
      R_schedule = c(R),
      draw_type = draw_type,
      seed = seed,
      n_starts = 1,
      start_jitter = 0.05,
      start_seed = 456,
      nm_maxit = 300,
      bfgs_maxit = 600,
      trace = trace,
      compute_hessian = FALSE,
      parallel_starts = FALSE
    )
  }

  fit_final <- compute_joint_msl_se_fix(
    fit0,
    hessian_method = hessian_method,
    use_robust = TRUE,
    eig_floor_abs = eig_floor_abs,
    eig_floor_rel = eig_floor_rel
  )

  out <- build_joint_msl_output_from_fit(
    fit_final,
    use_robust = TRUE
  )

  out$run_settings <- list(
    used_saved_best = has_saved_theta,
    theta_path = theta_path,
    R = R,
    draw_type = draw_type,
    seed = seed,
    hessian_method = hessian_method,
    eig_floor_abs = eig_floor_abs,
    eig_floor_rel = eig_floor_rel
  )

  out
}

# ---------------------------------------------------------
# RECOMMENDED MAIN RUN
# ---------------------------------------------------------
# If you already have cw2_best_theta.rds from a previous full run,
# this will be the fastest computation that stays closest to the
# same coefficients while still re-estimating properly.
#
# main_final <- run_cw2_main_final_one_run(
#   data_path = "prepdata.csv",
#   theta_path = "cw2_best_theta.rds",
#   use_saved_best = TRUE,
#   R = 500,
#   draw_type = "halton",
#   seed = 125,
#   trace = 1,
#   hessian_method = "optimHess",
#   eig_floor_abs = 1e-6,
#   eig_floor_rel = 1e-8
# )
#
# print(main_final$summary$optimization)
# print(main_final$summary$hessian_diagnostics)
# print(main_final$summary$selection_scales, n = Inf)
# print(main_final$summary$participation, n = Inf)
# print(main_final$summary$wage, n = Inf)
# print(main_final$coefficients_structural, n = Inf)
# print(cor(main_final$bargaining_index$mu_post,
#           main_final$bargaining_index$s_m,
#           use = "complete.obs"))
#
# saveRDS(main_final, "cw2_main_final_one_run_results.rds")
# write.csv(main_final$coefficients_structural,
#           "cw2_main_final_one_run_coefficients.csv",
#           row.names = FALSE)
# write.csv(main_final$bargaining_index,
#           "cw2_main_final_one_run_bargaining_index.csv",
#           row.names = FALSE)
