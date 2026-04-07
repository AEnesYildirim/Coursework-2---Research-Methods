# =========================================================
# Standard-error helper for the joint-MSL estimator
#
# Purpose:
#   - diagnose why Hessian-based SEs fail
#   - repair a nearly indefinite Hessian by eigenvalue clipping
#   - optionally build a regularised sandwich / robust VCOV
#
# IMPORTANT:
#   This is an inference helper. It does NOT change point estimates.
#   If the Hessian is repaired, note this clearly in the report.
# =========================================================

source("cw2_joint_msl_main_msl_upgrade.R")

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

compute_joint_msl_se_fix <- function(object,
                                     hessian_method = c("optimHess", "numDeriv"),
                                     use_robust = TRUE,
                                     eig_floor_abs = 1e-6,
                                     eig_floor_rel = 1e-8) {
  if (!inherits(object, "cw2_joint_msl")) {
    stop("object must inherit from 'cw2_joint_msl'.")
  }

  hessian_method <- match.arg(hessian_method)
  H <- compute_hessian_joint_msl(object, hessian_method = hessian_method)

  fix_obj <- make_pd_inverse(
    H,
    eig_floor_abs = eig_floor_abs,
    eig_floor_rel = eig_floor_rel
  )

  out <- object
  out$hessian <- H
  out$hessian_pd <- fix_obj$hessian_pd
  out$vcov <- fix_obj$vcov_pd
  out$hessian_diagnostics <- tibble::tibble(
    min_eigenvalue_original = min(fix_obj$eigen_original, na.rm = TRUE),
    min_eigenvalue_fixed = min(fix_obj$eigen_fixed, na.rm = TRUE),
    max_eigenvalue_original = max(fix_obj$eigen_original, na.rm = TRUE),
    clipping_floor = fix_obj$floor_val,
    n_eigenvalues_clipped = fix_obj$n_clipped,
    has_non_finite = any(!is.finite(H))
  )

  if (use_robust) {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      warning("numDeriv not installed; robust vcov not computed.")
      return(out)
    }

    score_mat <- numDeriv::jacobian(
      func = function(th) joint_loglik_obs(th, prepped = out$prepped, draws = out$draws),
      x = out$fit$par
    )
    meat <- crossprod(score_mat)
    out$score_matrix <- score_mat
    out$vcov_robust <- out$vcov %*% meat %*% out$vcov
  }

  out
}

build_joint_msl_output_from_fit <- function(fit_object,
                                            use_robust = FALSE) {
  list(
    fit = fit_object,
    summary = summarize_joint_msl(fit_object),
    coefficients_raw = tidy_joint_msl(fit_object, transformed = FALSE, use_robust = FALSE),
    coefficients_structural = tidy_joint_msl(fit_object, transformed = TRUE, use_robust = use_robust),
    bargaining_index = append_posterior_mu(fit_object)
  )
}

# -----------------------------
# Example usage
# -----------------------------
# source("cw2_joint_msl_se_fix.R")
#
# # Starting from an already estimated object: main_res
# fit_fixed <- compute_joint_msl_se_fix(
#   main_res$fit,
#   hessian_method = "optimHess",
#   use_robust = TRUE,
#   eig_floor_abs = 1e-6,
#   eig_floor_rel = 1e-8
# )
#
# fixed_res <- build_joint_msl_output_from_fit(fit_fixed, use_robust = TRUE)
#
# print(fixed_res$summary$hessian_diagnostics)
# print(fixed_res$coefficients_structural, n = Inf)
