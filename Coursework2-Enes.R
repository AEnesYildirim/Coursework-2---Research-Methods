# =========================================================
# COURSEWORK SCRIPT: Household bargaining, employment, wages
# =========================================================

# -----------------------------
# 0. Packages
# -----------------------------
library(tidyverse)
library(broom)
library(sampleSelection)

# =========================================================
# STEP 1: Load data, define variables, and inspect missingness
# =========================================================

# 1. Read the data
df <- read_csv("prepdata.csv", show_col_types = FALSE)

# 2. Make sure hid exists BEFORE anything else
if (!("hid" %in% names(df))) {
  df <- df %>% mutate(hid = row_number())
}

# 3. Define variable groups from the assignment
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

z_vars <- c(
  "duration", "abouttobeparents", "lr_inc",
  "d_educ", "d_age", "both_m"
)

# corrected participation controls: exclude lr_inc
z_vars_participation <- c(
  "duration", "abouttobeparents",
  "d_educ", "d_age", "both_m"
)

female_outcome_vars <- c("age.f", "educ.f", "income.f", "hours.f")

# 4. Check required columns
required_vars <- c("hid", "s_m", measurement_vars, trait_vars, z_vars, female_outcome_vars)

missing_columns <- setdiff(required_vars, names(df))

if (length(missing_columns) > 0) {
  stop("These required columns are missing from the data: ",
       paste(missing_columns, collapse = ", "))
}

# 5. Create baseline labour variables
df <- df %>%
  mutate(
    both_m = as.integer(both_m),
    D_f = case_when(
      is.na(income.f) ~ NA_integer_,
      income.f > 0    ~ 1L,
      TRUE            ~ 0L
    ),
    W_f = case_when(
      !is.na(income.f) & income.f > 0 ~ log(income.f),
      TRUE                            ~ NA_real_
    )
  )

# 6. Build analysis frame
analysis_df <- df %>%
  select(
    hid,
    s_m,
    all_of(measurement_vars),
    all_of(trait_vars),
    all_of(z_vars),
    age.f, educ.f, income.f, hours.f,
    D_f, W_f
  )

# 7. Missingness table
missingness_tbl <- analysis_df %>%
  summarise(across(everything(), ~ mean(is.na(.)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "missing_share"
  ) %>%
  arrange(desc(missing_share))

cat("\n=== Missingness table ===\n")
print(missingness_tbl, n = nrow(missingness_tbl))

# 8. Sample sizes
sample_sizes <- tibble(
  stage = c(
    "Raw data",
    "Measurement model sample: s_m + 12 measurement vars",
    "Bargaining determinants sample: s_m + measurement vars + trait vars + Z",
    "Participation model sample (corrected): D_f + mu_hat + controls without lr_inc",
    "Wage model sample: W_f + age.f + educ.f"
  ),
  n = c(
    nrow(df),
    analysis_df %>%
      filter(if_all(c("s_m", all_of(measurement_vars)), ~ !is.na(.))) %>%
      nrow(),
    analysis_df %>%
      filter(if_all(c("s_m", all_of(measurement_vars), all_of(trait_vars), all_of(z_vars)), ~ !is.na(.))) %>%
      nrow(),
    analysis_df %>%
      filter(if_all(c("D_f", "age.f", "educ.f", all_of(z_vars_participation)), ~ !is.na(.))) %>%
      nrow(),
    analysis_df %>%
      filter(!is.na(W_f), !is.na(age.f), !is.na(educ.f)) %>%
      nrow()
  )
)

cat("\n=== Sample sizes ===\n")
print(sample_sizes)

# 9. Basic descriptives
desc_tbl <- analysis_df %>%
  summarise(
    n = n(),
    s_m_mean = mean(s_m, na.rm = TRUE),
    s_m_sd   = sd(s_m, na.rm = TRUE),
    female_employment_rate = mean(D_f, na.rm = TRUE),
    female_income_mean = mean(income.f, na.rm = TRUE),
    female_income_median = median(income.f, na.rm = TRUE),
    female_wage_mean = mean(W_f, na.rm = TRUE),
    female_age_mean = mean(age.f, na.rm = TRUE),
    female_educ_mean = mean(educ.f, na.rm = TRUE)
  )

cat("\n=== Descriptives ===\n")
print(desc_tbl)

# 10. Quick plots
ggplot(analysis_df, aes(x = s_m)) +
  geom_histogram(bins = 30) +
  labs(title = "Distribution of male private consumption share (s_m)")

ggplot(analysis_df, aes(x = factor(D_f))) +
  geom_bar() +
  labs(title = "Female employment indicator", x = "D_f", y = "Count")

# =========================================================
# STEP 1B: Sanity checks before factor analysis
# =========================================================

range_tbl <- analysis_df %>%
  select(s_m, all_of(measurement_vars)) %>%
  summarise(across(
    everything(),
    list(
      min = ~ min(., na.rm = TRUE),
      p25 = ~ quantile(., 0.25, na.rm = TRUE),
      median = ~ median(., na.rm = TRUE),
      mean = ~ mean(., na.rm = TRUE),
      p75 = ~ quantile(., 0.75, na.rm = TRUE),
      max = ~ max(., na.rm = TRUE)
    )
  )) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_sep = "_(?=[^_]+$)"
  )

cat("\n=== Range table for measurement variables ===\n")
print(range_tbl, n = nrow(range_tbl))

corr_mat <- analysis_df %>%
  select(s_m, all_of(measurement_vars)) %>%
  cor(use = "pairwise.complete.obs")

cat("\n=== Correlation matrix: s_m and 12 measurement variables ===\n")
print(round(corr_mat, 2))

# =========================================================
# STEP 2: Measurement model for bargaining
# =========================================================

# Remove old versions if rerunning
df <- df %>% select(-any_of(c("mu_hat_13", "mu_hat_12")))
analysis_df <- analysis_df %>% select(-any_of(c("mu_hat_13", "mu_hat_12")))

measurement_13 <- c("s_m", measurement_vars)
measurement_12 <- measurement_vars

# 2.1 Complete-case samples
fa_sample_13 <- df %>%
  select(hid, all_of(measurement_13)) %>%
  drop_na()

fa_sample_12 <- df %>%
  select(hid, all_of(measurement_12)) %>%
  drop_na()

cat("\n13-indicator sample size:", nrow(fa_sample_13), "\n")
cat("12-indicator sample size:", nrow(fa_sample_12), "\n")

# 2.2 Standardise
X13 <- fa_sample_13 %>%
  select(-hid) %>%
  mutate(across(everything(), ~ as.numeric(scale(.x))))

X12 <- fa_sample_12 %>%
  select(-hid) %>%
  mutate(across(everything(), ~ as.numeric(scale(.x))))

# 2.3 Factor models
fa13 <- factanal(
  x = as.matrix(X13),
  factors = 1,
  scores = "regression"
)

fa12 <- factanal(
  x = as.matrix(X12),
  factors = 1,
  scores = "regression"
)

cat("\n=== 13-indicator factor model ===\n")
print(fa13, digits = 3, cutoff = 0)

cat("\n=== 12-indicator factor model ===\n")
print(fa12, digits = 3, cutoff = 0)

# 2.4 Loading tables
load13_raw <- tibble(
  variable = rownames(fa13$loadings),
  loading = as.numeric(fa13$loadings[, 1]),
  uniqueness = as.numeric(fa13$uniquenesses)
)

load12_raw <- tibble(
  variable = rownames(fa12$loadings),
  loading = as.numeric(fa12$loadings[, 1]),
  uniqueness = as.numeric(fa12$uniquenesses)
)

# 2.5 Sign normalisation using strongest-loading variable
ref13 <- load13_raw %>% slice_max(order_by = abs(loading), n = 1)
flip13 <- ifelse(ref13$loading < 0, -1, 1)

load13 <- load13_raw %>%
  mutate(loading = loading * flip13) %>%
  arrange(desc(abs(loading)))

scores13 <- tibble(
  hid = fa_sample_13$hid,
  mu_hat_13 = as.numeric(fa13$scores[, 1]) * flip13
)

ref12 <- load12_raw %>% slice_max(order_by = abs(loading), n = 1)
flip12 <- ifelse(ref12$loading < 0, -1, 1)

load12 <- load12_raw %>%
  mutate(loading = loading * flip12) %>%
  arrange(desc(abs(loading)))

scores12 <- tibble(
  hid = fa_sample_12$hid,
  mu_hat_12 = as.numeric(fa12$scores[, 1]) * flip12
)

# 2.6 Merge scores back
df <- df %>%
  left_join(scores13, by = "hid") %>%
  left_join(scores12, by = "hid")

analysis_df <- analysis_df %>%
  left_join(scores13, by = "hid") %>%
  left_join(scores12, by = "hid")

# 2.7 Diagnostics
diag13 <- tibble(
  model = "13 indicators: s_m + 12 projections",
  n = nrow(fa_sample_13),
  strongest_loading_var = ref13$variable,
  strongest_loading = abs(ref13$loading),
  sm_loading = load13_raw %>% filter(variable == "s_m") %>% pull(loading),
  sm_uniqueness = load13_raw %>% filter(variable == "s_m") %>% pull(uniqueness),
  corr_mu13_sm = cor(scores13$mu_hat_13, fa_sample_13$s_m, use = "complete.obs"),
  mu13_mean = mean(scores13$mu_hat_13),
  mu13_sd = sd(scores13$mu_hat_13)
)

diag12 <- tibble(
  model = "12 indicators: survey projections only",
  n = nrow(fa_sample_12),
  strongest_loading_var = ref12$variable,
  strongest_loading = abs(ref12$loading),
  sm_loading = NA_real_,
  sm_uniqueness = NA_real_,
  corr_mu13_sm = NA_real_,
  mu13_mean = mean(scores12$mu_hat_12),
  mu13_sd = sd(scores12$mu_hat_12)
)

score_compare <- scores13 %>%
  inner_join(scores12, by = "hid") %>%
  summarise(
    n_overlap = n(),
    corr_mu13_mu12 = cor(mu_hat_13, mu_hat_12, use = "complete.obs")
  )

cat("\n=== Diagnostics: 13-indicator model ===\n")
print(diag13)

cat("\n=== Diagnostics: 12-indicator model ===\n")
print(diag12)

cat("\n=== Comparison of factor scores ===\n")
print(score_compare)

cat("\n=== Loadings: 13-indicator model ===\n")
print(load13, n = nrow(load13))

cat("\n=== Loadings: 12-indicator model ===\n")
print(load12, n = nrow(load12))

sm_flag <- diag13 %>%
  mutate(
    sm_weak = abs(sm_loading) < 0.10 | sm_uniqueness > 0.95
  )

cat("\n=== Is s_m weak in the common factor? ===\n")
print(sm_flag)

# 2.8 Plots
ggplot(scores13, aes(x = mu_hat_13)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Latent bargaining factor: 13-indicator model",
    x = "mu_hat_13",
    y = "Count"
  )

ggplot(load13, aes(x = reorder(variable, loading), y = loading)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Factor loadings: 13-indicator model",
    x = "",
    y = "Loading"
  )

ggplot(scores12, aes(x = mu_hat_12)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Latent bargaining factor: 12-indicator survey-only model",
    x = "mu_hat_12",
    y = "Count"
  )

ggplot(load12, aes(x = reorder(variable, loading), y = loading)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Factor loadings: 12-indicator survey-only model",
    x = "",
    y = "Loading"
  )

# =========================================================
# STEP 3: What drives bargaining?
# =========================================================

determinants_vars <- c(trait_vars, z_vars)

# 3.1 Main model: mu_hat_13
mu_reg_sample_13 <- df %>%
  select(hid, mu_hat_13, all_of(determinants_vars)) %>%
  drop_na()

cat("\nDeterminants sample size (mu_hat_13):", nrow(mu_reg_sample_13), "\n")

mu_formula_13 <- as.formula(
  paste("mu_hat_13 ~", paste(determinants_vars, collapse = " + "))
)

mu_lm_13 <- lm(mu_formula_13, data = mu_reg_sample_13)

cat("\n=== Main regression: determinants of mu_hat_13 ===\n")
print(summary(mu_lm_13))

coef_tbl_13 <- broom::tidy(mu_lm_13, conf.int = TRUE) %>%
  arrange(p.value)

cat("\n=== Coefficients: mu_hat_13 ===\n")
print(coef_tbl_13, n = nrow(coef_tbl_13))

fit_tbl_13 <- broom::glance(mu_lm_13) %>%
  select(r.squared, adj.r.squared, sigma, statistic, p.value, df, df.residual)

cat("\n=== Fit statistics: mu_hat_13 ===\n")
print(fit_tbl_13)

# 3.2 Robustness model: mu_hat_12
mu_reg_sample_12 <- df %>%
  select(hid, mu_hat_12, all_of(determinants_vars)) %>%
  drop_na()

cat("\nDeterminants sample size (mu_hat_12):", nrow(mu_reg_sample_12), "\n")

mu_formula_12 <- as.formula(
  paste("mu_hat_12 ~", paste(determinants_vars, collapse = " + "))
)

mu_lm_12 <- lm(mu_formula_12, data = mu_reg_sample_12)

cat("\n=== Robustness regression: determinants of mu_hat_12 ===\n")
print(summary(mu_lm_12))

coef_tbl_12 <- broom::tidy(mu_lm_12, conf.int = TRUE) %>%
  arrange(p.value)

cat("\n=== Coefficients: mu_hat_12 ===\n")
print(coef_tbl_12, n = nrow(coef_tbl_12))

fit_tbl_12 <- broom::glance(mu_lm_12) %>%
  select(r.squared, adj.r.squared, sigma, statistic, p.value, df, df.residual)

cat("\n=== Fit statistics: mu_hat_12 ===\n")
print(fit_tbl_12)

# 3.3 Compare the two versions
coef_compare <- coef_tbl_13 %>%
  select(term, estimate_13 = estimate, p_13 = p.value) %>%
  inner_join(
    coef_tbl_12 %>%
      select(term, estimate_12 = estimate, p_12 = p.value),
    by = "term"
  ) %>%
  mutate(
    same_sign = sign(estimate_13) == sign(estimate_12),
    abs_diff = abs(estimate_13 - estimate_12)
  ) %>%
  arrange(desc(abs_diff))

cat("\n=== Coefficient comparison: mu_hat_13 vs mu_hat_12 ===\n")
print(coef_compare, n = nrow(coef_compare))

top_coef_13 <- coef_tbl_13 %>%
  filter(term != "(Intercept)") %>%
  mutate(abs_estimate = abs(estimate)) %>%
  arrange(desc(abs_estimate)) %>%
  slice_head(n = 15)

cat("\n=== Largest coefficients in main model ===\n")
print(top_coef_13)

ggplot(top_coef_13, aes(x = reorder(term, estimate), y = estimate)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Largest coefficients: determinants of mu_hat_13",
    x = "",
    y = "Coefficient"
  )

# =========================================================
# STEP 4: Female employment model (corrected)
# =========================================================

# 4.1 Main participation model
participation_vars_13 <- c(
  "mu_hat_13",
  "age.f", "educ.f",
  "duration", "abouttobeparents",
  "d_educ", "d_age", "both_m"
)

part_sample_13 <- df %>%
  select(hid, D_f, all_of(participation_vars_13)) %>%
  drop_na()

cat("\nParticipation sample size (mu_hat_13):", nrow(part_sample_13), "\n")

part_formula_13 <- as.formula(
  paste("D_f ~", paste(participation_vars_13, collapse = " + "))
)

part_probit_13 <- glm(
  formula = part_formula_13,
  data = part_sample_13,
  family = binomial(link = "probit")
)

cat("\n=== Probit: female employment on mu_hat_13 ===\n")
print(summary(part_probit_13))

part_coef_tbl_13 <- broom::tidy(part_probit_13, conf.int = TRUE) %>%
  arrange(p.value)

cat("\n=== Coefficients: participation model with mu_hat_13 ===\n")
print(part_coef_tbl_13, n = nrow(part_coef_tbl_13))

part_fit_tbl_13 <- broom::glance(part_probit_13) %>%
  select(null.deviance, df.null, deviance, df.residual, AIC)

cat("\n=== Fit statistics: participation model with mu_hat_13 ===\n")
print(part_fit_tbl_13)

mean_pred_prob_13 <- mean(predict(part_probit_13, type = "response"))
cat("\nAverage predicted female employment probability (mu_hat_13):",
    mean_pred_prob_13, "\n")

# 4.2 Robustness participation model
participation_vars_12 <- c(
  "mu_hat_12",
  "age.f", "educ.f",
  "duration", "abouttobeparents",
  "d_educ", "d_age", "both_m"
)

part_sample_12 <- df %>%
  select(hid, D_f, all_of(participation_vars_12)) %>%
  drop_na()

cat("\nParticipation sample size (mu_hat_12):", nrow(part_sample_12), "\n")

part_formula_12 <- as.formula(
  paste("D_f ~", paste(participation_vars_12, collapse = " + "))
)

part_probit_12 <- glm(
  formula = part_formula_12,
  data = part_sample_12,
  family = binomial(link = "probit")
)

cat("\n=== Probit: female employment on mu_hat_12 ===\n")
print(summary(part_probit_12))

part_coef_tbl_12 <- broom::tidy(part_probit_12, conf.int = TRUE) %>%
  arrange(p.value)

cat("\n=== Coefficients: participation model with mu_hat_12 ===\n")
print(part_coef_tbl_12, n = nrow(part_coef_tbl_12))

part_fit_tbl_12 <- broom::glance(part_probit_12) %>%
  select(null.deviance, df.null, deviance, df.residual, AIC)

cat("\n=== Fit statistics: participation model with mu_hat_12 ===\n")
print(part_fit_tbl_12)

mean_pred_prob_12 <- mean(predict(part_probit_12, type = "response"))
cat("\nAverage predicted female employment probability (mu_hat_12):",
    mean_pred_prob_12, "\n")

# 4.3 Compare the two versions
part_compare <- part_coef_tbl_13 %>%
  select(term, estimate_13 = estimate, p_13 = p.value) %>%
  inner_join(
    part_coef_tbl_12 %>%
      select(term, estimate_12 = estimate, p_12 = p.value),
    by = "term"
  ) %>%
  mutate(
    same_sign = sign(estimate_13) == sign(estimate_12),
    abs_diff = abs(estimate_13 - estimate_12)
  ) %>%
  arrange(desc(abs_diff))

cat("\n=== Comparison: participation model with mu_hat_13 vs mu_hat_12 ===\n")
print(part_compare, n = nrow(part_compare))

# 4.4 Approximate marginal effects at sample means
xb_bar_13 <- mean(predict(part_probit_13, type = "link"))
me_mu_13 <- dnorm(xb_bar_13) * coef(part_probit_13)["mu_hat_13"]

cat("\nApproximate marginal effect of mu_hat_13 at sample means:",
    me_mu_13, "\n")

xb_bar_12 <- mean(predict(part_probit_12, type = "link"))
me_mu_12 <- dnorm(xb_bar_12) * coef(part_probit_12)["mu_hat_12"]

cat("Approximate marginal effect of mu_hat_12 at sample means:",
    me_mu_12, "\n")

plot_df_13 <- part_sample_13 %>%
  mutate(pred_prob = predict(part_probit_13, type = "response"))

ggplot(plot_df_13, aes(x = mu_hat_13, y = pred_prob)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE) +
  labs(
    title = "Predicted female employment probability vs bargaining factor",
    x = "mu_hat_13",
    y = "Predicted probability of employment"
  )

# =========================================================
# STEP 5: Female wages and selection into work
# =========================================================

# 5.1 Baseline wage equation (workers only)
wage_vars <- c("W_f", "age.f", "educ.f")

wage_sample <- df %>%
  select(hid, all_of(wage_vars)) %>%
  drop_na()

cat("\nWage sample size:", nrow(wage_sample), "\n")

wage_lm <- lm(W_f ~ age.f + educ.f, data = wage_sample)

cat("\n=== OLS wage equation (workers only) ===\n")
print(summary(wage_lm))

wage_coef_tbl <- broom::tidy(wage_lm, conf.int = TRUE) %>%
  arrange(p.value)

cat("\n=== Wage equation coefficients ===\n")
print(wage_coef_tbl, n = nrow(wage_coef_tbl))

wage_fit_tbl <- broom::glance(wage_lm) %>%
  select(r.squared, adj.r.squared, sigma, statistic, p.value, df, df.residual)

cat("\n=== Wage equation fit statistics ===\n")
print(wage_fit_tbl)

# 5.2 Heckman-style selection model
heckman_vars <- c(
  "D_f", "W_f",
  "mu_hat_13",
  "age.f", "educ.f",
  "duration", "abouttobeparents",
  "d_educ", "d_age", "both_m"
)

heckman_sample <- df %>%
  select(hid, all_of(heckman_vars)) %>%
  mutate(D_f = as.integer(D_f)) %>%
  drop_na(D_f, age.f, educ.f, mu_hat_13, duration, abouttobeparents, d_educ, d_age, both_m)

cat("\nHeckman sample size (selection vars observed):", nrow(heckman_sample), "\n")
cat("Observed wages among these:", sum(!is.na(heckman_sample$W_f)), "\n")

heckman_fit <- selection(
  selection = D_f ~ mu_hat_13 + age.f + educ.f + duration +
    abouttobeparents + d_educ + d_age + both_m,
  outcome = W_f ~ age.f + educ.f,
  data = heckman_sample,
  method = "ml"
)

cat("\n=== Heckman selection model (ML) ===\n")
print(summary(heckman_fit))

# 5.3 Correct extraction of coefficients by equation
est_mat <- as.data.frame(summary(heckman_fit)$estimate)
est_mat$term <- rownames(est_mat)
rownames(est_mat) <- NULL

# Fixed block lengths for THIS model specification:
# selection: intercept + 8 regressors = 9
# outcome: intercept + 2 regressors = 3
# error terms: sigma, rho = 2
n_sel <- 9
n_out <- 3

est_mat$equation <- c(
  rep("selection", n_sel),
  rep("outcome", n_out),
  rep("error", nrow(est_mat) - n_sel - n_out)
)

# Reorder columns if available
if (all(c("Estimate", "Std. Error", "t value", "Pr(>|t|)") %in% names(est_mat))) {
  est_mat <- est_mat %>%
    select(equation, term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`)
}

cat("\n=== Heckman coefficients with equation labels ===\n")
print(est_mat, row.names = FALSE)

# 5.4 Pull out selection parameters
rho_row <- est_mat %>% filter(term == "rho")
sigma_row <- est_mat %>% filter(term == "sigma")

cat("\n=== Selection-on-unobservables parameters ===\n")
print(rho_row, row.names = FALSE)
print(sigma_row, row.names = FALSE)

# 5.5 Correct OLS vs Heckman wage-equation comparison
heckman_outcome <- est_mat %>%
  filter(equation == "outcome") %>%
  select(term, heckman_estimate = Estimate)

coef_compare_wage <- wage_coef_tbl %>%
  select(term, ols_estimate = estimate, ols_p = p.value) %>%
  left_join(heckman_outcome, by = "term") %>%
  mutate(abs_diff = abs(ols_estimate - heckman_estimate))

cat("\n=== OLS vs Heckman wage-coefficient comparison ===\n")
print(coef_compare_wage, n = nrow(coef_compare_wage))

# 5.6 Optional fitted-vs-observed wage plot
wage_sample_pred <- wage_sample %>%
  mutate(pred_wage_ols = predict(wage_lm, newdata = wage_sample))

ggplot(wage_sample_pred, aes(x = pred_wage_ols, y = W_f)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE) +
  labs(
    title = "Observed vs fitted log female wages (OLS, workers only)",
    x = "Fitted log wage",
    y = "Observed log wage"
  )



# =========================================================
# STEP 6: Robustness checks for labour outcomes
# Use mu_hat_12 instead of mu_hat_13
# =========================================================

library(tidyverse)
library(broom)
library(sampleSelection)

# -----------------------------
# 6.1 Robustness participation model with mu_hat_12
# -----------------------------
part_sample_robust <- df %>%
  select(
    hid, D_f, mu_hat_12,
    age.f, educ.f,
    duration, abouttobeparents,
    d_educ, d_age, both_m
  ) %>%
  drop_na()

cat("\nRobustness participation sample size (mu_hat_12):", nrow(part_sample_robust), "\n")

part_probit_robust <- glm(
  D_f ~ mu_hat_12 + age.f + educ.f + duration +
    abouttobeparents + d_educ + d_age + both_m,
  data = part_sample_robust,
  family = binomial(link = "probit")
)

cat("\n=== Robustness probit: female employment on mu_hat_12 ===\n")
print(summary(part_probit_robust))

part_coef_tbl_robust <- broom::tidy(part_probit_robust, conf.int = TRUE) %>%
  arrange(p.value)

cat("\n=== Robustness participation coefficients ===\n")
print(part_coef_tbl_robust, n = nrow(part_coef_tbl_robust))

part_fit_tbl_robust <- broom::glance(part_probit_robust) %>%
  select(null.deviance, df.null, deviance, df.residual, AIC)

cat("\n=== Robustness participation fit statistics ===\n")
print(part_fit_tbl_robust)

# Compare main vs robustness participation coefficients
part_compare_main_robust <- part_coef_tbl_13 %>%
  select(term, main_estimate = estimate, main_p = p.value) %>%
  inner_join(
    part_coef_tbl_robust %>%
      select(term, robust_estimate = estimate, robust_p = p.value),
    by = "term"
  ) %>%
  mutate(
    same_sign = sign(main_estimate) == sign(robust_estimate),
    abs_diff = abs(main_estimate - robust_estimate)
  ) %>%
  arrange(desc(abs_diff))

cat("\n=== Main vs robustness participation comparison ===\n")
print(part_compare_main_robust, n = nrow(part_compare_main_robust))

# -----------------------------
# 6.2 Robustness Heckman selection model with mu_hat_12
# -----------------------------
heckman_sample_robust <- df %>%
  select(
    hid, D_f, W_f,
    mu_hat_12,
    age.f, educ.f,
    duration, abouttobeparents,
    d_educ, d_age, both_m
  ) %>%
  mutate(D_f = as.integer(D_f)) %>%
  drop_na(D_f, age.f, educ.f, mu_hat_12, duration, abouttobeparents, d_educ, d_age, both_m)

cat("\nRobustness Heckman sample size:", nrow(heckman_sample_robust), "\n")
cat("Observed wages in robustness Heckman sample:", sum(!is.na(heckman_sample_robust$W_f)), "\n")

heckman_fit_robust <- selection(
  selection = D_f ~ mu_hat_12 + age.f + educ.f + duration +
    abouttobeparents + d_educ + d_age + both_m,
  outcome = W_f ~ age.f + educ.f,
  data = heckman_sample_robust,
  method = "ml"
)

cat("\n=== Robustness Heckman selection model (mu_hat_12) ===\n")
print(summary(heckman_fit_robust))

# -----------------------------
# 6.3 Extract robustness Heckman coefficients by equation
# -----------------------------
est_mat_robust <- as.data.frame(summary(heckman_fit_robust)$estimate)
est_mat_robust$term <- rownames(est_mat_robust)
rownames(est_mat_robust) <- NULL

# Same block lengths as main model:
# selection: intercept + 8 regressors = 9
# outcome: intercept + 2 regressors = 3
n_sel <- 9
n_out <- 3

est_mat_robust$equation <- c(
  rep("selection", n_sel),
  rep("outcome", n_out),
  rep("error", nrow(est_mat_robust) - n_sel - n_out)
)

if (all(c("Estimate", "Std. Error", "t value", "Pr(>|t|)") %in% names(est_mat_robust))) {
  est_mat_robust <- est_mat_robust %>%
    select(equation, term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`)
}

cat("\n=== Robustness Heckman coefficients with equation labels ===\n")
print(est_mat_robust)

rho_row_robust <- est_mat_robust %>% filter(term == "rho")
sigma_row_robust <- est_mat_robust %>% filter(term == "sigma")

cat("\n=== Robustness selection-on-unobservables parameters ===\n")
print(rho_row_robust)
print(sigma_row_robust)

# -----------------------------
# 6.4 Compare main vs robustness Heckman rho
# -----------------------------
rho_compare <- bind_rows(
  rho_row %>%
    mutate(model = "main_mu_hat_13") %>%
    select(model, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`),
  rho_row_robust %>%
    mutate(model = "robust_mu_hat_12") %>%
    select(model, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`)
)

cat("\n=== Main vs robustness rho comparison ===\n")
print(rho_compare)

# -----------------------------
# 6.5 Compare main vs robustness wage-equation coefficients
# -----------------------------
heckman_outcome_main <- est_mat %>%
  filter(equation == "outcome") %>%
  select(term, main_heckman = Estimate)

heckman_outcome_robust <- est_mat_robust %>%
  filter(equation == "outcome") %>%
  select(term, robust_heckman = Estimate)

wage_compare_main_robust <- heckman_outcome_main %>%
  inner_join(heckman_outcome_robust, by = "term") %>%
  mutate(
    abs_diff = abs(main_heckman - robust_heckman)
  )

cat("\n=== Main vs robustness Heckman wage-equation comparison ===\n")
print(wage_compare_main_robust)

# -----------------------------
# 6.6 Compact interpretation helpers
# -----------------------------
robust_summary_tbl <- tibble(
  model = c("Participation main (mu_hat_13)",
            "Participation robust (mu_hat_12)",
            "Heckman rho main",
            "Heckman rho robust"),
  estimate = c(
    part_coef_tbl_13 %>% filter(term == "mu_hat_13") %>% pull(estimate),
    part_coef_tbl_robust %>% filter(term == "mu_hat_12") %>% pull(estimate),
    rho_row %>% pull(Estimate),
    rho_row_robust %>% pull(Estimate)
  ),
  p_value = c(
    part_coef_tbl_13 %>% filter(term == "mu_hat_13") %>% pull(p.value),
    part_coef_tbl_robust %>% filter(term == "mu_hat_12") %>% pull(p.value),
    rho_row %>% pull(`Pr(>|t|)`),
    rho_row_robust %>% pull(`Pr(>|t|)`)
  )
)

cat("\n=== Compact robustness summary ===\n")
print(robust_summary_tbl)


# =========================================================
# STEP 7: Final tables for the report
# =========================================================

library(tidyverse)

# -----------------------------
# TABLE 1: Descriptives and sample sizes
# -----------------------------
table1_desc <- tibble(
  statistic = c(
    "N raw couples",
    "Measurement sample (13 indicators)",
    "Determinants sample",
    "Participation sample",
    "Wage sample",
    "Mean s_m",
    "SD s_m",
    "Female employment rate",
    "Mean female income",
    "Median female income",
    "Mean log female income"
  ),
  value = c(
    nrow(df),
    nrow(fa_sample_13),
    nrow(mu_reg_sample_13),
    nrow(part_sample_13),
    nrow(wage_sample),
    mean(df$s_m, na.rm = TRUE),
    sd(df$s_m, na.rm = TRUE),
    mean(df$D_f, na.rm = TRUE),
    mean(df$income.f, na.rm = TRUE),
    median(df$income.f, na.rm = TRUE),
    mean(df$W_f, na.rm = TRUE)
  )
)

cat("\n=== TABLE 1: Descriptives and sample sizes ===\n")
print(table1_desc, n = nrow(table1_desc))

# -----------------------------
# TABLE 2: Measurement model
# -----------------------------
table2_measurement <- load13 %>%
  mutate(
    abs_loading = abs(loading),
    weak_indicator = uniqueness > 0.95
  ) %>%
  select(variable, loading, uniqueness, weak_indicator)

cat("\n=== TABLE 2: Measurement model loadings (13-indicator factor) ===\n")
print(table2_measurement, n = nrow(table2_measurement))

table2_diag <- diag13 %>%
  select(
    n,
    strongest_loading_var,
    strongest_loading,
    sm_loading,
    sm_uniqueness,
    corr_mu13_sm
  )

cat("\n=== TABLE 2A: Key diagnostics for measurement model ===\n")
print(table2_diag)

# -----------------------------
# TABLE 3: Determinants + participation
# -----------------------------
table3_determinants <- coef_tbl_13 %>%
  filter(term != "(Intercept)") %>%
  mutate(
    sig_10 = p.value < 0.10,
    sig_05 = p.value < 0.05
  ) %>%
  arrange(p.value)

cat("\n=== TABLE 3A: Determinants of bargaining (mu_hat_13) ===\n")
print(table3_determinants, n = 15)

table3_participation <- part_coef_tbl_13 %>%
  mutate(
    sig_10 = p.value < 0.10,
    sig_05 = p.value < 0.05
  ) %>%
  select(term, estimate, std.error, p.value, conf.low, conf.high, sig_10, sig_05)

cat("\n=== TABLE 3B: Female employment probit ===\n")
print(table3_participation, n = nrow(table3_participation))

# -----------------------------
# TABLE 4: Wages and selection
# -----------------------------
table4_wage_ols <- wage_coef_tbl %>%
  select(term, ols_estimate = estimate, ols_se = std.error, ols_p = p.value)

table4_wage_heckman <- est_mat %>%
  filter(equation == "outcome") %>%
  select(term, heckman_estimate = Estimate, heckman_se = `Std. Error`, heckman_p = `Pr(>|t|)`)

table4_wage_compare <- table4_wage_ols %>%
  full_join(table4_wage_heckman, by = "term")

cat("\n=== TABLE 4A: Wage equation comparison (OLS vs Heckman outcome) ===\n")
print(table4_wage_compare, n = nrow(table4_wage_compare))

table4_selection <- bind_rows(
  rho_row %>% mutate(parameter = "rho"),
  sigma_row %>% mutate(parameter = "sigma")
) %>%
  select(parameter, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`)

cat("\n=== TABLE 4B: Selection parameters ===\n")
table4_selection %>% as_tibble() %>% print(n = nrow(.))

# -----------------------------
# OPTIONAL: export as csv files
# -----------------------------
write_csv(table1_desc, "table1_descriptives.csv")
write_csv(table2_measurement, "table2_measurement.csv")
write_csv(table3_determinants, "table3a_determinants.csv")
write_csv(table3_participation, "table3b_participation.csv")
write_csv(table4_wage_compare, "table4a_wages.csv")
write_csv(table4_selection, "table4b_selection.csv")

# =========================================================
# FAST HYBRID REDUCED MSL - PIECE 1
# Fix measurement block, estimate bargaining + participation
# =========================================================

library(tidyverse)

# Same variable order as before, but normalise on a strong-loading indicator
Y_names <- c(
  "yieldfor.proj.m",
  "s_m",
  setdiff(measurement_vars, "yieldfor.proj.m")
)

X_names <- c(trait_vars, z_vars)

Zp_names <- c(
  "age.f", "educ.f",
  "duration", "abouttobeparents",
  "d_educ", "d_age", "both_m"
)

msl_df <- df %>%
  select(hid, D_f, all_of(Y_names), all_of(X_names), all_of(Zp_names)) %>%
  drop_na()

cat("Hybrid MSL sample size:", nrow(msl_df), "\n")

Y_raw <- as.matrix(msl_df[, Y_names])
Y <- scale(Y_raw)

X <- as.matrix(msl_df[, X_names])
Zp <- as.matrix(msl_df[, Zp_names])
D <- as.numeric(msl_df$D_f)

N  <- nrow(Y)
J  <- ncol(Y)
Kx <- ncol(X)
Kz <- ncol(Zp)

cat("N =", N, " J =", J, " Kx =", Kx, " Kz =", Kz, "\n")

# Re-estimate a one-factor model on THIS exact MSL sample
fa_msl <- factanal(
  x = Y,
  factors = 1,
  scores = "none"
)

load_msl <- tibble(
  variable = rownames(fa_msl$loadings),
  loading = as.numeric(fa_msl$loadings[, 1]),
  uniqueness = as.numeric(fa_msl$uniquenesses)
)

# orient sign so yieldfor.proj.m is positive
ref_loading <- load_msl %>%
  filter(variable == "yieldfor.proj.m") %>%
  pull(loading)

flip <- ifelse(ref_loading < 0, -1, 1)

load_msl <- load_msl %>%
  mutate(
    loading = loading * flip,
    sigma_y = sqrt(uniqueness)
  )

print(load_msl, n = nrow(load_msl))

# Fixed measurement parameters
beta_fix <- rep(0, J)  # because Y is standardised
lambda_fix <- load_msl$loading[match(Y_names, load_msl$variable)]
sigma_y_fix <- load_msl$sigma_y[match(Y_names, load_msl$variable)]

cat("\nFixed measurement loadings:\n")
print(tibble(variable = Y_names, lambda_fix = lambda_fix, sigma_y_fix = sigma_y_fix), n = J)



# =========================================================
# FAST HYBRID REDUCED MSL - PIECE 2
# =========================================================

set.seed(123)
R <- 20
draws <- matrix(rnorm(N * R), nrow = N, ncol = R)

log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

# Parameters to estimate now:
# gamma      : Kx bargaining-equation coefficients
# log_sigma_u: sd of latent bargaining shock
# a0         : participation intercept
# a_mu       : effect of latent bargaining on participation
# a_z        : Kz participation coefficients
unpack_theta_fast <- function(theta, Kx, Kz) {
  idx <- 1
  
  gamma <- theta[idx:(idx + Kx - 1)]
  idx <- idx + Kx
  
  log_sigma_u <- theta[idx]
  sigma_u <- exp(log_sigma_u)
  idx <- idx + 1
  
  a0 <- theta[idx]
  idx <- idx + 1
  
  a_mu <- theta[idx]
  idx <- idx + 1
  
  a_z <- theta[idx:(idx + Kz - 1)]
  
  list(
    gamma = gamma,
    sigma_u = sigma_u,
    a0 = a0,
    a_mu = a_mu,
    a_z = a_z
  )
}

neg_loglik_fast <- function(theta, Y, X, Zp, D, draws, beta_fix, lambda_fix, sigma_y_fix) {
  N  <- nrow(Y)
  R  <- ncol(draws)
  Kx <- ncol(X)
  Kz <- ncol(Zp)
  
  p <- unpack_theta_fast(theta, Kx, Kz)
  
  mu_mean <- as.vector(X %*% p$gamma)
  MU <- mu_mean + p$sigma_u * draws
  
  ll_i <- numeric(N)
  
  for (i in 1:N) {
    log_terms <- numeric(R)
    
    for (r in 1:R) {
      mu_ir <- MU[i, r]
      
      # fixed measurement block
      mean_y <- beta_fix + lambda_fix * mu_ir
      log_fy <- sum(dnorm(Y[i, ], mean = mean_y, sd = sigma_y_fix, log = TRUE))
      
      # participation block
      index_ir <- p$a0 + p$a_mu * mu_ir + sum(Zp[i, ] * p$a_z)
      p_work <- pnorm(index_ir)
      p_work <- min(max(p_work, 1e-12), 1 - 1e-12)
      
      log_gd <- ifelse(D[i] == 1, log(p_work), log(1 - p_work))
      
      log_terms[r] <- log_fy + log_gd
    }
    
    ll_i[i] <- log_mean_exp(log_terms)
  }
  
  -sum(ll_i)
}




# =========================================================
# FAST HYBRID REDUCED MSL - PIECE 3
# =========================================================

theta_start_fast <- c(
  rep(0, Kx),   # gamma
  log(1),       # log_sigma_u
  0,            # a0
  0.1,          # a_mu
  rep(0, Kz)    # a_z
)

cat("Number of free parameters in fast version:", length(theta_start_fast), "\n")

msl_fit_fast <- optim(
  par = theta_start_fast,
  fn = neg_loglik_fast,
  Y = Y,
  X = X,
  Zp = Zp,
  D = D,
  draws = draws,
  beta_fix = beta_fix,
  lambda_fix = lambda_fix,
  sigma_y_fix = sigma_y_fix,
  method = "BFGS",
  control = list(maxit = 120, trace = 1, REPORT = 5),
  hessian = FALSE
)

cat("\nConvergence code:", msl_fit_fast$convergence, "\n")
cat("Negative log-likelihood:", msl_fit_fast$value, "\n")

p_fast <- unpack_theta_fast(msl_fit_fast$par, Kx, Kz)

participation_estimates_fast <- tibble(
  term = c("(Intercept)", "mu", Zp_names),
  estimate = c(p_fast$a0, p_fast$a_mu, p_fast$a_z)
)

bargaining_estimates_fast <- tibble(
  term = X_names,
  gamma_hat = p_fast$gamma
)

cat("\n=== Fast hybrid MSL: participation estimates ===\n")
print(participation_estimates_fast, n = nrow(participation_estimates_fast))

cat("\n=== Fast hybrid MSL: bargaining estimates (first 15) ===\n")
print(bargaining_estimates_fast, n = 15)

# =========================================================
# STEP 8B: Approximate standard errors for reduced MSL
# =========================================================

# We use the converged parameter vector from msl_fit_fast
theta_hat_fast <- msl_fit_fast$par

# 1. Numerical Hessian at the optimum
# Base R version first
t0_hess <- Sys.time()

hess_fast <- optimHess(
  par = theta_hat_fast,
  fn = neg_loglik_fast,
  Y = Y,
  X = X,
  Zp = Zp,
  D = D,
  draws = draws,
  beta_fix = beta_fix,
  lambda_fix = lambda_fix,
  sigma_y_fix = sigma_y_fix
)

t1_hess <- Sys.time()

cat("\n=== Hessian runtime ===\n")
print(t1_hess - t0_hess)

# 2. Try to invert Hessian
vcov_fast <- tryCatch(
  solve(hess_fast),
  error = function(e) NULL
)

if (is.null(vcov_fast)) {
  cat("\nHessian inversion failed.\n")
} else {
  cat("\nHessian inversion succeeded.\n")
}

# 3. Standard errors
if (!is.null(vcov_fast)) {
  se_fast <- sqrt(diag(vcov_fast))
} else {
  se_fast <- rep(NA_real_, length(theta_hat_fast))
}

# 4. Build a clean parameter table
param_names_fast <- c(
  paste0("gamma_", X_names),
  "log_sigma_u",
  "a0",
  "a_mu",
  paste0("a_z_", Zp_names)
)

msl_se_table <- tibble(
  parameter = param_names_fast,
  estimate = theta_hat_fast,
  std_error = se_fast
) %>%
  mutate(
    z_value = estimate / std_error,
    p_value = 2 * (1 - pnorm(abs(z_value)))
  )

cat("\n=== Approximate SE table for reduced MSL ===\n")
print(msl_se_table, n = nrow(msl_se_table))

# 5. Participation coefficients only
participation_se_fast <- msl_se_table %>%
  filter(parameter %in% c("a0", "a_mu", paste0("a_z_", Zp_names))) %>%
  mutate(
    term = c("(Intercept)", "mu", Zp_names)
  ) %>%
  select(term, estimate, std_error, z_value, p_value)

cat("\n=== Reduced MSL participation coefficients with approximate SEs ===\n")
print(participation_se_fast, n = nrow(participation_se_fast))

# 6. Bargaining-equation coefficients only
bargaining_se_fast <- msl_se_table %>%
  filter(grepl("^gamma_", parameter)) %>%
  mutate(
    term = sub("^gamma_", "", parameter)
  ) %>%
  select(term, estimate, std_error, z_value, p_value)

cat("\n=== Reduced MSL bargaining coefficients with approximate SEs ===\n")
print(bargaining_se_fast, n = nrow(bargaining_se_fast))

# 7. Delta-method SE for sigma_u = exp(log_sigma_u)
sigma_u_row <- msl_se_table %>%
  filter(parameter == "log_sigma_u")

if (nrow(sigma_u_row) == 1 && !is.na(sigma_u_row$std_error)) {
  sigma_u_hat <- exp(sigma_u_row$estimate)
  sigma_u_se <- sigma_u_hat * sigma_u_row$std_error
  
  sigma_u_table <- tibble(
    term = "sigma_u",
    estimate = sigma_u_hat,
    std_error = sigma_u_se,
    z_value = sigma_u_hat / sigma_u_se,
    p_value = 2 * (1 - pnorm(abs(z_value)))
  )
  
  cat("\n=== Delta-method SE for sigma_u ===\n")
  print(sigma_u_table)
}




