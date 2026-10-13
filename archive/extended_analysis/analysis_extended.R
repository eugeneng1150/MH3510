# MH3510: modelling average annual daily traffic (AADT)
# Run from this project folder: Rscript analysis.R
# Uses base R and MASS (a recommended package distributed with R).
# All attempted models, including unsuccessful ones, are retained.
# Outputs: results/analysis_output.txt, results/interpretation.txt, tables,
# fitted_models.rds, and a diagnostic figure for EVERY candidate model.

main <- function() {
  options(width = 110, digits = 6)
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Install the recommended MASS package before running this script.")
  }
  dir.create("results", showWarnings = FALSE)
  dir.create("results/plots", showWarnings = FALSE)
  sink("results/analysis_output.txt", split = TRUE)
  on.exit(sink(), add = TRUE)
  heading <- function(x) cat("\n\n", x, "\n", strrep("=", nchar(x)), "\n", sep = "")
  save_plot <- function(name, draw, width = 1500, height = 1000) {
    png(file.path("results/plots", name), width = width, height = height, res = 140)
    on.exit(dev.off())
    draw()
  }

  heading("1. DATA AND GRAPHICAL EXPLORATION")
  raw <- read.table("aadt.txt", header = FALSE)
  stopifnot(ncol(raw) >= 5, nrow(raw) > 5)
  d <- raw[, 1:5]
  names(d) <- c("aadt", "population", "lanes", "width", "control_original")
  stopifnot(!anyNA(d), all(is.finite(as.matrix(d))),
            all(d$aadt > 0), all(d$population > 0), all(d$lanes > 0),
            all(d$control_original %in% c(1, 2)))
  d$access <- as.integer(d$control_original == 1)
  cat("Observations:", nrow(d), "\nOnly original columns 1-5 are used.\n")
  cat("access = 1 means access control; access = 0 means no access control.\n")
  print(summary(d))
  print(with(d, table(lanes, control_original)))
  print(cor(d[, c("aadt", "population", "lanes", "width", "access")]))
  write.csv(d, "results/analysis_data.csv", row.names = FALSE)
  save_plot("01_scatterplot_matrix.png", function() {
    pairs(d[, c("aadt", "population", "lanes", "width", "access")],
          pch = 19, cex = 0.55, col = adjustcolor("#176b87", 0.55),
          main = "AADT and the four candidate predictors")
  }, 1400, 1400)
  save_plot("02_transformed_scatterplot_matrix.png", function() {
    pairs(data.frame(fifth_root_aadt = d$aadt^0.2,
                     log_population = log(d$population), log_lanes = log(d$lanes),
                     width = d$width, access = d$access),
          pch = 19, cex = 0.55, col = adjustcolor("#176b87", 0.55),
          main = "Relationships after transformation (natural logarithms)")
  }, 1400, 1400)

  heading("2. ALL CANDIDATE MODELS: PRESERVED EXPLORATION")
  # These candidates reproduce the exploratory sequence. Log population and
  # log lanes allow diminishing changes on the transformed response scale.
  specifications <- list(
    M1_raw = aadt ~ population + lanes + width + access,
    M2_log_y = log(aadt) ~ population + lanes + width + access,
    M3_log_y_log_population = log(aadt) ~ log(population) + lanes + width + access,
    M4_log_y_log_population_lanes = log(aadt) ~ log(population) + log(lanes) + width + access,
    M5_log_y_quadratic_lanes = log(aadt) ~ log(population) + lanes + I(lanes^2) + width + access,
    M6_log_y_reduced = log(aadt) ~ log(population) + log(lanes) + access,
    M7_fourth_root_full = I(aadt^0.25) ~ log(population) + log(lanes) + width + access,
    M8_square_root_full = I(aadt^0.5) ~ log(population) + log(lanes) + width + access,
    M9_fourth_root_reduced = I(aadt^0.25) ~ log(population) + log(lanes) + access,
    M10_fifth_root_full = I(aadt^0.2) ~ log(population) + log(lanes) + width + access,
    M11_fifth_root_reduced = I(aadt^0.2) ~ log(population) + log(lanes) + access
  )
  scales <- c("raw", rep("log", 5), "fourth_root", "square_root", "fourth_root",
              "fifth_root", "fifth_root")
  names(scales) <- names(specifications)
  models <- lapply(specifications, function(f) lm(formula = f, data = d))
  for (id in names(models)) {
    cat("\n", id, "\n", sep = "")
    print(summary(models[[id]]))
  }
  saveRDS(models, "results/fitted_models.rds")

  heading("3. RESIDUAL CHECKS FOR EVERY MODEL")
  # Studentized/Koenker BP diagnostic: n*R^2 from squared residuals regressed
  # on the actual model predictors; asymptotic chi-square reference.
  variance_check <- function(m) {
    x <- model.matrix(m)
    auxiliary <- lm.fit(x, residuals(m)^2)
    r2 <- 1 - sum(auxiliary$residuals^2) /
      sum((residuals(m)^2 - mean(residuals(m)^2))^2)
    statistic <- nobs(m) * r2
    c(statistic = statistic, df = ncol(x) - 1,
      p_value = pchisq(statistic, ncol(x) - 1, lower.tail = FALSE))
  }
  comparison <- do.call(rbind, lapply(names(models), function(id) {
    m <- models[[id]]
    data.frame(model = id, response_scale = scales[[id]], parameters = length(coef(m)),
               r_squared = summary(m)$r.squared,
               adjusted_r_squared = summary(m)$adj.r.squared,
               residual_sd = sigma(m), shapiro_p = shapiro.test(residuals(m))$p.value,
               variance_test_p = unname(variance_check(m)["p_value"]),
               max_cooks_distance = max(cooks.distance(m)))
  }))
  for (id in names(models)) {
    m <- models[[id]]
    cat("\n", id, "\n", sep = "")
    print(shapiro.test(residuals(m)))
    print(variance_check(m))
    save_plot(paste0(id, "_diagnostics.png"), function() {
      par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
      plot(m, which = 1, ask = FALSE, sub.caption = "")
      plot(m, which = 2, ask = FALSE, sub.caption = "")
      plot(m, which = 3, ask = FALSE, sub.caption = "")
      plot(log(d$population), residuals(m), pch = 19, cex = 0.6,
           xlab = "Log population", ylab = "Residual", main = "Residuals vs population")
      abline(h = 0, lty = 2)
      plot(d$lanes, residuals(m), pch = 19, cex = 0.6,
           xlab = "Lanes", ylab = "Residual", main = "Residuals vs lanes")
      abline(h = 0, lty = 2)
      plot(m, which = 4, ask = FALSE, sub.caption = "")
      mtext(id, outer = TRUE, font = 2)
    })
  }
  print(comparison, row.names = FALSE)
  cat("R-squared and residual SD are NOT directly comparable across response transformations.\n")
  cat("Diagnostic p-values are supporting evidence, not proof that assumptions hold.\n")
  cat("No time or spatial ordering is documented. Row number is not assumed to be time.\n")
  cat("A Durbin-Watson significance test would not establish independence here.\n")

  heading("4. TRANSFORMATION, MODEL COMPARISON AND SELECTION")
  # Box-Cox includes the transformation Jacobian; it can compare response
  # powers, unlike ordinary AIC values from different response scales.
  # Hold the predictor specification fixed while choosing a response power.
  boxcox_base <- lm(aadt ~ log(population) + log(lanes) + width + access,
                   data = d, y = TRUE, qr = TRUE)
  bc <- MASS::boxcox(boxcox_base, lambda = seq(-0.5, 1, by = 0.01), plotit = FALSE)
  supported <- bc$x[bc$y >= max(bc$y) - qchisq(0.95, 1) / 2]
  cat("Box-Cox maximizing lambda:", bc$x[which.max(bc$y)], "\n")
  cat("Approximate 95% profile support (grid):", range(supported), "\n")
  cat("Compare powers 0.20 and 0.25 using residuals and original-scale prediction error.\n")
  cat("Power 0.25 is just outside the 95% likelihood support; Box-Cox alone favors 0.20.\n")
  write.csv(data.frame(lambda = bc$x, profile_log_likelihood = bc$y),
            "results/boxcox_profile.csv", row.names = FALSE)
  save_plot("03_boxcox_profile.png", function() {
    plot(bc$x, bc$y, type = "l", lwd = 2, xlab = "Response power (lambda)",
         ylab = "Profile log likelihood", main = "Box-Cox: response transformation")
    abline(v = c(0, 0.2, 0.25, 0.5), col = c("grey50", "#176b87", "grey50", "grey50"), lty = 2)
    abline(h = max(bc$y) - qchisq(0.95, 1) / 2, lty = 3)
    legend("bottom", "Chosen convenient power: 0.20", lty = 2, col = "#176b87", bty = "n")
  })
  cat("\nWidth deletion test on fourth-root response:\n")
  width_test <- anova(models$M9_fourth_root_reduced, models$M7_fourth_root_full)
  print(width_test)
  cat("\nWidth deletion test for the earlier log-response candidate:\n")
  print(anova(models$M6_log_y_reduced, models$M4_log_y_log_population_lanes))
  cat("\nWidth deletion test on fifth-root response:\n")
  print(anova(models$M11_fifth_root_reduced, models$M10_fifth_root_full))
  cat("\nQuadratic lanes test within log-response models:\n")
  print(anova(models$M3_log_y_log_population, models$M5_log_y_quadratic_lanes))
  final <- models$M11_fifth_root_reduced
  curvature <- update(final, . ~ . + I(log(population)^2) + I(log(lanes)^2))
  cat("\nAdditional curvature check for selected fifth-root model:\n")
  print(anova(final, curvature))
  saveRDS(list(boxcox_base = boxcox_base, curvature_check = curvature),
          "results/supplementary_models.rds")

  # Compare predictions on the SAME original AADT scale. Retransformation
  # corrections are estimated ONLY from each training fold's residuals.
  original_mean <- function(m, newdata, scale) {
    mu <- as.numeric(predict(m, newdata = newdata))
    e <- residuals(m)
    if (scale == "raw") return(mu)
    if (scale == "log") return(exp(mu) * mean(exp(e)))
    if (scale == "square_root") return(vapply(mu, function(z) mean((z + e)^2), numeric(1)))
    if (scale == "fifth_root") return(vapply(mu, function(z) mean((z + e)^5), numeric(1)))
    vapply(mu, function(z) mean((z + e)^4), numeric(1))
  }
  set.seed(3510)
  cv_rows <- list()
  for (repeat_id in 1:10) {
    folds <- sample(rep(1:10, length.out = nrow(d)))
    for (id in names(models)) {
      p <- numeric(nrow(d))
      for (fold in 1:10) {
        train <- folds != fold
        fit <- lm(specifications[[id]], data = d[train, ])
        p[!train] <- original_mean(fit, d[!train, ], scales[[id]])
      }
      cv_rows[[length(cv_rows) + 1]] <- data.frame(
        repetition = repeat_id, model = id,
        rmse = sqrt(mean((p - d$aadt)^2)), mae = mean(abs(p - d$aadt)))
    }
  }
  cv <- do.call(rbind, cv_rows)
  write.csv(cv, "results/cross_validation_repetitions.csv", row.names = FALSE)
  cv_summary <- do.call(rbind, lapply(names(models), function(id) {
    rows <- cv[cv$model == id, ]
    data.frame(model = id, cv_rmse = mean(rows$rmse), cv_rmse_sd = sd(rows$rmse),
               cv_mae = mean(rows$mae), cv_mae_sd = sd(rows$mae))
  }))
  comparison <- merge(comparison, cv_summary, by = "model", sort = FALSE)
  write.csv(comparison, "results/model_comparison.csv", row.names = FALSE)
  print(comparison, row.names = FALSE)
  cat("CV: 10 repetitions of 10-fold CV, identical splits for all candidates.\n")
  cat("SD measures variation between repetitions, not a confidence interval.\n")
  cat("Transformations were explored using these data; this is internal validation,\n")
  cat("not an independent test set or an unbiased assessment of the selection process.\n")
  cat("Final selection balances original-scale CV, diagnostics, and simplicity; see interpretation.txt.\n")
  save_plot("04_prediction_error_comparison.png", function() {
    par(mfrow = c(1, 2), mar = c(5, 12, 4, 1))
    labels <- sub("^M[0-9]+_", "", comparison$model)
    barplot(comparison$cv_rmse, names.arg = labels, horiz = TRUE, las = 1,
            cex.names = 0.7, col = "#176b87", xlab = "Vehicles/day", main = "Cross-validated RMSE")
    barplot(comparison$cv_mae, names.arg = labels, horiz = TRUE, las = 1,
            cex.names = 0.7, col = "#419d78", xlab = "Vehicles/day", main = "Cross-validated MAE")
  }, 1800, 1000)

  heading("FINAL MODEL: INFLUENCE, COLLINEARITY AND SENSITIVITY")
  print(summary(final))
  print(confint(final)) # Ordinary coefficient intervals, not Bonferroni limits.
  x <- model.matrix(final)[, -1, drop = FALSE]
  vif <- sapply(seq_len(ncol(x)), function(j) {
    1 / (1 - summary(lm(x[, j] ~ x[, -j, drop = FALSE]))$r.squared)
  })
  names(vif) <- colnames(x)
  print(vif)
  influence <- data.frame(row = seq_len(nrow(d)), d,
                          studentized_residual = rstudent(final),
                          leverage = hatvalues(final), cooks_distance = cooks.distance(final))
  influence$flag <- influence$cooks_distance > 4 / nrow(d) |
    influence$leverage > 2 * length(coef(final)) / nrow(d) |
    abs(influence$studentized_residual) > 3
  write.csv(influence, "results/influence_diagnostics.csv", row.names = FALSE)
  print(influence[influence$flag, ], row.names = FALSE)
  largest <- which.max(cooks.distance(final))
  sensitivity <- lm(formula(final), data = d[-largest, ])
  cat("Sensitivity only: refit omitting highest Cook's-distance row", largest, "\n")
  print(cbind(all_observations = coef(final), sensitivity_only = coef(sensitivity)))
  cat("All observations remain in the chosen model; flags alone do not justify deletion.\n")
  # HC3 offers an additional check on coefficient inference if variance is imperfect.
  design <- model.matrix(final)
  bread <- solve(crossprod(design))
  weighted_design <- design * as.numeric(residuals(final) / (1 - hatvalues(final)))
  hc3 <- bread %*% crossprod(weighted_design) %*% bread
  hc3_se <- sqrt(diag(hc3))
  print(data.frame(estimate = coef(final), hc3_se = hc3_se,
                   approximate_t_p = 2 * pt(-abs(coef(final) / hc3_se), df.residual(final))))
  save_plot("05_final_residuals_all_predictors.png", function() {
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
    for (v in c("population", "lanes", "width")) {
      plot(d[[v]], residuals(final), xlab = v, ylab = "Residual", pch = 19, cex = 0.6,
           main = paste("Residuals vs", v))
      abline(h = 0, lty = 2)
    }
    boxplot(residuals(final) ~ d$access, xlab = "Access control (0=no, 1=yes)",
            ylab = "Residual", main = "Residuals by access control")
    abline(h = 0, lty = 2)
    plot(residuals(final), xlab = "File row (not known time order)", ylab = "Residual",
         pch = 19, cex = 0.6, main = "Row order: descriptive only")
    abline(h = 0, lty = 2)
    plot(fitted(final), residuals(final), xlab = "Fitted fifth-root AADT", ylab = "Residual",
         pch = 19, cex = 0.6, main = "Residuals vs fitted")
    abline(h = 0, lty = 2)
  })

  heading("5. REQUESTED PREDICTION AND INTERVALS")
  new_road <- data.frame(population = 50000, lanes = 3, width = 60,
                         control_original = 2, access = 0)
  print(new_road)
  cat("Original x4=2 is recoded to access=0. Width is omitted by the final model.\n")
  ci_transformed <- predict(final, new_road, interval = "confidence", level = 0.95)
  pi_transformed <- predict(final, new_road, interval = "prediction", level = 0.95)
  print(ci_transformed)
  print(pi_transformed)
  stopifnot(all(ci_transformed > 0), all(pi_transformed > 0))
  typical <- ci_transformed^5
  prediction_interval <- pi_transformed[, c("lwr", "upr"), drop = FALSE]^5
  mean_estimate <- original_mean(final, new_road, "fifth_root")
  # Inverse-transformed mean on the root scale is not mean AADT (Jensen's
  # inequality). Average over residuals to estimate original-scale mean AADT.
  # Pairs bootstrap gives an approximate CI for that corrected mean, conditional
  # on the selected formula. It does not repeat model/transform selection.
  set.seed(3511)
  boot_means <- replicate(2000, {
    sample_rows <- sample.int(nrow(d), replace = TRUE)
    m <- lm(formula(final), data = d[sample_rows, ])
    original_mean(m, new_road, "fifth_root")
  })
  mean_ci <- quantile(boot_means, c(0.025, 0.975))
  prediction <- data.frame(
    quantity = c("Mean AADT (residual-smearing correction)",
                 "Typical AADT (inverse transformed fitted value)",
                 "Individual road AADT prediction interval"),
    estimate = c(mean_estimate, typical[1, "fit"], NA_real_),
    lower_95 = c(mean_ci[1], typical[1, "lwr"], prediction_interval[1, "lwr"]),
    upper_95 = c(mean_ci[2], typical[1, "upr"], prediction_interval[1, "upr"]))
  print(prediction, row.names = FALSE)
  write.csv(prediction, "results/prediction.csv", row.names = FALSE)
  write.csv(data.frame(mean_aadt = boot_means), "results/bootstrap_mean_predictions.csv", row.names = FALSE)
  cat("The typical-value CI is NOT a confidence interval for mean AADT.\n")
  cat("The prediction interval uses approximately normal, constant-variance root-scale errors.\n")
  cat("The mean CI uses 2000 pairs-bootstrap resamples and a common residual distribution.\n")
  cat("Both are conditional on the selected model and assumed independent road sections.\n")
  cat("Three lanes is unobserved: the prediction interpolates between two and four lanes.\n")
  cat("Sensitivity-only mean prediction:", original_mean(sensitivity, new_road, "fifth_root"), "\n")

  heading("REPRODUCIBILITY")
  print(sessionInfo())
  cat("\nAll outputs were written under results/. Source data were not modified.\n")
}

main()
