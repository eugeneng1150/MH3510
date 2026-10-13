# MH3510 AADT project - course-aligned analysis, using only base R.
# Run: Rscript analysis.R (with this project folder as the working directory).
# Read analysis_notes.txt for explanations, lecture references and model decisions.
# Results go to coursework_results/. Every candidate model is retained.
# Earlier extended methods are preserved in archive/extended_analysis/.

run_analysis <- function() {
  options(width = 110, digits = 6, useFancyQuotes = FALSE)
  dir.create("coursework_results", showWarnings = FALSE)
  dir.create("coursework_results/plots", showWarnings = FALSE)
  sink("coursework_results/R_output.txt", split = TRUE)
  on.exit(sink(), add = TRUE)

  # These helpers only format output and save plots; they do not fit models.
  section <- function(title) cat("\n\n", title, "\n", strrep("=", nchar(title)), "\n", sep = "")
  save_plot <- function(filename, draw, width = 1600, height = 1100) {
    png(file.path("coursework_results/plots", filename), width, height, res = 140)
    on.exit(dev.off())
    draw()
  }

  section("STEP 1: LOAD AND EXPLORE THE DATA")
  raw <- read.table("aadt.txt", header = FALSE)
  d <- raw[, 1:5]
  names(d) <- c("aadt", "population", "lanes", "width", "control_original")
  stopifnot(!anyNA(d), all(is.finite(as.matrix(d))), all(d$aadt > 0),
            all(d$population > 0), all(d$lanes > 0),
            all(d$control_original %in% c(1, 2)))
  d$access <- ifelse(d$control_original == 1, 1, 0)
  cat("Number of road sections:", nrow(d), "\n")
  cat("access=1: access control; access=0: no access control.\n")
  print(summary(d))
  print(with(d, table(lanes, control_original)))
  print(cor(d[, c("aadt", "population", "lanes", "width", "access")]))
  write.csv(d, "coursework_results/analysis_data.csv", row.names = FALSE)
  save_plot("01_scatterplot_matrix.png", function() {
    pairs(d[, c("aadt", "population", "lanes", "width", "access")],
          pch = 19, cex = 0.6, col = adjustcolor("#176b87", 0.55),
          main = "AADT and the four candidate predictors")
  }, 1400, 1400)

  section("STEP 2: FIT THE INITIAL MULTIPLE REGRESSION")
  M1 <- lm(aadt ~ population + lanes + width + access, data = d)
  print(summary(M1))
  cat("Negative fitted AADT values:", sum(fitted(M1) < 0), "\n")
  cat("A significant overall F-test does not guarantee satisfactory residuals.\n")

  section("STEP 3: CHECK RESIDUALS, TRY TRANSFORMATIONS, AND RECHECK")
  # Chapter 2, slides 51-59: transform x for curvature and y for uneven spread.
  # Each model below is kept, even if it does not improve the diagnostics.
  M2 <- lm(log(aadt) ~ population + lanes + width + access, data = d)
  M3 <- lm(log(aadt) ~ log(population) + lanes + width + access, data = d)
  M4 <- lm(log(aadt) ~ log(population) + log(lanes) + width + access, data = d)
  # Chapter 4 Part II: polynomial terms are still linear regression coefficients.
  M5 <- lm(log(aadt) ~ log(population) + lanes + I(lanes^2) + width + access, data = d)
  M6 <- lm(log(aadt) ~ log(population) + log(lanes) + access, data = d)
  # A fourth root is sqrt(sqrt(y)); a fifth root is y^(1/5).
  M7 <- lm(I(aadt^0.25) ~ log(population) + log(lanes) + width + access, data = d)
  M8 <- lm(sqrt(aadt) ~ log(population) + log(lanes) + width + access, data = d)
  M9 <- lm(I(aadt^0.25) ~ log(population) + log(lanes) + access, data = d)
  M10 <- lm(I(aadt^0.2) ~ log(population) + log(lanes) + width + access, data = d)
  M11 <- lm(I(aadt^0.2) ~ log(population) + log(lanes) + access, data = d)
  models <- list(M1_raw = M1, M2_log_y = M2, M3_log_y_log_population = M3,
                 M4_log_y_log_population_lanes = M4, M5_log_y_quadratic_lanes = M5,
                 M6_log_y_reduced = M6, M7_fourth_root_full = M7,
                 M8_square_root_full = M8, M9_fourth_root_reduced = M9,
                 M10_fifth_root_full = M10, M11_fifth_root_reduced = M11)
  response_scales <- c("AADT", rep("log(AADT)", 5), "AADT^0.25", "sqrt(AADT)",
                       "AADT^0.25", "AADT^0.2", "AADT^0.2")
  fit_table <- data.frame()
  for (i in seq_along(models)) {
    m <- models[[i]]
    id <- names(models)[i]
    cat("\n", id, "\n", sep = "")
    print(formula(m))
    print(summary(m))
    fit_table <- rbind(fit_table, data.frame(
      model = id, response = response_scales[i], number_of_coefficients = length(coef(m)),
      r_squared = summary(m)$r.squared, adjusted_r_squared = summary(m)$adj.r.squared,
      residual_sd = sigma(m)))
    save_plot(paste0(id, "_residuals.png"), function() {
      par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
      plot(m, which = 1, sub.caption = "", ask = FALSE)
      plot(m, which = 2, sub.caption = "", ask = FALSE)
      for (v in c("population", "lanes", "width")) {
        plot(d[[v]], residuals(m), xlab = v, ylab = "Residual", pch = 19, cex = 0.6,
             main = paste("Residuals vs", v))
        abline(h = 0, lty = 2)
      }
      boxplot(residuals(m) ~ d$access, xlab = "Access control (0=no, 1=yes)",
              ylab = "Residual", main = "Residuals by access control")
      abline(h = 0, lty = 2)
      mtext(id, outer = TRUE, font = 2)
    })
  }
  print(fit_table, row.names = FALSE)
  write.csv(fit_table, "coursework_results/model_fit_summary.csv", row.names = FALSE)
  saveRDS(models, "coursework_results/all_candidate_models.rds")
  cat("Compare R-squared, adjusted R-squared and residual SD only on the SAME response scale.\n")
  cat("Across different response scales, judge residual patterns and explain the transformation.\n")

  section("STEP 4: FULL/REDUCED MODEL F-TESTS AND FINAL MODEL")
  cat("Log model: test whether road width is needed, H0: beta_width=0.\n")
  print(anova(M6, M4))
  cat("Log model: test whether the quadratic lane term is needed.\n")
  print(anova(M3, M5))
  cat("Fourth-root model: test whether road width is needed.\n")
  print(anova(M9, M7))
  cat("Fifth-root model: test whether road width is needed.\n")
  print(anova(M11, M10))
  # Never use anova() to compare different transformed responses.
  # Fifth-root residuals have a relatively flat trend, much more stable spread
  # than the raw/square-root models, and a reasonable Q-Q plot.
  # Fourth-root models are also reasonable alternatives, not failed fits.
  final_model <- M11
  cat("\nSelected model:\n")
  print(formula(final_model))
  print(summary(final_model))
  cat("Ordinary individual coefficient confidence intervals:\n")
  print(confint(final_model))
  cat("No confidence region or Bonferroni limits are calculated.\n")
  cat("Optional curvature check (Chapter 4 Part II polynomial regression):\n")
  curved_model <- lm(I(aadt^0.2) ~ log(population) + log(lanes) + access +
                      I(log(population)^2) + I(log(lanes)^2), data = d)
  print(anova(final_model, curved_model))

  # Chapter 2 slides 49 and 53: inspect influence and report sensitivity.
  # Do not remove a road merely because its residual is large.
  influential_row <- which.max(cooks.distance(final_model))
  cat("Largest Cook's-distance row:", influential_row, "\n")
  print(d[influential_row, ])
  sensitivity_model <- lm(formula(final_model), data = d[-influential_row, ])
  print(cbind(all_rows = coef(final_model), omit_one_for_sensitivity = coef(sensitivity_model)))
  cat("Final model retains ALL road sections. Sensitivity fit is only a check.\n")
  write.csv(data.frame(row = seq_len(nrow(d)),
                       standardized_residual = residuals(final_model) / sigma(final_model),
                       cooks_distance = cooks.distance(final_model)),
            "coursework_results/influence_check.csv", row.names = FALSE)
  saveRDS(list(curvature_check = curved_model, sensitivity_only = sensitivity_model),
          "coursework_results/additional_checks.rds")
  save_plot("02_final_influence_and_order.png", function() {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    plot(final_model, which = 4, sub.caption = "", ask = FALSE)
    plot(residuals(final_model), xlab = "File row (not known time order)",
         ylab = "Residual", pch = 19, main = "Residuals in file order")
    abline(h = 0, lty = 2)
  }, 1500, 650)
  # Reference tutorial gives this statistic; order is not documented here.
  e <- residuals(final_model)
  durbin_watson <- sum(diff(e)^2) / sum(e^2)
  cat("Durbin-Watson statistic in FILE order:", durbin_watson, "\n")
  cat("No time/spatial order is provided, so no serial-correlation conclusion is drawn.\n")
  cat("Independence remains an assumption; a row-order plot cannot establish it.\n")

  section("STEP 5: PREDICTION FOR THE REQUESTED ROAD")
  new_road <- data.frame(population = 50000, lanes = 3, width = 60, access = 0)
  cat("Original control code x4=2 means NO access control, hence access=0.\n")
  cat("Width=60 is recorded, but width is excluded from the selected model.\n")
  print(new_road)
  # Standard predict() intervals, as in the reference tutorial and Chapter 4.
  mean_root_interval <- predict(final_model, new_road, interval = "confidence", level = 0.95)
  individual_root_interval <- predict(final_model, new_road, interval = "prediction", level = 0.95)
  cat("95% confidence interval for mean FIFTH-ROOT AADT:\n")
  print(mean_root_interval)
  cat("95% prediction interval for one road, on the FIFTH-ROOT scale:\n")
  print(individual_root_interval)
  # Fifth power converts the individual prediction interval to AADT units.
  # It does NOT turn a mean-root confidence interval into a mean-AADT interval.
  point_aadt <- mean_root_interval[1, "fit"]^5
  individual_aadt_interval <- individual_root_interval[1, c("lwr", "upr")]^5
  cat("Back-transformed fitted value (typical AADT):", point_aadt, "vehicles/day\n")
  cat("95% individual-road prediction interval in vehicles/day:\n")
  print(individual_aadt_interval)
  cat("This point is a typical/approximate median prediction, NOT the arithmetic mean AADT.\n")
  cat("A mean-AADT CI requires an extra retransformation method; it is not claimed here.\n")
  cat("There are no 3-lane observations; 3 lanes interpolates between 2 and 4 lanes.\n")
  cat("Sensitivity-only back-transformed point:", as.numeric(predict(sensitivity_model, new_road))^5, "\n")
  cat("Fourth-root alternative's back-transformed point:", as.numeric(predict(M9, new_road))^4, "\n")
  predictions <- data.frame(
    quantity = c("Mean fifth-root AADT: confidence interval",
                 "Typical AADT point with individual-road prediction interval"),
    estimate = c(mean_root_interval[1, "fit"], point_aadt),
    lower_95 = c(mean_root_interval[1, "lwr"], individual_aadt_interval["lwr"]),
    upper_95 = c(mean_root_interval[1, "upr"], individual_aadt_interval["upr"]),
    units = c("(vehicles/day)^0.2", "vehicles/day"))
  write.csv(predictions, "coursework_results/prediction.csv", row.names = FALSE)
  print(predictions, row.names = FALSE)
  save_plot("03_final_transformed_scatterplots.png", function() {
    pairs(data.frame(fifth_root_aadt = d$aadt^0.2, log_population = log(d$population),
                     log_lanes = log(d$lanes), width = d$width, access = d$access),
          pch = 19, cex = 0.6, col = adjustcolor("#176b87", 0.55),
          main = "Transformed response and predictors (natural logarithms)")
  }, 1400, 1400)
  section("SESSION INFORMATION")
  print(sessionInfo())
  cat("\nRead analysis_notes.txt for the explanation of every candidate model.\n")
}

run_analysis()
