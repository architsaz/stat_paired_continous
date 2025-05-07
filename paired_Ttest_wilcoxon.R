# Load the CSV file
data <- read.csv("./table/aneu/continous_fields.csv")
print(names(data))

fields <- c("von", "eval_max", "eval_ratio")
stat_metrics <- c("mean", "mean_UQ", "max")
# Remove duplicate "area_ratio"
other_fields <- c("concen_von", "concen_eval_max", "force_ratio", "area_ratio", "area_udir")

# Define conditions
conditions <- c("homo", "hete")

# Create an empty list to hold field name pairs
pairs <- list()

# Create field pairs from main fields and stat metrics
for (field in fields) {
  for (metric in stat_metrics) {
    field_base <- paste0(metric, "_", field)
    pairs[[length(pairs) + 1]] <- c(paste0(field_base, "_homo"), paste0(field_base, "_hete"))
  }
}

# Add other field pairs
for (field in other_fields) {
  pairs[[length(pairs) + 1]] <- c(paste0(field, "_homo"), paste0(field, "_hete"))
}

# Prepare results storage
results <- data.frame(
  test_var = character(),
  method = character(),
  statistic = numeric(),
  p_value_ttest = numeric(),
  p_value_logis = numeric(),
  odds_ratio = numeric(),
  lower_ci = numeric(),
  upper_ci = numeric(),
  stringsAsFactors = FALSE
)

# Process each pair
for (pair in pairs) {
  cat("Working on:", pair[1], "vs", pair[2], "\n")

  # Skip if missing
  if (!(pair[1] %in% names(data)) || !(pair[2] %in% names(data))) {
    cat("Skipping missing column:", pair[1], "or", pair[2], "\n")
    next
  }

  homo_data <- data[[pair[1]]]
  hete_data <- data[[pair[2]]]
  diff_data <- hete_data - homo_data

  norm_homo <- shapiro.test(homo_data)$p.value > 0.05 
  norm_hete <- shapiro.test(hete_data)$p.value > 0.05 
  norm_diff <- shapiro.test(diff_data)$p.value > 0.05 

  if (norm_homo && norm_hete && norm_diff) {
    test <- t.test(homo_data, hete_data, paired = TRUE)
  } else {
    test <- wilcox.test(homo_data, hete_data, paired = TRUE)
  }

  # Logistic regression
  combined_stress <- c(homo_data, hete_data)
  condition <- c(rep(0, length(homo_data)), rep(1, length(hete_data)))  # 0 = homo, 1 = hete

  # Scale stress input
  stress_scaled <- scale(combined_stress)

  df <- data.frame(
    stress = stress_scaled,
    condition = as.factor(condition)
  )

  # Fit logistic model with control
  model <- glm(condition ~ stress, data = df, family = binomial, control = glm.control(maxit = 100))

  # Extract statistics
  odds_ratio <- exp(coef(model)["stress"])
  se <- summary(model)$coefficients["stress", "Std. Error"]
  lower_ci <- exp(coef(model)["stress"] - 1.96 * se)
  upper_ci <- exp(coef(model)["stress"] + 1.96 * se)

  # Append to results
  results <- rbind(
    results,
    data.frame(
      test_var = sub("_homo$", "", pair[1]),
      method = test$method,
      statistic = test$statistic,
      p_value_ttest = test$p.value,
      p_value_logis = summary(model)$coefficients["stress", "Pr(>|z|)"],
      odds_ratio = odds_ratio,
      lower_ci = lower_ci,
      upper_ci = upper_ci,
      stringsAsFactors = FALSE
    )
  )
}

# Save results
write.csv(results, "test_results_scaled.csv", row.names = FALSE)
