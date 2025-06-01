library(caper)
library(lmtest)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)
# Add your data as dataframes with columns 'Species' and 'dN/dS'

# Load or create your phylogeny with matching tip labels
myPhylogeny <- read.tree("bat_tags_noen.nwk")

# Data object creation
comp_dat <- comparative.data(
  phy = myPhylogeny,
  data = as.data.frame(wide_data),  # wide_data must be a data.frame (not just tibble). this is your dataframe
  names.col = "Species",           # which column has species names?
  vcv = TRUE,                      # whether to store the variance–covariance matrix
  warn.dropped = TRUE             # whether to warn if some species are dropped
)

# PGLS model. Choose the genes that should be used in the model.
model_pgls <- pgls(
  formula = LQ ~ P56704.2 + Q3B726.1 + Q9BT67.1 + P15692.3 + Q9Y6J9.1 + NP_004322.1 + Q9H9B1.4,
  data = comp_dat
)

summary(model_pgls) # Info about the model fit
bptest(model_pgls) # BP test for heteroscedasticity

# PLOTS

# Add prediction column into your main dataframe. In my case it is 'wide_data_df'
wide_data_df$pred_model <- predict(model_pgls, newdata = wide_data_df)

# Plot using ggplot2
ggplot(wide_data_df, aes(x = pred_model, y = LQ)) +
  geom_point(color = "darkred") +
  geom_line(aes(y = pred_model), color = "blue", size = 1) +
  geom_text(aes(label = Species), vjust = -0.5, hjust = 0.5, size =3) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  labs(
    title = "PGLS: LQ ~ Genes with Pos. Sel.",
    x = "Predicted LQ",
    y = "LQ"
  ) +
  theme_minimal() + coord_cartesian(clip = "off")
summary(model_pgls)

# LEAVE ONE OUT CROSS VALIDATION 

# Extract species names
species_list <- wide_data_df$Species
n_species <- length(species_list)
n_predictors <- 7  # Number of predictors in the model

# Create a data frame to store the results
results_df <- data.frame(
  Species = character(n_species),
  Absolute_Error = numeric(n_species),
  Test_R2 = numeric(n_species),
  Train_R2 = numeric(n_species),
  Adjusted_Train_R2 = numeric(n_species),
  Log_Loss = numeric(n_species),
  stringsAsFactors = FALSE
)

# Run LOOCV
for (i in 1:n_species) {
  # Exclude one species for validation
  train_data <- wide_data_df[-i, ]
  test_data <- wide_data_df[i, , drop = FALSE]
  
  # Prepare the training data for PGLS
  comp_dat_train <- comparative.data(
    phy = myPhylogeny,
    data = train_data,
    names.col = "Species",
    vcv = TRUE,
    warn.dropped = TRUE
  )
  
  # Fit the PGLS model
  pgls_model <- pgls(
    formula = LQ ~ P56704.2 + Q3B726.1 + Q9BT67.1 + P15692.3 + Q9Y6J9.1 + NP_004322.1 + Q9H9B1.4,
    data = comp_dat_train
  )
  
  # Predict for the left-out species
  pred_test <- predict(pgls_model, newdata = test_data)
  
  # Predict for the training species
  pred_train <- predict(pgls_model, newdata = train_data)
  
  # Calculate MAE for the test species
  mae <- abs(test_data$LQ - pred_test)
  
  # Calculate R^2 for the test species
  actual_test <- test_data$LQ
  mean_actual_test <- mean(train_data$LQ)
  ss_total_test <- sum((actual_test - mean_actual_test)^2)
  ss_residual_test <- sum((actual_test - pred_test)^2)
  r_squared_test <- 1 - (ss_residual_test / ss_total_test)
  
  # Calculate R^2 for the training species
  actual_train <- train_data$LQ
  mean_actual_train <- mean(actual_train)
  ss_total_train <- sum((actual_train - mean_actual_train)^2)
  ss_residual_train <- sum((actual_train - pred_train)^2)
  r_squared_train <- 1 - (ss_residual_train / ss_total_train)
  
  # Calculate Adjusted R^2 for the training species
  n_train <- nrow(train_data)
  adj_r_squared_train <- 1 - ( (1 - r_squared_train) * (n_train - 1) / (n_train - n_predictors - 1) )
  
  # Calculate log-loss (assuming normal distribution)
  epsilon <- 1e-15  # small value to avoid log(0)
  log_loss <- -sum(log(pmax(epsilon, dnorm(actual_test, mean = pred_test, sd = sd(train_data$LQ)))))
  
  # Store the results
  results_df[i, ] <- c(
    as.character(test_data$Species),
    mae,
    r_squared_test,
    r_squared_train,
    adj_r_squared_train,
    log_loss
  )
}

# Display the results
library(knitr)
kable(results_df, format = "markdown")

kable(results1_df, format="markdown")
filtered_results_df <- results_df[
  !(results_df$Species %in% c("Myotis_brandtii", "Molossus_molossus")),
  c("Species", "Absolute_Error", "Adjusted_Train_R2")
]
rownames(filtered_results_df) <- NULL
kable(filtered_results_df, format="markdown")
filtered_results_df$Adjusted_Train_R2 <- as.numeric(filtered_results_df$Adjusted_Train_R2)
filtered_results_df$Absolute_Error <- as.numeric(filtered_results_df$Absolute_Error)
mean(filtered_results_df$Absolute_Error)
summary(filtered_results_df)