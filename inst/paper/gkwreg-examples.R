## -------------------------------------------------------------------------- ##
## gkwreg: Generalized Kumaraswamy Regression Models for Bounded Data
## Author: Lopes, J. E. // 2025
## -------------------------------------------------------------------------- ##

if(!require(gkwreg)) install.packages("gkwreg")
if(!require(skimr)) install.packages("skimr")

## -------------------------------------------------------------------------- ##
## Extra functions for analysis
## -------------------------------------------------------------------------- ##

# Wrapper function to compare different gkw models for bounded data
fit_and_compare_models_2p <- function(
    formula, data, method = "BFGS",
    families = c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"),
    use_betareg = TRUE) {
  
  # Map family codes to descriptive names
  family_names <- c(
    gkw = "Generalized Kumaraswamy", bkw = "Beta-Kumaraswamy",
    kkw = "Kumaraswamy-Kumaraswamy", ekw = "Exponential Kumaraswamy",
    mc = "McDonald/Beta Power", kw = "Kumaraswamy", beta = "Beta (gkwreg)"
  )
  
  # Number of parameters for each family
  family_params <- list(
    gkw = 5, bkw = 4, kkw = 4, ekw = 3, mc = 3, kw = 2, beta = 2
  )
  
  results <- list()
  
  # Fit gkwreg models
  for (fam in families) {
    num_links <- family_params[[fam]]
    
    # Try to fit model with basic error handling
    fit <- try(gkwreg(formula, data = data, family = fam,
                      link = rep("log", num_links),
                      method = method, silent = TRUE), silent = TRUE)
    
    # Extract model statistics if fit successful
    if (!inherits(fit, "try-error")) {
      results[[fam]] <- data.frame(
        Model = family_names[[fam]],
        Code = fam,
        Parameters = length(unlist(coef(fit))),
        logLik = as.numeric(logLik(fit)),
        AIC = AIC(fit),
        BIC = BIC(fit),
        Conv = paste0(as.numeric(fit$convergence), collapse = ";"),
        row.names = fam
      )
    }
  }
  
  # Fit betareg model if requested
  if (use_betareg) {
    
    # Fit standard beta regression
    br_fit <- try(
      betareg::betareg(formula, data = data,
                       link = "log", link.phi = "log"),
      silent = TRUE)
    
    # Add betareg results if successful
    if (!inherits(br_fit, "try-error")) {
      results[["betareg"]] <- data.frame(
        Model = "Beta (betareg)",
        Code = "betareg",
        Parameters = length(coef(br_fit)),
        logLik = as.numeric(logLik(br_fit)),
        AIC = AIC(br_fit),
        BIC = BIC(br_fit),
        # Conv = max(as.numeric(br_fit$converged)),
        Conv = paste0(as.numeric(br_fit$converged), collapse = ";"),
        row.names = NULL
      )
    }
  }
  
  # Combine all results and sort by AIC
  results_df <- do.call(rbind, results)
  results_df <- results_df[order(results_df$AIC), ]
  rownames(results_df) <- NULL
  
  return(results_df)
}

# Function to fit and compare models with 3+ parameters
fit_and_compare_models_3p <- function(formula, data, method = "nlminb",
                               families = c("gkw", "bkw", "kkw", "ekw", "mc")) {
  # Map family codes to descriptive names  
  family_names <- c(
    gkw = "Generalized Kumaraswamy", bkw = "Beta-Kumaraswamy",
    kkw = "Kumaraswamy-Kumaraswamy", ekw = "Exponential Kumaraswamy",
    mc = "McDonald/Beta Power"
  )
  
  # Define number of parameters for each family
  family_params <- list(gkw = 5, bkw = 4, kkw = 4, ekw = 3, mc = 3)
  
  results <- list()
  
  # Fit each model family
  for (fam in families) {
    # Try to fit the model
    fit <- try(gkwreg(formula, 
                      data = data, 
                      family = fam,
                      link = rep("log", family_params[[fam]]),
                      method = method, 
                      silent = TRUE), 
               silent = TRUE)
    
    # Extract statistics if fit successful
    if (!inherits(fit, "try-error")) {
      results[[fam]] <- data.frame(
        Model = family_names[[fam]],
        Family = fam,
        Parameters = length(unlist(coef(fit))),
        logLik = as.numeric(logLik(fit)),
        AIC = AIC(fit),
        BIC = BIC(fit),
        Conv = max(as.numeric(fit$convergence)),
        row.names = NULL
      )
    }
  }
  
  # Combine results and sort by AIC
  results_df <- do.call(rbind, results)
  results_df <- results_df[order(results_df$AIC), ]
  rownames(results_df) <- NULL
  
  return(results_df)
}


## -------------------------------------------------------------------------- ##
## Ex-01 - FoodExpenditure data
## -------------------------------------------------------------------------- ##

# Get FoodExpenditure data and create 'y' as the response
food_data <- get_bounded_datasets("FoodExpenditure")
food_data <- within(food_data, {y = food / income})
skimr::skim(food_data)

# Define the formula: y depends on 'persons'
# We'll model only gamma and delta for Beta, keeping other parameters constant
formu_fe <- y ~ persons | income

# Fit all families
results_food_data <- fit_and_compare_models_2p(formu_fe, food_data, method = "nlminb")

# Best model
kw_model <- gkwreg(formu_fe, food_data, family = "kw",
                   link = rep("log", 2), method = "nlminb")

summary(kw_model)
plot(kw_model, use_ggplot = TRUE, arrange_plots = TRUE, sub.caption = "")

## -------------------------------------------------------------------------- ##
## Ex-02 - GasolineYield data
## -------------------------------------------------------------------------- ##

# Load GasolineYield data
gasoline_data <- get_bounded_datasets("GasolineYield")

# Formula: yield ~ batch + temp | temp
# First part (for alpha/gamma) includes batch and temp
# Second part (for beta/delta/phi) includes only temp
formu_gy <- yield ~ batch + temp | temp

# Load the GasolineYield data
gasoline_data <- get_bounded_datasets("GasolineYield")
skimr::skim(gasoline_data)

# Define the formula: yield ~ batch + temp | temp
# The first part (for alpha/gamma) includes batch and temp
# The second part (for beta/delta/phi) includes only temp
formu_gy <- yield ~ batch + temp | temp

# Execute the comparative analysis
results_gasoline <- fit_and_compare_models_2p(formu_gy, gasoline_data, method = "BFGS")

# Best model
kw_model_gas <- gkwreg(formu_gy, gasoline_data, family = "kw",
                       link = rep("log", 2), method = "BFGS")
summary(kw_model_gas)
plot(kw_model_gas, use_ggplot = TRUE, arrange_plots = TRUE, sub.caption = "")

## -------------------------------------------------------------------------- ##
## Ex-03 - sdac data
## -------------------------------------------------------------------------- ##

# Load the sdac data
sdac_data <- get_bounded_datasets("sdac")
skimr::skim(sdac_data)

# Formula: rcd ~ ageadj | chemo
formu_sd <- rcd ~ ageadj + chemo

# Compare different families in gkwreg
results_sdac <- fit_and_compare_models_2p(formu_sd, sdac_data,
                                          method = "nlminb", use_betareg = TRUE)

# Best model
ekw_model_gas <- gkwreg(formu_sd, sdac_data, family = "ekw", method = "BFGS")
summary(ekw_model_gas)
plot(ekw_model_gas, use_ggplot = TRUE, arrange_plots = TRUE, sub.caption = "")

## -------------------------------------------------------------------------- ##
## Ex-04 - retinal data
## -------------------------------------------------------------------------- ##

# Load the retinal data
retinal_data <- get_bounded_datasets("retinal")
skimr::skim(retinal_data)

# Formula modeling alpha, beta, gamma (for EKw/Mc) or more (KKw/BKw/GKw)
# alpha ~ LogT + LogT2 + Level
# beta  ~ LogT + Level
# gamma ~ Time (or constant if Kw/Beta)
formu_rt <- Gas ~ LogT + LogT2 + Level | LogT + Level | Time

# Fit model
results_retinal <- fit_and_compare_models_3p(formu_rt, retinal_data, method = "nlminb")

# Best model
ekw_model_ret <- gkwreg(formu_rt, retinal_data, family = "ekw", method = "nlminb")
summary(ekw_model_ret)
plot(ekw_model_ret, use_ggplot = TRUE, arrange_plots = TRUE, sub.caption = "")


## -------------------------------------------------------------------------- ##
## Ex-05 - WeatherTask Data
## -------------------------------------------------------------------------- ##

# Geta the "WeatherTask" data
df_weather <- get_bounded_datasets("WeatherTask")
skimr::skim(df_weather)

# Fit all seven distribution families to the 'agreement' data
# We'll use gkwfitall, which is a convenient wrapper for this purpose
fitall_weather <- gkwfitall(df_weather$agreement, method = "BFGS")

# Analyze the comparative results
summary(fitall_weather) # Displays the comparison table

# Obtain the name of the best family (based on the lowest AIC, for example)
# The table is already sorted by AIC
best_family_code <- fitall_weather$comparison$Family[1]

# Refit the best model to obtain additional details and plots
fit_best_weather <- gkwfit(
  df_weather$agreement, family = best_family_code,
  method = "BFGS", profile = TRUE, plot = TRUE, silent = TRUE)
# set plot = FALSE here if you want to generate later

# Generate the Goodness-of-Fit (GoF) report for the best model
gof_report <- gkwgof(
  fit_best_weather, theme = ggplot2::theme_classic(),
  plot = TRUE, print_summary = FALSE, verbose = FALSE)
# summary(gof_report) # Displays GoF statistics

# Extract GoF and fitting statistics for all the families
# fitted by gkwfitall
results_weathertask_df <- do.call(rbind, lapply(fitall_weather$fits, function(f){
  extract_gof_stats(gkwgof(f, plot = FALSE,
                           print_summary = FALSE, verbose = FALSE))
}))
results_weathertask_df <- results_weathertask_df[order(results_weathertask_df$AIC), ]
row.names(results_weathertask_df) <- NULL

# Generate the fit plot for the best model and save it as a PDF
# (The plot generated by gkwgof includes a histogram, empirical and 
# fitted density curves, P-P, and Q-Q plots)
plot(gkwgof(fit_best_weather, theme = ggplot2::theme_classic()), title = "")

# Display the formatted table with comparative results
# Select and rename columns for display
results_weathertask_display <- results_weathertask_df[,
                                c("family", "n_params", "logLik", "AIC", "BIC",
                                  "KS", "AD", "RMSE", "pseudo_R2")]
