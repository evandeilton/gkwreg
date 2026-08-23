# =============================================================================
# Precomputed Monte Carlo results for vignettes/gkwreg-vs-betareg.Rmd
# =============================================================================
#
# CRAN asked us to cut the vignette build time (Uwe Ligges, 2026-08-23:
# "checking re-building of vignette outputs ... [11m]"). The simulation study
# is 3 scenarios x 200 replications x 4 models = 2400 model fits, which is
# essentially all of that time. We follow the third option he offered --
# "providing precomputed results for the most lengthy parts" -- so the study
# keeps its 200 replications while the vignette only loads the summaries.
#
# Run this script to regenerate inst/extdata/vignette-simulations.rds:
#
#     Rscript data-raw/vignette-simulations.R
#
# Then re-knit the vignette and update any figures quoted in its prose, since
# the numbers change with the RNG stream.
#
# The data-generating processes below are copies of the `scenario*_dgp` chunks
# of the vignette, which displays them to the reader. Keep the two in sync.
# =============================================================================

library(gkwreg)
library(betareg)

MODEL_NAMES <- c(
  betareg = "Beta (betareg)",
  gkw_beta = "Beta (gkwreg)",
  gkw_kw = "Kumaraswamy",
  gkw_ekw = "Exp. Kumaraswamy"
)

# --- simulation infrastructure (was the vignette's `helper_functions` chunk) --

# ==============================================================================
# SIMULATION INFRASTRUCTURE
# ==============================================================================

# Unified simulation function (refactored to eliminate code duplication)
run_single_sim <- function(n, dgp_fun, params, formula, test_prop = 0.2,
                           return_fits = FALSE) {
  data_full <- tryCatch(
    {
      dgp_fun(n, params)
    },
    error = function(e) NULL
  )

  if (is.null(data_full)) {
    return(list(
      betareg = list(success = FALSE),
      gkw_beta = list(success = FALSE),
      gkw_kw = list(success = FALSE),
      gkw_ekw = list(success = FALSE)
    ))
  }

  n_test <- floor(n * test_prop)
  test_idx <- sample(seq_len(n), n_test)
  train_data <- data_full[-test_idx, ]
  test_data <- data_full[test_idx, ]

  results <- list()

  safe_fit <- function(fit_fun, name) {
    tryCatch(
      {
        t0 <- proc.time()[3]
        fit <- fit_fun()
        time_elapsed <- proc.time()[3] - t0

        if (name == "betareg") {
          loglik <- as.numeric(logLik(fit))
          aic <- AIC(fit)
          bic <- BIC(fit)
          converged <- isTRUE(fit$converged)
          se <- sqrt(diag(vcov(fit)))
        } else {
          loglik <- fit$loglik
          aic <- fit$aic
          bic <- fit$bic
          converged <- isTRUE(fit$convergence)
          se <- fit$se
        }

        pred <- predict(fit, newdata = test_data, type = "response")

        if (any(is.na(pred)) || any(is.infinite(pred))) {
          return(list(success = FALSE))
        }

        rmse <- sqrt(mean((test_data$y - pred)^2))

        if (!is.finite(loglik)) loglik <- NA_real_
        if (!is.finite(aic)) aic <- NA_real_
        if (!is.finite(bic)) bic <- NA_real_
        if (!is.finite(rmse)) rmse <- NA_real_
        if (!is.finite(time_elapsed)) time_elapsed <- NA_real_

        result <- list(
          loglik = as.numeric(loglik),
          aic = as.numeric(aic),
          bic = as.numeric(bic),
          rmse = as.numeric(rmse),
          time = as.numeric(time_elapsed),
          converged = converged,
          success = TRUE
        )

        # Optionally include full fit object and coefficients
        if (return_fits) {
          result$fit <- fit
          result$coef <- coef(fit)
          result$se <- se
        }

        result
      },
      error = function(e) {
        list(success = FALSE)
      }
    )
  }

  results$betareg <- safe_fit(
    function() betareg(formula, data = train_data),
    "betareg"
  )

  results$gkw_beta <- safe_fit(
    function() gkwreg(formula, data = train_data, family = "beta"),
    "gkwreg"
  )

  results$gkw_kw <- safe_fit(
    function() gkwreg(formula, data = train_data, family = "kw"),
    "gkwreg"
  )

  results$gkw_ekw <- safe_fit(
    function() gkwreg(formula, data = train_data, family = "ekw"),
    "gkwreg"
  )

  results
}

run_full_simulation <- function(n_sim, n, dgp_fun, params, formula, verbose = FALSE) {
  all_results <- vector("list", n_sim)

  for (i in seq_len(n_sim)) {
    all_results[[i]] <- run_single_sim(n, dgp_fun, params, formula)
  }

  models <- c("betareg", "gkw_beta", "gkw_kw", "gkw_ekw")
  summary_stats <- list()

  for (model in models) {
    success_idx <- sapply(all_results, function(x) {
      !is.null(x[[model]]) && isTRUE(x[[model]]$success)
    })

    n_success <- sum(success_idx)

    if (n_success == 0) {
      summary_stats[[model]] <- list(
        n_success = 0,
        conv_rate = 0,
        loglik_mean = NA_real_,
        aic_mean = NA_real_,
        bic_mean = NA_real_,
        rmse_mean = NA_real_,
        time_mean = NA_real_
      )
      next
    }

    successful <- all_results[success_idx]

    extract_metric <- function(metric_name) {
      vals <- sapply(successful, function(x) {
        val <- x[[model]][[metric_name]]
        if (is.null(val) || !is.numeric(val)) {
          return(NA_real_)
        }
        as.numeric(val)
      })
      if (all(is.na(vals))) {
        return(NA_real_)
      }
      vals
    }

    converged_vals <- sapply(successful, function(x) {
      val <- x[[model]]$converged
      if (is.null(val)) {
        return(FALSE)
      }
      if (is.logical(val)) {
        return(val)
      }
      if (is.numeric(val)) {
        return(val > 0)
      }
      FALSE
    })

    summary_stats[[model]] <- list(
      n_success = n_success,
      conv_rate = mean(converged_vals, na.rm = TRUE),
      loglik_mean = mean(extract_metric("loglik"), na.rm = TRUE),
      aic_mean = mean(extract_metric("aic"), na.rm = TRUE),
      bic_mean = mean(extract_metric("bic"), na.rm = TRUE),
      rmse_mean = mean(extract_metric("rmse"), na.rm = TRUE),
      time_mean = mean(extract_metric("time"), na.rm = TRUE)
    )
  }

  list(
    raw = all_results,
    summary = summary_stats,
    n = n,
    n_sim = n_sim,
    formula = formula
  )
}

make_comparison_table <- function(sim_results) {
  models <- names(sim_results$summary)

  safe_format <- function(x, digits = 2) {
    if (is.null(x) || is.na(x)) {
      return(NA_real_)
    }
    round(as.numeric(x), digits)
  }

  df_list <- lapply(models, function(model) {
    s <- sim_results$summary[[model]]

    if (is.null(s) || s$n_success == 0) {
      return(data.frame(
        Model = MODEL_NAMES[model],
        N_Success = 0,
        Conv_Rate = NA_real_,
        AIC = NA_real_,
        RMSE = NA_real_,
        Time = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      Model = MODEL_NAMES[model],
      N_Success = s$n_success,
      Conv_Rate = safe_format(s$conv_rate * 100, 1),
      AIC = safe_format(s$aic_mean, 2),
      RMSE = safe_format(s$rmse_mean, 4),
      Time = safe_format(s$time_mean, 4),
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, df_list)
  rownames(df) <- NULL
  df
}

# --- data-generating processes (mirror the vignette's `scenario*_dgp` chunks) -

dgp_beta <- function(n, params) {
  x1 <- rnorm(n, 0, 1)
  x2 <- runif(n, -1, 1)

  eta_mu <- params$beta_mu[1] + params$beta_mu[2] * x1 + params$beta_mu[3] * x2
  eta_phi <- params$beta_phi[1] + params$beta_phi[2] * x1

  mu <- plogis(eta_mu)
  phi <- exp(eta_phi)

  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  y <- pmax(pmin(y, 0.999), 0.001)

  data.frame(y = y, x1 = x1, x2 = x2, mu = mu, phi = phi)
}

dgp_heavy_tails <- function(n, params) {
  x1 <- rnorm(n, 0, 1)
  x2 <- rbinom(n, 1, 0.5)

  eta_alpha <- params$beta_alpha[1] + params$beta_alpha[2] * x1
  eta_beta <- params$beta_beta[1] + params$beta_beta[2] * x2
  eta_lambda <- params$beta_lambda[1]

  alpha <- exp(eta_alpha)
  beta <- exp(eta_beta)
  lambda <- exp(eta_lambda)

  u <- runif(n)
  y <- (1 - (1 - u^(1 / beta))^(1 / alpha))^(1 / lambda)
  y <- pmax(pmin(y, 0.9999), 0.0001)

  data.frame(y = y, x1 = x1, x2 = factor(x2))
}

dgp_extreme <- function(n, params) {
  x1 <- rnorm(n, 0, 1)
  group <- sample(c("J", "U"), n, replace = TRUE)

  alpha <- ifelse(
    group == "J",
    exp(params$alpha_J[1] + params$alpha_J[2] * x1),
    exp(params$alpha_U[1] + params$alpha_U[2] * x1)
  )

  beta <- ifelse(group == "J", exp(params$beta_J), exp(params$beta_U))

  u <- runif(n)
  y <- (1 - (1 - u)^(1 / beta))^(1 / alpha)
  y <- pmax(pmin(y, 0.9999), 0.0001)

  data.frame(y = y, x1 = x1, group = factor(group))
}

# --- warm up the TMB models -------------------------------------------------
#
# gkwreg compiles each TMB family to a shared library on first use, which costs
# roughly 30 seconds per family. `safe_fit()` times the fitting call, so without
# this warm-up that compilation lands inside replicate 1 and is then averaged
# over all replicates: at 200 replicates it inflates every reported gkwreg
# timing by ~0.15s, which is an order of magnitude more than a fit actually
# takes. Compiling here, outside the timed loops, makes the reported times
# measure estimation rather than compilation.

message("Warming up TMB models (compilation is not part of the timings) ...")
local({
  set.seed(1)
  n <- 200
  x1 <- rnorm(n)
  x2 <- runif(n, -1, 1)
  mu <- plogis(0.5 - 0.8 * x1 + 0.6 * x2)
  phi <- exp(1.5 + 0.4 * x1)
  y <- pmax(pmin(rbeta(n, mu * phi, (1 - mu) * phi), 0.999), 0.001)
  warm <- data.frame(y = y, x1 = x1, x2 = x2)
  for (fam in c("beta", "kw", "ekw")) {
    invisible(gkwreg(y ~ x1 + x2 | x1, data = warm, family = fam))
  }
})

# --- run the study ----------------------------------------------------------

N_SIM <- 200

message("Scenario 1: well-specified Beta ...")
set.seed(20260823)
results_s1 <- run_full_simulation(
  n_sim = N_SIM,
  n = 300,
  dgp_fun = dgp_beta,
  params = list(
    beta_mu = c(0.5, -0.8, 0.6),
    beta_phi = c(1.5, 0.4)
  ),
  formula = y ~ x1 + x2 | x1
)

message("Scenario 2: heavy tails ...")
set.seed(20260824)
results_s2 <- run_full_simulation(
  n_sim = N_SIM,
  n = 300,
  dgp_fun = dgp_heavy_tails,
  params = list(
    beta_alpha = c(0.8, -0.5),
    beta_beta = c(0.3, 0.4),
    beta_lambda = c(0.6)
  ),
  formula = y ~ x1 | x2
)

message("Scenario 3: extreme boundary concentration ...")
set.seed(20260825)
results_s3 <- run_full_simulation(
  n_sim = N_SIM,
  n = 400,
  dgp_fun = dgp_extreme,
  params = list(
    alpha_J = c(-1.5, 0.2),
    beta_J = 2.0,
    alpha_U = c(-1.8, 0.1),
    beta_U = -0.8
  ),
  formula = y ~ x1 * group | group
)

# Store only the summary tables. The raw per-replication fits are large and the
# vignette never touches them.
vignette_simulations <- list(
  comp_s1 = make_comparison_table(results_s1),
  comp_s2 = make_comparison_table(results_s2),
  comp_s3 = make_comparison_table(results_s3),
  n_sim = N_SIM,
  n_obs = c(s1 = 300, s2 = 300, s3 = 400),
  generated_on = Sys.Date(),
  r_version = R.version.string,
  gkwreg_version = as.character(utils::packageVersion("gkwreg")),
  betareg_version = as.character(utils::packageVersion("betareg"))
)

saveRDS(
  vignette_simulations,
  file.path("inst", "extdata", "vignette-simulations.rds"),
  compress = "xz"
)

message("Wrote inst/extdata/vignette-simulations.rds")
print(vignette_simulations$comp_s1)
print(vignette_simulations$comp_s2)
print(vignette_simulations$comp_s3)
