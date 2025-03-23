#=====================================================================
# Estudo Avançado da Distribuição Generalizada Kumaraswamy (GKw)
# Análise de seus casos especiais e propriedades
#=====================================================================

# Carregar os pacotes necessários
library(gkwreg)    # Contém: dgkw, pgkw, qgkw, rgkw, llgkw, grgkw, hsgkw, nrgkw
library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape2)
library(latex2exp)  # Para equações LaTeX nos gráficos

#-----------------------------------------------------------
# 1. Implementação dos Casos Especiais da Distribuição GKw
#-----------------------------------------------------------

# Mapeamento de parâmetros para casos especiais conforme o artigo
# Função para verificar se um conjunto de parâmetros pertence a um caso especial
identify_special_case <- function(alpha, beta, gamma, delta, lambda) {
  # Tolerância para comparação de números decimais
  tol <- 1e-6

  if (abs(lambda - 1) < tol && abs(gamma - 1) < tol && abs(delta) < tol) {
    return("Kumaraswamy (Kw)")
  } else if (abs(alpha - 1) < tol && abs(beta - 1) < tol) {
    return("McDonald (Mc)")
  } else if (abs(alpha - 1) < tol && abs(beta - 1) < tol && abs(lambda - 1) < tol) {
    return("Beta")
  } else if (abs(lambda - 1) < tol) {
    return("Beta Kumaraswamy (BKw)")
  } else if (abs(gamma - 1) < tol) {
    return("Kumaraswamy-Kumaraswamy (KwKw)")
  } else if (abs(delta) < tol && abs(gamma - 1) < tol) {
    return("Kumaraswamy Exponenciada (EKw)")
  } else if (abs(alpha - 1) < tol && abs(beta - 1) < tol) {
    return("Beta Power (BP)")
  } else {
    return("Generalizada Kumaraswamy (GKw)")
  }
}

# Funções específicas para cada caso especial
# Casos especiais da distribuição GKw conforme artigo (páginas 4-8):
# 1. Distribuição Kumaraswamy (Kw): λ = γ = 1, δ = 0
# 2. Distribuição McDonald (Mc): α = β = 1
# 3. Distribuição Beta: α = β = λ = 1
# 4. Distribuição Beta Kumaraswamy (BKw): λ = 1
# 5. Distribuição Kumaraswamy-Kumaraswamy (KwKw): γ = 1
# 6. Distribuição Kumaraswamy Exponenciada (EKw): δ = 0, γ = 1
# 7. Distribuição Beta Power (BP): α = 1, β = 1

# Vamos testar a generalização da função dgkw() para esses casos especiais

#-----------------------------------------------------------
# 2. Função Aprimorada para ajuste do modelo via MLE
#-----------------------------------------------------------
fit_gkw_model <- function(data, start_params = NULL, verbose = FALSE, model = "GKw") {
  # Definir parâmetros iniciais baseados no modelo escolhido
  if (is.null(start_params)) {
    if (model == "GKw") {
      start_params <- c(alpha = 2, beta = 2, gamma = 2, delta = 1, lambda = 2)
    } else if (model == "Kw") {
      start_params <- c(alpha = 2, beta = 2, gamma = 1, delta = 0, lambda = 1)
    } else if (model == "BKw") {
      start_params <- c(alpha = 2, beta = 2, gamma = 2, delta = 1, lambda = 1)
    } else if (model == "KwKw") {
      start_params <- c(alpha = 2, beta = 2, gamma = 1, delta = 1, lambda = 2)
    } else if (model == "EKw") {
      start_params <- c(alpha = 2, beta = 2, gamma = 1, delta = 0, lambda = 2)
    } else if (model == "Mc") {
      start_params <- c(alpha = 1, beta = 1, gamma = 2, delta = 1, lambda = 2)
    } else if (model == "Beta") {
      start_params <- c(alpha = 1, beta = 1, gamma = 2, delta = 1, lambda = 1)
    } else if (model == "BP") {
      start_params <- c(alpha = 1, beta = 1, gamma = 2, delta = 1, lambda = 2)
    }
  }

  # Aplicar restrições conforme o modelo - modificando os valores iniciais
  if (model == "Kw") {
    start_params["gamma"] <- 1
    start_params["delta"] <- 0
    start_params["lambda"] <- 1
  } else if (model == "BKw") {
    start_params["lambda"] <- 1
  } else if (model == "KwKw") {
    start_params["gamma"] <- 1
  } else if (model == "EKw") {
    start_params["gamma"] <- 1
    start_params["delta"] <- 0
  } else if (model == "Mc") {
    start_params["alpha"] <- 1
    start_params["beta"] <- 1
  } else if (model == "Beta") {
    start_params["alpha"] <- 1
    start_params["beta"] <- 1
    start_params["lambda"] <- 1
  } else if (model == "BP") {
    start_params["alpha"] <- 1
    start_params["beta"] <- 1
  }

  # Função para aplicar restrições após cada iteração do algoritmo
  # Esta abordagem não pode ser usada pois nrgkw não tem o argumento constrain
  # Em vez disso, vamos usar uma abordagem de estimação em duas etapas:

  # Etapa 1: Ajuste com os parâmetros iniciais adequados
  fit <- try(nrgkw(start_params, data, verbose = verbose), silent = TRUE)

  # Etapa 2: Se o modelo não for GKw, aplicamos as restrições nos parâmetros estimados
  if (!inherits(fit, "try-error") && is.list(fit) && model != "GKw") {
    estimated_params <- fit$parameters

    if (model == "Kw") {
      estimated_params["gamma"] <- 1
      estimated_params["delta"] <- 0
      estimated_params["lambda"] <- 1
    } else if (model == "BKw") {
      estimated_params["lambda"] <- 1
    } else if (model == "KwKw") {
      estimated_params["gamma"] <- 1
    } else if (model == "EKw") {
      estimated_params["gamma"] <- 1
      estimated_params["delta"] <- 0
    } else if (model == "Mc") {
      estimated_params["alpha"] <- 1
      estimated_params["beta"] <- 1
    } else if (model == "Beta") {
      estimated_params["alpha"] <- 1
      estimated_params["beta"] <- 1
      estimated_params["lambda"] <- 1
    } else if (model == "BP") {
      estimated_params["alpha"] <- 1
      estimated_params["beta"] <- 1
    }

    # Recalcular log-verossimilhança com os parâmetros restritos
    loglik <- llgkw(estimated_params, data)

    # Atualizar resultados
    fit$parameters <- estimated_params
    fit$loglik <- loglik

    # Calcular critérios de informação
    n <- length(data)
    k <- switch(model,
                "GKw" = 5,
                "Kw" = 2,
                "BKw" = 4,
                "KwKw" = 4,
                "EKw" = 3,
                "Mc" = 3,
                "Beta" = 2,
                "BP" = 3)

    fit$AIC <- -2 * loglik + 2 * k
    fit$BIC <- -2 * loglik + log(n) * k
  } else if (!inherits(fit, "try-error") && is.list(fit)) {
    # Calcular AIC e BIC para o modelo GKw
    n <- length(data)
    fit$AIC <- -2 * fit$loglik + 2 * 5  # GKw tem 5 parâmetros
    fit$BIC <- -2 * fit$loglik + log(n) * 5
  }

  # Adicionar informações extras ao resultado
  if (!inherits(fit, "try-error") && is.list(fit)) {
    fit$model <- model
  }

  return(fit)
}

#-----------------------------------------------------------
# 3. Função Melhorada para Gráficos Diagnósticos
#-----------------------------------------------------------
diagnostic_plots <- function(data, fits, plot_quantiles = TRUE) {
  # fits deve ser uma lista de ajustes de diferentes modelos

  plot_list <- list()

  # Histograma com curvas de densidade ajustadas sobrepostas
  hist_data <- data.frame(x = data)
  x_grid <- seq(min(data), max(data), length.out = 200)

  dens_df <- data.frame(x = x_grid)

  for (i in 1:length(fits)) {
    fit <- fits[[i]]
    model <- fit$model
    params <- fit$parameters

    # Calcular densidade para o grid usando sempre dgkw() com 5 parâmetros
    # conforme as restrições de cada caso especial
    dens <- dgkw(x_grid,
                 alpha = params["alpha"],
                 beta = params["beta"],
                 gamma = params["gamma"],
                 delta = params["delta"],
                 lambda = params["lambda"])

    dens_df[[model]] <- dens
  }

  # Preparar dados para ggplot
  dens_df_long <- reshape2::melt(dens_df, id.vars = "x", variable.name = "Model", value.name = "Density")

  # Criar paleta de cores
  n_models <- length(fits)
  colors <- rainbow(n_models)

  # Histograma e curvas de densidade
  p1 <- ggplot(hist_data, aes(x = x)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightgray", color = "black", alpha = 0.6) +
    geom_line(data = dens_df_long, aes(x = x, y = Density, color = Model, linetype = Model), size = 1) +
    scale_color_manual(values = colors) +
    labs(title = "Histograma e Curvas de Densidade Ajustadas",
         x = "Valores", y = "Densidade") +
    theme_minimal() +
    theme(legend.position = "bottom")

  plot_list[[1]] <- p1

  # QQ-Plot para cada modelo
  if (plot_quantiles) {
    for (i in 1:length(fits)) {
      fit <- fits[[i]]
      model <- fit$model
      params <- fit$parameters

      n <- length(data)
      p_seq <- (1:n) / (n + 1)

      # Calcular quantis teóricos usando sempre qgkw() com 5 parâmetros
      theo_q <- qgkw(p_seq,
                     alpha = params["alpha"],
                     beta = params["beta"],
                     gamma = params["gamma"],
                     delta = params["delta"],
                     lambda = params["lambda"])

      df_qq <- data.frame(empirico = sort(data), teorico = theo_q)

      p_qq <- ggplot(df_qq) +
        geom_point(aes(x = teorico, y = empirico), color = colors[i]) +
        geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
        labs(title = paste("QQ-Plot:", model),
             x = "Quantis Teóricos", y = "Quantis Empíricos") +
        theme_minimal()

      plot_list[[i+1]] <- p_qq
    }
  }

  # Tabela comparativa de AIC e BIC
  model_names <- sapply(fits, function(fit) fit$model)
  loglik_values <- sapply(fits, function(fit) fit$loglik)
  aic_values <- sapply(fits, function(fit) fit$AIC)
  bic_values <- sapply(fits, function(fit) fit$BIC)

  comparison_table <- data.frame(
    Model = model_names,
    LogLik = round(loglik_values, 2),
    AIC = round(aic_values, 2),
    BIC = round(bic_values, 2)
  )

  # Ordenar por AIC
  comparison_table <- comparison_table[order(comparison_table$AIC), ]

  # Imprimir a tabela
  cat("\nComparação de Modelos:\n")
  print(comparison_table, row.names = FALSE)

  # Arranjar gráficos em grid
  if (length(plot_list) <= 2) {
    grid.arrange(grobs = plot_list, ncol = 1)
  } else if (length(plot_list) <= 5) {
    grid.arrange(grobs = plot_list, ncol = 2)
  } else {
    grid.arrange(grobs = plot_list[1:6], ncol = 2)
    if (length(plot_list) > 6) {
      grid.arrange(grobs = plot_list[7:min(12, length(plot_list))], ncol = 2)
    }
  }

  return(invisible(comparison_table))
}

#-----------------------------------------------------------
# 4. Função para Visualizar as Formas da Densidade GKw e Casos Especiais
#-----------------------------------------------------------
visualize_gkw_shapes <- function() {
  # Grid de valores para x
  x_grid <- seq(0.001, 0.999, length.out = 200)

  # Exemplo 1: Variação de alpha
  alphas <- c(0.5, 2.0, 5.0, 15)
  dens_alpha <- sapply(alphas, function(a) {
    dgkw(x_grid, alpha = a, beta = 3.5, gamma = 1.5, delta = 2.5, lambda = 0.5)
  })
  df_alpha <- data.frame(x = x_grid, dens_alpha)
  names(df_alpha) <- c("x", paste("alpha =", alphas))

  # Exemplo 2: Variação de beta
  betas <- c(0.1, 0.32, 1.0, 6.0)
  dens_beta <- sapply(betas, function(b) {
    dgkw(x_grid, alpha = 3.5, beta = b, gamma = 1.5, delta = 2.5, lambda = 0.5)
  })
  df_beta <- data.frame(x = x_grid, dens_beta)
  names(df_beta) <- c("x", paste("beta =", betas))

  # Exemplo 3: Variação de gamma
  gammas <- c(0.5, 3.5, 15, 30)
  dens_gamma <- sapply(gammas, function(g) {
    dgkw(x_grid, alpha = 1.0, beta = 1.5, gamma = g, delta = 2.5, lambda = 0.5)
  })
  df_gamma <- data.frame(x = x_grid, dens_gamma)
  names(df_gamma) <- c("x", paste("gamma =", gammas))

  # Exemplo 4: Variação de delta
  deltas <- c(0.1, 0.5, 1.0, 1.5)
  dens_delta <- sapply(deltas, function(d) {
    dgkw(x_grid, alpha = 1.0, beta = 1.5, gamma = 2.5, delta = d, lambda = 0.5)
  })
  df_delta <- data.frame(x = x_grid, dens_delta)
  names(df_delta) <- c("x", paste("delta =", deltas))

  # Exemplo 5: Variação de lambda
  lambdas <- c(0.5, 1.0, 5.0, 10.0)
  dens_lambda <- sapply(lambdas, function(l) {
    dgkw(x_grid, alpha = 0.5, beta = 0.7, gamma = 0.1, delta = 3.0, lambda = l)
  })
  df_lambda <- data.frame(x = x_grid, dens_lambda)
  names(df_lambda) <- c("x", paste("lambda =", lambdas))

  # Preparar dados para ggplot
  df_alpha_long <- reshape2::melt(df_alpha, id.vars = "x", variable.name = "alpha", value.name = "Densidade")
  df_beta_long <- reshape2::melt(df_beta, id.vars = "x", variable.name = "beta", value.name = "Densidade")
  df_gamma_long <- reshape2::melt(df_gamma, id.vars = "x", variable.name = "gamma", value.name = "Densidade")
  df_delta_long <- reshape2::melt(df_delta, id.vars = "x", variable.name = "delta", value.name = "Densidade")
  df_lambda_long <- reshape2::melt(df_lambda, id.vars = "x", variable.name = "lambda", value.name = "Densidade")

  # Criar gráficos
  p1 <- ggplot(df_alpha_long, aes(x = x, y = Densidade, color = alpha)) +
    geom_line() +
    labs(title = "Variação de α") +
    theme_minimal()

  p2 <- ggplot(df_beta_long, aes(x = x, y = Densidade, color = beta)) +
    geom_line() +
    labs(title = "Variação de β") +
    theme_minimal()

  p3 <- ggplot(df_gamma_long, aes(x = x, y = Densidade, color = gamma)) +
    geom_line() +
    labs(title = "Variação de γ") +
    theme_minimal()

  p4 <- ggplot(df_delta_long, aes(x = x, y = Densidade, color = delta)) +
    geom_line() +
    labs(title = "Variação de δ") +
    theme_minimal()

  p5 <- ggplot(df_lambda_long, aes(x = x, y = Densidade, color = lambda)) +
    geom_line() +
    labs(title = "Variação de λ") +
    theme_minimal()

  # Casos especiais interessantes
  # Exemplo de J-shape
  j_shape <- dgkw(x_grid, alpha = 0.5, beta = 0.5, gamma = 0.5, delta = 0.1, lambda = 0.5)
  # Exemplo de U-shape
  u_shape <- dgkw(x_grid, alpha = 0.5, beta = 0.5, gamma = 0.1, delta = 0.1, lambda = 0.5)
  # Exemplo de forma de sino
  bell_shape <- dgkw(x_grid, alpha = 2.0, beta = 2.0, gamma = 2.0, delta = 0.5, lambda = 1.0)
  # Exemplo de forma exponencial
  exp_shape <- dgkw(x_grid, alpha = 1.0, beta = 1.0, gamma = 0.5, delta = 0.5, lambda = 0.5)

  df_special <- data.frame(
    x = x_grid,
    "J-shape" = j_shape,
    "U-shape" = u_shape,
    "Bell shape" = bell_shape,
    "Exponential" = exp_shape
  )

  df_special_long <- reshape2::melt(df_special, id.vars = "x", variable.name = "shape", value.name = "Densidade")

  p6 <- ggplot(df_special_long, aes(x = x, y = Densidade, color = shape)) +
    geom_line(size = 1) +
    labs(title = "Formas Especiais da GKw") +
    theme_minimal()

  # Casos especiais conforme artigo
  # 1. Kumaraswamy (Kw): λ = γ = 1, δ = 0
  kw_density <- dgkw(x_grid, alpha = 2.0, beta = 3.0, gamma = 1.0, delta = 0.0, lambda = 1.0)
  # 2. McDonald (Mc): α = β = 1
  mc_density <- dgkw(x_grid, alpha = 1.0, beta = 1.0, gamma = 2.0, delta = 1.5, lambda = 2.0)
  # 3. Beta: α = β = λ = 1
  beta_density <- dgkw(x_grid, alpha = 1.0, beta = 1.0, gamma = 2.0, delta = 3.0, lambda = 1.0)
  # 4. Beta Kumaraswamy (BKw): λ = 1
  bkw_density <- dgkw(x_grid, alpha = 2.0, beta = 2.0, gamma = 2.0, delta = 1.0, lambda = 1.0)
  # 5. Kumaraswamy-Kumaraswamy (KwKw): γ = 1
  kwkw_density <- dgkw(x_grid, alpha = 2.0, beta = 2.0, gamma = 1.0, delta = 1.0, lambda = 2.0)
  # 6. Kumaraswamy Exponenciada (EKw): δ = 0, γ = 1
  ekw_density <- dgkw(x_grid, alpha = 2.0, beta = 2.0, gamma = 1.0, delta = 0.0, lambda = 2.0)
  # 7. Beta Power (BP): α = 1, β = 1
  bp_density <- dgkw(x_grid, alpha = 1.0, beta = 1.0, gamma = 2.0, delta = 1.0, lambda = 2.0)

  df_special_cases <- data.frame(
    x = x_grid,
    "Kw" = kw_density,
    "Mc" = mc_density,
    "Beta" = beta_density,
    "BKw" = bkw_density,
    "KwKw" = kwkw_density,
    "EKw" = ekw_density,
    "BP" = bp_density
  )

  df_special_cases_long <- reshape2::melt(df_special_cases, id.vars = "x",
                                          variable.name = "Distribution", value.name = "Densidade")

  p7 <- ggplot(df_special_cases_long, aes(x = x, y = Densidade, color = Distribution)) +
    geom_line(size = 1) +
    labs(title = "Casos Especiais da GKw",
         subtitle = "Usando dgkw() com parâmetros específicos") +
    theme_minimal()

  # Exibir os gráficos
  grid.arrange(p1, p2, p3, p4, p5, p6, p7, ncol = 2)
}

#-----------------------------------------------------------
# 5. Função Aprimorada para Estudo de Simulação
#-----------------------------------------------------------
simulation_study <- function(n_sim = 100, n_samples = 500, true_params, models = c("GKw", "Kw", "BKw", "KwKw", "EKw"), verbose = FALSE) {
  # Nome dos parâmetros
  param_names <- c("alpha", "beta", "gamma", "delta", "lambda")

  # Criar lista para armazenar resultados por modelo
  results <- list()

  # Para cada modelo
  for (model in models) {
    # Matrizes para armazenar estimativas e erros padrão
    estimates <- matrix(NA, nrow = n_sim, ncol = length(param_names))
    colnames(estimates) <- param_names

    std_errors <- matrix(NA, nrow = n_sim, ncol = length(param_names))
    colnames(std_errors) <- param_names

    # Vetor para armazenar indicadores de convergência (1 = convergiu, 0 = não convergiu)
    convergence <- numeric(n_sim)

    # Vetores para armazenar métricas de ajuste
    loglik_values <- numeric(n_sim)
    aic_values <- numeric(n_sim)
    bic_values <- numeric(n_sim)

    # Loop de simulação
    for(i in 1:n_sim) {
      # Simula dados a partir da gKw
      sim_data <- rgkw(n_samples,
                       alpha = true_params["alpha"],
                       beta = true_params["beta"],
                       gamma = true_params["gamma"],
                       delta = true_params["delta"],
                       lambda = true_params["lambda"])

      # Ajuste do modelo
      fit <- try(fit_gkw_model(sim_data, verbose = verbose, model = model), silent = TRUE)

      # Se o ajuste for bem-sucedido e convergir
      if (!inherits(fit, "try-error") && is.list(fit) && fit$convergence == 0) {
        estimates[i, ] <- fit$parameters
        std_errors[i, ] <- fit$std_errors
        loglik_values[i] <- fit$loglik
        aic_values[i] <- fit$AIC
        bic_values[i] <- fit$BIC
        convergence[i] <- 1
      } else {
        convergence[i] <- 0
      }
    }

    # Armazenar resultados para este modelo
    results[[model]] <- list(
      estimates = estimates,
      std_errors = std_errors,
      convergence = convergence,
      logLik = loglik_values,
      AIC = aic_values,
      BIC = bic_values
    )
  }

  # Resultados globais
  global_results <- list(
    results = results,
    n_sim = n_sim,
    n_samples = n_samples,
    true_params = true_params,
    models = models
  )

  return(global_results)
}

#-----------------------------------------------------------
# 6. Função Aprimorada para Sumarizar os Resultados
#-----------------------------------------------------------
summarize_simulation <- function(sim_results) {
  results <- sim_results$results
  true_params <- sim_results$true_params
  n_sim <- sim_results$n_sim
  models <- sim_results$models

  # Criar dataframe para tabela de resumo
  summary_df <- data.frame()

  # Para cada modelo
  for (model in models) {
    model_results <- results[[model]]
    estimates <- model_results$estimates

    # Cálculo das estatísticas para cada parâmetro
    for (param in colnames(estimates)) {
      param_estimates <- estimates[, param]
      param_estimates <- param_estimates[!is.na(param_estimates)]

      # Se não há estimativas válidas, continue
      if (length(param_estimates) == 0) next

      mean_est <- mean(param_estimates)
      sd_est <- sd(param_estimates)
      bias <- mean_est - true_params[param]
      rmse <- sqrt(mean((param_estimates - true_params[param])^2))

      # Adicionar linha ao dataframe
      summary_df <- rbind(summary_df, data.frame(
        Model = model,
        Parameter = param,
        True = true_params[param],
        Mean = mean_est,
        SD = sd_est,
        Bias = bias,
        RMSE = rmse
      ))
    }

    # Taxa de convergência
    conv_rate <- mean(model_results$convergence) * 100

    # Métricas de ajuste médias
    mean_loglik <- mean(model_results$logLik, na.rm = TRUE)
    mean_aic <- mean(model_results$AIC, na.rm = TRUE)
    mean_bic <- mean(model_results$BIC, na.rm = TRUE)

    # Adicionar ao dataframe
    summary_df <- rbind(summary_df, data.frame(
      Model = model,
      Parameter = "Convergência (%)",
      True = NA,
      Mean = conv_rate,
      SD = NA,
      Bias = NA,
      RMSE = NA
    ))

    summary_df <- rbind(summary_df, data.frame(
      Model = model,
      Parameter = "LogLik Médio",
      True = NA,
      Mean = mean_loglik,
      SD = NA,
      Bias = NA,
      RMSE = NA
    ))

    summary_df <- rbind(summary_df, data.frame(
      Model = model,
      Parameter = "AIC Médio",
      True = NA,
      Mean = mean_aic,
      SD = NA,
      Bias = NA,
      RMSE = NA
    ))

    summary_df <- rbind(summary_df, data.frame(
      Model = model,
      Parameter = "BIC Médio",
      True = NA,
      Mean = mean_bic,
      SD = NA,
      Bias = NA,
      RMSE = NA
    ))
  }

  return(summary_df)
}

#-----------------------------------------------------------
# 7. Função para Visualizar Resultados da Simulação
#-----------------------------------------------------------
plot_simulation_results <- function(sim_results) {
  results <- sim_results$results
  true_params <- sim_results$true_params
  models <- sim_results$models

  # Para cada parâmetro
  for (param in names(true_params)) {
    # Criar dataframe para ggplot
    df <- data.frame()

    # Para cada modelo
    for (model in models) {
      model_results <- results[[model]]
      estimates <- model_results$estimates[, param]

      # Adicionar ao dataframe
      df <- rbind(df, data.frame(
        Model = model,
        Estimate = estimates
      ))
    }

    # Remover NA
    df <- df[!is.na(df$Estimate), ]

    # Verificar se há dados
    if (nrow(df) == 0) next

    # Boxplot
    p <- ggplot(df, aes(x = Model, y = Estimate, fill = Model)) +
      geom_boxplot() +
      geom_hline(yintercept = true_params[param], linetype = "dashed", color = "red", size = 1) +
      labs(title = paste("Distribuição das Estimativas de", param),
           subtitle = paste("Valor verdadeiro:", true_params[param]),
           y = param) +
      theme_minimal() +
      theme(legend.position = "none")

    print(p)
  }

  # Gráfico para métricas de ajuste (AIC e BIC)
  aic_df <- data.frame()
  bic_df <- data.frame()

  for (model in models) {
    model_results <- results[[model]]

    aic_df <- rbind(aic_df, data.frame(
      Model = model,
      AIC = model_results$AIC
    ))

    bic_df <- rbind(bic_df, data.frame(
      Model = model,
      BIC = model_results$BIC
    ))
  }

  # Remover NA
  aic_df <- aic_df[!is.na(aic_df$AIC), ]
  bic_df <- bic_df[!is.na(bic_df$BIC), ]

  # Verificar se há dados
  if (nrow(aic_df) > 0) {
    p_aic <- ggplot(aic_df, aes(x = Model, y = AIC, fill = Model)) +
      geom_boxplot() +
      labs(title = "Distribuição dos Valores de AIC",
           y = "AIC") +
      theme_minimal() +
      theme(legend.position = "none")

    print(p_aic)
  }

  if (nrow(bic_df) > 0) {
    p_bic <- ggplot(bic_df, aes(x = Model, y = BIC, fill = Model)) +
      geom_boxplot() +
      labs(title = "Distribuição dos Valores de BIC",
           y = "BIC") +
      theme_minimal() +
      theme(legend.position = "none")

    print(p_bic)
  }
}

#-----------------------------------------------------------
# 8. Teste dos Casos Especiais da GKw conforme o artigo
#-----------------------------------------------------------

# Definir nossa grid para visualização
x_grid <- seq(0.001, 0.999, length.out = 200)

# Parâmetros base para GKw
base_params <- c(alpha = 2.0, beta = 2.5, gamma = 1.5, delta = 1.0, lambda = 2.0)

#-----------------------------------------------------------
# 8.1 Comparação direta dos casos especiais da GKw
#-----------------------------------------------------------

# Testar cada caso especial mencionado no artigo (páginas 4-8)

# 1. Distribuição Kumaraswamy (Kw): λ = γ = 1, δ = 0
cat("\n\n1. Testando distribuição Kumaraswamy (Kw): λ = γ = 1, δ = 0\n")
kw_params <- base_params
kw_params["gamma"] <- 1
kw_params["delta"] <- 0
kw_params["lambda"] <- 1
kw_density <- dgkw(x_grid,
                   alpha = kw_params["alpha"],
                   beta = kw_params["beta"],
                   gamma = kw_params["gamma"],
                   delta = kw_params["delta"],
                   lambda = kw_params["lambda"])

# Verificação: Para Kw, a densidade deve ser: f(x; α, β) = αβx^(α-1)(1-x^α)^(β-1)
# Implementar manualmente para verificar
kw_density_manual <- kw_params["alpha"] * kw_params["beta"] *
  x_grid^(kw_params["alpha"]-1) *
  (1-x_grid^kw_params["alpha"])^(kw_params["beta"]-1)

# Comparar resultados
plot(x_grid, kw_density, type="l", col="blue", lwd=2,
     main="Teste Kumaraswamy (Kw)", xlab="x", ylab="Densidade")
lines(x_grid, kw_density_manual, col="red", lty=2, lwd=2)
legend("topright", legend=c("dgkw", "Manual"), col=c("blue", "red"), lty=c(1,2), lwd=2)
cat("Diferença máxima entre dgkw e cálculo manual:", max(abs(kw_density - kw_density_manual)), "\n")

# 2. Distribuição McDonald (Mc): α = β = 1
cat("\n\n2. Testando distribuição McDonald (Mc): α = β = 1\n")
mc_params <- base_params
mc_params["alpha"] <- 1
mc_params["beta"] <- 1
mc_density <- dgkw(x_grid,
                   alpha = mc_params["alpha"],
                   beta = mc_params["beta"],
                   gamma = mc_params["gamma"],
                   delta = mc_params["delta"],
                   lambda = mc_params["lambda"])

# Verificação: Para Mc, a densidade é uma generalização da beta
plot(x_grid, mc_density, type="l", col="blue", lwd=2,
     main="Teste McDonald (Mc)", xlab="x", ylab="Densidade")

# 3. Distribuição Beta: α = β = λ = 1
cat("\n\n3. Testando distribuição Beta: α = β = λ = 1\n")
beta_params <- base_params
beta_params["alpha"] <- 1
beta_params["beta"] <- 1
beta_params["lambda"] <- 1
beta_density <- dgkw(x_grid,
                     alpha = beta_params["alpha"],
                     beta = beta_params["beta"],
                     gamma = beta_params["gamma"],
                     delta = beta_params["delta"],
                     lambda = beta_params["lambda"])

# Verificação: Para Beta, compare com dbeta
beta_density_manual <- dbeta(x_grid, shape1 = beta_params["gamma"], shape2 = beta_params["delta"] + 1)

plot(x_grid, beta_density, type="l", col="blue", lwd=2,
     main="Teste Beta", xlab="x", ylab="Densidade")
lines(x_grid, beta_density_manual, col="red", lty=2, lwd=2)
legend("topright", legend=c("dgkw", "dbeta"), col=c("blue", "red"), lty=c(1,2), lwd=2)
cat("Diferença máxima entre dgkw e dbeta:", max(abs(beta_density - beta_density_manual)), "\n")

# 4. Distribuição Beta Kumaraswamy (BKw): λ = 1
cat("\n\n4. Testando distribuição Beta Kumaraswamy (BKw): λ = 1\n")
bkw_params <- base_params
bkw_params["lambda"] <- 1
bkw_density <- dgkw(x_grid,
                    alpha = bkw_params["alpha"],
                    beta = bkw_params["beta"],
                    gamma = bkw_params["gamma"],
                    delta = bkw_params["delta"],
                    lambda = bkw_params["lambda"])

# Verificação: Para BKw, a densidade é:
# f(x; α, β, γ, δ, 1) = (αβx^(α-1))/(B(γ, δ+1)) * (1-x^α)^(β(δ+1)-1) * [1-(1-x^α)^β]^(γ-1)
bkw_density_manual <- (bkw_params["alpha"] * bkw_params["beta"] * x_grid^(bkw_params["alpha"]-1) *
                         (1-x_grid^bkw_params["alpha"])^(bkw_params["beta"]*(bkw_params["delta"]+1)-1) *
                         (1-(1-x_grid^bkw_params["alpha"])^bkw_params["beta"])^(bkw_params["gamma"]-1)) /
  beta(bkw_params["gamma"], bkw_params["delta"]+1)

plot(x_grid, bkw_density, type="l", col="blue", lwd=2,
     main="Teste Beta Kumaraswamy (BKw)", xlab="x", ylab="Densidade")
lines(x_grid, bkw_density_manual, col="red", lty=2, lwd=2)
legend("topright", legend=c("dgkw", "Manual"), col=c("blue", "red"), lty=c(1,2), lwd=2)
cat("Diferença máxima entre dgkw e cálculo manual:", max(abs(bkw_density - bkw_density_manual)), "\n")

# 5. Distribuição Kumaraswamy-Kumaraswamy (KwKw): γ = 1
cat("\n\n5. Testando distribuição Kumaraswamy-Kumaraswamy (KwKw): γ = 1\n")
kwkw_params <- base_params
kwkw_params["gamma"] <- 1
kwkw_density <- dgkw(x_grid,
                     alpha = kwkw_params["alpha"],
                     beta = kwkw_params["beta"],
                     gamma = kwkw_params["gamma"],
                     delta = kwkw_params["delta"],
                     lambda = kwkw_params["lambda"])

# Verificação: Para KwKw, a densidade é:
# f(x; α, β, 1, δ, λ) = λαβ(δ+1)x^(α-1)(1-x^α)^(β-1)[1-(1-x^α)^β]^(λ-1){1-[1-(1-x^α)^β]^λ}^δ
kwkw_density_manual <- kwkw_params["lambda"] * kwkw_params["alpha"] * kwkw_params["beta"] *
  x_grid^(kwkw_params["alpha"]-1) *
  (1-x_grid^kwkw_params["alpha"])^(kwkw_params["beta"]-1) *
  (1-(1-x_grid^kwkw_params["alpha"])^kwkw_params["beta"])^(kwkw_params["lambda"]-1) *
  (1-(1-(1-x_grid^kwkw_params["alpha"])^kwkw_params["beta"])^kwkw_params["lambda"])^kwkw_params["delta"]

plot(x_grid, kwkw_density, type="l", col="blue", lwd=2,
     main="Teste Kumaraswamy-Kumaraswamy (KwKw)", xlab="x", ylab="Densidade")
lines(x_grid, kwkw_density_manual, col="red", lty=2, lwd=2)
legend("topright", legend=c("dgkw", "Manual"), col=c("blue", "red"), lty=c(1,2), lwd=2)
cat("Diferença máxima entre dgkw e cálculo manual:", max(abs(kwkw_density - kwkw_density_manual)), "\n")

# 6. Distribuição Kumaraswamy Exponenciada (EKw): δ = 0, γ = 1
cat("\n\n6. Testando distribuição Kumaraswamy Exponenciada (EKw): δ = 0, γ = 1\n")
ekw_params <- base_params
ekw_params["gamma"] <- 1
ekw_params["delta"] <- 0
ekw_density <- dgkw(x_grid,
                    alpha = ekw_params["alpha"],
                    beta = ekw_params["beta"],
                    gamma = ekw_params["gamma"],
                    delta = ekw_params["delta"],
                    lambda = ekw_params["lambda"])

# Verificação: Para EKw, a densidade é:
# f(x; α, β, 1, 0, λ) = λαβx^(α-1)(1-x^α)^(β-1)[1-(1-x^α)^β]^(λ-1)
ekw_density_manual <- ekw_params["lambda"] * ekw_params["alpha"] * ekw_params["beta"] *
  x_grid^(ekw_params["alpha"]-1) *
  (1-x_grid^ekw_params["alpha"])^(ekw_params["beta"]-1) *
  (1-(1-x_grid^ekw_params["alpha"])^ekw_params["beta"])^(ekw_params["lambda"]-1)

plot(x_grid, ekw_density, type="l", col="blue", lwd=2,
     main="Teste Kumaraswamy Exponenciada (EKw)", xlab="x", ylab="Densidade")
lines(x_grid, ekw_density_manual, col="red", lty=2, lwd=2)
legend("topright", legend=c("dgkw", "Manual"), col=c("blue", "red"), lty=c(1,2), lwd=2)
cat("Diferença máxima entre dgkw e cálculo manual:", max(abs(ekw_density - ekw_density_manual)), "\n")

# 7. Distribuição Beta Power (BP): α = 1, β = 1
cat("\n\n7. Testando distribuição Beta Power (BP): α = 1, β = 1\n")
bp_params <- base_params
bp_params["alpha"] <- 1
bp_params["beta"] <- 1
bp_density <- dgkw(x_grid,
                   alpha = bp_params["alpha"],
                   beta = bp_params["beta"],
                   gamma = bp_params["gamma"],
                   delta = bp_params["delta"],
                   lambda = bp_params["lambda"])

# Verificação: Para BP, a densidade é:
# f(x; 1, 1, γ, δ, λ) = (λ)/(B(γ, δ+1))x^(γλ-1)(1-x^λ)^δ
bp_density_manual <- (bp_params["lambda"] * x_grid^(bp_params["gamma"]*bp_params["lambda"]-1) *
                        (1-x_grid^bp_params["lambda"])^bp_params["delta"]) /
  beta(bp_params["gamma"], bp_params["delta"]+1)

plot(x_grid, bp_density, type="l", col="blue", lwd=2,
     main="Teste Beta Power (BP)", xlab="x", ylab="Densidade")
lines(x_grid, bp_density_manual, col="red", lty=2, lwd=2)
legend("topright", legend=c("dgkw", "Manual"), col=c("blue", "red"), lty=c(1,2), lwd=2)
cat("Diferença máxima entre dgkw e cálculo manual:", max(abs(bp_density - bp_density_manual)), "\n")

#-----------------------------------------------------------
# 8.2 Gerar dados de cada distribuição especial e ajustar
#-----------------------------------------------------------

# Definir tamanho da amostra
sample_size <- 1000

# 1. Gerar dados Kumaraswamy (Kw)
set.seed(123)
data_kw <- rgkw(sample_size,
                alpha = kw_params["alpha"],
                beta = kw_params["beta"],
                gamma = kw_params["gamma"],
                delta = kw_params["delta"],
                lambda = kw_params["lambda"])

# 2. Gerar dados Beta
set.seed(456)
data_beta <- rgkw(sample_size,
                  alpha = beta_params["alpha"],
                  beta = beta_params["beta"],
                  gamma = beta_params["gamma"],
                  delta = beta_params["delta"],
                  lambda = beta_params["lambda"])

# 3. Gerar dados Beta Kumaraswamy (BKw)
set.seed(789)
data_bkw <- rgkw(sample_size,
                 alpha = bkw_params["alpha"],
                 beta = bkw_params["beta"],
                 gamma = bkw_params["gamma"],
                 delta = bkw_params["delta"],
                 lambda = bkw_params["lambda"])

# 4. Gerar dados Kumaraswamy-Kumaraswamy (KwKw)
set.seed(101)
data_kwkw <- rgkw(sample_size,
                  alpha = kwkw_params["alpha"],
                  beta = kwkw_params["beta"],
                  gamma = kwkw_params["gamma"],
                  delta = kwkw_params["delta"],
                  lambda = kwkw_params["lambda"])

# 5. Gerar dados Kumaraswamy Exponenciada (EKw)
set.seed(112)
data_ekw <- rgkw(sample_size,
                 alpha = ekw_params["alpha"],
                 beta = ekw_params["beta"],
                 gamma = ekw_params["gamma"],
                 delta = ekw_params["delta"],
                 lambda = ekw_params["lambda"])

# 6. Gerar dados Beta Power (BP)
set.seed(131)
data_bp <- rgkw(sample_size,
                alpha = bp_params["alpha"],
                beta = bp_params["beta"],
                gamma = bp_params["gamma"],
                delta = bp_params["delta"],
                lambda = bp_params["lambda"])

# Ajustar modelo GKw completo a cada conjunto de dados e verificar recuperação dos parâmetros
cat("\n\nAjuste do modelo GKw a dados de casos especiais:\n")

fit_kw <- fit_gkw_model(data_kw, model = "GKw")
cat("\nDados Kumaraswamy (Kw) - parâmetros verdadeiros vs. estimados:\n")
print(data.frame(
  Parameter = names(kw_params),
  True = kw_params,
  Estimated = round(fit_kw$parameters, 4)
))

fit_beta <- fit_gkw_model(data_beta, model = "GKw")
cat("\nDados Beta - parâmetros verdadeiros vs. estimados:\n")
print(data.frame(
  Parameter = names(beta_params),
  True = beta_params,
  Estimated = round(fit_beta$parameters, 4)
))

fit_bkw <- fit_gkw_model(data_bkw, model = "GKw")
cat("\nDados Beta Kumaraswamy (BKw) - parâmetros verdadeiros vs. estimados:\n")
print(data.frame(
  Parameter = names(bkw_params),
  True = bkw_params,
  Estimated = round(fit_bkw$parameters, 4)
))

fit_kwkw <- fit_gkw_model(data_kwkw, model = "GKw")
cat("\nDados Kumaraswamy-Kumaraswamy (KwKw) - parâmetros verdadeiros vs. estimados:\n")
print(data.frame(
  Parameter = names(kwkw_params),
  True = kwkw_params,
  Estimated = round(fit_kwkw$parameters, 4)
))

fit_ekw <- fit_gkw_model(data_ekw, model = "GKw")
cat("\nDados Kumaraswamy Exponenciada (EKw) - parâmetros verdadeiros vs. estimados:\n")
print(data.frame(
  Parameter = names(ekw_params),
  True = ekw_params,
  Estimated = round(fit_ekw$parameters, 4)
))

fit_bp <- fit_gkw_model(data_bp, model = "GKw")
cat("\nDados Beta Power (BP) - parâmetros verdadeiros vs. estimados:\n")
print(data.frame(
  Parameter = names(bp_params),
  True = bp_params,
  Estimated = round(fit_bp$parameters, 4)
))

#-----------------------------------------------------------
# 8.3 Visualização conjunta de todos os casos especiais
#-----------------------------------------------------------

# Criar um dataframe para todos os casos especiais
special_cases_df <- data.frame(
  x = rep(x_grid, 7),
  Density = c(kw_density, beta_density, bkw_density,
              kwkw_density, ekw_density, bp_density,
              dgkw(x_grid, alpha = base_params["alpha"],
                   beta = base_params["beta"],
                   gamma = base_params["gamma"],
                   delta = base_params["delta"],
                   lambda = base_params["lambda"])),
  Distribution = factor(rep(c("Kw", "Beta", "BKw", "KwKw", "EKw", "BP", "GKw"),
                            each = length(x_grid)))
)

# Gráfico conjunto
ggplot(special_cases_df, aes(x = x, y = Density, color = Distribution)) +
  geom_line(size = 1) +
  labs(title = "Todos os Casos Especiais da Distribuição GKw",
       subtitle = "Usando dgkw() com parâmetros específicos",
       x = "x", y = "Densidade") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Visualizar formas da densidade GKw e casos especiais
visualize_gkw_shapes()
