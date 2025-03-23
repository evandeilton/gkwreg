rm(list = ls())
library(gkwreg)

# Função para aplicar a função de ligação inversa com base no código do link
apply_inverse_link <- function(eta, link_code, scale_factor = 1) {
  switch(as.character(link_code),
    "1" = exp(eta), # log
    "2" = scale_factor * plogis(eta), # logit
    "3" = scale_factor * pnorm(eta), # probit
    "4" = scale_factor * (atan(eta) / pi + 0.5), # cauchy
    "5" = scale_factor * (1 - exp(-exp(eta))), # cloglog
    "6" = pmax(eta, 0.001), # identity (garantindo positividade)
    stop("Link code não reconhecido: ", link_code)
  )
}

# Função para simulação do modelo de regressão Kumaraswamy Generalizado
simulate_gkw_regression <- function(n = 2000,
                                    links = c(1, 1, 1, 1, 1),
                                    scale_factors = c(1, 1, 1, 1, 1),
                                    seed = 123,
                                    true_coefs = NULL) {
  # Definir semente para reprodutibilidade
  set.seed(seed)

  # Nomes das funções de ligação para impressão
  link_names <- c("log", "logit", "probit", "cauchy", "cloglog", "identity")

  # Gerar covariáveis
  xa1 <- rnorm(n)
  xa2 <- rbinom(n, 1, 0.5)
  xb1 <- runif(n)
  xb2 <- rnorm(n, 2, 1)
  xg1 <- rpois(n, 2.0)
  xd1 <- runif(n)
  xl1 <- rpois(n, 1.5)

  # Definir os coeficientes verdadeiros se não fornecidos
  if (is.null(true_coefs)) {
    true_coefs <- list(
      alpha = c(1.0, 0.5, -1.0), # intercept, xa1, xa2
      beta = c(0.5, -0.3, 0.2), # intercept, xb1, xb2
      gamma = c(0.2, -0.2), # intercept, xg1
      delta = c(0.8, -0.2), # intercept, xd1
      lambda = c(0.3, -0.4) # intercept, xl1
    )
  }

  # Calcular os preditores lineares
  eta_alpha <- true_coefs$alpha[1] + true_coefs$alpha[2] * xa1 + true_coefs$alpha[3] * xa2
  eta_beta <- true_coefs$beta[1] + true_coefs$beta[2] * xb1 + true_coefs$beta[3] * xb2
  eta_gamma <- true_coefs$gamma[1] + true_coefs$gamma[2] * xg1
  eta_delta <- true_coefs$delta[1] + true_coefs$delta[2] * xd1
  eta_lambda <- true_coefs$lambda[1] + true_coefs$lambda[2] * xl1

  # Aplicar funções de ligação inversas com base no vetor de links
  alpha <- apply_inverse_link(eta_alpha, links[1], scale_factors[1])
  beta <- apply_inverse_link(eta_beta, links[2], scale_factors[2])
  gamma <- apply_inverse_link(eta_gamma, links[3], scale_factors[3])
  delta <- apply_inverse_link(eta_delta, links[4], scale_factors[4])
  lambda <- apply_inverse_link(eta_lambda, links[5], scale_factors[5])

  # Mostrar quais funções de ligação estão sendo usadas
  cat("Funções de ligação utilizadas:\n")
  cat("Alpha: ", link_names[links[1]], " (escala: ", scale_factors[1], ")\n", sep = "")
  cat("Beta: ", link_names[links[2]], " (escala: ", scale_factors[2], ")\n", sep = "")
  cat("Gamma: ", link_names[links[3]], " (escala: ", scale_factors[3], ")\n", sep = "")
  cat("Delta: ", link_names[links[4]], " (escala: ", scale_factors[4], ")\n", sep = "")
  cat("Lambda: ", link_names[links[5]], " (escala: ", scale_factors[5], ")\n\n", sep = "")

  # Intervalo esperado dos parâmetros
  cat("Resumo dos parâmetros após aplicação das funções de ligação:\n")
  print(summary(alpha))
  print(summary(beta))
  print(summary(gamma))
  print(summary(delta))
  print(summary(lambda))

  # Simulação dos dados
  cat("\nSimulando dados da distribuição GKw...\n")
  y <- c()
  for (i in seq_len(n)) {
    y[[i]] <- rgkw(1, alpha[i], beta[i], gamma[i], delta[i], lambda[i])
  }
  y <- unlist(y)

  # Histograma dos dados simulados
  hist(y, main = "Histograma dos dados simulados", xlab = "y", breaks = 30)

  # Criar data frame para ajuste
  data_sim <- data.frame(
    y = y,
    xa1 = xa1,
    xa2 = xa2,
    xb1 = xb1,
    xb2 = xb2,
    xg1 = xg1,
    xd1 = xd1,
    xl1 = xl1
  )

  # Ajustar o modelo
  cat("Ajustando o modelo de regressão GKw...\n")
  formula_model <- formula(y ~ xa1 + xa2 | xb1 + xb2 | xg1 | xd1 | xl1)

  fit <- tryCatch(
    {
      gkwreg(formula_model, data = data_sim, link = links, scale = scale_factors, silent = TRUE)
    },
    error = function(e) {
      cat("Erro no ajuste do modelo:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(fit)) {
    cat("O ajuste do modelo falhou. Verifique os parâmetros e dados.\n")
    return(NULL)
  }

  # Extrair e comparar coeficientes
  est_coefs <- fit$coefficients

  # Calcular erros relativos
  cat("\nComparação entre coeficientes verdadeiros e estimados:\n")
  for (param in names(true_coefs)) {
    cat("Parâmetro:", param, "- Link:", link_names[links[which(c("alpha", "beta", "gamma", "delta", "lambda") == param)]], "\n")
    true_val <- true_coefs[[param]]
    est_val <- est_coefs[[param]]
    err_abs <- est_val - true_val
    err_pct <- 100 * err_abs / true_val

    # Criar tabela para melhor visualização
    result_table <- data.frame(
      Coeficiente = paste0("Coef_", 1:length(true_val)),
      Verdadeiro = round(true_val, 4),
      Estimado = round(est_val, 4),
      Erro_Abs = round(err_abs, 4),
      Erro_Pct = round(err_pct, 2)
    )
    print(result_table)
    cat("\n")
  }

  # Análise de diagnóstico
  cat("Diagnóstico do ajuste:\n")
  cat("Convergência:", ifelse(fit$opt$convergence == 0, "Sucesso", "Falha"), "\n")
  cat("Log-verossimilhança:", round(fit$loglik, 4), "\n")
  cat("AIC:", round(fit$aic, 4), "\n")
  cat("BIC:", round(fit$bic, 4), "\n")

  # Retornar resultados para uso posterior
  return(list(
    fit = fit,
    true_coefs = true_coefs,
    est_coefs = est_coefs,
    data = data_sim,
    parameters = list(
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      delta = delta,
      lambda = lambda
    ),
    links = links,
    scale_factors = scale_factors
  ))
}

# Função para executar múltiplas simulações com diferentes combinações de links
run_multiple_simulations <- function(n_sim = 100, n = 2000, links_list, scale_factors_list) {
  results <- list()

  for (i in 1:length(links_list)) {
    cat("\n==== Simulação com combinação de links", i, "====\n")
    cat("Links:", links_list[[i]], "\n")
    cat("Scale factors:", scale_factors_list[[i]], "\n\n")

    sim_result <- simulate_gkw_regression(
      n = n,
      links = links_list[[i]],
      scale_factors = scale_factors_list[[i]],
      seed = 123 + i
    )

    results[[i]] <- sim_result
  }

  return(results)
}

# 1. Simulação com todos os links iguais a log (padrão)
sim1 <- simulate_gkw_regression(n = 2000, links = c(2, 2, 2, 2, 2), scale_factors = c(10, 10, 10, 10, 10))


# 2. Simulação com links diferentes
sim2 <- simulate_gkw_regression(
  n = 2000,
  links = c(1, 2, 3, 4, 5), # log, logit, probit, cauchy, cloglog
  scale_factors = c(1, 10, 10, 10, 10) # Ajustar escalas para garantir valores positivos
)

# 3. Executar múltiplas simulações
links_list <- list(
  c(1, 1, 1, 1, 1), # todos log
  c(1, 2, 1, 1, 1), # beta com logit
  c(1, 1, 3, 1, 1) # gamma com probit
)
scale_factors_list <- list(
  c(1, 1, 1, 1, 1),
  c(1, 10, 1, 1, 1),
  c(1, 1, 10, 1, 1)
)
multiple_results <- run_multiple_simulations(n_sim = 3, links_list = links_list, scale_factors_list = scale_factors_list)

#####
require(betareg)

## Section 4 from Ferrari and Cribari-Neto (2004)
data("GasolineYield", package = "betareg")
data("FoodExpenditure", package = "betareg")

## Table 1
gy <- betareg(yield ~ batch + temp, data = GasolineYield)
# y ~ alpha_terms | beta_terms | gamma_terms | delta_terms | lambda_terms
f1 <- gkwreg(yield ~ batch + temp, data = GasolineYield, link = c(2, 2, 1, 1, 1), scale = c(10, 10, 1, 1, 1))
gy$loglik
gy$coefficients
summary(gy)
f1$loglik
f1$coefficients
f1$sdreport

gy2 <- betareg(yield ~ batch + temp, data = GasolineYield, link = "logit", link.phi = "log");summary(gy2)
f2 <- gkwreg(yield ~ 1 | 1 | batch + temp | temp,
             data = GasolineYield,
             link = c(1, 1, 2, 1, 1),
             scale = c(1, 1, 10, 1, 1),
             );summary(f2)

gy2$loglik
gy2$coefficients
summary(gy2)
f2$loglik
f2$coefficients
f2$sdreport


summary(gy2)
summary(f2)


## Table 2
fe_beta <- betareg(I(food / income) ~ income + persons, data = FoodExpenditure)
summary(fe_beta)

## nested model comparisons via Wald and LR tests
fe_beta2 <- betareg(I(food / income) ~ income, data = FoodExpenditure)
summary(fe_beta2)


f1 <- gkwreg(yield ~ batch + temp, data = GasolineYield, link = c(2, 2, 1, 1, 1), scale = c(10, 10, 1, 1, 1))
round(f1$model$fn(), 2)
round(f1$model$gr(), 2)
round(f1$model$he(), 2)
