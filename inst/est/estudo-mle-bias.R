# ===========================================================================
# Estudo de Simulação para Estimação MLE da Distribuição GKw via Newton-Raphson
# Análise de Viés e Performance do Estimador
# ===========================================================================

# Carregar pacotes necessários
library(gkwreg)    # Pacote com funções dgkw, pgkw, qgkw, rgkw, nrgkw
library(ggplot2)
library(gridExtra)
library(dplyr)
library(reshape2)
library(pbapply)
library(parallel)
library(viridis)

# ===========================================================================
# 1. Configurações do estudo de simulação
# ===========================================================================
set.seed(8271364)  # Para reprodutibilidade

# Parâmetros da simulação
n_replicas <- 500  # Número de réplicas para cada configuração
tamanhos_amostra <- c(50, 100, 250, 500)  # Tamanhos amostrais

# ===========================================================================
# 2. Cenários para os parâmetros verdadeiros
# ===========================================================================

# Caso base: valores centrais para todos os parâmetros
base_params <- c(alpha = 2.0, beta = 2.0, gamma = 1.5, delta = 1.0, lambda = 1.5)

# Valores a serem testados para cada parâmetro
alpha_values <- c(0.5, 1.0, 2.0, 5.0)
beta_values <- c(0.5, 1.0, 2.0, 5.0)
gamma_values <- c(0.5, 1.0, 1.5, 3.0)
delta_values <- c(0.1, 0.5, 1.0, 2.0)
lambda_values <- c(0.5, 1.0, 2.0, 5.0)

# ===========================================================================
# 3. Função para gerar combinações de parâmetros
# ===========================================================================
gerar_cenarios <- function() {
  # Cenário base
  cenarios <- list(base = base_params)

  # Cenários variando cada parâmetro individualmente
  for (a in alpha_values) {
    params <- base_params
    params["alpha"] <- a
    cenarios[[paste0("alpha_", a)]] <- params
  }

  for (b in beta_values) {
    params <- base_params
    params["beta"] <- b
    cenarios[[paste0("beta_", b)]] <- params
  }

  for (g in gamma_values) {
    params <- base_params
    params["gamma"] <- g
    cenarios[[paste0("gamma_", g)]] <- params
  }

  for (d in delta_values) {
    params <- base_params
    params["delta"] <- d
    cenarios[[paste0("delta_", d)]] <- params
  }

  for (l in lambda_values) {
    params <- base_params
    params["lambda"] <- l
    cenarios[[paste0("lambda_", l)]] <- params
  }

  # Cenários especiais para casos importantes da GKw
  # Kumaraswamy (Kw): λ = γ = 1, δ = 0
  cenarios[["kw"]] <- c(alpha = 2.0, beta = 2.0, gamma = 1.0, delta = 0.0, lambda = 1.0)

  # Beta: α = β = λ = 1
  cenarios[["beta"]] <- c(alpha = 1.0, beta = 1.0, gamma = 2.0, delta = 1.0, lambda = 1.0)

  # Beta Kumaraswamy (BKw): λ = 1
  cenarios[["bkw"]] <- c(alpha = 2.0, beta = 2.0, gamma = 2.0, delta = 1.0, lambda = 1.0)

  # Kumaraswamy Exponenciada (EKw): δ = 0, γ = 1
  cenarios[["ekw"]] <- c(alpha = 2.0, beta = 2.0, gamma = 1.0, delta = 0.0, lambda = 2.0)

  return(cenarios)
}

# ===========================================================================
# 4. Função de réplica de simulação
# ===========================================================================
replicar_simulacao <- function(i, n, params_true, max_tentativas = 3) {
  tryCatch({
    # Criar cópia local dos parâmetros para evitar problemas de referência
    params_true <- as.numeric(params_true)
    names(params_true) <- c("alpha", "beta", "gamma", "delta", "lambda")

    # Definir semente diferente para cada réplica para evitar resultados idênticos
    set.seed(123456 + i)

    # Gerar dados da GKw
    dados <- rgkw(
      n = n,
      alpha = params_true["alpha"],
      beta = params_true["beta"],
      gamma = params_true["gamma"],
      delta = params_true["delta"],
      lambda = params_true["lambda"]
    )

    # Valores iniciais aleatórios próximos aos verdadeiros
    params_iniciais <- params_true * runif(5, 0.8, 1.2)
    params_iniciais[params_iniciais <= 0] <- 0.1  # Garantir positividade

    # Ajuste do modelo com múltiplas tentativas
    for (tentativa in 1:max_tentativas) {
      fit <- try(nrgkw(
        start_params = params_iniciais,
        data = dados,
        verbose = FALSE,
        max_iter = 200
      ), silent = TRUE)

      # Se convergiu, retornar resultados
      if (!inherits(fit, "try-error") && !is.null(fit$convergence) && fit$convergence == 0) {
        return(list(
          converged = TRUE,
          iterations = fit$iterations,
          params_est = fit$parameters,
          std_errors = fit$std_errors,
          loglik = fit$loglik,
          params_true = params_true
        ))
      }

      # Se não convergiu, tentar novos valores iniciais
      params_iniciais <- params_true * runif(5, 0.5, 1.5)
      params_iniciais[params_iniciais <= 0] <- 0.1
    }

    # Se todas as tentativas falharam
    return(list(
      converged = FALSE,
      iterations = NA,
      params_est = rep(NA, 5),
      std_errors = rep(NA, 5),
      loglik = NA,
      params_true = params_true
    ))

  }, error = function(e) {
    # Em caso de erro, retornar falha e imprimir mensagem
    message(paste("Erro na réplica", i, ":", conditionMessage(e)))
    return(list(
      converged = FALSE,
      iterations = NA,
      params_est = rep(NA, 5),
      std_errors = rep(NA, 5),
      loglik = NA,
      params_true = if(exists("params_true")) params_true else NULL,
      error_msg = conditionMessage(e)
    ))
  })
}

# ===========================================================================
# 5. Função para executar simulação completa para um cenário
# ===========================================================================
executar_simulacao <- function(cenario_nome, params_true, n_replicas, tamanhos_amostra) {
  resultados <- list()

  for (n in tamanhos_amostra) {
    cat(sprintf("\nExecutando simulação para cenário '%s', n = %d\n", cenario_nome, n))

    # Detectar número de cores disponíveis
    num_cores <- min(parallel::detectCores() - 1, 8)

    # Configurar cluster para processamento paralelo
    cl <- makeCluster(num_cores)
    on.exit(stopCluster(cl))

    # Carregar pacotes necessários em cada nó
    clusterEvalQ(cl, {
      library(gkwreg)
      NULL  # Retornar NULL para evitar saída desnecessária
    })

    # Exportar funções E DADOS necessários para o cluster
    # Passamos explicitamente o params_true para cada nó
    current_params_true <- params_true  # Criar uma cópia local
    clusterExport(cl, c("replicar_simulacao"), envir = environment())

    # Executar réplicas em paralelo com as variáveis corretas passadas diretamente
    replicas <- pblapply(1:n_replicas, function(i, n_atual, p_true) {
      replicar_simulacao(i, n_atual, p_true)
    }, n_atual = n, p_true = current_params_true, cl = cl)

    resultados[[paste0("n", n)]] <- replicas
  }

  return(resultados)
}

# ===========================================================================
# 6. Função para calcular estatísticas a partir dos resultados
# ===========================================================================
calcular_estatisticas <- function(resultados) {
  param_names <- c("alpha", "beta", "gamma", "delta", "lambda")
  estatisticas <- list()

  for (n_key in names(resultados)) {
    replicas <- resultados[[n_key]]
    n <- as.numeric(gsub("n", "", n_key))

    # Filtrar apenas réplicas convergentes
    convergentes <- which(sapply(replicas, function(r) r$converged))
    n_convergentes <- length(convergentes)

    if (n_convergentes > 0) {
      # Valores verdadeiros (iguais para todas as réplicas)
      valores_verdadeiros <- replicas[[1]]$params_true

      # Coletar estimativas
      estimativas <- t(sapply(replicas[convergentes], function(r) r$params_est))
      colnames(estimativas) <- param_names

      # Coletar erros padrão
      erros_padrao <- t(sapply(replicas[convergentes], function(r) r$std_errors))
      colnames(erros_padrao) <- param_names

      # Calcular viés
      vies <- colMeans(estimativas, na.rm = TRUE) - valores_verdadeiros

      # Calcular desvio padrão empírico
      sd_estimativas <- apply(estimativas, 2, sd, na.rm = TRUE)

      # Calcular erro quadrático médio (MSE)
      mse <- apply((estimativas - matrix(valores_verdadeiros, nrow = n_convergentes,
                                         ncol = 5, byrow = TRUE))^2,
                   2, mean, na.rm = TRUE)

      # Calcular cobertura dos intervalos de confiança
      cobertura_ic <- sapply(1:5, function(j) {
        ci_lower <- estimativas[, j] - 1.96 * erros_padrao[, j]
        ci_upper <- estimativas[, j] + 1.96 * erros_padrao[, j]
        mean(ci_lower <= valores_verdadeiros[j] & valores_verdadeiros[j] <= ci_upper,
             na.rm = TRUE)
      })
      names(cobertura_ic) <- param_names

      # Itens (média)
      mean_iterations <- mean(sapply(replicas[convergentes],
                                     function(r) r$iterations), na.rm = TRUE)

      estatisticas[[n_key]] <- list(
        n = n,
        n_replicas = length(replicas),
        n_convergentes = n_convergentes,
        taxa_convergencia = n_convergentes / length(replicas),
        valores_verdadeiros = valores_verdadeiros,
        estimativas_medias = colMeans(estimativas, na.rm = TRUE),
        vies = vies,
        sd = sd_estimativas,
        mse = mse,
        cobertura_ic = cobertura_ic,
        media_iteracoes = mean_iterations
      )
    } else {
      # Caso não haja réplicas convergentes
      estatisticas[[n_key]] <- list(
        n = n,
        n_replicas = length(replicas),
        n_convergentes = 0,
        taxa_convergencia = 0,
        valores_verdadeiros = replicas[[1]]$params_true,
        estimativas_medias = rep(NA, 5),
        vies = rep(NA, 5),
        sd = rep(NA, 5),
        mse = rep(NA, 5),
        cobertura_ic = rep(NA, 5),
        media_iteracoes = NA
      )
    }
  }

  return(estatisticas)
}

# ===========================================================================
# 7. Funções de visualização
# ===========================================================================

# 7.1 Visualizar a densidade GKw para diferentes parâmetros
visualizar_densidades <- function(plot_grid = TRUE) {
  x <- seq(0.001, 0.999, length.out = 200)
  param_names <- c("alpha", "beta", "gamma", "delta", "lambda")
  param_values <- list(
    alpha = alpha_values,
    beta = beta_values,
    gamma = gamma_values,
    delta = delta_values,
    lambda = lambda_values
  )

  plots <- list()

  # Criar um gráfico para cada parâmetro
  for (i in 1:5) {
    param <- param_names[i]
    values <- param_values[[i]]

    densities <- matrix(0, nrow = length(x), ncol = length(values))

    for (j in 1:length(values)) {
      # Criar cópia dos parâmetros base
      params <- base_params
      # Substituir o parâmetro atual pelo valor do loop
      params[param] <- values[j]

      # Calcular densidade
      densities[, j] <- dgkw(
        x,
        alpha = params["alpha"],
        beta = params["beta"],
        gamma = params["gamma"],
        delta = params["delta"],
        lambda = params["lambda"]
      )
    }

    # Criar dataframe para ggplot
    df <- data.frame(x = rep(x, length(values)),
                     y = as.vector(densities),
                     group = factor(rep(paste0(param, " = ", values),
                                        each = length(x))))

    # Criar plot
    p <- ggplot(df, aes(x = x, y = y, color = group)) +
      geom_line(linewidth = 1) +
      labs(title = paste0("Variação de ", param),
           x = "y", y = "Densidade") +
      theme_minimal() +
      theme(legend.title = element_blank(),
            legend.position = "bottom") +
      scale_color_viridis_d()

    plots[[param]] <- p
  }

  # Casos especiais da distribuição
  special_cases <- list(
    # Kumaraswamy (Kw): λ = γ = 1, δ = 0
    kw = c(alpha = 2.0, beta = 2.0, gamma = 1.0, delta = 0.0, lambda = 1.0),
    # Beta: α = β = λ = 1
    beta = c(alpha = 1.0, beta = 1.0, gamma = 2.0, delta = 1.0, lambda = 1.0),
    # Beta Kumaraswamy (BKw): λ = 1
    bkw = c(alpha = 2.0, beta = 2.0, gamma = 2.0, delta = 1.0, lambda = 1.0),
    # Kumaraswamy-Kumaraswamy (KwKw): γ = 1
    kwkw = c(alpha = 2.0, beta = 2.0, gamma = 1.0, delta = 1.0, lambda = 2.0),
    # Kumaraswamy Exponenciada (EKw): δ = 0, γ = 1
    ekw = c(alpha = 2.0, beta = 2.0, gamma = 1.0, delta = 0.0, lambda = 2.0)
  )

  # Calcular densidades para casos especiais
  special_densities <- matrix(0, nrow = length(x), ncol = length(special_cases))

  for (j in 1:length(special_cases)) {
    params <- special_cases[[j]]

    special_densities[, j] <- dgkw(
      x,
      alpha = params["alpha"],
      beta = params["beta"],
      gamma = params["gamma"],
      delta = params["delta"],
      lambda = params["lambda"]
    )
  }

  # Criar dataframe para casos especiais
  df_special <- data.frame(
    x = rep(x, length(special_cases)),
    y = as.vector(special_densities),
    group = factor(rep(names(special_cases), each = length(x)))
  )

  # Plot para casos especiais
  p_special <- ggplot(df_special, aes(x = x, y = y, color = group)) +
    geom_line(linewidth = 1) +
    labs(title = "Casos Especiais da GKw",
         x = "y", y = "Densidade") +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom") +
    scale_color_viridis_d()

  plots[["special"]] <- p_special

  # Mostrar grid de gráficos
  if (plot_grid) {
    do.call(grid.arrange, c(plots, ncol = 2))
  }

  # Retornar a lista de plots (invisível)
  invisible(plots)
}

# 7.2 Visualizar CDFs para diferentes parâmetros
visualizar_cdfs <- function(plot_grid = TRUE) {
  x <- seq(0.001, 0.999, length.out = 200)

  # Similar à função visualizar_densidades, mas usando pgkw
  # ...

  # Placeholder para implementação completa
  cat("Função de visualização de CDFs não implementada completamente\n")
}

# 7.3 Criar gráficos de viés ± SD
criar_graficos_vies <- function(estatisticas) {
  param_names <- c("alpha", "beta", "gamma", "delta", "lambda")
  plots <- list()

  # Preparar dados para os gráficos
  for (param in param_names) {
    df <- data.frame(
      n = sapply(estatisticas, function(e) e$n),
      vies = sapply(estatisticas, function(e) e$vies[param]),
      sd = sapply(estatisticas, function(e) e$sd[param]),
      ci_lower = sapply(estatisticas, function(e) e$vies[param] - 1.96 * e$sd[param] / sqrt(e$n_convergentes)),
      ci_upper = sapply(estatisticas, function(e) e$vies[param] + 1.96 * e$sd[param] / sqrt(e$n_convergentes))
    )

    # Ordenar por tamanho amostral
    df <- df[order(df$n), ]

    # Criar gráfico
    p <- ggplot(df, aes(x = factor(n), y = vies)) +
      geom_point(linewidth = 3) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste0("Viés ± SD para ", param),
           x = "Tamanho da amostra (n)",
           y = "Viés ± SD") +
      theme_minimal()

    plots[[param]] <- p
  }

  # Mostrar grid de gráficos
  do.call(grid.arrange, c(plots, ncol = 2))

  # Retornar a lista de plots (invisível)
  invisible(plots)
}

# 7.4 Criar gráfico de taxa de convergência
criar_grafico_convergencia <- function(estatisticas) {
  df <- data.frame(
    n = sapply(estatisticas, function(e) e$n),
    taxa = sapply(estatisticas, function(e) e$taxa_convergencia)
  )

  # Ordenar por tamanho amostral
  df <- df[order(df$n), ]

  # Criar gráfico
  p <- ggplot(df, aes(x = factor(n), y = taxa)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = sprintf("%.1f%%", taxa * 100)), vjust = -0.5) +
    labs(title = "Taxa de Convergência",
         x = "Tamanho da amostra (n)",
         y = "Taxa de Convergência") +
    ylim(0, 1) +
    theme_minimal()

  print(p)

  # Retornar o gráfico (invisível)
  invisible(p)
}

# 7.5 Criar gráfico de cobertura dos intervalos de confiança
criar_grafico_cobertura <- function(estatisticas) {
  param_names <- c("alpha", "beta", "gamma", "delta", "lambda")

  # Preparar dados para o gráfico
  df <- data.frame()

  for (param in param_names) {
    df_param <- data.frame(
      n = sapply(estatisticas, function(e) e$n),
      parametro = param,
      cobertura = sapply(estatisticas, function(e) e$cobertura_ic[param])
    )
    df <- rbind(df, df_param)
  }

  # Criar gráfico
  p <- ggplot(df, aes(x = factor(n), y = cobertura, fill = parametro)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    labs(title = "Cobertura dos Intervalos de Confiança de 95%",
         x = "Tamanho da amostra (n)",
         y = "Cobertura",
         fill = "Parâmetro") +
    ylim(0, 1) +
    theme_minimal() +
    scale_fill_viridis_d()

  print(p)

  # Retornar o gráfico (invisível)
  invisible(p)
}

# 7.6 Criar gráfico comparativo de viés para diferentes valores de um parâmetro
criar_grafico_comparativo_vies <- function(resultados_por_valor, param_name) {
  # Extrair valores do parâmetro sendo variado
  param_values <- as.numeric(gsub(paste0(param_name, "_"), "", names(resultados_por_valor)))

  # Criar dataframe para o gráfico
  df <- data.frame()

  for (i in seq_along(resultados_por_valor)) {
    valor <- param_values[i]
    estatisticas <- resultados_por_valor[[i]]

    for (n_key in names(estatisticas)) {
      est <- estatisticas[[n_key]]

      df_est <- data.frame(
        n = est$n,
        param_valor = valor,
        vies_alpha = est$vies["alpha"],
        vies_beta = est$vies["beta"],
        vies_gamma = est$vies["gamma"],
        vies_delta = est$vies["delta"],
        vies_lambda = est$vies["lambda"]
      )

      df <- rbind(df, df_est)
    }
  }

  # Transformar para formato long
  df_long <- reshape2::melt(df,
                            id.vars = c("n", "param_valor"),
                            measure.vars = paste0("vies_", c("alpha", "beta", "gamma", "delta", "lambda")),
                            variable.name = "parametro",
                            value.name = "vies")

  # Limpar nomes dos parâmetros
  df_long$parametro <- gsub("vies_", "", df_long$parametro)

  # Criar gráfico para cada parâmetro
  plots <- list()

  for (p in c("alpha", "beta", "gamma", "delta", "lambda")) {
    df_param <- df_long[df_long$parametro == p, ]

    p_plot <- ggplot(df_param, aes(x = factor(n), y = vies, color = factor(param_valor), group = param_valor)) +
      geom_point(linewidth = 2) +
      geom_line() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      labs(title = paste0("Viés de ", p, " para diferentes valores de ", param_name),
           x = "Tamanho da amostra (n)",
           y = "Viés",
           color = paste0(param_name, " =")) +
      theme_minimal() +
      scale_color_viridis_d()

    plots[[p]] <- p_plot
  }

  # Mostrar grid de gráficos
  do.call(grid.arrange, c(plots, ncol = 2))

  # Retornar a lista de plots (invisível)
  invisible(plots)
}

# ===========================================================================
# 8. Execução da simulação principal
# ===========================================================================

# Para fins de teste e desenvolvimento, reduzir o número de réplicas e tamanhos amostrais
# apenas durante a depuração
n_replicas_teste <- 10  # Valor menor para testes
tamanhos_amostra_teste <- c(50, 100)  # Apenas dois tamanhos para testes

# Visualização inicial das densidades da GKw
g1 <- visualizar_densidades()

# Gerar cenários
cenarios <- gerar_cenarios()

# Lista para armazenar resultados
resultados_todos <- list()
estatisticas_todos <- list()

# Verificar ambiente e configurações antes de executar
cat("\n=== VERIFICANDO AMBIENTE DE EXECUÇÃO ===\n")
cat("Número de cores disponíveis:", parallel::detectCores(), "\n")
cat("Número de réplicas configurado:", n_replicas, "\n")
cat("Tamanhos amostrais:", paste(tamanhos_amostra, collapse=", "), "\n")
cat("Pacote gkwreg carregado:", requireNamespace("gkwreg", quietly=TRUE), "\n")
cat("Parâmetros do caso base:", paste(names(cenarios$base), "=", round(cenarios$base, 3), collapse=", "), "\n")

# Executar uma única réplica para teste
cat("\n=== TESTE DE RÉPLICA ÚNICA ===\n")
teste_replica <- replicar_simulacao(1, 50, cenarios[["base"]])
cat("Teste de réplica concluído com sucesso:", !is.null(teste_replica), "\n")
cat("Convergência da réplica teste:", teste_replica$converged, "\n")

# Executar simulação para o caso base com valores de teste
cat("\n=== EXECUTANDO SIMULAÇÃO PARA O CASO BASE (TESTE) ===\n")
resultados_teste <- executar_simulacao("base_teste", cenarios[["base"]],
                                       n_replicas_teste, tamanhos_amostra_teste)
estatisticas_teste <- calcular_estatisticas(resultados_teste)

# Se o teste foi bem-sucedido, prosseguir com a simulação completa
cat("\n=== EXECUTANDO SIMULAÇÃO COMPLETA PARA O CASO BASE ===\n")
resultados_todos[["base"]] <- executar_simulacao("base", cenarios[["base"]],
                                                 n_replicas, tamanhos_amostra)
estatisticas_todos[["base"]] <- calcular_estatisticas(resultados_todos[["base"]])

# Visualizar resultados do caso base
cat("\n=== RESULTADOS PARA O CASO BASE ===\n")
criar_graficos_vies(estatisticas_todos[["base"]])
criar_grafico_convergencia(estatisticas_todos[["base"]])
criar_grafico_cobertura(estatisticas_todos[["base"]])

# Executar simulação para variações de alpha
resultados_alpha <- list()
estatisticas_alpha <- list()

for (a in alpha_values) {
  cenario_nome <- paste0("alpha_", a)
  cat(sprintf("\n=== EXECUTANDO SIMULAÇÃO PARA %s ===\n", cenario_nome))
  resultados_alpha[[cenario_nome]] <- executar_simulacao(
    cenario_nome, cenarios[[cenario_nome]], n_replicas, tamanhos_amostra
  )
  estatisticas_alpha[[cenario_nome]] <- calcular_estatisticas(resultados_alpha[[cenario_nome]])
}

# Visualizar resultados comparativos para alpha
cat("\n=== RESULTADOS COMPARATIVOS PARA ALPHA ===\n")
criar_grafico_comparativo_vies(estatisticas_alpha, "alpha")

# NOTA: O código acima pode ser repetido para os outros parâmetros (beta, gamma, delta, lambda)
# Porém, como o tempo de execução seria muito longo, focamos apenas em alpha neste exemplo.

# ===========================================================================
# 9. Salvar resultados
# ===========================================================================
save(
  resultados_todos, estatisticas_todos,
  resultados_alpha, estatisticas_alpha,
  cenarios,
  file = "resultados_simulacao_gkw.RData"
)

# ===========================================================================
# 10. Exemplo de um caso especial - Distribuição Kumaraswamy (Kw)
# ===========================================================================
cat("\n=== EXECUTANDO SIMULAÇÃO PARA O CASO KUMARASWAMY ===\n")
resultados_todos[["kw"]] <- executar_simulacao("kw", cenarios[["kw"]],
                                               n_replicas, tamanhos_amostra)
estatisticas_todos[["kw"]] <- calcular_estatisticas(resultados_todos[["kw"]])

cat("\n=== RESULTADOS PARA O CASO KUMARASWAMY ===\n")
criar_graficos_vies(estatisticas_todos[["kw"]])
criar_grafico_convergencia(estatisticas_todos[["kw"]])

# Atualizar arquivo de resultados
save(
  resultados_todos, estatisticas_todos,
  resultados_alpha, estatisticas_alpha,
  cenarios,
  file = "resultados_simulacao_gkw.RData"
)

cat("\n=== SIMULAÇÃO CONCLUÍDA ===\n")
