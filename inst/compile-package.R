# usethis::use_package("reshape2")
# usethis::use_mit_license()
# usethis::use_git()
# usethis::use_rcpp()
# usethis::use_namespace()
# usethis::use_readme_rmd()
# usethis::use_vignette("gkwreg-package", "Generalized Kumaraswamy Regression Models for Bounded Data")
# ?usethis::use_github_action_check_standard()
# usethis::use_testthat()

# Defina a flag de roxygen2 antes de documentar
# usethis::use_version()
options(roxygen2.running = TRUE)


## -------------------------------------------------------------------------- ##

devtools::clean_dll()
styler::style_pkg()
devtools::load_all(recompile = TRUE, quiet = FALSE)
Rcpp::compileAttributes()
roxygen2::roxygenize()
devtools::document()
# Verifique como o CRAN faria
rcmdcheck::rcmdcheck(args = c("--as-cran"))

# devtools::check(cran = TRUE)
# usethis::use_version(which = "dev")
devtools::test()

## -------------------------------------------------------------------------- ##
devtools::build()
devtools::build_readme()
devtools::build_manual()

devtools::install(force = TRUE)


# devtools::build_vignettes()

## -------------------------------------------------------------------------- ##

# pkgdown::init_site()
# pkgdown::build_site()
# pkgdown::build_site_github_pages()


## -------------------------------------------------------------------------- ##
## R-CRAN
## -------------------------------------------------------------------------- ##

# Instale as ferramentas necessárias
# install.packages(c("rcmdcheck", "spelling"))

# Verifique como o CRAN faria
rcmdcheck::rcmdcheck(args = c("--as-cran"))

# Verifique erros de ortografia
# spelling::spell_check_package()




#' Extrai Exemplos de Código de Arquivos Rd de um Pacote
#'
#' Percorre o diretório 'man/' de um pacote R (ou projeto com estrutura similar),
#' lê cada arquivo .Rd e extrai o código contido nas seções \examples{}
#' usando tools::Rd2ex.
#'
#' @param pkg_path O caminho para o diretório raiz do código fonte do pacote
#'   (o diretório que contém 'man/', 'R/', 'DESCRIPTION', etc.).
#'   O padrão é o diretório de trabalho atual (".").
#'
#' @return Uma lista nomeada onde cada nome é o nome do tópico (sem a extensão .Rd)
#'   e cada valor é uma string de caracteres contendo o código de exemplo
#'   extraído para aquele tópico. Retorna uma lista vazia se nenhum arquivo .Rd
#'   ou nenhum exemplo for encontrado.
#'
#' @importFrom tools Rd2ex file_path_sans_ext
#' @export # Opcional: exportar se estiver dentro de um pacote
#'
#' @examples
#' \dontrun{
#' # Exemplo: Extrair exemplos de um pacote local (substitua pelo caminho real)
#' # caminho_pacote <- "/caminho/para/codigo/fonte/meu_pacote"
#' # lista_exemplos <- extract_examples_from_package(caminho_pacote)
#'
#' # Ver exemplos de um tópico específico (ex: "minha_funcao")
#' # if ("minha_funcao" %in% names(lista_exemplos)) {
#' #   cat(lista_exemplos[["minha_funcao"]])
#' # }
#'
#' # Extrair exemplos do diretório atual (se for a raiz de um pacote)
#' # exemplos_aqui <- extract_examples_from_package(".")
#' # print(names(exemplos_aqui))
#' }
# extract_examples_from_package <- function(pkg_path = ".") {
#
#   # Carregar o pacote tools se não estiver carregado (geralmente está)
#   # requireNamespace("tools", quietly = TRUE) # Rd2ex está em tools
#
#   # 1. Encontrar arquivos Rd
#   rd_dir <- file.path(pkg_path, "man")
#   if (!dir.exists(rd_dir)) {
#     stop("Diretório 'man/' não encontrado em '", pkg_path, "'. ",
#          "Tem certeza que este é um diretório raiz de pacote R?")
#   }
#
#   rd_files <- list.files(rd_dir, pattern = "\\.Rd$", full.names = TRUE, ignore.case = TRUE)
#
#   if (length(rd_files) == 0) {
#     message("Nenhum arquivo .Rd encontrado em ", rd_dir)
#     return(list())
#   }
#
#   # 2. Inicializar lista para armazenar resultados
#   all_examples <- list()
#   message("Extraindo exemplos de ", length(rd_files), " arquivo(s) Rd...")
#
#   # 3. Iterar, extrair, ler e armazenar
#   for (rd_file in rd_files) {
#     topic_name <- tools::file_path_sans_ext(basename(rd_file))
#     # Usar arquivo temporário para a saída de Rd2ex
#     temp_out_file <- tempfile(fileext = ".R")
#     extracted_code_lines <- NULL # Para armazenar as linhas lidas
#
#     # Usar tryCatch para robustez caso Rd2ex falhe em algum arquivo
#     outcome <- tryCatch({
#       # Rd2ex escreve o código R executável no arquivo temp_out_file
#       tools::Rd2ex(rd_file, out = temp_out_file, commentDontrun = FALSE, commentDonttest = FALSE)
#
#       # Verificar se o arquivo foi criado e tem conteúdo
#       if (file.exists(temp_out_file) && file.info(temp_out_file)$size > 0) {
#         extracted_code_lines <- readLines(temp_out_file, warn = FALSE)
#       } else {
#         # Arquivo vazio ou não criado = sem exemplos (ou Rd2ex não extraiu nada)
#         extracted_code_lines <- character(0)
#       }
#       TRUE # Sucesso
#     }, error = function(e) {
#       warning("Falha ao processar exemplos de '", basename(rd_file), "': ",
#               conditionMessage(e), call. = FALSE)
#       FALSE # Falha
#     })
#
#     # Limpar arquivo temporário sempre
#     if (file.exists(temp_out_file)) {
#       unlink(temp_out_file)
#     }
#
#     # Armazenar se a extração funcionou e se há código
#     if (outcome && length(extracted_code_lines) > 0) {
#       # Juntar as linhas em uma única string com quebras de linha
#       all_examples[[topic_name]] <- paste(extracted_code_lines, collapse = "\n")
#     }
#   } # Fim do loop
#
#   message("Extração concluída. ", length(all_examples), " tópico(s) com exemplos encontrados.")
#   return(all_examples)
# }
#
# lista_exemplos <- extract_examples_from_package(".")
#
# sapply(lista_exemplos, cat)
#
# cat(lista_exemplos$summary.gkwfit)



#
# library(showtext)
# library(hexSticker)
#
# # Carregar a biblioteca ggplot2 (se ainda não estiver carregada)
# library(ggplot2)
#
# # Criar o gráfico usando expression()
# g <- ggplot() +
#   annotate(geom = "text",
#            x = 0.5,
#            y = 0.5,
#            label = expression( "(0,1)"^7 ),
#            size = 30,    # Tamanho grande
#            color = "black") +
#
#   # Definir limites (opcional, mas útil)
#   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
#
#   # Usar theme_void() para o fundo invisível
#   theme_void()
#
# sticker(g,
#         package="gkwreg", p_size=30, s_x=1.05, s_y=.8, s_width=2, s_height=1.5,
#         filename="inst/figures/gkwreg.png")

