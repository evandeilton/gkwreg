#' Check and Compile TMB Model Code with Persistent Cache
#'
#' @description
#' This utility function verifies whether the TMB model shared object (\code{.so/.dll}) file
#' has already been compiled for a specified DLL name. If not, it compiles the
#' corresponding C++ file and caches it in a persistent directory across R sessions.
#'
#' @details
#' The function works through the following steps:
#' \enumerate{
#'   \item Creates a persistent cache directory for storing compiled TMB models.
#'   \item Checks if a compiled file for the specified DLL already exists in the
#'         cache directory and whether it's up-to-date compared to the source code.
#'   \item If a valid compiled file exists, it loads it directly.
#'   \item If not, the function locates the corresponding C++ file inside the
#'         package, compiles it, and stores the result in the cache directory.
#'   \item Provides diagnostic messages regarding compilation status and exported symbols.
#' }
#'
#' @param dll_name A character string specifying the base name of the C++ file
#'   and the resulting DLL. The function assumes the code file is \code{dll_name.cpp}
#'   located in the \code{inst/tmb/} directory of the package.
#' @param pkg_name A character string specifying the package name. Defaults to "gkwreg".
#' @param force_recompile Logical; if \code{TRUE}, forces recompilation even if
#'   a valid compiled file exists (default is \code{FALSE}).
#' @param verbose Logical; if \code{TRUE}, prints detailed status messages
#'   (default is \code{TRUE}).
#'
#' @return
#' Returns (invisibly) a list with information about the compiled model, including path,
#' normalized path, name, and compilation status.
#' If any step fails, an error is thrown.
#'
#' @examples
#' \dontrun{
#' # Check if "gkwreg" DLL is compiled or else compile it:
#' .check_and_compile_TMB_code("gkwreg")
#'
#' # Compile a different model:
#' .check_and_compile_TMB_code("another_model")
#'
#' # Force recompilation:
#' .check_and_compile_TMB_code("gkwreg", force_recompile = TRUE)
#'
#' # Quiet mode:
#' .check_and_compile_TMB_code("gkwreg", verbose = FALSE)
#' }
#'
#' @importFrom TMB compile dynlib
#' @importFrom tools file_path_sans_ext
#' @keywords internal
.check_and_compile_TMB_code <- function(dll_name,
                                        pkg_name = "gkwreg",
                                        force_recompile = FALSE,
                                        verbose = TRUE) {
  # Helper function for logging
  log_msg <- function(...) {
    if (verbose) cat(...)
  }

  log_msg("Checking TMB model status for ", dll_name, "...\n")

  # Ensure .gkwreg_env exists
  if (!exists(".gkwreg_env")) {
    .gkwreg_env <- new.env(parent = emptyenv())
  }

  # 1. Create a persistent cache directory for TMB compilation
  # Use rappdirs if available, otherwise create directory in user's home folder
  if (requireNamespace("rappdirs", quietly = TRUE)) {
    cache_dir <- file.path(rappdirs::user_cache_dir(pkg_name), "tmb_cache")
  } else {
    # Fallback to a directory in the user's home folder
    cache_dir <- file.path(path.expand("~"), ".R", paste0(pkg_name, "_tmb_cache"))
  }

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
    log_msg("Created persistent cache directory: ", cache_dir, "\n")
  }

  # 2. Create a subdirectory for this specific R version and package version
  # This ensures recompilation when R or the package is updated
  r_version <- paste0("R", getRversion()[1, 1:2])
  pkg_version <- tryCatch(
    {
      as.character(packageVersion(pkg_name))
    },
    error = function(e) {
      "dev" # Use "dev" for package development mode
    }
  )

  version_dir <- file.path(cache_dir, r_version, pkg_version)
  if (!dir.exists(version_dir)) {
    dir.create(version_dir, recursive = TRUE)
    log_msg("Created version-specific cache directory: ", version_dir, "\n")
  }

  # 3. Define the full path for the shared object based on platform
  dll_ext <- .Platform$dynlib.ext # ".so" on Linux/Mac, ".dll" on Windows
  compiled_name <- paste0(dll_name, dll_ext)
  cached_dll_path <- file.path(version_dir, compiled_name)

  # Normalizar o caminho para evitar problemas com o til (~)
  cached_dll_path <- normalizePath(cached_dll_path, mustWork = FALSE)

  # 4. Locate the original C++ source file - try multiple standard locations
  cpp_file <- system.file("tmb", paste0(dll_name, ".cpp"), package = pkg_name)

  # If not found in "tmb", try "src/tmb"
  if (cpp_file == "") {
    cpp_file <- system.file("src/tmb", paste0(dll_name, ".cpp"), package = pkg_name)
  }

  # If still not found, try "src"
  if (cpp_file == "") {
    cpp_file <- system.file("src", paste0(dll_name, ".cpp"), package = pkg_name)
  }

  # If still not found, check standard paths for package development
  if (cpp_file == "") {
    potential_paths <- c(
      file.path("inst/tmb", paste0(dll_name, ".cpp"))
      # file.path("src/tmb", paste0(dll_name, ".cpp")),
      # file.path("src", paste0(dll_name, ".cpp"))
    )

    for (path in potential_paths) {
      if (file.exists(path)) {
        cpp_file <- normalizePath(path)
        break
      }
    }
  }

  if (cpp_file == "" || !file.exists(cpp_file)) {
    stop(
      "Could not find the TMB source file ", dll_name, ".cpp in any standard location. ",
      "Searched in inst/tmb/, src/tmb/, and src/ directories."
    )
  }

  log_msg("Found TMB source file: ", cpp_file, "\n")

  # 5. Check if the cached shared object exists and is up-to-date
  need_compile <- force_recompile

  if (file.exists(cached_dll_path) && !force_recompile) {
    log_msg("Found cached TMB model: ", cached_dll_path, "\n")

    # Check if source code is newer than compiled version
    if (file.info(cpp_file)$mtime > file.info(cached_dll_path)$mtime) {
      log_msg("Source code is newer than compiled version. Recompiling...\n")
      need_compile <- TRUE
    } else {
      # Try to load the existing library
      load_result <- try(
        {
          # Usar o caminho absoluto normalizado para evitar problemas com o til (~)
          norm_path <- normalizePath(cached_dll_path, mustWork = TRUE)

          # Check if this DLL is already loaded
          dll_loaded <- FALSE
          loaded_dlls <- getLoadedDLLs()
          for (i in seq_along(loaded_dlls)) {
            if (loaded_dlls[[i]][["path"]] == norm_path) {
              dll_loaded <- TRUE
              break
            }
          }

          # Only load if not already loaded
          if (!dll_loaded) {
            dyn.load(norm_path)
          }

          # Armazenar o caminho e nome da DLL em uma variável no ambiente para referência futura
          assign(paste0("TMB_", dll_name, "_PATH"), norm_path, envir = .gkwreg_env)
        },
        silent = TRUE
      )

      if (!inherits(load_result, "try-error")) {
        log_msg("Successfully loaded cached TMB model.\n")

        if (verbose) {
          # List registered symbols only in verbose mode
          symbols <- getDLLRegisteredRoutines(cached_dll_path)
          log_msg("Symbols registered: ", paste(names(symbols$.C), collapse = ", "), "\n")
        }

        # Create dll_info to return
        dll_info <- list(
          path = cached_dll_path,
          normalized_path = normalizePath(cached_dll_path, mustWork = FALSE),
          name = basename(tools::file_path_sans_ext(cached_dll_path)),
          compiled = FALSE
        )

        # Store in the environment
        assign(paste0("TMB_", dll_name, "_INFO"), dll_info, envir = .gkwreg_env)

        return(invisible(dll_info))
      } else {
        log_msg("Cached model exists but could not be loaded. Recompiling...\n")
        try(dyn.unload(cached_dll_path), silent = TRUE)
        need_compile <- TRUE
      }
    }
  } else if (!file.exists(cached_dll_path)) {
    log_msg("No cached model found. Compiling...\n")
    need_compile <- TRUE
  }

  # 6. Compile if needed
  if (need_compile) {
    # Copy the C++ file to a temporary location for compilation
    temp_compile_dir <- file.path(tempdir(), paste0(dll_name, "_compile_", format(Sys.time(), "%Y%m%d%H%M%S")))
    dir.create(temp_compile_dir, recursive = TRUE, showWarnings = FALSE)

    temp_cpp_file <- file.path(temp_compile_dir, basename(cpp_file))
    file.copy(cpp_file, temp_cpp_file, overwrite = TRUE)

    # Set current directory to temp compilation location
    old_dir <- getwd()
    setwd(temp_compile_dir)
    on.exit(setwd(old_dir), add = TRUE)

    log_msg("Compiling TMB model from ", temp_cpp_file, "...\n")

    # O TMB::compile requer o caminho completo ou o nome do arquivo no diretório atual
    # Vamos garantir que estamos usando o arquivo corretamente
    compile_result <- try(
      {
        # Use o arquivo com caminho completo para garantir que o TMB possa encontrá-lo
        TMB::compile(temp_cpp_file,
          safebounds = FALSE,
          safeunload = FALSE,
          verbose = verbose
        )
      },
      silent = !verbose
    )

    if (inherits(compile_result, "try-error")) {
      stop("TMB model compilation failed. Check the C++ code and dependencies.")
    }

    # Check if compilation was successful
    temp_dll_path <- file.path(temp_compile_dir, TMB::dynlib(dll_name))
    if (!file.exists(temp_dll_path)) {
      stop("Compiled file was not generated. Check permissions and file paths.")
    }

    # Copy the compiled file to the cache directory
    file.copy(temp_dll_path, cached_dll_path, overwrite = TRUE)
    log_msg("Copied compiled model to cache: ", cached_dll_path, "\n")

    # Clean up temporary files
    unlink(temp_compile_dir, recursive = TRUE)

    # Load the newly compiled model
    load_result <- try(
      {
        norm_path <- normalizePath(cached_dll_path, mustWork = TRUE)
        dyn.load(norm_path)

        # Armazenar o caminho em variável no ambiente
        assign(paste0("TMB_", dll_name, "_PATH"), norm_path, envir = .gkwreg_env)
      },
      silent = !verbose
    )

    if (inherits(load_result, "try-error")) {
      stop(
        "Failed to load compiled model. Error: ",
        attr(load_result, "condition")$message
      )
    }

    log_msg("Successfully compiled and loaded TMB model.\n")

    if (verbose) {
      # List the registered symbols
      symbols <- getDLLRegisteredRoutines(cached_dll_path)
      log_msg("Symbols registered: ", paste(names(symbols$.C), collapse = ", "), "\n")
    }
  }

  # Return o caminho da DLL e outras informações relevantes
  dll_info <- list(
    path = cached_dll_path,
    normalized_path = normalizePath(cached_dll_path, mustWork = FALSE),
    name = basename(tools::file_path_sans_ext(cached_dll_path)),
    compiled = need_compile
  )

  # Associar estas informações ao ambiente R para referência futura
  assign(paste0("TMB_", dll_name, "_INFO"), dll_info, envir = .gkwreg_env)

  invisible(dll_info)
}

#' Get TMB model path
#'
#' @param dll_name The name of the TMB model
#' @return The path to the compiled TMB model
#' @keywords internal
get_tmb_path <- function(dll_name) {
  if (!exists(".gkwreg_env")) {
    return(NULL)
  }
  get0(paste0("TMB_", dll_name, "_PATH"), envir = .gkwreg_env)
}

#' Get TMB model information
#'
#' @param dll_name The name of the TMB model
#' @return A list with information about the compiled TMB model
#' @keywords internal
get_tmb_info <- function(dll_name) {
  if (!exists(".gkwreg_env")) {
    return(NULL)
  }
  get0(paste0("TMB_", dll_name, "_INFO"), envir = .gkwreg_env)
}

#' List all available TMB model files
#'
#' @param pkg_name The package name (default: "gkwreg")
#' @return A character vector with the names of available TMB model files (without extension)
#' @export
list_tmb_models <- function(pkg_name = "gkwreg") {
  # Look for TMB models in the standard directories
  tmb_dir <- system.file("tmb", package = pkg_name)
  src_tmb_dir <- system.file("src/tmb", package = pkg_name)
  src_dir <- system.file("src", package = pkg_name)

  # Find all .cpp files
  cpp_files <- c(
    list.files(tmb_dir, pattern = "\\.cpp$", full.names = TRUE),
    list.files(src_tmb_dir, pattern = "\\.cpp$", full.names = TRUE),
    list.files(src_dir, pattern = "\\.cpp$", full.names = TRUE)
  )

  # Extract base names without extension
  model_names <- unique(tools::file_path_sans_ext(basename(cpp_files)))

  return(model_names)
}

#' Compile all available TMB models
#'
#' @param pkg_name The package name (default: "gkwreg")
#' @param force_recompile Logical; if TRUE, forces recompilation (default: FALSE)
#' @param verbose Logical; if TRUE, prints status messages (default: TRUE)
#' @return Invisibly returns a list with information about all compiled models
#' @export
compile_all_tmb_models <- function(pkg_name = "gkwreg",
                                   force_recompile = FALSE,
                                   verbose = TRUE) {
  # Get all available models
  models <- list_tmb_models(pkg_name)

  if (length(models) == 0) {
    warning("No TMB model files found in package ", pkg_name)
    return(invisible(list()))
  }

  if (verbose) {
    cat("Found", length(models), "TMB model(s):", paste(models, collapse = ", "), "\n")
  }

  # Compile each model
  results <- list()
  for (model in models) {
    if (verbose) {
      cat("\nProcessing model:", model, "\n")
      cat(paste(rep("-", nchar(model) + 16), collapse = ""), "\n")
    }

    results[[model]] <- try(
      .check_and_compile_TMB_code(
        dll_name = model,
        pkg_name = pkg_name,
        force_recompile = force_recompile,
        verbose = verbose
      ),
      silent = !verbose
    )

    if (inherits(results[[model]], "try-error")) {
      warning("Failed to compile model: ", model)
    }
  }

  if (verbose) {
    cat("\nSummary:\n")
    for (model in names(results)) {
      if (inherits(results[[model]], "try-error")) {
        cat("  ", model, ": FAILED\n")
      } else if (results[[model]]$compiled) {
        cat("  ", model, ": Compiled and loaded\n")
      } else {
        cat("  ", model, ": Loaded from cache\n")
      }
    }
  }

  invisible(results)
}


















#' #' Check and Compile TMB Model Code with Persistent Cache
#' #'
#' #' @description
#' #' This utility function verifies whether the TMB model shared object (\code{.so/.dll}) file
#' #' has already been compiled for a specified DLL name. If not, it compiles the
#' #' corresponding C++ file and caches it in a persistent directory across R sessions.
#' #'
#' #' @details
#' #' The function works through the following steps:
#' #' \enumerate{
#' #'   \item Creates a persistent cache directory for storing compiled TMB models.
#' #'   \item Checks if a compiled file for the specified DLL already exists in the
#' #'         cache directory and whether it's up-to-date compared to the source code.
#' #'   \item If a valid compiled file exists, it loads it directly.
#' #'   \item If not, the function locates the corresponding C++ file inside the
#' #'         package, compiles it, and stores the result in the cache directory.
#' #'   \item Provides diagnostic messages regarding compilation status and exported symbols.
#' #' }
#' #'
#' #' @param dll_name A character string specifying the base name of the C++ file
#' #'   and the resulting DLL (default is \code{"gkwreg"}). The function assumes the
#' #'   code file is \code{dll_name.cpp} located in the \code{inst/tmb/} directory
#' #'   of the \code{gkwreg} package.
#' #' @param force_recompile Logical; if \code{TRUE}, forces recompilation even if
#' #'   a valid compiled file exists (default is \code{FALSE}).
#' #' @param verbose Logical; if \code{TRUE}, prints detailed status messages
#' #'   (default is \code{TRUE}).
#' #'
#' #' @return
#' #' Returns (invisibly) \code{TRUE} if compilation and loading were successful.
#' #' If any step fails, an error is thrown.
#' #'
#' #' @examples
#' #' \dontrun{
#' #' # Check if "gkwreg" DLL is compiled or else compile it:
#' #' .check_and_compile_TMB_code("gkwreg")
#' #'
#' #' # Force recompilation:
#' #' .check_and_compile_TMB_code("gkwreg", force_recompile = TRUE)
#' #'
#' #' # Quiet mode:
#' #' .check_and_compile_TMB_code("gkwreg", verbose = FALSE)
#' #' }
#' #'
#' #' @importFrom TMB compile dynlib
#' #' @importFrom tools file_path_sans_ext
#' #' @keywords internal
#' .check_and_compile_TMB_code <- function(dll_name = "gkwreg",
#'                                         force_recompile = FALSE,
#'                                         verbose = TRUE) {
#'   # Helper function for logging
#'   log_msg <- function(...) {
#'     if (verbose) cat(...)
#'   }
#'
#'   log_msg("Checking TMB model status...\n")
#'
#'   # 1. Create a persistent cache directory for TMB compilation
#'   # Use rappdirs if available, otherwise create directory in user's home folder
#'   if (requireNamespace("rappdirs", quietly = TRUE)) {
#'     cache_dir <- file.path(rappdirs::user_cache_dir("gkwreg"), "tmb_cache")
#'   } else {
#'     # Fallback to a directory in the user's home folder
#'     cache_dir <- file.path(path.expand("~"), ".R", "gkwreg_tmb_cache")
#'   }
#'
#'   if (!dir.exists(cache_dir)) {
#'     dir.create(cache_dir, recursive = TRUE)
#'     log_msg("Created persistent cache directory: ", cache_dir, "\n")
#'   }
#'
#'   # 2. Create a subdirectory for this specific R version and package version
#'   # This ensures recompilation when R or the package is updated
#'   r_version <- paste0("R", getRversion()[1, 1:2])
#'   pkg_version <- tryCatch(
#'     {
#'       as.character(packageVersion("gkwreg"))
#'     },
#'     error = function(e) {
#'       "dev" # Use "dev" for package development mode
#'     }
#'   )
#'
#'   version_dir <- file.path(cache_dir, r_version, pkg_version)
#'   if (!dir.exists(version_dir)) {
#'     dir.create(version_dir, recursive = TRUE)
#'     log_msg("Created version-specific cache directory: ", version_dir, "\n")
#'   }
#'
#'   # 3. Define the full path for the shared object based on platform
#'   dll_ext <- .Platform$dynlib.ext # ".so" on Linux/Mac, ".dll" on Windows
#'   compiled_name <- paste0(dll_name, dll_ext)
#'   cached_dll_path <- file.path(version_dir, compiled_name)
#'
#'   # Normalizar o caminho para evitar problemas com o til (~)
#'   cached_dll_path <- normalizePath(cached_dll_path, mustWork = FALSE)
#'
#'   # 4. Locate the original C++ source file
#'   cpp_file <- system.file("tmb", paste0(dll_name, ".cpp"), package = "gkwreg")
#'   if (cpp_file == "") {
#'     stop("Could not find the file ", dll_name, ".cpp using system.file.")
#'   }
#'
#'   # 5. Check if the cached shared object exists and is up-to-date
#'   need_compile <- force_recompile
#'
#'   if (file.exists(cached_dll_path) && !force_recompile) {
#'     log_msg("Found cached TMB model: ", cached_dll_path, "\n")
#'
#'     # Check if source code is newer than compiled version
#'     if (file.info(cpp_file)$mtime > file.info(cached_dll_path)$mtime) {
#'       log_msg("Source code is newer than compiled version. Recompiling...\n")
#'       need_compile <- TRUE
#'     } else {
#'       # Try to load the existing library
#'       load_result <- try(
#'         {
#'           # Usar o caminho absoluto normalizado para evitar problemas com o til (~)
#'           norm_path <- normalizePath(cached_dll_path, mustWork = TRUE)
#'           dyn.load(norm_path)
#'
#'           # Armazenar o caminho e nome da DLL em uma variável global para referência futura
#'           # Isso permite que o código que usa esta DLL a encontre corretamente
#'           assign(paste0("TMB_", dll_name, "_PATH"), norm_path, envir = .gkwreg_env)
#'
#'           # Não verificamos por MakeADFun, pois é uma função R do pacote TMB,
#'           # não uma função C++ exportada pela DLL do modelo
#'         },
#'         silent = TRUE
#'       )
#'
#'       if (!inherits(load_result, "try-error")) {
#'         log_msg("Successfully loaded cached TMB model.\n")
#'
#'         if (verbose) {
#'           # List registered symbols only in verbose mode
#'           symbols <- getDLLRegisteredRoutines(cached_dll_path)
#'           log_msg("Symbols registered: ", paste(names(symbols$.C), collapse = ", "), "\n")
#'         }
#'
#'         return(invisible(TRUE))
#'       } else {
#'         log_msg("Cached model exists but could not be loaded. Recompiling...\n")
#'         try(dyn.unload(cached_dll_path), silent = TRUE)
#'         need_compile <- TRUE
#'       }
#'     }
#'   } else if (!file.exists(cached_dll_path)) {
#'     log_msg("No cached model found. Compiling...\n")
#'     need_compile <- TRUE
#'   }
#'
#'   # 6. Compile if needed
#'   if (need_compile) {
#'     # Copy the C++ file to a temporary location for compilation
#'     temp_compile_dir <- file.path(tempdir(), paste0(dll_name, "_compile_", format(Sys.time(), "%Y%m%d%H%M%S")))
#'     dir.create(temp_compile_dir, recursive = TRUE, showWarnings = FALSE)
#'
#'     temp_cpp_file <- file.path(temp_compile_dir, basename(cpp_file))
#'     file.copy(cpp_file, temp_cpp_file, overwrite = TRUE)
#'
#'     # Set current directory to temp compilation location
#'     old_dir <- getwd()
#'     setwd(temp_compile_dir)
#'     on.exit(setwd(old_dir), add = TRUE)
#'
#'     log_msg("Compiling TMB model from ", temp_cpp_file, "...\n")
#'
#'     # O TMB::compile requer o caminho completo ou o nome do arquivo no diretório atual
#'     # Vamos garantir que estamos usando o arquivo corretamente
#'     compile_result <- try(
#'       {
#'         # Use o arquivo com caminho completo para garantir que o TMB possa encontrá-lo
#'         TMB::compile(temp_cpp_file,
#'           safebounds = FALSE,
#'           safeunload = FALSE,
#'           verbose = verbose
#'         )
#'       },
#'       silent = !verbose
#'     )
#'
#'     if (inherits(compile_result, "try-error")) {
#'       stop("TMB model compilation failed. Check the C++ code and dependencies.")
#'     }
#'
#'     # Check if compilation was successful
#'     temp_dll_path <- file.path(temp_compile_dir, TMB::dynlib(dll_name))
#'     if (!file.exists(temp_dll_path)) {
#'       stop("Compiled file was not generated. Check permissions and file paths.")
#'     }
#'
#'     # Copy the compiled file to the cache directory
#'     file.copy(temp_dll_path, cached_dll_path, overwrite = TRUE)
#'     log_msg("Copied compiled model to cache: ", cached_dll_path, "\n")
#'
#'     # Clean up temporary files
#'     unlink(temp_compile_dir, recursive = TRUE)
#'
#'     # Load the newly compiled model
#'     load_result <- try(
#'       {
#'         norm_path <- normalizePath(cached_dll_path, mustWork = TRUE)
#'         dyn.load(norm_path)
#'
#'         # Armazenar o caminho em variável global
#'         assign(paste0("TMB_", dll_name, "_PATH"), norm_path, envir = .gkwreg_env)
#'       },
#'       silent = !verbose
#'     )
#'
#'     if (inherits(load_result, "try-error")) {
#'       stop(
#'         "Failed to load compiled model. Error: ",
#'         attr(load_result, "condition")$message
#'       )
#'     }
#'
#'     log_msg("Successfully compiled and loaded TMB model.\n")
#'
#'     if (verbose) {
#'       # List the registered symbols
#'       symbols <- getDLLRegisteredRoutines(cached_dll_path)
#'       log_msg("Symbols registered: ", paste(names(symbols$.C), collapse = ", "), "\n")
#'     }
#'   }
#'
#'   # Return o caminho da DLL e outras informações relevantes
#'   dll_info <- list(
#'     path = cached_dll_path,
#'     normalized_path = normalizePath(cached_dll_path, mustWork = FALSE),
#'     name = basename(tools::file_path_sans_ext(cached_dll_path)),
#'     compiled = need_compile
#'   )
#'
#'   # Associar estas informações ao ambiente R para referência futura
#'   assign(paste0("TMB_", dll_name, "_INFO"), dll_info, envir = .gkwreg_env)
#'
#'   invisible(dll_info)
#' }
#'
#' #' Get TMB model path
#' #'
#' #' @param dll_name The name of the TMB model
#' #' @return The path to the compiled TMB model
#' #' @keywords internal
#' get_tmb_path <- function(dll_name = "gkwreg") {
#'   get0(paste0("TMB_", dll_name, "_PATH"), envir = .gkwreg_env)
#' }
#'
#' #' Get TMB model information
#' #'
#' #' @param dll_name The name of the TMB model
#' #' @return A list with information about the compiled TMB model
#' #' @keywords internal
#' get_tmb_info <- function(dll_name = "gkwreg") {
#'   get0(paste0("TMB_", dll_name, "_INFO"), envir = .gkwreg_env)
#' }



#' #' Check and Compile TMB Model Code
#' #'
#' #' @description
#' #' This utility function verifies whether the TMB model shared object (\code{.so}) file
#' #' has already been compiled for a specified DLL name. If not, it compiles the
#' #' corresponding C++ file (located via \code{system.file()}) within a dedicated
#' #' temporary directory, generates the \code{.so} file, and then loads it.
#' #'
#' #' @details
#' #' The function works through the following steps:
#' #' \enumerate{
#' #'   \item Creates a unique temporary directory for compiling the TMB code, named
#' #'         after the DLL (e.g. \code{"gkwreg_compile"}).
#' #'   \item Checks if a \code{.so} file for the specified DLL already exists in the
#' #'         temporary directory. If so, it attempts to load it directly.
#' #'   \item If the \code{.so} does not exist or cannot be loaded, the function locates
#' #'         the corresponding C++ file (e.g. \code{dll_name.cpp}) inside the
#' #'         \code{gkwreg} package using \code{system.file()}.
#' #'   \item Copies the \code{.cpp} file to the temporary directory, compiles it using
#' #'         \code{TMB::compile()}, and then attempts to load the newly generated
#' #'         \code{.so} file.
#' #'   \item Provides some diagnostic messages regarding compilation success and the
#' #'         symbols exported by the shared object.
#' #' }
#' #'
#' #' @param dll_name A character string specifying the base name of the C++ file
#' #'   and the resulting DLL (default is \code{"gkwreg"}). The function assumes the
#' #'   code file is \code{dll_name.cpp} located in the \code{inst/tmb/} directory
#' #'   of the \code{gkwreg} package.
#' #'
#' #' @return
#' #' Returns (invisibly) \code{TRUE} if compilation and loading were successful.
#' #' If any step fails, an error is thrown.
#' #'
#' #' @examples
#' #' \dontrun{
#' #' # Check if "gkwreg" DLL is compiled or else compile it:
#' #' .check_and_compile_TMB_code("gkwreg")
#' #' }
#' #'
#' #' @importFrom TMB compile dynlib
#' #' @importFrom utils getDLLRegisteredRoutines
#' #' @keywords internal
#' .check_and_compile_TMB_code <- function(dll_name = "gkwreg") {
#'   cat("Checking if the TMB model is already compiled...\n")
#'
#'   # 1. Create a dedicated temporary directory for compilation
#'   tmp_dir <- file.path(tempdir(), paste0(dll_name, "_compile"))
#'   if (!dir.exists(tmp_dir)) {
#'     dir.create(tmp_dir, recursive = TRUE)
#'   }
#'
#'   # 2. Define the full path for the .so to be generated in the temporary directory
#'   tmp_dll_path <- file.path(tmp_dir, TMB::dynlib(dll_name))
#'
#'   # 3. If the .so file already exists, attempt to load it
#'   if (file.exists(tmp_dll_path)) {
#'     load_result <- try(dyn.load(tmp_dll_path), silent = TRUE)
#'     if (!inherits(load_result, "try-error")) {
#'       cat("Model already compiled and successfully loaded from", tmp_dll_path, "\n")
#'       return(invisible(TRUE))
#'     } else {
#'       cat("A .so file exists but could not be loaded. Recompiling...\n")
#'       try(dyn.unload(tmp_dll_path), silent = TRUE)
#'       file.remove(tmp_dll_path)
#'     }
#'   }
#'
#'   # 4. Locate the .cpp file using system.file within the "gkwreg" package
#'   cpp_file <- system.file("tmb", paste0(dll_name, ".cpp"), package = "gkwreg")
#'   if (cpp_file == "") {
#'     stop("Could not find the file ", dll_name, ".cpp using system.file.")
#'   }
#'
#'   # 5. Copy the .cpp file to the temporary directory
#'   tmp_cpp_file <- file.path(tmp_dir, basename(cpp_file))
#'   if (!file.copy(cpp_file, tmp_cpp_file, overwrite = TRUE)) {
#'     stop("Failed to copy the file ", cpp_file, " to the temporary directory ", tmp_dir)
#'   }
#'
#'   cat("Compiling the model from", tmp_cpp_file, "...\n")
#'
#'   # 6. Compile the .cpp file using TMB::compile in the temporary directory
#'   compile_result <- try(
#'     TMB::compile(tmp_cpp_file,
#'       safebounds = FALSE,
#'       safeunload = FALSE,
#'       verbose = TRUE
#'     ),
#'     silent = FALSE
#'   )
#'
#'   if (inherits(compile_result, "try-error")) {
#'     stop("Model compilation failed. Check the C++ code and dependencies.")
#'   }
#'
#'   # 7. Check if the .so file was generated
#'   if (!file.exists(tmp_dll_path)) {
#'     stop(".so file was not generated after compilation. Check permissions and file paths.")
#'   }
#'
#'   # 8. Attempt to load the DLL from the temporary directory
#'   load_result <- try(dyn.load(tmp_dll_path), silent = FALSE)
#'   if (inherits(load_result, "try-error")) {
#'     stop(
#'       "Failed to load the .so file after compilation. Error: ",
#'       attr(load_result, "condition")$message
#'     )
#'   }
#'
#'   cat("Model successfully compiled and loaded from", tmp_dll_path, "\n")
#'
#'   # 9. List the registered symbols for debugging
#'   symbols <- getDLLRegisteredRoutines(tmp_dll_path)
#'   cat("Symbols registered in the DLL:", paste(names(symbols$.C), collapse = ", "), "\n")
#'
#'   invisible(TRUE)
#' }
