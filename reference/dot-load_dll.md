# Load TMB DLL

Loads a TMB DLL and registers it with TMB's symbol table. CRITICAL:
Never calls dyn.unload() to avoid finalization errors.

## Usage

``` r
.load_dll(dll_path, dll_name, verbose = FALSE)
```

## Arguments

- dll_path:

  Full path to DLL file

- dll_name:

  Base name of DLL

- verbose:

  Logical

## Value

Logical, TRUE if successful
