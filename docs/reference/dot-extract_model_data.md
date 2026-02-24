# Extract Model Data for GKw Regression

Extract Model Data for GKw Regression

## Usage

``` r
.extract_model_data(
  formula_list,
  data,
  subset,
  weights,
  na.action,
  offset,
  contrasts,
  original_call
)
```

## Arguments

- formula_list:

  List of formulas for each parameter.

- data:

  Data frame containing the variables.

- subset:

  Optional subset specification.

- weights:

  Optional weights.

- na.action:

  Function to handle missing values.

- offset:

  Optional offset.

- contrasts:

  List of contrasts for factors.

- original_call:

  The original function call.

## Value

A list of model data including frames, matrices, etc.
