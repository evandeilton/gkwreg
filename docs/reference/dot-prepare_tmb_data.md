# Prepare TMB Data for GKw Regression

Prepare TMB Data for GKw Regression

## Usage

``` r
.prepare_tmb_data(
  model_data,
  family,
  param_names,
  fixed,
  link_ints,
  link_scale_list,
  y,
  param_positions
)
```

## Arguments

- model_data:

  List of model data.

- family:

  Family name.

- param_names:

  Names of parameters.

- fixed:

  List of fixed parameters and coefficients.

- link_ints:

  List of link function integers.

- link_scale_list:

  List of link scale values.

- y:

  Response variable.

- param_positions:

  Parameter position mapping for the family.

## Value

A list with TMB data.
