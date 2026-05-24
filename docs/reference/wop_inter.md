# Calculation of weight of partitions in pooled solution parameters for intermediate solution

`wop_inter` calculates the weight of partitions in the pooled solution
parameters (consistency, coverage) for the intermediate solution.

## Usage

``` r
wop_inter(
  dataset,
  units,
  time,
  cond,
  out,
  n_cut,
  incl_cut,
  intermediate,
  amb_selector
)
```

## Arguments

- dataset:

  Calibrated pooled dataset for partitioning and minimization

- units:

  Units defining the within-dimension of data (time series)

- time:

  Periods defining the between-dimension of data (cross sections)

- cond:

  Conditions used for the pooled analysis

- out:

  Outcome used for the pooled analysis

- n_cut:

  Frequency cut-off for designating truth table rows as observed

- incl_cut:

  Inclusion cut-off for designating truth table rows as consistent

- intermediate:

  A vector of directional expectations to derive the intermediate
  solutions

- amb_selector:

  Numerical value for selecting a single model in the presence of model
  ambiguity. Models are numbered according to their order produced by
  `minimize` by the `QCA` package.

## Value

A dataframe with information about the weight of the partitions for
pooled consistency and coverage scores and the following columns:

- `type`: The type of the partition. `between` stands for
  cross-sections; `within` stands for time series. `pooled` stands
  information about the pooled data.

- `partition`: Type of partition. For between-dimension, the unit
  identifiers are listed (argument `units`). For the within-dimension,
  the time identifiers are listed (argument `time`). The entry is `-`
  for the pooled data.

- `denom_cons`: Denominator of the consistency formula. It is the sum
  over the cases' membership in the solution.

- `num_cons`: Numerator of the consistency formula. It is the sum over
  the minimum of the cases' membership in the solution and the outcome.

- `denom_cov`: Denominator of the coverage formula. It is the sum over
  the cases' membership in the outcome.

- `num_cov`: Numerator of the coverage formula. It is the sum over the
  minimum of the cases' membership in the solution and the outcome.
  (identical to `num_cons`)

## Examples

``` r
# load data from Schwarz (2016; see data documentation)

data(Schwarz2016)
```
