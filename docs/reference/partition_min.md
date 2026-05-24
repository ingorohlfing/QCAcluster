# Generation of conservative or parsimonious solution for individual partitions

`partition_min` decomposes clustered data into individual partitions.
For panel data, for example, these can be cross sections, time series or
both. The function derives an individual solution for each partition and
the pooled data to assess the robustness of the solutions in a
comparative perspective.

## Usage

``` r
partition_min(
  dataset,
  units,
  time,
  cond,
  out,
  n_cut,
  incl_cut,
  solution,
  BE_cons,
  WI_cons,
  BE_ncut,
  WI_ncut
)
```

## Arguments

- dataset:

  Calibrated pooled dataset that is partitioned and minimized for
  deriving the pooled solution.

- units:

  Units defining the within-dimension of data (time series). If no units
  are specified, the data is assumed to lack a dimension and be
  hierarchical.

- time:

  Periods defining the between-dimension of data (cross sections). This
  should be specified because it does not make sense to partition a time
  series into individual data points.

- cond:

  Conditions used for minimization

- out:

  Outcome used for minimization

- n_cut:

  Frequency cut-off for designating truth table rows as observed as
  opposed to designating them as remainders for the *pooled* data.

- incl_cut:

  Inclusion (a.k.a. consistency) cut-off for designating truth table
  rows as consistent for the *pooled* data.

- solution:

  A character specifying the type of solution that should be derived.
  `C` produces the conservative (or complex) solution, `P` for the
  parsimonious solution. See `partition_min_inter` for a separate
  function for the intermediate solution.

- BE_cons:

  Inclusion thresholds for creating an individual truth table for each
  cross section. They must be specified as a numeric vector. Its length
  should be equal the number of cross sections. The order of thresholds
  corresponds to the order of the cross sections in the data defined by
  the cross-section ID in the dataset (such as years in ascending
  order).

- WI_cons:

  Inclusion thresholds for creating an individual truth table for each
  time series. They must be specified as a numeric vector. Its length
  should be equal the number of time series. The order of thresholds
  corresponds to the order of the of the time-series (unit) ID in the
  dataset (such as countries in alphabetical order).

- BE_ncut:

  For *cross sections*, the minimum number of members needed for
  declaring a truth table row as relevant as opposed to designating it
  as a remainder. Must be specified as a numeric vector. Its length
  should be equal the number of cross sections. The order of thresholds
  corresponds to the order of the cross sections in the data defined by
  the cross-section ID in the dataset (such as years in ascending
  order).

- WI_ncut:

  For *time series*, the minimum number of members needed for declaring
  a truth table row as relevant as opposed to designating it as a
  remainder. Must be specified as a numeric vector. Its length should be
  equal the number of time series. The order of thresholds corresponds
  to the order of the of the time-series (unit) ID in the dataset (such
  as countries in alphabetical order).

## Value

A dataframe summarizing the partition-specific and pooled solutions with
the following columns:

- `type`: The type of the partition. `pooled` are rows with information
  on the pooled data; `between` is for cross-section partitions;
  `within` is for time-series partitions.

- `partition`: Specific dimension of the partition at hand. For
  between-dimension, the unit identifiers are included here (argument
  `units`). For the within-dimension, the time identifier are listed
  (argument `time`). The entry is `-` for the pooled data without
  partitions.

- `solution`: The solution derived for the partition or the pooled data.
  Absence of a condition is denoted by the `~` sign.

- `model`: Running ID for models. In the presence of model ambiguity,
  each model has its own row with its individual solution and
  parameters. The rest of the information in the row is duplicated, for
  example by having two rows for the within-partition 1996. The column
  `model` highlights the presence of model ambiguity by numbering all
  models belonging to the same solution. For example, if three
  consecutive rows are numbered 1, 2 and 3, then these rows belong to
  the same solution and represent model ambiguity. If a 1 in a row is
  followed by another 1, then there is no model ambiguity.

- `consistency`: The consistency score (a.k.a. inclusion score) for the
  partition of the data or the pooled data.

- `coverage`: The coverage score for the partition of the data or the
  pooled data.

## Examples

``` r
# \donttest{
# load data from Thiem (2011; see data documentation)

data(Thiem2011)

Thiem_pars <- partition_min(
  dataset = Thiem2011,
  units = "country", time = "year",
  cond = c("fedismfs", "homogtyfs", "powdifffs", "comptvnsfs", "pubsupfs", "ecodpcefs"),
  out = "memberfs",
  n_cut = 1, incl_cut = 0.8,
  solution = "P",
  BE_cons = c(0.9, 0.8, 0.7, 0.8, 0.6, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8),
  WI_cons = c(0.5, 0.8, 0.7, 0.8, 0.6, rep(0.8, 10)))
 # }
 
```
