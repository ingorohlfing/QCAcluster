# Diversity of cases belonging to the same partition of the pooled data

`partition_div` calculates the diversity of cases that belong to the
same partition of the clustered data (a time series; a cross section;
etc.). Diversity is measured by the number of truth table rows that the
cases of a partition cover. `partition_div` calculates the partition
diversity for all truth table rows and for the subsets of consistent and
inconsistent rows.

## Usage

``` r
partition_div(dataset, units, time, cond, out, n_cut, incl_cut)
```

## Arguments

- dataset:

  Calibrated pooled dataset that is partitioned and minimized for
  deriving the pooled solution.

- units:

  Units defining the within-dimension of data (time series)

- time:

  Periods defining the between-dimension of data (cross sections)

- cond:

  Conditions used for the pooled analysis

- out:

  Outcome used for the pooled analysis

- n_cut:

  Frequency cut-off for designating truth table rows as observed in the
  pooled data

- incl_cut:

  Inclusion cut-off for designating truth table rows as consistent in
  the pooled data

## Value

A dataframe presenting the diversity of cases belonging to the same
partition with the following columns:

- `type`: The type of the partition. `pooled` are rows with information
  on the pooled data; `between` is for cross-section partitions;
  `within` is for time-series partitions.

- `partition`: Specific dimension of the partition at hand. For
  between-dimension, the unit identifiers are included here (argument
  `units`). For the within-dimension, the time identifier are listed
  (argument `time`). The entry is `-` for the pooled data without
  partitions.

- `diversity`: Count of all truth table rows with at least one member
  belonging to a partition.

- `diversity_1`: Count of consistent truth table rows with at least one
  member belonging to a partition.

- `diversity_0`: Count of inconsistent truth table rows with at least
  one member belonging to a partition.

- `diversity_per`: Ratio of the value for `diversity` and the total
  number of truth table rows from pooled data (`diversity` value for
  pooled data).

- `diversity_per_1`: Ratio of the value for `diversity_1` and the total
  number of consistent truth table rows from pooled data (`diversity_1`
  value for pooled data).

- `diversity_per_0`: Ratio of the value for `diversity_0` and the total
  number of inconsistent truth table rows from pooled data
  (`diversity_0` value for pooled data).

## Examples

``` r
# \donttest{
# load data from Schwarz (2016; see data documentation)

data(Schwarz2016)
Schwarz_diversity <- partition_div(Schwarz2016, 
units = "country", time = "year", 
cond = c("poltrans", "ecotrans", "reform", "conflict", "attention"), 
out = "enlarge", 1, 0.8)
# }
```
