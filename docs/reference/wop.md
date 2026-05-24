# Weight of partitions for pooled solution parameters for conservative or parsimonious solution

`wop` calculates the contribution or weight of partitions for the pooled
solution parameters of consistency and coverage for the conservative or
parsimonious solution.

## Usage

``` r
wop(dataset, units, time, cond, out, n_cut, incl_cut, solution, amb_selector)
```

## Arguments

- dataset:

  Calibrated pooled dataset for partitioning and minimization of pooled
  solution.

- units:

  Units that define the within-dimension of data (time series).

- time:

  Periods that define the between-dimension of data (cross sections).

- cond:

  Conditions used for the pooled analysis.

- out:

  Outcome used for the pooled analysis.

- n_cut:

  Frequency cut-off for designating truth table rows as observed in the
  pooled analysis.

- incl_cut:

  Inclusion cut-off for designating truth table rows as consistent in
  the pooled analysis.

- solution:

  A character specifying the type of solution that should be derived.
  `C` produces the conservative (or complex) solution, `P` the
  parsimonious solution. See `wop_inter` for deriving intermediate
  solution.

- amb_selector:

  Numerical value for selecting a single model in the presence of model
  ambiguity. Models are numbered according to their order produced by
  `minimize` by the `QCA` package.

## Value

A dataframe with information about the weight of the partitions with the
following columns:

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
# load data from Thiem (EPSR, 2011; see data documentation)

data(Thiem2011)
wop_pars <- wop(
  dataset = Thiem2011,
  units = "country", time = "year",
  cond = c("fedismfs", "homogtyfs", "powdifffs", "comptvnsfs", "pubsupfs", "ecodpcefs"),
  out = "memberfs",
  n_cut = 6, incl_cut = 0.8,
  solution = "P",
  amb_selector = 1)
wop_pars
#>       type partition denom_cons num_cons denom_cov num_cov
#> 1   pooled         -      80.64    72.39    101.76   72.39
#> 2  between      1996       7.46     5.49      7.95    5.49
#> 3  between      1997       6.77     5.25      7.95    5.25
#> 4  between      1998       6.79     5.34      7.29    5.34
#> 5  between      1999       7.31     6.31      8.29    6.31
#> 6  between      2000       6.66     6.45      9.91    6.45
#> 7  between      2001       6.80     6.57      9.91    6.57
#> 8  between      2002       6.85     6.59      9.91    6.59
#> 9  between      2003       7.24     7.02     10.13    7.02
#> 10 between      2004       7.73     7.68     11.04    7.68
#> 11 between      2005       8.22     8.09     11.04    8.09
#> 12 between      2006       8.81     7.60      8.34    7.60
#> 13  within        AT       4.32     1.89      2.28    1.89
#> 14  within        BE       9.99     7.71      7.86    7.71
#> 15  within        DE       9.65     9.65     10.73    9.65
#> 16  within        DK       1.73     1.59      7.46    1.59
#> 17  within        ES       8.41     7.95      9.83    7.95
#> 18  within        FI       0.60     0.60      5.26    0.60
#> 19  within        FR       9.88     9.86     10.87    9.86
#> 20  within        GR       1.46     1.29      3.80    1.29
#> 21  within        IE       0.73     0.21      0.21    0.21
#> 22  within        IT       9.21     9.21     10.77    9.21
#> 23  within        LU       0.33     0.33      3.80    0.33
#> 24  within        NL       6.26     5.73      7.06    5.73
#> 25  within        PT       0.81     0.81      3.80    0.81
#> 26  within        SE       6.94     5.27      7.16    5.27
#> 27  within        UK      10.32    10.29     10.87   10.29
```
