# Aggregation of individual conditions over partition-specific models

Models that have been derived for individual partitions are first
decomposed into conditions, that is single conditions or conditions that
are INUS (insufficient conditions that are necessary parts of a
conjunction that is unnecessary and sufficient). The individual
conditions are aggregated using UpSet plots to determine how frequent
they are individually and in combination.

## Usage

``` r
upset_conditions(df, nsets)
```

## Arguments

- df:

  Dataframe created with `partition_min` or `partition_min_inter`.

- nsets:

  Number of sets to include in plot (default is 5).

## Value

An UpSet plot produced with
[`upset`](https://rdrr.io/pkg/UpSetR/man/upset.html).

## Examples

``` r
# \donttest{
# load data from Grauvogel (2014; see data documentation)

data(Grauvogel2014)
GS_pars <- partition_min(
 dataset = Grauvogel2014,
 units = "Sender",
 cond = c("Comprehensiveness", "Linkage", "Vulnerability",
          "Repression", "Claims"),
 out = "Persistence",
 n_cut = 1, incl_cut = 0.75,
 solution = "P",
 BE_cons = rep(0.75, 3),
 BE_ncut = rep(1, 3))
upset_conditions(GS_pars, nsets = 5)
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue to the authors.
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue to the authors.
#> Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
#> ℹ Please use the `linewidth` argument instead.
#> ℹ The deprecated feature was likely used in the UpSetR package.
#>   Please report the issue to the authors.

# }
```
