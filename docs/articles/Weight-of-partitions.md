# Weight of partitions

``` r

library(QCAcluster)
library(knitr) # nicer html tables
```

### Conservative and parsimonious solution

We use the data from [Thiem
(2011)](https://doi.org/10.1017/S1755773910000251) for illustrating how
the function [`wop()`](../reference/wop.md) calculates the weights of
partitions. The weight of a partition is defined on the level of
individual models and can be calculated for the *consistency* and
*coverage* value of a model that has been derived from the pooled data.
The weight of a partition for the consistency value of the pooled
solution is calculated by applying the consistency formula only to the
cases that belong to a partition. The weight of partition is calculated
in *absolute* terms by calculating separately its contribution to the
numerator and denominator of the formula. When one divides the
partition-specific absolute contribution to the numerator by the
contribution to the denominator, then one receives the
partition-specific consistency or coverage score (depending on the type
of formula).

The arguments of the functions are:

- `n_cut`: Frequency threshold for pooled data
- `incl_cut`: Inclusion threshold (a.k.a. consistency threshold) for
  pooled data
- `solution`: Either `C` for conservative solution (a.k.a. complex
  solution) or `P` for parsimonious solution
- `amb_selector`: Numerical value for selecting a single model in the
  presence of model ambiguity. Models are numbered according to their
  order produced by minimize by the QCA package.

``` r

# load data (see data description for details)
data("Thiem2011")
# calculate weight of partitions
wop_pars <- wop(
  dataset = Thiem2011,
  units = "country", time = "year",
  cond = c("fedismfs", "homogtyfs", "powdifffs", "comptvnsfs", "pubsupfs", "ecodpcefs"),
  out = "memberfs",
  n_cut = 6, incl_cut = 0.8,
  solution = "P",
  amb_selector = 1)
kable(wop_pars)
```

| type    | partition | denom_cons | num_cons | denom_cov | num_cov |
|:--------|:----------|-----------:|---------:|----------:|--------:|
| pooled  | \-        |      80.64 |    72.39 |    101.76 |   72.39 |
| between | 1996      |       7.46 |     5.49 |      7.95 |    5.49 |
| between | 1997      |       6.77 |     5.25 |      7.95 |    5.25 |
| between | 1998      |       6.79 |     5.34 |      7.29 |    5.34 |
| between | 1999      |       7.31 |     6.31 |      8.29 |    6.31 |
| between | 2000      |       6.66 |     6.45 |      9.91 |    6.45 |
| between | 2001      |       6.80 |     6.57 |      9.91 |    6.57 |
| between | 2002      |       6.85 |     6.59 |      9.91 |    6.59 |
| between | 2003      |       7.24 |     7.02 |     10.13 |    7.02 |
| between | 2004      |       7.73 |     7.68 |     11.04 |    7.68 |
| between | 2005      |       8.22 |     8.09 |     11.04 |    8.09 |
| between | 2006      |       8.81 |     7.60 |      8.34 |    7.60 |
| within  | AT        |       4.32 |     1.89 |      2.28 |    1.89 |
| within  | BE        |       9.99 |     7.71 |      7.86 |    7.71 |
| within  | DE        |       9.65 |     9.65 |     10.73 |    9.65 |
| within  | DK        |       1.73 |     1.59 |      7.46 |    1.59 |
| within  | ES        |       8.41 |     7.95 |      9.83 |    7.95 |
| within  | FI        |       0.60 |     0.60 |      5.26 |    0.60 |
| within  | FR        |       9.88 |     9.86 |     10.87 |    9.86 |
| within  | GR        |       1.46 |     1.29 |      3.80 |    1.29 |
| within  | IE        |       0.73 |     0.21 |      0.21 |    0.21 |
| within  | IT        |       9.21 |     9.21 |     10.77 |    9.21 |
| within  | LU        |       0.33 |     0.33 |      3.80 |    0.33 |
| within  | NL        |       6.26 |     5.73 |      7.06 |    5.73 |
| within  | PT        |       0.81 |     0.81 |      3.80 |    0.81 |
| within  | SE        |       6.94 |     5.27 |      7.16 |    5.27 |
| within  | UK        |      10.32 |    10.29 |     10.87 |   10.29 |

When one aggregates the partition-specific absolute weights for the
between-dimension or within-dimension, one gets the absolute value for
the pooled solution. We illustrate this with the following chunk

``` r

# sum over all cross-sections for consistency denominator
sum(wop_pars[wop_pars$type == "between", ]$denom_cons)
#> [1] 80.64
# sum over all time series for coverage  numerator
sum(wop_pars[wop_pars$type == "within", ]$num_cov)
#> [1] 72.39
```

On the basis of the absolute weights, one can calculate the *relative
weight* of a partition by dividing its absolute contribution by the
corresponding value for the pooled solution.

``` r

# relative contribution of cross sections to denominator for consistency
wop_between  <- wop_pars[wop_pars$type == "between", ]
wop_between$rel_denom_cons <- round(wop_between$denom_cons / 
  sum(wop_between$denom_cons), digits = 2)
kable(wop_between)
```

|     | type    | partition | denom_cons | num_cons | denom_cov | num_cov | rel_denom_cons |
|:----|:--------|:----------|-----------:|---------:|----------:|--------:|---------------:|
| 2   | between | 1996      |       7.46 |     5.49 |      7.95 |    5.49 |           0.09 |
| 3   | between | 1997      |       6.77 |     5.25 |      7.95 |    5.25 |           0.08 |
| 4   | between | 1998      |       6.79 |     5.34 |      7.29 |    5.34 |           0.08 |
| 5   | between | 1999      |       7.31 |     6.31 |      8.29 |    6.31 |           0.09 |
| 6   | between | 2000      |       6.66 |     6.45 |      9.91 |    6.45 |           0.08 |
| 7   | between | 2001      |       6.80 |     6.57 |      9.91 |    6.57 |           0.08 |
| 8   | between | 2002      |       6.85 |     6.59 |      9.91 |    6.59 |           0.08 |
| 9   | between | 2003      |       7.24 |     7.02 |     10.13 |    7.02 |           0.09 |
| 10  | between | 2004      |       7.73 |     7.68 |     11.04 |    7.68 |           0.10 |
| 11  | between | 2005      |       8.22 |     8.09 |     11.04 |    8.09 |           0.10 |
| 12  | between | 2006      |       8.81 |     7.60 |      8.34 |    7.60 |           0.11 |

### Intermediate solution

The weight of partitions for intermediate solutions is produced with
[`wop_inter()`](../reference/wop_inter.md). We use data from [Schwarz
2016](https://doi.org/10.1080/07036337.2016.1203309) to illustrate the
function.

``` r

# load data (see data description for details)
data("Schwarz2016")
# calculating weight of partitions
Schwarz_wop_inter <- partition_min_inter(
  Schwarz2016,
  units = "country", time = "year",
  cond = c("poltrans", "ecotrans", "reform", "conflict", "attention"),
  out = "enlarge",
  n_cut = 1, incl_cut = 0.8,
  intermediate = c("1", "1", "1", "1", "1"))
kable(Schwarz_wop_inter)
```

### Other packages used in this vignette

Yihui Xie (2021): *knitr: A General-Purpose Package for Dynamic Report
Generation in R.* R package version 1.33.

Yihui Xie (2015): *Dynamic Documents with R and knitr.* 2nd edition.
Chapman and Hall/CRC. ISBN 978-1498716963

Yihui Xie (2014): *knitr: A Comprehensive Tool for Reproducible Research
in R.* In Victoria Stodden, Friedrich Leisch and Roger D. Peng, editors,
Implementing Reproducible Computational Research. Chapman and Hall/CRC.
ISBN 978-1466561595
