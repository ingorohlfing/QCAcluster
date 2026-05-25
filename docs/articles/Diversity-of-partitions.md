# Diversity of partitions

``` r

library(QCAcluster)
library(knitr) # nicer html tables
```

For illustration, we use data from [Schwarz
2016](https://doi.org/10.1080/07036337.2016.1203309). The data structure
is an unbalanced panel with eight countries, ten years and 74
observations in total.
[`partition_div()`](https://ingorohlfing.github.io/QCAcluster/reference/partition_div.md)
requires as input only parameters for the calculation of the pooled
solution plus identifiers for the units (`units`) and periods (`time`).

``` r

# load data (see data description for details)
data(Schwarz2016)
Schwarz_div <- partition_div(Schwarz2016, 
                             units = "country", time = "year", 
                             cond = c("poltrans", "ecotrans", "reform", "conflict", "attention"), 
                             out = "enlarge", 1, 0.8)
kable(Schwarz_div)
```

| type | partition | diversity | diversity_1 | diversity_0 | diversity_per | diversity_per_1 | diversity_per_0 |
|:---|:---|---:|---:|---:|---:|---:|---:|
| pooled | \- | 13 | 5 | 6 | 1.0000000 | 0.3846154 | 0.4615385 |
| between | 04 | 5 | 0 | 3 | 0.3846154 | 0.0000000 | 0.2307692 |
| between | 05 | 4 | 0 | 3 | 0.3076923 | 0.0000000 | 0.2307692 |
| between | 06 | 4 | 0 | 3 | 0.3076923 | 0.0000000 | 0.2307692 |
| between | 07 | 6 | 1 | 4 | 0.4615385 | 0.0769231 | 0.3076923 |
| between | 08 | 5 | 0 | 4 | 0.3846154 | 0.0000000 | 0.3076923 |
| between | 09 | 4 | 2 | 0 | 0.3076923 | 0.1538462 | 0.0000000 |
| between | 10 | 4 | 3 | 1 | 0.3076923 | 0.2307692 | 0.0769231 |
| between | 11 | 5 | 4 | 1 | 0.3846154 | 0.3076923 | 0.0769231 |
| between | 12 | 5 | 4 | 1 | 0.3846154 | 0.3076923 | 0.0769231 |
| between | 13 | 5 | 1 | 0 | 0.3846154 | 0.0769231 | 0.0000000 |
| within | AL | 4 | 1 | 3 | 0.3076923 | 0.0769231 | 0.2307692 |
| within | BA | 5 | 0 | 3 | 0.3846154 | 0.0000000 | 0.2307692 |
| within | HR | 4 | 2 | 0 | 0.3076923 | 0.1538462 | 0.0000000 |
| within | KS | 2 | 0 | 2 | 0.1538462 | 0.0000000 | 0.1538462 |
| within | ME | 4 | 3 | 0 | 0.3076923 | 0.2307692 | 0.0000000 |
| within | MK | 7 | 5 | 0 | 0.5384615 | 0.3846154 | 0.0000000 |
| within | RS | 6 | 0 | 4 | 0.4615385 | 0.0000000 | 0.3076923 |
| within | TR | 5 | 4 | 0 | 0.3846154 | 0.3076923 | 0.0000000 |

The dataframe shows how the cases are distributed across truth table
rows. The information is presented in absolute numbers and relative
terms and for all truth table rows and the subset of consistent and
inconsistent rows.

The table shows that while the pooled data covers 11 truth table rows,
the maximum number of rows that a partition covers is 5, which equals
about 45% of all rows. The minimum number of diversity is represented by
the partition for the year 2013 because it covers only one row.

The output of the function can also be used to see whether all cases of
a partition fall into consistent or inconsistent rows. Whenever the
value for `diversity_1` or `diversity_0` is 0, there is no variation in
the type of row for the partition. In this example, this concerns 13
partitions.

### Other packages used in this vignette

Yihui Xie (2021): *knitr: A General-Purpose Package for Dynamic Report
Generation in R.* R package version 1.33.

Yihui Xie (2015): *Dynamic Documents with R and knitr.* 2nd edition.
Chapman and Hall/CRC. ISBN 978-1498716963

Yihui Xie (2014): *knitr: A Comprehensive Tool for Reproducible Research
in R.* In Victoria Stodden, Friedrich Leisch and Roger D. Peng, editors,
Implementing Reproducible Computational Research. Chapman and Hall/CRC.
ISBN 978-1466561595
