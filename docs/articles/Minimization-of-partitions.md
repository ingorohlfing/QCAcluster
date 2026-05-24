# Minimization of partitions

``` r

library(QCAcluster)
library(knitr) # nicer html tables
```

Two functions allow empirical researchers to partition clustered data on
one or two dimensions and to derive solutions for the pooled data and
for each partition.

- [`partition_min()`](../reference/partition_min.md) is available for
  producing *conservative* or *parsimonious* models;  
- [`partition_min_inter()`](../reference/partition_min_inter.md) should
  be used for *intermediate* models. For programming purposes, we opted
  for a separate function for the intermediate solution.

## Panel data: Minimization of cross sections and time series

We first illustrate how one can decompose panel data on two dimensions.
In a *between-unit* perspective, the panel is partitioned into multiple
cross sections with the `time` argument that specifies the cross section
ID. In a *within-unit* perspective, the data is decomposed into multiple
time series with the `units` argument that specifies the unit (or time
series) ID. The arguments of the functions are:

- `n_cut`: Frequency threshold for pooled data
- `incl_cut`: Inclusion threshold (a.k.a. consistency threshold) for
  pooled data
- `solution` (only for
  [`partition_min()`](../reference/partition_min.md)): Either `C` for
  conservative solution (a.k.a. complex solution) or `P` for
  parsimonious solution
- `BE_cons` and `WI_cons`: Inclusion thresholds for cross sections and
  time series. The length of the numeric vector should equal the number
  of units and time series.
- `BE_ncut` and `WI_ncut`: Frequency thresholds for the cross sections
  and time series. The length of the numeric vector should equal the
  number of units and time series.

### Conservative and parsimonious solution

We first illustrate the parsimonious solution with dataset from [Thiem
(2011)](https://doi.org/10.1017/S1755773910000251).

``` r

# load data (see data description for details)
data(Thiem2011)
# partition data into time series (within-unit) and cross sections (between-unit)
Thiem_pars <- partition_min(
  dataset = Thiem2011,
  units = "country", time = "year",
  cond = c("fedismfs", "homogtyfs", "powdifffs", "comptvnsfs", 
           "pubsupfs", "ecodpcefs"),
  out = "memberfs",
  n_cut = 6, incl_cut = 0.8,
  solution = "P", # parsimonious solution
  BE_cons = c(0.9, 0.8, 0.7, 0.8, 0.85, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8),
  BE_ncut = rep(1, 11),
  WI_cons = c(0.75, 0.8, 0.9, 0.8, 0.85, rep(0.75, 10)),
  WI_ncut = rep(1, 15))
kable(Thiem_pars)
```

| type | partition | solution | model | consistency | coverage |
|:---|:---|:---|:---|---:|---:|
| pooled | \- | comptvnsfs+fedismfs \* pubsupfs | 1 | 0.8976935 | 0.7113797 |
| pooled | \- | comptvnsfs+fedismfs \* ecodpcefs | 2 | 0.8949502 | 0.7158019 |
| pooled | \- | comptvnsfs+homogtyfs \* pubsupfs | 3 | 0.8780259 | 0.7342767 |
| between | 1996 | fedismfs \* comptvnsfs | 1 | 0.9030303 | 0.3748428 |
| between | 1996 | comptvnsfs \* pubsupfs | 2 | 0.9885057 | 0.4327044 |
| between | 1997 | ~powdifffs | 1 | 0.9064748 | 0.6339623 |
| between | 1997 | comptvnsfs | 2 | 0.8910675 | 0.5144654 |
| between | 1997 | pubsupfs \* ~ecodpcefs | 3 | 0.8672769 | 0.4767296 |
| between | 1998 | comptvnsfs | 1 | 0.9288703 | 0.6090535 |
| between | 1999 | ~powdifffs+fedismfs \* ecodpcefs | 1 | 0.8876404 | 0.7623643 |
| between | 1999 | ~powdifffs+fedismfs \* ~homogtyfs+homogtyfs \* pubsupfs \* ecodpcefs | 2 | 0.8961039 | 0.7490953 |
| between | 2000 | comptvnsfs+fedismfs \* pubsupfs | 1 | 0.9684685 | 0.6508577 |
| between | 2000 | comptvnsfs+fedismfs \* ecodpcefs | 2 | 0.9417476 | 0.6851665 |
| between | 2000 | comptvnsfs+fedismfs \* ~homogtyfs+homogtyfs \* pubsupfs | 3 | 0.9708333 | 0.7053481 |
| between | 2001 | fedismfs+comptvnsfs | 1 | 0.9028436 | 0.7689203 |
| between | 2002 | fedismfs+~powdifffs+pubsupfs | 1 | 0.8149780 | 0.7467205 |
| between | 2002 | fedismfs+comptvnsfs+pubsupfs | 2 | 0.8214665 | 0.7800202 |
| between | 2003 | pubsupfs+~ecodpcefs | 1 | 0.7985213 | 0.8529121 |
| between | 2004 | fedismfs+~ecodpcefs | 1 | 0.9184290 | 0.8260870 |
| between | 2004 | pubsupfs+~ecodpcefs | 2 | 0.9081726 | 0.8958333 |
| between | 2005 | pubsupfs+~ecodpcefs | 1 | 0.9002695 | 0.9076087 |
| between | 2005 | fedismfs+_(homogtyfs+)ecodpcefs | 2 | 0.8868101 | 0.8586957 |
| between | 2006 | comptvnsfs+~pubsupfs | 1 | 0.8982118 | 0.7829736 |
| between | 2006 | ~pubsupfs+fedismfs \* ~ecodpcefs | 2 | 0.8335725 | 0.6966427 |
| within | AT | All truth table rows inconsistent | \- | NA | NA |
| within | BE | No variation in all conditions | \- | NA | NA |
| within | DE | All truth table rows consistent | \- | NA | NA |
| within | DK | ~pubsupfs | 1 | 0.8297389 | 0.9798928 |
| within | DK | ~ecodpcefs | 2 | 0.9469154 | 0.8847185 |
| within | ES | All truth table rows consistent | \- | NA | NA |
| within | FI | No variation in all conditions | \- | NA | NA |
| within | FR | All truth table rows consistent | \- | NA | NA |
| within | GR | All truth table rows inconsistent | \- | NA | NA |
| within | IE | All truth table rows inconsistent | \- | NA | NA |
| within | IT | No variation in all conditions | \- | NA | NA |
| within | LU | homogtyfs | 1 | 0.7629630 | 0.8131579 |
| within | NL | All truth table rows consistent | \- | NA | NA |
| within | PT | All truth table rows inconsistent | \- | NA | NA |
| within | SE | All truth table rows inconsistent | \- | NA | NA |
| within | UK | All truth table rows consistent | \- | NA | NA |

The output of [`partition_min()`](../reference/partition_min.md) is a
dataframe summarizing the solutions for the pooled data and the
partitions and the consistency and coverage values for the solution. The
column `model` shows whether model ambiguity is given for the pooled
data or individual partitions *if* one can derive any model from the
data in the first place.

There are different reasons why one might not be able to derive a
partition-specific solution:

- All rows could be consistent
- All rows could be inconsistent
- There is no variation across cases of a partition and all cases belong
  to the same truth table row.

When one the reason applies, it is listed in the column `solution`.

### Intermediate solution

The intermediate solution is derived with
[`partition_min_inter()`](../reference/partition_min_inter.md). The only
command that is new compared to
[`partition_min()`](../reference/partition_min.md) is `intermediate`
that is available for specifying the *directional expectations*. The
data structure for [Schwarz
2016](https://doi.org/10.1080/07036337.2016.1203309) is an unbalanced
panel with eight countries, ten years and 74 observations in total. We
assume that one is only interested in the between-unit dimension and
wants to derive one solution per cross section. For this reason, the
argument for the within-unit dimension (`unit`) is not specified.

``` r

# load data (see data description for details)
data(Schwarz2016)
# partition data into cross sections
Schwarz_inter <- partition_min_inter(
  Schwarz2016, 
  time = "year", 
  cond = c("poltrans", "ecotrans", "reform", "conflict", "attention"), 
  out = "enlarge", 
  n_cut = 1, incl_cut = 0.8, 
  WI_cons = rep(0.8, 8), BE_cons = c(0.75, 0.75, 0.75, 0.75, 0.75,
                                     0.8, 0.8, 0.8, 0.8, 0.8),
  WI_ncut = rep(1, 8), BE_ncut = rep(1, 10),
  intermediate = c("1", "1", "1", "1", "1"))
kable(Schwarz_inter)
```

| type | partition | solution | model | consistency | coverage |
|:---|:---|:---|:---|---:|---:|
| pooled | \- | poltrans \* ecotrans \* reform+poltrans \* reform \* conflict \* attention | 1 | 0.8008497 | 0.7783001 |
| between | 04 | All inconsistent | \- | NA | NA |
| between | 05 | All inconsistent | \- | NA | NA |
| between | 06 | All inconsistent | \- | NA | NA |
| between | 07 | poltrans \* ecotrans \* reform+poltrans \* reform \* ~conflict | 1 | 0.7552752 | 0.8692104 |
| between | 08 | poltrans \* ecotrans \* reform \* conflict | 1 | 0.7626173 | 0.8482275 |
| between | 09 | All consistent | \- | NA | NA |
| between | 10 | poltrans \* ecotrans \* reform \* attention+poltrans \* reform \* conflict \* attention | 1 | 0.8760953 | 0.8125806 |
| between | 11 | poltrans \* conflict \* attention+poltrans \* ecotrans \* reform \* attention | 1 | 0.8195671 | 0.9566749 |
| between | 12 | poltrans \* conflict+poltrans \* ecotrans \* reform | 1 | 0.8411864 | 0.8865839 |
| between | 13 | All consistent | \- | NA | NA |

## Multilevel data

Clustered data can be partitioned on a single dimension if there is only
one dimension as an in multilevel data where lower-level units are
nested in higher-level units. The analysis is then similar to the
partition of panel data along one dimension. We use the dataset by
[Grauvogel and von Soest
(2014)](https://doi.org/10.1017/S1755773910000251) for illustrating the
analysis of multilevel data. The study analyzes the effect of sanctions
on authoritarian regimes. The data distinguishes between the source of
the sanction (`Sender`) and the target country (`Target`). All sanctions
have been imposed by the EU, UN or US, which means that target countries
are nested in three different senders. We partition the data on the
dimension of senders to see how solutions differ across senders.

``` r

# load data (see data description for details)
data(Grauvogel2014)
# partition data by sender country (higher-level unit)
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
kable(GS_pars)
```

| type | partition | solution | model | consistency | coverage |
|:---|:---|:---|---:|---:|---:|
| pooled | \- | ~Comprehensiveness \* Claims+~Linkage \* Claims+~Repression \* Claims+~Comprehensiveness \* ~Linkage \* ~Repression+Comprehensiveness \* ~Vulnerability \* Repression | 1 | 0.7758164 | 0.7336208 |
| pooled | \- | ~Comprehensiveness \* Claims+~Linkage \* Claims+~Repression \* Claims+~Comprehensiveness \* ~Linkage \* ~Repression+Linkage \* ~Vulnerability \* Repression | 2 | 0.7776948 | 0.7245792 |
| within | EU | Vulnerability+~Comprehensiveness \* ~Repression+Linkage \* Claims | 1 | 0.6293355 | 0.9145825 |
| within | EU | Vulnerability+~Linkage \* ~Repression+Linkage \* Claims | 2 | 0.6327684 | 0.9049634 |
| within | EU | Vulnerability+Linkage \* Repression+~Repression \* Claims | 3 | 0.6310549 | 0.9022701 |
| within | EU | Vulnerability+Linkage \* Claims+~Repression \* Claims | 4 | 0.6320277 | 0.9126587 |
| within | EU | Vulnerability+~Comprehensiveness \* Linkage+~Comprehensiveness \* ~Repression+Comprehensiveness \* Claims | 5 | 0.6273610 | 0.8945748 |
| within | EU | Vulnerability+~Comprehensiveness \* Linkage+Comprehensiveness \* Repression+~Repression \* Claims | 6 | 0.6238361 | 0.9022701 |
| within | EU | Vulnerability+~Comprehensiveness \* Linkage+Comprehensiveness \* Claims+~Linkage \* ~Repression | 7 | 0.6253391 | 0.8868796 |
| within | EU | Vulnerability+~Comprehensiveness \* Linkage+Comprehensiveness \* Claims+~Repression \* Claims | 8 | 0.6261682 | 0.9022701 |
| within | EU | Vulnerability+~Comprehensiveness \* ~Repression+Comprehensiveness \* Claims+Linkage \* Repression | 9 | 0.6303763 | 0.9022701 |
| within | UN | Comprehensiveness+Linkage+Claims | 1 | 0.7458176 | 0.9207195 |
| within | US | Comprehensiveness \* ~Linkage \* ~Vulnerability+Linkage \* ~Repression \* Claims+~Vulnerability \* Repression \* Claims | 1 | 0.7864914 | 0.5616704 |
| within | US | Comprehensiveness \* ~Vulnerability \* Claims+Linkage \* ~Repression \* Claims+~Vulnerability \* Repression \* Claims | 2 | 0.7866918 | 0.5396568 |
| within | US | Comprehensiveness \* ~Linkage \* ~Vulnerability+Linkage \* ~Vulnerability \* Claims+~Vulnerability \* Repression \* Claims+Comprehensiveness \* Linkage \* Vulnerability \* ~Repression | 3 | 0.7960289 | 0.5710586 |
| within | US | Comprehensiveness \* ~Vulnerability \* Claims+Linkage \* ~Vulnerability \* Claims+~Vulnerability \* Repression \* Claims+Comprehensiveness \* Linkage \* Vulnerability \* ~Repression | 4 | 0.7975690 | 0.5522823 |

### Other packages used in this vignette

Yihui Xie (2021): *knitr: A General-Purpose Package for Dynamic Report
Generation in R.* R package version 1.33.

Yihui Xie (2015): *Dynamic Documents with R and knitr.* 2nd edition.
Chapman and Hall/CRC. ISBN 978-1498716963

Yihui Xie (2014): *knitr: A Comprehensive Tool for Reproducible Research
in R.* In Victoria Stodden, Friedrich Leisch and Roger D. Peng, editors,
Implementing Reproducible Computational Research. Chapman and Hall/CRC.
ISBN 978-1466561595
