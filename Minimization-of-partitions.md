---
title: "Minimization of partitions"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Minimization of partitions}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(QCAcluster)
```

Two functions allow empirical researchers to partition clustered data on one or two dimensions and to derive solutions for the pooled data and for each partition. `partition_min()` is available for producing conservative or parsimonious solutions. `partition_min_inter()` has to be used for the intermediate solution. (For programming purposes, we opted for separate functions for the intermediate solution and the other two.)

## Panel data: Minimization for cross sections and time series
We first illustrate how one can dempose panel data on two dimensions. In a *between-unit* perspective, the panel is partitioned into multiple cross sections with the `time` argument that specifies the cross section ID (here: years). In a *within-unit* perspective, the data is decomposed into multiple time series with the `units` argument that specifies the time series ID (here: countries). In short, the other parameters are:

- `n_cut`: Frequency threshold for pooled data
- `incl_cut`: Inclusion threshold (a.k.a. consistency threshold) for pooled data
- `solution`: Either `C` for conservative solution (a.k.a. complex solution) or `P` for parsimonious solution
- `BE_cons` and `WI_cons`: Inclusion thresholds for the cross sections and time series. The length of the numeric vector should equal the number of units and time series.
- `BE_ncut` and `WI_ncut`: Frequency thresholds for the cross sections and time series. The length of the numeric vector should equal the number of units and time series.

### Conservative or parsimonious solution
We first illustrate the parsimonious solution with dataset from [Thiem (2011)](https://doi.org/10.1017/S1755773910000251).

```r
# loading panel data (see data description for emails)
data(Thiem2011)
# partitioning data by countries (within-unit) and years (between-unit)
Thiem_pars <- partition_min(
  dataset = Thiem2011,
  units = "country", time = "year",
  cond = c("fedismfs", "homogtyfs", "powdifffs", "comptvnsfs", "pubsupfs", "ecodpcefs"),
  out = "memberfs",
  n_cut = 6, incl_cut = 0.8,
  solution = "P",
  BE_cons = c(0.9, 0.8, 0.7, 0.8, 0.85, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8),
  BE_ncut = rep(1, 11),
  WI_cons = c(0.75, 0.8, 0.9, 0.8, 0.85, rep(0.75, 10)),
  WI_ncut = rep(1, 15))
Thiem_pars
#>       type partition                                                    solution model consistency  coverage
#> 1   pooled         -                                comptvnsfs+fedismfs*pubsupfs     1   0.8976935 0.7113797
#> 2   pooled         -                               comptvnsfs+fedismfs*ecodpcefs     2   0.8949502 0.7158019
#> 3   pooled         -                               comptvnsfs+homogtyfs*pubsupfs     3   0.8780259 0.7342767
#> 4  between      1996                                         fedismfs*comptvnsfs     1   0.9030303 0.3748428
#> 5  between      1996                                         comptvnsfs*pubsupfs     2   0.9885057 0.4327044
#> 6  between      1997                                                  ~powdifffs     1   0.9064748 0.6339623
#> 7  between      1997                                                  comptvnsfs     2   0.8910675 0.5144654
#> 8  between      1997                                         pubsupfs*~ecodpcefs     3   0.8672769 0.4767296
#> 9  between      1998                                                  comptvnsfs     1   0.9288703 0.6090535
#> 10 between      1999                               ~powdifffs+fedismfs*ecodpcefs     1   0.8876404 0.7623643
#> 11 between      1999 ~powdifffs+fedismfs*~homogtyfs+homogtyfs*pubsupfs*ecodpcefs     2   0.8961039 0.7490953
#> 12 between      2000                                comptvnsfs+fedismfs*pubsupfs     1   0.9684685 0.6508577
#> 13 between      2000                               comptvnsfs+fedismfs*ecodpcefs     2   0.9417476 0.6851665
#> 14 between      2000           comptvnsfs+fedismfs*~homogtyfs+homogtyfs*pubsupfs     3   0.9708333 0.7053481
#> 15 between      2001                                         fedismfs+comptvnsfs     1   0.9028436 0.7689203
#> 16 between      2002                                fedismfs+~powdifffs+pubsupfs     1   0.8149780 0.7467205
#> 17 between      2002                                fedismfs+comptvnsfs+pubsupfs     2   0.8214665 0.7800202
#> 18 between      2003                                         pubsupfs+~ecodpcefs     1   0.7985213 0.8529121
#> 19 between      2004                                         fedismfs+~ecodpcefs     1   0.9184290 0.8260870
#> 20 between      2004                                         pubsupfs+~ecodpcefs     2   0.9081726 0.8958333
#> 21 between      2005                                         pubsupfs+~ecodpcefs     1   0.9002695 0.9076087
#> 22 between      2005                              fedismfs+~homogtyfs+~ecodpcefs     2   0.8868101 0.8586957
#> 23 between      2006                                        comptvnsfs+~pubsupfs     1   0.8982118 0.7829736
#> 24 between      2006                               ~pubsupfs+fedismfs*~ecodpcefs     2   0.8335725 0.6966427
#> 25  within        AT                           All truth table rows inconsistent     -          NA        NA
#> 26  within        BE                              No variation in all conditions     -          NA        NA
#> 27  within        DE                             All truth table rows consistent     -          NA        NA
#> 28  within        DK                                                   ~pubsupfs     1   0.8297389 0.9798928
#> 29  within        DK                                                  ~ecodpcefs     2   0.9469154 0.8847185
#> 30  within        ES                             All truth table rows consistent     -          NA        NA
#> 31  within        FI                              No variation in all conditions     -          NA        NA
#> 32  within        FR                             All truth table rows consistent     -          NA        NA
#> 33  within        GR                           All truth table rows inconsistent     -          NA        NA
#> 34  within        IE                           All truth table rows inconsistent     -          NA        NA
#> 35  within        IT                              No variation in all conditions     -          NA        NA
#> 36  within        LU                                                   homogtyfs     1   0.7629630 0.8131579
#> 37  within        NL                             All truth table rows consistent     -          NA        NA
#> 38  within        PT                           All truth table rows inconsistent     -          NA        NA
#> 39  within        SE                           All truth table rows inconsistent     -          NA        NA
#> 40  within        UK                             All truth table rows consistent     -          NA        NA
```

The output of `partition_min()` is a dataframe summarizing the solutions for the pooled data and the partitions and the consistency and coverage values for the solution. There are different reasons why one might not be able to derive a partition-specific solution: All rows could be consistent; all rows could be inconsistent; there is no variation and all cases belong to the same truth table row. If one the reason applies, it is specifically listed in the column `solution`.

The function produces the information than one sees in the dataframe. It can serve as a basis for spotting interesting insights such as partitions with no variation in conditions that can be explored further in a manual analysis of the data.

### Intermediate solution
The intermediate solution is derived with `partition_min_inter()`. (We decided to have a separate function for the intermediate solution for computational purposes.) For the intermediate solution, one has to specify the directional expectations. 

```r
data(Schwarz2016)
Schwarz_inter_1 <- partition_min_inter(
  Schwarz2016, 
  units = "country", time = "year", 
  cond = c("poltrans", "ecotrans", "reform", "conflict", "attention"), 
  out = "enlarge", 
  n_cut = 1, incl_cut = 0.8, 
  intermediate = c("1", "1", "1", "1", "1"))
#> Error in has_error(susu <- try(suppressWarnings(truthTable(x, outcome = out, : konnte Funktion "has_error" nicht finden
```



## Multilevel data
Clustered data can be partitioned on a single dimension if the second dimension is not of interest or if there is only one dimension such as an in multilevel data where lower-level units are nested in higher-level units. We use the dataset by [Grauvogel and von Soest (2014)](https://doi.org/10.1017/S1755773910000251) for illustrating the analysis of multilevel data. The study analyzes the effect of sanctions on authoritarian regimes. The data distinguishes between the source of the sanction (`Sender`) and the target country (`Target). All sanctions have been imposed by the EU, UN or US, which means that target countries are nested in three different senders. We partition the data on the dimension of senders to see how solutions differ across senders.


```r
# loading multilevel data (see data description for emails)
data(Grauvogel2014)
# partitioning data by sender country (higher-level unit)
GS_pars <- partition_min(
  dataset = Grauvogel2014,
  units = "Sender",
  cond = c("Comprehensiveness", "Linkage", "Vulnerability","Repression", "Claims"),
  out = "Persistence",
  n_cut = 1, incl_cut = 0.75,
  solution = "P",
  BE_cons = rep(0.75, 3),
  BE_ncut = rep(1, 3))
GS_pars
#>      type partition
#> 1  pooled         -
#> 2  pooled         -
#> 3  within        EU
#> 4  within        EU
#> 5  within        EU
#> 6  within        EU
#> 7  within        EU
#> 8  within        EU
#> 9  within        EU
#> 10 within        EU
#> 11 within        EU
#> 12 within        UN
#> 13 within        US
#> 14 within        US
#> 15 within        US
#> 16 within        US
#>                                                                                                                                                        solution
#> 1              ~Comprehensiveness*Claims+~Linkage*Claims+~Repression*Claims+~Comprehensiveness*~Linkage*~Repression+Comprehensiveness*~Vulnerability*Repression
#> 2                        ~Comprehensiveness*Claims+~Linkage*Claims+~Repression*Claims+~Comprehensiveness*~Linkage*~Repression+Linkage*~Vulnerability*Repression
#> 3                                                                                                   Vulnerability+~Comprehensiveness*~Repression+Linkage*Claims
#> 4                                                                                                             Vulnerability+~Linkage*~Repression+Linkage*Claims
#> 5                                                                                                           Vulnerability+Linkage*Repression+~Repression*Claims
#> 6                                                                                                               Vulnerability+Linkage*Claims+~Repression*Claims
#> 7                                                              Vulnerability+~Comprehensiveness*Linkage+~Comprehensiveness*~Repression+Comprehensiveness*Claims
#> 8                                                                      Vulnerability+~Comprehensiveness*Linkage+Comprehensiveness*Repression+~Repression*Claims
#> 9                                                                        Vulnerability+~Comprehensiveness*Linkage+Comprehensiveness*Claims+~Linkage*~Repression
#> 10                                                                         Vulnerability+~Comprehensiveness*Linkage+Comprehensiveness*Claims+~Repression*Claims
#> 11                                                                     Vulnerability+~Comprehensiveness*~Repression+Comprehensiveness*Claims+Linkage*Repression
#> 12                                                                                                                             Comprehensiveness+Linkage+Claims
#> 13                                                        Comprehensiveness*~Linkage*~Vulnerability+Linkage*~Repression*Claims+~Vulnerability*Repression*Claims
#> 14                                                          Comprehensiveness*~Vulnerability*Claims+Linkage*~Repression*Claims+~Vulnerability*Repression*Claims
#> 15 Comprehensiveness*~Linkage*~Vulnerability+Linkage*~Vulnerability*Claims+~Vulnerability*Repression*Claims+Comprehensiveness*Linkage*Vulnerability*~Repression
#> 16   Comprehensiveness*~Vulnerability*Claims+Linkage*~Vulnerability*Claims+~Vulnerability*Repression*Claims+Comprehensiveness*Linkage*Vulnerability*~Repression
#>    model consistency  coverage
#> 1      1   0.7758164 0.7336208
#> 2      2   0.7776948 0.7245792
#> 3      1   0.6293355 0.9145825
#> 4      2   0.6327684 0.9049634
#> 5      3   0.6310549 0.9022701
#> 6      4   0.6320277 0.9126587
#> 7      5   0.6273610 0.8945748
#> 8      6   0.6238361 0.9022701
#> 9      7   0.6253391 0.8868796
#> 10     8   0.6261682 0.9022701
#> 11     9   0.6303763 0.9022701
#> 12     1   0.7458176 0.9207195
#> 13     1   0.7864914 0.5616704
#> 14     2   0.7866918 0.5396568
#> 15     3   0.7960289 0.5710586
#> 16     4   0.7975690 0.5522823
```

