% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/upset_conditions.R
\name{upset_conditions}
\alias{upset_conditions}
\title{Aggregation of individual conditions over partition-specific models}
\usage{
upset_conditions(df, nsets)
}
\arguments{
\item{df}{Dataframe created with \code{\link{partition_min}} or
\code{\link{partition_min_inter}}.}

\item{nsets}{Number of sets to include in plot (default is 5).}
}
\value{
An UpSet plot produced with \code{\link[UpSetR]{upset}}.
}
\description{
Models that have been derived for individual partitions are first
decomposed into conditions, that is single conditions or
conditions that are INUS (insufficient conditions that are necessary
parts of a conjunction that is unnecessary and sufficient). The individual
conditions are aggregated using UpSet plots to determine how frequent
they are individually and in combination.
}
\examples{
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

}
