#' Aggregation of individual configurations  over partition-specific models
#'
#' @description 
#' Models that have been derived for individual partitions are first 
#' decomposed into sufficient terms, that is single sufficient conditions or 
#' configurations. The individual terms are aggregated using UpSet plots to 
#' determine how frequent they are individually and in combination.
#'
#' @importFrom magrittr %>%
#' @importFrom stringi stri_trim stri_unique stri_split_fixed
#' @importFrom purrr map
#' @import UpSetR
#'
#' @param df Dataframe created with \code{\link{partition_min}} or
#' \code{\link{partition_min_inter}}.
#' @param nsets Number of sets to include in plot (default is 5).
#' @md
#'
#' @return An UpSet plot produced with \code{\link[UpSetR]{upset}}.
#' @md
#' 
#' @examples
#' data(Grauvogel2014)
#' GS_cons <- partition_min(
#'  dataset = Grauvogel2014,
#'  units = "Sender",
#'  cond = c("Comprehensiveness", "Linkage", "Vulnerability",
#'           "Repression", "Claims"),
#'  out = "Persistence",
#'  n_cut = 1, incl_cut = 0.75,
#'  solution = "C",
#'  BE_cons = rep(0.75, 3),
#'  BE_ncut = rep(1, 3))
#' upset_conf(GS_pars, nsets = 6)
#'
#' @export
upset_configurations <- function(df, nsets) {
  temp1 <- unlist(df$solution) # is this needed? df$solution is a column of a dataframe
  temp1 <- purrr::map(temp1, function(x) stringi::stri_trim(x))
  all_values <- stringi::stri_unique(unlist(temp1))
  finl <- mdetection_upset(temp1, all_values)
  colnames(finl) <- all_values
  UpSetR::upset(finl, order.by = "freq", nsets = nsets)
}
