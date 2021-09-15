#' Aggregation of individual configurations  over partition-specific models
#'
#' @description 
#' Models that have been derived for individual partitions are first 
#' decomposed into sufficient terms, that is single sufficient conditions or 
#' configurations. The individual terms are aggregated using UpSet plots to 
#' determine how frequent they are individually and in combination.
#'
#' @importFrom magrittr %>%
#' @importFrom stringi stri_trim stri_unique
#' @importFrom purrr map
#' @import UpSetR
#'
#'
#' @param df Dataframe created with \code{\link{partition_min}} or
#' \code{\link{partition_min_inter}}.
#' @param nsets Number of sets to include in plot (default is 5).
#' @md
#'
#' @return An UpSet plot produced with \code{\link[UpSetR]{upset}}.
#' @md
#'
#' @export
upset_conf <- function(df, nsets) {
  

  temp1 <- unlist(df$solution)
  temp1 <- purrr::map(temp1, function(x) stringi::stri_trim(x))
  temp1 <- purrr::map(temp1, function(x) stringi::stri_split_fixed(x, "+"))
  all_values <- stringi::stri_unique(unlist(temp1))
  all_values
  
  finl <- mdetection_upset(temp1, all_values)
  colnames(finl) <- all_values
  UpSetR::upset(finl, order.by = "freq", nsets = nsets)
}
