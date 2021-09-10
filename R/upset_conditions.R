#' Aggregation of individual conditions partition-specific 
#' models
#'  
#' Models that have been derived for individual partitions are first 
#' decomposed into individual conditions, that is single conditions or 
#' conditions that are INUS (insufficient conditions that are necessary
#' parts of a conjunction that is unnecessary and sufficient).
#' The individual conditions are aggregated using UpSet plots to determine
#' how frequent they are.
#'  
#' @importFrom magrittr %>% 
#' @importFrom stringi stri_split_fixed stri_unique
#' @importFrom plyr ldply
#' @importFrom purrr map
#' @import UpSetR
#'
#' @param df Dataframe created with \code{\link{partition_min}} or
#' \code{\link{partition_min_inter}}.
#' @md
#' @param nsets Number of terms to include in plot
#' 
#' @examples 
#' data(Grauvogel2014)
#' GS_pars <- partition_min(
#'  dataset = Grauvogel2014,
#'  units = "Sender",
#'  cond = c("Comprehensiveness", "Linkage", "Vulnerability",
#'           "Repression", "Claims"),
#'  out = "Persistence",
#'  n_cut = 1, incl_cut = 0.75,
#'  solution = "P",
#'  BE_cons = rep(0.75, 3),
#'  BE_ncut = rep(1, 3))
#' upset_plot(GS_pars, nsets = 5)
#' 
#' @export
upset_plot <- function(df, nsets) {

temp1 <- purrr::map(unlist(df$solution), function(x)stringi::stri_split_fixed(x, "*") %>% 
                      unlist())
temp1 <- purrr::map(temp1, function(x)
    stringi::stri_split_fixed(x, "+") %>% unlist())
  all_values <- stringi::stri_unique(unlist(temp1))
  final_matrix <- plyr::ldply(temp1, function(y)
    comparison(x = all_values, y = y, num = T))
  final_matrix$.id <- NULL
  colnames(final_matrix) <- all_values
  UpSetR::upset(final_matrix, order.by = "freq", nsets = nsets)

}