#' UpSetR for configurations
#'
#'
#' @importFrom magrittr %>%
#' @importFrom stringi stri_trim stri_unique
#' @importFrom purrr map
#' @import UpSetR
#'
#'
#' @param df Dataframe created with \code{\link{partition_min}} or
#' \code{\link{partition_min_inter}}.
#' 
#' @param nsets Specifies number of sets to be plotted. Argument
#' imported from the \code{\link{upset}} function from \pkg{UpSetR}.
#'
#' @return A plot presenting the frequency of individual
#' terms and their cooccurrences across QCA solutions.
#'
#' @export
upset_conf <- function(df, nsets) {
  

  temp1 <- unlist(df$solution)
  temp1 <- purrr::map(temp1, function(x) stringi::stri_trim(x))
  temp1 <- purrr::map(temp1, function(x) stringi::stri_split_fixed(x, "+"))
  all_values <- stringi::stri_unique(unlist(temp1))
  all_values
  
  finl <- detection(temp1, all_values)
  colnames(finl) <- all_values
  UpSetR::upset(finl, order.by = "freq", nsets = nsets)
}
