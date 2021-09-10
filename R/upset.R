#'  Aggregation of partition-specific models
#'
#' @param df Dataframe created with \code{\link{partition_min}} or
#' \code{\link{partition_min_inter}}.
#' @param nsets Number of terms to include in plot
#' @md
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