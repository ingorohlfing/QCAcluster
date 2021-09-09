# UpSetR
# Documentation required
# Input: data frame that has a column with solutions
# Output: an Upset Plot
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