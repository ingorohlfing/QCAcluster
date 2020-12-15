#' select.fnc
#' 
#' A function that selects element pertaining
#' to each part of ambiguous solutions
#' @importFrom purrr map
#'
#' @param sol QCA solution
#' @param elmnt is specific part of the solution that
#' @param overall a logical argument that allows to move between different 
#' levels of ambiguous solution. When FALSE returns elements per each model
#' within an ambiguous solution. When TRUE returns overall values
#' is needed (as saved by QCA package)
#' @return returns coverage and consistency information
#' for each part of the ambiguous solution

select.fnc <- function(sol, elmnt, overall = FALSE){
  
  if (overall == TRUE) {
    tmp <- sol$IC[[elmnt]]
  } else {
    tmp <- purrr::map(sol$IC$individual, function(x) x[[elmnt]])
  }
  return(tmp)
}
