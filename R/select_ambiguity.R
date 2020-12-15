#' select.fnc
#' 
#' A function that selects element pertaining
#' to each part of ambiguous solutions
#'
#'
#' @param sol QCA solution
#' @param elmnt is specific part of the solution that 
#' is needed (as saved by QCA package)
#' @return returns coverage and consistency information
#' for each part of the ambiguous solution

select.fnc <- function(sol, elmnt){
  tmp <- map(sol$IC$individual, function(x) x[[elmnt]])
  return(tmp)
}
