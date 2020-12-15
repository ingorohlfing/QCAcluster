#' ambig
#' 
#' A function that detects if a QCA solution is ambiguous
# 1. ambig: returns information about ambiguity
#'
#'
#' subsets of consistent and inconsistent rows.
#'
#' @param sol QCA solution
#' 
#' @return returns FALSE if the solution is not ambiguous
#' if the solution is ambiguous, the function prints out
#' the length of the ambiguity

ambig <- function(sol){
  
  a <- length(sol$IC$individual)
  
  if (a < 2) {
    print(FALSE)
  } else {
    print(a)
  }
}
