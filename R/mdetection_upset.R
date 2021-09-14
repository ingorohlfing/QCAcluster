#' Internal comparison function for configurations.
#'
#' @importFrom  stringi stri_detect_fixed
#'
#' @param ls List of QCA solutions 
#' @noRd
#'
#' @return A list counting the individual solutions
#' or configurations.
mdetection_upset <- function(ls, x){
  mtr <- NULL
  for(l in x) {
    vctr <- as.numeric(stringi::stri_detect_fixed(ls, l))
    mtr <- cbind(mtr, vctr)
    mtr <- as.data.frame(mtr)
  }

  return(mtr)
}
