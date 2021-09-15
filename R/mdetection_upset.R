#' Internal comparison function for aggregation over configurations.
#'
#' @importFrom  stringi stri_detect_fixed
#'
#' @param ls List of QCA models 
#' @noRd
#'
#' @return A list counting the individual models or configurations.
mdetection_upset <- function(ls, x) {
  mtr <- NULL
  for(l in x) {
    vctr <- as.numeric(stringi::stri_detect_fixed(ls, l))
    mtr <- cbind(mtr, vctr)
    mtr <- as.data.frame(mtr)
  }

  return(mtr)
}
