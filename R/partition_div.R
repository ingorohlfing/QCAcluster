#' Title
#'
#' @param x 
#' @param units 
#' @param time 
#' @param cond 
#' @param out 
#' @param n_cut 
#' @param incl_cut 
#' @param BE_cons 
#' @param WI_cons 
#'
#' @return
#' @export
#'
#' @examples
diversity <- function(x, units, time, cond, out, n_cut, incl_cut, BE_cons, WI_cons) {
  
  #### Splitting the data ####
  options(warn = -1)
  if (missing(units)) {
    colnames(x)[which(names(x) == time)] <- "time"
    x <- x[with(x, order(time)), ]
    xB <- x
    xB$consis <- rep(incl_cut, times = nrow(x))
    BE_list <- split(xB, xB[, "time"])
    xxx <- 1
  } else if (missing(time)) {
    colnames(x)[which(names(x) == units)] <- "units"
    x <- x[with(x, order(units)), ]
    xW <- x
    xW$consis <- rep(incl_cut, times = nrow(x))
    WI_list <- split(xW, xW[, "units"])
    xxx <- 2
    
  } else {
    xxx <- 3
    colnames(x)[which(names(x) == time)] <- "time"
    colnames(x)[which(names(x) == units)] <- "units"
    
    x <- x[with(x, order(time)), ]
    
    lengtht <- length(unique(x$time))
    lengthu <- length(unique(x$units))
    
    xB <- x
    xW <- x
    
    
    if (missing(BE_cons)) {
      BE_cons <- rep(incl_cut, times = lengtht)
      xB$consis <- rep(incl_cut, times = nrow(x))
      
    } else {
      BE_cons <- BE_cons
      xB$consis <- rep(BE_cons, each = lengthu)
    }
    
    if (missing(WI_cons)) {
      WI_cons <- rep(incl_cut, times = lengthu)
      xW$consis <- rep(incl_cut, times = nrow(x))
    } else {
      WI_cons <- WI_cons
      xW$consis <- rep(WI_cons, times = lengtht)
    }
    
    
    BE_list <- split(xB, xB[, "time"])
    WI_list <- split(xW, xW[, "units"])
    
  }
  x$consis <- incl_cut
  PO_list <- list(x)
  
  paster <- function(x) {
    x <- paste(x, collapse = "+")
    x
  }
  
  #### Function for Between and Within Solutions ####
  pqmcc <- function(x) {
    
    if (xxx == 1) {
      part <- as.character(x$time[1])
      type <- "between"
    } else if (xxx == 2) {
      part <- as.character(x$units[1])
      type <- "within"
    } else {
      partition <- unlist(x$time)
      
      if (partition[1] == partition[2]) {
        part <- as.character(x$time[1])
        type <- "between"
      } else {
        part <- as.character(x$units[1])
        type <- "within"
      }
    }
    check <- x[cond]
    check[check < 0.5] <- 0
    check[check > 0.5] <- 1
    check2 <- as.data.frame(colMeans(check))
    check2[check2 == 1] <- 0
    check3 <- as.numeric(colMeans(check2))
    
    if (check3 == 0) {
      
      zz <- as.data.frame(part)
      zz$type <- type
      zz$diversity <- "No variation in all coniditions"
      zz$diversity_per <- "-"
      zz <- zz[!duplicated(zz), ]
      colnames(zz)[1] <- "partition"
      
    } else {
      
      
      # s <- has_error(truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[,ncol(x)][1], n.cut = n_cut))
      s <- has_error(susu <- try(truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)][1], n.cut = n_cut), silent = TRUE))
      
      if (s == F) {
        x1 <- try(truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)][1], n.cut = n_cut), silent = TRUE)
        
        zz <- as.data.frame(part)
        zz$type <- type
        zz$diversity <- as.numeric(length(x1$indexes))
        zz$diversity_1 <- as.numeric(sum(x1$tt$OUT == 1))
        zz$diversity_0 <- as.numeric(sum(x1$tt$OUT == 0))
        zz$diversity_per <- "???"
        zz <- zz[!duplicated(zz), ]
        colnames(zz)[1] <- "partition"
        
      } else {
        
        zz <- as.data.frame(part)
        zz$type <- type
        zz$DIV <- "no combinations at this frequency cutoff"
        zz$DIV_1 <- "-"
        zz$DIV_0 <- "-"
        zz$DIV_per <- "-"
        zz <- zz[!duplicated(zz), ]
        colnames(zz)[1] <- "part"
        
      }
      
    }
    zz
  }
  
  
  #### Application of Function ####
  
  if (missing(time)) {
    WI_list1 <- lapply(WI_list, pqmcc)
    PO_list1 <- lapply(PO_list, pqmcc)
    dff2 <- ldply(WI_list1)[, -1]
    dff3 <- ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    
    total <- rbind(dff3, dff2)
  } else if (missing(units)) {
    BE_list1 <- lapply(BE_list, pqmcc)
    PO_list1 <- lapply(PO_list, pqmcc)
    
    dff1 <- ldply(BE_list1)[, -1]
    dff3 <- ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    
    total <- rbind(dff3, dff1)
    
  } else {
    BE_list1 <- lapply(BE_list, pqmcc)
    WI_list1 <- lapply(WI_list, pqmcc)
    PO_list1 <- lapply(PO_list, pqmcc)
    
    dff1 <- ldply(BE_list1)[, -1]
    dff2 <- ldply(WI_list1)[, -1]
    dff3 <- ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    
    total <- rbind(dff3, dff1, dff2)
  }
  
  total$diversity_old <- total$diversity
  total$diversity[total$diversity == "No variation in all coniditions"] <- NA
  total$diversity[total$diversity == "no combinations at this frequency cutoff"] <- NA
  total$diversity_1[total$diversity_1 == "-"] <- NA
  total$diversity_0[total$diversity_0 == "-"] <- NA
  total$diversity <- as.numeric(total$diversity)
  total$diversity_1 <- as.numeric(total$diversity_1)
  total$diversity_0 <- as.numeric(total$diversity_0)
  
  y <- as.numeric(max(total$diversity, na.rm = T))
  
  total$diversity_per <- ifelse(is.na(total$diversity), NA, total$diversity/y)
  total$diversity_per_1 <- ifelse(is.na(total$diversity_1), NA, total$diversity_1/y)
  total$diversity_per_0 <- ifelse(is.na(total$diversity_0), NA, total$diversity_0/y)
  total$diversity <- total$diversity_old
  total$diversity_old <- NULL
  return(total)
  
}
