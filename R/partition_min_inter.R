#' Title
#'
#' @param x 
#' @param units 
#' @param time 
#' @param cond 
#' @param out 
#' @param n_cut 
#' @param incl_cut 
#' @param intermediate 
#' @param BE_cons 
#' @param WI_cons 
#'
#' @return
#' @export
#'
#' @examples
partition_min_inter <- function(x, units, time, cond, out, n_cut, incl_cut, intermediate, BE_cons, WI_cons) {
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
  
  
  intersol <- function(x) {
    neux <- lapply(x$solution, paster)
    neuxx <- unlist(neux)
    zz <- as.data.frame(neuxx)
    zz
  }
  intercons1 <- function(x) {
    
    zz <- x$IC$sol.incl.cov[1]
    
    if (is.null(zz)) {
      BSP2 <- x$IC
      BSP3 <- BSP2$individual
      neu <- sapply(BSP3, function(x) x[2])
      zz <- sapply(neu, function(x) x[1])
    }
    zz <- as.data.frame(zz)
    zz
  }
  
  intercov1 <- function(x) {
    zz <- x$IC$sol.incl.cov[3]
    if (is.null(zz)) {
      BSP2 <- x$IC
      BSP3 <- BSP2$individual
      neu <- sapply(BSP3, function(x) x[2])
      zz <- sapply(neu, function(x) x[3])
    }
    zz <- as.data.frame(zz)
    zz
  }
  
  intersol2 <- function(x) {
    solu <- as.data.frame(x[1])
    solu <- unlist(solu)
    solu
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
      
      x$CONS <- "-"
      zz <- as.data.frame(x$CONS)
      zz$coverage <- "-"
      zz$solution <- "No variation in all coniditions"
      zz$model <- "-"
      zz$partition <- part
      zz$type <- type
      zz <- zz[!duplicated(zz), ]
      colnames(zz)[1] <- "consistency"
      
    } else {
      
      s <- has_error(susu <- try(truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)][1], n.cut = n_cut), silent = TRUE))
      
      if (s == F) {
        x1 <- try(truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)][1], n.cut = n_cut), silent = TRUE)
        
        x2 <- x1$tt$OUT
        x2[x2 == "?"] <- NA
        x2 <- as.numeric(x2)
        # x2[is.na(x2)] <- 0.5
        x2 <- na.omit(x2)
        x2 <- mean(x2)
        
        if (x2 == 0) {
          x$CONS <- "-"
          zz <- as.data.frame(x$CONS)
          zz$coverage <- "-"
          zz$solution <- "All inconsistent"
          zz$model <- "-"
          zz$partition <- part
          zz$type <- type
          zz <- zz[!duplicated(zz), ]
          colnames(zz)[1] <- "consistency"
          
        } else if (x2 == 1) {
          
          x$consistency <- "-"
          zz <- as.data.frame(x$consistency)
          zz$coverage <- "-"
          zz$solution <- "All consistent"
          zz$model <- "-"
          zz$partition <- part
          zz$type <- type
          zz <- zz[!duplicated(zz), ]
          colnames(zz)[1] <- "consistency"
          
        } else {
          t <- has_error(susu <- try(minimize(x1, explain = "1", dir.exp = intermediate, include = "?", details = T, show.cases = T, 
                                              all.sol = T, row.dom = F), silent = TRUE))
          
          if (t == F) {
            
            x <- minimize(x1, explain = "1", dir.exp = intermediate, include = "?", details = T, show.cases = T, all.sol = T, 
                          row.dom = F)
            
            
            a <- x$i.sol
            ININ <- lapply(a, intersol)
            ININ1 <- lapply(ININ, intersol2)
            zz <- unlist(ININ1)
            zz <- as.data.frame(zz)
            ININCONS <- lapply(a, intercons1)
            ININCOV <- lapply(a, intercov1)
            zzcons <- unlist(ININCONS)
            zzcov <- unlist(ININCOV)
            zz$consistency <- zzcons
            zz$coverage <- zzcov
            zz <- zz[!duplicated(zz), ]
            
            rownames(zz) <- c()
            colnames(zz)[1] <- "solution"
            rownames(zz) <- c()
            zz$model <- as.numeric(rownames(zz))
            zz$partition <- part
            zz$type <- type
            zz <- zz[, c(2, 3, 1, 4, 5, 6)]
          } else {
            x$CONS <- "-"
            zz <- as.data.frame(x$CONS)
            zz$coverage <- "-"
            zz$solution <- "Values specified in the directional expectations do not appear in the data"
            zz$model <- "-"
            zz$partition <- part
            zz$type <- type
            zz <- zz[!duplicated(zz), ]
            colnames(zz)[1] <- "consistency"
          }
        }
        
      } else {
        x$CONS <- "-"
        zz <- as.data.frame(x$CONS)
        zz$coverage <- "-"
        zz$solution <- "no combinations at this frequency cutoff"
        zz$model <- "-"
        zz$partition <- part
        zz$type <- type
        zz <- zz[!duplicated(zz), ]
        colnames(zz)[1] <- "consistency"
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
  
  
  #### Rounding ####
  total$consistency[total$model == "-"] <- NA
  total$coverage[total$model == "-"] <- NA
  total$consistency <- as.numeric(total$consistency)
  total$coverage <- as.numeric(total$coverage)
  
  return(total)
  
}
