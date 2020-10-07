#' Calculates the weight of partitions for intermediate solution
#'
#' \code{wop} calculates the weight of partitions for each 
#' individual solution derived from the \code{\link{partition_min_inter}}
#' function for pooled solution parameters. It derives a a score for the
#' numerator and denimonator for each solution of the clustered data.
#' 
#' @param x Calibrated pooled dataset for partitioning and minimization
#' @param units Units defining the within-dimension of data (time series)
#' @param time Periods defining the between-dimension of data (cross sections)
#' @param cond Conditions used for the pooled analysis
#' @param out Outcome used for the pooled analysis
#' @param n_cut Frequency cut-off for designating truth table rows as observed
#' @param incl_cut Inclusion cut-off for designating truth table rows as
#' consistent
#' @param intermediate A vector of directional expectations to derive intermediate solutions
#' @param BE_cons Inclusion (or consistency) thresholds for cross sections. 
#' Must be specified as a numeric vector with length equaling the number of
#' cross sections. Numbers correspond to the order of the cross section ID
#' in the data (such as years in ascending order).
#' @param WI_cons Inclusion (or consistency) thresholds for time series. 
#' Must be specified as a numeric vector with length equaling the number of
#' time series. Numbers correspond to the order of the time series (unit) ID
#' in the data (such as countries in alphabetical order).
#' @param BE_ncut The minimum number of cases under which a truth table 
#' row is declared as a remainder for creating an individual truth 
#' table for each cross section.
#' They must be specified as a numeric vector. Its length should be
#' equal the number of cross sections. The order of thresholds corresponds
#' to the order of the cross sections in the data defined by the cross-section
#' ID in the dataset (such as years in ascending order).
#' @param WI_ncut The minimum number of cases under which a truth table
#' row is declared as a remainder for creating an individual truth 
#' table for each time series. 
#' They must be specified as a numeric vector. Its length should be
#' equal the number of time series. The order of thresholds corresponds
#' to the order of the of the time-series (unit) ID
#' in the dataset (such as countries in alphabetical order).
#' 
#' @return A dataframe with the weight of the partitions for pooled consistency
#' scores.
#' 
#' @examples
#' data(schwarz2016)
#' Schwarz_WOP_Inter <- wop_intermediate(schwarz2016, units = "country", time = "year", cond = c("poltrans", "ecotrans", "reform", "conflict", 
#' "attention"), out = "enlarge", 1, 0.8, intermediate = c("1", "1", "1", "1", "1"))
#'
#' @export
wop_intermediate <- function(dataset, units, time, cond, out, n_cut, incl_cut, 
                             intermediate, BE_cons, WI_cons,
                             BE_ncut,WI_ncut) {
  
  # turning of warnings
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 
  #### Splitting the data ####
  x <- dataset
  if (missing(units)) {
    colnames(x)[which(names(x) == time)] <- "time"
    x <- x[with(x, order(time)), ]
    xB <- x
    if (missing(BE_cons)) {
      BE_cons <- rep(incl_cut, times = lengtht)
      xB$consis <- rep(incl_cut, times = nrow(x))
    } else {
      BE_cons <- BE_cons
      xB$consis <- rep(BE_cons, each = lengthu)
    }
    if (missing(BE_ncut)) {
      BE_ncut <- rep(n_cut, times = lengtht)
      xB$ncut <- rep(n_cut, times = nrow(x))
    } else {
      BE_ncut <- BE_ncut
      xB$ncut <- rep(BE_ncut, each = lengthu)
    }
    xxx <- 1
    BE_list <- split(xB, xB[, "time"])
  } else if (missing(time)) {
    colnames(x)[which(names(x) == units)] <- "units"
    x <- x[with(x, order(units)), ]
    xW <- x
    if (missing(WI_cons)) {
      WI_cons <- rep(incl_cut, times = lengthu)
      xW$consis <- rep(incl_cut, times = nrow(x))
    } else {
      WI_cons <- WI_cons
      xW$consis <- rep(WI_cons, times = lengtht)
    }
    if (missing(WI_ncut)) {
      WI_ncut <- rep(n_cut, times = lengthu)
      xW$ncut <- rep(n_cut, times = nrow(x))
    } else {
      WI_ncut <- WI_ncut
      xW$ncut <- rep(WI_ncut, each = lengtht)
    }
    xxx <- 2
    WI_list <- split(xW, xW[, "units"])
    
  } else {
    xxx <- 3
    colnames(x)[which(names(x) == time)] <- "time"
    colnames(x)[which(names(x) == units)] <- "units"
    
    x <- x[with(x, order(time)), ]
    
    lengtht <- length(unique(x$time))
    lengthu <- length(unique(x$units))
    
    xB <- x
    xW <- x
    
    # Assigning individual consistency thresholds if available
    if (missing(BE_cons)) {
      BE_cons <- rep(incl_cut, times = lengtht)
      xB$consis <- rep(incl_cut, times = nrow(x))
    } else {
      BE_cons <- BE_cons
      xB$consis <- rep(BE_cons, each = lengthu)
    }
    if (missing(BE_ncut)) {
      BE_ncut <- rep(n_cut, times = lengtht)
      xB$ncut <- rep(n_cut, times = nrow(x))
    } else {
      BE_ncut <- BE_ncut
      xB$ncut <- rep(BE_ncut, each = lengthu)
    }
    
    if (missing(WI_cons)) {
      WI_cons <- rep(incl_cut, times = lengthu)
      xW$consis <- rep(incl_cut, times = nrow(x))
    } else {
      WI_cons <- WI_cons
      xW$consis <- rep(WI_cons, times = lengtht)
    }
    if (missing(WI_ncut)) {
      WI_ncut <- rep(n_cut, times = lengthu)
      xW$ncut <- rep(n_cut, times = nrow(x))
    } else {
      WI_ncut <- WI_ncut
      xW$ncut <- rep(WI_ncut, each = lengtht)
    }
    
    BE_list <- split(xB, xB[, "time"])
    WI_list <- split(xW, xW[, "units"])
  }
  
  
  x$consis <- incl_cut
  x$ncut <- n_cut
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

      SOL <- "No variation in all coniditions"
      zz <- as.data.frame(SOL)
      zz$model <- "-"
      zz$partition <- part
      zz$type <- type
      zz$denom <- "-"
      zz$num <- "-"
      zz <- zz[!duplicated(zz), ]
      colnames(zz)[1] <- "solution"
      
    } else {
      
      s <- has_error(susu <- try(suppressWarnings(truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)-1][1], n.cut = x[, ncol(x)][1])), silent = TRUE))
      
      if (s == F) {
        x1 <- try(suppressWarnings(truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)-1][1], n.cut = x[, ncol(x)][1])), silent = TRUE)
        
        x2 <- x1$tt$OUT
        x2[x2 == "?"] <- NA
        x2 <- as.numeric(x2)
        # x2[is.na(x2)] <- 0.5
        x2 <- na.omit(x2)
        x2 <- mean(x2)
        
        if (x2 == 0) {
          
          SOL <- "All inconsistent"
          zz <- as.data.frame(SOL)
          zz$model <- "-"
          zz$partition <- part
          zz$type <- type
          zz$denom <- "-"
          zz$num <- "-"
          zz <- zz[!duplicated(zz), ]
          colnames(zz)[1] <- "solution"
          
        } else if (x2 == 1) {
          
          SOL <- "All consistent"
          zz <- as.data.frame(SOL)
          zz$model <- "-"
          zz$partition <- part
          zz$type <- type
          zz$denom <- "-"
          zz$num <- "-"
          zz <- zz[!duplicated(zz), ]
          colnames(zz)[1] <- "solution"
          
        } else {
          t <- has_error(susu <- try(minimize(x1, explain = "1", dir.exp = intermediate, include = "?", details = T, show.cases = T, 
                                              all.sol = T, row.dom = F), silent = TRUE))
          
          if (t == F) {
            
            x3 <- minimize(x1, explain = "1", dir.exp = intermediate, include = "?", details = T, show.cases = T, all.sol = T, 
                           row.dom = F)
            
            a <- x3$i.sol
            ININ <- lapply(a, intersol)
            ININ1 <- lapply(ININ, intersol2)
            zz <- unlist(ININ1)
            zz <- as.data.frame(zz)
            
            pim <- x3$pims
            pimlength <- as.numeric(ncol(pim))
            pim$max <- do.call(pmax, pim[1:pimlength])
            denom <- sum(pim$max)
            
            pim1 <- x3$pims
            pimlength1 <- as.numeric(ncol(pim1))
            pim1$max <- do.call(pmax, pim1[1:pimlength])
            pim1$out <- unlist(x[out])
            pim1$min <- with(pim1, pmin(max, out))
            num <- sum(pim1$min)
            
            zz$denom <- denom
            zz$num <- num
            rownames(zz) <- c()
            colnames(zz)[1] <- "solution"
            zz <- zz[!duplicated(zz), ]
            numberrows <- nrow(zz)
            if (numberrows == 1) {
              zz$model <- 1
            } else {
              zz$model <- as.numeric(rownames(zz))
            }
            zz$partition <- part
            zz$type <- type
            
          } else {
            SOL <- "Values specified in the directional expectations do not appear in the data"
            zz <- as.data.frame(SOL)
            zz$model <- "-"
            zz$partition <- part
            zz$type <- type
            zz$denom <- "-"
            zz$num <- "-"
            zz <- zz[!duplicated(zz), ]
            colnames(zz)[1] <- "solution"
            
          }
        }
      } else {
        
        SOL <- "no combinations at this frequency cutoff"
        zz <- as.data.frame(SOL)
        zz$model <- "-"
        zz$partition <- part
        zz$type <- type
        zz$denom <- "-"
        zz$num <- "-"
        zz <- zz[!duplicated(zz), ]
        colnames(zz)[1] <- "solution"
        
      }
    }
    zz
  }
  
  if (missing(time)) {
    WI_list1 <- quiet(lapply(WI_list, pqmcc))
    PO_list1 <- quiet(lapply(PO_list, pqmcc))
    dff2 <- ldply(WI_list1)[, -1]
    dff3 <- ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    ntot <- as.numeric(mean(dff3$num))
    dtot <- as.numeric(mean(dff3$denom))
    
    total <- rbind(dff3, dff2)
  } else if (missing(units)) {
    BE_list1 <- quiet(lapply(BE_list, pqmcc))
    PO_list1 <- quiet(lapply(PO_list, pqmcc))
    
    dff1 <- ldply(BE_list1)[, -1]
    dff3 <- ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    ntot <- as.numeric(mean(dff3$num))
    dtot <- as.numeric(mean(dff3$denom))
    
    total <- rbind(dff3, dff1)
    
  } else {
    BE_list1 <- quiet(lapply(BE_list, pqmcc))
    WI_list1 <- quiet(lapply(WI_list, pqmcc))
    PO_list1 <- quiet(lapply(PO_list, pqmcc))
    
    dff1 <- ldply(BE_list1)[, -1]
    dff2 <- ldply(WI_list1)[, -1]
    dff3 <- ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    ntot <- as.numeric(mean(dff3$num))
    dtot <- as.numeric(mean(dff3$denom))
    
    total <- rbind(dff3, dff1, dff2)
  }
  total <- total[, c(6, 5, 1, 4, 3, 2)]
  return(total)
}