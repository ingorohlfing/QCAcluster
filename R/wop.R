#' Calculation of weight of partitions in pooled solution parameters
#' for conservative or parsimonious solution
#' 
#' \code{wop} calculates the weight of partitions in the pooled
#' solution parameters (consistecy, coverage) for the conservative
#' and parsimonious solution. 
#'
#' @importFrom plyr ldply
#' @importFrom testit has_error
#' @import QCA
#' 
#' @param x Calibrated pooled dataset for partitioning and minimization
#' @param units Units defining the within-dimension of data (time series)
#' @param time Periods defining the between-dimension of data (cross sections)
#' @param cond Conditions used for the pooled analysis
#' @param out Outcome used for the pooled analysis
#' @param n_cut Frequency cut-off for designating truth table rows as observed
#' @param incl_cut Inclusion cut-off for designating truth table rows as
#' consistent
#' @param solution A character specifying the type of solution that should
#' be derived. "C" produces the conservative (or complex) solution, "P" the
#' parsimonious solution. See \code{\link{wop_intermediate}} for the
#' intermediate solution.
#' @param BE_cons Inclusion (or consistency) thresholds for cross sections. 
#' Must be specified as a numeric vector with length equaling the number of
#' cross sections. Numbers correspond to the order of the cross section ID
#' in the data (such as years in ascending order).
#' @param WI_cons Inclusion (or consistency) thresholds for time series. 
#' Must be specified as a numeric vector with length equaling the number of
#' time series. Numbers correspond to the order of the time series (unit) ID
#' in the data (such as countries in alphabetical order).
#' @param BE_ncut For *cross sections*, the minimum number of members needed
#' for declaring a truth table row 
#' as relevant as opposed to designating it as a remainder.
#' Must be specified as a numeric vector. Its length should be
#' equal the number of cross sections. The order of thresholds corresponds
#' to the order of the cross sections in the data defined by the cross-section
#' ID in the dataset (such as years in ascending order).
#' @param WI_ncut For *time series*, the minimum number of members needed
#' for declaring a truth table row 
#' as relevant as opposed to designating it as a remainder. 
#' Must be specified as a numeric vector. Its length should be
#' equal the number of time series. The order of thresholds corresponds
#' to the order of the of the time-series (unit) ID
#' in the dataset (such as countries in alphabetical order).
#' 
#' @return A dataframe with the weight of the partitions for pooled consistency
#' and coverage scores.
#' 
#' @examples
#' data(Thiem2016)
#' Thiem_wop_pars <- wop(
#'   dataset = Thiem2011,
#'   units = "country", time = "year",
#'   cond = c("fedismfs", "homogtyfs", "powdifffs", "comptvnsfs", "pubsupfs", "ecodpcefs"),
#'   out = "memberfs",
#'   n_cut = 6, incl_cut = 0.8,
#'   solution = "P",
#'   BE_cons = c(0.9, 0.8, 0.7, 0.8, 0.85, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8),
#'   BE_ncut = rep(1, 11),
#'   WI_cons = c(0.75, 0.8, 0.9, 0.8, 0.85, rep(0.75, 10)),
#'   WI_ncut = rep(1, 15))
#' Thiem_wop_pars

#' @export
wop <- function(dataset, units, time, cond, out, n_cut, incl_cut, 
                solution, BE_cons, WI_cons, BE_ncut, WI_ncut) {
  
  # turning of warnings
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 
  
  x1 <- dataset
  x2 <- dataset
  if (!missing(BE_cons) & !missing(time)) {
    colnames(x1)[which(names(x1) == time)] <- "time"
    beconstest <- all(length(unique(x1$time)) == length(BE_cons))  
  } else{
    beconstest <- T
  } 
  if (!missing(WI_cons) & !missing(units)) {
    colnames(x1)[which(names(x1) == units)] <- "units"
    wiconstest <- all(length(unique(x1$units)) == length(WI_cons))  
  } else{
    wiconstest <- T
  } 
  if (!missing(BE_ncut) & !missing(time)) {
    colnames(x2)[which(names(x2) == time)] <- "time"
    bencuttest <- all(length(unique(x2$time)) == length(BE_ncut))  
  } else{
    bencuttest <- T
  } 
  if (!missing(WI_ncut) & !missing(units)) {
    colnames(x2)[which(names(x2) == units)] <- "units"
    wincuttest <- all(length(unique(x2$units)) == length(WI_ncut))  
  } else{
    wincuttest <- T
  } 
  
  if(beconstest == F) {
    stop("The number of BE_cons values does not match the number of periods")
  } else {
    if(wiconstest == F) {
      stop("The number of WI_cons values does not match the number of units")
    } else {
      if(bencuttest == F) {
        stop("The number of BE_ncut values does not match the number of periods")
      } else {  
        if(wincuttest == F) {
          stop("The number of WI_ncut values does not match the number of units")
        } else {
  
  # splitting the data if time and unit values are available
  x <- dataset
  if (missing(units)) {
    colnames(x)[which(names(x) == time)] <- "time"
    xB <- x
    
    if (missing(BE_cons)) {
      xB$consis <- rep(incl_cut, times = nrow(x))
    } else {
      
      time_data <- as.data.frame(unique(x$time))
      names(time_data)[1] <- "time"
      time_data$consis <- BE_cons
      
      xB <- merge(xB, time_data, by = "time", all.x = T)
    }
    if (missing(BE_ncut)) {
      xB$ncut <- rep(n_cut, times = nrow(x))
    } else {
      
      time_data1 <- as.data.frame(unique(x$time))
      names(time_data1)[1] <- "time"
      time_data1$ncut <- BE_ncut
      
      xB <- merge(xB, time_data1, by = "time", all.x = T)
    }
    xxx <- 1
    BE_list <- split(xB, xB[, "time"])
    
  } else if (missing(time)) {
    
    colnames(x)[which(names(x) == units)] <- "units"
    xW <- x
    
    if (missing(WI_cons)) {
      xW$consis <- rep(incl_cut, times = nrow(x))
    } else {
      
      unit_data <- as.data.frame(unique(x$units))
      names(unit_data)[1] <- "units"
      unit_data$consis <- WI_cons
      
      xW <- merge(xW, unit_data, by = "units", all.x = T)
    }
    if (missing(WI_ncut)) {
      xW$ncut <- rep(n_cut, times = nrow(x))
    } else {
      
      unit_data1 <- as.data.frame(unique(x$units))
      names(unit_data1)[1] <- "units"
      unit_data1$ncut <- WI_ncut
      
      xW <- merge(xW, unit_data1, by = "units", all.x = T)
    }
    xxx <- 2
    WI_list <- split(xW, xW[, "units"])
    
  } else {
    xxx <- 3
    colnames(x)[which(names(x) == time)] <- "time"
    colnames(x)[which(names(x) == units)] <- "units"
    
    xB <- x
    xW <- x
    
    # Assigning individual consistency thresholds if available
    if (missing(BE_cons)) {
      xB$consis <- rep(incl_cut, times = nrow(x))
    } else {
      
      time_data <- as.data.frame(unique(x$time))
      names(time_data)[1] <- "time"
      time_data$consis <- BE_cons
      
      xB <- merge(xB, time_data, by = "time", all.x = T)
    }
    if (missing(BE_ncut)) {
      xB$ncut <- rep(n_cut, times = nrow(x))
    } else {
      
      time_data1 <- as.data.frame(unique(x$time))
      names(time_data1)[1] <- "time"
      time_data1$ncut <- BE_ncut
      
      xB <- merge(xB, time_data1, by = "time", all.x = T)
    }
    
    if (missing(WI_cons)) {
      xW$consis <- rep(incl_cut, times = nrow(x))
    } else {
      
      unit_data <- as.data.frame(unique(x$units))
      names(unit_data)[1] <- "units"
      unit_data$consis <- WI_cons
      
      xW <- merge(xW, unit_data, by = "units", all.x = T)
    }
    if (missing(WI_ncut)) {
      xW$ncut <- rep(n_cut, times = nrow(x))
    } else {
      
      unit_data1 <- as.data.frame(unique(x$units))
      names(unit_data1)[1] <- "units"
      unit_data1$ncut <- WI_ncut
      
      xW <- merge(xW, unit_data1, by = "units", all.x = T)
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
      
      s <- testit::has_error(susu <- try(suppressWarnings(QCA::truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)-1][1], n.cut = x[, ncol(x)][1])), silent = TRUE))
      
      if (s == F) {
        x1 <- try(suppressWarnings(QCA::truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)-1][1], n.cut = x[, ncol(x)][1])), silent = TRUE)
        
        
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
          
        } else if (x2 == 1 & solution == "P") {
          
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
          
          if (solution == "C") {
            
            x3 <- QCA::minimize(x1, explain = "1", include = "1", details = T, show.cases = T, all.sol = T, row.dom = F)
            
          } else if (solution == "P") {
            
            x3 <- QCA::minimize(x1, explain = "1", include = "?", details = T, show.cases = T, all.sol = T, row.dom = F)
            
          } else {
            
            x3 <- QCA::minimize(x1, explain = "1", include = "1", details = T, show.cases = T, all.sol = T, row.dom = F)
          }
          
          SOL <- x3$solution[]
          tete <- list(cons = SOL)
          neux <- lapply(tete$cons, paster)
          neuxx <- unlist(neux)
          zz <- as.data.frame(neuxx)
          
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
          numberrows <- nrow(zz)
          if (numberrows == 1) {
            zz$model <- 1
          } else {
            zz$model <- as.numeric(rownames(zz))
          }
          zz$partition <- part
          zz$type <- type
          colnames(zz)[1] <- "solution"
          
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
  
  #### Application of Function ####
  
  if (missing(time)) {
    WI_list1 <- quiet(lapply(WI_list, pqmcc))
    PO_list1 <- quiet(lapply(PO_list, pqmcc))
    dff2 <- plyr::ldply(WI_list1)[, -1]
    dff3 <- plyr::ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    ntot <- as.numeric(mean(dff3$num))
    dtot <- as.numeric(mean(dff3$denom))
    
    total <- rbind(dff3, dff2)
  } else if (missing(units)) {
    BE_list1 <- quiet(lapply(BE_list, pqmcc))
    PO_list1 <- quiet(lapply(PO_list, pqmcc))
    
    dff1 <- plyr::ldply(BE_list1)[, -1]
    dff3 <- plyr::ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    ntot <- as.numeric(mean(dff3$num))
    dtot <- as.numeric(mean(dff3$denom))
    
    total <- rbind(dff3, dff1)
    
  } else {
    BE_list1 <- quiet(lapply(BE_list, pqmcc))
    WI_list1 <- quiet(lapply(WI_list, pqmcc))
    PO_list1 <- quiet(lapply(PO_list, pqmcc))
    
    dff1 <- plyr::ldply(BE_list1)[, -1]
    dff2 <- plyr::ldply(WI_list1)[, -1]
    dff3 <- plyr::ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    ntot <- as.numeric(mean(dff3$num))
    dtot <- as.numeric(mean(dff3$denom))
    
    total <- rbind(dff3, dff1, dff2)
  }
  total <- total[, c(6, 5, 1, 4, 3, 2)]
  total$solution <- gsub("\\*"," * ",total$solution)
  return(total)}}}}
}
