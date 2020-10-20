#' Diversity of cases belonging to the same partition of the pooled data
#' 
#' \code{diversity} calculates the diversity of cases that belong to the same
#' partition of the clustered data (time series, cross section etc.).
#' Diversity is measured by the number of truth table rows that the cases
#' cover. It calculates diversity across all truth table rows and for the
#' subsets of consistent and inconsistent rows.
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
#' @return A dataframe presenting the diversity of cases belonging to the
#' same partition. 
#'
#' @examples
#' data(schwarz2016)
#' Schwarz_diversity <- diversity(schwarz2016, units = "country", time = "year", 
#' cond = c("poltrans", "ecotrans", "reform", "conflict", "attention"), 
#' out = "enlarge", 1, 0.8)
#' 
#' @export
diversity <- function(dataset, 
                      units, time, 
                      cond, out, 
                      n_cut, incl_cut, 
                      BE_cons, WI_cons, 
                      BE_ncut, WI_ncut) {
  
  # turning of warnings
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 
  
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
  #x <- BE_list[1]
  #x <- as.data.frame(BE_list[1])
  #x <- Thiem2011[Thiem2011$time == 1996,]
  #x$consis <- 0.8
  #x$ncut <- 2
  
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
      s <- testit::has_error(susu <- try(suppressWarnings(QCA::truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)-1][1], n.cut = x[, ncol(x)][1])), silent = TRUE))
      
      if (s == F) {
        x1 <- try(suppressWarnings(QCA::truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)-1][1], n.cut = x[, ncol(x)][1])), silent = TRUE)
        
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
        zz$diversity <- "no combinations at this frequency cutoff"
        zz$diversity_per <- "-"
        zz <- zz[!duplicated(zz), ]
        colnames(zz)[1] <- "partition"
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
    
    total <- rbind(dff3, dff2)
  } else if (missing(units)) {
    BE_list1 <- quiet(lapply(BE_list, pqmcc))
    PO_list1 <- quiet(lapply(PO_list, pqmcc))
    
    dff1 <- plyr::ldply(BE_list1)[, -1]
    dff3 <- plyr::ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    
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
  total <- total[, c(2, 1, 3, 4, 5, 6, 7, 8)]
  return(total)
}
