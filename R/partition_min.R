#' Generation of conservative or parsimonious solution for individual 
#' partitions
#' 
#' \code{partition_min} decomposes clustered data into individual
#' partitions. For panel data, for example, these can be cross sections,
#' time series or both. The function derives an individual solution for
#' each partition and the pooled data to assess the robustness of the 
#' solutions in a comparative perspective.
#' 
#' @importFrom plyr ldply
#' @importFrom testit has_error
#' @importFrom stats na.omit
#' @import QCA
#'
#' @param dataset Calibrated pooled dataset that is partitioned and minimized for
#' deriving the pooled solution.
#' @param units Units defining the within-dimension of data (time series). If no
#' units are specified, the data is assumed to lack a dimension and be
#' hierarchical.
#' @param time Periods defining the between-dimension of data (cross sections). 
#' This should be specified because it does not make sense to partition a 
#' time series into individual data points.
#' @param cond Conditions used for minimization
#' @param out Outcome used for minimization
#' @param n_cut Frequency cut-off for designating truth table rows as observed
#' as opposed to designating them as remainders for the *pooled* data.
#' @param incl_cut Inclusion (a.k.a. consistency) cut-off for designating 
#' truth table rows as consistent for the *pooled* data.
#' @param solution A character specifying the type of solution that should
#' be derived. \code{C} produces the conservative (or complex) solution, 
#' \code{P} for the parsimonious solution. See \code{\link{partition_min_inter}} 
#' for a separate function for the intermediate solution.
#' @param BE_cons Inclusion thresholds for creating an individual truth table
#' for each cross section.
#' They must be specified as a numeric vector. Its length should be
#' equal the number of cross sections. The order of thresholds corresponds
#' to the order of the cross sections in the data defined by the cross-section
#' ID in the dataset (such as years in ascending order).
#' @param WI_cons Inclusion thresholds for creating an individual truth table
#' for each time series. 
#' They must be specified as a numeric vector. Its length should be
#' equal the number of time series. The order of thresholds corresponds
#' to the order of the of the time-series (unit) ID
#' in the dataset (such as countries in alphabetical order).
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
#' @return A dataframe summarizing the partition-specific and pooled solutions
#' with the following columns:
#' 
#' * \code{type}: The type of the partition. \code{pooled} are rows with information
#' on the pooled data; \code{between} is for cross-section partitions;
#' \code{within} is for time-series partitions.
#' * \code{partition}: Specific dimension of the partition at hand. For 
#' between-dimension, the unit identifiers are included here  (argument \code{units}).
#' For the within-dimension, the time identifier are listed (argument \code{time}).
#' The entry is \code{-} for the pooled data without partitions.
#' * \code{solution}: The solution derived for the partition or the pooled data.
#' Absence of a condition is denoted by the \code{~} sign.
#' * \code{model}: Running ID for models. In the presence of model ambiguity, each
#' model has its own row with its individual solution and parameters. The rest of
#' the information in the row is duplicated, for example by having two rows for
#' the within-partition 1996. The column \code{model} highlights the presence of
#' model ambiguity by numbering all models belonging to the same solution. For 
#' example, if three consecutive rows are numbered 1, 2 and 3, then these rows
#' belong to the same solution and represent model ambiguity. If a 1 in a row
#' is followed by another 1, then there is no model ambiguity.
#' * \code{consistency}: The consistency score (a.k.a. inclusion score) 
#' for the partition of the data or the pooled data.
#' * \code{coverage}: The coverage score for the partition of the data 
#' or the pooled data.
#' @md
#'
#' @examples
#' # loading data from Thiem (EPSR, 2011; see data documentation)
#' data(Thiem2011)
#' 
#' # running function for parsimonious solution
#' Thiem_pars <- partition_min(
#'   dataset = Thiem2011,
#'   units = "country", time = "year",
#'   cond = c("fedismfs", "homogtyfs", "powdifffs", "comptvnsfs", "pubsupfs", "ecodpcefs"),
#'   out = "memberfs",
#'   n_cut = 1, incl_cut = 0.8,
#'   solution = "P",
#'   BE_cons = c(0.9, 0.8, 0.7, 0.8, 0.6, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8),
#'   WI_cons = c(0.5, 0.8, 0.7, 0.8, 0.6, rep(0.8, 10)))
#' 
#' @export
partition_min <- function(dataset, 
                          units, time, 
                          cond, out, 
                          n_cut, incl_cut, 
                          solution, 
                          BE_cons, WI_cons,
                          BE_ncut, WI_ncut) {
  
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
  
  # Creating a function for the transformation of QCA results
  paster <- function(x) {
    x <- paste(x, collapse = "+")
    x
  }
  
  #### Function for between and within solutions ####
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
      x$consistency <- "-"
      zz <- as.data.frame(x$consistency)
      zz$coverage <- "-"
      zz$solution <- "No variation in all conditions"
      zz$model <- "-"
      zz$partition <- part
      zz$type <- type
      zz <- zz[!duplicated(zz), ]
      colnames(zz)[1] <- "consistency"
    } else {
      
      s <- testit::has_error(susu <- try(suppressWarnings(QCA::truthTable(
        x, 
        outcome = out, 
        conditions = cond, 
        incl.cut1 = x[, ncol(x)-1][1],
        n.cut = x[, ncol(x)][1])), silent = TRUE))
      
      if (s == F) {
        x1 <- try(suppressWarnings(QCA::truthTable(x, outcome = out, 
                                                   conditions = cond, 
                                                   incl.cut1 = x[, ncol(x)-1][1], 
                                                   n.cut = x[, ncol(x)][1])), silent = TRUE)
        
        x2 <- x1$tt$OUT
        x2[x2 == "?"] <- NA
        x2 <- as.numeric(x2)
        # x2[is.na(x2)] <- 0.5
        x2 <- na.omit(x2)
        x2 <- mean(x2)
        
        if (x2 == 0) {
          
          x$consistency <- "-"
          zz <- as.data.frame(x$consistency)
          zz$coverage <- "-"
          zz$solution <- "All truth table rows inconsistent"
          zz$model <- "-"
          zz$partition <- part
          zz$type <- type
          zz <- zz[!duplicated(zz), ]
          colnames(zz)[1] <- "consistency"
          
        } else if (x2 == 1 & solution == "P") {
          
          x$consistency <- "-"
          zz <- as.data.frame(x$consistency)
          zz$coverage <- "-"
          zz$solution <- "All truth table rows consistent"
          zz$model <- "-"
          zz$partition <- part
          zz$type <- type
          zz <- zz[!duplicated(zz), ]
          colnames(zz)[1] <- "consistency"
          
        } else {
          
          if (solution == "C") {
            x <- QCA::minimize(x1, explain = "1", include = "1", 
                               details = T, show.cases = T, 
                               all.sol = T, row.dom = F)
          } else if (solution == "P") {
            x <- QCA::minimize(x1, explain = "1", include = "?", 
                               details = T, show.cases = T, 
                               all.sol = T, row.dom = F)
          } else {
            x <- QCA::minimize(x1, explain = "1", include = "1", 
                               details = T, show.cases = T, 
                               all.sol = T, row.dom = F)
          }
          
          # x$CONS <- x$IC$incl.cov$incl[1]
          x$consistency <- x$IC$sol.incl.cov[1]
          
          if (is.null(x$consistency)) {
            BSP2 <- x$IC
            BSP3 <- BSP2$individual
            neu <- sapply(BSP3, function(x) x[2])
            x$consistency <- sapply(neu, function(x) x[1])
            x$consistency <- unlist(x$consistency)
          }
          
          zz <- as.data.frame(x$consistency)
          x$coverage <- x$IC$sol.incl.cov[3]
          
          if (is.null(x$coverage)) {
            BSP2 <- x$IC
            BSP3 <- BSP2$individual
            neu <- sapply(BSP3, function(x) x[2])
            x$coverage <- sapply(neu, function(x) x[3])
            x$coverage <- unlist(x$coverage)
          }
          
          zz$coverage <- as.numeric(x$coverage)
          x$solution <- x$solution[]
          tete <- list(cons = x$solution[])
          neux <- lapply(tete$cons, paster)
          neuxx <- unlist(neux)
          zz$solution <- neuxx
          numberrows <- nrow(zz)
          if (numberrows == 1) {
            zz$model <- 1
          } else {
            zz$model <- as.numeric(rownames(zz))
          }
          zz$partition <- part
          zz$type <- type
          colnames(zz)[1] <- "consistency"
        }
        
      } else {
        x$consistency <- "-"
        zz <- as.data.frame(x$consistency)
        zz$coverage <- "-"
        zz$solution <- "All rows remainders (n < frequency cutoff)"
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
  
  #### Rounding ####
  total$consistency[total$model == "-"] <- NA
  total$coverage[total$model == "-"] <- NA
  total$consistency <- as.numeric(total$consistency)
  total$coverage <- as.numeric(total$coverage)
  total <- total[, c(6, 5, 3, 4, 1, 2)]
  total$solution <- gsub("\\*"," * ",total$solution)
  return(total)}}}}
}