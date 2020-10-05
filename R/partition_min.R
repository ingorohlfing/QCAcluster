#' Generation of conservative or parsimonious solution for individual 
#' partitions
#' 
#' \code{partition_min} decomposes clustered data into individual
#' partitions. For panel data, for example, these can be cross sections,
#' time series or both. The function derives an individual solution for
#' each partition and the pooled data to assess the robustness of the 
#' solutions in a comparative perspective.
#'
#' @param dataset Calibrated pooled dataset that is partitioned and minimized for
#' deriving the pooled solution.
#' @param units Units defining the within-dimension of data (time series)
#' @param time Periods defining the between-dimension of data (cross sections)
#' @param cond Conditions used for minimization
#' @param out Outcome used for minimization
#' @param n_cut Frequency cut-off for designating truth table rows as observed
#' as opposed to designating them as remainders.
#' @param incl_cut Inclusion (a.k.a. consistency) cut-off for designating 
#' truth table rows as consistent.
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
#'
#' @return A dataframe summarizing the partition-specific and pooled solutions
#' with the following columns:
#' 
#' * \code{type}: The type of the partition. \code{pooled} are rows with information
#' on the pooled data; \code{between} is for cross-section partitions;
#' \code{within} is for time-series partitions
#' * \code{partition}: Specific dimension of the partition at hand. For 
#' between-dimension, the unit identifiers are included here  (argument \code{units}).
#' For the within-dimension, the time identifier are listed (argument \code{time}).
#' The entry is \code{-} for the pooled data without partitions.
#' * consistency
#' @md
#'
#' @examples
#' # loading data from Thiem (EPSR, 2011; see data documentation)
#' data(Thiem2011)
#' 
#' # running function for parsimonious solution
# Thiem_pars_1 <- partition_min(
#   dataset = Thiem2011,
#   units = "units", time = "time",
#   cond = c("fedismfs", "homogtyfs", "powdifffs", "comptvnsfs", "pubsupfs", "ecodpcefs"),
#   out = "memberfs",
#   n_cut = 6, incl_cut = 0.8,
#   solution = "P",
#   BE_cons = c(0.9, 0.8, 0.7, 0.8, 0.6, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8),
#   WI_cons = c(0.5, 0.8, 0.7, 0.8, 0.6, rep(0.8, 10)))
#' 
#' @export
partition_min <- function(dataset, 
                          units, time, 
                          cond, out, 
                          n_cut, incl_cut, 
                          solution, 
                          BE_cons, WI_cons) {
  
  # turning of warnings
  # options(warn = -1)
  
  # splitting the data if time and unit values are available
  x <- dataset
  if (missing(units)) {
    colnames(x)[which(names(x) == time)] <- "time"
    x <- x[with(x, order(time)), ]
    xB <- x
    xB$consis <- rep(incl_cut, times = nrow(x))
    xxx <- 1
    BE_list <- split(xB, xB[, "time"])
  } else if (missing(time)) {
    colnames(x)[which(names(x) == units)] <- "units"
    x <- x[with(x, order(units)), ]
    xW <- x
    xW$consis <- rep(incl_cut, times = nrow(x))
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
      
      # s <- has_error(truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[,ncol(x)][1], n.cut = n_cut))
      s <- testit::has_error(susu <- try(QCA::truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)][1], n.cut = n_cut), silent = TRUE))
      
      if (s == F) {
        x1 <- try(QCA::truthTable(x, outcome = out, conditions = cond, incl.cut1 = x[, ncol(x)][1], n.cut = n_cut), silent = TRUE)
        
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
            x <- QCA::minimize(x1, explain = "1", include = "1", details = T, show.cases = T, all.sol = T, row.dom = F)
          } else if (solution == "P") {
            x <- QCA::minimize(x1, explain = "1", include = "?", details = T, show.cases = T, all.sol = T, row.dom = F)
          } else {
            x <- QCA::minimize(x1, explain = "1", include = "1", details = T, show.cases = T, all.sol = T, row.dom = F)
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
    WI_list1 <- lapply(WI_list, pqmcc)
    PO_list1 <- lapply(PO_list, pqmcc)
    dff2 <- plyr::ldply(WI_list1)[, -1]
    dff3 <- plyr::ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    total <- rbind(dff3, dff2)
    
  } else if (missing(units)) {
    BE_list1 <- lapply(BE_list, pqmcc)
    PO_list1 <- lapply(PO_list, pqmcc)
    dff1 <- plyr::ldply(BE_list1)[, -1]
    dff3 <- plyr::ldply(PO_list1)[, ]
    dff3$type <- "pooled"
    dff3$partition <- "-"
    total <- rbind(dff3, dff1)
    
  } else {
    BE_list1 <- lapply(BE_list, pqmcc)
    WI_list1 <- lapply(WI_list, pqmcc)
    PO_list1 <- lapply(PO_list, pqmcc)
    
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
  return(total)
}
