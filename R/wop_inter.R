#' Calculation of weight of partitions in pooled solution parameters
#' for intermediate solution
#' 
#' \code{wop_inter} calculates the weight of partitions in the pooled
#' solution parameters (consistency, coverage) for the conservative
#' and parsimonious solution. 
#' 
#' @importFrom plyr ldply
#' @importFrom testit has_error
#' @importFrom purrr map
#' @import QCA
#' 
#' @param Dataset Calibrated pooled dataset for partitioning and minimization
#' @param units Units defining the within-dimension of data (time series)
#' @param time Periods defining the between-dimension of data (cross sections)
#' @param cond Conditions used for the pooled analysis
#' @param out Outcome used for the pooled analysis
#' @param n_cut Frequency cut-off for designating truth table rows as observed
#' @param incl_cut Inclusion cut-off for designating truth table rows as
#' consistent
#' @param intermediate A vector of directional expectations to derive 
#' intermediate solutions
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
#' @return A dataframe with information about the weight of the partitions 
#' for pooled consistency and coverage scores and the following columns:
#' 
#' * \code{type}: The type of the partition. \code{pooled} are rows with information
#' on the pooled data; \code{between} is for cross-section partitions;
#' \code{within} is for time-series partitions
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
#' * \code{num}: The contribution of the data to the numerator of the consistency
#' formula. For the pooled data, the value draws on all cases. For
#' partitions, the value represents the weight of the partition. The larger the
#' value, the larger the weight. The sum over all partition-specific weights is 
#' equal to the value for the pooled data.
#' * \code{denom}: The contribution of the data to the denominator of the 
#' consistency formula. For the pooled data, the value draws on all cases.
#' For partitions, the value represents the weight of the partition. 
#' The larger the value, the larger the weight. The sum over all 
#' partition-specific weights is equal to the value for the pooled data.
#' @md
#' 
#' @examples
#' data(schwarz2016)
# Schwarz_wop_inter <- partition_min_inter(
#   Schwarz2016,
#   units = "country", time = "year",
#   cond = c("poltrans", "ecotrans", "reform", "conflict", "attention"),
#   out = "enlarge",
#   n_cut = 1, incl_cut = 0.8,
#   intermediate = c("1", "1", "1", "1", "1"))
#'
#' @export
wop_inter <- function(dataset, units, time, cond, out, n_cut, incl_cut, 
                      intermediate, BE_cons, WI_cons,
                      BE_ncut,WI_ncut) {
  
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
    #datapart <- x[out]
    #datapart <- unlist(datapart[1])
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
      
      s <- testit::has_error(susu <- 
                               try(suppressWarnings(QCA::truthTable(x,
                                                                    outcome = out, 
                                                                    conditions = cond, 
                                                                    incl.cut1 = x[, ncol(x)-1][1], 
                                                                    n.cut = x[, ncol(x)][1])), 
                                   silent = TRUE))
      
      if (s == F) {
        x1 <- try(suppressWarnings(QCA::truthTable(x, outcome = out, 
                                                   conditions = cond, 
                                                   incl.cut1 = x[, ncol(x)-1][1], 
                                                   n.cut = x[, ncol(x)][1])),
                  silent = TRUE)
        
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
          t <- testit::has_error(susu <- try(QCA::minimize(x1, 
                                                           explain = "1",
                                                           dir.exp = intermediate, 
                                                           include = "?", 
                                                           details = T, 
                                                           show.cases = T, 
                                              all.sol = T, row.dom = F), 
                                             silent = TRUE))
          
          if (t == F) {
            
            x3 <- QCA::minimize(x1, explain = "1", 
                                dir.exp = intermediate, 
                                include = "?", details = T, 
                                show.cases = T, all.sol = T, 
                           row.dom = F)
            
            a <- x3$i.sol
            ININ <- lapply(a, intersol)
            ININ1 <- lapply(ININ, intersol2)
            zz <- unlist(ININ1)
            zz <- as.data.frame(zz)
            
            
            denomfunc <- function(x) {
              pimlength <- as.numeric(ncol(x))
              x$max <- do.call(pmax, x[1:pimlength])
              deno <- sum(x$max)
              deno
            }
            
            numfunc2 <- function(x) {
              pimlength <- as.numeric(ncol(x))
              x$max <- do.call(pmax, x[1:pimlength])
              outcomedata <- as.numeric(unlist(x3$tt$initial.data[out]))
              x$out <- outcomedata
              x$min <- with(x, pmin(max, out))
              nu <- sum(x$min)
              nu
            }
            
            getpims <- function(x){
              if(length(x$solution)<2){
                x <- x$pims
              }else{
                allpi <- x$IC$individual
                x <- purrr::map(allpi,3)
              }
              x
            }
            
            pimsisi <- lapply(a, getpims)
            
            getpims2 <- function(x){
              if(is.data.frame(x)==T){
                pimlength <- as.numeric(ncol(x))
                x$max <- do.call(pmax, x[1:pimlength])
                x <- sum(x$max)
              }else{
                x <- unlist(lapply(x, denomfunc))
              }
              x
            }
            
            denom <- unlist(lapply(pimsisi, getpims2))
            
            getpims3 <- function(x, x3, out){
              if(is.data.frame(x)==T){
                pimlength <- as.numeric(ncol(x))
                x$max <- do.call(pmax, x[1:pimlength])
                outcomedata <- as.numeric(unlist(x3$tt$initial.data[out]))
                x$out <- outcomedata
                x$min <- with(x, pmin(max, out))
                x <- sum(x$min)
              }else{
                x <- unlist(lapply(x, numfunc2))
              }
              x
            }
            
            num <- unlist(lapply(pimsisi, x3=x3, out = out, getpims3))
            
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
        
        SOL <- "No combinations at this frequency cutoff"
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
