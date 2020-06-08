

library(QCA)
library(data.table)
library(plyr)
library(testit)
library(devtools)
library(roxygen2)
library(testthat)
library(knitr)
has_devel()


### Functions for Panel results, diversity and partition weights
#' @export

partition_params_inter <- function(x, units, time, cond, out, n_cut, incl_cut, intermediate, BE_cons, WI_cons) {
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

wop <- function(x, units, time, cond, out, n_cut, incl_cut, solution, BE_cons, WI_cons) {
    
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
                    
                    x3 <- minimize(x1, explain = "1", include = "1", details = T, show.cases = T, all.sol = T, row.dom = F)
                    
                  } else if (solution == "P") {
                    
                    x3 <- minimize(x1, explain = "1", include = "?", details = T, show.cases = T, all.sol = T, row.dom = F)
                    
                  } else {
                    
                    x3 <- minimize(x1, explain = "1", include = "1", details = T, show.cases = T, all.sol = T, row.dom = F)
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
        WI_list1 <- lapply(WI_list, pqmcc)
        PO_list1 <- lapply(PO_list, pqmcc)
        dff2 <- ldply(WI_list1)[, -1]
        dff3 <- ldply(PO_list1)[, ]
        dff3$type <- "pooled"
        dff3$partition <- "-"
        ntot <- as.numeric(mean(dff3$num))
        dtot <- as.numeric(mean(dff3$denom))
        
        total <- rbind(dff3, dff2)
    } else if (missing(units)) {
        BE_list1 <- lapply(BE_list, pqmcc)
        PO_list1 <- lapply(PO_list, pqmcc)
        
        dff1 <- ldply(BE_list1)[, -1]
        dff3 <- ldply(PO_list1)[, ]
        dff3$type <- "pooled"
        dff3$partition <- "-"
        ntot <- as.numeric(mean(dff3$num))
        dtot <- as.numeric(mean(dff3$denom))
        
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
        ntot <- as.numeric(mean(dff3$num))
        dtot <- as.numeric(mean(dff3$denom))
        
        total <- rbind(dff3, dff1, dff2)
    }
    
    return(total)
    
}


Thiem2011 <- read.table(file = "ThiemPanel051.csv", fileEncoding = "WINDOWS-1252",
                        stringsAsFactors = FALSE, header = TRUE, sep = ";")


Thiem_panelpars <- PWBfunctionCP(Thiem2011, units = "units", time = "time", cond = c("fedismfs", "homogtyfs", "powdifffs", "comptvnsfs", 
    "pubsupfs", "ecodpcefs"), out = "memberfs", 6, 0.8, solution = "P", 
    BE_cons = c(0.9, 0.8, 0.7, 0.8, 0.6, 0.8, 0.8, 0.8, 0.8, 0.8, 
    0.8), 
    WI_cons = c(0.5, 0.8, 0.7, 0.8, 0.6, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8))

Thiem_panelpars <- PWBfunctionCP(Thiem2011, units = "units", time = "time", cond = c("fedismfs", "homogtyfs", "powdifffs", "comptvnsfs", 
                                                                                     "pubsupfs", "ecodpcefs"), out = "memberfs", 2, 0.8, solution = "P")
### Schwarz 2016

setwd("D:/QCA")

schwarz2016 <- read.csv("schwarz2016.csv", as.is = T, sep = ";")

# calibration
schwarz2016$enlarge <- calibrate(schwarz2016$enlarge, type = "fuzzy", logistic = T, thresholds = c(0, 9, 21))

schwarz2016$poltrans <- calibrate(schwarz2016$poltrans, type = "fuzzy", logistic = T, thresholds = c(4, 6, 10))

schwarz2016$ecotrans <- calibrate(schwarz2016$ecotrans, type = "fuzzy", logistic = T, thresholds = c(3, 7, 10))

schwarz2016$reform <- calibrate(schwarz2016$reform, type = "fuzzy", logistic = T, thresholds = c(3, 5.6, 10))

schwarz2016$conflict <- calibrate(schwarz2016$conflict, type = "fuzzy", logistic = T, thresholds = c(5, 2, 0))

schwarz2016$attention <- calibrate(schwarz2016$attention, type = "fuzzy", logistic = T, thresholds = c(-0.96, 0, 2.93))

schwarz2016$year <- schwarz2016$Case.ID
schwarz2016$year <- gsub("[A-Z]", "", schwarz2016$year)

schwarz2016$country <- schwarz2016$Case.ID
schwarz2016$country <- gsub("[0-9]", "", schwarz2016$country)

Schwarz_panelinter_1 <- PWBfunctionINTER(schwarz2016, units = "country", time = "year", cond = c("poltrans", "ecotrans", "reform", "conflict", 
    "attention"), out = "enlarge", 1, 0.8, intermediate = c("1", "1", "1", "1", "1"))
Schwarz_panelinter_2 <- PWBfunctionINTER(schwarz2016, units = "country", cond = c("poltrans", "ecotrans", "reform", "conflict", "attention"), 
    out = "enlarge", n_cut = 1, incl_cut = 0.8, intermediate = c("1", "1", "1", "1", "1"))

Schwarz_panelcons <- PWBfunctionCP(schwarz2016, units = "country", time = "year", cond = c("poltrans", "ecotrans", "reform", "conflict", 
    "attention"), out = "enlarge", 1, 0.8, solution = "P")

Schwarz_DIV <- DIVfunction(schwarz2016, units = "country", time = "year", cond = c("poltrans", "ecotrans", "reform", "conflict", "attention"), 
    out = "enlarge", 1, 0.8)

Schwarz_WOP_Inter <- WOPfunctionInter(schwarz2016, units = "country", time = "year", cond = c("poltrans", "ecotrans", "reform", "conflict", 
    "attention"), out = "enlarge", 1, 0.8, intermediate = c("1", "1", "1", "1", "1"))


### Grauvogel 2014

grauvogel2014 <- read.csv("grauvogel2014.csv", as.is = T, sep = ";")
grauvogel2014$Sender <- gsub("EU ", "EU", grauvogel2014$Sender, fixed = FALSE)
grauvogel2014$Sender <- gsub("UN ", "UN", grauvogel2014$Sender, fixed = FALSE)


grauvogel_panelinter <- PWBfunctionINTER(grauvogel2014, units = "Sender", cond = c("Comprehensiveness", "Linkage", "Vulnerability", "Repression", 
    "Claims"), out = "Persistence", n_cut = 1, incl_cut = 0.75, intermediate = c("1", "1", "1", "1", "1"))

grauvogel_panelp <- PWBfunctionCP(grauvogel2014, units = "Sender", cond = c("Comprehensiveness", "Linkage", "Vulnerability", "Repression", 
    "Claims"), out = "Persistence", n_cut = 1, incl_cut = 0.75, solution = "C")

grauvogel_DIV <- DIVfunction(grauvogel2014, units = "Sender", cond = c("Comprehensiveness", "Linkage", "Vulnerability", "Repression", 
    "Claims"), out = "Persistence", n_cut = 1, incl_cut = 0.75)

grauvogel_WOP <- WOPfunctionInter(grauvogel2014, units = "Sender", cond = c("Comprehensiveness", "Linkage", "Vulnerability", "Repression", 
    "Claims"), out = "Persistence", n_cut = 1, incl_cut = 0.75, intermediate = c("1", "1", "1", "1", "1"))


### Schneider 2014

schneider2014 <- read.csv("schneider2014.csv", as.is = T, sep = ";")

schneider2014$year <- schneider2014$X
schneider2014$year <- gsub("[A-Z]", "", schneider2014$year)
unique(schneider2014$year)

schneider2014_panelp <- PWBfunctionCP(schneider2014, time = "year", cond = c("high_lmx_l2", "high_wc_l2", "high_ud_l2", "high_epl_l2"), 
    out = "low_fde_l2", n_cut = 1, incl_cut = 0.8, solution = "P")

schneider2014_DIV <- DIVfunction(schneider2014, time = "year", cond = c("high_lmx_l2", "high_wc_l2", "high_ud_l2", "high_epl_l2"), out = "low_fde_l2", 
    n_cut = 1, incl_cut = 0.8)

schneider2014_WOP <- WOPfunctionCP(schneider2014, time = "year", cond = c("high_lmx_l2", "high_wc_l2", "high_ud_l2", "high_epl_l2"), 
    out = "low_fde_l2", n_cut = 9, incl_cut = 0.8, solution = "P")
