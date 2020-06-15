

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

schwarz2016 <- read.csv("Schwarz2016.csv", as.is = T, sep = ";")

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
