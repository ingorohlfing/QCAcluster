

library(QCA)
library(data.table)
library(plyr)
library(testit)
library(devtools)
library(testthat)
has_devel()

### Schwarz 2016


Schwarz_WOP_Inter <- WOPfunctionInter(schwarz2016, units = "country", time = "year", cond = c("poltrans", "ecotrans", "reform", "conflict", 
    "attention"), out = "enlarge", 1, 0.8, intermediate = c("1", "1", "1", "1", "1"))


### Grauvogel 2014


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
