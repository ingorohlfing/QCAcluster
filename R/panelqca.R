

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

