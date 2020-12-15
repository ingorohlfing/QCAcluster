# Functions
# 1. ambig: returns information about ambiguity

# the function prints out FALSE if the solution is not ambiguous
# if it is ambiguous, then prints out the length of the ambiguity
ambig <- function(sol){
  
  a <- length(sol$IC$individual)
  
  if (a < 2) {
    print(FALSE)
  } else {
    print(a)
  }
}


# 2. select.fnc: Selecting the details for each of the ambiguous models
# sol: the solution , elmnt: element that needs to be extracted (as it is saved by the QCA package)
select.fnc <- function(sol, elmnt){
  tmp <- map(sol$IC$individual, function(x) x[[elmnt]])
  return(tmp)
}



# checking the functions (with the reproduction article by Csergo2017)
# Ambiguity
ambig(parspos1)
ambig(parspos2) 
ambig(poscons2)

# Extracting details 
select.fnc(parspos2, "incl.cov")
select.fnc(parspos1, "incl.cov")

