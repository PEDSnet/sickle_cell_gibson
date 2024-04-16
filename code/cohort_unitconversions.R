# to convert to cells/uL
neutrophils_unitconversion <- function(data) {
    
  if (unit == "10^9/L") {
    return(neutrophils)
  } else if (unit == "10^3/uL") {
    return(neutrophils * 1000)
  } else {
    stop("Invalid unit")
  }
}