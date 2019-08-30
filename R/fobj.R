source("R/michal.R")

## You can replace the following by your own functions
# A d-dimensional objective function
fobj <- function (u) {
  ## d-dimensional sphere function sum_j=1^d (u_j-1)^2. 
  #  with a minimum at (1,1, ...., 1); 
  return(michal(u))
}
