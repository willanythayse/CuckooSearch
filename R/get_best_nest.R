
source("R/fobj.R")

## Find the current best nest
get_best_nest <- function (nest,newnest,fitness){
  # Evaluating all new solutions
  for (j in 1:dim(nest)[1]) {
    fnew <- fobj(newnest[j,])
    if (fnew <= fitness[j]) {
      fitness[j] <- fnew
      nest[j,] <- newnest[j,]
    }
  }
  # Find the current best
  fmin <- min(fitness)
  best <- nest[which.min(fitness)]
  return(list('fmin' = fmin, 'best' = best, 'nest' = nest, 'fitness' = fitness))
}