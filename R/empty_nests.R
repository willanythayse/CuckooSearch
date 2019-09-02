
source("R/simplesbounds.R")

## Replace some nests by constructing new solutions/nests
empty_nests <- function (nest,Lb,Ub,pa){
  # A fraction of worse nests are discovered with a probability pa
  #  n <- dim(nest)[1,2]
  n <- dim(nest)[1]
  # Discovered or not -- a status vector
  K <- matrix(runif(n*dim(nest)[2]), n, dim(nest)[2] ) > pa
  
  # In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
  # this cuckoo's egg is less likely to be discovered, thus the fitness should 
  # be related to the difference in solutions.  Therefore, it is a good idea 
  # to do a random walk in a biased way with some random step sizes.  
  ## New solution by biased/selective random walks
  stepsize <- runif(1)*(nest[sample(n),]-nest[sample(n),])
  new_nest <- nest+stepsize*K
  #for (j in 1:dim(new_nest)[1]) {
  #  s <- new_nest
  #  new_nest = simplebounds(s,Lb,Ub)  
  #}
  #
  return(new_nest)
}