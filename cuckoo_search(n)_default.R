#####################################
##      Willany Thayse 04/23/19    ##
#####################################
##     https://bit.ly/2uAr4LU      ##
#####################################

# -----------------------------------------------------------------
# Cuckoo Search (CS) algorithm by Xin-She Yang and Suash Deb      #
# Programmed by Xin-She Yang at Cambridge University              #
# Programming dates: Nov 2008 to June 2009                        #
# Last revised: Dec  2009   (simplified version for demo only)    #
# -----------------------------------------------------------------
# Papers -- Citation Details:
# 1) X.-S. Yang, S. Deb, Cuckoo search via Levy flights,
# in: Proc. of World Congress on Nature & Biologically Inspired
# Computing (NaBIC 2009), December 2009, India,
# IEEE Publications, USA,  pp. 210-214 (2009).
# http://arxiv.org/PS_cache/arxiv/pdf/1003/1003.1594v1.pdf 
# 2) X.-S. Yang, S. Deb, Engineering optimization by cuckoo search,
# Int. J. Mathematical Modelling and Numerical Optimisation, 
# Vol. 1, No. 4, 330-343 (2010). 
# http://arxiv.org/PS_cache/arxiv/pdf/1005/1005.2908v2.pdf
# ----------------------------------------------------------------#
# This demo program only implements a standard version of         #
# Cuckoo Search (CS), as the Levy flights and generation of       #
# new solutions may use slightly different methods.               #
# The pseudo code was given sequentially (select a cuckoo etc),   #
# but the implementation here uses Matlab's vector capability,    #
# which results in neater/better codes and shorter running time.  # 
# This implementation is different and more efficient than the    #
# the demo code provided in the book by 
#    "Yang X. S., Nature-Inspired Metaheuristic Algoirthms,       # 
#     2nd Edition, Luniver Press, (2010).                 "       #
# --------------------------------------------------------------- #
# =============================================================== #
# Notes:                                                          #
# Different implementations may lead to slightly different        #
# behavour and/or results, but there is nothing wrong with it,    #
# as this is the nature of random walks and all metaheuristics.   #
# -----------------------------------------------------------------
cuckoo_search <- function (n = 25, maxIter = 10^5, pa = 0.25, Tol = 1.0e-5, nd = 2, lb, ub) {
  # Max Iteration
#  maxIter = 10^5
#  # Number of nests (or different solutions)
#  n <- 25
#  # Discovery rate of alien eggs/solutions
#  pa <- 0.25
#  ## Change this if you want to get better results
#  # Tolerance
#  Tol <- 1.0e-5
#  ## Simple bounds of the search domain
#  # Lower bounds
#  nd <- 2   #Dimension
  Lb <- matrix(lb, 1, nd)
  # Upper bounds
  Ub <- matrix(ub, 1, nd)
  # Random initial solutions
  nest <- matrix(0, n, nd)
  for (i in 1:n) {
    nest[i,] <- Lb+(Ub-Lb)* runif(nd)
  }
  # Get the current best
  fitness <- 10^10 * matrix(1, n, 1) 
  current <- get_best_nest(nest, nest, fitness)
  fmin <- current$fmin
  bestnest <- current$best
  nest <- current$nest
  fitness <- current$fitness
  N_iter <- 0
  N_iter_rand <- n
  
  while (N_iter < maxIter) {
   # Generate new solutions (but keep the current best)
   new_nest <- get_cuckoos(nest,bestnest,Lb,Ub)
   new_best <- get_best_nest(nest,new_nest,fitness)
   fnew <- new_best$fmin 
   best <- new_best$best
   nest <- new_best$nest
   fitness <- new_best$fitness
   
   # Update the counter
   N_iter <- N_iter + n

   # Discovery and randomization
   new_nest <- empty_nests(nest,Lb,Ub,pa) 
   
   # Evaluate this set of solutions
   new_best <- get_best_nest(nest,new_nest,fitness)
   fnew <- new_best$fmin 
   best <- new_best$best
   nest <- new_best$nest
   fitness <- new_best$fitness
   
   # Update the counter again
    N_iter_rand <- N_iter * 2

   # Find the best objective so far  
   if (fnew < fmin) {
     fmin <- fnew
     bestnest <- best
   }
   print(cat('Iteration number:', N_iter, 'fitness:', fmin, 'Total number of iteration:', N_iter_rand))
 }  ## End of iterations
  ## Post-optimization processing
  ## Display all the nests
  print(cat('Number loop iterations=', N_iter))
  print(cat('Total number of iteration= ', N_iter_rand))
  return(list('fmin' = fmin, 'bestnest' = bestnest))
}

# --------------- All subfunctions are list below ------------------
# Get cuckoos by ramdom walk
get_cuckoos <- function (nest,best,Lb,Ub) {
  ## Levy flights
  n <- dim(nest)[1]
  # Levy exponent and coefficient
  # For details, see equation (2.21), Page 16 (chapter 2) of the book
  # X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
  beta <- 3/2
  sigma <- (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta)
  for (j in 1:n) {
    s <- nest[j,]
    size <- dim(nest)[2]
    # This is a simple way of implementing Levy flights
    # For standard random walks, use step=1;
    ## Levy flights by Mantegna's algorithm
    u <- runif(size)*sigma
    v <- runif(size)
    step <- (u/abs(v))^(1/beta)
    # I#n the next equation, the difference factor (s-best) means that 
    # when the solution is the best solution, it remains unchanged.     
    stepsize <- 0.01*step*(s-best)
    # Here the factor 0.01 comes from the fact that L/100 should the typical
    # step size of walks/flights where L is the typical lenghtscale; 
    # otherwise, Levy flights may become too aggresive/efficient, 
    # which makes new solutions (even) jump out side of the design domain 
    # (and thus wasting evaluations).
    # Now the actual random walks or flights
    s <- s + (stepsize*rnorm(size))
    # Apply simple bounds/limits
    nest[j,]=simplebounds(s,Lb,Ub)
  }
  return(nest)
}

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

#### Application of simple constraints ####
simplebounds <- function (s,Lb,Ub) {
  ns_tmp <- s
  ## Apply the lower bound 
  i <- ns_tmp < Lb
  ns_tmp[i] <- Lb[i]
  ## Apply the upper bounds 
  j <- ns_tmp > Ub 
  ns_tmp[j] <- Ub[j]
  # Update this new move
  s <- ns_tmp
  return(s)
}

## You can replace the following by your own functions
# A d-dimensional objective function
fobj <- function (u) {
  ## d-dimensional sphere function sum_j=1^d (u_j-1)^2. 
  #  with a minimum at (1,1, ...., 1); 
  return(cs_aux(u))
}


#cuckoo_search(n = 25, maxIter = 10^5, pa = 0.25, Tol = 1.0e-5, nd = 2, lb, ub)
cuckoo_search(25, 10^5, 0.25, 1.0e-5, 2, lb, ub)