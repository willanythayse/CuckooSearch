source("R/simplesbounds.R")

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