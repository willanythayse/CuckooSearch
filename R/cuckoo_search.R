# http://r-pkgs.had.co.nz/ Some useful keyboard shortcuts for package authoring: Install Package: 'Ctrl + Shift + B' Check
# Package: 'Ctrl + Shift + E' Test Package: 'Ctrl + Shift + T'

##################################### Willany Thayse 04/23/19 ## https://bit.ly/2uAr4LU ##

# ----------------------------------------------------------------- Cuckoo Search (CS) algorithm by Xin-She Yang and Suash
# Deb # Programmed by Xin-She Yang at Cambridge University # Programming dates: Nov 2008 to June 2009 # Last revised: Dec
# 2009 (simplified version for demo only) # ----------------------------------------------------------------- Papers --
# Citation Details: 1) X.-S. Yang, S. Deb, Cuckoo search via Levy flights, in: Proc. of World Congress on Nature &
# Biologically Inspired Computing (NaBIC 2009), December 2009, India, IEEE Publications, USA, pp. 210-214 (2009).
# http://arxiv.org/PS_cache/arxiv/pdf/1003/1003.1594v1.pdf 2) X.-S. Yang, S. Deb, Engineering optimization by cuckoo
# search, Int. J. Mathematical Modelling and Numerical Optimisation, Vol. 1, No. 4, 330-343 (2010).
# http://arxiv.org/PS_cache/arxiv/pdf/1005/1005.2908v2.pdf
# ----------------------------------------------------------------# This demo program only implements a standard version
# of # Cuckoo Search (CS), as the Levy flights and generation of # new solutions may use slightly different methods.  #
# The pseudo code was given sequentially (select a cuckoo etc), # but the implementation here uses Matlab's vector
# capability, # which results in neater/better codes and shorter running time.  # This implementation is different and
# more efficient than the # the demo code provided in the book by 'Yang X. S., Nature-Inspired Metaheuristic Algoirthms, #
# 2nd Edition, Luniver Press, (2010).  ' # --------------------------------------------------------------- #
# =============================================================== # Notes: # Different implementations may lead to
# slightly different # behavour and/or results, but there is nothing wrong with it, # as this is the nature of random
# walks and all metaheuristics.  # -----------------------------------------------------------------

# This Function should be find the best nest 'global minimum'
# The function of optimization (aux function) should be defined with your Lower bound 'lb' and Upper bound 'ub'

source("get_best_nest.R")
source("get_cuckoos.R")
source("empty_nests.R")


cuckoo_search <- function(n = 25, maxIter = 10^5, pa = 0.25, Tol = 1e-05, nd = 2, lb = 0, ub = pi) {
    # Max Iteration maxIter = 10^5 # Number of nests (or different solutions) n <- 25 # Discovery rate of alien eggs/solutions
    # pa <- 0.25 ## Change this if you want to get better results # Tolerance Tol <- 1.0e-5 ## Simple bounds of the search
    # domain # Lower bounds nd <- 2 #Dimension
    Lb <- matrix(lb, 1, nd)
    # Upper bounds
    Ub <- matrix(ub, 1, nd)
    # Random initial solutions
    nest <- matrix(0, n, nd)
    for (i in 1:n) {
        nest[i, ] <- Lb + (Ub - Lb) * runif(nd)
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
        new_nest <- get_cuckoos(nest, bestnest, Lb, Ub)
        new_best <- get_best_nest(nest, new_nest, fitness)
        fnew <- new_best$fmin
        best <- new_best$best
        nest <- new_best$nest
        fitness <- new_best$fitness

        # Update the counter
        N_iter <- N_iter + n

        # Discovery and randomization
        new_nest <- empty_nests(nest, Lb, Ub, pa)

        # Evaluate this set of solutions
        new_best <- get_best_nest(nest, new_nest, fitness)
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
        print(cat("Iteration number:", N_iter, "fitness:", fmin, "Total number of iteration:", N_iter_rand))
    }  ## End of iterations
    ## Post-optimization processing Display all the nests
    print(cat("Number loop iterations=", N_iter))
    print(cat("Total number of iteration= ", N_iter_rand))
    return(list(fmin = fmin, bestnest = bestnest))
}
