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