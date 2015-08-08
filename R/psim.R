#------------------------------
# PVALUE SIMULATION
#------------------------------

simulate <- function(n) {
  p.sample <- runif(n=n, min=0, max=1)
  return(tukey(p.sample))
}

generateDistributions <- function(max, k) {
  distributions <- list(1:max)
  for (i in 1:max) {
    distributions[[i]] <- replicate(k, simulate(n=i))
  }
  return(distributions)
}

psim <- function(p, list) {
  list <- unlist(list) 
  return(sum(list <= p) / length(list))
}
