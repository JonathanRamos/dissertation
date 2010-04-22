## This variant of forest metrics is based on the Laplacian matrix.
## It should work for directed graphs also.

forest.metrics <- function(L, alpha = 1){

  n <- nrow(L)
  Q = solve(eye(n) + alpha*L)
  H = kappa(Q)

  return(H)
} 

