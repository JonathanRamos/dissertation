
## Diffusion distances takes as input a probability transition matrix P

diffusion.distance <- function(P, t){
  
  n <- nrow(P)
  decomposition <- eigen(t(P),symmetric=FALSE)
  w <- abs(decomposition$vectors[,1])
  w <- w/sum(w)

  P.t <- mtx.exp(P,t)
  D.squared <- kappa(P.t %*% diag(1/w) %*% t(P.t))
  D <- sqrt(D.squared)

  return(D)
}

diffusion.map <- function(P, t, k=2){

    n <- nrow(P)
    decomposition <- eigen(t(P),symmetric=FALSE)
    w <- abs(decomposition$vectors[,1])
    w <- w/sum(w)
    
    Pi.sqrt <- diag(sqrt(w))
    Pi.sqrt.inv <- diag(1/sqrt(w))
    A <- Pi.sqrt %*% P %*% Pi.sqrt.inv
    
    decomp <- eigen(A)
    eigen.vals <- decomp$values[2:(k+1)]
    eigen.vects <- decomp$vectors[,2:(k+1)]
    eigen.vals.powered <- eigen.vals^t

    Psi <- Pi.sqrt.inv %*% eigen.vects %*% diag(eigen.vals.powered)
    return(Psi)
}
