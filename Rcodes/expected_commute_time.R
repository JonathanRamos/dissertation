ect <- function(P){
  
    n <- nrow(P)
    decomposition <- eigen(t(P),symmetric=FALSE)
    w <- abs(decomposition$vectors[,1])
    w <- w/sum(w)
    Q <- seq(1,1,length.out=n) %*% t(w)

    Z <- solve(diag(n) - P + Q)
    H <- kappa(Z %*% diag(1/w))
    D <- sqrt(H)
    
    return(D)
}

ect.map <- function(P,k=2){
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
    eigen.vals.transformed <- 1/sqrt(1 - eigen.vals)

    Psi <- Pi.sqrt.inv %*% eigen.vects %*% diag(eigen.vals.transformed)

    return(Psi)
}

ect.sqrd <- function(P){
    n <- nrow(P)
    decomposition <- eigen(t(P),symmetric=FALSE)
    w <- abs(decomposition$vectors[,1])
    w <- w/sum(w)
    Q <- seq(1,1,length.out=n) %*% t(w)

    Z <- solve(diag(n) - P + Q)
    Z.sqrd <- Z%*%Z
    H <- kappa(Z.sqrd %*% diag(1/w))
    D <- sqrt(H)
    
    return(D)
}

ect.cube <- function(P){
    n <- nrow(P)
    decomposition <- eigen(t(P),symmetric=FALSE)
    w <- abs(decomposition$vectors[,1])
    w <- w/sum(w)
    Q <- seq(1,1,length.out=n) %*% t(w)

    Z <- solve(diag(n) - P + Q)
    Z.cube <- Z%*%Z%*%Z
    H <- kappa(Z.cube %*% diag(1/w))
    D <- sqrt(H)
    
    return(D)
}

ect.pow <- function(P,t){
    n <- nrow(P)
    decomposition <- eigen(t(P),symmetric=FALSE)
    w <- abs(decomposition$vectors[,1])
    w <- w/sum(w)
    Q <- seq(1,1,length.out=n) %*% t(w)

    Z <- solve(diag(n) - P + Q)
    Z.powered = mtx.exp(Z,t)
    H <- kappa(Z.powered %*% diag(1/w))
    D <- sqrt(H)
    
    return(D)
}


## k = 40 
## Psi <- M$map[,1:k]
## H <- sqrt(M$dist)
## D1 <- as.matrix(dist(Psi, method="euclidean"))
## norm1 <- norm(D1 - H, type ="F")
## Z <- cmdscale(H,eig = TRUE, k)
## D2 <- as.matrix(dist(Z$points,method="euclidean"))
## norm2 <- norm(D2 - H, type = "F")
