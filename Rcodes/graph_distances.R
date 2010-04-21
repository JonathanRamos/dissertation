require("Matrix")

gaussian.similarity <- function(dat, epsilon){
    n <- nrow(dat)
    W <- diag(n)
    
    W <- kappa(dat %*% t(dat))
    W = exp(-W/epsilon)

    diag(W) <- 0

    return(W)
}

## Simulate Gaussian similarity for directed graphs. The simulation is done
## by setting entries of the Gaussian similarity matrix to 0 at random. 

gaussian.similarity.directed <- function(dat, epsilon, sparsity = 0.5){
    
    W <- gaussian.similarity(dat,epsilon)
    n <- nrow(W)

    for(i in 1:n){
        t <- runif(n)
        W[i,t < sparsity] <- 0
    }

    return(W)
}
        
laplacian <- function(W){
    S <- apply(W,1,sum)
    L <- diag(S) - W
    return(L)
}

## Returns the kappa transform of a square matrix X
kappa <- function(X){
    n <- nrow(X)
    e <- seq(1,1,length.out = n)
    x <- diag(X)
    D <- x%*%t(e) - X - t(X) + e%*%t(x)
    return(D)
}

## Returns the tau transform of a square matrix D
tau <- function(D){
    n <-  nrow(D)
    P <- projection.matrix(n)
    X <- (-1/2)*P%*%D%*%P
    
    return(X)
}

projection.matrix <- function(n){
    e <- seq(1,1,length.out=n)
    J <- (e%*%t(e))/n
    P <- diag(n) - J
}

# Return the transition matrix given a similarity measure
transition.matrix <- function(W){
    S.inv <- 1/apply(W,1,sum)
    P <- diag(S.inv)%*%W
    return(P)
}

transition.matrix.star <- function(P){
    n <- nrow(P)
    decomposition <- eigen(t(P),symmetric=FALSE)
    w <- abs(decomposition$vectors[,1])
    w <- w/sum(w)
    P.star <- diag(1/w) %*% t(P) %*% diag(w)

    return(P.star)
}
    

## Construct a new similarity matrix corresponding to the addition of a vertex
## to a graph. The new vertex will have similarity 1 with respect to the other vertices.

augmented.similarity <- function(W){
    n <- nrow(W)
    s = seq(1,1,length.out = n)
    W.new <- rbind(c(0,s),cbind(1,W))
    return(W.new)
}

# Compute the exponent of matrix
mtx.exp <- function(X,n){
    if(n != round(n)){
        n <- round(n)
        warning("rounding exponent `n` to", n)
    }
    phi <- diag(nrow = nrow(X))
    pot <- X
    while(n > 0){
        if(n %% 2)
            phi <- phi %*% pot
        n <- n %/% 2
        pot <- pot %*% pot
    }
    return(phi)
}
            
## Expected commute time

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

## Diffusion distances

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

## Forest metrics
forest.metrics <- function(L, alpha){

    n = nrow(L)
    Q = solve(eye(n) + alpha*L)
    H = kappa(Q)

    return(H)
} 

## Exponential distance
exp.dist <- function(P){
    
    decomposition <- eigen(t(P),symmetric=FALSE)
    w <- abs(decomposition$vectors[,1])
    w <- w/sum(w)
    
    H = expm(P)%*%diag(1/w)
    D = kappa(H)
}

## data <- clusters2.data(50)
## W <- gaussian.similarity.directed(data, 1, sparsity = 0.1)
## P <- transition.matrix(W)
## D <- exp.dist(P)
## B <-  tau(D)
## L <- eigen(B,only.values = TRUE)
## L$values[L$values < 0]
## k = 40 
## Psi <- M$map[,1:k]
## H <- sqrt(M$dist)
## D1 <- as.matrix(dist(Psi, method="euclidean"))
## norm1 <- norm(D1 - H, type ="F")
## Z <- cmdscale(H,eig = TRUE, k)
## D2 <- as.matrix(dist(Z$points,method="euclidean"))
## norm2 <- norm(D2 - H, type = "F")

