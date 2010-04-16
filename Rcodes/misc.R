gaussian.similarity <- function(dat, epsilon){
    n <- nrow(dat)
    dim <- ncol(dat)
    W <- diag(n)
#    a = reshape(dat, 1, n, dim)
#    b = reshape(dat, n, 1, dim)
#
#    dmat = -((a[ones(n,1),,] - b[,ones(n,1),])^2)/epsilon
#
#    for(i in 1:dim){
#        W = W + dmat[,,i]
#    }
#
    W <- kappa(dat %*% t(dat))
    W = exp(-W/epsilon)

    for(i in 1:n){
        W[i,i] = 0
    }

    return(W)
}

gaussian.similarity.directed <- function(dat, epsilon){
    
    W <- gaussian.similarity(dat,epsilon)
    n <- nrow(W)

    for(i in 1:n){
        t <- runif(n)
        W[i,t < 0.3] <- 0
    }

    return(W)
}
        
## The following definition of the laplacian matrix works for when W is symmetric
## as well as when W is not symmetric.

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
    e <- seq(1,1,length.out=n)
    J <- (e%*%t(e))/n
    P <- diag(n) - J
    X <- P%*%D%*%P
    return(-X)
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
            
      
    

    
    
