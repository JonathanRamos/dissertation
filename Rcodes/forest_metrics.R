require("MASS")
require("cluster")

forest.metrics <- function(dat, epsilon){

    W = gaussian.similarity(dat,epsilon)
    
    ## L = -W
    ## n = nrow(dat)
    ## for(i in 1:n){
    ##     L[i,i] = -sum(L[i,])
    ## }

    L <- laplacian(W)

    Q = solve(eye(n) + L)
    H = kappa(Q)

    return(H)
} 

# This variant of forest metrics is based on the directed Laplacian matrix 

forest.metrics.directed1 <- function(dat, epsilon){
    W <- gaussian.similarity.directed(dat, epsilon)
    L <- laplacian(W)
    n <- nrow(W)

    Q = solve(diag(n) + 10*L)
    H = kappa(Q)
    return(H)
}
    
