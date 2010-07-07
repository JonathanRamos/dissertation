require("MASS")
require("Matrix")
## B.list is a list of fallible inner product matrices
## dim.max is the maximum dimension of the group space G
## dim.list is the list of dimension for each projected space

max0 <- function(x){
    if(x >= 0)
      return (x)
    else
      return (0)
}

initialize.mds <- function(B.list, dim.max){

    N <- nrow(B.list[[1]])
    B.group <- diag(N)
    diag(B.group) <- 0

    for( i in 1:length(B.list)){
        B.group <- B.group + B.list[[i]]
    }

    B.group <- B.group/length(B.list)

    G.list <- list()
    for( i in 1:length(B.list)){
        B.i <- B.list[[i]]
        B.i.eigen <- eigen(B.i, only.values = FALSE)
        truncated.eigen = sapply(B.i.eigen$values, max0)
        U <- B.i.eigen$vectors[, c(1:dim.max)]
        G.list[[i]] <- U%*%diag(sqrt(truncated.eigen[c(1:dim.max)]))
    }
        

    ## B.group.eigen <- eigen(B.group, only.values = FALSE)
            
    ## U <- B.group.eigen$vectors
    ## U.d <- U[,c(1:dim.max)]

    ## G <- U.d
    ## G <- cmdscale(B.group, k = dim.max)
    objective.value <- 0
    min.objective.value <- Inf
    i = 1
    for(i in 1:length(B.list)){
        B.i <- B.list[[i]]
        L.i <- G.list[[i]] %*% t(G.list[[i]])
        objective.value <- norm(B.i - L.i, type = "F") 
        if(objective.value < min.objective.value){
          min.index = i
          min.objective.value <- objective.value
      }
    }

    G = G.list[[min.index]]    

    objective.value <- 0
    for(i in 1:length(B.list)){
        B.i <- B.list[[i]]
        L <- G %*% t(G)
        objective.value <- objective.value + norm(B.i - L, type = "F") 
    }

    return (list(group.space = G, score = objective.value))
}

three.way.mds.projected <- function(B.list, dim.max, dim.list, max.iter, tol = 1e-6){

    P.list <- list()
    B.num <- length(B.list)
    N <- nrow(B.list[[1]])

    ## Initiate each projection matrix be a truncated
    ## identity matrix of rank dim.list[i] 
    
    for( i in 1:B.num){
        P.list[[i]] <- diag(dim.max)
        diag(P.list[[i]]) <- c(seq(1,1,length.out = dim.list[i]),
                               seq(0,0, length.out = dim.max - dim.list[i]))
    }

    ## Initiate the group space to be the cmdscale solution of
    ## the sum of the fallible inner product matrices in B.list

    G.struct <- initialize.mds(B.list, dim.max)
    G <- G.struct$group.space
    objective.value <- G.struct$score

    iter <- 1

    while(TRUE){

        ## Optimize over P.i for each i, keeping G fixed
        for(i in 1:B.num){
            B.i <- B.list[[i]]
            X.i <- t(G) %*% B.i %*% G
            X.i.eigen <- eigen(X.i, only.values = FALSE)
            
            d <- dim.list[i]
            U <- X.i.eigen$vectors
            U.d <- U[,c(1:d)]

            P.i <- U.d %*% t(U.d)
            
            P.list[[i]] <- P.i
        }

        ## Optimize over G, keeping the P.i fixed 
        ## The optimization is done similar to the IDIOSCAL/CANDECOMP algorithm

        ## G.tmp.1 <- matrix(seq(0,0,length.out = N*dim.max),
        ##                   nrow = N, ncol = dim.max)
        ## G.tmp.2 <- matrix(seq(0,0,length.out = dim.max*dim.max),
        ##                   nrow = dim.max, ncol = dim.max)

        ## for(i in 1:B.num){
        ##     Z.i <- G %*% P.list[[i]]
        ##     G.tmp.1 <- G.tmp.1 + B.list[[i]] %*% Z.i
        ##     G.tmp.2 <- G.tmp.2 + t(Z.i) %*% Z.i
        ## }

        ## G.tmp.2.inv <- solve(G.tmp.2)
        ## G.new <- G.tmp.1 %*% G.tmp.2.inv

        G.new <- candecomp(G, B.list, P.list, dim.max)

        objective.value.new <- strain(G.new, B.list, P.list)

        if(abs(objective.value - objective.value.new) < tol || iter >= max.iter)
          break

        iter <- iter + 1
        objective.value <- objective.value.new
        G <- G.new
    }

    return(list(group.space = G, project.list = P.list))
}

strain <- function(G, B.list, P.list){
    
    B.num <- length(B.list)
    objective.value <- 0
    
    for(i in 1:B.num){
            
        P.i <- P.list[[i]]
        B.i <- B.list[[i]]
            
        G.i<- G %*% P.i
        L.i<- G.i%*% t(G.i)
        objective.value <- objective.value + norm(B.i - L.i, type = "F") 
    }

    objective.value
}
    
candecomp <- function(G, B.list, P.list, dim.max, max.iter=5, tol=1e-6){ 
    
    B.num <- length(B.list)
    N <- nrow(B.list[[1]])

    X <- G
    iter <- 1
    while(TRUE){
        old.strain <- strain(X, B.list, P.list)
        
        G.tmp.1 <- matrix(seq(0,0,length.out = N*dim.max),
                          nrow = N, ncol = dim.max)
        G.tmp.2 <- matrix(seq(0,0,length.out = dim.max*dim.max),
                          nrow = dim.max, ncol = dim.max)
        for(i in 1:B.num){
            Z.i <- X %*% P.list[[i]]
            G.tmp.1 <- G.tmp.1 + B.list[[i]] %*% Z.i
            G.tmp.2 <- G.tmp.2 + t(Z.i) %*% Z.i
        }

        G.tmp.2.inv <- ginv(G.tmp.2)
        Y <- G.tmp.1 %*% G.tmp.2.inv


        G.tmp.3 <- matrix(seq(0,0,length.out = N*dim.max),
                          nrow = N, ncol = dim.max)
        G.tmp.4 <- matrix(seq(0,0,length.out = dim.max*dim.max),
                          nrow = dim.max, ncol = dim.max)

        for(i in 1:B.num){
            Z.i <- Y %*% P.list[[i]]
            G.tmp.3 <- G.tmp.3 + B.list[[i]] %*% Z.i
            G.tmp.4 <- G.tmp.4 + t(Z.i) %*% Z.i
        }
        
        G.tmp.4.inv <- ginv(G.tmp.4)
        X <- G.tmp.3 %*% G.tmp.4.inv

        new.strain <- strain(X, B.list, P.list)

        if(abs(new.strain - old.strain) < tol)
          break

        if(iter >= max.iter)
          break

        iter <- iter + 1
    }

    G.new <- X
    G.new
}

        

        
