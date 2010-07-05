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

    break.loop <- FALSE
    iter <- 1

    while( break.loop == FALSE){

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
        ## The optimization is done similar to the IDIOSCAL algorithm

        G.tmp.1 <- matrix(seq(0,0,length.out = N*dim.max),
                          nrow = N, ncol = dim.max)
        G.tmp.2 <- matrix(seq(0,0,length.out = dim.max*dim.max),
                          nrow = dim.max, ncol = dim.max)

        for(i in 1:B.num){
            Z.i <- G %*% P.list[[i]]
            G.tmp.1 <- G.tmp.1 + B.list[[i]] %*% Z.i
            G.tmp.2 <- G.tmp.2 + t(Z.i) %*% Z.i
        }

        G.tmp.2.inv <- solve(G.tmp.2)
        G.new <- G.tmp.1 %*% G.tmp.2.inv

        objective.value.new = 0
        for(i in 1:B.num){
            
            P.i <- P.list[[i]]
            B.i <- B.list[[i]]
            
            G.i.new <- G.new %*% P.i
            L.i.new <- G.i.new %*% t(G.i.new)
            objective.value.new <- objective.value + norm(B.i - L.i.new, type = "F") 
        }

        if(abs(objective.value - objective.value.new) < tol || iter >= max.iter)
          break.loop = TRUE

        iter <- iter + 1
        objective.value <- objective.value.new
        G <- G.new
    }

    return(list(group.space = G, project.list = P.list))
}

        
