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

## Is it important to normalize the B.i in B.list to have
## Frobenius norm = dim.max ?

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

    G <- G.list[[min.index]]    

    G
}

three.way.mds.projected <- function(B.list, dim.max, dim.list, weights.list, max.iter, tol = 1e-6){

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

    G <- initialize.mds(B.list, dim.max)
    objective.value <- strain(G, B.list, P.list, weights.list) 
    print(objective.value)

    sigma.vals <- sigma.values(B.list)

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

        ## G.new <- candecomp(G, B.list, P.list, dim.max)

        G.new <- majorize.G(G, B.list, P.list, weights.list, sigma.vals, max.iter = 20, tol = 1e-6)

        objective.value.new <- strain(G.new, B.list, P.list, weights.list)
        print(objective.value.new)

        if(abs(objective.value - objective.value.new) < tol || iter >= max.iter)
          break

        iter <- iter + 1
        objective.value <- objective.value.new
        G <- G.new
    }

    return(list(group.space = G, project.list = P.list))
}

strain <- function(G, B.list, P.list, weights.list){
    
    B.num <- length(B.list)
    objective.value <- 0
    
    for(i in 1:B.num){
            
        P.i <- P.list[[i]]
        B.i <- B.list[[i]]
            
        G.i<- G %*% P.i
        L.i<- G.i%*% t(G.i)
        objective.value <- objective.value + weights.list[i]*norm(B.i - L.i, type = "F") 
    }

    objective.value
}
    
candecomp <- function(G, B.list, P.list, dim.max, max.iter=5, tol=1e-6){ 
    
    B.num <- length(B.list)
    N <- nrow(B.list[[1]])

    X <- G
    iter <- 1
    while(TRUE){
        old.strain <- strain(X, B.list, P.list, weights.list)
        
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

        new.strain <- strain(X, B.list, P.list, weights.list)

        if(abs(new.strain - old.strain) < tol)
          break

        if(iter >= max.iter)
          break

        iter <- iter + 1
    }

    G.new <- X
    G.new
}

sigma.values <- function(B.list){

  B.num <- length(B.list)

  sigma.vals <- seq(0,0,length.out = B.num)

  for(i in 1:B.num){
    tmp <- svd(B.list[[i]], nu = 0, nv = 0)
    sigma.vals[i] <- tmp$d[1]
  }

  sigma.vals
}

majorize.G <- function(G, B.list, P.list, weights.list, sigma.vals, max.iter, tol=1e-6){

  iter <- 1
  C <- 2/sum(sigma.vals)
  B.num <- length(B.list)

##  while(TRUE){
    G.tmp <- matrix(seq(0,0,length.out = nrow(G)*ncol(G)),
                          nrow = nrow(G), ncol = ncol(G))

    for(i in 1:B.num){
      G.tmp <- G.tmp + weights.list[i]*B.list[[i]]%*%G%*%P.list[[i]]
    }

    G.tmp <- G + C*G.tmp

    G.tmp.svd <- svd(G.tmp)

    U <- G.tmp.svd$u
    V <- G.tmp.svd$v

    G.new <- U %*% t(V)

  ## B.list.norm <- seq(0,0,length.out = B.num)
  
  ## for(i in 1:B.num){
  ##      B.list.norm[i] <- norm(B.list[[i]], type = "F")
  ##  }
  
  ##    G.new <- G.new * sqrt(max(B.list.norm))/sqrt(sqrt(ncol(G.new)))
  
  ##     objective.value <- 0
  ##     for(i in 1:B.num){
  ##       objective.value <- objective.value + (-2)*sum(diag(B.list[[i]] %*% G %*% P.list[[i]] %*% t(G)))
  ##     }
  ## ##    print(objective.value)

  ##     objective.value.new <- 0
  ##     for(i in 1:B.num){
  ##       objective.value.new <- objective.value.new + (-2)*sum(diag(B.list[[i]] %*% G.new %*% P.list[[i]] %*% t(G.new)))
  ##     }

  ## G <- G.new
    
  ##    print(objective.value.new)

  ## if(iter >= max.iter)
  ##   break
  ## if(abs(objective.value - objective.value.new) < tol)
  ##   break

  ## iter <- iter + 1
  ##  }

  G.new
}

upper.triangular <- function(N){
    U <- diag(N)
    for(i in 1:N){
        for(j in i:N){
            U[i,j] <- 1
        }
    }
    U
}

lower.triangular <- function(N){
    L <- diag(N)
    for(i in 1:N){
        for(j in 1:i){
            L[i,j] <- 1
        }
    }
    diag(L) <- 0
    L
}

S.upper <- function(S){
    N <- nrow(S)
    U <- upper.triangular(N)
    L <- lower.triangular(N)
    S.u <- S*U + t(S)*L
}

S.lower <- function(S){
    
    N <- nrow(S)
    U <- upper.triangular(N)
    L <- lower.triangular(N)
    S.l <- S*L + t(S)*U
}
