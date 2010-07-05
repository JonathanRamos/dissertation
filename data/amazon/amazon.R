## Read the Amazon data file
read.dat <- function(fname){
    ff <- file(fname, open="rt")

    ## Read Header
    l1 <- readLines(ff, n = 1)
    ## N = number of objects
    N <- strtoi(strsplit(l1, '\\t')[[1]][1])

    S <- matrix(seq(0,0,length.out = N*N), nrow = N, ncol = N)

    ## Skip description of data
    tmp <- readLines(ff, n = 1)

    while(TRUE){
        line <- readLines(ff, n = 1)
        if(length(line) == 0){
            break
        }
        tmp1 <- strsplit(line, '\\t')
        tmp2 <- tmp1[[1]][1]
        tmp3 <- strsplit(tmp2, ' ')
        i <- strtoi(tmp3[[1]][1]) + 1
        j <- strtoi(tmp3[[1]][2]) + 1

        s <- as.numeric(tmp1[[1]][2])

        S[i,j] <- s
    }
    close(ff)
    S
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


    
    
