data.setup <- function(){
    spam.data <- read.csv("spambase.data", header = FALSE, sep=",")
    spam.data <- as.matrix(spam.data)
    spam.sd <- sd(spam.data)
    spam.mean <- apply(spam.data, 2, mean)

    spam.normalized = matrix(nrow = nrow(spam.data), ncol = ncol(spam.data))

    for(i in 1:nrow(spam.data)){
        spam.normalized[i,] = (spam.data[i,] - spam.mean)/spam.sd
    }

#    spam.normalized[,55:57] = spam.normalized[,55:57]*0.1

    labels = spam.data[,ncol(spam.data)]

    return(list(l = labels, d = spam.normalized[,-ncol(spam.normalized)]))
}

spam.similarities <- function(dat, k){

    dist.mat <- as.matrix(dist(dat))
    W <- matrix(nrow = nrow(dist.mat), ncol = ncol(dist.mat))
    
    for(i in 1:nrow(dist.mat)){
        W[i,] <- c(1/(200*k))
        tmp <- sort(dist.mat[i,], index.return = TRUE)
        index = tmp$ix[2:(k+1)]
        vals = 1/seq(1:k)
        W[i,index] = vals
    }

    diag(W) <- 0

    return(W)
}
        
small.spam <- function(pct, k){
    spam.frame <- data.setup()
    spam.l <- spam.frame$l
    spam.d <- spam.frame$d
    
    u <- runif(nrow(spam.d))
    data <- spam.d[u < pct,]
    labels <- spam.l[u < pct] 

    W <- spam.similarities(data, k)
    P <- transition.matrix(W)
    spam.ect <- ect(P)
    map.ect <- cmdscale(spam.ect, k = 3)
    spheres3d(map.ect[,1], map.ect[,2], map.ect[,3], col = labels + 2)
#    map.ect <- ect.map(P)
#    map.diff <- diffusion.map(P, t = 2, 3)
#    plot(map.ect[,1],map.ect[,2],col = labels + 2)
#    plot(map.diff[,1],map.diff[,2],col = labels + 2)
#    spheres3d(map.diff[,1], map.diff[,2], map.diff[,3], radius = 1, col = labels + 2)
}

    
