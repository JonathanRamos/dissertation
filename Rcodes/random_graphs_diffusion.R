require("igraph")
source("diffusion_distances.R")

random.graph <- function(data, p){
    n <- nrow(data)
    W <- zeros(n,n)

    for(i in 1:n){
        sample <- runif(n)
        W[i,] = (sample <= p)
    }

    diag(W) <- 0
    ##    B <- graph.adjacency(W,mode="undirected",diag=FALSE)
    ##    plot.igraph(B, layout=layout.fruchterman.reingold,vertex.size=10)
    return(W)
}

diffusion.random.graph <- function(n,p,t,iterations){

    for(i in 1:iterations){
        G <- erdos.renyi.game(n,p)
        W <- get.adjacency(G,type="both")
    
        diffusion.struct <- diffusion.distance(W, t)
        dist <- as.matrix(diffusion.struct$dist)
    
        kNN <- 5

        for(j in 1:n){
            dist.j <- dist[j,]
            tmp <- sort(dist.j, index.return = TRUE)
            dist[j,] <- 0
            dist[j,tmp$ix[1:kNN]] <- 1
        }

        A <- graph.adjacency(W, mode="undirected",weighted=NULL, diag=FALSE) 
        B <- graph.adjacency(dist,mode="undirected",weighted=NULL,diag=FALSE)

        B.all.pairs.shortest.path <- shortest.paths(B,mode="all")
        inf_idx <- find(B.all.pairs.shortest.path == Inf)
        B.all.pairs.shortest.path[inf_idx] <- 10000.0
        idx <- find(W > 0)
        tmp <- table(B.all.pairs.shortest.path[idx])
        if(i == 1){
            tsukue <- tmp
            ## all <- B.all.pairs.shortest.path
        }
        else{
            z <- c(tsukue,tmp)
            ## all <- c(all,B.all.pairs.shortest.path)
            tsukue <- tapply(z, names(z), sum)
        }
    }
    ## par(mfrow=c(1,1))
    ## plot.igraph(A,layout=layout.fruchterman.reingold,vertex.size=10)
    ## plot.igraph(B,layout=layout.fruchterman.reingold,vertex.size=10)
    ## plot(density(B.all.pairs.shortest.path[idx]),\
    ## main="Shortest path distance between original k-NN pairs",col="red")
    ## plot(density(as.data.frame(tsukue),
    ## main="Shortest path distance between original k-NN pairs", col = "green"))
    ## hist(all,main=sprintf("Distribution of shortest path distance between original k-NN pairs at time scale %d",t),freq=FALSE)
    ## lines(density(tsukue),col="red")
    barplot(tsukue,main=sprintf("Distribution of shortest path distance between original k-NN pairs at time scale %d",t))
    ## return(B)
}

## num.clust <-  20
## X <- 10*rnorm(50)
## Y <- 10*rnorm(50)
## data <- cbind(X,Y)
## labels <- seq(1,by=1,length.out = 50)

## W <- gaussian.similarity(data, epsilon = 4.0)
## length(find(W > 0.05))
## diffusion <- diffusion.distance(W, t = 5, k = 2, sparsity = FALSE, cutoff = 0)
## d.diffusion = diffusion$dist
## map <- cmdscale(d.diffusion, k = 2)
## ## tree.diffusion = hclust(as.dist(d.diffusion))
## ## cluster.diffusion = cutree(tree.diffusion, k = num.clust)
## ## data.table <- as.table(cbind(data, cluster.diffusion))
## clusplot(data, cluster.diffusion, color = TRUE, shade = FALSE, labels = 0, lines = 0)
## par(mfrow = c(1,2))
## plot(X,Y,pch=19,col=labels)
## #plot(X,Y, col = data.table[,3],pch=19, xlab = "X", ylab = "Y")
## plot(map[,1],map[,2],pch=19,col=labels)
## #text(map[,1],map[,2] + 0.05,labels)

## data <- clusters2.data(n = 40)
## df <- data.frame(x=data[,1],y=data[,2])
## W <- gaussian.similarity(data, epsilon = 0.57)
## diffusion.struct <- diffusion.distance(W, t = 5)
## df1 <- data.frame(x1=diffusion.struct$map[,2],y1=diffusion.struct$map[,3])

## W <- gaussian.similarity(data, epsilon = 0.58)
## diffusion.struct <- diffusion.distance(W, t = 5)
## df2 <- data.frame(x2=diffusion.struct$map[,2],y2=diffusion.struct$map[,3])

## W <- gaussian.similarity(data, epsilon = 0.60)
## diffusion.struct <- diffusion.distance(W, t = 5)
## df3 <- data.frame(x3=diffusion.struct$map[,2],y3=diffusion.struct$map[,3])

## W <- gaussian.similarity(data, epsilon = 0.62)
## diffusion.struct <- diffusion.distance(W, t = 5)
## df4 <- data.frame(x4=diffusion.struct$map[,2],y4=diffusion.struct$map[,3])

## W <- gaussian.similarity(data, epsilon = 0.64)
## diffusion.struct <- diffusion.distance(W, t = 5)
## df5 <- data.frame(x5=diffusion.struct$map[,2],y5=diffusion.struct$map[,3])

## dframe <- cbind(df,df1,df2,df3,df4,df5)

## g <- ggobi(dframe)
