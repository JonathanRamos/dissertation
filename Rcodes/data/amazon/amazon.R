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

X <- three.way.mds.projected(B.list, 3, dim.list = c(2,2), weights.list <- c(1,10), max.iter = 500, tol = 1e-6)
coord1 <- X$group.space %*% X$project.list[[1]]
coord2 <- X$group.space %*% X$project.list[[2]]

##coords <- coord2
coords <- rbind(coord1, coord2)
pch.list <- c(seq(17,17,length.out = nrow(coord1)), seq(19,19,length.out = nrow(coord2)))
colors <- c(seq("red","red", length.out = nrow(coord1)), seq("blue", "blue", length.out = nrow(coord2)))

tmp1 <- seq(65, by=1, length.out = 26)
tmp2 <- sapply(tmp1, intToUtf8)
tmp3 <- seq(1, by = 1, length.out = 9)
labels <- c(tmp2, as.character(tmp3), as.character(0))
#labels <- c(seq(1,by=1,length.out = nrow(coord1)),seq(1,by=1,length.out = nrow(coord1)))

pdf("dge_morsecode_bimension12.pdf", useDingbats = FALSE)
plot(coords[,1], coords[,2], col = colors, pch = pch.list, asp = 1)
text(coords[,1], coords[,2], pos = 3, labels)

dev.off()

pdf("dge_morsecode_bimension13.pdf", useDingbats = FALSE)
plot(coords[,1], coords[,3], col = colors, pch = pch.list, asp = 1)
text(coords[,1], coords[,3], pos = 3, labels)

dev.off()

pdf("dge_morsecode_bimension23.pdf", useDingbats = FALSE)
plot(coords[,2], coords[,3], col = colors, pch = pch.list, asp = 1)
text(coords[,2], coords[,3], pos = 3, labels)
dev.off()
    
dev.new()
xxx.coord <- rbind(xxx1,xxx2)
plot(xxx.coord, col = colors, pch = pch.list, asp = 1)
text(xxx.coord[,1], xxx.coord[,2], pos = 3, labels)
    

x <- matrix(rnorm(50), nrow = 25)
dist1 <- as.matrix(dist(x))
dist2 <- as.matrix(dist(x[,1]))

B1 <-  tau(dist1^2)
B2 <-  tau(dist2^2)

B.list.other <- list(B1,B2)

X <- three.way.mds.projected(B.list.other, 2, dim.list = c(2,1), weights.list <- c(1,1), max.iter = 500, tol = 1e-6)

coord1 <- X$group.space %*% X$project.list[[1]]
coord2 <- X$group.space %*% X$project.list[[2]]

coords <- rbind(coord1, coord2)
pch.list <- c(seq(17,17,length.out = nrow(coord1)), seq(19,19,length.out = nrow(coord2)))
colors <- c(seq("red","red", length.out = nrow(coord1)), seq("blue", "blue", length.out = nrow(coord2)))

labels <- seq(1, by = 1, length.out = nrow(coord1))
labels <- c(labels, labels)

pdf("dge_toy1_projected.pdf", useDingbats = FALSE)
plot(coords[,1], coords[,2], xlab = "", ylab = "", col = colors, pch = pch.list, asp = 1)
text(coords[,1], coords[,2], pos = 3, labels)
dev.off()

zzz1 <- cmdscale(dist1, 2)
zzz2 <- cmdscale(dist2, 2)
zzz <- rbind(zzz1, zzz2)

pdf("dge_toy1_separate.pdf", useDingbats = FALSE)
plot(zzz[,1], zzz[,2], col = colors, xlab = "", ylab = "", pch = pch.list, asp = 1)
text(zzz[,1], zzz[,2], pos = 3, labels)
dev.off()

lll.upper <- S.upper(lll)
lll.lower <- S.lower(lll)
lll.symmetric <- (lll + t(lll))/2
lll.difference <- abs(lll.upper - lll.lower)/2
B.upper <- tau(lll.upper^2)
B.lower <- tau(lll.lower^2)
B.symmetric <- tau(lll.symmetric)
B.difference <- tau(lll.difference^2)
B.symmetric <- B.symmetric/norm(B.symmetric, type = "F")*sqrt(2)
B.difference <- B.difference/norm(B.difference, type = "F")*sqrt(2)
lll.list <- list(B.symmetric, B.difference)
X <- three.way.mds.projected(lll.list, 3, dim.list = c(2,2), weights.list <- c(1,1), max.iter = 500, tol = 1e-6)

coord1 <- X$group.space %*% X$project.list[[1]]
coord2 <- X$group.space %*% X$project.list[[2]]

coords <- rbind(coord1, coord2)
pch.list <- c(seq(17,17,length.out = nrow(coord1)), seq(19,19,length.out = nrow(coord2)))
colors <- c(seq("red","red", length.out = nrow(coord1)), seq("blue", "blue", length.out = nrow(coord2)))

labels <- seq(1, by = 1, length.out = nrow(coord1))
labels <- c(labels, labels)

pdf("cola_switching.pdf", useDingbats = FALSE, width=10)
plot(coords[,1], coords[,2], xlab = "", ylab = "", col = colors, pch = pch.list, asp = 1)
smartlegend(x="left", y = "top", inset = 0,
            c("1: Coke decaf", "2: Coke diet decaf", "3: Pepsi diet decaf",
              "4: Pepsi decaf", "5: Canfield", "6: Coke", "7: Coke classic",
              "8: Coke diet", "9: Pepsi diet", "10: RC diet", "11: Rite diet",
              "12: Pepsi", "13: Private", "14: RC", "15: Wildwood"))
text(coords[,1], coords[,2], pos = 3, labels)
dev.off()
dev.new()
plot(coords[,1], coords[,3], xlab = "", ylab = "", col = colors, pch = pch.list, asp = 1)
text(coords[,1], coords[,3], pos = 3, labels)
dev.new()
plot(coords[,2], coords[,3], xlab = "", ylab = "", col = colors, pch = pch.list, asp = 1)
text(coords[,2], coords[,3], pos = 3, labels)
