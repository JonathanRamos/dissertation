data.setup <- function(){
    abalone.data <- read.csv("abalone.data", header = FALSE, sep = ",")
    sex <- abalone.data[,1]

    abalone.mat <- as.matrix(abalone.data[,-c(1,3,6,7,9)])

    abalone.male <- abalone.mat[sex == "M",]
    abalone.female <- abalone.mat[sex == "F",]
    abalone.infant <- abalone.mat[sex == "I",]

    male.mean <- apply(abalone.male, 2, mean)
    male.sd <- sd(abalone.male)

    female.mean <- apply(abalone.female, 2, mean)
    female.sd <- sd(abalone.female)
    
    infant.mean <- apply(abalone.infant, 2, mean)
    infant.sd <- sd(abalone.infant)

    male.mean.mat <- outer(seq(1,1,length.out = nrow(abalone.male)), male.mean)
    female.mean.mat <- outer(seq(1,1,length.out = nrow(abalone.female)), female.mean)
    infant.mean.mat <- outer(seq(1,1,length.out = nrow(abalone.infant)), infant.mean)
    
    male.sd.mat <- outer(seq(1,1,length.out = nrow(abalone.male)), male.sd)
    female.sd.mat <- outer(seq(1,1,length.out = nrow(abalone.female)), female.sd)
    infant.sd.mat <- outer(seq(1,1,length.out = nrow(abalone.infant)), infant.sd)

    male.mat <- (abalone.male - male.mean.mat)/male.sd.mat
    female.mat <- (abalone.female - female.mean.mat)/female.sd.mat
    infant.mat <- (abalone.infant - infant.mean.mat)/infant.sd.mat

    male.labels = abalone.data[sex == "M",9]
    female.labels = abalone.data[sex == "F",9]
    infant.labels = abalone.data[sex == "I",9]

    return(list(male=list(data = male.mat, labels = male.labels), female = list(data = female.mat, labels = female.labels), infant = list(data = infant.mat, labels = infant.labels)))
}

## abalone.struct <- data.setup()

## abalone.male <- abalone.struct$data[abalone.struct$sex == "M",]
## abalone.female <- abalone.struct$data[abalone.struct$sex == "F",]
## abalone.infant <- abalone.struct$data[abalone.struct$sex == "I",]
## labels.male <- abalone.struct$labels[abalone.struct$sex == "M"]
## labels.female <- abalone.struct$labels[abalone.struct$sex == "F"]
## labels.infant <- abalone.struct$labels[abalone.struct$sex == "I"]

## abalone.male.small <- abalone.male[1:500,]
## labels.male.small <- labels.male[1:500]
## W <- gaussian.similarity(abalone.male.small, epsilon = 1)
## W.sparse <- sparsify.similarity(W, 50)
## vol <- sum(W.sparse)
## P <- transition.matrix(W.sparse)
## abalone.male.ect <- ect(P)
## #abalone.male.diff <- diffusion.distance(P, t = 20)
## abalone.male.map <- cmdscale(abalone.male.ect/sqrt(vol),k=3)
## #plot(abalone.male.map[,1], abalone.male.map[,2], col = labels.male.small)
## spheres3d(abalone.male.map[,1], abalone.male.map[,2], abalone.male.map[3,], col = labels.male.small)

abalone.infant.fun <- function(n,k,epsilon,d){
    abalone.struct <- data.setup()
    data.infant <- abalone.struct$infant$data
    labels.infant <- abalone.struct$infant$labels
    data.infant.small <- data.infant[1:n,]
    labels.infant.small <- labels.infant[1:n]
    W <- gaussian.similarity(data.infant.small, epsilon)
    P <- transition.matrix(W)
    infant.ect <- ect(P)
#abalone.infant.ect[is.na(abalone.infant.ect)] <- 100000
#abalone.male.diff <- diffusion.distance(P, t = 20)
    infant.map <- cmdscale(infant.ect, 2)
#    infant.map <- laplacian.map(W.sparse,d)
    plot(infant.map[,1], infant.map[,2], col = floor(labels.infant.small)/4 + 2)
#spheres3d(infant.map[,1], infant.map[,2], infant.map[,3], col = labels.infant.small)
}


    
    
