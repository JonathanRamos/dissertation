library("R.matlab")
require("gplots")
mnist.dat <- readMat("mnist_all.mat")

cutoff = 0.1
one.train.small <- mnist.dat$train1[runif(dim(mnist.dat$train1)[1]) < cutoff,]
two.train.small <- mnist.dat$train2[runif(dim(mnist.dat$train2)[1]) < cutoff,]
three.train.small <- mnist.dat$train3[runif(dim(mnist.dat$train3)[1]) < cutoff,]
four.train.small <- mnist.dat$train4[runif(dim(mnist.dat$train4)[1]) < cutoff,]
five.train.small <- mnist.dat$train5[runif(dim(mnist.dat$train5)[1]) < cutoff,]
six.train.small <- mnist.dat$train6[runif(dim(mnist.dat$train6)[1]) < cutoff,]
seven.train.small <- mnist.dat$train7[runif(dim(mnist.dat$train7)[1]) < cutoff,]
eight.train.small <- mnist.dat$train8[runif(dim(mnist.dat$train8)[1]) < cutoff,]
nine.train.small <- mnist.dat$train9[runif(dim(mnist.dat$train9)[1]) < cutoff,]
zero.train.small <- mnist.dat$train0[runif(dim(mnist.dat$train0)[1]) < cutoff,]

digits.small <- rbind(one.train.small,two.train.small,three.train.small,
                      four.train.small, five.train.small,six.train.small,
                      seven.train.small,eight.train.small,
                      nine.train.small,zero.train.small)

zero.eight.small <- rbind(zero.train.small,eight.train.small)
three.nine.small <- rbind(three.train.small,nine.train.small)
zero.one.small <-  rbind(zero.train.small, one.train.small)
four.five.small <- rbind(four.train.small, five.train.small)
two.six.small <- rbind(two.train.small, six.train.small)
one.seven.small <- rbind(one.train.small, seven.train.small)

zero.eight.labels <- c(seq(10,10,length.out = nrow(zero.train.small)),
                       seq(8,8,length.out = nrow(eight.train.small)))

three.nine.labels <- c(seq(3,3,length.out = nrow(three.train.small)),
                       seq(9,9,length.out = nrow(nine.train.small)))
zero.one.labels <- c(seq(10,10,length.out = nrow(zero.train.small)),
                     seq(1,1,length.out = nrow(one.train.small)))

four.five.labels <- c(seq(4,4,length.out = nrow(four.train.small)),
                      seq(5,5,length.out = nrow(five.train.small)))
two.six.labels <- c(seq(2,2,length.out = nrow(two.train.small)),
                    seq(6,6,length.out = nrow(six.train.small)))

one.seven.labels <- c(seq(1,1,length.out = nrow(one.train.small)),
                      seq(7,7,length.out = nrow(seven.train.small)))

digits.small <- rbind(one.train.small,two.train.small,three.train.small,four.train.small,five.train.small,six.train.small,
                      seven.train.small,eight.train.small, nine.train.small,zero.train.small)
labels <- c(seq("red","red",length.out = nrow(one.train.small)),
            seq("green","green",length.out = nrow(two.train.small)),
            seq("cyan","cyan",length.out = nrow(three.train.small)),
            seq("blue","blue",length.out = nrow(four.train.small)),
            seq("orange","orange",length.out = nrow(five.train.small)),
            seq("magenta","magenta",length.out = nrow(six.train.small)),
            seq("gray","gray",length.out = nrow(seven.train.small)),
            seq("yellow","yellow",length.out = nrow(eight.train.small)),
            seq("brown","brown",length.out = nrow(nine.train.small)),
            seq("violet","violet",length.out = nrow(zero.train.small)))
#digits.small <- rbind(one.train.small,five.train.small,seven.train.small,eight.train.small)

#labels <- c(seq(1,1,length.out = nrow(one.train.small)),
#seq(5,5,length.out = nrow(five.train.small)),
#seq(7,7,length.out = nrow(seven.train.small)), seq(8,8,length.out = nrow(eight.train.small)))

epsilon = 500000

W <- gaussian.similarity(digits.small, epsilon)
W.sparse <- sparsify.similarity(W, 50)
P <- transition.matrix(W.sparse)
digits.map <- ect.map(P)

plot(digits.map, col = labels, xlab="", ylab="", pch = 21)


W <- gaussian.similarity(four.five.small, epsilon)
P <- transition.matrix(W)

digits.map.ect <- ect.map(P)
pdf("mnist45_small_ect_map.pdf", useDingbats = FALSE)
plot(digits.map.ect, col = four.five.labels,xlab="",ylab="",pch=21)
smartlegend(x = "left", y="top", inset = 0,
            c("Digit 4", "Digit 5"), fill = c(4,5))
dev.off()

digits.map.diffusion <- diffusion.map(P, t=10)
pdf("mnist45_small_diffusion_map.pdf", useDingbats = FALSE)
plot(-digits.map.diffusion[,1], digits.map.diffusion[,2], col = four.five.labels,xlab="",ylab="",pch=21)

smartlegend(x = "left", y="top", inset = 0,
            c("Digit 4", "Digit 5"), fill = c(4,5))
dev.off()

W <- gaussian.similarity(100000)
W.sparse <- sparsify.similarity.directed(W,5)
P <- transition.matrix(W.sparse)
dist.ect <- ect(P)
ect.cmds <- cmdscale(dist.ect, k = 2)
#pdf("mnist45_small_ect_cmds.pdf", useDingbats = FALSE)
plot(ect.cmds, col = four.five.labels,xlab="",ylab="",pch=21)
smartlegend(x = "right", y="bottom", inset = 0,
            c("Digit 4", "Digit 5"), fill = c(4,5))
dev.off()


W <- gaussian.similarity(zero.one.small, epsilon)
P <- transition.matrix(W)
digits.map <- ect.map(P)
pdf("mnist01_small.pdf", useDingbats = FALSE)
plot(digits.map, col = zero.one.labels,xlab="",ylab="",pch=21)
smartlegend(x = "left", y="top", inset = 0,
            c("Digit 0", "Digit 1"), fill = c(10,1))
dev.off()

W <- gaussian.similarity(zero.eight.small, epsilon)
P <- transition.matrix(W)
digits.map <- ect.map(P)

pdf("mnist08_small.pdf", useDingbats = FALSE)
plot(digits.map, col = zero.eight.labels,xlab="",ylab="",pch=21)
smartlegend(x = "left", y="top", inset = 0,
            c("Digit 0", "Digit 8"), fill = c(10,8))
dev.off()

W <- gaussian.similarity(two.six.small, epsilon)
P <- transition.matrix(W)
digits.map <- ect.map(P)

pdf("mnist26_small.pdf", useDingbats = FALSE)
plot(digits.map, col = two.six.labels,xlab="",ylab="",pch=21)
smartlegend(x = "left", y="top", inset = 0,
            c("Digit 2", "Digit 6"), fill = c(2,6))
dev.off()

W <- gaussian.similarity(three.nine.small, epsilon)
P <- transition.matrix(W)
digits.map <- ect.map(P)

pdf("mnist39_small.pdf", useDingbats = FALSE)
plot(digits.map, col = three.nine.labels,xlab="",ylab="",pch=21)
smartlegend(x = "left", y="top", inset = 0,
            c("Digit 3", "Digit 9"), fill = c(3,9))
dev.off()


W <- gaussian.similarity(one.seven.small, epsilon)
P <- transition.matrix(W)
digits.map <- ect.map(P)

pdf("mnist17_small.pdf", useDingbats = FALSE)
plot(digits.map, col = one.seven.labels,xlab="",ylab="",pch=21)
smartlegend(x = "left", y="top", inset = 0,
            c("Digit 1", "Digit 7"), fill = c(1,7))
dev.off()

four.test.small = mnist.dat$test4[runif(nrow(mnist.dat$test4)) < cutoff,]
five.test.small = mnist.dat$test5[runif(nrow(mnist.dat$test5)) < cutoff,]

four.n = nrow(four.test.small)
five.n = nrow(five.test.small)

four.k = matrix(nrow = four.n, ncol = nrow(four.five.small))
five.k = matrix(nrow = five.n, ncol = nrow(four.five.small))

for(i in (1:four.n)){
    x = outer(seq(1,1,length.out = nrow(four.five.small)), four.test.small[i,])
    l = apply((x - four.five.small)^2, 1, sum)
    four.k[i,] = exp(-l/epsilon)
}
    
for(i in (1:five.n)){
    x = outer(seq(1,1,length.out = nrow(four.five.small)), five.test.small[i,])
    l = apply((x - four.five.small)^2, 1, sum)
    five.k[i,] = exp(-l/epsilon)
}

P.four <- transition.matrix(four.k)
P.five <- transition.matrix(five.k)

four.kern = matrix(nrow = four.n, ncol = nrow(four.five.small))
five.kern = matrix(nrow = five.n, ncol = nrow(four.five.small))
for(i in (1:four.n)){
    four.kern[i,] = out.of.sample.ect(P,P.four[i,])
}
    
for(i in (1:five.n)){
    five.kern[i,] = out.of.sample.ect(P,P.five[i,])
}

four.map = matrix(nrow = four.n, ncol = ncol(digits.map))
five.map = matrix(nrow = five.n, ncol = ncol(digits.map))

for(i in (1:four.n)){
    four.map[i,1] = sum(four.kern[i,] * digits.map[,1])/sum(four.kern[i,])
    four.map[i,2] = sum(four.kern[i,] * digits.map[,2])/sum(four.kern[i,])
}

for(i in (1:five.n)){
    five.map[i,1] = sum(five.kern[i,] * digits.map[,1])/sum(five.kern[i,])
    five.map[i,2] = sum(five.kern[i,] * digits.map[,2])/sum(five.kern[i,])
}

new.points = rbind(four.map, five.map)
new.labels = c(seq(10,10,length.out = four.n), seq(9,9,length.out = five.n))

plot(digits.map, col = four.five.labels,xlab="",ylab="",pch=21)
points(new.points, pch = 17, cex=1, col = new.labels)

X <- prcomp(four.five.small)
pdf("four_five_pca.pdf", useDingbats = FALSE)
plot(X$x, col = four.five.labels)
smartlegend(x = "left", y="top", inset = 0,
            c("Digit 4", "Digit 5"), fill = c(4,5))
dev.off()

X <- prcomp(zero.eight.small)
pdf("zero_eight_pca.pdf", useDingbats = FALSE)
plot(X$x, col = zero.eight.labels)
dev.off()

X <- prcomp(one.seven.small)
pdf("one_seven_pca.pdf", useDingbats = FALSE)
plot(X$x, col = one.seven.labels)
dev.off()

X <- prcomp(three.nine.small)
pdf("three_nine_pca.pdf", useDingbats = FALSE)
plot(X$x, col = three.nine.labels)
dev.off()

X <- prcomp(two.six.small)
pdf("two_six_pca.pdf", useDingbats = FALSE)
plot(X$x, col = two.six.labels)
dev.off()

X <- prcomp(zero.one.small)
pdf("zero_one_pca.pdf", useDingbats = FALSE)
plot(X$x, col = zero.one.labels)
dev.off()
