library("R.matlab")
require("gplots")
cutoff = 0.03

mnist.dat <- readMat("mnist_all.mat")
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

digits.small <- rbind(one.train.small,two.train.small,three.train.small,four.train.small,five.train.small,six.train.small,seven.train.small,eight.train.small, nine.train.small,zero.train.small)

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

digits.small <- rbind(one.train.small,two.train.small,three.train.small,four.train.small,five.train.small,six.train.small,seven.train.small,eight.train.small, nine.train.small,zero.train.small)
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

#labels <- c(seq(1,1,length.out = nrow(one.train.small)), seq(5,5,length.out = nrow(five.train.small)), seq(7,7,length.out = nrow(seven.train.small)), seq(8,8,length.out = nrow(eight.train.small)))

epsilon = 400000

W <- gaussian.similarity(digits.small, epsilon)
W.sparse <- sparsify.similarity(W, 20)
P <- transition.matrix(W.sparse)
digits.map <- ect.map(P)

plot(digits.map, col = labels, xlab="", ylab="", pch = 21)


W <- gaussian.similarity(four.five.small, epsilon)
P <- transition.matrix(W)
digits.map <- ect.map(P)

pdf("mnist45_small.pdf", useDingbats = FALSE)
plot(digits.map, col = fo
> ur.five.labels,xlab="",ylab="",pch=21)

smartlegend(x = "left", y="top", inset = 0,
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



