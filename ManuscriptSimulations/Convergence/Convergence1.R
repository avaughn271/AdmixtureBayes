library("coda") # We load the coda library

#We format the plots
par(mfrow=c(3,3))
 # bottom  left   top right
par(mar = c(4.2, 4, 1.5, 1))   # changed bottom from 2

#We assume that the output of the 3 independently run chains are
#chain1.txt, chain2.txt, and chain3.txt
Chain1 = read.csv("chain1.csv")
Chain2 = read.csv("chain2.csv")
Chain3 = read.csv("chain3.csv")

#We extract the total branch length
Branches1 = Chain1[,which(names(Chain1) == "total_branch_length")]
Branches2 = Chain2[,which(names(Chain2) == "total_branch_length")]
Branches3 = Chain3[,which(names(Chain3) == "total_branch_length")]

#We do the same thing for the posterior
Post1 = Chain1[,which(names(Chain1) == "posterior")]
Post2 = Chain2[,which(names(Chain2) == "posterior")]
Post3 = Chain3[,which(names(Chain3) == "posterior")]

#And for the number of admixture events
admixes1 = Chain1[,which(names(Chain1) == "no_admixes")]
admixes2 = Chain2[,which(names(Chain2) == "no_admixes")]
admixes3 = Chain3[,which(names(Chain3) == "no_admixes")]

#We plot our trace plots (these include all samples, including those in the burn-in period)
posteriormax = max(c(Post1, Post2, Post3))
posteriormin =  min(c(Post1[-c(1:1000)], Post2[-c(1:1000)], Post3[-c(1:1000)]))
plot(Post1, type = "l", col ="steelblue4", xlab = "", ylab = "log-Posterior", main = "Chain 1", ylim = c(posteriormin, 130))
plot(Post2, type = "l", col ="steelblue4", xlab = "", ylab = "", main = "Chain 2", ylim = c(posteriormin, 130))
plot(Post3, type = "l", col ="steelblue4", xlab = "", ylab = "", main = "Chain 3", ylim = c(posteriormin, 130))

plot(Branches1, type = "l", col ="limegreen",  xlab = "", ylab = "Total Branch Length", ylim = c(2,12))
plot(Branches2, type = "l", col ="limegreen",  xlab = "", ylab = "", ylim = c(2,12))
plot(Branches3, type = "l", col ="limegreen",  xlab = "", ylab = "", ylim = c(2,12))

plot(admixes1, pch = 19, col ="orange", xlab = "Iteration", ylab = "Number of Admixture Events", cex = 0.48, ylim = c(0,  20))
plot(admixes2, pch = 19, col ="orange", xlab = "Iteration", ylab = "", cex = 0.48, ylim = c(0,  20))
plot(admixes3, pch = 19, col ="orange", xlab = "Iteration", ylab = "", cex = 0.48, ylim = c(0,  20))
