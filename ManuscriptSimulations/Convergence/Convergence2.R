library("coda") # We load the coda library

#This is our fraction of samples to discard as burn-in. This can be changed to
#decide what a good burn-in is
BurnInFraction = 0.35

#We assume that the output of the 3 independently run chains are
#chain1.txt, chain2.txt, and chain3.txt
Chain1 = read.csv("chain1.csv")
Chain2 = read.csv("chain2.csv")
Chain3 = read.csv("chain3.csv")

#We extract the total branch length
Branches1 = Chain1[,which(names(Chain1) == "total_branch_length")]
Branches2 = Chain2[,which(names(Chain2) == "total_branch_length")]
Branches3 = Chain3[,which(names(Chain3) == "total_branch_length")]

#And our total branch length after a burn-in
Branches1Burn = Branches1[-(1:round(BurnInFraction*length(Branches1)))]
Branches2Burn = Branches2[-(1:round(BurnInFraction*length(Branches2)))]
Branches3Burn = Branches3[-(1:round(BurnInFraction*length(Branches3)))]

#We do the same thing for the posterior
Post1 = Chain1[,which(names(Chain1) == "posterior")]
Post2 = Chain2[,which(names(Chain2) == "posterior")]
Post3 = Chain3[,which(names(Chain3) == "posterior")]
Post1Burn = Post1[-(1:round(BurnInFraction*length(Post1)))]
Post2Burn = Post2[-(1:round(BurnInFraction*length(Post2)))]
Post3Burn = Post3[-(1:round(BurnInFraction*length(Post3)))]

#And for the number of admixture events
admixes1 = Chain1[,which(names(Chain1) == "no_admixes")]
admixes2 = Chain2[,which(names(Chain2) == "no_admixes")]
admixes3 = Chain3[,which(names(Chain3) == "no_admixes")]
admixes1Burn = admixes1[-(1:round(BurnInFraction*length(admixes1)))]
admixes2Burn = admixes2[-(1:round(BurnInFraction*length(admixes2)))]
admixes3Burn = admixes3[-(1:round(BurnInFraction*length(admixes3)))]

#We reformat the plotting space for the Gelman-Rubin Plots
par(mfrow=c(1,1))

#This computes and plots the median of the scale reduction factor for total branch length
#after the burn-in period for a cumulative increasing number of samples
y = c()
x = c()
for ( i in 1:50) {
  seq = 1:(round(length(Branches1Burn)/50*i))
  x = c(x, (round(length(Branches1Burn)/50*i)))
 a = gelman.diag(mcmc.list(list(mcmc(Branches1Burn[seq]), mcmc(Branches2Burn[seq]), mcmc(Branches3Burn[seq]))), autoburnin = F)
y = c(y, (a[[1]][1]) )
 }
plot(x, y, type = "l", ylab = "Gelman-Rubin Scale Reduction Factor", main = "Gelman-Rubin Convergence Diagnostics",
     xlab = paste0("Cumulative Total Samples (After Burn-in of ",toString(round(BurnInFraction*length(admixes1))) , " Samples)"),
     col ="steelblue4", lwd = 2, ylim = c(0.98, 1.5))

#We do the same thing for the number of admixture events
y = c()
x = c()
for ( i in 1:50) {
  seq = 1:(round(length(admixes1Burn)/50*i))
  x = c(x, (round(length(admixes1Burn)/50*i)))
  a = gelman.diag(mcmc.list(list(mcmc(admixes1Burn[seq]), mcmc(admixes2Burn[seq]), mcmc(admixes3Burn[seq]))), autoburnin = F)
  y = c(y, (a[[1]][1]) )
}

lines(x, y, col ="orange", lwd = 2)

#And the posterior
y = c()
x = c()
for ( i in 1:50) {
  seq = 1:(round(length(Post1)/50*i))
  x = c(x, (round(length(Post1)/50*i)))
  a = gelman.diag(mcmc.list(list(mcmc(Post1Burn[seq]), mcmc(Post2Burn[seq]), mcmc(Post3Burn[seq]))), autoburnin = F)
  y = c(y, (a[[1]][1]))
}

lines(x, y, col ="limegreen", lwd = 2)

#We plot a horizontal line at 1 for reference as well as a legend
lines(c(-100, 1e10), c(1,1), col = "red", lty = 3)

legend("topright", legend=c("Posterior", "Total Branch Length", "Number of Admixtures"),
       col=c("steelblue4", "limegreen", "orange"), lwd = 2)
