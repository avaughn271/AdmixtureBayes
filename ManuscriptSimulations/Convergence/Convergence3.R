Chain1 = read.csv("chain1.csv")

Branches1 = Chain1[,which(names(Chain1) == "total_branch_length")]

Post1 = Chain1[,which(names(Chain1) == "posterior")]

par(mfcol=c(2,1))
acf(Post1[-(1:round(0.35*length(Post1)))],lag.max = 5000 , ci = -1, col ="steelblue4", ylim = c(-0.2,1), main = "Autocorrelation of Posterior (Chain 1)", ylab = "Autocorrelation")
acf(Branches1[-(1:round(0.35*length(Branches1)))],lag.max = 5000 , ci = -1, col ="limegreen", ylim = c(-0.2,1), main = "Autocorrelation of Total Branch Length (Chain 1)", ylab = "Autocorrelation")