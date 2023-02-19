
par(mfrow=c(3,1))
B = read.table("AdBayesCov.txt",header = F)

boxplot(B, ylab = "Covariance Distance", xaxt = "n",
        cex.axis = 1.3 ,  cex.lab = 1.5,
       main = "Method Comparison for Graph G1",
       ylim = c(0, covarmax), col = "steelblue3", range = 0)

B = read.table("AdBayesSet.txt",header = F)

boxplot(B, ylab = "Set Distance", xaxt = "n",
        cex.axis = 1.3 , cex.lab = 1.5,range=0, col = "steelblue3")

B = read.table("AdBayesTop.txt", header = F)

boxplot(B, ylab = "Topology Equality",
        cex.lab = 1.5, ylim = c(0,1),  cex.axis = 1.3,
        col = "steelblue3" , range=0)
