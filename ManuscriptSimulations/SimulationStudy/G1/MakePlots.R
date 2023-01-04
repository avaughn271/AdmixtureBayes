

pdf(file = "Plot.pdf",     width = 9*0.65,  height = 10*0.65)

par(mfrow=c(3,1))
par(mar = c(1, 4, 1.3, 0.3))  
covarmax = 0.5
A = read.table("TreemixCov.txt",header = F)
A = rbind(A, c(min(A[[1]]), min(A[[2]]) ))
B = read.table("AdBayesCov.txt",header = F)

df = data.frame(A[[2]], A[[1]], B[[1]] , B[[2]])
boxplot(df, ylab = "Covariance Distance", xaxt = "n", 
       main = "Method Comparison for Graph G1", ylim = c(0, covarmax), col = c("red", "green"  ,"steelblue3" , "lightskyblue"), range = 0)
par(mar = c(2, 4, 1, 0.3))

A = read.table("TreemixSet.txt",header = F)
A = rbind(A, c(min(A[[1]]), min(A[[2]]) ))

B = read.table("AdBayesSet.txt",header = F)

df = data.frame(A[[2]], A[[1]], B[[1]] , B[[2]])

boxplot(df, ylab = "Set Distance", xaxt = "n" , range=0, col = c("red", "green"  ,"steelblue3" , "lightskyblue"))
par(mar = c(3, 4, 0, 0.3))

A = read.table("TreemixTop.txt",header = F)
A = rbind(A, c(min(A[[1]]), min(A[[2]]) ))
B = read.table("AdBayesTop.txt",header = F)

df = data.frame(A[[2]], A[[1]], B[[1]] , B[[2]])
names(df) = c("TreeMix (20)", "OrientAGraph (20)", "AdmixtureBayes MAP", "AdmixtureBayes Mean")

boxplot(df, ylab = "Topology Equality", ylim = c(0,1), col = c("red", "green"  ,"steelblue3" , "lightskyblue"), range=0)






dev.off()
