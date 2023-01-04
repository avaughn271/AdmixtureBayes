pdf(file = "Plot.pdf",     width = 9*0.65,  height = 10*0.65)

par(mfrow=c(3,1))
par(mar = c(1, 4, 1.3, 0.3))  
covarmax = 1
A = read.table("TreemixCov.txt",header = F)

B = data.frame(matrix(c(0.249948555,	0.248737205,
                        0.264279468,	0.268954394,
                        0.222525683,	0.219334758,
                        0.277909792,	0.279420573,
                        0.228400199,	0.220768156,
                        0.24465227,	0.25013655,
                        0.256573229,	0.254464029,
                        0.296461224,	0.28947583,
                        0.208463023,	0.211412666,
                        0.262552207,	0.262942794,
                        0.27156618,	0.273033464,
                        0.220934065,	0.229304025,
                        0.260434998,	0.262876629,
                        0.267127613,	0.270224305,
                        0.252486084,	0.254642533), ncol  = 2, byrow = 2))


df = data.frame(rep(A[[2]]) , rep(A[[1]]) ,  rep(B[[1]],4) , rep(B[[2]],4) )
boxplot(df, ylab = "Covariance Distance", xaxt = "n", 
       main = "Method Comparison for Graph G3", ylim = c(0, covarmax), col = c("red", "green"  ,"steelblue3" , "lightskyblue"), range = 0)
par(mar = c(2, 4, 1, 0.3))

A = read.table("TreemixSet.txt",header = F)


B = data.frame(matrix(c(0	,0,
                        0	,0,
                        0	,0,
                        0	,0,
                        0	,0,
                        0	,0,
                        0,	0,
                        0,	0,
                        0,	0.06,
                        0,	0,
                        0,	0,
                        0,	0.03,
                        0,	0,
                        0,	0,
                        0,	0), ncol  = 2, byrow = 2))


df = data.frame(  A[[2]] , A[[1]] , rep(B[[1]],4) , rep(B[[2]],4) )

boxplot(df, ylab = "Set Distance", xaxt = "n", ylim = c(0,6), range=0, col = c("red", "green"  ,"steelblue3" , "lightskyblue"))
par(mar = c(3, 4, 0, 0.3))

A = read.table("TreemixTop.txt",header = F)

B = data.frame(matrix(c(1,	1,
                        1,	0.99,
                        1,	1,
                        1,	1,
                        1,	0.99,
                        1,	0.99,
                        1,	1,
                        1,	0.99,
                        1,	0.98,
                        1,	1,
                        1,	1,
                        1,	0.99,
                        1,	1,
                        1,	0.99,
                        1,	0.99), ncol  = 2, byrow = 2))


df = data.frame(   mean(A[[2]] ),mean(A[[1]]) ,  mean( rep(B[[1]],4)) , rep(B[[2]],4) )
names(df) = c("TreeMix (3)", "OrientAGraph (20)", "AdmixtureBayes MAP", "AdmixtureBayes Mean")

boxplot(df, ylab = "Topology Equality", ylim = c(0,1), col = c("red", "green"  ,"steelblue3" , "lightskyblue"), range=0)

###don't plot treemix, as it does not work. orientagraph does very poorly, and only includes 6 examples
dev.off()
