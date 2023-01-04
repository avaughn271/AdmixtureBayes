pdf(file = "Plot.pdf",     width = 9*0.65,  height = 10*0.65)

par(mfrow=c(3,1))
par(mar = c(1, 4, 1.3, 0.3))  
covarmax = 1
A = read.table("TreemixCov.txt",header = F)
B = data.frame(matrix(c(0.42331234,	0.4125072,
0.512455486,	0.518382022,
0.460252227,	0.460711294,
0.508685644,	0.512046658,
0.47414573,	0.474331157,
0.488320154	,0.494007502,
0.405355171	,0.400461537,
0.424235415	,0.425095789,
0.500478547	,0.497633966,
0.455246635,	0.45519789,
0.450010888,	0.438622798,
0.47791136,	0.486261583,
0.347333764,	0.352724506,
0.41839384,	0.419397172,
0.489177204	,0.485933159), ncol  = 2, byrow = 2))

df = data.frame( rep(1000,15) ,  rep(A[[2]],15) ,rep(B[[1]],6) , rep(B[[2]],6) )
boxplot(df, ylab = "Covariance Distance", xaxt = "n", 
       main = "Method Comparison for Graph G4", ylim = c(0, covarmax), col = c("red", "green"  ,"steelblue3" , "lightskyblue"), range = 0)
par(mar = c(2, 4, 1, 0.3))

A = read.table("TreemixSet.txt",header = F)

B = data.frame(matrix(c(0,	0.62,
0,	0.08,
0	,0.35,
2,	1.75,
1,	1.21,
1,	0.65,
1,	0.57,
0,	1.01,
1,	0.78,
0,	0.38,
0,	0.6,
0,	0.63,
1,	0.97,
0,	0.54,
1,	1.5), ncol  = 2, byrow = 2))





df = data.frame(  rep(1000,15) ,rep(A[[1]],15) , rep(B[[1]],6) , rep(B[[2]],6) )

boxplot(df, ylab = "Set Distance", xaxt = "n", ylim = c(0,10), range=0, col = c("red", "green"  ,"steelblue3" , "lightskyblue"))
par(mar = c(3, 4, 0, 0.3))

A = read.table("TreemixTop.txt",header = F)

B = data.frame(matrix(c(1,	0.48,
1,	0.96,
1,	0.71,
0,	0.06,
0,	0,
0,	0.46,
0,	0.55,
1,	0.31,
0,	0.3,
1,	0.7,
1,	0.65,
1,	0.56,
0,	0.25,
1,	0.6,
0,	0), ncol  = 2, byrow = 2))


df = data.frame( rep(1000,15) , rep(A[[2]],15), rep(mean(B[[1]]),6) , rep(B[[2]],6) )
names(df) = c("TreeMix (0)", "OrientAGraph (6)", "AdmixtureBayes MAP", "AdmixtureBayes Mean")

boxplot(df, ylab = "Topology Equality", ylim = c(0,1), col = c("red", "green"  ,"steelblue3" , "lightskyblue"), range=0)

###don't plot treemix, as it does not work. orientagraph does very poorly, and only includes 6 examples
dev.off()
