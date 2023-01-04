TotalList = list()
Type1 = c()
Type2 = c()
Type3 = c()
#              4        0.3
par(mar = c(1, 4, 2, 0.3))  

layout(matrix(c(1,2), byrow = TRUE), heights=c(7,8))

for (i in 1:100) {
  datatable = read.csv(paste0("Small/" ,toString(i), "rankings.txt") , header = F)
  
  temp = c(Type1, (datatable[[2]])[which(datatable[[1]] == "c.0.w-c.0")])
  if (length(temp) == length(Type1)) {Type1  = c(Type1 , 0.0)}
  else {Type1 = temp}
  temp = c(Type2, (datatable[[2]])[which(datatable[[1]] == "c.w.0-c.0")])
  
  if (length(temp) == length(Type2)) {Type2  = c(Type2 , 0.0)}
  else {Type2 = temp}
  temp = c(Type3, (datatable[[2]])[which(datatable[[1]] == "w.c.1-c.0")])
  
  if (length(temp) == length(Type3)) {Type3  = c(Type3 , 0.0)}
  else {Type3 = temp}
}

#Type 1 is the true value

boxplot(data.frame(Type1, Type2, Type3 )  ,  #xlab = "Topology", 
        ylab = "Posterior Density" , ylim = c(0,1), range =  0,  
        col = "lightskyblue", main = "4 Haplotypes in Each Population", xaxt = "n")
#axis(1, at = 1:3, labels=c("((pop1,pop2),pop3)", "((pop2,pop3),pop1)", "((pop1,pop3),pop2)"))

lines(c(-100,100), c(1/3, 1/3), col  ="red", lwd = 2,lty = 2)

print(sum(Type1) + sum(Type2) + sum(Type3))



TotalList = list()
Type1 = c()
Type2 = c()
Type3 = c()

for (i in 1:100) {
  datatable = read.csv(paste0("Large/" ,toString(i), "rankings.txt") , header = F)
  
  temp = c(Type1, (datatable[[2]])[which(datatable[[1]] == "c.0.w-c.0")])
  if (length(temp) == length(Type1)) {Type1  = c(Type1 , 0.0)}
  else {Type1 = temp}
  temp = c(Type2, (datatable[[2]])[which(datatable[[1]] == "c.w.0-c.0")])
  
  if (length(temp) == length(Type2)) {Type2  = c(Type2 , 0.0)}
  else {Type2 = temp}
  temp = c(Type3, (datatable[[2]])[which(datatable[[1]] == "w.c.1-c.0")])
  
  if (length(temp) == length(Type3)) {Type3  = c(Type3 , 0.0)}
  else {Type3 = temp}
}

#Type 1 is the true value
par(mar = c(4.2, 4, 2, 0.3))  

boxplot(data.frame(Type1, Type2, Type3 )  , xlab = "Topology", 
        ylab = "Posterior Density" ,  main = "40 Haplotypes in Each Population", ylim = c(0,1), range =  0,  
        col = "lightskyblue", xaxt = "n")
axis(1, at = 1:3, labels=c("((pop1,pop2),pop3)", "((pop2,pop3),pop1)", "((pop1,pop3),pop2)"))

lines(c(-100,100), c(1/3, 1/3), col  ="red", lwd = 2,lty = 2)

print(sum(Type1) + sum(Type2) + sum(Type3))