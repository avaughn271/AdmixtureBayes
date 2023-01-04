convertoutput = function(edgessss, verticesss, outputtt, outgroupppp) {
EDGES = read.table(edgessss)
VERTICES = read.table(verticesss)[,1:10]
numberss = VERTICES[[1]]
nodess = VERTICES[[2]]

A = matrix("", nrow = 0, ncol = 4)

for (edge in 1:nrow(EDGES)) {
  if (is.na(nodess[(which(EDGES[edge,2] == numberss))])) {
    A = rbind(A, c( EDGES[edge,2], EDGES[edge,1], EDGES[edge,3], EDGES[edge,4] ) )
  }
  else {  A = rbind(A, c( nodess[which(EDGES[edge,2] == numberss)] ,EDGES[edge,1], EDGES[edge,3], EDGES[edge,4])) }
}
a = A[which(A[,1] == commandArgs(trailingOnly=TRUE)[1]) , 2] # removes both of the outgroup lines

#a = A[which(A[,1] ==  "pop5") , 2]

if (length(a)  > 1) {print("admixtire involving the outgroup. Problem")}
outgroupdistance = as.numeric((A[A[,2] == a ,])[1,3]) * 2
A = A[ A[,2] != a   , ]

#It is here where we split the nodes in the appropriate way
problemnodes = c()
for (i in A[,1]) {
  if (sum(A[,1] == i) > 1){problemnodes = c(problemnodes, i)} 
}
if (length(unique(problemnodes)) * 2 != length(problemnodes) ) {print("problem")}
problemnodes = unique(problemnodes)
# for each of them, split if they are a leaf node or a drivergence node.
#perform a check that parents is only 2.
for (i in problemnodes) {
  if (i %in% A[,2]) {
    for (j in 1:length(A[,2])) {
      if (A[j,2] == i) {A[j,2] = paste0(i, "_split")}
    }
    A = rbind(A, c( paste0(i, "_split"), i , 0.000000000001, 1))
  }
  else {
    A[A == i] = paste0(i, "_split")
    A = rbind(A, c(i , paste0(i, "_split"), 0.000000000001, 1))
  }
}
write.table(A, outputtt, quote = F, row.names = F, col.names = F)
write.table(outgroupdistance, outgroupppp, quote = F, row.names = F, col.names = F)
}


for (jj in 1:20) {
convertoutput(
 paste0("TemporaryFiles/treeoutput",toString(jj) , ".edges" ),
 paste0("TemporaryFiles/treeoutput",toString(jj) , ".vertices" ),
 paste0("TemporaryFiles/TREEMIXoutput",toString(jj) , ".txt" ),
 paste0("TemporaryFiles/TREEMIXoutgroup",toString(jj) , ".txt" ))
}

for (jj in 1:20) {
convertoutput(
 paste0("TemporaryFiles/orientoutput",toString(jj) , ".edges" ),
 paste0("TemporaryFiles/orientoutput",toString(jj) , ".vertices" ),
 paste0("TemporaryFiles/ORIENToutput",toString(jj) , ".txt" ),
 paste0("TemporaryFiles/ORIENToutgroup",toString(jj) , ".txt" ))

}
