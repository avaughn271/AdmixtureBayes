convertoutput = function(edgessss, verticesss, outputtt, outgroupppp) {
EDGES = read.table(edgessss)
VERTICES = read.table(verticesss)[,1:10]
outgroupnode = VERTICES[[1]][which(VERTICES[[2]] == "pop5")]
rootnode = VERTICES[[1]][which(VERTICES[[3]] == "ROOT")]

if (sum(EDGES[[2]] == outgroupnode) > 1) {return("BAD")}

  if (EDGES[which(EDGES[[2]] == outgroupnode),1] == rootnode) {
    return("GOOD")
}
return("BAD")
}
args = 
 j  = as.numeric(commandArgs(trailingOnly=TRUE)[2]) 

 TREEEE = c()
 ORIENTTTT = c()
for (ii in 1:j) {
  TREEEE = c(TREEEE,
convertoutput(
 paste0("TemporaryFiles/treeoutput",toString(ii) , ".edges" ),
 paste0("TemporaryFiles/treeoutput",toString(ii) , ".vertices" ),
 paste0("TemporaryFiles/TREEMIXoutput",toString(ii) , ".txt" ),
 paste0("TemporaryFiles/TREEMIXoutgroup",toString(ii) , ".txt" )))
  ORIENTTTT = c(ORIENTTTT, convertoutput(
 paste0("TemporaryFiles/orientoutput",toString(ii) , ".edges" ),
 paste0("TemporaryFiles/orientoutput",toString(ii) , ".vertices" ),
 paste0("TemporaryFiles/ORIENToutput",toString(ii) , ".txt" ),
 paste0("TemporaryFiles/ORIENToutgroup",toString(ii) , ".txt" )))

}
 TREEEE
 ORIENTTTT
 
 sum(TREEEE == "GOOD")
 sum(ORIENTTTT == "GOOD")