#tree2 will have the ones fixed

checkequality = function(tree1, tree2, tree2orig) {
  

  if (length(tree1) == 2 & length(tree2) == 2 & tree1[1] == tree2[1] &  tree1[2] == tree2[2]) {
    return(1)
  }
#  if (length(tree1) == 2 & length(tree2) == 2  ) { ### recently added
#    return(0)
#  }

  if (nrow(tree1) != nrow(tree2)) {return(0)}
  
  for (i in 1:nrow(tree2)) {
    if (!(tree2[i,1] %in% tree2[,2]) & sum(tree2[i,1] == tree2) == 1) {
      leaff2 = tree2[i,1]
      parent2 = tree2[i,2]
      if (sum(leaff2 == tree1) != 1) {return(0)}
      if (sum(leaff2 == tree1[,1]) != 1) {return(0)}
      j = which(leaff2 == tree1[,1])
      leaff1 = leaff2
      parent1 = tree1[j,2]
      if (substr(parent1, nchar(parent1) - 2, nchar(parent1)) == "_xx" &  parent1 != parent2) {return(0)}
      tree1[tree1 == parent1] = parent2
      return(checkequality(tree1[-j,], tree2[-i,], tree2orig))
    }
  }
  for (i in 1:nrow(tree1)) {
    for (j in 1:nrow(tree2)) {
      if (tree1[i,1] == tree2[j,1] & tree1[i,2] == tree2[j,2]) {
        return(checkequality(tree1[-i,], tree2[-j,], tree2orig))}
  }
  }

  for (i in 1:nrow(tree1)) {
      if (substr(tree1[i,1], nchar(tree1[i,1]) - 2, nchar(tree1[i,1])) == "_xx" &
      substr(tree1[i,2], nchar(tree1[i,2]) - 2, nchar(tree1[i,2])) == "_xx") {
if (  !(tree1[i,2] %in%   (tree2orig[which(tree1[i,1] == tree2orig[,1]), 2])  )){return(0)     }
  }
}

  print("PROBLEM")
return(-1000000)
}

EDGES1 = cbind(read.table("TemporaryFiles/MAPtree.txt")[,1],read.table("TemporaryFiles/MAPtree.txt")[,2])   
EDGES2 = cbind(read.table("TemporaryFiles/TrueTree.txt")[,1], read.table("TemporaryFiles/TrueTree.txt")[,2])   
for (i in 1:length(EDGES2)) {
  if (EDGES2[i] %in% EDGES2[,2]) EDGES2[i] = paste0(EDGES2[i], "_xx")
}
MAP = checkequality(EDGES1 , EDGES2, EDGES2)

write.table(MAP, file = "TemporaryFiles/TopEq.txt", row.names = F, quote = F)