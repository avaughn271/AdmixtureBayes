gettopsetnode = function(Table, currentnode) {
  if (!(currentnode %in% Table[[2]])) {
    return(currentnode)
  }
  else if (sum(currentnode ==  Table[[2]]) == 1){
    return(gettopsetnode(Table, (Table[[1]])[which(Table[[2]] == currentnode)]))
  }
  else {
    return(union(
      gettopsetnode(Table, (Table[[1]])[which(Table[[2]] == currentnode)[1]]),
      gettopsetnode(Table, (Table[[1]])[which(Table[[2]] == currentnode)[2]])
    ))
  }
}

caltopsets = function(Table){
  LISTT = list()
  indexx = 1
  for (i in unique(c(Table[[1]], Table[[2]]))) {
    LISTT[[indexx]] = gettopsetnode(Table, i)
    indexx = indexx + 1
  }
  for (i in 1:length(LISTT)) {LISTT[[i]] = sort(LISTT[[i]])}
  
  for (i in length(LISTT):2) {
    for (j in (i-1):1) {
      if (setequal(LISTT[[i]],  LISTT[[j]]  )) {
        LISTT[[i]] = NULL
        break
      }
    }
  }
  for (i in length(LISTT):1) {
    if (length(LISTT[[i]]) > 1) {
      temp = (LISTT[[i]])[1]
      for (jj in 2:length(LISTT[[i]])) {
        temp = paste0(temp, (LISTT[[i]])[jj])
      }
      LISTT[[i]] = temp
    }
  }
  return(unlist(LISTT))
}

args = as.numeric(commandArgs(trailingOnly=TRUE)[1]) 
for (j in 1:args)  {
  jj = toString(j)
  maptree = read.table(paste0("TemporaryFiles/MAPtree" , jj , ".txt"))
  truetree = read.table(paste0("TemporaryFiles/Truetree.txt"))
  
MAP = (length(setdiff(caltopsets(maptree), caltopsets(truetree)))) + (length(setdiff(caltopsets(truetree), caltopsets(maptree))))

write.table(MAP, file = paste0("TemporaryFiles/SetDistance" , jj, ".txt"), row.names = F, quote = F)

}