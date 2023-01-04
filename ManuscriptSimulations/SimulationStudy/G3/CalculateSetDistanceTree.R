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
print(caltopsets(read.table("TemporaryFiles/TrueTree.txt")))
TreeMIXResult = (length(setdiff(caltopsets(read.table("TemporaryFiles/TREEMIXoutput.txt")), 
                                caltopsets(read.table("TemporaryFiles/TrueTree.txt")))))+ (length(setdiff(caltopsets(read.table("TemporaryFiles/TrueTree.txt")), 
                  caltopsets(read.table("TemporaryFiles/TREEMIXoutput.txt")))))

OrientResult = (length(setdiff(caltopsets(read.table("TemporaryFiles/ORIENToutput.txt")), 
                                caltopsets(read.table("TemporaryFiles/TrueTree.txt")))))+ (length(setdiff(caltopsets(read.table("TemporaryFiles/TrueTree.txt")), 
                 caltopsets(read.table("TemporaryFiles/ORIENToutput.txt")))))
 
write.table(data.frame( OrientResult ,TreeMIXResult), file = "TemporaryFiles/SetDistance.txt",
            row.names = F, quote = F)