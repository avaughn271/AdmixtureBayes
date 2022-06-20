args <- commandArgs(TRUE)
filename=args[1]
proportion=as.numeric(args[2])
outname=args[3]
verbose_level=args[4]
summaries=args[5:length(args)]

suppressMessages(library('coda', quietly = T))
suppressMessages(library('rwty', quietly = T))

summaries_without_trees=c()
tree_summaries=c()
subset_summary=FALSE

for(summary in summaries){
	if(grepl('Ntree', summary)){
		tree_summaries=c(tree_summaries, summary)
	}
  else if(grepl('descendant_sets', summary)){
    subset_summary=TRUE
  }
	else{
		summaries_without_trees=c(summaries_without_trees, summary)
	}
}

get_focal_element=function(listi){
  return(listi[sample(length(listi),1)])
}

convert_to_sets=function(listi){
  return(strsplit(listi, split = '-',fixed = T))
}

distance=function(a,v){
  return(length(c(setdiff(a,v),setdiff(v,a))))
}


df=read.csv(filename, header=T)

if( 'layer' %in% colnames(df)){
dfa=subset(df, layer==0)} else{
  dfa=df
}

dfa=dfa[(floor(proportion*nrow(dfa))):nrow(dfa),]
dfb=as.data.frame(dfa[,summaries_without_trees])
colnames(dfb) <- summaries_without_trees
df2=apply(dfb,c(1,2),as.numeric)

if(subset_summary){
  l=as.character(dfa$descendant_sets)
  l2=convert_to_sets(l)
  focal_set=get_focal_element(l2)[[1]]
  set_dists=sapply(X = l2, FUN = distance, v=focal_set)
  old_cols=colnames(df2)
  df2=cbind(df2,as.numeric(set_dists))
  colnames(df2) <- c(old_cols, 'descendant_sets')
  summaries_without_trees=c(summaries_without_trees, 'descendant_sets')
}

all_nums=function(df){
  mcmcobj=mcmc(df)
  res_df=data.frame(summary=colnames(df))
  return(effectiveSize(mcmcobj))
}

tree_nums=function(df){
	ids=floor(seq(1,nrow(df), length.out = min(nrow(df),1000)))
	df=df[ids,,drop=F]
	res=c()
	for(tree_summary in tree_summaries){
		write(paste0('\t','tree gen.', ids,' = [&U] ',as.character(df[,tree_summary]),';'),'trees_tmp.txt')
		invisible(capture.output(chain <- load.trees('trees_tmp.txt', type='newick')))
		invisible(capture.output(a <- topological.pseudo.ess(chain,n=1)))
		if(a>0.5){#if all trees are identical, it will cause a ess of 0 - and not 500 as it should. Therefore, I put in the recognizable value, 777.77777
			res=c(res, a[1,1])
		}
		else{
			res=c(res,777.77777)
		}
	}
	return(res)
}

esss=all_nums(df2)

dfc=as.data.frame(dfa[,tree_summaries])
colnames(dfc) <- tree_summaries
tree_esss=tree_nums(dfc)


to_print=cbind(c(summaries_without_trees,tree_summaries),c(esss,tree_esss))

write.table(to_print, file= outname, quote=FALSE)
