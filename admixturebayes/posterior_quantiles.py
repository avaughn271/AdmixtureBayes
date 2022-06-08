
import pandas as pd

def get_numeric(string_tree):
    branch_lengths_string, admixture_proportion_string= string_tree.split(';')[1:]
    branch_lengths=list(map(float,branch_lengths_string.split('-')))
    if admixture_proportion_string:
        admixture_proportions=list(map(float, admixture_proportion_string.split('-')))
    else:
        admixture_proportions=[]
    return branch_lengths, admixture_proportions


def branch_and_proportion_quantiles(list_of_string_trees):
    ''' extracts the branch lengths and admixture proportions and returns them as tuples of four. the list should not be empty'''
    branches=[]
    admixtures=[]
    for string_tree in list_of_string_trees:
        b,a=get_numeric(string_tree)
        branches.append(b)
        admixtures.append(a)
    bdat=pd.DataFrame.from_records(branches)
    adat=pd.DataFrame.from_records(admixtures)
    bresults=[]
    aresults=[]
    for n,(lower,mean, upper) in enumerate(zip(bdat.quantile(0.025),bdat.mean(),bdat.quantile(0.975))):
        bresults.append(('c'+str(n+1), lower,mean,upper))
    if len(admixtures[0])>0:
        for n,(lower,mean, upper) in enumerate(zip(adat.quantile(0.025),adat.mean(),adat.quantile(0.975))):
            aresults.append(('ax'+str(n+1), lower,mean,upper))
    return bresults, aresults
