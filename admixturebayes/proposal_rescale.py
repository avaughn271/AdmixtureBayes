from copy import deepcopy
from numpy.random import normal

sigma=0.01
def rescale(tree,pks={}):
    
    cop=deepcopy(tree)#only changes the memory in this line
    
    for branch in cop:
        #rescale branchlengths
        branch[2::2]=branch[2::2]+normal(0,sigma,size=len(branch[2::2]))
        
        #rescale admixtureproportions
        for admixevent in branch[3::2]:
            if admixevent[2] is not None:
                admixevent[2]=admixevent[2]+normal(0,sigma,size=1)[0]               
    return cop, 1,1,1