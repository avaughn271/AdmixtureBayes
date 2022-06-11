from Rtree_operations import update_all_admixtures, get_number_of_admixes
from copy import deepcopy
from numpy.random import normal
from math import sqrt

class updater(object):
    
    def __init__(self, sigma):
        self.sigma=sigma

    def __call__(self):
        return normal(scale=self.sigma)

def rescale_admixtures(tree, sigma=0.01, pks={}):
    k=get_number_of_admixes(tree)
    pks['rescale_admixtures_adap_param']=sigma
    new_tree=deepcopy(tree)
    if k>0:
        updat=updater(sigma/sqrt(k))
        new_tree=update_all_admixtures(new_tree, updat)
        if new_tree is None:
            return tree,1,0 #rejecting by setting backward jump probability to 0.
    else:
        return new_tree,1,0.234 #not to have to deal with the admix=0 case, I return 0.234 such that the adaptive parameter is not affected by these special cases.
    return new_tree ,1,1

class rescale_admixtures_class(object):
    new_nodes=0
    proposal_name='rescale_admixtures'
    adaption=True
    input='tree'
    require_admixture=1
    reverse_require_admixture=1
    admixture_change=0
    reverse='rescale_admixtures'
    
    def __call__(self, *args, **kwargs):
        return rescale_admixtures(*args, **kwargs)

def rescale(add, sigma=0.01, pks={}):
    pks['rescale_add_adap_param']=sigma
    new_add=add+normal()*sigma
    if new_add<0:
        return add,1,0 #rejecting by setting backward jump probability to 0.
    return new_add,1,1

class rescale_add_class(object):
    new_nodes=0
    proposal_name='rescale_add'
    adaption=True
    input='add'
    require_admixture=0
    reverse_require_admixture=0
    reverse='rescale_add'
    admixture_change=0
    
    def __call__(self,*args, **kwargs):
        return rescale(*args, **kwargs)

    