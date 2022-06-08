from Rtree_operations import update_all_branches, get_number_of_leaves, get_number_of_admixes
from copy import deepcopy
from numpy.random import normal
from math import sqrt

class updater(object):
    
    def __init__(self, sigma):
        self.sigma=sigma

    def __call__(self):
        return normal(scale=self.sigma)

def rescale(tree, sigma=0.01, pks={}):
    n=get_number_of_leaves(tree)
    k=get_number_of_admixes(tree)
    pks['rescale_adap_param']=sigma
    new_tree=deepcopy(tree)
    updat=updater(sigma/sqrt(2*n-2+4*k))
    new_tree=update_all_branches(new_tree, updat)
    if new_tree is None:
        return tree,1,0 #rejecting by setting backward jump probability to 0.
    return new_tree ,1,1

class rescale_class(object):
    new_nodes=0
    proposal_name='rescale'
    adaption=True
    input='tree'
    require_admixture=0
    reverse_require_admixture=0
    reverse='rescale'
    admixture_change=0
    
    def __call__(self,*args, **kwargs):
        return rescale(*args, **kwargs)
