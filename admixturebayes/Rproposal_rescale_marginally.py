from Rtree_operations import update_node
from copy import deepcopy
from numpy.random import normal, choice

class updater(object):
    
    def __init__(self, sigma):
        self.sigma=sigma

    def __call__(self):
        return normal(scale=self.sigma)
    
def get_random_key(tree):
    return choice(list(tree.keys()), size=1)[0]

def rescale(tree, sigma=0.01, pks={}):
    pks['rescale_marginally_adap_param']=sigma
    new_tree=deepcopy(tree)
    updat=updater(sigma)
    update_key=get_random_key(tree)
    new_tree=update_node(new_tree,update_key, updat)
    if new_tree is None:
        return tree,1,0 #rejecting by setting backward jump probability to 0.
    return new_tree ,1,1

class rescale_marginally_class(object):
    new_nodes=0
    proposal_name='rescale_marginally'
    adaption=True
    input='tree'
    require_admixture=0
    reverse_require_admixture=0
    admixture_change=1
    reverse='rescale_marginally'
    
    def __call__(self,*args, **kwargs):
        return rescale(*args, **kwargs)
    
