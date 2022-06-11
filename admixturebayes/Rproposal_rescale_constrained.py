from Rtree_operations import update_specific_branch_lengths
from copy import deepcopy
from numpy.random import normal
from Rtree_to_coefficient_matrix import get_orthogonal_branch_space

    
def reverse_dic_to_list(dic):
    l_determined=[None]*len(dic)
    for e,i in list(dic.items()):
        l_determined[i]=e
    return l_determined

def get_added_branch_pieces(org, param):
    return org.dot(normal(scale=param, size=org.shape[1]))

def rescale(x, sigma=0.01, pks={}, update_add=True):
    tree, add=x
    pks['rescale_constrained_adap_param']=sigma
    new_tree=deepcopy(tree)
    orgs, bi,_= get_orthogonal_branch_space(new_tree, add_one_column=update_add)
    #removedprin norm(orgs, axis=0)
    branches=reverse_dic_to_list(bi)
    branch_pieces= get_added_branch_pieces(orgs, sigma)
    if update_add:
        b=branch_pieces[:-1]
        #removedprin 'ADD!'
    else:
        #removedprin 'NO ADD!'
        b=branch_pieces
    new_tree= update_specific_branch_lengths(new_tree, branches, b, add=True)
    if new_tree is None:
        return x,1,0 #rejecting by setting backward jump probability to 0.
    if update_add:
        new_add=add+branch_pieces[-1]
    else:
        new_add=0
    if new_add<0:
        return x,1,0 #rejecting by setting backward jump probability to 0.
    return  (new_tree, new_add),1,1

class rescale_constrained_class(object):
    new_nodes=0
    proposal_name='rescale_constrained'
    adaption=True
    input='both'
    require_admixture=0
    reverse_require_admixture=0
    admixture_change=0
    reverse='rescale_constrained'
    
    
    def __init__(self, **kwargs):
        self.kwargs=kwargs
    
    def __call__(self,*args, **kwargs):
        kwargs.update(self.kwargs)
        #removedprin kwargs
        return rescale(*args, **kwargs)

