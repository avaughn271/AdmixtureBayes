from Rtree_operations import update_all_admixtures, get_number_of_admixes, get_specific_branch_lengths, get_leaf_keys, update_specific_branch_lengths
from copy import deepcopy
from numpy.random import normal
from math import sqrt, exp
from Rtree_to_coefficient_matrix import make_coefficient_matrix
from scipy.stats import norm
from numpy import identity, array

class updater(object):
    
    def __init__(self, sigma):
        self.sigma=sigma

    def __call__(self):
        return normal(scale=self.sigma)

def rescale_admix_correction(tree, sigma=0.01, pks={}, make_correction=True, return_without_correction=False):
    k=get_number_of_admixes(tree)
    pks['rescale_admix_correction_adap_param']=sigma
    new_tree=deepcopy(tree)
    if k>0:
        updat=updater(sigma/sqrt(k))
        new_tree=update_all_admixtures(new_tree, updat)
        if new_tree is None:
            return tree,1,0 #rejecting by setting backward jump probability to 0.
    else:
        return new_tree,1,0.234 #not to have to deal with the admix=0 case, I return 0.234 such that the adaptive parameter is not affected by these special cases.
    if make_correction:
        untouched_tree=deepcopy(new_tree)
        new_tree, qforward, qbackward = getcorrection(tree, new_tree,sigma/20)
    else:
        qforward=qbackward=1.0
    if new_tree is None:
        return tree, 1,0
    if return_without_correction:
        return new_tree ,qforward, qbackward, untouched_tree
    return new_tree ,qforward, qbackward#, untouched_tree

def mm(U,L,initial_value, tol=0.00000001):
    x_old=deepcopy(initial_value)
    #removedprin 'x_old', x_old
    x_new=U/(L.dot(x_old))*x_old
    #removedprin 'multipliers', U/(L.dot(x_old))
    #removedprin 'x_new', x_new
    while sum((x_new-x_old)**2)>tol:
        x_old=x_new
        x_new=U/(L.dot(x_old))*x_old
    return x_new

def getcorrection(old_tree, new_tree,sigma):
    
    node_keys=sorted(get_leaf_keys(old_tree))
    
    B,_,bi1= make_coefficient_matrix(old_tree, node_keys=node_keys)
    A,_,bi2= make_coefficient_matrix(new_tree, node_keys=node_keys)
    
    branches=reverse_dic_to_list(bi1)
    
    x_A=array(get_specific_branch_lengths(old_tree, branches))
    x_B=array(get_specific_branch_lengths(new_tree, branches))
    x_old=deepcopy(x_A)
   # print x_A
    #removedprin x_B
    #removedprin x_old
    
    upper=x_A.dot(B.T.dot(A)+identity(len(branches)))
    
    
    lower_first=A.T.dot(A)+identity(len(branches))
    
    mu_new=mm(U=upper, L=lower_first, initial_value=x_B)
    
    x_new=mu_new+norm.rvs(scale=sigma, size= len(mu_new))
    
    q_forward=sum(norm.logpdf(mu_new-x_new, scale=sigma))
    
    upper_reverse=x_new.dot((A.T.dot(B)+identity(len(branches))))
    lower_first_reverse=B.T.dot(B)+identity(len(branches))
    
    mu_reverse=mm(U=upper_reverse, L=lower_first_reverse, initial_value=array(x_new))
        
    q_backward=sum(norm.logpdf(mu_reverse-x_A, scale=sigma))
        
    new_tree=update_specific_branch_lengths(new_tree,branches, x_new)


    return new_tree, 1.0, exp(q_backward-q_forward)
    
    
def reverse_dic_to_list(dic):
    l_determined=[None]*len(dic)
    for e,i in list(dic.items()):
        l_determined[i]=e
    return l_determined

class rescale_admix_correction_class(object):
    new_nodes=0
    proposal_name='rescale_admix_correction'
    adaption=True
    input='tree'
    require_admixture=1
    reverse_require_admixture=1
    admixture_change=0
    reverse='rescale_admix_correction'
    
    def __call__(self, *args, **kwargs):
        return rescale_admix_correction(*args, **kwargs)
