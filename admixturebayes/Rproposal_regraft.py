from copy import deepcopy
from numpy.random import choice
from scipy.stats import expon, uniform
from Rtree_operations import (node_is_non_admixture,
remove_parent_attachment, graft, get_real_parents, halfbrother_is_uncle,
get_branch_length, get_all_branch_descendants_and_rest)
from random import getrandbits

def get_possible_regrafters(tree):
    res=[]
    for key in tree:
        parents=get_real_parents(tree[key])
        #removedprin parents
        for branch, parent in enumerate(parents):
            #removedprin key,(not is_root(parent)),node_is_non_admixture(tree[parent]),(not has_child_admixture(tree, parent))
            #if (not is_root(parent)) and node_is_non_admixture(tree[parent]) and (not has_child_admixture(tree, parent)):
            if parent=='r' or (node_is_non_admixture(tree[parent]) and not halfbrother_is_uncle(tree, key, parent)):
                res.append((key,branch))
    return res

class regraft_class(object):
    
    new_nodes=1
    proposal_name='regraft'
    input='tree'
    require_admixture=0
    reverse_require_admixture=0
    adaption=False
    reverse='regraft'
    admixture_change=0
    
    def __call__(self,*args, **kwargs):
        return make_regraft(*args, **kwargs)

def make_regraft(tree, new_node=None, pks={}):
    possible_nodes= get_possible_regrafters(tree)
    
    assert len(possible_nodes)>0, 'There were no regraft locations possible, which is strange because the root is regraftable and always look the same.'
    
    new_tree= deepcopy(tree)
    regraft_key, regraft_branch= possible_nodes[choice(len(possible_nodes), 1)[0]]
    pks['regraft_key']=regraft_key
    pks['regraft_branch']=regraft_branch
    new_tree, remove_distrub, remove_val, remove_par = remove_parent_attachment(new_tree, regraft_key, regraft_branch)
    q_backward=back_density(remove_distrub, remove_val, remove_par)
    children, other= get_all_branch_descendants_and_rest(new_tree, regraft_key, regraft_branch)
    candidates=thin_out_sibling(new_tree, other, regraft_key)+[('r',0)]
    ch= choice(len(candidates),1)[0]
    recipient_key, recipient_branch=candidates[ch]
    new_tree, q_forward= regraft(new_tree, regraft_key, regraft_branch, recipient_key, new_node=new_node, which_branch=recipient_branch)

    return new_tree, q_forward, q_backward

def thin_out_sibling(tree, branches, key):
    return [(r,w) for r,w in branches if (r!='r' and tree[r][3+w]!='closed_branch' and r!=key)]

def back_density(distrub, val, par):
    if distrub=='r':
        return expon.pdf(val)
    if distrub=='u':#this looks kind of contraintuitive
        return uniform.pdf(val, scale=par)
    
def simulate_and_forward_density(distrub, par=None):
    if distrub == 'r':
        insertion_spot=expon.rvs()
        q=expon.pdf(insertion_spot)
    else:
        insertion_spot=uniform.rvs()
        branch_length=par
        q=uniform.pdf(insertion_spot*branch_length,scale=branch_length)
    return insertion_spot, q

def regraft(tree, remove_key,remove_branch, add_to_branch, new_node=None,which_branch=0):
    
    if add_to_branch=='r':
        insertion_spot, q=simulate_and_forward_density('r')
    else:
        branch_length=get_branch_length(tree, add_to_branch,which_branch)
        #removedprin branch_length
        insertion_spot, q=simulate_and_forward_density('u', branch_length)
    if new_node is None:
        new_node=str(getrandbits(68)).strip()
    tree=graft(tree, remove_key, add_to_branch, insertion_spot, new_node, which_branch, remove_branch=remove_branch)
    return tree,q
