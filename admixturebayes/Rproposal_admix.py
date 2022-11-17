from numpy.random import choice, random
from copy import deepcopy
from Rtree_operations import (node_is_admixture, insert_admixture_node_halfly, 
                              graft, node_is_non_admixture,
                              parent_is_spouse, halfbrother_is_uncle,
                              parent_is_sibling, other_branch, get_branch_length, change_admixture,
                              get_all_branches, get_all_branch_descendants_and_rest, remove_admix2,
                              readjust_length,get_keys_and_branches_from_children,
                              update_branch_length)
from random import getrandbits
from scipy.stats import expon

import warnings

class addadmix_class(object):
    
    new_nodes=2
    proposal_name='addadmix'
    adaption=False
    input='tree'
    require_admixture=0
    admixture_change=1
    reverse='deladmix'
    
    def __init__(self, **kwargs):
        self.kwargs=kwargs
    
    def __call__(self,*args, **kwargs):
        kwargs.update(self.kwargs)
        return addadmix(*args, **kwargs)
    
class deladmix_class(object):
    
    new_nodes=0
    adaption=False
    proposal_name='deladmix'
    input='tree'
    require_admixture=1
    admixture_change=-1
    reverse='addadmix'
    
    def __init__(self, **kwargs):
        self.kwargs=kwargs
    
    def __call__(self,*args, **kwargs):
        kwargs.update(self.kwargs)
        return deladmix(*args, **kwargs)

def addadmix(tree,new_node_names=None,pks={}, fixed_sink_source=None, new_branch_length=None, new_to_root_length=None, preserve_root_distance=True):
    '''
    This proposal adds an admixture to the tree. There are a lot of free parameters but only 5 are in play here:
        c1: the branch length of the source population
        c2: the branch length of the population genes are migrating into. 
        u1: the position of the admixture source on the source population branch
        u2: the position of the admixture destination on the sink population branch
        w: the admixture proportion.
        The connecting function, h, (see Green & Hastie 2009) is
            h(c1,c2,u1,u2,w)=(c1*u1, c1*(1-u1), c2*u2, c2*(1-u2), 0.5*w)
    '''
    
    possible_nodes=get_all_branches(tree)
        
    new_tree= deepcopy(tree)
    #removedprin possible_nodes
    sink_key, sink_branch=possible_nodes[choice(len(possible_nodes), 1)[0]]
    if fixed_sink_source is not None:
        sink_key,sink_branch,source_key,source_branch = fixed_sink_source
    other= get_all_branch_descendants_and_rest(tree, sink_key, sink_branch)
    candidates=other+[('r',0)]
    ch= choice(len(candidates),1)[0]
    if fixed_sink_source is None:
        source_key, source_branch=candidates[ch]

    if fixed_sink_source is not None:
        new_tree, forward_density, backward_density, multip= insert_admix(new_tree, source_key, source_branch, sink_key, sink_branch, pks=pks, new_branch_length=new_branch_length, new_to_root_length=new_to_root_length, preserve_root_distance=preserve_root_distance)
    elif new_node_names is None:
        new_tree, forward_density, backward_density, multip= insert_admix(new_tree, source_key, source_branch, sink_key, sink_branch, pks=pks, preserve_root_distance=preserve_root_distance)
    else:
        new_tree, forward_density ,backward_density, multip= insert_admix(new_tree, source_key, source_branch, sink_key, sink_branch, pks=pks, source_name=new_node_names[0], sink_name=new_node_names[1], preserve_root_distance=preserve_root_distance)
    
    choices_forward=float(len(possible_nodes)*len(candidates))*2
    choices_backward=float(len(_get_removable_admixture_branches(new_tree)))
    
    return new_tree,forward_density/choices_forward, backward_density/choices_backward*multip 

    
def get_admixture_branch_length(x=None):
    if x is None:
        x=expon.rvs()
        return x, expon.pdf(x)
    else:
        return expon.pdf(x)
    
def get_root_branch_length(x=None):
    return get_admixture_branch_length(x)
    
def get_admixture_proportion(x=None):
    if x is None:
        return random(),1
    else: 
        return 1
    
def get_insertion_spot(x=None, length=1.0):
    if x is None:
        if length<1e-8:
            warnings.warn('length of branch was too close to 0 - incorrect calculation follows')
            return random(), 1.0
        return random(), 1.0/length
    else:
        if length<1e-8:
            warnings.warn('length of branch was too close to 0 - incorrect calculation follows')
            return 1.0
        return 1.0/length
    
def insert_admix(tree, source_key, source_branch, sink_key, sink_branch, source_name=None, sink_name=None, pks={}, new_branch_length=None, new_to_root_length=None, preserve_root_distance=False):
    if source_key=='r':
        u1,q1=get_root_branch_length()
        if new_to_root_length is not None:
            u1,q1 = new_to_root_length, get_root_branch_length(new_to_root_length)
    else:
        u1,q1=get_insertion_spot(length=get_branch_length(tree,source_key,source_branch))
    u2,q2=get_insertion_spot(length=get_branch_length(tree,sink_key,sink_branch))
    if new_branch_length is not None:
        t4,q4= new_branch_length, get_admixture_branch_length(new_branch_length)
    else:
        t4,q4=get_admixture_branch_length()
    u3,q3=get_admixture_proportion()
    if sink_name is None:
        sink_name=str(getrandbits(128)).strip()
    if source_name is None:
        source_name=str(getrandbits(128)).strip()
    tree=insert_admixture_node_halfly(tree, sink_key, sink_branch, u2, admix_b_length=t4, new_node_name=sink_name, admixture_proportion= u3)
    tree=graft(tree, sink_name, source_key, u1, source_name, source_branch, remove_branch=1)

    if preserve_root_distance:
        tree[sink_name], multip=readjust_length(tree[sink_name])
    else:
        multip=1.0
    
    if random()<0.5:
        tree[sink_name]=change_admixture(tree[sink_name])
    return tree,q1*q2*q3*q4,1, multip


def deladmix(tree,pks={}, fixed_remove=None, preserve_root_distance=True):
    '''
    Reversible Jump MCMC transition which removes a random admixture branch. This is the reverse of the other proposal in this file. 
    '''
    
    #making copy that we can erase branches from. 
    cop=deepcopy(tree)
    
    candidates=_get_removable_admixture_branches(cop)
    #removedprin candidates
    if len(candidates)==0:
        return tree,1,1
    if fixed_remove is None:
        remove_key, remove_branch = candidates[choice(len(candidates),1)[0]]
    else:
        remove_key, remove_branch = fixed_remove
    
    new_tree, (t1,t2,t3,t4,t5), alpha = remove_admix2(cop, remove_key, remove_branch, pks=pks)
    
    if preserve_root_distance:
        #removedprin t1
        multip=(alpha**2+(1.0-alpha)**2)
        old_length=t2+multip*t1
        t1=old_length-t2
        #removedprin old_length, t1,t2
        child_key, child_branch= get_keys_and_branches_from_children(tree, remove_key)[0]
        update_branch_length(new_tree, child_key, child_branch, old_length)
    else:
        multip=1.0
    backward_density= get_backward_remove_density(t1,t2,t3,t4,t5, alpha)
    forward_density= 1.0
    
    forward_choices=float(len(candidates))
    backward_choices=float(get_possible_admixture_adds(new_tree, pks['orphanota_key'], pks['orphanota_branch']))*2
    
    return new_tree, forward_density/forward_choices, backward_density/backward_choices*multip

def get_backward_remove_density(t1,t2,t3,t4,t5, alpha):
    '''
    remembering this ugly sketch:
    
                parent_key          sparent_key
                    |                |
                    |t_1             | t_4
                    |   __---- source_key
                  rkey/   t_5       \
                    |                \t_3
                    |t_2          sorphanota_key  
                orphanota_key   

    we want to get the density of all the choices. If 't_4 is None', it is because source_key was the root. So we want to find the density of the insertion spot on the
    parent_key-orphanota_key branch (=u2) the insertion spot on the sorphanota_key-sparent_key branch (=u1) (which could be exponentially distributed. We also want the density of t5
    and the admixture proportion, alpha
    '''
    if t4 is None:
        u1=t3
        q1=get_root_branch_length(u1)
    else:
        q1=get_insertion_spot(t3, t3+t4)
    q2=get_insertion_spot(t2, t1+t2)
    q3=get_admixture_branch_length(t5)
    q4=get_admixture_proportion(alpha)
    
    return q1*q2*q3*q4
    
def get_possible_admixture_adds(tree, sink_key,sink_branch):
    possible_nodes=get_all_branches(tree)
    other= get_all_branch_descendants_and_rest(tree, sink_key, sink_branch)
    candidates=other+[('r',0)]
    return len(possible_nodes)*len(candidates)

def _check_node(tree,key,direction):
    parent_key=tree[key][direction]
    return ((parent_key=='r' or node_is_non_admixture(tree[parent_key])) and 
            not parent_is_spouse(tree, key, other_branch(direction)) and
            (parent_key=='r' or not halfbrother_is_uncle(tree, key, parent_key)) and
            (parent_key=='r' or not (parent_is_spouse(tree,key,direction) and parent_is_sibling(tree, key, direction))))
        
def _get_removable_admixture_branches(tree):
    res=[]
    for key, node in list(tree.items()):
        if node_is_admixture(node):
            if _check_node(tree, key, 0):
                res.append((key,0))
            if _check_node(tree, key, 1):
                res.append((key, 1))
    return res
