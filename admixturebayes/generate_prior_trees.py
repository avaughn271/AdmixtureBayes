from numpy.random import choice
from Rtree_operations import rename_root, get_number_of_leaves
from collections import Counter
from scipy.stats import uniform, expon, geom


def rvss(fro=0.0, to=1.0):
    guess=uniform.rvs(fro, to-fro)
    while uniform.rvs() > (to-(guess-fro))/(to-fro):
        guess=uniform.rvs(fro, to-fro)
    return guess

def set_outgoing_branch(node, parent_name, branch, length):
    node[branch]=parent_name
    node[branch+3]=length
    return node

factor_to_index_number={'parent_key':0, 
                        'parent_key2':1,
                        'admixture_proportion':2,
                        'branch_length':3,
                        'branch_length2':4,
                        'child_key':5,
                        'child_key2':6}

def create_node(**kwargs):
    return update_node([None]*7, **kwargs)

def update_node(node, **kwargs):
    for factor, value in list(kwargs.items()):
        node[factor_to_index_number[factor]]=value
    return node

def generate_admix_topology(size, admixes, leaf_nodes=None):
    if leaf_nodes is None:
        leaf_nodes = [ 's'+str(i+1) for i in range(size)]
    free_admixes=admixes
    no_totally_free_coalescences=size-1+admixes
    
    ready_lineages=[(leaf_node,0) for leaf_node in leaf_nodes]
    tree={key:[None]*7 for key in leaf_nodes}
    halfly_free_coalescences=[]
    
    node_name=_get_node_name()
    
    while True:
        ready_lineages, tree, no_totally_free_coalescences, halfly_free_coalescences, free_admixes = simulate_generation(no_totally_free_coalescences, 
                                                                                                                         halfly_free_coalescences, 
                                                                                                                         free_admixes, 
                                                                                                                         ready_lineages, 
                                                                                                                         tree, 
                                                                                                                         node_name)
        if no_totally_free_coalescences+len(halfly_free_coalescences)+free_admixes==0:
            break
        
    for key, node in list(tree.items()):
        if node[0] is None and node[1] is None:
            del tree[key]
            tree=rename_root(tree,key)
    return tree

def simulate_number_of_admixture_events(p=0.5):
    return geom.rvs(p=p)-1

def get_admixture_factor(n,k):
    return float(2*n-2+3*k)/float(2*n-2)

def generate_phylogeny(size,admixes=None, p=0.5, leaf_nodes=None, skewed_admixture_prior=False):
    if admixes is None:
        admixes=simulate_number_of_admixture_events(p)
    tree=generate_admix_topology(size, admixes, leaf_nodes)
    n=get_number_of_leaves(tree)
    factor=get_admixture_factor(n, admixes)
    for node in list(tree.values()):
        node=_resimulate(node, factor, skewed_admixture_prior)
    return tree

def generate_add():
    return expon.rvs()


def _resimulate(node, factor=1.0, skewed_admixture_prior=False):
    if node[2] is not None:
        if skewed_admixture_prior:
            node[2]=rvss()
        else:
            node[2]=uniform.rvs()
    if node[3] is not None:
        node[3]=expon.rvs()/factor
    if node[4] is not None:
        node[4]=expon.rvs()/factor
    return node

def _allowed_generation(chosen_indexes, no_totally_free, no_halfly_frees, no_admixes, illegal_indexes):
    
    #checking if there are double bands
    for i1,i2 in illegal_indexes:
        c1,c2=chosen_indexes[i1], chosen_indexes[i2]
        if c2==c1+1 and c1%2==0 and c2<no_totally_free:
            return False
        if c1==c2+1 and c2%2==0 and c1<no_totally_free:
            return False
    
    no_doubles=0
    tmp=sorted([c for c in chosen_indexes if c<no_totally_free])
    if len(tmp)>=2:
        for c1,c2 in zip(tmp[:-1],tmp[1:]):
            if (c2==c1+1 and c1%2==0):
                no_doubles+=1
    no_admixtures=sum(chosen_index>=no_totally_free+no_halfly_frees for chosen_index in chosen_indexes)
    no_singles=sum(chosen_index>=no_totally_free for chosen_index in chosen_indexes)-no_admixtures
    if no_admixtures+no_singles+no_doubles==0:
        return False
    if no_admixes==1 and no_admixtures==0 and ((no_halfly_frees==0 and no_totally_free==4 and no_doubles==1) or
                                               (no_halfly_frees==1 and no_totally_free==2 and no_singles==1)):
        return False
    
    return True

class _get_node_name(object):
    
    def __init__(self, letter='n', admixture_letter='a'):
        self.letter=letter
        self.admixture_letter=admixture_letter
        self.counter=0
        
    def __call__(self, admixture=False):
        self.counter+=1
        if admixture:
            return self.admixture_letter+str(self.counter)
        return self.letter+str(self.counter)
    
def _pair_everyhting_up_nicely(indexes, no_totally_free, halfly_frees, no_admixes, node_name):
    parents=[]
    types=[]
    seen_parents={}
    for n,index in enumerate(indexes):
        type_of_index= _classify_type(index, no_totally_free*2, len(halfly_frees), no_admixes)
        if type_of_index=='free':
            if (index-1 in seen_parents and index%2==1 and index<no_totally_free*2): 
                parents.append(seen_parents[index-1])
                types.append('second_coalescence')
            elif (index+1 in seen_parents and index%2==0 and index+1<no_totally_free*2):
                parents.append(seen_parents[index+1])
                types.append('second_coalescence')
            else:
                parent_key=node_name()
                parents.append(parent_key)
                types.append('first_coalescence')
                seen_parents[index]=parent_key
        elif type_of_index=='admix':
            parents.append(node_name(admixture=True))
            types.append('admixture')
        else:
            which_suitor=index-no_totally_free*2
            parents.append(halfly_frees[which_suitor])
            types.append('second_coalescence')
    return parents, types

def _get_illegal_indexes(lineages):
    keys, _ =list(zip(*lineages))
    dic= Counter(keys)
    res=[]
    for key in keys:
        if dic[key]==2:
           res.append((lineages.index((key,0)), lineages.index((key,1))))
    return res 
            
def simulate_generation(no_totally_free, halfly_frees, no_admixes, lineages, tree, node_name):
    no_halfly_frees=len(halfly_frees)
    new_lineages=[]
    indexes=choice(no_totally_free*2+no_halfly_frees+no_admixes, size = len(lineages), replace=False )
    illegal_indexes=_get_illegal_indexes(lineages)
    COUNT=0
    while (not _allowed_generation(indexes, no_totally_free*2, no_halfly_frees, no_admixes, illegal_indexes) and COUNT<100):
        indexes=choice(no_totally_free*2+no_halfly_frees+no_admixes, size = len(lineages), replace=False)
    parent_keys, types= _pair_everyhting_up_nicely(indexes, no_totally_free, halfly_frees, no_admixes, node_name)
    new_lineages=[]
    for (key,branch), parent_key, typ in zip(lineages, parent_keys, types):
        if typ=='first_coalescence':
            tree[key]=set_outgoing_branch(tree[key], parent_key, branch, 0.14)
            no_totally_free-=1
            halfly_frees.append(parent_key)
            tree[parent_key]=create_node(child_key=key)
        if typ=='second_coalescence':
            tree[key]=set_outgoing_branch(tree[key], parent_key, branch, 0.13)
            halfly_frees.remove(parent_key)
            tree[parent_key]=update_node(tree[parent_key], child_key2=key)
            new_lineages.append((parent_key,0))
        if typ=='admixture':
            no_admixes-=1
            tree[key]=set_outgoing_branch(tree[key], parent_key, branch, 0.12)
            tree[parent_key]=create_node(admixture_proportion=0.6, child_key=key)
            new_lineages.append((parent_key,0))
            new_lineages.append((parent_key,1))
    return new_lineages, tree, no_totally_free, halfly_frees, no_admixes
                    
def _classify_type(index, n_frees, n_halfs, n_admixs):
    if index<n_frees:
        return 'free'
    elif index<n_halfs+n_frees:
        return 'half'
    return 'admix'
