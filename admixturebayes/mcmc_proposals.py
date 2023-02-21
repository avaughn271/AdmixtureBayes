from numpy.random import choice, random, normal
from copy import deepcopy
from Rtree_operations import (node_is_admixture, insert_admixture_node_halfly,  find_rooted_nodes, get_number_of_leaves,
                              graft, node_is_non_admixture, get_real_parents, move_node, get_number_of_admixes,
                              parent_is_spouse, halfbrother_is_uncle, get_parent_of_branch, update_specific_branch_lengths,
                              parent_is_sibling, other_branch, get_branch_length, change_admixture,
                              get_all_branches, get_all_branch_descendants_and_rest, remove_admix2,
                              readjust_length,get_keys_and_branches_from_children, update_all_branches,
                              update_branch_length, get_real_children, mother_or_father, get_leaf_keys)
from random import getrandbits

import warnings
from math import sqrt
from scipy.stats import chi2, uniform, expon

from numpy import zeros, insert
from Rtree_to_covariance_matrix import Population, _add_to_waiting, _thin_out_dic
from scipy.linalg import svd

def rescale(add, sigma=0.01, pks={}):
    new_add=add+normal()*sigma
    if new_add < 0:
        new_add = (-1) * new_add
    return new_add,1,1

class rescale_add_class(object):
    new_nodes=0
    proposal_name='rescale_add'
    adaption=True
    input='add'
    require_admixture=0
    reverse='rescale_add'
    admixture_change=0
    
    def __call__(self,*args, **kwargs):
        return rescale(*args, **kwargs)

class Coefficient_Matrix():
    
    def __init__(self, nodes_to_index, branch_to_index, custom_list={}):
        self.ni=nodes_to_index
        self.bi=branch_to_index
        self.custom_list=custom_list
        if custom_list:
            self.cofmatr=zeros((len(custom_list), len(branch_to_index)))
        else:
            self.cofmatr=zeros((len(nodes_to_index), len(branch_to_index)))
        
    def update(self, branch, population):
        j=self.bi[branch]
        if self.custom_list:
            for pop1_index in range(len(population.members)):
                for pop2_index in range(pop1_index, len(population.members)):
                    pop1=population.members[pop1_index]
                    pop2=population.members[pop2_index]
                    key=None
                    if (pop1,pop2) in self.custom_list:
                        key=(pop1,pop2)
                    elif (pop2,pop1) in self.custom_list:
                        key=(pop2,pop1)
                    if key is not None:
                        self.cofmatr[self.custom_list[key],j]=population.weights[pop1_index]*population.weights[pop2_index]
        else:
            for pop, w in zip(population.members, population.weights):
                self.cofmatr[self.ni[pop],j]=w**2

def nullspace(A):
    u, s, vh = svd(A)
    nnz = (s >= 1e-13).sum()
    ns = vh[nnz:].conj().T
    return ns

def get_orthogonal_branch_space(tree, add_one_column=True):
    cof,bi= make_coefficient_matrix(tree)
    if add_one_column:
        cof=insert(cof, cof.shape[1], 1, axis=1)
    ad=nullspace(cof)
    return ad, bi
    
def make_coefficient_matrix(tree):
    '''
    Instead of constructing the covariance matrix, this function calculates the coefficient matrix, C, to solve
    
    w=Cx
    
    where w is the diagonal of the covariance matrix and x is the vector of branch lengths. Hence, C depends on the admixture proportions and the topology.
    '''
    node_keys=sorted(get_leaf_keys(tree))
    branch_keys=get_all_branches(tree)
    pops=[Population([1.0],[node]) for node in node_keys]
    ready_nodes=list(zip(node_keys,pops))
    ni={node_key:n for n,node_key in enumerate(node_keys)}
    bi={branch:n for n,branch in enumerate(branch_keys)}
    cofmat=Coefficient_Matrix(ni,bi, get_all_pairs(node_keys))
    waiting_nodes={}
    taken_nodes=[]
    while True:
        for key,pop in ready_nodes:
            upds=leave_node(key, tree[key], pop, cofmat)
            for upd in upds:
                waiting_nodes=_add_to_waiting(waiting_nodes, upd,tree)
            taken_nodes.append(key)
        waiting_nodes,ready_nodes=_thin_out_dic(waiting_nodes, taken_nodes[:])
        if len(ready_nodes)==0:
            return None
        if len(ready_nodes)==1 and ready_nodes[0][0]=="r":
            break

    return cofmat.cofmatr,bi

def leave_node(key, node, population, cofmat):
    if node_is_non_admixture(node):
        return [follow_branch(parent_key=node[0],branch=(key,0), population=population, cofmat=cofmat)]
    else:
        new_pop=population.remove_partition(1.0-node[2])
        return [follow_branch(parent_key=node[0],branch=(key,0), population=population, cofmat=cofmat, dependent='none'), #changed dependent='none' to go to most loose restriction that still makes sense. To go back,put dependent=node[1
                follow_branch(parent_key=node[1],branch=(key,1), population=new_pop, cofmat=cofmat, dependent='none')]

def follow_branch(parent_key, branch, population, cofmat, dependent="none"):
    cofmat.update(branch, population)
    return parent_key, population, dependent

def get_all_pairs(nodes):
    all_nodes={}
    count=0
    for node_index in range(len(nodes)):
        for node_index2 in range(node_index, len(nodes)):
            node1=nodes[node_index]
            node2=nodes[node_index2]
            all_nodes[(node1, node2)]=count
            count+=1
    return all_nodes
    
def reverse_dic_to_list(dic):
    l_determined=[None]*len(dic)
    for e,i in list(dic.items()):
        l_determined[i]=e
    return l_determined

def get_added_branch_pieces(org, param):
    return org.dot(normal(scale=param, size=org.shape[1]))

def rescale_constrained(x, sigma=0.01, pks={}, update_add=True):
    tree, add=x
    new_tree=deepcopy(tree)
    orgs, bi= get_orthogonal_branch_space(new_tree, add_one_column=update_add)
    branches=reverse_dic_to_list(bi)
    branch_pieces= get_added_branch_pieces(orgs, sigma)
    if update_add:
        b=branch_pieces[:-1]
    else:
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
    admixture_change=0
    reverse='rescale_constrained'
    
    def __init__(self, **kwargs):
        self.kwargs=kwargs
    
    def __call__(self,*args, **kwargs):
        kwargs.update(self.kwargs)
        return rescale_constrained(*args, **kwargs)

def rescale_admixtures(tree, sigma=0.01, pks={}):
    new_tree=deepcopy(tree)
    if get_number_of_admixes(tree) > 0:
        for key, node in list(new_tree.items()):
            if node_is_admixture(node):
                node[2]=uniform.rvs()
    else:
        return new_tree,1,0.234 #not to have to deal with the admix=0 case, I return 0.234 such that the adaptive parameter is not affected by these special cases.
    return new_tree ,1,1

class rescale_admixtures_class(object):
    new_nodes=0
    proposal_name='rescale_admixtures'
    adaption=True
    input='tree'
    require_admixture=1
    admixture_change=0
    reverse='rescale_admixtures'
    
    def __call__(self, *args, **kwargs):
        return rescale_admixtures(*args, **kwargs)

class updater(object):
    
    def __init__(self, sigma):
        self.sigma=sigma

    def __call__(self):
        return normal(scale=self.sigma)

def rescale_normal(tree, sigma=0.01, pks={}):
    n=get_number_of_leaves(tree)
    k=get_number_of_admixes(tree)
    new_tree=deepcopy(tree)
    updat=updater(sigma/sqrt(2*n-2+4*k))
    new_tree=update_all_branches(new_tree, updat)
    return new_tree ,1,1

class rescale_class(object):
    new_nodes=0
    proposal_name='rescale'
    adaption=True
    input='tree'
    require_admixture=0
    reverse='rescale'
    admixture_change=0
    
    def __call__(self,*args, **kwargs):
        return rescale_normal(*args, **kwargs)

def get_possible_regrafters(tree):
    res=[]
    for key in tree:
        parents=get_real_parents(tree[key])
        for branch, parent in enumerate(parents):
            if parent=='r' or (node_is_non_admixture(tree[parent]) and not halfbrother_is_uncle(tree, key, parent)):
                res.append((key,branch))
    return res

def thin_out_sibling(tree, branches, key):
    return [(r,w) for r,w in branches if (r!='r' and tree[r][3+w]!='closed_branch' and r!=key)]

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
        insertion_spot, q=simulate_and_forward_density('u', branch_length)
    if new_node is None:
        new_node=str(getrandbits(68)).strip()
    tree=graft(tree, remove_key, add_to_branch, insertion_spot, new_node, which_branch, remove_branch=remove_branch)
    return tree,q

class piece(object):
    
    def __init__(self, start_lattitude, end_lattitude, start_distance, end_distance, child_key, child_branch, parent_key):
        self.start_lattitude=start_lattitude
        self.end_lattitude=end_lattitude
        self.start_distance=start_distance
        self.end_distance=end_distance
        self.child_key=child_key
        self.child_branch=child_branch
        self.parent_key=parent_key
        if end_lattitude is not None and start_lattitude>end_lattitude:
            self.direction='to_leaves'
        else:
            self.direction='to_root'
        
    def __str__(self):
        return ', '.join(map(str,[(self.child_key,self.child_branch, self.parent_key), self.start_lattitude, self.end_lattitude, self.start_distance, self.end_distance, self.direction]))
    
    def get_branch_key(self):
        return (self.child_key, self.child_branch)
    
    def contains_distance(self, distance):
        if self.end_distance is None:
            return distance>=self.start_distance
        return distance<=self.end_distance and distance>=self.start_distance
    
    def get_leaf_and_root_sided_length(self, distance):
        if self.end_distance is None:
            return distance-self.start_distance, None
        if self.direction=='to_leaves':
            return self.end_distance-distance, distance-self.start_distance
        else:
            return distance-self.start_distance, self.end_distance-distance

class lineage(object):
    
    def __init__(self, key,  distance=0, lattitude=0, topological_distance=False):
        self.key,self.distance, self.lattitude= key, distance, lattitude
        if topological_distance:
            self.get_branch_length=1
        else:
            self.get_branch_length=get_branch_length
        self.topological_distance=topological_distance
        
    def follow(self, tree, visited_keys=[]):
        new_lineages=[]
        pieces=[]
        for key in get_real_children(tree[self.key]):
            if key not in visited_keys:
                branch=mother_or_father(tree, child_key=key, parent_key=self.key)
                l=get_branch_length(tree, key, branch)
                pieces.append(piece(self.lattitude, self.lattitude-l, self.distance, self.distance+l, key, branch, self.key))
                new_lineages.append(lineage(key,self.distance+l, self.lattitude-l, topological_distance=self.topological_distance))
        for key in get_real_parents(tree[self.key]):
            if key not in visited_keys:
                branch=mother_or_father(tree, child_key=self.key, parent_key=key)
                l=get_branch_length(tree, self.key, branch)
                pieces.append(piece(self.lattitude, self.lattitude+l, self.distance, self.distance+l, self.key, branch, key))
                new_lineages.append(lineage(key,self.distance+l, self.lattitude+l, topological_distance=self.topological_distance))
        if self.key=='r':#add the very long piece
            pieces.append(piece(self.lattitude, None, self.distance, None,'r',0,None))
        return new_lineages, pieces

    def under_cap(self, cap):
        return self.distance<cap

def insert_root_in_tree(tree):
    (child_key1,_,_),(child_key2,_,_)=find_rooted_nodes(tree)
    tree['r']=[None,None,None,None,None,child_key1, child_key2]
    
def distanced_branch_lengths(tree, start_key, visited_keys=[], upper_cap=float('inf'), topological_distance=False):
    insert_root_in_tree(tree)
    pieces=[]
    lineages=[lineage(start_key, 0, 0, topological_distance=topological_distance)]
    while lineages:
        lineages.sort(key=lambda x: x.distance)
        lin=lineages.pop(0)
        if lin.key not in visited_keys:
            visited_keys.append(lin.key)
            new_lineages, new_pieces= lin.follow(tree, visited_keys)
            pieces.extend(new_pieces)
            lineages.extend([new_lineage for new_lineage in new_lineages if new_lineage.under_cap(upper_cap)])
    del tree['r']
    return pieces

class sliding_regraft_class(object):
    new_nodes=1
    input='tree'
    require_admixture=0
    admixture_change=0
    proposal_name='sliding_regraft'
    adaption=True
    reverse='sliding_regraft'
    
    def __init__(self, param=False):
        self.param=param
    
    def __call__(self,*args, **kwargs):
        return make_sliding_regraft(*args, resimulate_moved_branch_length=self.param, **kwargs)

def get_thinned_pieces(tree,regraft_key, regraft_branch, distance_to_regraft, parent_key):
    pieces= distanced_branch_lengths(tree, parent_key, visited_keys=[regraft_key], upper_cap=distance_to_regraft)
    other= get_all_branch_descendants_and_rest(tree, regraft_key, regraft_branch)
    candidates= thin_out_sibling(tree, other, regraft_key)+[('r',0)]
    thinned_pieces=[piece for piece in pieces if (piece.get_branch_key() in candidates and piece.contains_distance(distance_to_regraft) )]
    return thinned_pieces
    
def make_sliding_regraft(tree, new_node=None, param=0.1, resimulate_moved_branch_length=False, pks={}):
    
    possible_nodes= get_possible_regrafters(tree)
        
    new_tree= deepcopy(tree)
    regraft_key, regraft_branch= possible_nodes[choice(len(possible_nodes), 1)[0]]
    
    distance_to_regraft = chi2.rvs(1)*param
    parent_key= get_parent_of_branch(tree, regraft_key, regraft_branch)
    thinned_pieces_forward=get_thinned_pieces(tree, regraft_key, regraft_branch, distance_to_regraft, parent_key)
        
    forward_choices=len(thinned_pieces_forward)
    chosen_piece=thinned_pieces_forward[choice(len(thinned_pieces_forward),1)[0]]
    
    if new_node is None:
        new_node=str(getrandbits(68)).strip()
    new_tree =move_node(new_tree, regraft_key, regraft_branch, parent_key, distance_to_regraft, chosen_piece, new_node_name=new_node)
    
    parent_key= get_parent_of_branch(new_tree, regraft_key, regraft_branch)
    thinned_pieces_backward=get_thinned_pieces(new_tree, regraft_key, regraft_branch, distance_to_regraft, parent_key)
    
    backward_choices=len(thinned_pieces_backward)
    
    return new_tree, 1.0/forward_choices, 1.0/backward_choices

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
