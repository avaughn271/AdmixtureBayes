from scipy.stats import chi2, uniform, expon
from numpy.random import choice
from copy import deepcopy

from Rtree_operations import (node_is_non_admixture, get_parent_of_branch, move_node, find_rooted_nodes, get_real_children, mother_or_father,
graft, get_real_parents, halfbrother_is_uncle,
get_branch_length, get_all_branch_descendants_and_rest)
from random import getrandbits

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
        node=tree[self.key]
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

class sliding_regraft_class_resimulate(object):
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
    
def simulate_regraft_distance(param):
    return chi2.rvs(1)*param

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
    
    distance_to_regraft= simulate_regraft_distance(param)
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