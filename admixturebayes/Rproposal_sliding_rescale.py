from Rtree_operations import get_parent_of_branch, get_all_branches
from nearest_branch_search import distanced_branch_lengths
from numpy.random import normal, choice
from copy import deepcopy

class sliding_rescale_class(object):
    
    new_nodes=0
    proposal_name='sliding_rescale'
    adaption=True
    input='tree'
    require_admixture=0
    reverse_require_admixture=0
    admixture_change=0
    reverse='sliding_rescale'
    
    def __init__(self):
        self.recently_called_number=None
    
    def __call__(self,*args, **kwargs):
        res, self.recently_called_number = make_sliding_rescale(*args, **kwargs)
        return res
    
def simulate_regraft_distance(max_top_dist):
    c=choice(max_top_dist+1)
    return c

def get_thinned_pieces(tree,rescale_key, rescale_branch, cutoff_distance, parent_key):
    pieces= distanced_branch_lengths(tree, parent_key, visited_keys=[], upper_cap=cutoff_distance, topological_distance=True)
    thinned_pieces=[piece for piece in pieces if  piece.within_distance(cutoff_distance)]
    return thinned_pieces    

def resimulate_pieces(tree, pieces, max_distance, param, admixture_multiplier=1.0, pks={}):
    visited_keys=[]
    for piece in pieces:
        if piece.child_key!='r':
            visited_keys.append(piece.parent_key)
            distance=piece.get_start_distance()
            if max_distance==0:
                d=1
            else:
                d=(max_distance-distance)/max_distance
    
            if piece.child_branch==1:
                tree[piece.child_key][2]+= normal()*param*admixture_multiplier*d
            tree[piece.child_key][3+piece.child_branch]+= normal()*param*d
    return tree
                
        
    
def make_sliding_rescale(tree, param=0.1,pks={}):
    '''
    
    '''
    
    possible_branches= get_all_branches(tree)
        
    new_tree= deepcopy(tree)
    rescale_key, rescale_branch= possible_branches[choice(len(possible_branches), 1)[0]]
    pks['rescale_key']=rescale_key
    pks['rescale_branch']=rescale_branch
    pks['sliding_rescale_adap_param']= param
    
    cutoff_distance= simulate_regraft_distance(3-1)
    parent_key= get_parent_of_branch(tree, rescale_key, rescale_branch)
    thinned_pieces_forward=get_thinned_pieces(tree, rescale_key, rescale_branch, cutoff_distance, parent_key)
    
    pks['cutoff_distance']= cutoff_distance
    pks['number_of_pieces']=len(thinned_pieces_forward)
    
    
    new_tree=resimulate_pieces(new_tree, thinned_pieces_forward, cutoff_distance, param, pks=pks)

    
    forward=1.0
    backward=1.0
    
    return (new_tree, forward, backward), cutoff_distance

    