from Rproposal_regraft import get_possible_regrafters, thin_out_sibling
from Rtree_operations import get_all_branch_descendants_and_rest, get_parent_of_branch, move_node, get_branch_length, update_branch_length
from nearest_branch_search import distanced_branch_lengths
from scipy.stats import chi2
from numpy.random import choice
from copy import deepcopy
from random import getrandbits
from math import exp
from scipy.stats import gamma, uniform

def conditional_gamma_rvs(mean,shape, upper_limit):
    shape,scale=transform_to_shape_scale(mean, shape=shape)
    max_U= gamma.cdf(upper_limit, a=shape, scale=scale)
    drawn_U= uniform.rvs()*max_U
    x=gamma.ppf(drawn_U, a=shape, scale=scale)        
    return gamma.ppf(drawn_U, a=shape, scale=scale)

def conditional_gamma_logpdf(value, mean, shape, upper_limit):
    shape,scale=transform_to_shape_scale(mean, shape=shape)
    logprob_conditional= gamma.logcdf(upper_limit, a=shape, scale=scale)
    return gamma.logpdf(value, a=shape, scale=scale) - logprob_conditional

def rvs(t,delta_L, shape=20):
    if delta_L<0:
        mean= -delta_L
        return t-conditional_gamma_rvs(mean, shape=shape, upper_limit=t)
    else:
        return t+clean_gamma_rvs(delta_L, shape=shape)

def logpdf(x,t,delta_L, shape=20):
    if delta_L < 0:
        mean=-delta_L
        simulated_part=t-x
        return conditional_gamma_logpdf(simulated_part, mean=mean, shape=shape, upper_limit=t)
    else:
        simulated_part=x-t
        return clean_gamma_logpdf(simulated_part, mean=delta_L, shape=shape)
        
    
def clean_gamma_rvs(mean, shape):
    shape,scale=transform_to_shape_scale(mean, shape=shape)
    return gamma.rvs(a=shape, scale=scale)

def clean_gamma_logpdf(value, mean, shape):
    shape,scale=transform_to_shape_scale(mean, shape=shape)
    return gamma.logpdf(value, a=shape, scale=scale)
    
def transform_to_shape_scale(mean, variance=None, shape=None):
    if shape is None:
        scale=variance/mean
        shape=mean/scale
    else:
        scale=mean/shape
    return shape,scale

class sliding_regraft_class_resimulate(object):
    new_nodes=1
    input='tree'
    require_admixture=0
    reverse_require_admixture=0
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
    children, other= get_all_branch_descendants_and_rest(tree, regraft_key, regraft_branch)
    candidates= thin_out_sibling(tree, other, regraft_key)+[('r',0)]
    thinned_pieces=[piece for piece in pieces if (piece.get_branch_key() in candidates and piece.contains_distance(distance_to_regraft) )]
    return thinned_pieces
    
    
def make_sliding_regraft(tree, new_node=None, param=0.1, resimulate_moved_branch_length=False, pks={}):
    '''
    
    '''
    
    possible_nodes= get_possible_regrafters(tree)
        
    new_tree= deepcopy(tree)
    regraft_key, regraft_branch= possible_nodes[choice(len(possible_nodes), 1)[0]]
    pks['regraft_key']=regraft_key
    pks['regraft_branch']=regraft_branch
    pks['sliding_regraft_adap_param']= param
    
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
    
    
    if resimulate_moved_branch_length:
        new_tree, forward, backward= resimulate_moved_branch(new_tree, regraft_key, regraft_branch, chosen_piece.get_lattitude(distance_to_regraft), resimulate_moved_branch_length)
    else:
        forward=1.0
        backward=1.0
    
    return new_tree, forward/forward_choices, backward/backward_choices

def resimulate_moved_branch(tree, key, branch, delta_L, alpha):
    old_length= get_branch_length(tree, key, branch)
    new_length= gamma_restricted.rvs(old_length, -delta_L, alpha)
    logpdf=gamma_restricted.logpdf(new_length, old_length, -delta_L, alpha)
    logpdf_back=gamma_restricted.logpdf(old_length, new_length, delta_L, alpha)
    backward=exp(-logpdf+logpdf_back) #for stability, both forward and backward are put in backward. 
    update_branch_length(tree, key, branch, new_length)
    return tree, 1.0, backward
