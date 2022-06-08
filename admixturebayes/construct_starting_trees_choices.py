from generate_prior_trees import generate_phylogeny, generate_add
from tree_statistics import identifier_to_tree_clean, generate_predefined_list_string
from Rtree_operations import create_trivial_tree, scale_tree, rename_leaves
from Rtree_to_covariance_matrix import make_covariance
from rescale_covariance import rescale_empirical_covariance
from copy import deepcopy
from warnings import warn


def get_starting_trees(inputs, 
                       no_chains, 
                       adds=[], 
                       nodes=None, 
                       pipeline=[],
                       multiplier=None,
                       scale_tree_factor=1.0,
                       start='trivial',
                       add_start=None,
                       prefix='',
                       starting_tree_scaling='trivial',
                       starting_tree_use_scale_tree_factor=False,
                       scale_goal='min', 
                       mscale_file=None,
                       no_add=False):
    add_vals=[]
    if adds:
        for add in adds:
            with open(add, 'r') as f:
                add_vals.append(float(f.readline()))
    trees=[]
    for input in inputs:
        trees.append(input_to_tree(input, nodes))
        
    if not trees:
        no_pops=len(nodes)#error if nodes is not specified
        if start=='trivial':
            trees=[]
            for _ in range(no_chains):
                tree=create_trivial_tree(no_pops,1.0)
                rename_leaves(tree, nodes)
                trees.append(tree)
        elif start=='random':
            no_pops=len(nodes) #error if nodes is not specified
            trees=[generate_phylogeny(no_pops, leaf_nodes=nodes, skewed_admixture_prior=False) for _ in range(no_chains)]
        elif start=='perfect':
            tree=get_perfect_tree(prefix, nodes)
            trees=[deepcopy(tree) for _ in range(no_chains)]
        else:
            assert False, 'wrong input to start'
    
    if not add_vals:
        if add_start is None:
            add_start=start
        elif start=='trivial' or no_add:
            add_vals=[0]
        elif start=='random':
            no_pops=len(nodes) #error if nodes is not specified
            add_vals=[generate_add() for _ in range(no_chains)]
        elif start=='perfect':
            add_vals=[get_perfect_add(prefix)]
        else:
            assert False, 'wrong input to add_start'
    
    xs=match_trees_and_adds(trees, add_vals)
    
    if len(xs)==1 and no_chains>1:
        tmp=[deepcopy(xs[0]) for _ in range(no_chains)]
        xs=tmp
        
    if mscale_file is None:
        mscale_file=prefix+'m_scale.txt'
    
    return scale_xs(xs, multiplier, scale_tree_factor, starting_tree_scaling, starting_tree_use_scale_tree_factor, scale_goal, mscale_file=mscale_file)

def retrieve_mscale(mscale_file):
    with open(mscale_file, 'r') as f:
        return float(f.readline())

def scale_xs(xs, multiplier, scale_tree_factor, starting_tree_scaling, starting_tree_use_scale_tree_factor, scale_goal, mscale_file):
    if starting_tree_scaling=='None':
        return xs
    elif starting_tree_scaling=='empirical_trace':
        assert multiplier is not None, 'the empirical trace was not communicated to the starting tree somehow'
        xs=[(scale_tree(tree, multiplier), multiplier*add) for tree,add in xs]
        if starting_tree_use_scale_tree_factor:
            xs=[(scale_tree(tree, scale_tree_factor), add*scale_tree_factor) for tree,add in xs]
        return xs
    elif starting_tree_scaling=='starting_tree_trace':
        new_xs=[]
        for tree, add in xs:
            cov=make_covariance(tree)+add
            _, multiplier = rescale_empirical_covariance(cov, normalizer=scale_goal)
            new_xs.append((scale_tree(tree,multiplier), add*multiplier))
        return new_xs
    elif starting_tree_scaling=='treemix_tree':
        mscale=retrieve_mscale(mscale_file)
        xs=[(scale_tree(tree, multiplier/mscale), multiplier*add/mscale) for tree,add in xs]#multiplying with the 
        return xs
    elif starting_tree_scaling=='scalar':
        assert starting_tree_use_scale_tree_factor, 'Illegal combination of options.'
        xs=[(scale_tree(tree, scale_tree_factor), add*scale_tree_factor) for tree,add in xs]
        return xs
    else:
        assert False, 'unknown input for scale_xs: '+str(starting_tree_scaling)
            

def get_perfect_tree(prefix, nodes):
    with open(prefix+'true_tree.txt', 'r') as f:
        first_string=f.readline().rstrip()
        if ';' in first_string:
            return identifier_to_tree_clean(first_string, leaves=generate_predefined_list_string(deepcopy(nodes)))
        nodes2=first_string.split()
        assert set(nodes)==set(nodes2), 'the perfect tree does not have the correct leaf nodes'
        return identifier_to_tree_clean(f.readline().rstrip(), leaves=generate_predefined_list_string(deepcopy(nodes)))
    
def get_perfect_add(prefix):
    with open(prefix+'true_add.txt', 'r') as f:
        a=f.readline().rstrip()
        return float(a)
        
        
def match_trees_and_adds(list_of_trees, list_of_adds):
    if len(list_of_adds)==0:
        return [(t,0) for t in list_of_trees]
    elif len(list_of_adds)==1:
        return [(t, list_of_adds[0]) for t in list_of_trees]
    elif len(list_of_adds)==len(list_of_trees):
        return [(t,a) for t,a in zip(list_of_trees, list_of_adds)]
    else:
        assert False, 'couldnt match adds and trees in starting_trees'
        
        
            
            
def fill_up_list(list_of_trees, no_chains):
    list_of_trees=[]
    if len(list_of_trees)>no_chains:
        warn('Some starting trees will not be used')
    trees=[]
    count=0
    for _ in range(no_chains):
        trees.append(list_of_trees[count])
        count+=1
        if count==len(list_of_trees):
            count=0
    return trees
        

            
def input_to_tree(input, nodes, skewed_admixture_prior=False):
    assert isinstance(input, str), 'unrecognized input:'+ str(input)+ ' has type '+ str(type(input))
    if ';' in input:
        return identifier_to_tree_clean(input, nodes=nodes)
    elif '.' in input:
        with open(input, 'r') as f:
            f.readline() #removing empty file
            return identifier_to_tree_clean(f.readline().rstrip(), leaves=generate_predefined_list_string(deepcopy(nodes)))
    elif ',' in input:
        first_part, last_part= input.split(',')
        n=int(first_part[1:])
        k=int(last_part[:-1])
        return generate_phylogeny(n,k, leaf_nodes=nodes, skewed_admixture_prior=skewed_admixture_prior)
    else:#tries to open file anyway
        with open(input, 'r') as f:
            f.readline() #removing empty file
            return identifier_to_tree_clean(f.readline().rstrip(), leaves=generate_predefined_list_string(deepcopy(nodes)))
