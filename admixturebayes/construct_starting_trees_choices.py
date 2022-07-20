from generate_prior_trees import generate_phylogeny, generate_add
from tree_statistics import identifier_to_tree_clean, generate_predefined_list_string
from Rtree_operations import create_trivial_tree, rename_leaves
from copy import deepcopy

def get_starting_trees(inputs, 
                       no_chains, 
                       adds=[], 
                       nodes=None, 
                       start='trivial',
                       add_start=None):
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
    
    if not add_vals:
        if add_start is None:
            add_start=start
        elif start=='trivial':
            add_vals=[0]
        elif start=='random':
            no_pops=len(nodes) #error if nodes is not specified
            add_vals=[generate_add() for _ in range(no_chains)]
        else:
            assert False, 'wrong input to add_start'
    
    xs=match_trees_and_adds(trees, add_vals)
    
    if len(xs)==1 and no_chains>1:
        tmp=[deepcopy(xs[0]) for _ in range(no_chains)]
        xs=tmp
    
    return xs

def match_trees_and_adds(list_of_trees, list_of_adds):
    if len(list_of_adds)==0:
        return [(t,0) for t in list_of_trees]
    elif len(list_of_adds)==1:
        return [(t, list_of_adds[0]) for t in list_of_trees]
    elif len(list_of_adds)==len(list_of_trees):
        return [(t,a) for t,a in zip(list_of_trees, list_of_adds)]
    else:
        assert False, 'couldnt match adds and trees in starting_trees'
 
def input_to_tree(input, nodes):
    with open(input, 'r') as f:
        f.readline() #removing empty file
        return identifier_to_tree_clean(f.readline().rstrip(), leaves=generate_predefined_list_string(deepcopy(nodes)))