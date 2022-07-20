from Rtree_to_covariance_matrix import leave_node, _thin_out_dic, Population, _add_to_waiting
from Rtree_operations import (node_is_non_admixture, get_leaf_keys, get_real_parents, get_real_children, rename_root, 
                              screen_and_prune_one_in_one_out, remove_non_mixing_admixtures)
from copy import deepcopy

def leave_node(key, node, population, target_nodes, follow_branch):
    if node_is_non_admixture(node): 
        return [follow_branch(parent_key=node[0],branch=0, population=population, target_nodes=target_nodes, child_key=key)]
    else:
        new_pop=population.remove_partition(1.0-node[2])
        return [follow_branch(parent_key=node[0],branch=0, population=population, target_nodes=target_nodes, child_key=key, dependent='none'), #changed dependent='none' to go to most loose restriction that still makes sense. To go back,put dependent=node[1
                follow_branch(parent_key=node[1],branch=1, population=new_pop, target_nodes=target_nodes, child_key=key, dependent='none')]
        
class follow_branch_class(object):
    
    def __init__(self, sub_graph_nodes):
        self.sub_graph_nodes=sub_graph_nodes
        self.seen_merging=False
        
    def __call__(self, parent_key, branch, population, target_nodes, child_key, dependent='none'):
        if self.seen_merging:
            return parent_key, population, dependent
        subset=population.subset_of_the_candidates(self.sub_graph_nodes)
        if subset=='partly':
            target_nodes.append((child_key, branch))
        elif subset=='all':
            self.seen_merging=True
        return parent_key, population, dependent

def get_branches_to_keep(tree, subgraph_keys):
    node_keys=get_leaf_keys(tree)
    pops=[Population([1.0],[node]) for node in node_keys]
    follow_branch=follow_branch_class(subgraph_keys)
    ready_nodes=list(zip(node_keys,pops))
    waiting_nodes={}
    taken_nodes=[]
    target_nodes=[]
    while True:
        for key,pop in ready_nodes:
        
            upds=leave_node(key, tree[key], pop, target_nodes, follow_branch)
            for upd in upds:
                waiting_nodes=_add_to_waiting(waiting_nodes, upd,tree)
            taken_nodes.append(key)
        waiting_nodes,ready_nodes=_thin_out_dic(waiting_nodes, taken_nodes[:])
        if len(ready_nodes)==0:
            return None
        if len(ready_nodes)==1 and ready_nodes[0][0]=="r":
            break

    return target_nodes

def find_root_name(tree):
    parents_seen=set()
    for k in tree:
        ps=get_real_parents(tree[k])
        for p in ps:
            parents_seen.add(p)
    rootset=parents_seen-set(tree.keys())
    return next(iter(rootset))

def prune_to_subtree(tree,branches_to_keep):
    sub_tree={}
    for key,b in branches_to_keep:
        sub_tree[key]=tree[key]
    root_name=find_root_name(sub_tree)
    sub_tree=rename_root(sub_tree, root_name)
    sub_tree=remove_empty_children(sub_tree)
    sub_tree=screen_and_prune_one_in_one_out(sub_tree)
    return sub_tree

def get_subtree(tree, subgraph_keys):
    tree2=remove_non_mixing_admixtures(deepcopy(tree))
    branches_to_keep=get_branches_to_keep(tree2, subgraph_keys)
    return prune_to_subtree(tree2, branches_to_keep)

def remove_empty_children(tree):
    for k in tree:
        child_keys=get_real_children(tree[k])
        children_to_keep=[]
        for child_key in child_keys:
            if child_key in tree:
                children_to_keep.append(child_key)
        if len(child_keys)!=len(children_to_keep):
            for n,ch in enumerate(children_to_keep):
                tree[k][5+n]=ch
            for n in range(n+1,2):
                tree[k][5+n]=None
    return tree
