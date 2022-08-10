from tree_statistics import (identifier_to_tree_clean, generate_predefined_list_string, unique_identifier_and_branch_lengths)
import pandas as pd
from copy import deepcopy
from Rtree_to_covariance_matrix import leave_node, _thin_out_dic, Population, _add_to_waiting, get_populations

from tree_statistics import admixture_sorted_unique_identifier

from Rtree_operations import (node_is_non_admixture, get_leaf_keys, get_real_parents, get_real_children, rename_root, 
                              screen_and_prune_one_in_one_out, remove_non_mixing_admixtures)

from collections import Counter

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

def get_list_of_turned_topologies(trees, true_tree):
    nodes=get_leaf_keys(true_tree)
    return [admixture_sorted_unique_identifier(tree, nodes) for tree in trees], admixture_sorted_unique_identifier(true_tree, nodes)

def identity(x):
    return x

def iterate_over_output_file(outfile, 
                             cols=[], 
                             pre_thin_data_set_function=identity, 
                             row_summarize_functions=[],
                             thinned_d_dic=identity,
                             **constant_kwargs):
    
    df= pd.read_csv(outfile, usecols=cols, dtype={'no_admixes':object})
    df = df[cols]
    df= pre_thin_data_set_function(df)
    all_results=[]
    
    for n,(i,r) in enumerate(df.iterrows()):
        cont=False
        d_dic={colname:r[k] for k, colname in enumerate(cols)}
        d_dic.update(constant_kwargs)

        for row_summarize_function in row_summarize_functions:
            add_dic, skip=row_summarize_function(**d_dic)
            d_dic.update(add_dic)
        all_results.append(thinned_d_dic(d_dic))
    return all_results

class make_Rtree(object):
    
    def __init__(self, nodes_to_be_sorted, remove_sadtrees=False, subnodes=[], outgroup_name=''):
        self.nodes=sorted(nodes_to_be_sorted)
        self.remove_sadtrees=remove_sadtrees
        self.subnodes=subnodes
        self.outgroup_name=outgroup_name
        
    def __call__(self, tree, **not_needed):
        first_level=tree.split('-')[0]
        no_pops=len(first_level.split('.'))
        if no_pops==len(self.nodes): #checking if
            Rtree=identifier_to_tree_clean(tree, leaves=generate_predefined_list_string(deepcopy(self.nodes)))
        elif no_pops==len(self.nodes)+1 and self.outgroup_name:
            self.nodes=sorted(self.nodes+[self.outgroup_name])
            Rtree = identifier_to_tree_clean(tree, leaves=generate_predefined_list_string(deepcopy(self.nodes)))
        else:
            assert False, 'Either the outgroup name was not specified or something is seriously wrong because' \
                          ' the number of nodes did not match the size of the trees'

        if self.subnodes:#DETTE TAGER IKKE ORDENTLIG HOJDE FOR KOVARIANSMATRICERNE SOM BLIVER FORKERTE
            try:
                Rtree=get_subtree(Rtree, self.subnodes)
            except AssertionError:
                from tree_plotting import plot_as_directed_graph
                plot_as_directed_graph(Rtree)
                print('input_tree', tree)
                print('nodes', self.nodes)
                print('subnodes', self.subnodes)
                assert False
        return {'Rtree':Rtree}, False
    
class make_full_tree(object):
    
    def __init__(self, outgroup_name='out', remove_sadtrees=False, subnodes=[], reroot_population='',
                 reroot_method='stop'):
        self.outgroup_name=outgroup_name
        self.remove_sadtrees=remove_sadtrees
        self.subnodes=subnodes
        self.reroot_population=reroot_population
        self.reroot_method=reroot_method
        
    def __call__(self, Rtree=None, add=None, **kwargs):
        full_tree=deepcopy(Rtree)
        if self.subnodes:
            full_tree=get_subtree(full_tree, self.subnodes)
        return {'full_tree':full_tree}, False

class make_string_tree(object):

    def __init__(self, nodes, tree_unifier=None):
        self.nodes=nodes
        self.tree_unifier=tree_unifier
        self.node_string='='.join(self.nodes)+'='

    def __call__(self, full_tree, **kwargs):
        stree=unique_identifier_and_branch_lengths(full_tree, leaf_order=self.nodes)
        if self.tree_unifier is not None:
            stree=self.tree_unifier(stree)
            
        string_tree=self.node_string+stree
        return {'string_tree':string_tree},  False
    
def get_subpops(pops, sub_graph_keys):
    ss_subgraph_keys=set(sub_graph_keys)
    new_pops=[]
    for pop in pops:
        new_pop=ss_subgraph_keys.intersection(pop.split('.'))
        new_pops.append('.'.join(new_pop))
    if '' in new_pops:
        new_pops.remove('')
    return '_'.join(sorted(list(set(new_pops))))

class topology(object):
    
    def __init__(self, nodes):
        self.nodes=nodes
        
    def __call__(self, Rtree=None, **kwargs):
        if 'string_tree' in kwargs:
            topology=kwargs['string_tree'].split('=')[-1].split(';')[0]
            return {'topology':topology},False
        if Rtree is None:
            Rtree=kwargs['full_tree']
        top=admixture_sorted_unique_identifier(Rtree, leaf_order=self.nodes, not_opposite=True)
        return {'topology':top}, False
    
class get_pops(object):
    
    def __init__(self, min_w=0.0, keys_to_include=None):
        self.min_w=min_w
        self.keys_to_include=keys_to_include
            
    def __call__(self, full_tree=None, **kwargs):
        if full_tree is None:
            tree=kwargs['Rtree']
        else:
            tree=full_tree
        pops=get_populations(tree, self.min_w, keys_to_include=self.keys_to_include)
        return {'pops':'-'.join(pops)}, False

class thinning(object):
    
    def __init__(self, burn_in_fraction=None, total=None, **values_to_filter_by):
        self.burn_in_fraction=burn_in_fraction
        self.total=total
        self.values_to_filter_by=values_to_filter_by
        
    def __call__(self, df):
        #first_removing burn-in
        n=len(df)
        print('Dataframe read with ' + str(n) + ' samples.')
        if self.burn_in_fraction is not None:
            df=df[int(n*self.burn_in_fraction):]
        print('Burn-in of ' + str(n-len(df)) + ' samples removed. There are now ' + str(len(df)) + ' samples.')
        for column,value in list(self.values_to_filter_by.items()):
            print('filtering on ', column,'==',value)
            df=df.loc[df[column]==value,:]
        if self.total is not None:
            n=len(df)
            stepsize = self.total
            df=df[::stepsize]
            print('Thinning every ' + str(stepsize) + ' samples complete. There are now ' + str(len(df)) + ' samples')
        return df

def mode(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]
    
