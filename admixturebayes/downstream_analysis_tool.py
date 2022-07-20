from Rtree_operations import add_outgroup
from tree_statistics import (identifier_to_tree_clean, generate_predefined_list_string, unique_identifier_and_branch_lengths)
from Rtree_to_covariance_matrix import get_populations
import pandas as pd
from copy import deepcopy
from tree_to_data import file_to_emp_cov

from tree_statistics import admixture_sorted_unique_identifier
from Rtree_operations import get_leaf_keys

from collections import Counter
from subgraphing import get_subtree

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
    
    def __init__(self, add_multiplier=1, outgroup_name='out', remove_sadtrees=False, subnodes=[], reroot_population='',
                 reroot_method='stop'):
        self.add_multiplier=add_multiplier
        self.outgroup_name=outgroup_name
        self.remove_sadtrees=remove_sadtrees
        self.subnodes=subnodes
        self.reroot_population=reroot_population
        self.reroot_method=reroot_method
        
    def __call__(self, Rtree=None, add=None, **kwargs):
        if Rtree is None:
            assert 'sfull_tree' in kwargs, 'sfull_tree not specified'
            nodes=sorted(kwargs['full_nodes'])
            sfull_tree=kwargs['sfull_tree']
            full_tree=identifier_to_tree_clean(sfull_tree, leaves=generate_predefined_list_string(deepcopy(nodes)))
            if self.subnodes:
                full_tree=get_subtree(full_tree, self.subnodes)
            return {'full_tree':full_tree}, False
        if self.outgroup_name and self.outgroup_name not in Rtree:
            full_tree= add_outgroup(deepcopy(Rtree),
                                    inner_node_name='new_node',
                                    to_new_root_length=float(add)*self.add_multiplier,
                                    to_outgroup_length=0,
                                    outgroup_name=self.outgroup_name)
        else:
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
    
def read_true_values(true_covariance_and_multiplier='',
                      true_variance_correction=None,
                      subnodes_wo_outgroup=[]):
    Rcovariance,multiplier,vc=file_to_emp_cov(true_covariance_and_multiplier, sort_nodes_alphabetically=True, vc=true_variance_correction, return_only_covariance=False, subnodes=sorted(subnodes_wo_outgroup))
    return (Rcovariance,multiplier)

def mode(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]
    
    
          
