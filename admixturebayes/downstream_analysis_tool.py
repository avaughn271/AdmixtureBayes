from tree_statistics import (identifier_to_tree_clean, generate_predefined_list_string, unique_identifier_and_branch_lengths)
from Rtree_to_covariance_matrix import get_populations
import pandas as pd
from copy import deepcopy
from tree_to_data import file_to_emp_cov

from tree_statistics import admixture_sorted_unique_identifier

from subgraphing import get_subtree

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
        d_dic={colname:r[k] for k, colname in enumerate(cols)}
        d_dic.update(constant_kwargs)

        for row_summarize_function in row_summarize_functions:
            add_dic, skip=row_summarize_function(**d_dic)
            d_dic.update(add_dic)
        all_results.append(thinned_d_dic(d_dic))
    return all_results

class make_Rtree(object):
    
    def __init__(self, nodes_to_be_sorted, subnodes=[]):
        self.nodes=sorted(nodes_to_be_sorted)
        self.subnodes=subnodes
        
    def __call__(self, tree):
        first_level=tree.split('-')[0]
        no_pops=len(first_level.split('.'))
        if no_pops==len(self.nodes):
            Rtree=identifier_to_tree_clean(tree, leaves=generate_predefined_list_string(deepcopy(self.nodes)))
        else:
            assert False, 'Either the outgroup name was not specified or something is seriously wrong because' \
                          ' the number of nodes did not match the size of the trees'

        if self.subnodes:
            Rtree=get_subtree(Rtree, self.subnodes)
        return {'Rtree':Rtree}, False
    
class make_full_tree(object):
    
    def __init__(self, add_multiplier=1, subnodes=[]):
        self.add_multiplier=add_multiplier
        self.subnodes=subnodes
        
    def __call__(self, Rtree=None, add=None, **kwargs):
        if Rtree is None:
            assert 'sfull_tree' in kwargs, 'sfull_tree not specified'
            nodes=sorted(kwargs['full_nodes'])
            sfull_tree=kwargs['sfull_tree']
            full_tree=identifier_to_tree_clean(sfull_tree, leaves=generate_predefined_list_string(deepcopy(nodes)))
            if self.subnodes:
                full_tree=get_subtree(full_tree, self.subnodes)
            return {'full_tree':full_tree}, False
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
    
    def __init__(self, burn_in_fraction=None, total=None):
        self.burn_in_fraction=burn_in_fraction
        self.total=total
        
    def __call__(self, df):
        n=len(df)
        print('Dataframe read with ' + str(n) + ' samples.')
        if self.burn_in_fraction is not None:
            df=df[int(n*self.burn_in_fraction):]
        print('Burn-in of ' + str(n-len(df)) + ' samples removed. There are now ' + str(len(df)) + ' samples.')
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
