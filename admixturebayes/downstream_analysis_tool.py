from reduce_covariance import reduce_covariance
from Rtree_operations import add_outgroup
from tree_statistics import (identifier_to_tree_clean, generate_predefined_list_string, unique_identifier_and_branch_lengths)
from Rtree_to_covariance_matrix import make_covariance, get_populations
import numpy as np
from construct_covariance_choices import read_one_line
import pandas as pd
from copy import deepcopy
from tree_to_data import file_to_emp_cov

from tree_statistics import admixture_sorted_unique_identifier
from Rtree_operations import get_leaf_keys

from collections import Counter
from subgraphing import get_subtree
from find_true_trees import get_unique_plottable_tree

def get_list_of_turned_topologies(trees, true_tree):
    nodes=get_leaf_keys(true_tree)
    return [admixture_sorted_unique_identifier(tree, nodes) for tree in trees], admixture_sorted_unique_identifier(true_tree, nodes)

def always_true(**args):
    return True

def identity(x):
    return x

def iterate_over_output_file(outfile, 
                             cols=[], 
                             pre_thin_data_set_function=identity, 
                             while_thin_data_set_function=always_true,
                             row_summarize_functions=[],
                             thinned_d_dic=identity,
                             full_summarize_functions=[],
                             **constant_kwargs):
    
    df= pd.read_csv(outfile, usecols=cols, dtype={'no_admixes':object})
    df = df[cols]
    df= pre_thin_data_set_function(df)
    full_summs=[full_summarize_function(df) for full_summarize_function in full_summarize_functions]
    
    all_results=[]
    
    for n,(i,r) in enumerate(df.iterrows()):
        cont=False
        d_dic={colname:r[k] for k, colname in enumerate(cols)}
        d_dic.update(constant_kwargs)
        if not while_thin_data_set_function(**d_dic):
            print('breaking from while_thin')
            continue

        for row_summarize_function in row_summarize_functions:
            #removedprin row_summarize_function
            #removedprin d_dic
            add_dic, skip=row_summarize_function(**d_dic)
            if skip:
                print('breaks because of the value', add_dic)
                cont=True
                break
            d_dic.update(add_dic)
        if cont:
            print('breaking because of cont')
            continue
        all_results.append(thinned_d_dic(d_dic))
    return all_results, full_summs

class make_Rtree(object):
    
    def __init__(self, nodes_to_be_sorted, remove_sadtrees=False, subnodes=[], outgroup_name=''):
        self.nodes=sorted(nodes_to_be_sorted)
        self.remove_sadtrees=remove_sadtrees
        self.subnodes=subnodes
        self.outgroup_name=outgroup_name
        
    def __call__(self, tree, **not_needed):
        #removedprin tree
        #removedprin not_needed
        #removedprin tree
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

class make_Rcovariance(object):
    
    def __init__(self, nodes, add_multiplier=1):
        self.nodes=nodestr
        self.add_multiplier=add_multiplier
        
    def __call__(self, Rtree=None, add=None, **kwargs):
        if Rtree is None:
            full_tree=kwargs['full_tree']
            outgroup_name=list(set(get_leaf_keys(full_tree))-set(self.nodes))[0]
            cov=make_covariance(full_tree, node_keys=[outgroup_name]+self.nodes)
            Rcov=reduce_covariance(cov, 0)
            return {'Rcov':Rcov}, False
        #removedprin self.nodes
        Rcov=make_covariance(Rtree, node_keys=self.nodes)+float(add)*self.add_multiplier
        
        return {'Rcov':Rcov}, False
    
class cov_truecov(object):
    
    def __init__(self, true_covariance):
        self.true_covariance=true_covariance
        
    def __call__(self, Rcov, **not_needed):
        dist=np.linalg.norm(Rcov-self.true_covariance)
        return {'cov_dist':dist}, False
    
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


class topology_without_outgroup(object):

    def __init__(self, outgroup_to_remove):
        self.outgroup_to_remove = outgroup_to_remove

    def __call__(self, Rtree=None, **kwargs):
        if Rtree is None:
            Rtree = kwargs['full_tree']
        print("SEVERE ERROR HAS OCCURRED")
        top = admixture_sorted_unique_identifier(tmp_tree, leaf_order=self.nodes, not_opposite=True)
        return {'topology': top}, False


    
class subgraph(object):
    
    def __init__(self, subgraph_keys, identifier='', max_num=10, total_probability=1.0, prefix='',**not_needed):
        self.subgraph_keys=subgraph_keys
        self.identifier=identifier
        self.max_num=max_num
        self.total_probability=total_probability
        self.prefix=prefix
        
    def __call__(self, full_tree=None, **kwargs):
        if full_tree is None:
            full_tree=kwargs['Rtree']
        full_tree=deepcopy(full_tree)
        subgraph=get_subtree(Rtree, self.subgraph_keys)
        sub_stree=get_unique_plottable_tree(subgraph)
        return {'subgraph_'+self.identifier:sub_stree}, False 

class subsets(object):
    
    def __init__(self, subgraph_keys, identifier='', prefix='', max_num=10, **not_needed):
        self.subgraph_keys=subgraph_keys
        self.identifier=identifier
        self.prefix=prefix
        self.max_num=max_num
        
    def __call__(self, pops, **kwargs):
        subsets=get_subpops(pops, self.subgraph_keys)
        return {'subsets_'+self.identifier:subsets}, False   
         
    
class topology_identity(object):
    
    def __init__(self, true_Rtree, nodes):
        self.full_scaled_turned_topology=admixture_sorted_unique_identifier(true_Rtree, leaf_order=nodes, not_opposite=True)
        
    def __call__(self, topology, **not_needed):
        #removedprin topology, self.full_scaled_turned_topology
        ident_top=(topology==self.full_scaled_turned_topology)
        return {'top_identity':ident_top}, False
    
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
    
class compare_pops(object):
    def __init__(self, true_tree, min_w=0.0, keys_to_include=None):
        self.true_pops=set(get_populations(true_tree, min_w, keys_to_include))
        
    def __call__(self, pops, **not_needed):
        diffs=set(pops).symmetric_difference(self.true_pops)
        return {'set_differences':len(diffs)}, False

class set_identity(object):
    def __init__(self):
        pass
    
    def __call__(self, set_differences,**not_needed):
        ident=(set_differences==0)
        return {'set_identity':ident}

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
        #print('after filtering operations there are now', len(df),'samples')
        if self.total is not None:
            n=len(df)
            #stepsize=max(n//self.total,1) #ANDREWDEBUG
            stepsize = self.total
            df=df[::stepsize]
            print('Thinning every ' + str(stepsize) + ' samples complete. There are now ' + str(len(df)) + ' samples')
        return df
    
def read_true_values(true_scaled_tree='', 
                      true_tree='',
                      true_add='',
                      true_covariance_reduced='',
                      true_covariance_and_multiplier='',
                      true_no_admix='',
                      true_m_scale='',
                      true_variance_correction=None,
                      true_df='',
                      subnodes_wo_outgroup=[],
                      subnodes_with_outgroup=[]):
    scaled_tree,tree,add,covariance_reduced,(Rcovariance,multiplier), no_admix,m_scale,vc,df=None,None,None,None,(None,None),None,None,None,None
    if true_covariance_and_multiplier:
        Rcovariance,multiplier,vc=file_to_emp_cov(true_covariance_and_multiplier, sort_nodes_alphabetically=True, vc=true_variance_correction, return_only_covariance=False, subnodes=sorted(subnodes_wo_outgroup))
    if true_m_scale:
        m_scale=float(read_one_line(true_m_scale))
    if true_df:
        df=float(read_one_line(true_df))
    return scaled_tree,tree,add,covariance_reduced,(Rcovariance,multiplier), no_admix,m_scale,vc,df

def mode(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]
    
    
          
