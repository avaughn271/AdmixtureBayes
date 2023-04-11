from argparse import ArgumentParser, SUPPRESS
import pandas as pd
import os
from copy import deepcopy
from collections import Counter
from graphviz import Digraph
from numpy.random import choice

from itertools import chain
from numpy import random
import numpy
from math import exp, log

import Rtree_to_covariance_matrix
import construct_starting_trees_choices
from construct_covariance_choices import get_covariance, estimate_degrees_of_freedom_scaled_fast
from posterior import posterior_class

from tree_statistics import (identifier_to_tree_clean, get_admixture_proportion_string, generate_predefined_list_string, admixture_sorted_unique_identifier, unique_identifier_and_branch_lengths,topological_identifier_to_tree_clean, identifier_to_tree)

from Rtree_operations import node_is_admixture, rename_key, get_admixture_proportion_from_key, get_all_admixture_origins, to_networkx_format,change_admixture, get_categories, get_number_of_ghost_populations, get_number_of_admixes, scale_tree_copy

from mcmc_proposals import addadmix_class, deladmix_class, sliding_regraft_class, rescale_class, rescale_admixtures_class, rescale_constrained_class, rescale_add_class

from multiprocessing import Queue, Process

class basic_chain_class_as_process(object):
    
    def __init__(self, basic_chain_class):
        self.chain=basic_chain_class
        self.process= Process(target=self.chain)
        self.process.start()
    
class basic_chain_class(object):
    
    def __init__(self, summaries, posterior_function, proposal, resxeed):
        self.summaries=summaries
        self.posterior_function=posterior_function
        self.proposal=proposal
        self.task_queue= Queue()
        self.response_queue = Queue()

    def __call__(self):
        while True:
            input = self.task_queue.get()
            self.response_queue.put(self.run_chain(input))
            
    def run_chain(self, p):
        start_tree, post, N, sample_verbose_scheme, overall_thinning, i_start_from, temperature, proposal_update, multiplier = p
        return basic_chain(start_tree,  self.summaries,  self.posterior_function, self.proposal,  post,  N,  sample_verbose_scheme, i_start_from, temperature, proposal_update, multiplier)
        
class basic_chain_pool(object):
    
    def __init__(self, summaries, posterior_function, proposals, seeds, posterior_function_list=[]): #SxEEDDEBUG
        self.group=[basic_chain_class_as_process(basic_chain_class(summaries, posterior_function, proposals[0],None))]

    def order_calculation(self, list_of_lists_of_arguments):
        for i in list_of_lists_of_arguments:
            listargs = i
        (self.group[0]).chain.task_queue.put(listargs)
        a = [(self.group[0]).chain.response_queue.get()]
        return a

def one_jump(x, post, temperature, posterior_function, proposal, pks={}):
    
    newx,g1,g2,Jh,j1,j2=proposal(x,pks)
    post_new=posterior_function(newx)

    likelihood_old, prior_old = post[:2]
    likelihood_new, prior_new = post_new[:2]
    
    TargetDensitynew = likelihood_new + prior_new
    TargetOld = prior_old + likelihood_old
    if TargetDensitynew == -float('inf'):
        return(x, post, 0, 1)
    logmhr = TargetDensitynew/temperature - TargetOld/temperature # maybe only change this and no MC^3 switches and no adaption.
    if logmhr>100:
        mhr=float('inf')
    else:
        mhr=exp(logmhr)
            
    u=random.random()
    if u<mhr:
        return newx,post_new,1,TargetOld>TargetDensitynew
    return x,post,0,TargetOld>TargetDensitynew

def basic_chain(start_x, summaries, posterior_function, proposal, post=None, N=10000, 
                sample_verbose_scheme=None, i_start_from=0, 
                temperature=1.0, proposal_update=None, multiplier=None):
    if proposal_update is not None:
        proposal.node_naming.n=proposal_update['n']
    
    x=start_x
    if post is None:
        post=posterior_function(x)
    
    iteration_summary=[]
    NumberOfAcceptedUphillMoves = 0
    NumberOfUphillProposals = 0
    NumberOfAcceptedTotal = 0
    NumberOfProposals = 0
    for i in range(i_start_from,i_start_from+N):
        proposal_knowledge_scraper={}
        new_x,new_post,acceptedd,wasuphill=one_jump(x, post, temperature, posterior_function, proposal, proposal_knowledge_scraper)
        NumberOfProposals = NumberOfProposals + 1
        NumberOfAcceptedTotal = NumberOfAcceptedTotal + acceptedd
        NumberOfUphillProposals = NumberOfUphillProposals + wasuphill
        NumberOfAcceptedUphillMoves = NumberOfAcceptedUphillMoves + wasuphill * acceptedd
        if i%40==0:
            iteration_summary.append(_calc_and_print_summaries(sample_verbose_scheme,
                                                               summaries,
                                                               tree=scale_tree_copy(new_x[0], 1.0/multiplier),
                                                               add=new_x[1],
                                                               posterior=new_post,
                                                               old_post=post,
                                                               old_tree=scale_tree_copy(x[0],1.0/multiplier),
                                                               iteration_number=i,**proposal_knowledge_scraper))
        x=new_x
        post=new_post
    if temperature < 1:
        print("Accepted " + str(NumberOfAcceptedUphillMoves) + " out of " + str(NumberOfUphillProposals) + " uphill moves (" +str(round(NumberOfAcceptedUphillMoves/NumberOfUphillProposals,4)*100)[0:5] + "%) at temp=" + '{:.3e}'.format(temperature) + ".", end =' ',flush=True)
    else:
        print("Accepted " + str(NumberOfAcceptedUphillMoves) + " out of " + str(NumberOfUphillProposals) + " uphill moves (" +str(round(NumberOfAcceptedUphillMoves/NumberOfUphillProposals,4)*100)[0:5] + "%) at temp=" + str(round(temperature, 3)) + ".", end =' ', flush = True)
    return x, post, list(zip(*iteration_summary)),proposal.get_exportable_state()

def _calc_and_print_summaries(sample_verbose_scheme,summaries,**kwargs):
    iteration=kwargs['iteration_number']
    res=[iteration]
    for s in summaries:
        save_num,print_num=sample_verbose_scheme.get(s.name, (0,0))
        save_bool = (save_num!=0) and (iteration % save_num==0)
        if save_bool:
            val=s(**kwargs)
            res.append(val)
        else:
            res.append(None)
    return res

def MCMCMC(coolerMultiple, starting_trees,    posterior_function, summaries, temperature_scheme,  printing_schemes, 
           iteration_scheme,  proposal_scheme, n_arg, verboseee, numpy_seeds=None, multiplier= None, result_file=None):
    '''
    this function runs a MC3 using the basic_chain_unpacker. Let no_chains=number of chains. The inputs are
        starting_trees: a list of one or more trees that the chains should be started with
        proposal_scheme: a list of instances of classes that handles proposals and updates of said proposals. It has the basic functions
                                    - prop(x,pks): proposes the next tree(and returns proposal densities and statistics in pks)
                                    - adapt(mhr): updates the proposals based on the mhr ratio, mhr                                    
    '''
    df_result=None
    xs = starting_trees

    pool = basic_chain_pool(summaries, posterior_function, proposal_scheme, numpy_seeds)
    posteriors = [posterior_function(x) for x in xs]

    proposal_updates=[proposal.get_exportable_state() for proposal in proposal_scheme]
    
    cum_iterations=0
    for no_iterations in iteration_scheme:
        numberremaining = int(len(iteration_scheme) - cum_iterations/no_iterations)
        if cum_iterations != 0: print(str(numberremaining) + " temps left.")
        #letting each chain run for no_iterations:
        iteration_object=_pack_everything(xs, posteriors, temperature_scheme, printing_schemes, no_iterations, cum_iterations, proposal_updates, multiplier)
        new_state = pool.order_calculation(iteration_object)
        xs, posteriors, df_add,proposal_updates = _unpack_everything(new_state, summaries)
        df_result=_update_results(df_result, df_add)
        if result_file is not None:
            if cum_iterations==0:
                start_data_frame(df_result, result_file)
            elif df_result.shape[0]>1000:
                add_to_data_frame(df_result, result_file)
                df_result=df_result[0:0]
        #making the mc3 flips and updating:
        temperature_scheme[0] = temperature_scheme[0] * coolerMultiple
        cum_iterations+=no_iterations
    print("Done!")
    for chain in pool.group:
            chain.process.terminate()

def start_data_frame(df, result_file):
    df=df.loc[df.layer==0,:]
    df.to_csv(result_file, header=True)

def add_to_data_frame(df_add, result_file):
    df_add=df_add.loc[df_add.layer==0,:]
    with open(result_file, 'a') as f:
        df_add.to_csv(f, header=False)
        
def _update_results(df_result, df_add):
    if df_result is None:
        df_result = df_add
    else:
        df_result = pd.concat([df_result, df_add])
    return df_result

def _pack_everything(xs, posteriors, temperature_scheme,printing_schemes,no_iterations,cum_iterations, proposal_updates=None, multiplier=None):
    return ([x, posterior,no_iterations, printing_scheme, 40, cum_iterations, temperature_scheme[i], proposal_update, multiplier] for i,(x,posterior,printing_scheme,proposal_update) in enumerate(zip(xs,posteriors,printing_schemes,proposal_updates)))

def _unpack_everything(new_state, summaries):
    xs,posteriors, summs, proposal_updates = list(zip(*new_state))
    list_of_smaller_data_frames=[]
    for summ_data, i in zip(summs, list(range(1))):
        iter_chain=chain((('iteration', summ_data[0]),),
                         ((summ_object.name,summ_col) for summ_object,summ_col in zip(summaries,summ_data[1:])))
        df=pd.DataFrame.from_dict(dict(iter_chain))
        df['layer']=i
        list_of_smaller_data_frames.append(df)
    df=pd.concat(list_of_smaller_data_frames)
    return list(xs), list(posteriors), df, list(proposal_updates)

def rename_root_process(tree, old_name):
    for _,node in list(tree.items()):
        if node[0]==old_name:
            node[0]='r'
        if (node[1] is not None and node[1]==old_name):
            node[1]='r'
    return tree

def update_parent_and_branch_length_process(tree, child_key, child_branch, new_parent, new_branch_length):
    tree[child_key][child_branch]=new_parent
    tree[child_key][child_branch+3]=new_branch_length
    return tree

def insert_children_in_tree_process(tree):
    children={key:[] for key in tree}
    for key in tree:
        parents = get_real_parents_process(tree[key])
        for parent in parents:
            if parent!='r':
                children[parent].append(key)
    for key in tree:
        tree[key]=_update_parents_process(tree[key], children[key])
    return tree

def _update_parents_process(node, new_parents):
    if len(new_parents)==1:
        res=node[:5]+[new_parents[0],None]
        return res
    if len(new_parents)==2:
        res=node[:5]+new_parents
        return res
    if len(new_parents)==0:
        res=node[:5]+[None]*2
        return res

def get_real_parents_process(node):
    ps=node[:2]
    return [p for p in ps if p is not None]

class generate_numbered_nodes_process(object):

    def __init__(self, prefix='n', admixture_prefix='a'):
        self.node_count=0
        self.admixture_count=0
        self.prefix='n'
        self.admixture_prefix='a'

    def __call__(self, admixture=False):
        if admixture:
            self.admixture_count+=1
            return self.admixture_prefix+str(self.admixture_count)
        self.node_count+=1
        return self.prefix+str(self.node_count)

class generate_predefined_list_string_process(object):

    def __init__(self, listi):
        self.listi=listi

    def __call__(self):
        return self.listi.pop(0)

def identifier_to_tree_process(identifier, leaves=None, inner_nodes=None, branch_lengths=None, admixture_proportions=None):
    '''
    Transforms an identifier of the form qwert-uio-asdfg-jk into a dictionary tree using the generators of leaves, inner_nodes, branch_lengths and admixture_proportions.
    '''
    levels=identifier.split('-')
    n_leaves=len(levels[0].split('.'))
    if leaves is None:
        leaf_values=sorted(['s'+str(n+1) for n in range(n_leaves)])
    else:
        leaf_values=[leaves() for _ in range(n_leaves)]
    tree={leaf:[None]*5 for leaf in leaf_values}
    trace_lineages=[(leaf,0) for leaf in leaf_values]
    if inner_nodes is None:
        inner_nodes=generate_numbered_nodes_process('n')
    if branch_lengths is None:
        def f():
            return 1.0
        branch_lengths= f
    if admixture_proportions is None:
        def g():
            return 0.4
        admixture_proportions=g
    for level in levels:
        identifier_lineages=level.split('.')
        parent_index={}
        indexes_to_be_removed=[]
        for n,identifier_lineage in enumerate(identifier_lineages):
            if identifier_lineage=='c':
                ##there is a coalecence for the n'th lineage, and it should be replaced by a new lineage
                new_key=inner_nodes()
                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length_process(tree, old_key, old_branch, new_key, new_branch_length)
                tree[new_key]=[None]*5
                parent_index[n]=new_key
                trace_lineages[n]=(new_key,0)
            elif identifier_lineage=='w':
                pass
            elif identifier_lineage=='a':
                new_key=inner_nodes(admixture=True)
                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length_process(tree, old_key, old_branch, new_key, new_branch_length)
                new_admixture_proportion=admixture_proportions()
                tree[new_key]=[None,None,new_admixture_proportion,None,None]
                trace_lineages[n]=(new_key,0)
                trace_lineages.append((new_key,1))
            else:
                try:
                    new_key=parent_index[int(identifier_lineage)]
                except KeyError as e:
                    print(e)

                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length_process(tree, old_key, old_branch, new_key, new_branch_length)
                indexes_to_be_removed.append(n)

        ##remove lineages
        trace_lineages=[trace_lineage for n,trace_lineage in enumerate(trace_lineages) if n not in indexes_to_be_removed]
    root_key=new_key
    del tree[root_key]
    tree=rename_root_process(tree, new_key)

    return insert_children_in_tree_process(tree)

def get_numeric_process(string_tree):
    branch_lengths_string, admixture_proportion_string= string_tree.split(';')[1:]
    branch_lengths=list(map(float,branch_lengths_string.split('-')))
    if admixture_proportion_string:
        admixture_proportions=list(map(float, admixture_proportion_string.split('-')))
    else:
        admixture_proportions=[]
    return branch_lengths, admixture_proportions

def branch_and_proportion_quantiles_process(list_of_string_trees):
    ''' extracts the branch lengths and admixture proportions and returns them as tuples of four. the list should not be empty'''
    branches=[]
    admixtures=[]
    for string_tree in list_of_string_trees:
        b,a=get_numeric_process(string_tree)
        branches.append(b)
        admixtures.append(a)
    bdat=pd.DataFrame.from_records(branches)
    adat=pd.DataFrame.from_records(admixtures)
    bresults=[]
    aresults=[]
    for n,(mean) in enumerate(bdat.mean()):
        bresults.append(('c'+str(n+1), mean))
    if len(admixtures[0])>0:
        for n,(mean) in enumerate(adat.mean()):
            aresults.append(('ax'+str(n+1), mean))
    return bresults, aresults

def mainnnoutput():
    df = pd.read_csv('mcmc_samples_annealing.csv')
    totalnummmofrows = df.shape[0]
    posteriorr = df['posterior'].values.tolist()
    modetree = ((df[["tree"]]).iloc[totalnummmofrows-1])[0]
    print(df.shape, "total number of rows", posteriorr[totalnummmofrows-1] )
    print(posteriorr.index(max(posteriorr)), "index of posterior mode",  max(posteriorr))
    ff = open("MAPadd.txt", "w")
    ff.write(str((df[["add"]]).iloc[totalnummmofrows-1][0] ) +  "\n")
    ff.close
    leaves = sorted(list(set([leaf for node in (df['descendant_sets'].tolist()[0].split('-')) for leaf in node.split('.')])))

    chosen_topology = modetree
    temptopology = chosen_topology.split(';')[0]
    for i in leaves[::-1]:
        chosen_topology = i + "=" + chosen_topology
    branches_intervals, admixture_proportion_intervals=branch_and_proportion_quantiles_process([chosen_topology])

    branch_names=[branches_interval[0] for branches_interval in branches_intervals]
    admixture_names=[ad[0] for ad in admixture_proportion_intervals]

    tree = identifier_to_tree_process(temptopology, leaves=generate_predefined_list_string_process(deepcopy(leaves)), branch_lengths=generate_predefined_list_string_process(deepcopy(branch_names)), admixture_proportions=generate_predefined_list_string_process(deepcopy(admixture_names)))

    f = open("MAPtree.txt", "w")

    for node in tree:
        if (tree[node][2] is not None):
            for adds in admixture_proportion_intervals:
                for branchh in branches_intervals:
                    if tree[node][2] == adds[0]:
                        if tree[node][3] == branchh[0]:
                            f.write( node + " " + tree[node][0] + " "+ str(branchh[1]) + " "+ str(adds[1])+"\n")
                        if tree[node][4] == branchh[0]:
                            f.write( node + " " + tree[node][1] + " "+ str(branchh[1]) + " "+ str(1-adds[1])+"\n")
        else:
            for branchh in branches_intervals:
                if tree[node][3] == branchh[0]:
                            f.write( node + " " + tree[node][0] + " "+ str(branchh[1]) + " 1.00"+"\n")
    f.close



class new_node_naming_policy(object):
    
    def __init__(self, n=0):
        self.n=0
        
    def next_nodes(self, no_nodes):
        if no_nodes==2:
            self.n+=1
            return ['x'+str(self.n)+a for a in ['a','b']]
        elif no_nodes==1:
            self.n+=1
            return 'x'+str(self.n)
        else:
            return ''

def initialize_proposals(proposals):
    all_props=[addadmix_class, deladmix_class, rescale_class, sliding_regraft_class, rescale_add_class, rescale_constrained_class,  rescale_admixtures_class]
    all_props_dic={cl.proposal_name:cl for cl in all_props}
    res=[]
    for proposal in proposals:
        res.append(all_props_dic[proposal]())
    return res
    
def draw_proposal(props, k, proportions):
    
    legal_indices=[i for i,prop in enumerate(props) if prop.require_admixture<=k]
    normaliser=sum([proportion for n,proportion in enumerate(proportions) if n in legal_indices])
    new_proportions=[float(proportion)/normaliser for n,proportion in enumerate(proportions) if n in legal_indices]
    
    chosen_index_i= choice(len(legal_indices), 1, p=new_proportions)[0]
    chosen_index=legal_indices[chosen_index_i]
    
    effect_of_chosen_index=props[chosen_index].admixture_change
    if effect_of_chosen_index!=0:
        legal_indices2=[i for i,prop in enumerate(props) if prop.require_admixture <= k+effect_of_chosen_index]    
        normaliser2=sum([proportion for n,proportion in enumerate(proportions) if n in legal_indices2])
        new_proportions2=[float(proportion)/normaliser2 for n,proportion in enumerate(proportions) if n in legal_indices2]
        reverse_type= props[chosen_index].reverse
        reverse_index= next((index for index, prop in enumerate(props) if prop.proposal_name==reverse_type))
        reverse_index_i= next((index_i for index_i, index in enumerate(legal_indices2) if index==reverse_index))
        return chosen_index, new_proportions[chosen_index_i], new_proportions2[reverse_index_i]
    else:
        return chosen_index, 1.96,1.96 #it is not really 1.96 and 1.96 but only the ratio between them matters and I like 1.96
    
def get_args2(names):
    args=[]
    if names:
        args.append(names)
    return args    

class simple_adaptive_proposal(object):
    
    def __init__(self, proposals, proportions):
        self.props=initialize_proposals(proposals)
        self.proportions=proportions
        self.node_naming=new_node_naming_policy()
        self.recently_called_index=None
        
    def __call__(self, x, pks={}):
        tree,add=x
        k=get_number_of_admixes(tree)
        index, jforward, jbackward = draw_proposal(self.props, k, self.proportions)
        
        names=self.node_naming.next_nodes(self.props[index].new_nodes)
        self.recently_called_index=index
        proposal_input= self.props[index].input
        args=get_args2(names)
        
        if proposal_input=='add':
            new_add, forward, backward = self.props[index](add, *args, pks=pks)
            return (tree, new_add), forward, backward, 1.0, jforward, jbackward
        if proposal_input=='tree':
            new_tree, forward, backward = self.props[index](tree, *args, pks=pks)
            return (new_tree, add), forward, backward, 1.0, jforward, jbackward
        else:
            new_x, forward, backward = self.props[index](x, *args, pks=pks)
            return new_x, forward, backward, 1.0, jforward, jbackward
        
    def get_exportable_state(self):
        information={}
        information['n']=self.node_naming.n
        return information
        



class make_read_name(object):
    
    def __init__(self):
        self.count=0
        
    def __call__(self, admix=False):
        self.count+=1
        if admix:
            return 'a'+str(self.count)
        else:
            return 'n'+str(self.count)

def node_structure_to_networkx(node_structure, node_freq_dictionary=None):
    edges=[]
    admixture_nodes=[]
    pure_leaves=[]
    admixture_leaves=[]
    root=[]
    coalescence_nodes=[]
    namer=make_read_name()
    string_names={}
    for key,node in list(node_structure.items()):
        if key==frozenset(['']):
            continue
        plotting_type=node.get_plotting_type()
        if plotting_type=='coalescence':
            string_names[key]=namer()
            coalescence_nodes.append(string_names[key])
        elif plotting_type=='admixture':
            string_names[key]=namer(admix=True)
            admixture_nodes.append(string_names[key])
        elif plotting_type=='leaf':
            string_names[key]=list(key)[0]
            pure_leaves.append(string_names[key])
        elif plotting_type=='admixture_leaf':
            string_names[key]=list(key)[0]
            admixture_leaves.append(string_names[key])
        elif plotting_type=='root':
            string_names[key]='r'
            root.append(string_names[key])
        else:
            assert False, 'unknown plotting type'+str(plotting_type)
        
    for key, node in list(node_structure.items()):
        parents=node.parents
        for parent in parents:
            edges.append((string_names[parent.name],string_names[key]))
    return pure_leaves, admixture_leaves, coalescence_nodes, admixture_nodes, root, edges

class Node():
    
    def __init__(self, name):
        self.name=name
        self.parents=[]
        self.children=[]
        
    def add_child(self, child_node):
        self.children.append(child_node)
        
    def add_parent(self, parent_node):
        self.parents.append(parent_node)
    
    def get_plotting_type(self):
        if self.children:
            if len(self.parents)==0:
                return 'root'
            if len(self.parents)==1:
                return 'coalescence'
            if len(self.parents)>1:
                return 'admixture'
        else:
            if len(self.parents)>1:
                return 'admixture_leaf'
            return 'leaf'
    
    def has_parent(self):
        return len(self.parents)>0
    
    def is_admixture(self):
        return len(self.children)==1 and len(self.parents)==2

def plot_node_structure_as_directed_graph(node_structure, drawing_name='tmp_.png', popup=True, node_dic=None, verbose=True):
    pure_leaves, admix_leaves, pure_coalescence, admix_coalescence, root, edges= node_structure_to_networkx(node_structure, node_dic)
    filename, image_format= drawing_name.split('.')
    image_format = "pdf" #ANDREWDEBUG
    G=Digraph('G', filename=filename)
    leaves_graph=Digraph('l')
    leaves_graph.node_attr.update(style='filled', color='cadetblue1')
    for leaf in pure_leaves:
        leaves_graph.node(leaf)
        
    aleaves_graph=Digraph('l')
    aleaves_graph.node_attr.update(style='filled', color='slateblue1')
    for leaf in admix_leaves:
        aleaves_graph.node(leaf)
        
    admixture_graph=Digraph()
    admixture_graph.node_attr.update(style='filled', fillcolor='coral1', shape='box')
    for adm_n in admix_coalescence:
        admixture_graph.node(adm_n, fillcolor='#%02x%02x%02x' % (0, 153, 0))
        
    coalescence_graph=Digraph()
    coalescence_graph.node_attr.update(style='filled', fillcolor='greenyellow')
    for cn in pure_coalescence:
        coalescence_graph.node(cn, fillcolor='#%02x%02x%02x' % (0, 153, 0))
        
    G.subgraph(leaves_graph)
    G.subgraph(aleaves_graph)
    G.subgraph(admixture_graph)
    G.subgraph(coalescence_graph)
    G.node('r', shape='egg', color='black', style='filled', fontcolor='white')
    G.edges(edges)
    G.format = image_format
    G.render(view=popup)
    if verbose:
        if os.path.exists(filename):
            os.remove(filename)
        else:
            print("The file does not exist")
        
def plot_as_directed_graph(tree, drawing_name='tmp.png', popup=True, plot_edge_lengths=False, verbose=True ,labeldict = {}, admixturegroup = {}, valtostring=[], ISLABELS = False, outgrooupp = "", disttt = 0, dashess = {}):

    leaves, admixture_nodes, coalescence_nodes, root, edges, edge_lengths= to_networkx_format(tree)
    filename, image_format= drawing_name.split('.')
    image_format = "pdf" #ANDREWDEBUG
    G=Digraph('G', filename=filename)
    
    leaves_graph=Digraph('l')
    leaves_graph.node_attr.update(style='filled', color='cadetblue1')
    for leaf in leaves:
        leaves_graph.node(leaf)
        
    admixture_graph=Digraph()
    admixture_graph.node_attr.update(style='filled', fillcolor='coral1', shape='box')
    for adm_n in admixture_nodes:
        if len(valtostring) == 0:
            admixture_graph.node(adm_n)
        else:
            for jjjj in valtostring:
                if jjjj[0] == adm_n:
                    admixture_graph.node(adm_n, label = str(round(   min(jjjj[1]  ,1-jjjj[1]  ), 3)))
        
    coalescence_graph=Digraph()
    coalescence_graph.node_attr.update(style='filled', fillcolor='greenyellow')
    for cn in coalescence_nodes:
        coalescence_graph.node(cn)
        
    G.subgraph(leaves_graph)
    G.subgraph(admixture_graph)
    G.subgraph(coalescence_graph)
    G.node('r', shape='egg', color='black', style='filled', fontcolor='white')
    if plot_edge_lengths:
        for (to_node, from_node),label in zip(edges, edge_lengths):
            if to_node[0] == "a" and to_node[1] == "x" and to_node[2] in ["0","1","2","3","4","5","6","7","8","9"]:
                tonodecheck = "a" + to_node[2:]
            else:
                tonodecheck = to_node
            if from_node[0] == "a" and from_node[1] == "x" and from_node[2] in ["0","1","2","3","4","5","6","7","8","9"]:
                from_nodecheck = "a" + from_node[2:]
            else:
                from_nodecheck = from_node
            if from_nodecheck in dashess and dashess[from_nodecheck][0] == tonodecheck and dashess[from_nodecheck][1] != dashess[from_nodecheck][0]:
                if dashess[from_nodecheck][2] < 0.0001:
                    G.edge(to_node, from_node, label=  f'{dashess[from_nodecheck][2]:.2e}', style = "dashed")
                else:
                    G.edge(to_node, from_node, label=  str(round(dashess[from_nodecheck][2], 4)), style = "dashed")
            else:
                if dashess[from_nodecheck][3] < 0.0001:
                    G.edge(to_node, from_node, label=  f'{dashess[from_nodecheck][3]:.2e}')
                else:
                    G.edge(to_node, from_node, label=  str(round(dashess[from_nodecheck][3], 4)))
    else:
        G.edges(edges)
    if ISLABELS:
        G.node(outgrooupp, shape='diamond', color='gray', style='filled')
        if disttt < 0.0001:
            G.edge( outgrooupp, 'r', dir="none", label = f'{disttt:.2e}')
        else:
            G.edge( outgrooupp, 'r', dir="none", label = str(round(disttt, 4)))
    G.format = image_format
    G.render(view=popup)
    if verbose:
        if drawing_name[0:10] == "topology_l":
            strlength = len(drawing_name)
            numm = (drawing_name[16:(strlength-4)]) #ANDREWDEBUG
        if os.path.exists(filename):
            os.remove(filename)
        else:
            print("The file does not exist")
            
def get_numeric(string_tree):
    branch_lengths_string, admixture_proportion_string= string_tree.split(';')[1:]
    branch_lengths=list(map(float,branch_lengths_string.split('-')))
    if admixture_proportion_string:
        admixture_proportions=list(map(float, admixture_proportion_string.split('-')))
    else:
        admixture_proportions=[]
    return branch_lengths, admixture_proportions

def branch_and_proportion_quantiles(list_of_string_trees):
    ''' extracts the branch lengths and admixture proportions and returns them as tuples of four. the list should not be empty'''
    branches=[]
    admixtures=[]
    for string_tree in list_of_string_trees:
        b,a=get_numeric(string_tree)
        branches.append(b)
        admixtures.append(a)
    bdat=pd.DataFrame.from_records(branches)
    adat=pd.DataFrame.from_records(admixtures)
    bresults=[]
    aresults=[]
    for n,mean in enumerate(zip(bdat.mean())):
        bresults.append(('c'+str(n+1),mean[0]))
    if len(admixtures[0])>0:
        for n,mean in enumerate(zip(adat.mean())):
            aresults.append(('ax'+str(n+1), mean[0]))
    return bresults, aresults

def main_plot(plottype, outtttname, outputprefix, disttooutgroup):

    def combine_nodes(node_structure, new_node, seen_sets):
        candidate=new_node.name
        seen=[]
        for lists_of_fixed_size in seen_sets[::-1]:
            for attached_branch in lists_of_fixed_size:
                if( attached_branch.issubset(candidate) and
                   ((not attached_branch.issubset(seen)) or (not node_structure[attached_branch].has_parent()))):
                    seen.extend(list(attached_branch))
                    new_node.add_child(node_structure[attached_branch])
                    node_structure[attached_branch].add_parent(new_node)
        return node_structure

    def node_combinations_to_node_structure(node_combinations):
        length_sorted={}
        for node_combination in node_combinations:
            leaves=frozenset(node_combination.split('.'))
            k=len(leaves)
            if k in length_sorted:
                length_sorted[k].append(leaves)
            else:
                length_sorted[k]=[leaves]
        length_sorted_list=[length_sorted.get(k,[]) for k in range(1,max(length_sorted.keys())+1)]
        node_structure={}
        for leaf_node in length_sorted_list[0]:
            node_structure[leaf_node]=Node(leaf_node)
        added_sets=[length_sorted_list[0]]
        for lists_of_fixed_size in length_sorted_list[1:]:
            for branch_set in lists_of_fixed_size:
                new_node=Node(branch_set)
                combine_nodes(node_structure, new_node, added_sets)
                node_structure[branch_set]=new_node
            added_sets.append(lists_of_fixed_size)
        return node_structure

    if plottype=='top_minimal_topologies':
        df = pd.read_csv("thinned_samples_annealing.csv", sep=',', usecols=['pops'])
        nodes_list = df['pops'].tolist()
        seen_combinations = {}
        for nodes in nodes_list:
            for node in nodes.split('-'):
                seen_combinations[node] = seen_combinations.get(node, 0) + 1
        N = len(nodes_list)
        if plottype=='top_minimal_topologies':
            c=Counter(nodes_list)
            to_plots=c.most_common(1)
            c=Counter(seen_combinations)
            node_count_dic={frozenset(key.split('.')):float(count)/N for key,count in c.most_common(1000)}
            for i, (to_plot,count) in enumerate(to_plots):
                node_structure = node_combinations_to_node_structure(to_plot.split('-'))
                plot_node_structure_as_directed_graph(node_structure, drawing_name= outputprefix + '_minimal_topology.png',
                                                        node_dic=node_count_dic,  popup=False)
    elif plottype=='top_trees':
        df = pd.read_csv("thinned_samples_annealing.csv", sep=',', usecols=['pops','topology'])
        trees_list = df['topology'].tolist()
        N=len(trees_list)
        c = Counter(trees_list)
        to_plots = c.most_common(1)

        nodes=df['pops'].tolist()[0].split('-')
        leaves=list(set([leaf for node in nodes for leaf in node.split('.')]))
        leaves=sorted(leaves)

        for i, (to_plot, count) in enumerate(to_plots):
            tree=topological_identifier_to_tree_clean(to_plot, leaves=generate_predefined_list_string(deepcopy(leaves)))
            plot_as_directed_graph(tree,drawing_name=outputprefix+'_topology_.png', popup=False )
    elif plottype=='estimates':
        try:
            df = pd.read_csv("thinned_samples_annealing.csv", sep=',', usecols=['string_tree', 'topology', 'pops'])
        except ValueError as e:
            raise Exception('Unexpected columns in the posterior_distribution file. Did you turn on the --faster flag in AdmixtureBayes posterior?')
        maptree = pd.read_csv("MAPtree.txt", sep=' ', header = None)
        reviseddashedspecifier = {}
        firstcol = maptree.iloc[:, 0].tolist()
        secondcol = maptree.iloc[:, 1].tolist()
        thirdcol = maptree.iloc[:, 2].tolist()
        fourthcol =  maptree.iloc[:, 3].tolist()

        for iiiiii in range(len(firstcol)):
            if firstcol.count(firstcol[iiiiii]) > 1:
                indicess = numpy.where(numpy.array(firstcol) == firstcol[iiiiii])[0]
                indexx1 = indicess[0]
                indexx2 = indicess[1]
                if fourthcol[indexx2] > 0.5: # 1 is lower
                    reviseddashedspecifier[firstcol[indexx1]] = (secondcol[indexx1], secondcol[indexx2], thirdcol[indexx1], thirdcol[indexx2] )
                else:
                    reviseddashedspecifier[firstcol[indexx1]] = (secondcol[indexx2], secondcol[indexx1], thirdcol[indexx2], thirdcol[indexx1] )
            else:
                reviseddashedspecifier[firstcol[iiiiii]] = (secondcol[iiiiii], secondcol[iiiiii], thirdcol[iiiiii], thirdcol[iiiiii] )


        #for index, row in maptree.iterrows():
        #    if maptree.loc[index,3] < 0.5:
        #        reviseddashedspecifier[maptree.loc[index,0]] =  maptree.loc[index,1]
        topologies_list = df['topology'].tolist()
        string_tree_list=df['string_tree'].tolist()
        c = Counter(topologies_list)
        to_plots = c.most_common(1)
        cleaned_topology_list=[d[0] for d in to_plots]

        nodes=df['pops'].tolist()[0].split('-')
        leaves=list(set([leaf for node in nodes for leaf in node.split('.')]))
        leaves=sorted(leaves)

        for i, to_plot in enumerate(cleaned_topology_list):
            relevant_string_trees=[]
            for string_tree, topology in zip(string_tree_list, topologies_list):
                if topology==to_plot:
                    relevant_string_trees.append(string_tree)
            branches_intervals, admixture_proportion_intervals=branch_and_proportion_quantiles(relevant_string_trees)
            branch_names=[branches_interval[0] for branches_interval in branches_intervals]

            admixture_names=[ad[0] for ad in admixture_proportion_intervals]
            tree = identifier_to_tree(to_plot,
                                      leaves=generate_predefined_list_string(deepcopy(leaves)),
                                      branch_lengths=generate_predefined_list_string(deepcopy(branch_names)),
                                      admixture_proportions=generate_predefined_list_string(deepcopy(admixture_names)))

            org_keys = list(tree.keys())
            for key in org_keys:
                node = tree[key]
                if node_is_admixture(node):
                    new_name = get_admixture_proportion_from_key(tree, key)
                    tree = rename_key(tree, key, new_name)
            adms=get_all_admixture_origins(tree)
            adm_interpretation={}
            LabelLengthDictionary = {}
            for iiiii in branches_intervals:
                LabelLengthDictionary[iiiii[0]] = iiiii[1]
            for key, (branch_name, node_destination) in list(adms.items()):
                adm_interpretation[key]='For the lineages that pass through {}, this is the proportion that follows branch {} to node {}'.format(key, branch_name,node_destination)
            plot_as_directed_graph(tree, drawing_name=outputprefix+'_topology_labels.png', plot_edge_lengths=True,  popup=False, labeldict = LabelLengthDictionary, admixturegroup = adms, valtostring= admixture_proportion_intervals, ISLABELS = True, outgrooupp = outtttname, disttt = disttooutgroup, dashess = reviseddashedspecifier)

def identity(x):
    return x

def iterate_over_output_file(outfile, 
                             cols=[], 
                             row_summarize_functions=[],
                             thinned_d_dic=identity,
                             **constant_kwargs):
    
    df= pd.read_csv(outfile, usecols=cols, dtype={'no_admixes':object})
    df = df[cols]
    df= df.tail(1)
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
    
    def __init__(self, nodes_to_be_sorted, subnodes=[], outgroup_name=''):
        self.nodes=sorted(nodes_to_be_sorted)
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

        return {'Rtree':Rtree}, False
    
class make_full_tree(object):
    
    def __init__(self, outgroup_name='out', subnodes=[]):
        self.outgroup_name=outgroup_name
        
    def __call__(self, Rtree=None, add=None, **kwargs):
        full_tree=deepcopy(Rtree)
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
    
    def __init__(self, keys_to_include=None):
        self.keys_to_include=keys_to_include
            
    def __call__(self, full_tree=None, **kwargs):
        if full_tree is None:
            tree=kwargs['Rtree']
        else:
            tree=full_tree
        pops=Rtree_to_covariance_matrix.get_populations(tree, keys_to_include=self.keys_to_include)
        return {'pops':'-'.join(pops)}, False

class tree_unifier(object):
    
    def __init__(self):
        '''
        key:(val1,val2, val3) where 
            key=lookup topology string, 
            val1 is unique topology string, 
            val2 is the permutation of branches and 
            val3 is the (signed) permutation of the admixture proportions.
        '''
        self.seen_trees={}
        
    def __call__(self, stree):
        topology,branches,admixtures=stree.split(';')
        update_dic= analyze_tree(topology, branches, admixtures)
        self.seen_trees.update(update_dic)
        target_topology, branch_permutation, admixture_permutation=self.seen_trees[topology]
        new_branch_string=make_branch_string(branches, branch_permutation)
        new_admixtures_string=make_admixture_string(admixtures, admixture_permutation)
        return ';'.join([target_topology, new_branch_string, new_admixtures_string])
    
def make_branch_string(branches, branch_permutations):
    branch_pieces=branches.split('-')
    return '-'.join([branch_pieces[branch_permutations[i]] for i in range(len(branch_pieces))])

def make_admixture_string(admixes, admixture_permutations):
    if not admixes:
        return ''
    res=[]
    admixture_pieces=admixes.split('-')
    for i in range(len(admixture_pieces)):
        target_admixture=admixture_permutations[i]
        if target_admixture<-0.5:
            res.append(1.0-float(admixture_pieces[abs(target_admixture)-1]))
        else:
            res.append(float(admixture_pieces[abs(target_admixture)-1]))
    return '-'.join(map(str,res))

def analyze_tree(topology, branches, admixtures):
    
    id_branches='-'.join(map(str, list(range(len(branches.split('-'))))))
    id_admixtures='-'.join(map(str, list(range(1,len(admixtures.split('-'))+1))))
    id_stree=';'.join([topology,id_branches, id_admixtures])
    no_leaves=len((id_stree.split('-')[0]).split('.'))
    id_tree=identifier_to_tree_clean(id_stree.strip())
    
    strees= sorted(get_possible_permutation_strees(id_tree))
    top_topology=strees[0].split(';')[0]
    res={}
    for stree in strees:
        lookup_topology, branches_sperm, admixtures_sperm= stree.split(';')
        rf=list(map(round, list(map(float, branches_sperm.split('-')))))
        branches_permutation= list(map(int, rf))
        admixtures_permutation=get_admixtures_permutation(admixtures_sperm)
        res[lookup_topology]=(top_topology, branches_permutation, admixtures_permutation)
    return res
    
def get_admixtures_permutation(admixtures):
    res=[]
    for a in admixtures.split('-'):
        if a:
            number=int(round(float(a)))
            if number>100000:
                res.append(-(number%100000))
            else:
                res.append(number)
    return res

def get_possible_permutation_strees(tree):
    leaves,admixture_keys=get_categories(tree)
    k=len(admixture_keys)
    format_code='{0:0'+str(k)+'b}'
    
    n_trees=[]
    for i in range(2**k):
        pruned_tree = deepcopy(tree)
        bina= format_code.format(i)
        for adm_key,str_bin in zip(admixture_keys, list(bina)):
            int_bin=int(str_bin)
            if int_bin==1:
                pruned_tree[adm_key]=change_admixture(pruned_tree[adm_key])
        n_tree= unique_identifier_and_branch_lengths(pruned_tree)
        n_trees.append(n_tree)
    return n_trees
    
def run_posterior_main():

    subnodes_with_outgroup=[]
    subnodes_wo_outgroup=[]

    totallist = []
    a = pd.read_csv("mcmc_samples_annealing.csv", nrows=3)
    stringg = (a.loc[0,["descendant_sets"]])[0]
    stringg = (stringg.split('-'))
    for i in stringg:
        totallist.extend(i.split('.'))
    nodes = []
    for i in totallist:
        if i not in nodes:
            nodes.append(i)
    nodes.sort()

    nodes_wo_outgroup = deepcopy(nodes)
    nodes_with_outgroup = deepcopy(nodes_wo_outgroup)
    nodes_with_outgroup.sort()
    nodes_wo_outgroup.sort()
    subnodes_with_outgroup.sort()
    subnodes_wo_outgroup.sort()

    row_sums=[]

    class pointers(object):

        def __init__(self):
            self.count=0
            self.dic={}

        def __call__(self, name):
            self.dic[name]=self.count
            self.count+=1

        def __getitem__(self, key):
            return self.dic[key]

    name_to_rowsum_index=pointers()

    row_sums.append(make_Rtree(deepcopy(nodes_wo_outgroup), subnodes=subnodes_wo_outgroup, outgroup_name=''))
    name_to_rowsum_index('Rtree')
    row_sums.append(make_full_tree(outgroup_name='', subnodes=[]))
    name_to_rowsum_index('full_tree')

    nodes=nodes_with_outgroup
    row_sums.append(make_string_tree(deepcopy(nodes), tree_unifier())) #calling make_string_tree
    name_to_rowsum_index('string_tree')
    row_sums.append(topology(nodes=nodes))
    name_to_rowsum_index('topology')
    row_sums.append(get_pops(keys_to_include=nodes))
    name_to_rowsum_index('pops')

    def save_thin_columns(d_dic):
        return {summ:d_dic[summ] for summ in list(set(['no_admixes', 'topology', 'pops','string_tree']+[]))}
    all_results=iterate_over_output_file("mcmc_samples_annealing.csv",
                                             cols=['tree', 'add', 'layer', 'no_admixes'],
                                             row_summarize_functions=row_sums,
                                             thinned_d_dic=save_thin_columns)

    if True:
        summaries=list(all_results[0].keys())
        with open("thinned_samples_annealing.csv", 'w') as f:
            f.write(','.join(summaries)+'\n')
            for row in all_results:
                s_summs=[str(row[summ]) for summ in summaries]
                f.write(','.join(s_summs)+ '\n')
    
def removefile(filename):
    if os.path.exists(filename):
        os.remove(filename)

def get_summary_scheme(no_chains=1):
    summaries=[construct_starting_trees_choices.s_posterior(),
               construct_starting_trees_choices.s_likelihood(),
               construct_starting_trees_choices.s_prior(),
               construct_starting_trees_choices.s_no_admixes(),
               construct_starting_trees_choices.s_variable('add', output='double'), 
               construct_starting_trees_choices.s_total_branch_length(),
               construct_starting_trees_choices.s_basic_tree_statistics(get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               construct_starting_trees_choices.s_basic_tree_statistics(Rtree_to_covariance_matrix.get_populations_string, 'descendant_sets', output='string'),
               construct_starting_trees_choices.s_basic_tree_statistics(unique_identifier_and_branch_lengths, 'tree', output='string'),
               construct_starting_trees_choices.s_basic_tree_statistics(get_admixture_proportion_string, 'admixtures', output='string')]
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    sample_verbose_scheme_first=deepcopy(sample_verbose_scheme)
    return [sample_verbose_scheme_first]+[{}]*(no_chains-1), summaries

def main(args):
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"

    parser = ArgumentParser(usage='pipeline for Admixturebayes')

    #input/output options
    parser.add_argument('--input_file', type=str, required=True, help='the input file of the pipeline. It should be of the same type as the treemix input file with a header of population names and each line representing a snp (unless --covariance_pipeline is altered).')
    parser.add_argument('--result_file', type=str, default='mcmc_samples_annealing.csv', help='file in which to save results.')
    parser.add_argument('--output_prefix', type=str, default='out', help='file in which to save results.')
    parser.add_argument('--outgroup', type=str, default='',
                        help='The name of the population that should be outgroup for the covariance matrix. If the covariance matrix is supplied at stage 8 , this argument is not needed.')
    parser.add_argument('--save_covariance', default=False, action='store_true', help='saving the covariance matrix')
    #Important arguments
    parser.add_argument('--bootstrap_blocksize', type=int, default=1000,
                        help='the size of the blocks to bootstrap in order to estimate the degrees of freedom in the wishart distribution')

    #annealing
    parser.add_argument('--starting_temp', type=str, required=True, help='the inputuu')
    parser.add_argument('--ending_temp', type=str, required=True, help='the inputppp')
    parser.add_argument('--temp_scaling', type=str, required=True, help='the inputnnn')
    parser.add_argument('--iter_per_temp', type=str, required=True, help='the inputxx')
    parser.add_argument('--dontcleanup', default = False, action='store_true', help='the inputxx')

    #convenience arguments
    parser.add_argument('--verbose_level', default='silent', choices=['normal', 'silent'],
                        help='this will set the amount of status out prints throughout running the program.')

    options=parser.parse_args(args)

    temporaryfoldername =  options.output_prefix
    os.mkdir(os.getcwd() + "/" + temporaryfoldername)

    
    assert not (any((i < 8 for i in [6,8,9])) and not options.outgroup), 'In the requested analysis, the outgroup needs to be specified by the --outgroup flag and it should match one of the populations'

    #Here is the only thing we should be changing.
    temp = pd.read_csv(options.input_file, sep ="\s+")
    colnames = list(temp.columns.values)
    for pop_name in colnames:
        assert pop_name.isalnum(), 'Population names can only contain alphanumeric characters (A-Z, a-z, and 0-9). Special characters such as commas, hyphens, and underscores are not allowed.'
    
    assert options.outgroup in colnames, 'The outgroup name is not in the given list of populations. Population names are case-sensitive.'
    
    temp = temp[sorted(colnames)]
    temp.to_csv(os.getcwd() +"/"+ temporaryfoldername + "/temp_input.txt", sep =" ", index = False)

    mp= [simple_adaptive_proposal(['deladmix', 'addadmix', 'rescale', 'rescale_add', 'rescale_admixtures', 'rescale_constrained', 'sliding_regraft'],
     [1, 1, 1, 1, 1, 1, 1]) for _ in range(1)]

    with open(os.getcwd() +"/"+ temporaryfoldername + "/temp_input.txt", 'r') as f:
        full_nodes = f.readline().rstrip().split()
    reduced_nodes=deepcopy(full_nodes)
    reduced_nodes.remove(options.outgroup)

    estimator_arguments=dict(reducer=options.outgroup, nodes=full_nodes, add_variance_correction_to_graph=True, save_variance_correction=True)
                             
    covariance=get_covariance(os.getcwd() +"/"+ temporaryfoldername + "/temp_input.txt", 
    varcovfilename = os.getcwd() +"/"+ temporaryfoldername + "/variance_correction.txt",
    full_nodes=full_nodes,
     reduce_covariance_node=options.outgroup,
     estimator_arguments=estimator_arguments,filename =  os.getcwd() +"/"+ temporaryfoldername + "/covariance_and_multiplier.txt")
    estimator_arguments['save_variance_correction']=False
    df=estimate_degrees_of_freedom_scaled_fast(os.getcwd() +"/"+ temporaryfoldername + "/temp_input.txt",
                                               varcovfilename = os.getcwd() +"/"+ temporaryfoldername + "/variance_correction.txt",
                                            bootstrap_blocksize=options.bootstrap_blocksize,
                                            cores=2,
                                            est=estimator_arguments, 
                                            verbose_level=options.verbose_level)
    multiplier=covariance[1]

    #NEWEDIT, change the thing below to False, then []
    starting_trees=construct_starting_trees_choices.get_starting_trees([], 1, adds=[], nodes=reduced_nodes)

    summary_verbose_scheme, summaries=get_summary_scheme(no_chains=1)

    posterior = posterior_class(emp_cov=covariance[0], M=df, multiplier=covariance[1], nodes=reduced_nodes, 
                                 varcovname=os.getcwd() +"/"+ temporaryfoldername + "/variance_correction.txt")
    
    removefile("variance_correction.txt")
    removefile("temp_starttree.txt")
    removefile("temp_start_tree.txt")
    removefile("temp_add.txt")

    if options.save_covariance:
        removefile(os.getcwd() + "/covariance_matrix.txt")
        Liness = open(os.getcwd() +"/"+ temporaryfoldername +"/covariance_and_multiplier.txt", 'r').readlines()
        covarfile = open(os.getcwd() + "/covariance_matrix.txt", "a")
        covarfile.writelines(Liness)
        covarfile.close()

    removefile("covariance_and_multiplier.txt")
    removefile(os.getcwd() + "/temp_input.txt")

    if os.path.exists(os.getcwd() + "/temp_adbayes"):
        os.rmdir(os.getcwd() + "/temp_adbayes")
    
    StartingTemp =  float(options.starting_temp) #    100
    EndingTemp = float(options.ending_temp) #  0.0001
    TempDecrease = float(options.temp_scaling) # 0.9
    NumberAtEach = int(options.iter_per_temp) # 5000

    MCMCMC(TempDecrease, starting_trees=starting_trees,
            posterior_function= posterior,
            summaries=summaries,
            temperature_scheme=[StartingTemp],
            printing_schemes=summary_verbose_scheme,
            iteration_scheme=[NumberAtEach]*int(log(EndingTemp / StartingTemp) / log(TempDecrease)),
            proposal_scheme= mp,
            multiplier=multiplier,
            result_file=options.result_file,
            n_arg=int(log(EndingTemp / StartingTemp) / log(TempDecrease)), verboseee=options.verbose_level)
    
    removefile(os.getcwd() +"/"+ temporaryfoldername + "/" + "covariance_and_multiplier.txt")
    removefile(os.getcwd() +"/"+ temporaryfoldername + "/" + "temp_input.txt")
    removefile(os.getcwd() +"/"+ temporaryfoldername + "/" + "variance_correction.txt")
    if os.path.exists(os.getcwd() +"/"+ temporaryfoldername + "/temp_adbayes" ):
        os.rmdir(os.getcwd() +"/"+ temporaryfoldername + "/temp_adbayes" )
    if os.path.exists(os.getcwd() +"/"+ temporaryfoldername):
        os.rmdir(os.getcwd() +"/"+ temporaryfoldername)
        
    return(options.outgroup, options.output_prefix, multiplier, options.dontcleanup)

if __name__=='__main__':
    import sys, os

    #We here run the main algorithm
    outtname, outputprefixx, multiplierrr, dontcleanup = main(sys.argv[1:])
    run_posterior_main()

    mainnnoutput()

    scaleddistancetooutgroup = (pd.read_csv("MAPadd.txt", header = None))[0][0]

    actualoutgroupdistance = scaleddistancetooutgroup / multiplierrr

    outputprefixx
    with open(outputprefixx + '_outgroup.txt', 'w') as ffff:
        ffff.write(str(actualoutgroupdistance) + "\n")

    main_plot("top_minimal_topologies", "out", outputprefixx, 0) ####file is right, picture is wrong, jsut read in map tree and edit it there!!!
    main_plot("top_trees", "out", outputprefixx, 0) #most of the admixture should go to the side of the tree with 2 leaf nodes
    main_plot("estimates", outtname, outputprefixx, actualoutgroupdistance)

    (pd.read_csv("MAPtree.txt", header = None)).to_csv(outputprefixx + '_tree.txt',  header=False, index = False)

    if os.path.exists("mcmc_samples_annealing.csv"):
        os.remove("mcmc_samples_annealing.csv")
    if os.path.exists("thinned_samples_annealing.csv"):
        os.remove("thinned_samples_annealing.csv")
    
    if not dontcleanup:
        if os.path.exists("MAPadd.txt"):
            os.remove("MAPadd.txt")
        if os.path.exists("MAPtree.txt"):
            os.remove("MAPtree.txt")
