from collections import Counter
import pandas as pd
from tree_statistics import generate_predefined_list_string, topological_identifier_to_tree_clean, identifier_to_tree
from copy import deepcopy
from Rtree_operations import node_is_admixture, rename_key, get_admixture_proportion_from_key, get_all_admixture_origins, to_networkx_format
import sys

from graphviz import Digraph
import os
    
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
        
def plot_as_directed_graph(tree, drawing_name='tmp.png', popup=True, plot_edge_lengths=False, verbose=True ,labeldict = {}, admixturegroup = {}, valtostring=[]):

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
    print(admixture_nodes)
    for adm_n in admixture_nodes:
        if len(valtostring) == 0:
            admixture_graph.node(adm_n)
        else:
            for jjjj in valtostring:
                if jjjj[0] == adm_n:
                    admixture_graph.node(adm_n, label = str(round(jjjj[1], 3)))
        
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
            if labeldict[label] < 0.0001:
                if from_node in admixturegroup and admixturegroup[from_node][1] == to_node:
                    G.edge(to_node, from_node, label=  f'{labeldict[label]:.2e}', style = "dashed")
                else:
                    G.edge(to_node, from_node, label=  f'{labeldict[label]:.2e}')
            else:
                if from_node in admixturegroup and admixturegroup[from_node][1] == to_node:
                    G.edge(to_node, from_node, label=  str(round(labeldict[label], 4)), style = "dashed")
                else:
                    G.edge(to_node, from_node, label=  str(round(labeldict[label], 4)))
    else:
        G.edges(edges)
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

def main(plottype):

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
        df = pd.read_csv("thinned_samples.csv", sep=',', usecols=['pops'])
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
                plot_node_structure_as_directed_graph(node_structure, drawing_name='minimal_topology_' +str(i+1)+'.png',
                                                        node_dic=node_count_dic,  popup=False)
    elif plottype=='top_trees':
        df = pd.read_csv("thinned_samples.csv", sep=',', usecols=['pops','topology'])
        trees_list = df['topology'].tolist()
        N=len(trees_list)
        c = Counter(trees_list)
        to_plots = c.most_common(1)

        nodes=df['pops'].tolist()[0].split('-')
        leaves=list(set([leaf for node in nodes for leaf in node.split('.')]))
        leaves=sorted(leaves)

        for i, (to_plot, count) in enumerate(to_plots):
            tree=topological_identifier_to_tree_clean(to_plot, leaves=generate_predefined_list_string(deepcopy(leaves)))
            plot_as_directed_graph(tree,drawing_name='topology_' + str(i + 1) + '.png', popup=False )
    elif plottype=='estimates':
        try:
            df = pd.read_csv("thinned_samples.csv", sep=',', usecols=['string_tree', 'topology', 'pops'])
        except ValueError as e:
            raise Exception('Unexpected columns in the posterior_distribution file. Did you turn on the --faster flag in AdmixtureBayes posterior?')

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
            

            plot_as_directed_graph(tree, drawing_name='topology_labels_' + str(i + 1) + '.png', plot_edge_lengths=True,  popup=False, labeldict = LabelLengthDictionary, admixturegroup = adms, valtostring= admixture_proportion_intervals)

if __name__=='__main__':
    main("top_minimal_topologies")
    main("top_trees")
    main("estimates")