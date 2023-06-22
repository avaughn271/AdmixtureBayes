from collections import Counter
import pandas as pd
from argparse import ArgumentParser
from tree_statistics import generate_predefined_list_string, topological_identifier_to_tree_clean, identifier_to_tree
from copy import deepcopy
from Rtree_operations import node_is_admixture, rename_key, get_admixture_proportion_from_key, get_all_admixture_origins, to_networkx_format
import sys
import numpy as np
from graphviz import Digraph
import os

from math import floor
    
class make_read_name(object):
    
    def __init__(self):
        self.count=0
        
    def __call__(self, admix=False):
        self.count+=1
        if admix:
            return 'a'+str(self.count)
        else:
            return 'n'+str(self.count)
        
def rgb_to_hex(r,g,b):
    return '#%02x%02x%02x' % (r, g, b)

def freq_to_white_green_hex(f):
    if f<0.6:
        r,g,b=int((0.6-f)*255.0/0.6),255, int((0.6-f)*255/0.6)
    else:
        r,g,b=0,255*(1-(f-0.6)),0
    return rgb_to_hex(r,int(g),b)

class dummy_dic(object):
    
    def __init__(self, val):
        self.val=val
        
    def __getitem__(self,key):
        return freq_to_white_green_hex(self.val)

class color_dic(object):
    
    def __init__(self, node_freq_dic):
        self.nfd=node_freq_dic
        
    def __getitem__(self, key):
        freq=self.nfd[key]
        return freq_to_white_green_hex(freq)

def node_structure_to_networkx(node_structure, node_freq_dictionary=None):
    if node_freq_dictionary is None:
        colors=dummy_dic(0.2)
    else:
        colors=color_dic(node_freq_dictionary)
    edges=[]
    admixture_nodes=[]
    pure_leaves=[]
    admixture_leaves=[]
    root=[]
    coalescence_nodes=[]
    namer=make_read_name()
    inner_node_colors={}
    string_names={}
    for key,node in list(node_structure.items()):
        if key==frozenset(['']):
            continue
        plotting_type=node.get_plotting_type()
        if plotting_type=='coalescence':
            string_names[key]=namer()
            if node_freq_dictionary is not None:
                string_names[key]+='('+str(int(floor(100*node_freq_dictionary[key])))+')'
            coalescence_nodes.append(string_names[key])
            inner_node_colors[string_names[key]]=colors[key]
        elif plotting_type=='admixture':
            string_names[key]=namer(admix=True)
            if node_freq_dictionary is not None:
                string_names[key]+='('+str(int(floor(100*node_freq_dictionary[key])))+')'
            admixture_nodes.append(string_names[key])
            inner_node_colors[string_names[key]]=colors[key]
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
    return pure_leaves, admixture_leaves, coalescence_nodes, admixture_nodes, root, edges, inner_node_colors

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
    pure_leaves, admix_leaves, pure_coalescence, admix_coalescence, root, edges, inner_node_colors= node_structure_to_networkx(node_structure, node_dic)
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
        admixture_graph.node(adm_n, fillcolor=inner_node_colors[adm_n])
        
    coalescence_graph=Digraph()
    coalescence_graph.node_attr.update(style='filled', fillcolor='greenyellow')
    for cn in pure_coalescence:
        coalescence_graph.node(cn, fillcolor=inner_node_colors[cn])
        
    G.subgraph(leaves_graph)
    G.subgraph(aleaves_graph)
    G.subgraph(admixture_graph)
    G.subgraph(coalescence_graph)
    G.node('r', shape='egg', color='black', style='filled', fontcolor='white')
    G.edges(edges)
    G.format = image_format
    G.render(view=popup)
    if verbose:
        print('written plot to file', drawing_name[0:(len(drawing_name)-4)] + ".pdf")
        if os.path.exists(filename):
            os.remove(filename)
        else:
            print("The file does not exist")
        
def plot_as_directed_graph(tree, drawing_name='tmp.png', popup=True, plot_edge_lengths=False, verbose=True):

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
        admixture_graph.node(adm_n)
        
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
            G.edge(to_node, from_node, label=label)
    else:
        G.edges(edges)
    G.format = image_format
    G.render(view=popup)
    if verbose:
        print('written plot to file', drawing_name[0:(len(drawing_name)-4)] + ".pdf")
        if drawing_name[0:10] == "topology_l":
            strlength = len(drawing_name)
            numm = (drawing_name[16:(strlength-4)]) #ANDREWDEBUG
            print("written branch lengths to file branch_estimates_" + numm + ".txt")
            print("written admixture proportions to file admixture_estimates_" + numm + ".txt")
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
    for n,(lower,mean, upper) in enumerate(zip(bdat.quantile(0.025),bdat.mean(),bdat.quantile(0.975))):
        bresults.append(('c'+str(n+1), lower,mean,upper))
    if len(admixtures[0])>0:
        for n,(lower,mean, upper) in enumerate(zip(adat.quantile(0.025),adat.mean(),adat.quantile(0.975))):
            aresults.append(('ax'+str(n+1), lower,mean,upper))
    return bresults, aresults

def checkequality(tree1, tree2, tree2orig):
  
  if tree1.shape[0] == 1 and tree2.shape[0] == 1 and tree1[0,0] == tree2[0,0] and tree1[0,1] and tree2[0,1]:
     return(1)
  
  if tree1.shape[0] != tree2.shape[0]:
      return(0)
  
  for i in range(tree2.shape[0]):
    if ( (tree2[i,0] not in tree2[:,1]) and np.sum(tree2[i,0] == tree2) == 1): #what is left is a leaf or divergence node
      leaff2 = tree2[i,0]
      parent2 = tree2[i,1]
      if np.sum(leaff2 == tree1) != 1:
        return(0)
      if np.sum(leaff2 == tree1[:,0]) != 1:
        return(0)
      j = int(np.where(leaff2 == tree1[:, 0])[0][0])
      parent1 = tree1[j,1]
      if parent1[-3:] == "_xx" and parent1 != parent2:
         return 0
      tree1[tree1 == parent1] = parent2
      return checkequality(np.delete(tree1, j, 0), np.delete(tree2, i, 0), tree2orig)
  for i in range(tree1.shape[0]): # remove rows that are entirely equivalent between the trees
        for j in range(tree2.shape[0]):
            if tree1[i, 0] == tree2[j, 0] and tree1[i, 1] == tree2[j, 1]:
                return checkequality(np.delete(tree1, i, 0), np.delete(tree2, j, 0), tree2orig)
  for i in range(tree1.shape[0]): #if everything is fixed in a row but does not match the original version.
        if tree1[i, 0][-3:] == "_xx" and tree1[i, 1][-3:] == "_xx":
            if tree1[i, 1] not in tree2orig[np.where(tree1[i, 0] == tree2orig[:, 0]), 1]:
                return 0
  for i in range(tree2.shape[0]):
          if (tree2[i,0] not in tree2[:,1]) and (tree2[i,0] in tree1[:,0]): #what is left is an admixture node, present in both
              admixturestring = tree2[i,0]
              tree1_index1 = (np.where(admixturestring== tree1[:, 0])[0][0])
              tree1_index2 = (np.where(admixturestring== tree1[:, 0])[0][1])
              tree2_index1 = (np.where(admixturestring== tree2[:, 0])[0][0])
              tree2_index2 = (np.where(admixturestring== tree2[:, 0])[0][1])
              tree1_new = np.copy(tree1)
              tree1_new2 = np.copy(tree1)
              tree1_new[tree1_new == tree1[tree1_index1,1]] = tree2[tree2_index1,1] # try both ways of assigning the admixture node parents.
              tree1_new[tree1_new == tree1[tree1_index2,1]] = tree2[tree2_index2,1]

              tree1_new2[tree1_new2 == tree1[tree1_index2,1]] = tree2[tree2_index1,1]
              tree1_new2[tree1_new2 == tree1[tree1_index1,1]] = tree2[tree2_index2,1]

              AA  = checkequality(np.delete(tree1_new, [tree1_index1,tree1_index2], 0), np.delete(tree2, [tree2_index1,tree2_index2], 0), tree2orig)
              BB  = checkequality(np.delete(tree1_new2, [tree1_index1,tree1_index2], 0), np.delete(tree2, [tree2_index1,tree2_index2], 0), tree2orig)

              return AA or BB
  print("Problem resolving topology equivalences. Contact Andrew Vaughn at ahv36@berkeley.edu.")
  return(-1000000)

    
def IsEquivalentTopology(tree1, tree2):
    Tree1List = []
    for node in tree1:
        if tree1[node][1] is not None:
            Tree1List.append([node,  tree1[node][0]] )
            Tree1List.append([node,  tree1[node][1]] )
        else:
            Tree1List.append([node,  tree1[node][0]] )
    Tree2List = []
    for node in tree2:
        if tree2[node][1] is not None:
            Tree2List.append([node,  tree2[node][0]] )
            Tree2List.append([node,  tree2[node][1]] )
        else:
            Tree2List.append([node,  tree2[node][0]] )
    
    EDGES1 = np.array(Tree1List, dtype='<U20')
    EDGES2 = np.array(Tree2List, dtype='<U20')

    edgess2 =  np.copy(EDGES2[:,1])
    for row in range(EDGES2.shape[0]):
        for col in range(EDGES2.shape[1]):
            if EDGES2[row,col] in edgess2:
                EDGES2[row,col] = EDGES2[row,col] + "_xx"
    return checkequality(EDGES1 , EDGES2, EDGES2)

def counterequivalence(originalcounter, Equivs):
    OriginalOrder = originalcounter
    UnpackList = np.array(sorted(originalcounter.elements()))

    while True:
        didsub = 0
        for i in range(len(Equivs)):
            firstel = Equivs[i][0]
            secondel = Equivs[i][1]
            maxindex = 0
            if OriginalOrder[secondel] > OriginalOrder[firstel]:
                maxindex = 1
            if Equivs[i][1 - maxindex] in UnpackList:
                didsub = didsub + 1
                UnpackList[UnpackList == Equivs[i][1 - maxindex]] = Equivs[i][ maxindex]
                break
        if didsub == 0:
            break

    return (Counter(UnpackList))


def main(args):
    parser = ArgumentParser(usage='pipeline for plotting posterior distribution summaries.')

    parser.add_argument('--posterior', required=True, type=str, help='The file containing posterior distributions from the "AdmixtureBayes posterior" command. It needs the two columns "pops" and topology.')
    parser.add_argument('--plot', choices=['consensus_trees', 'top_minimal_topologies', 'top_trees','estimates'], required=True,
                        help='The type of plot to make. Choose between: 1) consensus_trees. '
                             'It plots an admixture graph based on all nodes that have a higher (marginal) posterior probability of X. '
                             'Different X\'s can be supplied with the command --consensus_threshold \n'
                             '2) top_minimal_topologies. It plots the X highest posterior combinations of node types '
                             'and creates the corresponding minimal topologies.  X can be supplied through the command --top_minimal_topologies_to_plot'
                             '3) top_trees. It plots the X highest posterior topologies. X can be supplied by the command --top_trees_to_plot.'
                             '4) estimates. It creates a table  with continuous parameters estimated from the posterior sample'
                             'It also plots the concerned topologies with labels. It does this for either the X highest posterior topologies ')
    parser.add_argument('--outgroup', default='outgroup', help='name of the outgroup to plot')
    parser.add_argument('--consensus_thresholds', default=[0.25, 0.5, 0.75, 0.9, 0.95, 0.99], type=float, nargs='+',
                        help='The posterior thresholds for which to draw different consensus trees.')
    parser.add_argument('--top_minimal_topologies_to_plot', type=int, default=3,
                        help='The number of node trees (or minimal topologies) to plot')
    parser.add_argument('--top_trees_to_plot', type=int, default=3,
                        help='The number of trees (or topologies) to plot ')
    parser.add_argument('--top_trees_to_estimate', type=int, default=3,
                        help='The number of trees (or topologies) to plot ')
    parser.add_argument('--write_estimates_to_file', default=[], type=str, nargs='+',
                        help='The file in which to put the tables when plotting estimates. ')
    parser.add_argument('--write_rankings', type=str, default='', help='if a file is supplied here, the natural rankings for each of the plots is written here.')
    parser.add_argument('--rankings_to_write_to_file', type=int, default=1000,
                        help='the number of rankings(nodes, min topology or topology depending on --plot) to write to the ranking file.')
    parser.add_argument('--popup', default=False, action='store_true')
    parser.add_argument('--output_prefix', type=str,  default="default")

    options= parser.parse_args(args)

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
        #length_sorted_list is of the form [[[A],[B],[C]],[[A,B],[B,C]],...,[[A,B,C]]]
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

    if options.plot=='consensus_trees' or options.plot=='top_minimal_topologies':
        df = pd.read_csv(options.posterior, sep=',', usecols=['pops'])
        nodes_list = df['pops'].tolist()
        seen_combinations = {}
        for nodes in nodes_list:
            for node in nodes.split('-'):
                seen_combinations[node] = seen_combinations.get(node, 0) + 1
        N = len(nodes_list)
        if options.plot=='consensus_trees':
            node_combinations = []
            for threshold in options.consensus_thresholds:
                total_threshold = int(N * threshold)
                final_node_combinations = [k for k, v in list(seen_combinations.items()) if v > total_threshold]
                node_combinations.append(final_node_combinations)
            node_count_dic={frozenset(k.split('.')):float(v)/N for k,v in list(seen_combinations.items())}
            for i, final_node_combinations in enumerate(node_combinations):
                final_node_structure = node_combinations_to_node_structure(final_node_combinations)
                drawingnamm = options.output_prefix + '_'
                if drawingnamm == "default_":
                    drawingnamm = 'consensus_'
                plot_node_structure_as_directed_graph(final_node_structure, drawing_name = drawingnamm + str(int(100*options.consensus_thresholds[i]))+'.png', node_dic=node_count_dic,  popup=options.popup)
            if options.write_rankings:
                with open(options.write_rankings, 'w') as f:
                    c = Counter(seen_combinations)
                    to_write = c.most_common(options.rankings_to_write_to_file)
                    for node, frequency in to_write:
                        f.write(node+','+str(float(frequency)/N)+'\n')
        elif options.plot=='top_minimal_topologies':
            c=Counter(nodes_list)
            to_plots=c.most_common(options.top_minimal_topologies_to_plot)
            if options.write_rankings:
                with open(options.write_rankings, 'w') as f:
                    for tree, frequency in c.most_common(options.rankings_to_write_to_file):
                        f.write(tree + ',' + str(float(frequency) / N) + '\n')
            c=Counter(seen_combinations)
            node_count_dic={frozenset(key.split('.')):float(count)/N for key,count in c.most_common(1000)}
            for i, (to_plot,count) in enumerate(to_plots):
                node_structure = node_combinations_to_node_structure(to_plot.split('-'))
                drawingnamm = options.output_prefix + '_'
                if drawingnamm == "default_":
                    drawingnamm = 'minimal_topology_'
                plot_node_structure_as_directed_graph(node_structure, drawing_name = drawingnamm +str(i+1)+'.png',
                                                        node_dic=node_count_dic,  popup=options.popup)
    elif options.plot=='top_trees':
        df = pd.read_csv(options.posterior, sep=',', usecols=['pops','topology'])
        trees_list = df['topology'].tolist()
        no_leaves=len(trees_list[0].split('-')[0].split('.'))
        N=len(trees_list)
        c = Counter(trees_list)

        nodes=df['pops'].tolist()[0].split('-')
        leaves=list(set([leaf for node in nodes for leaf in node.split('.')]))
        if len(leaves)==no_leaves:
            pass #everything is good
        elif len(leaves)==no_leaves-1:
            leaves.append(options.outgroup)
        leaves=sorted(leaves)

        AllEquivs = []

        for i, to_plot1 in enumerate(c):
            for j, to_plot2 in enumerate(c):
                if j > i:
                    TREE1 = topological_identifier_to_tree_clean(to_plot1, leaves=generate_predefined_list_string(deepcopy(leaves)))
                    TREE2 = topological_identifier_to_tree_clean(to_plot2, leaves=generate_predefined_list_string(deepcopy(leaves)))
                    if IsEquivalentTopology(TREE1, TREE2):
                        AllEquivs.append([to_plot1, to_plot2])
        c = counterequivalence(c, AllEquivs)
        to_plots = c.most_common(options.top_trees_to_plot)

        if options.write_rankings:
            with open(options.write_rankings, 'w') as f:
                for tree, frequency in c.most_common(options.rankings_to_write_to_file):
                    f.write(tree + ',' + str(float(frequency) / N) + '\n')

        for i, (to_plot, count) in enumerate(to_plots):
            tree=topological_identifier_to_tree_clean(to_plot, leaves=generate_predefined_list_string(deepcopy(leaves)))
            drawingnamm = options.output_prefix + '_'
            if drawingnamm == "default_":
                drawingnamm = 'topology_'
            plot_as_directed_graph(tree,drawing_name = drawingnamm + str(i + 1) + '.png', popup=options.popup)
    elif options.plot=='estimates':
        try:
            df = pd.read_csv(options.posterior, sep=',', usecols=['string_tree', 'topology', 'pops'])
        except ValueError as e:
            raise Exception('Unexpected columns in the posterior_distribution file. Did you turn on the --faster flag in AdmixtureBayes posterior?')

        topologies_list = df['topology'].tolist()
        string_tree_list=df['string_tree'].tolist()
        c = Counter(topologies_list)

        no_leaves=len(topologies_list[0].split('-')[0].split('.'))

        nodes=df['pops'].tolist()[0].split('-')
        leaves=list(set([leaf for node in nodes for leaf in node.split('.')]))
        if len(leaves)==no_leaves:
            pass #everything is good
        elif len(leaves)==no_leaves-1:
            leaves.append(options.outgroup)
        leaves=sorted(leaves)

        AllEquivs = []

        for i, to_plot1 in enumerate(c):
            for j, to_plot2 in enumerate(c):
                if j > i:
                    TREE1 = topological_identifier_to_tree_clean(to_plot1, leaves=generate_predefined_list_string(deepcopy(leaves)))
                    TREE2 = topological_identifier_to_tree_clean(to_plot2, leaves=generate_predefined_list_string(deepcopy(leaves)))
                    if IsEquivalentTopology(TREE1, TREE2):
                        AllEquivs.append([to_plot1, to_plot2])
        c = counterequivalence(c, AllEquivs)






        to_plots = c.most_common(options.top_trees_to_estimate)
        cleaned_topology_list=[d[0] for d in to_plots]
        no_leaves = len(topologies_list[0].split('-')[0].split('.'))

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
            for key, (branch_name, node_destination) in list(adms.items()):
                adm_interpretation[key]='For the lineages that pass through {}, this is the proportion that follows branch {} to node {}'.format(key, branch_name,node_destination)
            drawingnamm = options.output_prefix + '_'
            if drawingnamm == "default_":
                drawingnamm = 'topology_labels_'
            plot_as_directed_graph(tree, drawing_name=drawingnamm + str(i + 1) + '.png', plot_edge_lengths=True,  popup=options.popup)
            if options.write_estimates_to_file:
                branch_file=options.write_estimates_to_file[i*2+0]
                admixtures_file=options.write_estimates_to_file[i*2+1]
            else:
                if options.output_prefix == "default":
                    branch_file='branch_estimates_'+str(i+1)+'.txt'
                    admixtures_file='admixture_estimates_'+str(i+1)+'.txt'
                else:
                    branch_file = options.output_prefix +'_branch_estimates_'+str(i+1)+'.txt'
                    admixtures_file = options.output_prefix   + '_admixture_estimates_'+str(i+1)+'.txt'
            with open(branch_file, 'w') as f:
                f.write(','.join(['branch label','lower 95%','mean','upper 95%'])+'\n')
                for v in branches_intervals:
                    f.write(','.join(map(str,v))+'\n')
            with open(admixtures_file, 'w') as f:
                f.write(','.join(['branch label','lower 95%', 'mean', 'upper 95%','interpretation'])+'\n')
                for v in admixture_proportion_intervals:
                    f.write(','.join(map(str,list(v)+[adm_interpretation[v[0]]]))+'\n')

    sys.exit()

if __name__=='__main__':
    main(sys.argv[1:])
