from csv import writer
from Rtree_operations import to_aarhus_admixture_graph, to_networkx_format, node_is_admixture, node_is_coalescence, node_is_leaf_node
from node_structure import node_structure_to_networkx
from subprocess import call
from graphviz import Digraph
import os

import imp

try:
    imp.find_module('PIL')
    PIL_found = True
    from PIL import Image
except ImportError:
    PIL_found = False

file_suffix=[s+'.csv' for s in ['leaves', 'inner_nodes','edges','adm_props']]

def plot_graph(*args, **kwargs):
    plot_as_directed_graph(*args, **kwargs)

def plot_as_admixture_tree(tree, file_prefix='', drawing_name='tmp.png', popup=True):
    aarhus_tree = to_aarhus_admixture_graph(tree)
    file_names=[file_prefix+s for s in file_suffix]
    write_aarhus_tree_to_files(aarhus_tree, file_names)
    make_R_draw_from_files(drawing_name, file_names)
    if popup and PIL_found:
        img=Image.open(drawing_name)
        img.show()
        
def plot_node_structure_as_directed_graph(node_structure, file_prefix='', drawing_name='tmp_.png', popup=True, node_dic=None, verbose=True):
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
        
def plot_as_directed_graph(tree, file_prefix='', drawing_name='tmp.png', popup=True, plot_edge_lengths=False, verbose=True):

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
        #with open("topology_labels_2_old") as f: contents = f.readlines()
        #contents[0] = contents[0] + "graph [ dpi = 300 ]; \n"
        #f = open("demofile2.txt", "a")
        #f.write("Now the file has more content!")
        #f.close()
        if os.path.exists(filename):
            os.remove(filename)
        else:
            print("The file does not exist")
    
    
def write_aarhus_tree_to_files(aarhus_tree, file_names):
    for rows, name in zip(aarhus_tree, file_names):
        with open(name, "wb") as f:
            writer2 = writer(f)
            writer2.writerows(rows)
            
def pretty_print(tree):
    keys,vals=list(tree.keys()),list(tree.values())
    chosen_indexes=[]
    print('{ ')
    for key,val in zip(keys,vals):
        if node_is_leaf_node(val):
            print('  '+key+': '+str(val))
    print('  ,')
    for key,val in zip(keys,vals):
        if node_is_admixture(val):
            print('  '+key+': '+str(val))
    print('  ,')
    for key,val in zip(keys,vals):
        if node_is_coalescence(val):
            print('  '+key+': '+str(val))
    print('}')
    
def pretty_string(tree):
    keys,vals=list(tree.keys()),list(tree.values())
    res=''
    res+='{ '+'\n'
    for key,val in zip(keys,vals):
        if node_is_leaf_node(val):
            res+='  '+key+': '+str(val)+'\n'
    res+='  ,'+'\n'
    for key,val in zip(keys,vals):
        if node_is_admixture(val):
            res+='  '+key+': '+str(val)+'\n'
    res+='  ,'+'\n'
    for key,val in zip(keys,vals):
        if node_is_coalescence(val):
            res+='  '+key+': '+str(val)+'\n'
    res+='}'
    return res
    
    
def make_R_draw_from_files(drawing_name, file_names):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    cmd=['Rscript', dir_path+os.path.sep+'make_drawing.R', drawing_name]+file_names
    print(cmd)
    call(cmd)
