from math import floor
def make_readable_leaf_name(leaf_node):
    return list(leaf_node)[0]
    
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
    root_set=None
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
            string_names[key]=make_readable_leaf_name(key)
            pure_leaves.append(string_names[key])
        elif plotting_type=='admixture_leaf':
            string_names[key]=make_readable_leaf_name(key)
            admixture_leaves.append(string_names[key])
        elif plotting_type=='root':
            string_names[key]='r'
            root.append(string_names[key])
        else:
            assert False, 'unknown plotting type'+str(plotting_type)
        
    for key, node in list(node_structure.items()):
        parents=node.get_parents()
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
    
    def get_parents(self):
        return self.parents
    
    def get_children(self):
        return self.children
    
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
