class Node():
    
    def __init__(self, name):
        self.name=name
        self.parents=[]
        self.children=[]
        
    def add_child(self, child_node):
        self.children.append(child_node)
        
    def add_parent(self, parent_node):
        self.parents.append(parent_node)
        
    def set_children(self, children):
        self.children=children
        
    def set_parents(self, parents):
        self.parents=parents
        
    def clear_parents(self):
        self.parents=[]
        
    def clear_children(self):
        self.children=[]
    
    def get_parents(self):
        return self.parents
    
    def get_children(self):
        return self.children
    
    def get_number_of_children(self):
        return len(self.children)
        
    def get_number_of_parents(self):
        return len(self.parents)
    
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
    
    def has_children(self):
        return len(self.children)>0
    
    def has_parent(self):
        return len(self.parents)>0
    
    def is_admixture(self):
        return len(self.children)==1 and len(self.parents)==2
    
    def is_wrong(self):
        leaf= len(self.children)==0 and len(self.parents)==1
        admixture= len(self.children)==1 and len(self.parents)==2
        coalescence= len(self.children)==2 and len(self.parents)==1
        root= len(self.children)==2 and len(self.parents)==0
        situation_0= len(self.children)==0 and len(self.parents)>1
        situation_1= len(self.children)==1 and len(self.parents)>2
        situation_2= len(self.children)==2 and len(self.parents)>1
        assert leaf or admixture or coalescence or root or situation_0 or situation_1 or situation_2, 'Illegal tree structure detected'
        return not (root or leaf or admixture or coalescence)
    
    def is_too_wrong(self):
        return len(self.children)>2
    
    def replace_child(self, old_child, new_child):
        self.children.remove(old_child)
        self.children.append(new_child)
        
    def replace_parent(self, old_parent, new_parent):
        self.parents.remove(old_parent)
        self.parents.append(new_parent)
        
    def remove_parent(self, parent):
        self.parents.remove(parent)

def initor(a):
    if not isinstance(a, str):
        return a[0]
    else:
        return a
    