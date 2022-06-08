from Rtree_operations import find_rooted_nodes, get_real_children, get_real_parents, get_branch_length, mother_or_father, pretty_print
from piece_of_tree import piece


def topological_branch_length(*args, **kwargs):
    return 1
    

class lineage(object):
    
    def __init__(self, key,  distance=0, lattitude=0, topological_distance=False):
        self.key,self.distance, self.lattitude= key, distance, lattitude
        if topological_distance:
            self.get_branch_length=topological_branch_length
        else:
            self.get_branch_length=get_branch_length
        self.topological_distance=topological_distance
        
    def follow(self, tree, visited_keys=[]):
        node=tree[self.key]
        new_lineages=[]
        pieces=[]
        new_keys_visited=[]
        for key in get_real_children(tree[self.key]):
            if key not in visited_keys:
                branch=mother_or_father(tree, child_key=key, parent_key=self.key)
                l=get_branch_length(tree, key, branch)
                pieces.append(piece(self.lattitude, self.lattitude-l, self.distance, self.distance+l, key, branch, self.key))
                new_lineages.append(lineage(key,self.distance+l, self.lattitude-l, topological_distance=self.topological_distance))
        for key in get_real_parents(tree[self.key]):
            if key not in visited_keys:
                branch=mother_or_father(tree, child_key=self.key, parent_key=key)
                l=get_branch_length(tree, self.key, branch)
                pieces.append(piece(self.lattitude, self.lattitude+l, self.distance, self.distance+l, self.key, branch, key))
                new_lineages.append(lineage(key,self.distance+l, self.lattitude+l, topological_distance=self.topological_distance))
        if self.key=='r':#add the very long piece
            pieces.append(piece(self.lattitude, None, self.distance, None,'r',0,None))
        return new_lineages, pieces

    def under_cap(self, cap):
        return self.distance<cap


def insert_root_in_tree(tree):
    (child_key1,_,_),(child_key2,_,_)=find_rooted_nodes(tree)
    tree['r']=[None,None,None,None,None,child_key1, child_key2]

def remove_root_from_tree(tree):
    del tree['r']
    

def distanced_branch_lengths(tree, start_key, visited_keys=[], upper_cap=float('inf'), topological_distance=False):
    insert_root_in_tree(tree)
    pieces=[]
    lineages=[lineage(start_key, 0, 0, topological_distance=topological_distance)]
    while lineages:
        lineages.sort(key=lambda x: x.distance)
        lin=lineages.pop(0)
        if lin.key not in visited_keys:
            visited_keys.append(lin.key)
            #removedprin lin.key
            new_lineages, new_pieces= lin.follow(tree, visited_keys)
            pieces.extend(new_pieces)
            lineages.extend([new_lineage for new_lineage in new_lineages if new_lineage.under_cap(upper_cap)])
    remove_root_from_tree(tree)
    return pieces

    