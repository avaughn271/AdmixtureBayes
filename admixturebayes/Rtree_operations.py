from copy import deepcopy
import warnings

def create_trivial_tree(size, total_height=1.0):
    '''
    constructs tree of the form (..((s1,s2),s3),s4)...)
    '''
    step_size=total_height/size
    tree={'s1':['n1',None,None,step_size,None, None,None],
          's2':['n1',None,None,step_size,None, None,None],
          'n1':['n2',None,None,step_size,None, 's1','s2']}
    nex_inner_node='n2'
    new_inner_node='n1'
    for k in range(3,size+1):
        old_inner_node='n'+str(k-2)
        new_inner_node='n'+str(k-1)
        nex_inner_node='n'+str(k)
        new_leaf='s'+str(k)
        tree[new_leaf]=[new_inner_node, None,None, step_size*(k-1),None, None,None]
        tree[new_inner_node]=[nex_inner_node, None,None, step_size, None, new_leaf,old_inner_node]
    del tree[new_inner_node]
    return rename_root(tree, new_inner_node)

def rename_leaves(tree, new_leaf_names):
    old_keys=get_leaf_keys(tree)
    assert len(new_leaf_names) == len(old_keys), 'number of renamed nodes did not match the actual number of nodes'
    unique_identifier= max(old_keys, key=len)+max(new_leaf_names, key=len)+'_'
    temporary_keys=[ok+unique_identifier for ok in old_keys]
    for old_key, new_key in zip(old_keys, temporary_keys):
        tree=rename_key(tree, old_key, new_key )
    for old_key, new_key in zip(temporary_keys, new_leaf_names):
        tree = rename_key(tree, old_key, new_key)
    return tree

def max_distance_to_leaf(tree,key, parent_key=None):
    if key=='r':
        (child_key1,_,_),(child_key2,_,_)=find_rooted_nodes(tree)
        return max(max_distance_to_leaf(tree, child_key1, key), max_distance_to_leaf(tree, child_key2,key))
    node=tree[key]
    if parent_key is not None:
        branch=mother_or_father(tree, key, parent_key)
        add=node[3+branch]
    else:
        add=0
    if node_is_leaf_node(node):
        return add
    if node_is_coalescence(node):
        return add+max(max_distance_to_leaf(tree, node[5], key), max_distance_to_leaf(tree, node[6],key))
    if node_is_admixture(node):
        return add+max_distance_to_leaf(tree, node[5], key)
    assert False, 'strange node caused no exit.'
    
def time_adjust_node(key, node, timed_events):
    p1,p2=node[:2]
    t1=timed_events[p1]-timed_events[key]
    node[3]=t1
    if p2 is not None:
        t2=timed_events[p2]-timed_events[key]
        node[4]=t2
    return node

def get_max_timing(tree):
    timed_events={'r':get_max_distance_to_root(tree)}
    for key in tree:
        timed_events[key]=max_distance_to_leaf(tree, key)
    return timed_events    
    
def time_adjust_tree(tree):
    timed_events=get_max_timing(tree)
    for key, node in list(tree.items()):
        tree[key]=time_adjust_node(key, node, timed_events)
    return tree

def create_trivial_equibranched_tree(size, height=1.0):
    '''
    constructs tree of the form (..((s1,s2),s3),s4)...)
    '''
    tree={'s1':['n1',None,None,height,None, None,None],
          's2':['n1',None,None,height,None, None,None],
          'n1':['n2',None,None,height,None, 's1','s2']}
    nex_inner_node='n2'
    new_inner_node='n1'
    for k in range(3,size+1):
        old_inner_node='n'+str(k-2)
        new_inner_node='n'+str(k-1)
        nex_inner_node='n'+str(k  )
        new_leaf='s'+str(k)
        tree[new_leaf]=[new_inner_node, None,None, height,None, None,None]
        tree[new_inner_node]=[nex_inner_node, None,None, height, None, new_leaf,old_inner_node]
    del tree[new_inner_node]
    return rename_root(tree, new_inner_node)

def add_outgroup(tree, inner_node_name='new', to_new_root_length=0.5, to_outgroup_length=0.5, outgroup_name=None):
    (child_key1, child_branch1,_),(child_key2, child_branch2, _)=find_rooted_nodes(tree)
    tree[inner_node_name]=['r', None, None, to_new_root_length, None, child_key1, child_key2]
    tree[child_key1][child_branch1]=inner_node_name
    tree[child_key2][child_branch2]=inner_node_name
    if outgroup_name is None:
        n=get_number_of_leaves(tree)
        outgroup_name='s'+str(n+1)
    tree[outgroup_name]=['r', None, None,to_outgroup_length, None, None, None]
    return tree

def non_admixture_path(tree, key):
    if key=='r':
        return True
    if node_is_admixture(tree[key]):
        return False
    parent_key=tree[key][0]
    return non_admixture_path(tree, parent_key)

def get_first_admixture_meeting(tree, key):
    if key=='r':
        return None
    if node_is_admixture(tree[key]):
        return key
    parent_key=tree[key][0]
    return get_first_admixture_meeting(tree, parent_key)

def get_all_admixture_meetings(tree, key):
    if key=='r':
        return []
    parent_key = tree[key][0]
    parent2_key=tree[key][1]
    if node_is_admixture(tree[key]):
        left_admixtures=get_all_admixture_meetings(tree, parent_key)
        right_admixtures=get_all_admixture_meetings(tree, parent2_key)
        uniq=list(set(left_admixtures+right_admixtures))
        return uniq+[key]
    return get_all_admixture_meetings(tree, parent_key)


def get_branches_to_reverse(tree, key, so_far=None):
    if so_far is None:
        so_far=[]
    if is_root(key):
        (key1, branch1, length1),(key2,branch2,length2)=find_rooted_nodes(tree)
        if key1==so_far[-1][0]:
            so_far.append((key2,tree[key2][branch2+3],tree[key2][branch2]))
        else:
            so_far.append((key1,tree[key1][branch1+3],tree[key1][branch1]))
        return so_far
    else:
        so_far.append((key,tree[key][3], tree[key][0]))
        return get_branches_to_reverse(tree, tree[key][0], so_far)
    
            

def rename_rootname(tree,old_name, new_name):
    for key,node in list(tree.items()):
        if node[0]==old_name:
            tree[key][0]=new_name
        if node[1]==old_name:
            tree[key][1]=new_name 
    return tree

def remove_children(tree):
    for key in tree:
        tree[key]=tree[key][:5]
    return tree

def rearrange_root_foolproof(tree, new_outgroup):
    '''
    Like rearrange_root this changes the outgroup by reversing branches. When the 
    '''
    if tree[new_outgroup][0]=='r' or tree[new_outgroup][1]=='r':
        warnings.warn('The root was already in the requested location so no rearranging performed', UserWarning)
        return tree
    while True:
        admixture_to_remove=get_first_admixture_meeting(tree, new_outgroup)
        print('tryna remove', admixture_to_remove)
        if admixture_to_remove is None:
            break
        pretty_print(tree)
        tree=remove_admix(tree, admixture_to_remove, 1)[0]
    return rearrange_root(tree, new_outgroup)

def rearrange_root(tree, new_outgroup):
    if tree[new_outgroup][0]=='r':
        warnings.warn('The root was already in the requested location so no rearranging performed', UserWarning)
        return tree
    assert non_admixture_path(tree, new_outgroup), 'There were admixtures on the path from the requested outgroup to the old root.'
    reversers=get_branches_to_reverse(tree, new_outgroup)
    #removedprin reversers[:-2]
    for child_key,length,parent_key in reversers[:-2]:
        #removedprin 'child,length,parent',(child_key,length,parent_key)
        tree[parent_key][0]=child_key
        tree[parent_key][3]=length
    (before_root, length1, _),(after_root,length2,_)=reversers[-2:]
    if is_root(tree[after_root][0]):
        tree[after_root][0]=before_root
        tree[after_root][3+0]=length1+length2
    else:
        tree[after_root][1]=before_root
        tree[after_root][3+1]=length1+length2
    tree[new_outgroup][0]='r'
    tree[new_outgroup][0+3]=0.0
    first_reverser_parent=reversers[0][2]
    tree[first_reverser_parent][0]='r'
    tree=insert_children_in_tree(tree)
    return tree
    
def reverse_node(tree, key, old_child, new_parent_key=None):
    if new_parent_key is None:
        new_parent_key=old_child
    new_branch_length=get_branch_length_from_parent(tree, old_child, key)
    forgetten_branch_length=get_branch_length(tree, key, branch=0) #we know that it is a non-admixture_node
    
    
    
def prune_double_nodes(tree):
    '''
    A double node is a node with two parents and two children
    '''    
    for key, node in list(tree.items()):
        if len(get_real_children(node))!=1 and len(get_real_parents(node))==2:
            tree=split_up_double_node(tree, key)
    return tree

def split_up_double_node(tree, key):
    old_node=tree[key]
    new_node=key+'pruned'
    for parent in get_real_parents(old_node):
        tree[parent]=_rename_child(tree[parent], old_name=key, new_name=new_node)
    tree[new_node]=deepcopy(old_node)
    tree[new_node][5]=key
    tree[new_node][6]=None
    tree[key]=[new_node,None,None,0,None,old_node[5], old_node[6]]
    return tree

def rename_key(tree, old_key_name, new_key_name):
    node=tree[old_key_name]
    tree[new_key_name]=node
    ps= get_real_parents(node)
    for p in ps:
        if not is_root(p):
            tree[p]=_rename_child(tree[p], old_key_name, new_key_name)
    cs=get_real_children(node)
    for c in cs:
        tree[c]=rename_parent(tree[c], old_key_name, new_key_name)
    del tree[old_key_name]
    return tree

def remove_outgroup(tree, remove_key='s1', return_add_distance=False):
    (child_key1, child_branch1,length1),(child_key2, child_branch2, length2)=find_rooted_nodes(tree)
    if remove_key==child_key1:
        root_key=child_key2
    elif remove_key== child_key2:
        root_key= child_key1
    else:
        assert remove_key==child_key1 or remove_key==child_key2, 'the removed key is not an outgroup'
    del tree[remove_key]
    del tree[root_key]
    tree= rename_root(tree, root_key)
    if return_add_distance:
        return tree, length1+length2
    return tree


def simple_reorder_the_leaves_after_removal_of_s1(tree):
    no_leaves=get_number_of_leaves(tree)
    for n in range(no_leaves):
        tree=rename_key(tree, 's'+str(n+2), 's'+str(n+1))
    return tree
        


def find_children(tree, parent_key):
    res=[]
    for key,node in list(tree.items()):
        if node[0]==parent_key:
            res.append(key)
        if node[1]==parent_key:
            res.append(key)
    while len(res)<2:
        res.append(None)
    return res

def create_balanced_tree(size, height=1.0):
    return finish_tree_with_coalescences({}, get_trivial_nodes(size), height)
        

def finish_tree_with_coalescences(tree, keys_to_finish, height=1.0):
    while len(keys_to_finish)>1:
        key1,key2=keys_to_finish[:2]
        if len(keys_to_finish)==2:
            new_key='r'
        else:
            new_key=key1+key2
        tree[key1]=[new_key,None,None,height,None]+find_children(tree,key1)
        tree[key2]=[new_key,None,None,height,None]+find_children(tree,key2)
        keys_to_finish=keys_to_finish[2:]
        keys_to_finish.append(new_key)
        #removedprin len(keys_to_finish)
    return tree
        

def create_burled_leaved_tree(size, height):
    res={}
    for i in range(size):
        res['s'+str(i+1)]=['a'+str(i+1), None, None, height, None,None,None]
        res['a'+str(i+1)]=['b'+str(i+1), 'm'+str(i+1), 0.5, height, height,'s'+str(i+1),None]
        res['b'+str(i+1)]=['n'+str(i+1), 'm'+str(i+1), 0.5, height, height,'a'+str(i+1),None]
        res['m'+str(i+1)]=['n'+str(i+1),None ,None, height, None,'a'+str(i+1),'b'+str(i+1)]
    return finish_tree_with_coalescences(res, ['n'+str(i+1) for i in range(size)], height)
    

def get_trivial_nodes(size):
    return ['s'+str(n+1) for n in range(size)]

def get_distance_to_root(tree, key, function=max):
    if key=='r':
        return 0.0
    node=tree[key]
    if node_is_admixture(node):
        return function(get_distance_to_root(tree, node[0], function=function)+node[3], 
                   get_distance_to_root(tree, node[1], function=function)+node[4], node[2])
    else:
        return get_distance_to_root(tree, node[0], function=function)+node[3]

def special_max(x,y,z=None):
    return max(x,y)

def special_min(x,y,z=None):
    return min(x,y)

def average_admixture_node(x,y,z):
    return z*x+(1-z)*y

def get_max_distance_to_root(tree):
    return max(get_leaf_distances_to_root(tree, function=special_max))

def get_min_distance_to_root(tree):
    return min(get_leaf_distances_to_root(tree, function=special_min))

def get_average_distance_to_root(tree):
    av_lengths=get_leaf_distances_to_root(tree, function=average_admixture_node)
    return float(sum(av_lengths))/len(av_lengths)

def get_leaf_distances_to_root(tree, function=max):
    res=[]
    for key, node in list(tree.items()):
        if node_is_leaf_node(node):
            res.append(get_distance_to_root(tree, key, function=function))
    return res

def get_number_of_ghost_populations(tree):
    count=0
    for key,node in list(tree.items()):
        if node_is_coalescence(node):
            child1,child2=get_children(node)
            if node_is_admixture(tree[child1]) and node_is_admixture(tree[child2]):
                count+=1
    return count

    
def update_all_branches(tree, updater):
    for key, node in list(tree.items()):
        if node_is_admixture(node):
            node[2]+=updater()
            node[3]+=updater()
            node[4]+=updater()
            if node[2]<0 or node[2]>1 or node[3]<0 or node[4]<0:
                return None
        else:
            node[3]+=updater()
            if node[3]<0:
                return None        
    return tree

def update_all_admixtures(tree, updater):
    for key, node in list(tree.items()):
        if node_is_admixture(node):
            node[2]+=updater()
            if node[2]<0 or node[2]>1:
                return None  
    return tree

def update_node(tree, key, updater, admixture_proportion_multiplier=1.0):
    node=tree[key]
    if node_is_admixture(node):
        node[2]+=updater()*admixture_proportion_multiplier
        node[3]+=updater()
        node[4]+=updater()
        if node[2]<0 or node[2]>1 or node[3]<0 or node[4]<0:
            return None
    else:
        node[3]+=updater()
        if node[3]<0:
            return None 
    tree[key]=node       
    return tree

def update_branch(tree, key, branch, updater):
    tree[key][branch+3]+= updater()
    return tree


def update_branch_length(tree,key,branch, new_length):
    tree[key][branch+3]=new_length

def extend_branch(node, pkey, grand_parent_key, p_to_gp):
    #removedprin node, pkey, grand_parent_key, p_to_gp
    if node[0]==pkey:
        node[0]=grand_parent_key
        u=node[3]/(node[3]+p_to_gp)
        node[3]+=p_to_gp
        return node,u,node[3]
    elif node[1]==pkey:
        node[1]=grand_parent_key
        u=node[4]/(node[4]+p_to_gp)
        node[4]+=p_to_gp
        return node,u,node[4]
    else:
        assert False, 'extension of branch was not possible'

    
def remove_parent_attachment(tree, orphanota_key, orphanota_branch):
    '''
    This takes the tree and removes the parent of orphanonte_key.
    '''
    pkey=get_parent_of_branch(tree, orphanota_key, orphanota_branch)
    
    if pkey=='r':
        return remove_root_attachment(tree, orphanota_key, orphanota_branch)
    grand_pkey=get_parents(tree[pkey])[0]
    child_of_parent=get_other_children(tree[pkey], orphanota_key)[0]
    sib_node=tree[child_of_parent]
    tree[child_of_parent],u,extended_branch_length=extend_branch(sib_node, pkey, grand_pkey, tree[pkey][3])
    del tree[pkey]
    if grand_pkey!='r':
        tree[grand_pkey]=_rename_child(tree[grand_pkey], pkey, child_of_parent)
    tree[orphanota_key][orphanota_branch]=None
    return tree,"u",u*extended_branch_length,extended_branch_length

def remove_root_attachment(tree, orphanota_key, orphanota_branch):
    '''
    The situation is different when the root is removed because of the special naming strategy.
    
                r
               / \
             /    \
      orphanota   new_root 
    Here a new root is born.
    
    '''
    rooted_keys=find_rooted_nodes(tree)
    for key,branch,len_to_root in rooted_keys:
        if key!=orphanota_key or orphanota_branch!=branch:
            if node_is_coalescence(tree[key]):
                tree=rename_root(tree, key)
                r=len_to_root
                del tree[key]
            else:
                tree[key],r=get_branch_length_and_reset(tree[key], 'r', 'closed_branch')
                #removedprin 'closed_branch!'
            tree[orphanota_key][orphanota_branch]=None
    return tree,'r', r,None
    
def get_branch_length_and_reset(node, parent_key, new_length,add=False):
    if node[0]==parent_key:
        old_length=node[3]
        if add:
            node[3]+=new_length
        else:
            node[3]=new_length
        return node, old_length
    elif node[1]==parent_key:
        old_length=node[4]
        if add:
            node[4]+=new_length
        else:
            node[4]=new_length
        return node, old_length
    else:
        assert False, 'could not give new length because the parent was unknown'
    
def rename_root(tree, old_name):
    for _,node in list(tree.items()):
        if node[0]==old_name:
            node[0]='r'
        if (node[1] is not None and node[1]==old_name):
            node[1]='r'
    return tree

def _rename_child(node, old_name, new_name):
    if node[5]==old_name:
        node[5]=new_name
    elif node[6]==old_name:
        node[6]=new_name
    else:
        assert False, "tried to rename a child that didn't exist in its parents documents."
    return node

def rename_parent(node, old_name, new_name):
    if node[0]==old_name:
        node[0]=new_name
    elif node[1]==old_name:
        node[1]=new_name
    else:
        assert False, 'tried to rename a parent that didnt exist in its childs document'
    return node
    
def find_rooted_nodes(tree):
    res=[]
    for key,node in list(tree.items()):
        if node[0]=='r' or (node[1] is not None and node[1]=='r'):
            if node[0]=='r':
                res.append((key,0,node[3]))
            else:
                res.append((key,1,node[4]))
    return res

def _find_rooted_branches(tree):
    res=[]
    for key,node in list(tree.items()):
        if node[0]=='r':
            res.append((key,0))
        if node[1] is not None and node[1]=='r':
            res.append((key,1))
    return res

def _get_root_sibling(tree, child_key, child_branch):
    root_keys=_find_rooted_branches(tree)
    if len(root_keys)==1:
        return child_key, child_branch
    else:
        if root_keys[0][0]==child_key and root_keys[0][1]==child_branch:
            return root_keys[1]
        else:
            return root_keys[0]

def node_is_non_admixture(node):
    return (node[1] is None)

def node_is_admixture(node):
    return (node[1] is not None)

def node_is_coalescence(node):
    return (node[1] is None and node[5] is not None)

def node_is_leaf_node(node):
    return (node[1] is None and node[5] is None)



def get_descendants_and_rest(tree, key):
    all_keys=list(tree.keys())
    #removedprin tree, key
    descendant_keys=_get_descendants(tree, key)
    return descendant_keys, list(set(all_keys)-set(descendant_keys))

def _get_descendants(tree, key):
    if tree[key][5] is None:
        return [key]
    else:
        ans=[key]+_get_descendants(tree, tree[key][5])
        if tree[key][6] is None:
            return ans
        return ans+_get_descendants(tree, tree[key][6])

def get_all_branches(tree):
    res=[]
    for key, node in list(tree.items()):
        if node_is_admixture(node):
            res.extend([(key, 0),(key,1)])
        else:
            res.append((key,0))
    return res

def get_specific_branch_lengths(tree, branches):
    res=[]
    for key, branch in branches:
        res.append(tree[key][branch+3])
    return res

def update_specific_branch_lengths(tree, branches, new_lengths, add=False):
    if add:
        for (key, branch), new_length in zip(branches,new_lengths):
            tree[key][branch+3]+=new_length
            if tree[key][branch+3]<0:
                return None
    else:
        for (key, branch), new_length in zip(branches,new_lengths):
            tree[key][branch+3]=new_length
            if tree[key][branch+3]<0:
                return None
    return tree


        

def get_all_branch_descendants_and_rest(tree, key,branch):
    all_branches=get_all_branches(tree)
    descendant_branches=_get_descendant_branches(tree, key,branch)
    return descendant_branches, list(set(all_branches)-set(descendant_branches))

def _get_descendant_branches(tree, key, branch):
    child_key1=tree[key][5]
    if child_key1 is None: #branch is leaf
        return [(key,0)]
    else:
        child_key1_branch=mother_or_father(tree, child_key1, key)
        ans=[(key, branch)]+_get_descendant_branches(tree, child_key1, child_key1_branch)
        child_key2=tree[key][6]
        if child_key2 is None:
            return ans
        return ans+_get_descendant_branches(tree, child_key2, mother_or_father(tree,child_key2, key))

def get_other_children(node, child_key):
    res=[]
    for n in get_children(node):
        if n is not None and n!=child_key:
            res.append(n)
    return res

def get_sibling_on_parent_side(tree, parent_key, child_key):
    if parent_key=='r':
        (child_key1, _,_), (child_key2,_,_) = find_rooted_nodes(tree)
        if child_key==child_key1:
            return child_key2
        elif child_key==child_key2:
            return child_key1
        else:
            assert False, 'Claimed child of root was not a child of root'
    else:
        return get_other_children(tree[parent_key], child_key)

def get_children(node):
    return node[5:7]

def _get_index_of_parent(node, parent):
    if node[0]==parent:
        return 0
    if node[1]==parent:
        return 1
    return -1

def has_child_admixture(tree, key):
    node=tree[key]
    child1,child2=get_children(node)
    if child1 is not None:
        if node_is_admixture(tree[child1]):
            return True
    if child2 is not None:
        return (node_is_admixture(tree[child2]))
    return False

def screen_and_prune_one_in_one_out(tree):
    keys_to_check=list(tree.keys())
    for key in keys_to_check:
        if key in tree:
            node=tree[key]
            children=get_real_children(node)
            parents=get_real_parents(node)
            if len(children)==1 and len(parents)==1:
                tree=remove_one_in_one_out(tree, children[0], key, parents[0])
    return tree
        
def remove_one_in_one_out(tree, child, key, parent_key):
    branch=mother_or_father(tree, key, parent_key)
    length=get_branch_length(tree, key, branch)
    child_branch=mother_or_father(tree, child, key)
    child_length=get_branch_length(tree, child, child_branch)
    tree=update_parent_and_branch_length(tree, child, child_branch, new_parent=parent_key, new_branch_length=child_length+length)
    if parent_key!='r':
        tree[parent_key]=_update_child(tree[parent_key], key, child)
    del tree[key]
    return tree
    
    
def get_leaf_keys(tree):
    res=[]
    for key, node in list(tree.items()):
        if node_is_leaf_node(node):
            res.append(key)
    return res

def get_no_leaves(tree):
    return len(get_leaf_keys(tree))

def get_categories(tree):
    leaves=[]
    admixture_nodes=[]
    coalescence_nodes=[]
    for key, node in list(tree.items()):
        if node_is_leaf_node(node):
            leaves.append(key)
        if node_is_coalescence(node):
            coalescence_nodes.append(key)
        if node_is_admixture(node):
            admixture_nodes.append(key)
    return leaves, coalescence_nodes, admixture_nodes

def get_parent_of_branch(tree, key, branch):
    assert key!='r', 'Tried to access the parent of the root branch'
    return tree[key][branch]

def get_branch_length(tree,key,branch):
    assert key!='r', 'Tried to access the length of the root branch'
    return tree[key][branch+3]


def get_branch_length_from_parent(tree, child_key, parent_key):
    return tree[child_key][3+mother_or_father(tree, child_key, parent_key)]

def get_admixture_proportion_from_key(tree, key):
    return tree[key][2]

def get_admixture_proportion(tree, child_key,child_branch):
    key=get_parent_of_branch(tree, child_key,child_branch)
    assert node_is_admixture(tree[key]), 'Tried to get the admixture proportion of a non-admixture node'
    return tree[key][2]

def update_parent_and_branch_length(tree, child_key, child_branch, new_parent, new_branch_length):
    assert child_key!='r', 'Tried to update the root branch'
    tree[child_key][child_branch]=new_parent
    tree[child_key][child_branch+3]=new_branch_length
    return tree

def get_admixture_keys_and_proportions(tree):
    keys=[]
    props=[]
    for key, node in list(tree.items()):
        if node_is_admixture(node):
            keys.append(key)
            props.append(node[2])
    return keys, props

def get_pruned_tree_and_add(tree, outgroup):
    assert tree[outgroup][0]=='r', 'can not remove outgroup which is not outgroup'
    (child1, branch1, a), (child2,branch2, b)=find_rooted_nodes(tree)
    del tree[child2]
    del tree[child1]
    if child1==outgroup:
        tree=rename_root(tree, child2)
    else:
        tree=rename_root(tree, child1)
    return tree, a+b
        
            
    


def get_destination_of_lineages(tree, ready_lineages):
    single_coalescences={} #list of tuples ((key,branch),(sister_key,sister_branch))
    double_coalescences=[]
    admixtures=[]
    for key, branch in ready_lineages:
        if (key,branch) in single_coalescences:
            double_coalescences.append(((key,branch),single_coalescences[(key,branch)]))
            del single_coalescences[(key,branch)]
            continue
        parent_key=tree[key][branch]
        if parent_key=='r':
            sister_key, sister_branch=_get_root_sibling(tree, key, branch)
            single_coalescences[(sister_key,sister_branch)]=(key,branch)
            continue
        parent=tree[parent_key]
        if node_is_coalescence(parent):
            sister_key, sister_branch=get_sister_branch(tree, parent,key, branch)
            single_coalescences[(sister_key,sister_branch)]=(key,branch)
        elif node_is_admixture(parent):
            admixtures.append((key,branch))
        else:
            assert False, 'the parent of a node was neither admixture nor coalescence'
    return double_coalescences, single_coalescences, admixtures

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

def pretty_print(tree):
    print(pretty_string(tree))

def get_sister_branch(tree, parent, key, branch):
    #removedprin parent,key
    #removedprin pretty_string(tree)
    if parent[5]==parent[6]:
        return key, other_branch(branch)
    else:
        if parent[5]==key:
            return parent[6], mother_or_father(tree, parent[6], tree[key][branch])
        elif parent[6]==key:
            return parent[5], mother_or_father(tree, parent[5], tree[key][branch])
        else:
            assert False, "the parent was not really a parent"+'\n'+pretty_string(tree)+'\n'+'parent,key,branch='+str(parent)+','+str(key)+','+str(branch)
        
        

def propagate_married(tree, list_of_pairs):
    res=[]
    for (key1,branch1),(_,_) in list_of_pairs:
        parent_key=get_parent_of_branch(tree, key1, branch1)
        res.append((parent_key,0))
    return res

def propagate_admixtures(tree, list_of_admixtures):
    res=[]
    for key,branch in list_of_admixtures:
        parent_key=get_parent_of_branch(tree, key, branch)
        res.append((parent_key,0))
        res.append((parent_key,1))
    return res
    
        
def mother_or_father(tree, child_key, parent_key):
    if tree[child_key][0]==parent_key:
        return 0
    elif tree[child_key][1]==parent_key:
        return 1
    assert False, 'The child did not match its parent'+\
                  '\n'+pretty_string(tree)+'\n\n'+\
                  'child_key,parent_key='+str(child_key)+\
                  ','+str(parent_key)
    
def insert_children_in_tree(tree):
    children={key:[] for key in tree}
    for key in tree:
        parents = get_real_parents(tree[key])
        for parent in parents:
            if parent!='r':
                children[parent].append(key)
    #removedprin children
    for key in tree:
        tree[key]=_update_parents(tree[key], children[key])
    return tree



def get_all_branch_lengths(tree):
    res=[]
    for key, node in list(tree.items()):
        if node_is_non_admixture(node):
            res.append(node[3])
        else:
            res.extend(node[3:5])
    return res
    
def get_all_admixture_proportions(tree):
    res=[]
    for key, node in list(tree.items()):
        if node_is_admixture(node):
            res.append(node[2])
    return res

def convert_to_vector(tree, keys=None):
    if keys is None:
        keys=list(tree.keys())
    res=[]
    for key in keys:
        node=tree[key]
        if node_is_non_admixture(node):
            res.append(node[3])
        else:
            res.extend(node[3:5])
            res.append(node[2])
    return res

def graft(tree, remove_key, add_to_branch, insertion_spot, new_node_code, which_branch, remove_branch=0):
    #if we are at the root, things are different.
    if add_to_branch=='r':
        return graft_onto_root(tree, insertion_spot, remove_key, new_name_for_old_root=new_node_code, remove_branch=remove_branch)
    
    #updating the grafted_branch. Easy.
    tree[remove_key][remove_branch]=new_node_code
    
    #dealing with the other child of the new node
    index_of_pop_name_to_change=which_branch
    index_of_branchl_to_change=which_branch+3
    #saving info:
    parent_of_branch=tree[add_to_branch][index_of_pop_name_to_change]
    length_of_piece_to_break_up=tree[add_to_branch][index_of_branchl_to_change]
    #updating node
    tree[add_to_branch][index_of_pop_name_to_change]=new_node_code
    tree[add_to_branch][index_of_branchl_to_change]=length_of_piece_to_break_up*insertion_spot
    
    #taking care of the new node
    tree[new_node_code]=[parent_of_branch, None, None, length_of_piece_to_break_up*(1.0-insertion_spot), None, add_to_branch, remove_key]

    #taking care of the parent of the new node
    if parent_of_branch != 'r':
        tree[parent_of_branch]=_update_child(tree[parent_of_branch], add_to_branch, new_node_code)
        
    return tree
    
def graft_onto_root(tree, insertion_spot, remove_key, new_name_for_old_root, remove_branch=0):    
    root_keys=find_rooted_nodes(tree)

    if len(root_keys)==1:#this is the special case where an admixture leads to the new root
        return graft_onto_rooted_admixture(tree, insertion_spot, remove_key, root_keys[0], remove_branch=remove_branch)
    
    #updating the grafted_branch. Easy.
    tree[remove_key][remove_branch]='r'
    
    #dealing with the other child of the new node, but since the new node is the root, the old root is the new node. If that makes sense.
    #removedprin 'root_keys', root_keys
    tree[new_name_for_old_root]=['r', None, None, insertion_spot, None,root_keys[0][0], root_keys[1][0]]
    
    #dealing with the children of the new node.
    for rkey, _, b_l in root_keys:
        tree[rkey]=_update_parent(tree[rkey], 'r', new_name_for_old_root)
    
    return tree

def materialize_root(tree, new_key):
    (child_key1, child_branch1,_),(child_key2, child_branch2,_) = find_rooted_nodes(tree)
    tree[new_key]=['r',None,None,None,None,child_key1, child_key2]
    tree[child_key1][child_branch1]=new_key
    tree[child_key2][child_branch2]=new_key
    return tree

def move_node(tree, regraft_key, regraft_branch, parent_key, distance_to_regraft, chosen_piece, new_node_name='x'):
    if parent_key=='r':
        sister_key,sister_branch= _get_root_sibling(tree, regraft_key, regraft_branch)
        if chosen_piece.child_key=='r':
            u1,_=chosen_piece.get_leaf_and_root_sided_length(distance_to_regraft)
            tree[sister_key][sister_branch+3]=tree[sister_key][sister_branch+3]+u1
        elif chosen_piece.child_key==sister_key and chosen_piece.child_branch==sister_branch:
            u1,_=chosen_piece.get_leaf_and_root_sided_length(distance_to_regraft)
            tree[sister_key][sister_branch+3]=u1
        else:
            tree=rename_root(tree, sister_key)
            del tree[sister_key]
            u1,u2 =chosen_piece.get_leaf_and_root_sided_length(distance_to_regraft)
            if chosen_piece.parent_key == sister_key:
                tree[new_node_name]=['r', None,None,u2,None,regraft_key,chosen_piece.child_key]
            else:
                tree[new_node_name]=[chosen_piece.parent_key, None,None,u2,None,regraft_key,chosen_piece.child_key]
            if chosen_piece.parent_key=='r':
                pass
            elif chosen_piece.parent_key != sister_key:
                tree[chosen_piece.parent_key]=_rename_child(tree[chosen_piece.parent_key], chosen_piece.child_key, new_node_name)
            tree=update_parent_and_branch_length(tree, chosen_piece.child_key, chosen_piece.child_branch, new_node_name, u1)
            tree[regraft_key][regraft_branch]=new_node_name
    else:
        if chosen_piece.child_key=='r':
            tree=materialize_root(tree, new_node_name)
            u1,_=chosen_piece.get_leaf_and_root_sided_length(distance_to_regraft)
            tree[new_node_name][3]=u1       
            tree[regraft_key][regraft_branch]='r'     
        else:
            u1,u2=chosen_piece.get_leaf_and_root_sided_length(distance_to_regraft)
            #removedprin u1,u2
            tree[new_node_name]=[chosen_piece.parent_key, None,None,u2,None,regraft_key,chosen_piece.child_key]
            if chosen_piece.parent_key=='r':
                pass
            else:
                tree[chosen_piece.parent_key]=_rename_child(tree[chosen_piece.parent_key], chosen_piece.child_key, new_node_name)
            tree=update_parent_and_branch_length(tree, chosen_piece.child_key, chosen_piece.child_branch, new_node_name, u1)
            tree[regraft_key][regraft_branch]=new_node_name
        sibling_key=get_other_children(tree[parent_key], regraft_key)[0]
        grandfather_key=tree[parent_key][0]
        _=get_branch_length_and_reset(tree[sibling_key], parent_key, tree[parent_key][3], add=True)
        tree[sibling_key]=rename_parent(tree[sibling_key], parent_key, grandfather_key)
        if grandfather_key!='r':
            tree[grandfather_key]=_rename_child(tree[grandfather_key], parent_key, sibling_key)
        del tree[parent_key]
    
    return tree
        
    

def graft_onto_rooted_admixture(tree, insertion_spot, remove_key, root_key, remove_branch=0):
    #removedprin 'undoing a closed branch', insertion_spot, remove_key, root_key
    tree[remove_key][remove_branch]='r'
    tree[root_key[0]],_=get_branch_length_and_reset(tree[root_key[0]], 'r', insertion_spot)
    return tree

def change_admixture(node):
    assert node_is_admixture(node), 'tried to change admixture branches for a node without admixture'
    new_node=[node[1],node[0],1.0-node[2],node[4],node[3]]+node[5:]
    return new_node

def readjust_length(node):
    multip= 1.0/((1.0-node[2])**2+node[2]**2)
    node[3]=node[3]*multip
    return node, multip

def get_number_of_admixes(tree):
    return sum((1 for node in list(tree.values()) if node_is_admixture(node)))

def get_number_of_leaves(tree):
    return sum((1 for node in list(tree.values()) if node_is_leaf_node(node)))
            
def remove_admix(tree, rkey, rbranch):
    '''
    removes an admixture. besides the smaller tree, also the disappeared branch lengths are returned.
    parent_key          sparent_key
        |                |
        |t_1             | t_4
        |   __---- source_key
      rkey/   t_5       \
        |                \t_3
        |t_2          sorphanota_key  
    orphanota_key   
    and alpha=admixture proportion. The function returns new_tree, (t1,t2,t3,t4,t5), alpha 

    Note that t_4 could be None if source key is root. The source_key node is not an admixture node by assumption.
    '''
    rnode= tree[rkey]
    orphanota_key= get_children(rnode)[0]
    parent_key= rnode[other_branch(rbranch)]
    source_key= rnode[rbranch]
    t1= rnode[3+other_branch(rbranch)]
    t5= rnode[3+rbranch]
    alpha=rnode[2]
    
    tree[orphanota_key],t2=get_branch_length_and_reset(tree[orphanota_key], rkey, t1, add=True)
    tree[orphanota_key]=_update_parent(tree[orphanota_key], rkey, parent_key)
    orphanota_branch=_get_index_of_parent(tree[orphanota_key], parent_key)
    
    if parent_key!='r':
        tree[parent_key]=_update_child(tree[parent_key], old_child=rkey, new_child=orphanota_key)
    if source_key=='r':
        tree,_,t3,_=remove_root_attachment(tree, rkey, other_branch(rbranch)) #now sorphanota_key is new root
        del tree[rkey]
        return tree, (t1,t2,t3,None,t5),alpha,(orphanota_key,orphanota_branch)
    del tree[rkey]
    source_node=tree[source_key]
    sorphanota_key=get_other_children(source_node, child_key=rkey)[0]
    sparent_key=source_node[0]
    t4=source_node[3]
    if sparent_key!='r':
        tree[sparent_key]=_update_child(tree[sparent_key], source_key, sorphanota_key)
    tree[sorphanota_key],t3=get_branch_length_and_reset(tree[sorphanota_key], source_key, t4, add=True)
    tree[sorphanota_key]=_update_parent(tree[sorphanota_key], source_key, sparent_key)
    del tree[source_key]
    return tree, (t1,t2,t3,t4,t5), alpha, (orphanota_key,orphanota_branch)

def remove_admix2(tree, rkey, rbranch, pks={}):
    '''
    removes an admixture. besides the smaller tree, also the disappeared branch lengths are returned.
    parent_key          sparent_key
        |                |
        |t_1             | t_4
        |   __---- source_key
      rkey/   t_5       \
        |                \t_3
        |t_2          sorphanota_key  
    orphanota_key   
    and alpha=admixture proportion. The function returns new_tree, (t1,t2,t3,t4,t5), alpha 

    Note that t_4 could be None if source key is root. The source_key node is not an admixture node by assumption.
    '''
    rnode= tree[rkey]
    orphanota_key= get_children(rnode)[0]
    parent_key= rnode[other_branch(rbranch)]
    source_key= rnode[rbranch]
    t1= rnode[3+other_branch(rbranch)]
    t5= rnode[3+rbranch]
    if rbranch == 0:
        alpha=1-rnode[2]
    else:
        alpha=rnode[2]
        
    tree[orphanota_key],t2=get_branch_length_and_reset(tree[orphanota_key], rkey, t1, add=True)
    tree[orphanota_key]=_update_parent(tree[orphanota_key], rkey, parent_key)
    orphanota_branch=_get_index_of_parent(tree[orphanota_key], parent_key)
    pks['orphanota_key']=orphanota_key
    pks['orphanota_branch']=orphanota_branch
    
    if parent_key!='r':
        tree[parent_key]=_update_child(tree[parent_key], old_child=rkey, new_child=orphanota_key)
    if source_key=='r':
        tree,_,t3,_=remove_root_attachment(tree, rkey, 0) #now sorphanota_key is new root
        del tree[rkey]
        pks['sorphanota_key']='r'
        pks['sorphanota_branch']=0
        return tree, (t1,t2,t3,None,t5),alpha
    del tree[rkey]
    source_node=tree[source_key]
    sorphanota_key=get_other_children(source_node, child_key=rkey)[0]
    sparent_key=source_node[0]
    t4=source_node[3]
    if sparent_key!='r':
        tree[sparent_key]=_update_child(tree[sparent_key], source_key, sorphanota_key)
    tree[sorphanota_key],t3=get_branch_length_and_reset(tree[sorphanota_key], source_key, t4, add=True)
    tree[sorphanota_key]=_update_parent(tree[sorphanota_key], source_key, sparent_key)
    pks['sorphanota_key']=sorphanota_key
    pks['sorphanota_branch']=mother_or_father(tree, sorphanota_key, sparent_key)
    del tree[source_key]
    return tree, (t1,t2,t3,t4,t5), alpha

def tree_to_0tree(tree):
    leaves,_,admixture_keys=get_categories(tree)
    pruned_tree = deepcopy(tree)
    for adm_key in admixture_keys:
        if adm_key in pruned_tree:
            #removedprin '------------------------------------------'
            #removedprin 'removing', (adm_key, int_bin) , 'from tree:'
            pruned_tree=remove_admixture(pruned_tree, adm_key, 1)
    return pruned_tree

def direct_all_admixtures(tree, smaller_than_half=True):
    for key, node in list(tree.items()):
        if node_is_admixture(node):
            alpha=get_admixture_proportion_from_key(tree, key)
            if (smaller_than_half and alpha<0.5) or (not smaller_than_half and alpha>0.5):
                tree[key]=change_admixture(node)
    return tree
    
def other_branch(branch):
    if branch==0:
        return 1
    elif branch==1:
        return 0
    else:
        assert False, 'illegal branch'
    

def _update_parents(node, new_parents):
    if len(new_parents)==1:
        res=node[:5]+[new_parents[0],None]
        return res
    if len(new_parents)==2:
        res=node[:5]+new_parents
        return res
    if len(new_parents)==0:
        res=node[:5]+[None]*2
        return res
    assert False, 'how many parents do you think you have?'
    
def _update_parent(node, old_parent, new_parent):
    if node[0]==old_parent:
        node[0]=new_parent
    elif node[1]==old_parent:
        node[1]=new_parent
    else:
        assert False, 'parent could not be updated'
    return node
        
def _update_child(node, old_child, new_child):
    if node[5]==old_child:
        node[5]=new_child
    elif node[6]==old_child:
        node[6]=new_child
    else:
        assert False, 'child could not be updated'
    return node

def remove_non_mixing_admixtures(tree, limit=1e-7):
    '''
    Trees can have admixture events with admixture proportion very close to 0 or 1. This will remove those admixture proportions.
    '''
    admixtures_to_remove=[]
    for key, node in list(tree.items()):
        if node_is_admixture(node):
            if node[2]<limit:
                admixtures_to_remove.append((key,0))
            elif node[2]>1.0-limit:
                admixtures_to_remove.append((key,1))
    for adm_key, adm_branch in admixtures_to_remove:
        if adm_key in tree:#it could have been removed by others
            tree=remove_admixture(tree, adm_key, adm_branch)
    return tree
        
    
        
def get_parents(node):
    return node[:2]

def insert_admixture_node_halfly(tree, source_key, source_branch, insertion_spot, admix_b_length, new_node_name, admixture_proportion=0.51):
    '''
    Source should be sink, I guess in hindsight.
    '''
    node=tree[source_key]
    old_parent=node[source_branch]
    old_branch_length=node[source_branch+3]
    node[source_branch]=new_node_name
    node[source_branch+3]=old_branch_length*insertion_spot
    tree[new_node_name]=[old_parent, None, admixture_proportion, old_branch_length*(1-insertion_spot), admix_b_length, source_key, None]
    if old_parent != 'r':
        tree[old_parent]=_update_child(tree[old_parent], old_child=source_key, new_child=new_node_name)
    return tree

def get_real_parents(node):
    ps=node[:2]
    return [p for p in ps if p is not None]

def get_real_children(node):
    cs=node[5:7]
    return [c for c in cs if c is not None]

def get_keys_and_branches_from_children(tree, key):
    '''
    the key has to be an admixture
    '''
    
    child_keys=get_real_children(tree[key])
    branches=[]
    for child_key in child_keys:
        branches.append(mother_or_father(tree, child_key, key))
    return list(zip(child_keys, branches))
    
    

def get_real_children_root(tree, key):
    if key=='r':
        a,b=find_rooted_nodes(tree)
        return [(a[0],a[1]), (b[0],b[1])]
    else:
        c_keys=get_real_children(tree[key])
        res=[]
        for c_key in c_keys:
            res.append((c_key, mother_or_father(tree, c_key, key)))
        return res
        

def get_other_parent(node, parent_key):
    if parent_key==node[0]:
        return node[1]
    elif parent_key==node[1]:
        return node[0]
    else:
        assert False, 'the shared parent was not a parent of the sibling.'
        
def halfbrother_is_uncle(tree, key, parent_key):
    '''
    parent_key is a non admixture node. This function checks if a sibling is an admixture that goes to the parent and the grandparent at the same time.
    If removed, there would be a loop where one person has two of the same parent.
    '''
    sibling_key=get_other_children(tree[parent_key], key)[0]
    if node_is_non_admixture(tree[sibling_key]):
        return False
    bonus_parent=get_other_parent(tree[sibling_key], parent_key)
    grand_parent_key=tree[parent_key][0]
    return bonus_parent==grand_parent_key

def remove_admixture(tree, key, branch):
    parent_key=tree[key][branch]
    
    while parent_key!='r' and node_is_admixture(tree[parent_key]):
        tree=remove_admixture(tree, parent_key, 1)
        parent_key=tree[key][branch]
        
    return remove_admix2(tree, key,branch)[0]

def scale_tree(tree, factor):
    for key,node in list(tree.items()):
        node[3]*=factor
        if node_is_admixture(node):
            node[4]*=factor
        tree[key]=node
    return tree

def scale_tree_copy(tree, factor):        
    cop=deepcopy(tree)
    for key,node in list(cop.items()):
        node[3]*=factor
        if node_is_admixture(node):
            node[4]*=factor
        cop[key]=node
    return cop

def parent_is_spouse(tree, key, direction):
    '''
    key is an admixture node. This function checks if the parent in the direction of 'direction' also has a child with the admixture node.
    '''
    parent_key=tree[key][direction]
    offspring_key=get_real_children(tree[key])[0] #there is only one because key is an admixture node
    if node_is_non_admixture(tree[offspring_key]):
        return False
    spouse_key=get_other_parent(tree[offspring_key], key)
    return spouse_key == parent_key

def parent_is_sibling(tree, key, direction):
    '''
    key is an admixture node. This function checks if the parent in the direction of 'direction' is also the child of the parent in the direction of 
    'other_branch(direction)'. 
    '''
    parent_key=tree[key][direction]
    other_parent_key=tree[key][other_branch(direction)] #there is only one because key is an admixture node
    return (parent_key=='r' or other_parent_key in get_real_parents(tree[parent_key]))

    
def get_admixture_branches(tree):
    res=[]
    for key,node in list(tree.items()):
        if node_is_admixture(node):
            res.append((key,0))
            res.append((key,1))

def get_all_admixture_origins(tree):
    res={}
    for key,node in list(tree.items()):
        if node_is_admixture(node):
            res[key]=(node[3], get_parents(node)[0])
    return res
        

def is_root(*keys):
    ad=[key=='r' for key in list(keys)]
    return any(ad)

def to_aarhus_admixture_graph(tree):
    leaves=[]
    inner_nodes=[]
    edges=[]
    admixture_proportions=[]
    for key,node in list(tree.items()):
        if node_is_leaf_node(node):
            leaves.append([key])
        else:
            inner_nodes.append([key])
            if node_is_admixture(node):
                admixture_proportions.append([key, node[0],node[1],node[2]])
        ps=get_real_parents(node)
        for p in ps:
            edges.append([key,p])
    return leaves, inner_nodes, edges, admixture_proportions

def to_networkx_format(tree):
    edges=[]
    admixture_nodes=[]
    leaves=[]
    root=['r']
    coalescence_nodes=[]
    edge_lengths=[]
    for key,node in list(tree.items()):
        if node_is_leaf_node(node):
            leaves.append(key)
        else:
            if node_is_coalescence(node):
                coalescence_nodes.append(key)
            if node_is_admixture(node):
                admixture_nodes.append(key)
        ps=get_real_parents(node)
        for p in ps:
            edges.append((p,key))
            branch=mother_or_father(tree, key, p)
            edge_lengths.append(get_branch_length(tree, key, branch))
    return leaves, admixture_nodes, coalescence_nodes, root, edges, edge_lengths
    

def make_consistency_checks(tree, leaf_nodes=None):
    key_to_parents_by_def={key:[] for key in list(tree.keys())}
    key_to_children_by_def={key:[] for key in list(tree.keys())}
    key_to_children_by_family={key:[] for key in list(tree.keys())}
    key_to_parents_by_family={key:[] for key in list(tree.keys())}
    rooted_nodes=[]
    doppel_bands=[]
    pseudo_nodes=[]
    child_is_parent=[]
    recorded_leaf_nodes=[]
    illegal_admixture_props=[]
    illegal_branch_lengths=[]
    for key,node in list(tree.items()):
        parents=[r for r in get_parents(node) if r is not None]
        children=[r for r in get_children(node) if r is not None]
        key_to_parents_by_def[key]+=[r for r in parents if r!='r']
        key_to_children_by_def[key]+=children
        if len(children)==2 and children[0]==children[1]:
            doppel_bands.append((key, ('children', children[0],children[1])))
        if len(parents)==2 and parents[0]==parents[1]:
            doppel_bands.append((key, ('parents', parents[0],parents[1])))
        if len(children)==1 and len(parents)==1:
            pseudo_nodes.append((key, ('child', children[0]), ('parent', parents[0])))
        if len(children)==0 and len(parents)==1:
            recorded_leaf_nodes.append(key)
        for r in children:
            key_to_parents_by_family[r]+=[key]
        for r in parents:
            if r!='r':
                key_to_children_by_family[r]+=[key]
            else:
                rooted_nodes.append(key)   
        if list(set(parents).intersection(children)):
            child_is_parent.append((key, ('parents', str(parents)), ('children', str(children))))
        if node_is_non_admixture(node):
            if node[3]<0:
                illegal_branch_lengths.append((key, (node[0], node[3])))
        else:
            if node[3]<0:
                illegal_branch_lengths.append((key, (node[0], node[3])))
            if node[4]<0:
                illegal_branch_lengths.append((key, (node[1], node[4])))
            if node[2]<0 or node[2]>1:
                illegal_admixture_props.append((key, node[2]))
                  
    
    def _transform_dic(dic):
        for key in list(dic.keys()):
            dic[key]=set(dic[key])
        return dic
    key_to_children_by_def=_transform_dic(key_to_children_by_def)
    key_to_parents_by_def=_transform_dic(key_to_parents_by_def)
    key_to_children_by_family=_transform_dic(key_to_children_by_family)
    key_to_parents_by_family=_transform_dic(key_to_parents_by_family)
    def _print_inconsistencies(dic_def,dic_fam, pref=''):
        res=''
        for key in list(dic_def.keys()):
            if dic_def[key]!=dic_fam[key]:
                res+=key+' '+pref +': '+str(dic_def[key])+'><'+str(dic_fam[key])+'  '
        return res
    
    family1=_print_inconsistencies(key_to_children_by_def, key_to_children_by_family, 'ch')
    family2=_print_inconsistencies(key_to_parents_by_def, key_to_parents_by_family, 'pa')
    
    bools=[]
    names=[]
    messages=[]
    
    
    consensus_bool=(key_to_children_by_def == key_to_children_by_family and key_to_parents_by_def == key_to_parents_by_family)
    consensus_message=family1+family2
    bools.append(consensus_bool)
    names.append('consensus')
    messages.append(consensus_message)
    
    roots_bool=(len(rooted_nodes)==2)
    roots_message=str(rooted_nodes)
    bools.append(roots_bool)
    names.append('roots')
    messages.append(roots_message)
    
    doppel_bands_bool=(len(doppel_bands)==0)
    doppel_bands_message=str(doppel_bands)
    bools.append(doppel_bands_bool)
    names.append('doppel_bands')
    messages.append(doppel_bands_message)
    
    pseudo_nodes_bool=(len(pseudo_nodes)==0)
    pseudo_nodes_message=str(pseudo_nodes)
    bools.append(pseudo_nodes_bool)
    names.append('pseudo_nodes')
    messages.append(pseudo_nodes_message)
    
    child_is_parent_bool=(len(child_is_parent)==0)
    child_is_parent_message=str(child_is_parent)
    bools.append(child_is_parent_bool)
    names.append('child_is_parent')
    messages.append(child_is_parent_message)
 
    if leaf_nodes is not None:
        leaf_nodes_bool=(set(leaf_nodes)==set(recorded_leaf_nodes))
    else:
        leaf_nodes_bool=True
    if leaf_nodes_bool:
        leaf_nodes_message=str(recorded_leaf_nodes)
    else:
        sl=set(leaf_nodes)
        srl=set(recorded_leaf_nodes)
        leaf_nodes_message=str(set(sl)-set(srl))+'><'+str(set(srl)-set(sl))
    bools.append(leaf_nodes_bool)
    names.append('leaf_nodes')
    messages.append(leaf_nodes_message) 
    
    illegal_branch_lengths_bool=(len(illegal_branch_lengths)==0)
    illegal_branch_lengths_message=str(illegal_branch_lengths)
    bools.append(illegal_branch_lengths_bool)
    names.append('illegal_branch_lengths')
    messages.append(illegal_branch_lengths_message)
    
    illegal_admixture_props_bool=(len(illegal_admixture_props)==0)
    illegal_admixture_props_message=str(illegal_admixture_props)
    bools.append(illegal_admixture_props_bool)
    names.append('illegal_admixture_props')
    messages.append(illegal_admixture_props_message)

    
    res_bool=all(bools)
    res_dic={name:(bool, message) for name,bool,message in zip(names,bools, messages)}
    
    return res_bool, res_dic
