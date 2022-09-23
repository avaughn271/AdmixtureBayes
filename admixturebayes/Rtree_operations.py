from copy import deepcopy

def get_trivial_nodes(size):
    return ['s'+str(n+1) for n in range(size)]

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

def update_branch_length(tree,key,branch, new_length):
    tree[key][branch+3]=new_length

def extend_branch(node, pkey, grand_parent_key, p_to_gp):
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

def get_all_branches(tree):
    res=[]
    for key, node in list(tree.items()):
        if node_is_admixture(node):
            res.extend([(key, 0),(key,1)])
        else:
            res.append((key,0))
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
    return list(set(all_branches)-set(descendant_branches))

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

def get_children(node):
    return node[5:7]

def _get_index_of_parent(node, parent):
    if node[0]==parent:
        return 0
    if node[1]==parent:
        return 1
    return -1

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

def get_categories(tree):
    leaves=[]
    admixture_nodes=[]
    for key, node in list(tree.items()):
        if node_is_leaf_node(node):
            leaves.append(key)
        if node_is_admixture(node):
            admixture_nodes.append(key)
    return leaves, admixture_nodes

def get_parent_of_branch(tree, key, branch):
    assert key!='r', 'Tried to access the parent of the root branch'
    return tree[key][branch]

def get_branch_length(tree,key,branch):
    assert key!='r', 'Tried to access the length of the root branch'
    return tree[key][branch+3]

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

def get_sister_branch(tree, parent, key, branch):
    if parent[5]==parent[6]:
        return key, other_branch(branch)
    else:
        if parent[5]==key:
            return parent[6], mother_or_father(tree, parent[6], tree[key][branch])
        elif parent[6]==key:
            return parent[5], mother_or_father(tree, parent[5], tree[key][branch])

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

def insert_children_in_tree(tree):
    children={key:[] for key in tree}
    for key in tree:
        parents = get_real_parents(tree[key])
        for parent in parents:
            if parent!='r':
                children[parent].append(key)
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
    del tree[source_key]
    return tree, (t1,t2,t3,t4,t5), alpha

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
