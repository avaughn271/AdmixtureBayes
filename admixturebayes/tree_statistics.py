from Rtree_operations import (get_categories, get_destination_of_lineages, propagate_married, 
                              propagate_admixtures, get_branch_length,update_parent_and_branch_length, 
                              insert_children_in_tree, rename_root,
                              get_admixture_proportion,
                              get_admixture_keys_and_proportions,
                              direct_all_admixtures)
from copy import deepcopy
from numpy.random import random

def matchmake(single_coalescences, coalescences_on_hold):
    happy_couples=[]
    continuing_singles=[]
    for key,branch in coalescences_on_hold:
        if (key,branch) in single_coalescences:
            happy_couples.append(((key,branch),single_coalescences[(key,branch)]))
            del single_coalescences[(key,branch)]
        else:
            continuing_singles.append((key,branch))
    return single_coalescences, happy_couples, continuing_singles
    
def make_dics_first_and_second(double_list):
    if double_list:
        firsts, seconds=list(map(list,list(zip(*double_list))))
        dic={a:b for a,b in zip(firsts+seconds, seconds+firsts)}
        return dic, firsts, seconds
    else:
        return {},[],[]
    
def unique_identifier(tree, leaf_order=None):
    leaves, coalescences_nodes, admixture_nodes=get_categories(tree)
    if leaf_order is not None:
        assert set(leaves)==set(leaf_order), 'the specified leaf order did not match the leaves of the tree.'
        leaves_ordered=leaf_order
    else:
        leaves_ordered=sorted(leaves)
    ready_lineages=[(key,0) for key in leaves_ordered]
    lineages=deepcopy(ready_lineages)
    gen_to_column=list(range(len(ready_lineages)))
    list_of_gens=[]
    coalescences_on_hold=[]
    gone=[]
    
    res=[]
    
    while True:
        sames, single_coalescences, admixtures=get_destination_of_lineages(tree, ready_lineages)
        waiting_coalescences, awaited_coalescences, still_on_hold = matchmake(single_coalescences, coalescences_on_hold)
        res.append((len(sames), len(waiting_coalescences), len(awaited_coalescences), len(admixtures)))
        
        sames_dic, first_sames, second_sames = make_dics_first_and_second(sames)
        awaited_dic, first_awaited, second_awaited = make_dics_first_and_second(awaited_coalescences)
        gen=[]
        for n,element in enumerate(lineages):
            if element in gone:
                gen.append('_')
            elif element in sames_dic:
                partner_index=lineages.index(sames_dic[element])
                if n<partner_index:
                    gen.append('c')
                else:
                    gen.append(partner_index)
            elif element in awaited_dic:
                partner_index=lineages.index(awaited_dic[element])
                if n<partner_index:
                    gen.append('c')
                else:
                    gen.append(partner_index)
            elif element in admixtures:
                gen.append('a')
            else:
                gen.append('w')
        #removedprin 'gen',gen
        #removedprin 'gone', gone
        list_of_gens,gone, lineages =update_lineages(list_of_gens,gen,gone, lineages, tree)
        for gon in gone:
            lineages.remove(gon)
        gone=[]      
                
    
        #updating lineages
        coalescences_on_hold=still_on_hold+list(waiting_coalescences.values())
        ready_lineages=propagate_married(tree, awaited_coalescences)
        ready_lineages.extend(propagate_married(tree, sames))
        ready_lineages.extend(propagate_admixtures(tree, admixtures))

        #stop criteria
        if len(ready_lineages)==1 and ready_lineages[0][0]=='r':
            break
    return _list_identifier_to_string(list_of_gens)

def admixture_sorted_unique_identifier(tree, leaf_order=None, not_opposite=True):
    '''
    because of a mis-implementation a good amount of calculations could be salvaged by the option not_opposite. It should in a bug-free environment be True.
    '''
    return unique_identifier(direct_all_admixtures(tree, not_opposite), leaf_order)

def _list_identifier_to_string(list_of_gens):
    return '-'.join(['.'.join(map(str,[c for c in l if c!='_'])) for l in list_of_gens])

def _list_double_to_string(list_of_doubles, digits=3):
    format_string='.'+str(digits)+'f'
    return '-'.join(map(str, [format(double_, format_string) for double_ in list_of_doubles]))

class generate_numbered_nodes(object):
    
    def __init__(self):
        self.node_count=0
        self.admixture_count=0
        
    def __call__(self, admixture=False):
        if admixture:
            self.admixture_count+=1
            return 'a' + str(self.admixture_count)
        self.node_count+=1
        return 'n' + str(self.node_count)
    
class generate_predefined_list(object):
    
    def __init__(self, listi):
        self.listi=listi
        
    def __call__(self):
        return float(self.listi.pop(0))
    
class generate_predefined_list_string(object):
    
    def __init__(self, listi):
        self.listi=listi
        
    def __call__(self):
        return self.listi.pop(0)
    
class same_number(object):
    
    def __init__(self, value):
        self.value=value
    
    def __call__(self):
        return self.value

class uniform_generator(object):
    
    def __call__(self):
        return random()
  
def identifier_to_tree(identifier, leaves=None, inner_nodes=None, branch_lengths=None, admixture_proportions=None):
    '''
    Transforms an identifier of the form qwert-uio-asdfg-jk into a dictionary tree using the generators of leaves, inner_nodes, branch_lengths and admixture_proportions.
    '''
    print("dddd")
    levels=identifier.split('-')
    n_leaves=len(levels[0].split('.'))
    leaf_values=[leaves() for _ in range(n_leaves)]
    tree={leaf:[None]*5 for leaf in leaf_values}
    trace_lineages=[(leaf,0) for leaf in leaf_values]
    print(leaf_values)
    if inner_nodes is None:
        inner_nodes=generate_numbered_nodes()
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
        assert len(trace_lineages)==len(identifier_lineages), 'the number of traced lineages did not match the number of lineages in the identifier '+\
                                                               '\n\n'+'trace_lineages:'+'\n'+str(trace_lineages)+\
                                                               '\n\n'+'identifier_lineages:'+'\n'+str(identifier_lineages)
        parent_index={}
        indexes_to_be_removed=[]
        for n,identifier_lineage in enumerate(identifier_lineages):
            if identifier_lineage=='c':
                ##there is a coalecence for the n'th lineage, and it should be replaced by a new lineage
                new_key=inner_nodes()
                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length(tree, old_key, old_branch, new_key, new_branch_length)
                tree[new_key]=[None]*5
                parent_index[n]=new_key
                trace_lineages[n]=(new_key,0)
            elif identifier_lineage=='w':
                pass
            elif identifier_lineage=='a':
                new_key=inner_nodes(admixture=True)
                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length(tree, old_key, old_branch, new_key, new_branch_length)
                new_admixture_proportion=admixture_proportions()
                tree[new_key]=[None,None,new_admixture_proportion,None,None]
                trace_lineages[n]=(new_key,0)
                trace_lineages.append((new_key,1))
            else:
                ##there is a coalescence but this lineage disappears
                try:
                    new_key=parent_index[int(identifier_lineage)]
                except KeyError as e:
                    print(e)
                
                old_key,old_branch=trace_lineages[n]
                new_branch_length=branch_lengths()
                tree=update_parent_and_branch_length(tree, old_key, old_branch, new_key, new_branch_length)
                indexes_to_be_removed.append(n)
        
        ##remove lineages
        trace_lineages=[trace_lineage for n,trace_lineage in enumerate(trace_lineages) if n not in indexes_to_be_removed]
    root_key=new_key
    del tree[root_key]
    tree=rename_root(tree, new_key)
    
    return insert_children_in_tree(tree)
              
def identifier_to_tree_clean(identifier, **kwargs):
    print("a")     
    ad2, branch_lengths, admixture_proportions=divide_triple_string(identifier)
    tree_good2= identifier_to_tree(ad2,
                                   branch_lengths=string_to_generator(branch_lengths), 
                                   admixture_proportions=string_to_generator(admixture_proportions),
                                   **kwargs)
    return tree_good2

def topological_identifier_to_tree_clean(identifier, **kwargs):
    print("b")
    tree_good2= identifier_to_tree(identifier,
                                   branch_lengths=uniform_generator(),
                                   admixture_proportions=uniform_generator(),
                                   **kwargs)
    return tree_good2

def unique_identifier_and_branch_lengths(tree, leaf_order=None):
    leaves, coalescences_nodes, admixture_nodes=get_categories(tree)
    if leaf_order is not None:
        assert set(leaves)==set(leaf_order), 'the specified leaf order did not match the leaves of the tree.'
        leaves_ordered=leaf_order
    else:
        leaves_ordered=sorted(leaves)
    ready_lineages=[(key,0) for key in leaves_ordered]
    lineages=deepcopy(ready_lineages)
    list_of_gens=[]
    coalescences_on_hold=[]
    gone=[]
    
    res=[]
    branch_lengths=[]
    admixture_proportions=[]
    
    while True:
        sames, single_coalescences, admixtures=get_destination_of_lineages(tree, ready_lineages)
        waiting_coalescences, awaited_coalescences, still_on_hold = matchmake(single_coalescences, coalescences_on_hold)
        res.append((len(sames), len(waiting_coalescences), len(awaited_coalescences), len(admixtures)))
        
        sames_dic, first_sames, second_sames = make_dics_first_and_second(sames)
        awaited_dic, first_awaited, second_awaited = make_dics_first_and_second(awaited_coalescences)
        gen=[]
        for n,element in enumerate(lineages):
            if element in gone:
                print('entered gone!')
                gen.append('_')
            elif element in sames_dic:
                partner_index=lineages.index(sames_dic[element])
                if n<partner_index:
                    gen.append('c')
                    branch_lengths.append(get_branch_length(tree, element[0],element[1]))
                else:
                    gen.append(partner_index)
                    branch_lengths.append(get_branch_length(tree, element[0],element[1]))
            elif element in awaited_dic:
                partner_index=lineages.index(awaited_dic[element])
                if n<partner_index:
                    gen.append('c')
                    branch_lengths.append(get_branch_length(tree, element[0],element[1]))
                else:
                    gen.append(partner_index)
                    branch_lengths.append(get_branch_length(tree, element[0],element[1]))
            elif element in admixtures:
                gen.append('a')
                branch_lengths.append(get_branch_length(tree, element[0],element[1]))
                admixture_proportions.append(get_admixture_proportion(tree, child_key=element[0],child_branch=element[1]))
            else:
                gen.append('w')
        #removedprin 'gone', gone
        list_of_gens,gone, lineages =update_lineages(list_of_gens,gen,gone, lineages, tree)
        for gon in gone:
            lineages.remove(gon)
        gone=[]
    
        #updating lineages
        coalescences_on_hold=still_on_hold+list(waiting_coalescences.values())
        ready_lineages=propagate_married(tree, awaited_coalescences)
        ready_lineages.extend(propagate_married(tree, sames))
        ready_lineages.extend(propagate_admixtures(tree, admixtures))

        #stop criteria
        if len(ready_lineages)==1 and ready_lineages[0][0]=='r':
            break
    return ';'.join([_list_identifier_to_string(list_of_gens),
                     _list_double_to_string(branch_lengths, 9),
                     _list_double_to_string(admixture_proportions, 3)])

def list_to_generator(listi):
    return generate_predefined_list(listi)
    
def string_to_generator(stringi, floats=True):
    if floats:
        return list_to_generator(list(map(str,stringi.split('-'))))
    return list_to_generator(stringi.split('-'))
    
def divide_triple_string(stringbig):
    return stringbig.split(';')
    
def update_lineages(lists, new, gone, lineages, tree):
    for n,element in enumerate(new):
        if element=='a':
            key_of_admix_child, branch_of_admix_child=lineages[n]
            admix_key=tree[key_of_admix_child][branch_of_admix_child]
            lineages.append((admix_key,1))
            lineages[n]=(admix_key,0)
        elif element=='c':
            key_of_child, branch_of_child=lineages[n]
            key=tree[key_of_child][branch_of_child]
            lineages[n]=(key,0)
        elif element=='w' or element=='_':
            pass
        else:
            gone.append(lineages[n])
    lists.append(new)
    return lists, gone, lineages

def get_admixture_proportion_string(tree):
    keys, props= get_admixture_keys_and_proportions(tree)
    return '-'.join(keys)+';'+_list_double_to_string(props)
