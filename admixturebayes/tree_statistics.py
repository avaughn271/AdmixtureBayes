from Rtree_operations import (get_categories, get_destination_of_lineages, propagate_married,
                              propagate_admixtures, get_branch_length,
                              get_admixture_proportion,
                              get_admixture_keys_and_proportions)
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
        return {a:b for a,b in zip(firsts+seconds, seconds+firsts)}
    else:
        return {}

def _list_identifier_to_string(list_of_gens):
    return '-'.join(['.'.join(map(str,[c for c in l if c!='_'])) for l in list_of_gens])

def _list_double_to_string(list_of_doubles, digits=3):
    format_string='.'+str(digits)+'f'
    return '-'.join(map(str, [format(double_, format_string) for double_ in list_of_doubles]))

class generate_numbered_nodes(object):

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

class uniform_generator(object):

    def __call__(self):
        return random()
 
def unique_identifier_and_branch_lengths(tree, leaf_order=None):
    leaves,admixture_nodes=get_categories(tree)
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

        sames_dic = make_dics_first_and_second(sames)
        awaited_dic= make_dics_first_and_second(awaited_coalescences)
        gen=[]
        for n,element in enumerate(lineages):
            if element in gone:
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
