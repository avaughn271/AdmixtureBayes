from argparse import ArgumentParser
from copy import deepcopy
import sys
import pandas as pd

from tree_statistics import (identifier_to_tree_clean, generate_predefined_list_string, admixture_sorted_unique_identifier, unique_identifier_and_branch_lengths)
from Rtree_to_covariance_matrix import get_populations

from Rtree_operations import (change_admixture, get_categories)

def identity(x):
    return x

def iterate_over_output_file(outfile, 
                             cols=[], 
                             row_summarize_functions=[],
                             thinned_d_dic=identity,
                             **constant_kwargs):
    
    df= pd.read_csv(outfile, usecols=cols, dtype={'no_admixes':object})
    df = df[cols]
    df= df.tail(1)
    all_results=[]
    
    for n,(i,r) in enumerate(df.iterrows()):
        d_dic={colname:r[k] for k, colname in enumerate(cols)}
        d_dic.update(constant_kwargs)

        for row_summarize_function in row_summarize_functions:
            add_dic, skip=row_summarize_function(**d_dic)
            d_dic.update(add_dic)
        all_results.append(thinned_d_dic(d_dic))
    return all_results

class make_Rtree(object):
    
    def __init__(self, nodes_to_be_sorted, subnodes=[], outgroup_name=''):
        self.nodes=sorted(nodes_to_be_sorted)
        self.outgroup_name=outgroup_name
        
    def __call__(self, tree, **not_needed):
        first_level=tree.split('-')[0]
        no_pops=len(first_level.split('.'))
        if no_pops==len(self.nodes): #checking if
            Rtree=identifier_to_tree_clean(tree, leaves=generate_predefined_list_string(deepcopy(self.nodes)))
        elif no_pops==len(self.nodes)+1 and self.outgroup_name:
            self.nodes=sorted(self.nodes+[self.outgroup_name])
            Rtree = identifier_to_tree_clean(tree, leaves=generate_predefined_list_string(deepcopy(self.nodes)))
        else:
            assert False, 'Either the outgroup name was not specified or something is seriously wrong because' \
                          ' the number of nodes did not match the size of the trees'

        return {'Rtree':Rtree}, False
    
class make_full_tree(object):
    
    def __init__(self, outgroup_name='out', subnodes=[]):
        self.outgroup_name=outgroup_name
        
    def __call__(self, Rtree=None, add=None, **kwargs):
        full_tree=deepcopy(Rtree)
        return {'full_tree':full_tree}, False

class make_string_tree(object):

    def __init__(self, nodes, tree_unifier=None):
        self.nodes=nodes
        self.tree_unifier=tree_unifier
        self.node_string='='.join(self.nodes)+'='

    def __call__(self, full_tree, **kwargs):
        stree=unique_identifier_and_branch_lengths(full_tree, leaf_order=self.nodes)
        if self.tree_unifier is not None:
            stree=self.tree_unifier(stree)
            
        string_tree=self.node_string+stree
        return {'string_tree':string_tree},  False

class topology(object):
    
    def __init__(self, nodes):
        self.nodes=nodes
        
    def __call__(self, Rtree=None, **kwargs):
        if 'string_tree' in kwargs:
            topology=kwargs['string_tree'].split('=')[-1].split(';')[0]
            return {'topology':topology},False
        if Rtree is None:
            Rtree=kwargs['full_tree']
        top=admixture_sorted_unique_identifier(Rtree, leaf_order=self.nodes, not_opposite=True)
        return {'topology':top}, False
    
class get_pops(object):
    
    def __init__(self, keys_to_include=None):
        self.keys_to_include=keys_to_include
            
    def __call__(self, full_tree=None, **kwargs):
        if full_tree is None:
            tree=kwargs['Rtree']
        else:
            tree=full_tree
        pops=get_populations(tree, keys_to_include=self.keys_to_include)
        return {'pops':'-'.join(pops)}, False

class tree_unifier(object):
    
    def __init__(self):
        '''
        key:(val1,val2, val3) where 
            key=lookup topology string, 
            val1 is unique topology string, 
            val2 is the permutation of branches and 
            val3 is the (signed) permutation of the admixture proportions.
        '''
        self.seen_trees={}
        
    def __call__(self, stree):
        topology,branches,admixtures=stree.split(';')
        update_dic= analyze_tree(topology, branches, admixtures)
        self.seen_trees.update(update_dic)
        target_topology, branch_permutation, admixture_permutation=self.seen_trees[topology]
        new_branch_string=make_branch_string(branches, branch_permutation)
        new_admixtures_string=make_admixture_string(admixtures, admixture_permutation)
        return ';'.join([target_topology, new_branch_string, new_admixtures_string])
    
def make_branch_string(branches, branch_permutations):
    branch_pieces=branches.split('-')
    return '-'.join([branch_pieces[branch_permutations[i]] for i in range(len(branch_pieces))])

def make_admixture_string(admixes, admixture_permutations):
    if not admixes:
        return ''
    res=[]
    admixture_pieces=admixes.split('-')
    for i in range(len(admixture_pieces)):
        target_admixture=admixture_permutations[i]
        if target_admixture<-0.5:
            res.append(1.0-float(admixture_pieces[abs(target_admixture)-1]))
        else:
            res.append(float(admixture_pieces[abs(target_admixture)-1]))
    return '-'.join(map(str,res))

def analyze_tree(topology, branches, admixtures):
    
    id_branches='-'.join(map(str, list(range(len(branches.split('-'))))))
    id_admixtures='-'.join(map(str, list(range(1,len(admixtures.split('-'))+1))))
    id_stree=';'.join([topology,id_branches, id_admixtures])
    no_leaves=len((id_stree.split('-')[0]).split('.'))
    id_tree=identifier_to_tree_clean(id_stree.strip())
    
    strees= sorted(get_possible_permutation_strees(id_tree))
    top_topology=strees[0].split(';')[0]
    res={}
    for stree in strees:
        lookup_topology, branches_sperm, admixtures_sperm= stree.split(';')
        rf=list(map(round, list(map(float, branches_sperm.split('-')))))
        branches_permutation= list(map(int, rf))
        admixtures_permutation=get_admixtures_permutation(admixtures_sperm)
        res[lookup_topology]=(top_topology, branches_permutation, admixtures_permutation)
    return res
    
def get_admixtures_permutation(admixtures):
    res=[]
    for a in admixtures.split('-'):
        if a:
            number=int(round(float(a)))
            if number>100000:
                res.append(-(number%100000))
            else:
                res.append(number)
    return res

def get_possible_permutation_strees(tree):
    leaves,admixture_keys=get_categories(tree)
    k=len(admixture_keys)
    format_code='{0:0'+str(k)+'b}'
    
    n_trees=[]
    for i in range(2**k):
        pruned_tree = deepcopy(tree)
        bina= format_code.format(i)
        for adm_key,str_bin in zip(admixture_keys, list(bina)):
            int_bin=int(str_bin)
            if int_bin==1:
                pruned_tree[adm_key]=change_admixture(pruned_tree[adm_key])
        n_tree= unique_identifier_and_branch_lengths(pruned_tree)
        n_trees.append(n_tree)
    return n_trees
    
def run_posterior_main(args):

    parser = ArgumentParser(usage='pipeline for post analysis')

    parser.add_argument('--mcmc_results', required=True, type=str, help='The output file from an AdmixtureBayes run.')
    parser.add_argument('--result_file', default='thinned_samples.csv', type=str,
                        help='The resulting file. It will be comma-separated and contain one column per summary plus a header.')

    options= parser.parse_args(args)


    subnodes_with_outgroup=[]
    subnodes_wo_outgroup=[]

    totallist = []
    a = pd.read_csv(options.mcmc_results, nrows=3)
    stringg = (a.loc[0,["descendant_sets"]])[0]
    stringg = (stringg.split('-'))
    for i in stringg:
        totallist.extend(i.split('.'))
    nodes = []
    for i in totallist:
        if i not in nodes:
            nodes.append(i)
    nodes.sort()

    nodes_wo_outgroup = deepcopy(nodes)
    nodes_with_outgroup = deepcopy(nodes_wo_outgroup)
    nodes_with_outgroup.sort()
    nodes_wo_outgroup.sort()
    subnodes_with_outgroup.sort()
    subnodes_wo_outgroup.sort()

    row_sums=[]

    class pointers(object):

        def __init__(self):
            self.count=0
            self.dic={}

        def __call__(self, name):
            self.dic[name]=self.count
            self.count+=1

        def __getitem__(self, key):
            return self.dic[key]

    name_to_rowsum_index=pointers()

    row_sums.append(make_Rtree(deepcopy(nodes_wo_outgroup), subnodes=subnodes_wo_outgroup, outgroup_name=''))
    name_to_rowsum_index('Rtree')
    row_sums.append(make_full_tree(outgroup_name='', subnodes=[]))
    name_to_rowsum_index('full_tree')

    nodes=nodes_with_outgroup
    row_sums.append(make_string_tree(deepcopy(nodes), tree_unifier())) #calling make_string_tree
    name_to_rowsum_index('string_tree')
    row_sums.append(topology(nodes=nodes))
    name_to_rowsum_index('topology')
    row_sums.append(get_pops(keys_to_include=nodes))
    name_to_rowsum_index('pops')

    def save_thin_columns(d_dic):
        return {summ:d_dic[summ] for summ in list(set(['no_admixes', 'topology', 'pops','string_tree']+[]))}
    all_results=iterate_over_output_file(options.mcmc_results,
                                             cols=['tree', 'add', 'layer', 'no_admixes'],
                                             row_summarize_functions=row_sums,
                                             thinned_d_dic=save_thin_columns)

    if True:
        summaries=list(all_results[0].keys())
        with open(options.result_file, 'w') as f:
            f.write(','.join(summaries)+'\n')
            for row in all_results:
                s_summs=[str(row[summ]) for summ in summaries]
                f.write(','.join(s_summs)+ '\n')
    
if __name__=='__main__':
    import sys
    run_posterior_main(sys.argv[1:])