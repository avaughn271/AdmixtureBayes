from argparse import ArgumentParser, SUPPRESS
from downstream_analysis_tool import (thinning, iterate_over_output_file, make_Rtree, make_full_tree,
                                    get_pops, topology, make_string_tree)
from copy import deepcopy
import sys
import pandas as pd
from tree_statistics import identifier_to_tree_clean, unique_identifier_and_branch_lengths
from Rtree_operations import change_admixture, get_categories

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
        if topology in self.seen_trees:
            target_topology, branch_permutation, admixture_permutation= self.seen_trees[topology]
        else:
            update_dic= analyze_tree(topology, branches, admixtures)
            if len(update_dic)>100: #computational reasons
                self.seen_trees[topology]=update_dic[topology]
            else:
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
    leaves,_,admixture_keys=get_categories(tree)
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

    possible_summaries={'Rtree': make_Rtree,
                        'full_tree':make_full_tree,
                        'string_tree':make_string_tree,
                        'topology':topology,
                        'pops':get_pops,
                        }
    parser = ArgumentParser(usage='pipeline for post analysis')

    parser.add_argument('--mcmc_results', required=True, type=str, help='The output file from an AdmixtureBayes run.')
    parser.add_argument('--subnodes', default=[], type=str, nargs='+',
                        help='The subset of populations to perform the analysis on. If not declared, the analysis will be done on the full dataset.')
    parser.add_argument('--result_file', default='thinned_samples.csv', type=str,
                        help='The resulting file. It will be comma-separated and contain one column per summary plus a header.')
    parser.add_argument('--slower', default=False,    action='store_true', #ANDREWDEBUG
                        help='This will make the program not calculate the string_tree summary which can be very slow when there are many admixture events. '
                             'As a consequence, the option "--plot estimates" can not be used by AdmixtureBayes plot.')

    parser.add_argument('--thinning_rate', default=10, type=int,
                        #help='an upper limit on the number of rows to reduce computational pressure')
                        help='thinning rate')
    parser.add_argument('--burn_in_fraction', default=0.5, type=float,
                        help='the proportion of the rows that are discarded as burn in period')
    parser.add_argument('--calculate_summaries', default=['Rtree', 'pops','full_tree','string_tree','topology'], choices=list(possible_summaries.keys()),
                        nargs='*', type=str, help='The summaries to calculate')
    parser.add_argument('--save_summaries', default=['no_admixes', 'topology', 'pops','string_tree'], nargs='*', type=str,
                        help='The list of summaries to save')
    parser.add_argument('--min_w', default=0.0, type=float,
                        help='a lower threshold of which descendants matter when the consensus_method is descendant_frequencies.')
    parser.add_argument('--use_cols', default=['tree', 'add', 'layer', 'no_admixes'], type=str, nargs='+',
                        help='The columns to load from the input file')
    parser.add_argument('--outgroup_name', default='', type=str, help='Name of the outgroup. By default this is argument is empty meaning that the outgroup will not be included in any summary.')

    options= parser.parse_args(args)

    if options.subnodes:
        if not options.outgroup_name:
            subnodes_wo_outgroup=options.subnodes
            subnodes_with_outgroup=options.subnodes
        elif options.outgroup_name in options.subnodes:
            subnodes_with_outgroup=options.subnodes
            subnodes_wo_outgroup=deepcopy(options.subnodes)
            subnodes_wo_outgroup.remove(options.outgroup_name)
        else:
            subnodes_with_outgroup=deepcopy(options.subnodes)+[options.outgroup_name]
            subnodes_wo_outgroup=options.subnodes
    else:
        subnodes_with_outgroup=[]
        subnodes_wo_outgroup=[]

    thinner=thinning(burn_in_fraction=options.burn_in_fraction, total=options.thinning_rate)

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
    if options.outgroup_name:
        assert options.outgroup_name not in nodes_wo_outgroup, 'The outgroup_name=' + options.outgroup_name + ' occured in the covariance which an outgroup should not.'
        nodes_with_outgroup = nodes_wo_outgroup + [options.outgroup_name]
    else:
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

    if 'Rtree' in options.calculate_summaries:
        row_sums.append(possible_summaries['Rtree'](deepcopy(nodes_wo_outgroup),False, subnodes=subnodes_wo_outgroup, outgroup_name=options.outgroup_name))
        name_to_rowsum_index('Rtree')
    if 'full_tree' in options.calculate_summaries:

        row_sums.append(possible_summaries['full_tree'](outgroup_name=options.outgroup_name,
                                                        remove_sadtrees=False,
                                                        subnodes=options.subnodes,
                                                        reroot_population='',
                                                        reroot_method='stop'))
        name_to_rowsum_index('full_tree')

    if options.subnodes:
        nodes=options.subnodes
        nodes_with_outgroup=subnodes_with_outgroup
        nodes_wo_outgroup=subnodes_wo_outgroup
    else:
        nodes=nodes_with_outgroup
    if 'string_tree' in options.calculate_summaries and options.slower:
        row_sums.append(possible_summaries['string_tree'](deepcopy(nodes), tree_unifier())) #calling make_string_tree
        name_to_rowsum_index('string_tree')
    if not options.slower:
        if 'string_tree' in options.calculate_summaries:
            options.calculate_summaries.remove('string_tree')

        if 'string_tree' in options.save_summaries:
            options.save_summaries.remove('string_tree')

    if 'topology' in options.calculate_summaries:
        row_sums.append(possible_summaries['topology'](nodes=nodes))
        name_to_rowsum_index('topology')
    if 'pops' in options.calculate_summaries:
        row_sums.append(possible_summaries['pops'](min_w=options.min_w, keys_to_include=nodes))
        name_to_rowsum_index('pops')

    def save_thin_columns(d_dic):
        return {summ:d_dic[summ] for summ in list(set(options.save_summaries+[]))}
    all_results=iterate_over_output_file(options.mcmc_results,
                                             cols=options.use_cols,
                                             pre_thin_data_set_function=thinner,
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