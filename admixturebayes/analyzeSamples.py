from argparse import ArgumentParser, SUPPRESS
from downstream_analysis_tool import (thinning, iterate_over_output_file, make_Rtree, make_full_tree,
                                    get_pops, topology, make_string_tree)
from find_true_trees import tree_unifier
from copy import deepcopy
import sys
import pandas as pd

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