from argparse import ArgumentParser
from downstream_analysis_tool import (thinning, iterate_over_output_file, make_Rtree, make_full_tree, read_true_values,
                                    get_pops, topology, make_string_tree)
from find_true_trees import tree_unifier
from construct_covariance_choices import read_one_line
from copy import deepcopy
import sys

def run_posterior_main(args):

    possible_summaries={'Rtree': make_Rtree,
                        'full_tree':make_full_tree,
                        'string_tree':make_string_tree,
                        'topology':topology,
                        'pops':get_pops
                        }
    parser = ArgumentParser(usage='pipeline for post analysis')

    parser.add_argument('--mcmc_results', required=True, type=str, help='The output file from an AdmixtureBayes run.')
    parser.add_argument('--covariance', required=True, type=str,
                        help='file containing the covariance matrix with a header with all the population names and a line with the multiplier. It has the ending covariance_and_multiplier.txt.')
    parser.add_argument('--subnodes', default=[], type=str, nargs='+',
                        help='The subset of populations to perform the analysis on. If not declared, the analysis will be done on the full dataset.')
    parser.add_argument('--result_file', default='thinned_samples.csv', type=str,
                        help='The resulting file. It will be comma-separated and contain one column per summary plus a header.')
    parser.add_argument('--slower', default=False,    action='store_true', #ANDREWDEBUG
                        help='This will make the program not calculate the string_tree summary which can be very slow when there are many admixture events. '
                             'As a consequence, the option "--plot estimates" can not be used by AdmixtureBayes plot.')

    parser.add_argument('--thinning_rate', default=10, type=int,
                        help='thinning rate')
    parser.add_argument('--burn_in_fraction', default=0.5, type=float,
                        help='the proportion of the rows that are discarded as burn in period')
    parser.add_argument('--save_summaries', default=['no_admixes', 'topology', 'pops','string_tree'], nargs='*', type=str,
                        help='The list of summaries to save')

    options= parser.parse_args(args)

    if options.subnodes:
        subnodes_wo_outgroup=options.subnodes
        subnodes_with_outgroup=options.subnodes
    else:
        subnodes_with_outgroup=[]
        subnodes_wo_outgroup=[]
    
    outp=read_true_values(true_covariance_and_multiplier=options.covariance,
                          subnodes_wo_outgroup=subnodes_wo_outgroup)
    (emp_covariance_scaled,multiplier)=outp

    thinner=thinning(burn_in_fraction=options.burn_in_fraction, total=options.thinning_rate)

    nodes=read_one_line(options.covariance).split() #this will not include any outgroup.
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

    row_sums.append(possible_summaries['Rtree'](deepcopy(nodes_wo_outgroup),False, subnodes=subnodes_wo_outgroup))
    name_to_rowsum_index('Rtree')
    if multiplier is None:
        add_multiplier=1.0
    else:
        add_multiplier=1.0/multiplier

    row_sums.append(possible_summaries['full_tree'](add_multiplier=add_multiplier,
                                                        subnodes=options.subnodes))
    name_to_rowsum_index('full_tree')

    if options.subnodes:
        nodes=options.subnodes
        nodes_with_outgroup=subnodes_with_outgroup
        nodes_wo_outgroup=subnodes_wo_outgroup
    else:
        nodes=nodes_with_outgroup
    if options.slower:
        row_sums.append(possible_summaries['string_tree'](deepcopy(nodes), tree_unifier()))
        name_to_rowsum_index('string_tree')
    if not options.slower:
        options.save_summaries.remove('string_tree')

    row_sums.append(possible_summaries['topology'](nodes=nodes))
    name_to_rowsum_index('topology')
    row_sums.append(possible_summaries['pops'](min_w=0.0, keys_to_include=nodes))
    name_to_rowsum_index('pops')

    def save_thin_columns(d_dic):
        return {summ:d_dic[summ] for summ in list(set(options.save_summaries+[]))}
    all_results=iterate_over_output_file(options.mcmc_results,
                                             cols=['tree', 'add', 'layer', 'no_admixes'],
                                             pre_thin_data_set_function=thinner,
                                             row_summarize_functions=row_sums,
                                             thinned_d_dic=save_thin_columns)

    summaries=list(all_results[0].keys())
    with open(options.result_file, 'w') as f:
        f.write(','.join(summaries)+'\n')
        for row in all_results:
            s_summs=[str(row[summ]) for summ in summaries]
            f.write(','.join(s_summs)+ '\n')
    
if __name__=='__main__':
    import sys
    run_posterior_main(sys.argv[1:])