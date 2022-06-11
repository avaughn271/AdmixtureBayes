from argparse import ArgumentParser, SUPPRESS
from downstream_analysis_tool import (thinning, iterate_over_output_file, always_true, make_Rtree, make_full_tree, read_true_values,
                                      make_Rcovariance, cov_truecov, topology_identity,get_pops,compare_pops,extract_number_of_sadmixes, 
                                      read_one_line,summarize_all_results, topology, float_mean,
                                       subgraph, subsets,
                                      make_string_tree, topology_without_outgroup)
from find_true_trees import tree_unifier
from copy import deepcopy
import sys
from custom_summary import all_custom_summaries

def run_posterior_main(args):

    possible_summaries={'Rtree': make_Rtree,
                        'full_tree':make_full_tree,
                        'string_tree':make_string_tree,
                        'Rcov':make_Rcovariance,
                        'cov_dist':cov_truecov,
                        'topology':topology,
                        'topology_without_outgroup':topology_without_outgroup,
                        'subgraph':subgraph,
                        'subsets':subsets,
                        'top_identity':topology_identity,
                        'pops':get_pops,
                        'set_differences':compare_pops,
                        'no_sadmixes':extract_number_of_sadmixes
                        }
    possible_summaries.update(all_custom_summaries())

    parser = ArgumentParser(usage='pipeline for post analysis')

    parser.add_argument('--mcmc_results', required=True, type=str, help='The output file from an AdmixtureBayes run.')
    parser.add_argument('--covariance', required=True, type=str,
                        help='file containing the covariance matrix with a header with all the population names and a line with the multiplier. It has the ending covariance_and_multiplier.txt.')
    parser.add_argument('--subnodes', default=[], type=str, nargs='+',
                        help='The subset of populations to perform the analysis on. If not declared, the analysis will be done on the full dataset.')
    parser.add_argument('--result_file', default='posterior_distributions.csv', type=str,
                        help='The resulting file. It will be comma-separated and contain one column per summary plus a header.')
    parser.add_argument('--faster', default=False, action='store_true',
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
    parser.add_argument('--custom_summaries', default=[], nargs='*', choices=list(possible_summaries.keys()),
                        help='This will add summaries (to both calculate_summaries and save_summaries). They are defined in the class custom_summary.py.')
    parser.add_argument('--reroot', default='', type=str,
                        help='If a population is given, this will reroot all graphs to the population given. '
                             'Supplying this flag will often result in error because some graphs can not be rerooted'
                             ' without becoming cyclic. See flag --reroot_erro to remedy this.')
    parser.add_argument('--reroot_error', default='stop', choices=['stop','force','ignore'],
                        help='If a rerooting can not be performed (due to admixture events with a sink in a path between '
                             'the reroot population and the graph root), this will either 1) "stop" the program'
                             'or 2) "force" the rerooting through by removing the problematic admixture events or '
                             '3) "ignore" rerooting by simply not doing the rerooting when not possible.')
    parser.add_argument('--summarize_posterior_distributions', default=False,
                        help='If set to true, the posterior distibutions will be summarized even further.')
    parser.add_argument('--min_w', default=0.0, type=float,
                        help='a lower threshold of which descendants matter when the consensus_method is descendant_frequencies.')
    parser.add_argument('--constrain_number_of_admixes', default='', type=str,
                        choices=['', 'true_val'] + list(map(str, list(range(21)))),
                        help='The number of admixture events that there are constrained on in the data set. If negative there are no constraints')
    parser.add_argument('--constrain_number_of_effective_admixes', default='',
                        choices=['', 'true_val'] + list(map(str, list(range(21)))), type=str,
                        help='The number of effective(visible)_admixture events that there are constrained on in the data set. If negative there are no constraints.')
    parser.add_argument('--constrain_sadmix_trees', default=False, action='store_true',
                        help='this will remove the graphs which has invisible admixtures. This will produce worse, but perhaps more easily interpretable results.')
    parser.add_argument('--no_sort', default=False, action='store_true',
                        help='often the tree is sorted according to the leaf names. no_sort willl assumed that they are not sorted according to this but sorted according to ')
    parser.add_argument('--use_cols', default=['tree', 'add', 'layer', 'no_admixes'], type=str, nargs='+',
                        help='The columns to load from the input file')
    parser.add_argument('--outgroup_name', default='', type=str, help='Name of the outgroup. By default this is argument is empty meaning that the outgroup will not be included in any summary.')
    parser.add_argument('--emp_m_scale', type=str, default='', help=SUPPRESS)
    parser.add_argument('--emp_variance_correction', type=str, default='', help=SUPPRESS)
    parser.add_argument('--emp_df', type=str, default='', help=SUPPRESS)
    parser.add_argument('--emp_covariance_and_multiplier', default='', type=str, help=SUPPRESS)
    parser.add_argument('--emp_covariance_reduced', default='', type=str, help=SUPPRESS)

    parser.add_argument('--choice_if_no_thinned_graphs', default='error', choices=['error', 'nearest_admixture_events'],
                        help='If the thinning leaves no graphs left, this is what will be done in stead. error will throw an error and nearest_admixture_events will expand the band of allowed number of admixture events(if the chain has been thinned on number of admixture events).')

    parser.add_argument('--summary_summaries', default=['mean'], nargs='*', type=str,
                        help=SUPPRESS)#'How each list is summarized as a single, numerical value. If it doesnt have the same length as save summaries the arguments will be repeated until it does')
    parser.add_argument('--number_of_top_pops', default=10, type=int,
                        help='if top_pops is added to summary_summaries this is the number of set topologies saved. negative values means all topologies are saved.')
    parser.add_argument('--true_scaled_tree',  type=str, default='',help=SUPPRESS)
    parser.add_argument('--true_tree',  type=str, default='',help=SUPPRESS)
    parser.add_argument('--true_add',  type=str, default='',help=SUPPRESS)
    parser.add_argument('--true_covariance_reduced',  type=str, default='',help=SUPPRESS)
    parser.add_argument('--true_covariance_and_multiplier',  type=str, default='',help=SUPPRESS)
    parser.add_argument('--true_no_admix',  type=str, default='',help=SUPPRESS)
    parser.add_argument('--treemix_post_analysis', action='store_true', default=False,
                        help=SUPPRESS)#'this will convert the treemix input fil ../../../../Dropbox/Bioinformatik/AdmixtureBayes/test_final_grid/ai_2_5true/_true_tree.txtes into a suitable csv file for ',help=SUPPRESS)
    parser.add_argument('--treemix_tree', default='', type=str,help=SUPPRESS)
    parser.add_argument('--treemix_add', default='', type=str, help=SUPPRESS)
    parser.add_argument('--treemix_full_tree', default='',help=SUPPRESS)
    parser.add_argument('--treemix_csv_output', default='treemix.csv', type=str,help=SUPPRESS)
    parser.add_argument('--subgraph_file', default='', type=str,
                        help='file where each line has a space separated list of leaf labels to calculate subtrees from. If a double underscore(__) occurs, it means that the following two arguments are max number of sub topologies and total posterior probability.')

    options= parser.parse_args(args)

    assert not ('string_tree' in options.calculate_summaries and not 'full_tree' in options.calculate_summaries), 'The full tree flag is needed for the string tree'

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

    outp=read_true_values(true_scaled_tree=options.true_scaled_tree,
                          true_tree=options.true_tree,
                          true_add=options.true_add,
                          true_covariance_reduced=options.true_covariance_reduced,
                          true_covariance_and_multiplier=options.true_covariance_and_multiplier,
                          true_no_admix=options.true_no_admix,
                          subnodes_with_outgroup=subnodes_with_outgroup,
                          subnodes_wo_outgroup=subnodes_wo_outgroup)
    true_scaled_tree, true_tree, true_add, true_covariance_reduced, (true_covariance_scaled,true_multiplier), true_no_admix, _, _, _=outp
    outp=read_true_values(true_covariance_reduced=options.emp_covariance_reduced,
                          true_covariance_and_multiplier=options.covariance,
                          true_m_scale=options.emp_m_scale,
                          subnodes_with_outgroup=subnodes_with_outgroup,
                          subnodes_wo_outgroup=subnodes_wo_outgroup)
    _, _, _, emp_covariance_reduced, (emp_covariance_scaled,multiplier), _, emp_m_scale, vc, df=outp


    thinner=thinning(burn_in_fraction=options.burn_in_fraction, total=options.thinning_rate)

    nodes=read_one_line(options.covariance).split() #this will not include any outgroup.
    nodes_wo_outgroup = deepcopy(nodes)
    if options.outgroup_name:
        assert options.outgroup_name not in nodes_wo_outgroup, 'The outgroup_name=' + options.outgroup_name + ' occured in the covariance which an outgroup should not.'
        nodes_with_outgroup = nodes_wo_outgroup + [options.outgroup_name]
    else:
        nodes_with_outgroup = deepcopy(nodes_wo_outgroup)
    if not options.no_sort:
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
    possible_summary_summaries={'mean':float_mean}

    special_summaries=['Rtree','full_tree','string_tree','subgraph','Rcov','cov_dist','topology','topology_without_outgroup','top_identity','pops','subsets', 'set_differences','no_admixes']
    if 'Rtree' in options.calculate_summaries:
        row_sums.append(possible_summaries['Rtree'](deepcopy(nodes_wo_outgroup),options.constrain_sadmix_trees, subnodes=subnodes_wo_outgroup, outgroup_name=options.outgroup_name))
        name_to_rowsum_index('Rtree')
    if 'full_tree' in options.calculate_summaries:
        if multiplier is None:
            add_multiplier=1.0
        else:
            add_multiplier=1.0/multiplier

        row_sums.append(possible_summaries['full_tree'](add_multiplier=add_multiplier,
                                                        outgroup_name=options.outgroup_name,
                                                        remove_sadtrees=options.constrain_sadmix_trees,
                                                        subnodes=options.subnodes,
                                                        reroot_population=options.reroot,
                                                        reroot_method=options.reroot_error))
        name_to_rowsum_index('full_tree')

    if options.subnodes:
        nodes=options.subnodes
        nodes_with_outgroup=subnodes_with_outgroup
        nodes_wo_outgroup=subnodes_wo_outgroup
    else:
        nodes=nodes_with_outgroup
    if 'string_tree' in options.calculate_summaries and not options.faster:
        row_sums.append(possible_summaries['string_tree'](deepcopy(nodes), tree_unifier())) #calling make_string_tree
        name_to_rowsum_index('string_tree')
    if options.faster:
        if 'string_tree' in options.calculate_summaries:
            options.calculate_summaries.remove('string_tree')

        if 'string_tree' in options.save_summaries:
            options.save_summaries.remove('string_tree')

    if 'topology' in options.calculate_summaries:
        row_sums.append(possible_summaries['topology'](nodes=nodes))
        name_to_rowsum_index('topology')
        assert  'topology_without_outgroup' not in options.calculate_summaries, 'not possible to have both topology_without_outgroup and topology'
    if 'pops' in options.calculate_summaries:
        row_sums.append(possible_summaries['pops'](min_w=options.min_w, keys_to_include=nodes))
        name_to_rowsum_index('pops')

    for summary in possible_summaries:
        if summary not in special_summaries:
            if summary in options.calculate_summaries or summary in options.custom_summaries:
                row_sums.append(possible_summaries[summary]())
                name_to_rowsum_index(summary)

    def save_thin_columns(d_dic):
        return {summ:d_dic[summ] for summ in list(set(options.save_summaries+options.custom_summaries))}

    all_results,_=iterate_over_output_file(options.mcmc_results,
                                             cols=options.use_cols,
                                             pre_thin_data_set_function=thinner,
                                             while_thin_data_set_function=always_true,
                                             row_summarize_functions=row_sums,
                                             thinned_d_dic=save_thin_columns,
                                             full_summarize_functions=[])

    if not options.summarize_posterior_distributions:
        summaries=list(all_results[0].keys())
        with open(options.result_file, 'w') as f:
            f.write(','.join(summaries)+'\n')
            for row in all_results:
                s_summs=[str(row[summ]) for summ in summaries]
                f.write(','.join(s_summs)+ '\n')
        sys.exit()

    n=len(options.save_summaries)
    summary_summaries=options.summary_summaries
    while len(summary_summaries)<n:#repeat arguments until the number of arguments is correct
        summary_summaries+=options.summary_summaries

    summary_summaries_functions=[possible_summary_summaries[summ] for summ in summary_summaries]


    summ_results=summarize_all_results(all_results, options.save_summaries, summary_summaries_functions)
    res=[]
    header=[]
    with open(options.result_file, 'w') as f:
        for n,(summ_func_name, summ_name) in enumerate(zip(summary_summaries, options.save_summaries)):
            res.append(summ_results[n])
            header.append(summ_name+'_'+summ_func_name)
        f.write(','.join(['mcmc_results']+header)+'\n')
        f.write(','.join([options.mcmc_results]+list(map(str,res))))
    
if __name__=='__main__':
    import sys
    run_posterior_main(sys.argv[1:])