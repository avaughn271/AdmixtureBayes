import summary
import Rtree_operations
import tree_statistics
import Rtree_to_covariance_matrix
from copy import deepcopy

def get_summary_scheme(light_newick_tree_summaries=False,
                       full_tree=True, 
                       proposals=None, 
                       acceptance_rate_information=False,
                       admixture_proportion_string=False,
                       priors=False,
                       no_chains=1,
                       nodes=None, 
                       verbose_level='normal', 
                       only_coldest_chain=True):
    
    if proposals is not None:
        props=proposals.props
        prop_names=[prop.proposal_name for prop in props]
        adaption=[prop.adaption for prop in props]
    
    summaries=[summary.s_posterior(),
               summary.s_likelihood(),
               summary.s_prior(),
               summary.s_no_admixes(),
               summary.s_variable('add', output='double'), 
               summary.s_average_branch_length(),
               summary.s_total_branch_length(),
               summary.s_basic_tree_statistics(Rtree_operations.get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               summary.s_basic_tree_statistics(Rtree_operations.get_max_distance_to_root, 'max_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_min_distance_to_root, 'min_root'),
               summary.s_basic_tree_statistics(Rtree_operations.get_average_distance_to_root, 'average_root'),
               summary.s_basic_tree_statistics(Rtree_to_covariance_matrix.get_populations_string, 'descendant_sets', output='string')]
    if full_tree:
        summaries.append(summary.s_basic_tree_statistics(tree_statistics.unique_identifier_and_branch_lengths, 'tree', output='string'))
    if admixture_proportion_string:
        summaries.append(summary.s_basic_tree_statistics(tree_statistics.get_admixture_proportion_string, 'admixtures', output='string'))
    if light_newick_tree_summaries:
        summaries.append(summary.s_basic_tree_statistics(tree_statistics.tree_to_0ntree, 'Zero_Ntree',output='string'))
        summaries.append(summary.s_basic_tree_statistics(tree_statistics.tree_to_random_ntree, 'random_Ntree',output='string'))
        summaries.append(summary.s_basic_tree_statistics(tree_statistics.tree_to_mode_ntree, 'mode_Ntree',output='string'))
    if acceptance_rate_information:
        summaries.append(summary.s_variable('mhr', output='double_missing'))
        summaries.append(summary.s_variable('proposal_type', output='string'))
        if proposals is not None:
            for prop_name,adapt in zip(prop_names, adaption):
                if adapt:
                    summaries.append(summary.s_variable(prop_name+"_adap_param", output= 'double_missing'))
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    sample_verbose_scheme_first=deepcopy(sample_verbose_scheme)
    if no_chains==1:
        return [sample_verbose_scheme_first], summaries
    elif only_coldest_chain:
        return [sample_verbose_scheme_first]+[{}]*(no_chains-1), summaries
    else:
        return [sample_verbose_scheme_first]+[sample_verbose_scheme]*(no_chains-1), summaries
