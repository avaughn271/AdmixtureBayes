from argparse import ArgumentParser, SUPPRESS

import construct_starting_trees_choices
from construct_covariance_choices import get_covariance, estimate_degrees_of_freedom_scaled_fast
from posterior import posterior_class
from MCMCMC import MCMCMC
import os
from numpy import random
import math
import pandas
from meta_proposal import simple_adaptive_proposal

import Rtree_operations
import tree_statistics
import Rtree_to_covariance_matrix
from copy import deepcopy

def removefile(filename):
    if os.path.exists(filename):
        os.remove(filename)

def get_summary_scheme(no_chains=1):
    summaries=[construct_starting_trees_choices.s_posterior(),
               construct_starting_trees_choices.s_likelihood(),
               construct_starting_trees_choices.s_prior(),
               construct_starting_trees_choices.s_no_admixes(),
               construct_starting_trees_choices.s_variable('add', output='double'), 
               construct_starting_trees_choices.s_total_branch_length(),
               construct_starting_trees_choices.s_basic_tree_statistics(Rtree_operations.get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               construct_starting_trees_choices.s_basic_tree_statistics(Rtree_to_covariance_matrix.get_populations_string, 'descendant_sets', output='string'),
               construct_starting_trees_choices.s_basic_tree_statistics(tree_statistics.unique_identifier_and_branch_lengths, 'tree', output='string'),
               construct_starting_trees_choices.s_basic_tree_statistics(tree_statistics.get_admixture_proportion_string, 'admixtures', output='string')]
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    sample_verbose_scheme_first=deepcopy(sample_verbose_scheme)
    return [sample_verbose_scheme_first]+[{}]*(no_chains-1), summaries

def main(args):
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"

    parser = ArgumentParser(usage='pipeline for Admixturebayes')

    #input/output options
    parser.add_argument('--input_file', type=str, required=True, help='the input file of the pipeline. It should be of the same type as the treemix input file with a header of population names and each line representing a snp (unless --covariance_pipeline is altered).')
    parser.add_argument('--result_file', type=str, default='mcmc_samples.csv', help='file in which to save results.')

    parser.add_argument('--outgroup', type=str, default='',
                        help='The name of the population that should be outgroup for the covariance matrix. If the covariance matrix is supplied at stage 8 , this argument is not needed.')
    parser.add_argument('--save_covariance', default=False, action='store_true', help='saving the covariance matrix')
    #Important arguments
    parser.add_argument('--bootstrap_blocksize', type=int, default=1000,
                        help='the size of the blocks to bootstrap in order to estimate the degrees of freedom in the wishart distribution')

    #convenience arguments
    parser.add_argument('--starting_temp', type=str, required=True, help='the inputuu')

    parser.add_argument('--ending_temp', type=str, required=True, help='the inputppp')

    parser.add_argument('--temp_scaling', type=str, required=True, help='the inputnnn')

    parser.add_argument('--iter_per_temp', type=str, required=True, help='the inputxx')

    #convenience arguments
    parser.add_argument('--verbose_level', default='normal', choices=['normal', 'silent'],
                        help='this will set the amount of status out prints throughout running the program.')


    options=parser.parse_args(args)
    assert not (any((i < 8 for i in [6,8,9])) and not options.outgroup), 'In the requested analysis, the outgroup needs to be specified by the --outgroup flag and it should match one of the populations'

    #Here is the only thing we should be changing.
    temp = pandas.read_csv(options.input_file, sep ="\s+")
    colnames = list(temp.columns.values)
    for pop_name in colnames:
        assert pop_name.isalnum(), 'Population names can only contain alphanumeric characters (A-Z, a-z, and 0-9). Special characters such as commas, hyphens, and underscores are not allowed.'
    
    assert options.outgroup in colnames, 'The outgroup name is not in the given list of populations. Population names are case-sensitive.'
    
    temp = temp[sorted(colnames)]
    temp.to_csv(os.getcwd() + "/temp_input.txt", sep =" ", index = False)

    mp= [simple_adaptive_proposal(['deladmix', 'addadmix', 'rescale', 'rescale_add', 'rescale_admixtures', 'rescale_constrained', 'sliding_regraft'],
     [1, 1, 1, 1, 1, 1, 1]) for _ in range(1)]

    with open(os.getcwd() + "/temp_input.txt", 'r') as f:
        full_nodes = f.readline().rstrip().split()
    reduced_nodes=deepcopy(full_nodes)
    reduced_nodes.remove(options.outgroup)

    estimator_arguments=dict(reducer=options.outgroup, nodes=full_nodes, add_variance_correction_to_graph=True, save_variance_correction=True)
                             
    covariance=get_covariance(os.getcwd() + "/temp_input.txt", 
    full_nodes=full_nodes,
     reduce_covariance_node=options.outgroup,
     estimator_arguments=estimator_arguments)

    estimator_arguments['save_variance_correction']=False
    df=estimate_degrees_of_freedom_scaled_fast(os.getcwd() + "/temp_input.txt",
                                            bootstrap_blocksize=options.bootstrap_blocksize,
                                            cores=1,
                                            est=estimator_arguments, 
                                            verbose_level=options.verbose_level)

    multiplier=covariance[1]
    #NEWEDIT, change the thing below to False, then []
    starting_trees=construct_starting_trees_choices.get_starting_trees([], 1, adds=[], nodes=reduced_nodes)

    summary_verbose_scheme, summaries=get_summary_scheme(no_chains=1)

    posterior = posterior_class(emp_cov=covariance[0], M=df, multiplier=covariance[1], nodes=reduced_nodes)

    removefile("covariance_without_reduce_name.txt")
    removefile("variance_correction.txt")
    removefile("temp_starttree.txt")
    removefile("temp_start_tree.txt")
    removefile("temp_add.txt")

    if options.save_covariance:
        removefile(os.getcwd() + "/covariance_matrix.txt")
        Liness = open("covariance_and_multiplier.txt", 'r').readlines()
        covarfile = open(os.getcwd() + "/covariance_matrix.txt", "a")
        covarfile.writelines(Liness)
        covarfile.close()

    removefile("covariance_and_multiplier.txt")
    removefile(os.getcwd() + "/temp_input.txt")

    if os.path.exists(os.getcwd() + "/temp_adbayes"):
        os.rmdir(os.getcwd() + "/temp_adbayes")


    
    StartingTemp =  float(options.starting_temp) #    100  
    EndingTemp = float(options.ending_temp) #  0.0001
    TempDecrease = float(options.temp_scaling) # 0.9
    NumberAtEach = int(options.iter_per_temp) # 5000

    MCMCMC(TempDecrease, starting_trees=starting_trees,
            posterior_function= posterior,
            summaries=summaries,
            temperature_scheme=[StartingTemp],
            printing_schemes=summary_verbose_scheme,
            iteration_scheme=[NumberAtEach]*int(math.log(EndingTemp / StartingTemp) / math.log(TempDecrease)),
            proposal_scheme= mp,
            multiplier=multiplier,
            result_file=options.result_file,
            n_arg=int(math.log(EndingTemp / StartingTemp) / math.log(TempDecrease)), verboseee=options.verbose_level)

if __name__=='__main__':
    import sys
    main(sys.argv[1:])
