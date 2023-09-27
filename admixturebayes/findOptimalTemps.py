from argparse import ArgumentParser, SUPPRESS

import construct_starting_trees_choices
from construct_covariance_choices import get_covariance, estimate_degrees_of_freedom_scaled_fast
from posterior import posterior_class
import os
import pandas
from meta_proposal import simple_adaptive_proposal
from numpy import linspace
import Rtree_operations
import tree_statistics
import Rtree_to_covariance_matrix
from copy import deepcopy
from pathos.multiprocessing import freeze_support
from MCMC import basic_chain_pool
from MCMCMC import flipping, _pack_everything, _unpack_everything
from math import log

def MCMCMC_findoptimal(starting_trees,    posterior_function, summaries, temperature_scheme, priortempscheme,  printing_schemes, 
           iteration_scheme,  proposal_scheme, n_arg, verboseee,
           no_chains=None, numpy_seeds=None, multiplier= None):
    '''
    this function runs a MC3 using the basic_chain_unpacker. Let no_chains=number of chains. The inputs are
        starting_trees: a list of one or more trees that the chains should be started with
        proposal_scheme: a list of instances of classes that handles proposals and updates of said proposals. It has the basic functions
                                    - prop(x,pks): proposes the next tree(and returns proposal densities and statistics in pks)
                                    - adapt(mhr): updates the proposals based on the mhr ratio, mhr                                    
    '''
    total_permutation=list(range(no_chains))
    xs = starting_trees

    freeze_support() # maybe delete this line.

    pool = basic_chain_pool(summaries, posterior_function, proposal_scheme, numpy_seeds)
    posteriors = [posterior_function(x) for x in xs]

    proposal_updates=[proposal.get_exportable_state() for proposal in proposal_scheme]
    flipacceptances = [0] * (no_chains - 1)
    flipprops= [0] * (no_chains - 1)
    
    cum_iterations=0
    for no_iterations in iteration_scheme:
        #letting each chain run for no_iterations:
        iteration_object=_pack_everything(xs, posteriors, temperature_scheme, priortempscheme, printing_schemes, no_iterations, cum_iterations, proposal_updates, multiplier)
        new_state = pool.order_calculation(iteration_object)
        xs, posteriors, df_add,proposal_updates = _unpack_everything(new_state, summaries, total_permutation)
        #making the mc3 flips and updating:
        xs, posteriors, permut, proposal_updates, flipprops, flipacceptances = flipping(xs, posteriors, temperature_scheme, priortempscheme ,  proposal_updates, flipprops, flipacceptances)
        total_permutation=[total_permutation[n] for n in permut]
        cum_iterations+=no_iterations
    for chain in pool.group:
            chain.process.terminate()
    flippropsstring = []
    for i in range(len(flipprops)):
        fracc = float(flipacceptances[i])  / flipprops[i] 
        flippropsstring.append(fracc)
    return(flippropsstring)

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
    parser.add_argument('--temperature_file', type=str, default='optimal_temperatures.txt', help='the input file of the pipeline. It should be of the same type as the treemix input file with a header of population names and each line representing a snp (unless --covariance_pipeline is altered).')

    parser.add_argument('--outgroup', type=str, default='',
                        help='The name of the population that should be outgroup for the covariance matrix. If the covariance matrix is supplied at stage 8 , this argument is not needed.')
    #Important arguments
    parser.add_argument('--bootstrap_blocksize', type=int, default=1000,
                        help='the size of the blocks to bootstrap in order to estimate the degrees of freedom in the wishart distribution')

    #convenience arguments
    parser.add_argument('--verbose_level', default='silent', choices=['normal', 'silent'],
                        help='this will set the amount of status out prints throughout running the program.')
    
    CURRENTTEMPERATURESCHEME = [1.0, 47.28708045015879, 100.0]
    numberofmcmcflipstodo = 4000
    PRIORTEMPS = linspace(1.0, 1.3, num=len(CURRENTTEMPERATURESCHEME)).tolist()

    upperbounds = [100.0]
    lowerbounds = [1.0]
    ComputedDF = False
    df = -1.0
    while True:
        flippropsstring = []
        for i in range(len(CURRENTTEMPERATURESCHEME)):
            flippropsstring.append(round(CURRENTTEMPERATURESCHEME[i] , 2) )

        print("Assessing mixing with temperature scheme of: ", flippropsstring)
        NUMBEROFCHAINS = len(CURRENTTEMPERATURESCHEME)

        options=parser.parse_args(args)

        temporaryfoldername = (options.temperature_file).replace('.', '') + "_tempfilefolder"
        if ComputedDF == False:
            os.mkdir(os.getcwd() + os.sep + temporaryfoldername)

        assert not (any((i < 8 for i in [6,8,9])) and not options.outgroup), 'In the requested analysis, the outgroup needs to be specified by the --outgroup flag and it should match one of the populations'

        #Here is the only thing we should be changing.
        temp = pandas.read_csv(options.input_file, sep ="\s+")
        colnames = list(temp.columns.values)
        for pop_name in colnames:
            assert pop_name.isalnum(), 'Population names can only contain alphanumeric characters (A-Z, a-z, and 0-9). Special characters such as commas, hyphens, and underscores are not allowed.'

        assert options.outgroup in colnames, 'The outgroup name is not in the given list of populations. Population names are case-sensitive.'

        temp = temp[sorted(colnames)]
        temp.to_csv(os.getcwd() + os.sep + temporaryfoldername + os.sep + "temp_input.txt", sep =" ", index = False)

        mp= [simple_adaptive_proposal(['deladmix', 'addadmix', 'rescale', 'rescale_add', 'rescale_admixtures', 'rescale_constrained', 'sliding_regraft'],
            [1, 1, 1, 1, 1, 1, 1]) for _ in range(NUMBEROFCHAINS)]

        with open(os.getcwd() +os.sep+ temporaryfoldername + os.sep + "temp_input.txt", 'r') as f:
            full_nodes = f.readline().rstrip().split()
        reduced_nodes=deepcopy(full_nodes)
        reduced_nodes.remove(options.outgroup)

        estimator_arguments=dict(reducer=options.outgroup, nodes=full_nodes, add_variance_correction_to_graph=True, save_variance_correction=True)
        if ComputedDF == False:
            covariance=get_covariance(os.getcwd() +os.sep+ temporaryfoldername + os.sep + "temp_input.txt", 
        varcovfilename = os.getcwd() +os.sep+ temporaryfoldername + os.sep + "variance_correction.txt",
        full_nodes=full_nodes,
            reduce_covariance_node=options.outgroup,
            estimator_arguments=estimator_arguments,filename =  os.getcwd() +os.sep + temporaryfoldername + os.sep  +"covariance_and_multiplier.txt")

        estimator_arguments['save_variance_correction']=False
        if ComputedDF == False:
            df=estimate_degrees_of_freedom_scaled_fast(os.getcwd() + os.sep + temporaryfoldername + os.sep  + "temp_input.txt",
                                                        varcovfilename = os.getcwd() +os.sep + temporaryfoldername + os.sep  +"variance_correction.txt",
                                                    bootstrap_blocksize=options.bootstrap_blocksize,
                                                    cores=NUMBEROFCHAINS,
                                                    est=estimator_arguments, 
                                                    verbose_level=options.verbose_level)
            ComputedDF = True

        multiplier=covariance[1]

        starting_trees=construct_starting_trees_choices.get_starting_trees([], NUMBEROFCHAINS, adds=[], nodes=reduced_nodes)

        summary_verbose_scheme, summaries=get_summary_scheme(no_chains=NUMBEROFCHAINS)

        posterior = posterior_class(emp_cov=covariance[0], M=df, multiplier=covariance[1], nodes=reduced_nodes, 
                                        varcovname=os.getcwd() +os.sep + temporaryfoldername + os.sep  + "variance_correction.txt")

        removefile(os.getcwd() + os.sep + temporaryfoldername + os.sep  +  "temp_starttree.txt")
        removefile(os.getcwd() + os.sep + temporaryfoldername + os.sep  +  "temp_start_tree.txt")
        removefile(os.getcwd() + os.sep + temporaryfoldername + os.sep  + "temp_add.txt")

        PRIORTEMPS = linspace(1.0, 1.3, num=len(CURRENTTEMPERATURESCHEME)).tolist()

        MixingRates = MCMCMC_findoptimal(starting_trees=starting_trees,
                posterior_function= posterior,
                summaries=summaries,
                temperature_scheme=CURRENTTEMPERATURESCHEME, priortempscheme = PRIORTEMPS, 
                printing_schemes=summary_verbose_scheme,
                iteration_scheme=[50]* numberofmcmcflipstodo ,
                proposal_scheme= mp,
                no_chains=NUMBEROFCHAINS,
                multiplier=multiplier,  #numpy_seeds = random_seeds,
                n_arg= numberofmcmcflipstodo, verboseee=options.verbose_level)
        currenttempindex = 1



        flippropsstring = []
        for i in range(len(MixingRates)):
            flippropsstring.append( str(round(MixingRates[i] * 100 , 2)  ) + "%" )
        print("Between chain mixing rates: ", flippropsstring)
        if MixingRates[0] > 0.15 and MixingRates[1] > 0.15:
            with open(options.temperature_file, 'w') as fp:
                for item in CURRENTTEMPERATURESCHEME:
                    # write each item on a new line
                    fp.write("%s\n" % item)
            print("Optimal temperature scheme found! Saving results to file: ", options.temperature_file)
            removefile(os.getcwd() +os.sep + temporaryfoldername +os.sep + "covariance_and_multiplier.txt")
            removefile(os.getcwd() +os.sep + temporaryfoldername + os.sep  + "temp_input.txt")
            removefile(os.getcwd() +os.sep + temporaryfoldername + os.sep + "variance_correction.txt")
            if os.path.exists(os.getcwd() +os.sep + temporaryfoldername + os.sep  +"temp_adbayes" ):
                os.rmdir(os.getcwd() +os.sep + temporaryfoldername + os.sep  + "temp_adbayes" )
            if os.path.exists(os.getcwd() +os.sep + temporaryfoldername):
                os.rmdir(os.getcwd() + os.sep + temporaryfoldername)
            break
        elif MixingRates[currenttempindex] > 0.4:
            upperbounds.append(CURRENTTEMPERATURESCHEME[currenttempindex])
            CURRENTTEMPERATURESCHEME[currenttempindex] = 10**((log(CURRENTTEMPERATURESCHEME[currenttempindex], 10) +  log(max(lowerbounds), 10) )/2 )
        elif MixingRates[currenttempindex] < 0.15:
            lowerbounds.append(CURRENTTEMPERATURESCHEME[currenttempindex])
            CURRENTTEMPERATURESCHEME[currenttempindex] = 10**((log(CURRENTTEMPERATURESCHEME[currenttempindex], 10) +  log(min(upperbounds), 10) )/2 )
        else:
            upperbounds = [CURRENTTEMPERATURESCHEME[1]]
            CURRENTTEMPERATURESCHEME.insert(1, (CURRENTTEMPERATURESCHEME[0] + CURRENTTEMPERATURESCHEME[1]) / 2)
            lowerbounds = [1.0]

    removefile(os.getcwd() +os.sep + temporaryfoldername +os.sep + "covariance_and_multiplier.txt")
    removefile(os.getcwd() +os.sep + temporaryfoldername + os.sep  + "temp_input.txt")
    removefile(os.getcwd() +os.sep + temporaryfoldername + os.sep + "variance_correction.txt")
    if os.path.exists(os.getcwd() +os.sep + temporaryfoldername + os.sep  +"temp_adbayes" ):
        os.rmdir(os.getcwd() +os.sep + temporaryfoldername + os.sep  + "temp_adbayes" )
    if os.path.exists(os.getcwd() +os.sep + temporaryfoldername):
        os.rmdir(os.getcwd() + os.sep + temporaryfoldername)

if __name__=='__main__':
    import sys
    main(sys.argv[1:])
