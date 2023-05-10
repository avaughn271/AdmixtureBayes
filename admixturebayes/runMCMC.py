from argparse import ArgumentParser, SUPPRESS

import construct_starting_trees_choices
from construct_covariance_choices import get_covariance, estimate_degrees_of_freedom_scaled_fast
from posterior import posterior_class
from MCMCMC import MCMCMC
import os
from numpy import random
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
    parser.add_argument('--MCMC_chains', type=int, default=8,
                        help='The number of chains to run the MCMCMC with. Optimally, the number of cores matches the number of chains.')
    parser.add_argument('--n', type=int, default=200, help='the number of MCMCMC flips throughout the chain.')
    parser.add_argument('--bootstrap_blocksize', type=int, default=1000,
                        help='the size of the blocks to bootstrap in order to estimate the degrees of freedom in the wishart distribution')

    #convenience arguments
    parser.add_argument('--verbose_level', default='normal', choices=['normal', 'silent'],
                        help='this will set the amount of status out prints throughout running the program.')
    #start arguments
    parser.add_argument('--continue_samples', type=str, nargs='+', default=[],
                        help='filenames of trees to start in. If empty, the trees will either be simulated with the flag --random_start')
    
    parser.add_argument('--maxtemp', type=float, default=1000.0, help='the max temp of the hottest chain')

    options=parser.parse_args(args)

    temporaryfoldername = (options.result_file).replace('.', '') + "_tempfilefolder"
    os.mkdir(os.getcwd() + os.sep + temporaryfoldername)

    assert options.MCMC_chains > 1, 'At least 2 chains must be run for the MCMCMC to work properly'
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
     [1, 1, 1, 1, 1, 1, 1]) for _ in range(options.MCMC_chains)]

    with open(os.getcwd() +os.sep+ temporaryfoldername + os.sep + "temp_input.txt", 'r') as f:
        full_nodes = f.readline().rstrip().split()
    reduced_nodes=deepcopy(full_nodes)
    reduced_nodes.remove(options.outgroup)

    estimator_arguments=dict(reducer=options.outgroup, nodes=full_nodes, add_variance_correction_to_graph=True, save_variance_correction=True)
                             
    covariance=get_covariance(os.getcwd() +os.sep+ temporaryfoldername + os.sep + "temp_input.txt", 
    varcovfilename = os.getcwd() +os.sep+ temporaryfoldername + os.sep + "variance_correction.txt",
    full_nodes=full_nodes,
     reduce_covariance_node=options.outgroup,
     estimator_arguments=estimator_arguments,filename =  os.getcwd() +os.sep + temporaryfoldername + os.sep  +"covariance_and_multiplier.txt")

    estimator_arguments['save_variance_correction']=False
    df=estimate_degrees_of_freedom_scaled_fast(os.getcwd() + os.sep + temporaryfoldername + os.sep  + "temp_input.txt",
                                               varcovfilename = os.getcwd() +os.sep + temporaryfoldername + os.sep  +"variance_correction.txt",
                                            bootstrap_blocksize=options.bootstrap_blocksize,
                                            cores=options.MCMC_chains,
                                            est=estimator_arguments, 
                                            verbose_level=options.verbose_level)

    multiplier=covariance[1]
    
    if options.continue_samples != []:
        #This is where the continuation is all happening.
        #We first save the tree to a temporary file
        removefile(os.getcwd() + os.sep +  "temp_start_tree.txt")
        temp = pandas.read_csv(os.getcwd() + os.sep + (options.continue_samples[0]))
        temp2 = (temp[["tree"]])
        f = open(os.getcwd() + os.sep  +  "temp_start_tree.txt", "a")
        f.write("\n")
        f.write(temp2.iloc[len(temp2.index) - 1, 0])
        f.write("\n")
        f.close()

        removefile(os.getcwd() + os.sep  + "temp_starttree.txt")
        gii = open(os.getcwd() + os.sep  + temporaryfoldername + os.sep  + "covariance_and_multiplier.txt", "r")
        g = gii.readlines()
        g = g[len(g)-1]
        g = g.split("=")
        multiplier = float(g[len(g)-1])
        gii.close()
        fff = open(os.getcwd() + os.sep  + "temp_start_tree.txt", "r")

        f = fff.readlines()
        secondline = f[1]
        splitted = secondline.split(";")

        FinalString = "\n" + splitted[0] + ";"

        relevantbranches = splitted[1]
        splitbranches = relevantbranches.split("-")

        #This part is new
        while(True):
            for i in range(len(splitbranches) - 1):
                currentstring = splitbranches[i]
                nexstring = splitbranches[i+1]
                if currentstring[len(currentstring) - 1] == "e":
                    splitbranches[i] = currentstring + "-" + nexstring
                    splitbranches.pop(i + 1)
                    break
            break
        for i in range(len(splitbranches)):
            splitbranches[i] = "{:.9f}".format(float(splitbranches[i]))

        for i in splitbranches:
            FinalString = FinalString + "{:.12f}".format(float(i) * multiplier) + "-"

        FinalString = FinalString[0:(len(FinalString)-1)]
        FinalString = FinalString + ";" + splitted[2]
        fff.close()
        gg = open(os.getcwd() +os.sep  +  "temp_starttree.txt", "a")
        gg.write(FinalString)
        gg.close()

        temp = pandas.read_csv(os.getcwd() + os.sep + (options.continue_samples[0]))
        addvalue = (temp[["add"]])
        addvalue = float(addvalue.iloc[len(addvalue.index) - 1, 0])

        removefile(os.getcwd() + os.sep + "temp_add.txt")
        f = open(os.getcwd() + os.sep  + "temp_add.txt", "a")
        f.write(str(addvalue) + "\n")
        f.close()

        #compute addfile as a file with a number
        starting_trees=construct_starting_trees_choices.get_starting_trees([os.getcwd() + os.sep +  "temp_starttree.txt"],
                                        options.MCMC_chains,
                                        adds=[os.getcwd() + os.sep  + "temp_add.txt"],
                                        nodes=reduced_nodes)
    else:
        starting_trees=construct_starting_trees_choices.get_starting_trees(options.continue_samples, options.MCMC_chains, adds=[], nodes=reduced_nodes)

    summary_verbose_scheme, summaries=get_summary_scheme(no_chains=options.MCMC_chains)

    posterior = posterior_class(emp_cov=covariance[0], M=df, multiplier=covariance[1], nodes=reduced_nodes, 
                                 varcovname=os.getcwd() +os.sep + temporaryfoldername + os.sep  + "variance_correction.txt")

    removefile("temp_starttree.txt")
    removefile("temp_start_tree.txt")
    removefile("temp_add.txt")

    if options.save_covariance:
        removefile(os.getcwd() + os.sep  +"covariance_matrix.txt")
        Liness = open(os.getcwd() +os.sep + temporaryfoldername +os.sep  + "covariance_and_multiplier.txt", 'r').readlines()
        covarfile = open(os.getcwd() + os.sep  +"covariance_matrix.txt", "a")
        covarfile.writelines(Liness)
        covarfile.close()

        #random_seeds = []
        #for i in range(options.MCMC_chains):
        #    random_seeds.append(givenseed + i)
        #print(random_seeds)
    MCMCMC(starting_trees=starting_trees,
            posterior_function= posterior,
            summaries=summaries,
            temperature_scheme=[(options.maxtemp)**(float(i)/(float(options.MCMC_chains)-1.0)) for i in range(options.MCMC_chains)],  ##[1.0 + 0.1 * float(i) for i in range(options.MCMC_chains)]
            printing_schemes=summary_verbose_scheme,
            iteration_scheme=[50]*options.n,
            proposal_scheme= mp,
            no_chains=options.MCMC_chains,
            multiplier=multiplier,  #numpy_seeds = random_seeds,
            result_file=options.result_file,
            n_arg=options.n, verboseee=options.verbose_level)

    if options.continue_samples != []:
        oldcsv = pandas.read_csv(os.getcwd() + os.sep  + (options.continue_samples[0]))
        newcsv = pandas.read_csv(os.getcwd() + os.sep + options.result_file)
        result = pandas.concat([oldcsv,newcsv])
        removefile(os.getcwd() + os.sep + options.result_file)
        result.to_csv(os.getcwd() +os.sep + options.result_file, index = False)
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
