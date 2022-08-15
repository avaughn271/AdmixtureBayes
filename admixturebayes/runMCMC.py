from argparse import ArgumentParser, SUPPRESS

from construct_starting_trees_choices import get_starting_trees
from construct_covariance_choices import get_covariance, estimate_degrees_of_freedom_scaled_fast
from posterior import posterior_class
from MCMCMC import MCMCMC
from MCMC import basic_chain
import os
from numpy import random
import pandas
from meta_proposal import simple_adaptive_proposal

import summary
import Rtree_operations
import tree_statistics
import Rtree_to_covariance_matrix
from copy import deepcopy

def get_summary_scheme(no_chains=1):
    
    summaries=[summary.s_posterior(),
               summary.s_likelihood(),
               summary.s_prior(),
               summary.s_no_admixes(),
               summary.s_variable('add', output='double'), 
               summary.s_total_branch_length(),
               summary.s_basic_tree_statistics(Rtree_operations.get_number_of_ghost_populations, 'ghost_pops', output='integer'),
               summary.s_basic_tree_statistics(Rtree_to_covariance_matrix.get_populations_string, 'descendant_sets', output='string')]
    summaries.append(summary.s_basic_tree_statistics(tree_statistics.unique_identifier_and_branch_lengths, 'tree', output='string'))
    summaries.append(summary.s_basic_tree_statistics(tree_statistics.get_admixture_proportion_string, 'admixtures', output='string'))
    sample_verbose_scheme={summary.name:(1,0) for summary in summaries}
    sample_verbose_scheme_first=deepcopy(sample_verbose_scheme)
    if no_chains==1:
        return [sample_verbose_scheme_first], summaries
    else:
        return [sample_verbose_scheme_first]+[{}]*(no_chains-1), summaries

class fixed_geometrical(object):

    def __init__(self, maxT, no_chains):
        if no_chains==1:
            self.temps=[1.0]
        else:
            self.temps=[maxT**(float(i)/(float(no_chains)-1.0)) for i in range(no_chains)]

    def get_temp(self,i):
        return self.temps[i]

def get_nodes(input_file, reduce_node):
    nodes=read_one_line(input_file)
    reduced_nodes=deepcopy(nodes)
    reduced_nodes.remove(reduce_node)
    return nodes, reduced_nodes

def read_one_line(filename):
    with open(filename, 'r') as f:
        return f.readline().rstrip().split()

def main(args):
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"

    parser = ArgumentParser(usage='pipeline for Admixturebayes')

    #input/output options
    parser.add_argument('--input_file', type=str, required=True, help='the input file of the pipeline. It should be of the same type as the treemix input file with a header of population names and each line representing a snp (unless --covariance_pipeline is altered).')
    parser.add_argument('--result_file', type=str, default='mcmc_samples.csv', help='file in which to save results.')

    parser.add_argument('--outgroup', type=str, default='',
                        help='The name of the population that should be outgroup for the covariance matrix. If the covariance matrix is supplied at stage 8 , this argument is not needed.')

    #Important arguments
    parser.add_argument('--MCMC_chains', type=int, default=8,
                        help='The number of chains to run the MCMCMC with. Optimally, the number of cores matches the number of chains.')
    parser.add_argument('--n', type=int, default=200, help='the number of MCMCMC flips throughout the chain.')
    parser.add_argument('--bootstrap_blocksize', type=int, default=1000,
                        help='the size of the blocks to bootstrap in order to estimate the degrees of freedom in the wishart distribution')

    #medium important convergence arguments
    parser.add_argument('--m', type=int, default=50, help='the number of MCMC steps before between each MCMCMC flip')

    #convenience arguments
    parser.add_argument('--verbose_level', default='normal', choices=['normal', 'silent'],
                        help='this will set the amount of status out prints throughout running the program.')
    parser.add_argument('--thinning_coef', type=int, default=40,
                        help=SUPPRESS)#'The number of MCMC steps between each saved instance. It has to be lower than --m.')

    #more obscure arguments
    parser.add_argument('--p', type=float, default=0.5,
                        help='the geometrical parameter in the prior. The formula is p**x(1-p)')
    parser.add_argument('--no_bootstrap_samples', type=int, default=100,
                        help='the number of bootstrap samples to make to estimate the degrees of freedom in the wishart distribution.')
    #start arguments
    parser.add_argument('--continue_samples', type=str, nargs='+', default=[],
                        help='filenames of trees to start in. If empty, the trees will either be simulated with the flag --random_start')

    options=parser.parse_args(args)

    assert not (any((i < 8 for i in [6,8,9])) and not options.outgroup), 'In the requested analysis, the outgroup needs to be specified by the --outgroup flag and it should match one of the populations'

    #Here is the only thing we should be changing.
    temp = pandas.read_csv(options.input_file, sep =" ")
    colnames = list(temp.columns.values)
    temp = temp[sorted(colnames)]
    temp.to_csv(os.getcwd() + "/temp_input.txt", sep =" ", index = False)

    mp= [simple_adaptive_proposal(['deladmix', 'addadmix', 'rescale', 'rescale_add', 'rescale_admixtures', 'rescale_constrained', 'sliding_regraft'],
     [1, 1, 1, 1, 1, 1, 1]) for _ in range(options.MCMC_chains)]

    full_nodes, reduced_nodes=get_nodes(os.getcwd() + "/temp_input.txt", options.outgroup)

    treemix_in_file=os.getcwd() + "/temp_input.txt"

    estimator_arguments=dict(reducer=options.outgroup, 
                             nodes=full_nodes,
                             Simulator_fixed_sxeed=True,
                             add_variance_correction_to_graph=True,
                             save_variance_correction=True)
                             
    covariance=get_covariance(os.getcwd() + "/temp_input.txt",
                              full_nodes=full_nodes,
                              reduce_covariance_node=options.outgroup,
                              estimator_arguments=estimator_arguments)

    estimator_arguments['save_variance_correction']=False
    df, boot_covs=estimate_degrees_of_freedom_scaled_fast(treemix_in_file,
                                            bootstrap_blocksize=options.bootstrap_blocksize,
                                            no_bootstrap_samples=options.no_bootstrap_samples,
                                            cores=options.MCMC_chains,
                                            est=estimator_arguments, 
                                            verbose_level=options.verbose_level)

    multiplier=covariance[1]
        
    tree_nodes=reduced_nodes
    if options.continue_samples != []:
        #This is where the continuation is all happening.
        #We first save the tree to a temporary file
        if os.path.exists(os.getcwd() + "/" +  "temp_start_tree.txt"):
            os.remove(os.getcwd() + "/" +  "temp_start_tree.txt")
        temp = pandas.read_csv(os.getcwd() + "/" + (options.continue_samples[0]))
        temp2 = (temp[["tree"]])
        f = open(os.getcwd() + "/" +  "temp_start_tree.txt", "a")
        f.write("\n")
        f.write(temp2.iloc[len(temp2.index) - 1, 0])
        f.write("\n")
        f.close()

        if os.path.exists(os.getcwd() + "/" +  "temp_starttree.txt"):
            os.remove(os.getcwd() + "/" +  "temp_starttree.txt")
        gii = open(os.getcwd() + "/" + "covariance_and_multiplier.txt", "r")
        g = gii.readlines()
        g = g[len(g)-1]
        g = g.split("=")
        multiplier = float(g[len(g)-1])
        gii.close()
        fff = open(os.getcwd() + "/" + "temp_start_tree.txt", "r")

        FinalString = "\n"

        f = fff.readlines()
        secondline = f[1]
        splitted = secondline.split(";")

        FinalString = FinalString + splitted[0] + ";"

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
        gg = open(os.getcwd() + "/" +  "temp_starttree.txt", "a")
        gg.write(FinalString)
        gg.close()

        temp = pandas.read_csv(os.getcwd() + "/" + (options.continue_samples[0]))
        addvalue = (temp[["add"]])
        addvalue = float(addvalue.iloc[len(addvalue.index) - 1, 0])

        if os.path.exists(os.getcwd() + "/" + "temp_add.txt"):
            os.remove(os.getcwd() + "/" + "temp_add.txt")
        f = open(os.getcwd() + "/" + "temp_add.txt", "a")
        f.write(str(addvalue) + "\n")
        f.close()

        #compute addfile as a file with a number
        starting_trees=get_starting_trees([os.getcwd() + "/" +  "temp_starttree.txt"],
                                        options.MCMC_chains,
                                        adds=[os.getcwd() + "/" + "temp_add.txt"],
                                        nodes=tree_nodes)
    else:
        starting_trees=get_starting_trees(options.continue_samples,
                                      options.MCMC_chains,
                                      adds=[],
                                      nodes=tree_nodes)

    summary_verbose_scheme, summaries=get_summary_scheme(no_chains=options.MCMC_chains)

    sim_lengths=[options.m]*options.n

    likelihood_nodes=reduced_nodes
    posterior = posterior_class(emp_cov=covariance[0],
                                M=df,
                                p=options.p,
                                multiplier=covariance[1],
                                nodes=likelihood_nodes)

    posterior_function_list=[]

    temperature_scheme=fixed_geometrical(1000,options.MCMC_chains)

    #####    ANDREW DEBUG   !!!!!!!!!!
    if os.path.exists("covariance_without_reduce_name.txt"):
        os.remove("covariance_without_reduce_name.txt")
    else:
        print("The file does not exist")

    if os.path.exists("variance_correction.txt"):
        os.remove("variance_correction.txt")
    else:
        print("The file does not exist")

    if os.path.exists("temp_starttree.txt"):
        os.remove("temp_starttree.txt")
    if os.path.exists("temp_start_tree.txt"):
        os.remove("temp_start_tree.txt")
    if os.path.exists("temp_add.txt"):
        os.remove("temp_add.txt")

    if os.path.exists(os.getcwd() + "/temp_adbayes"):
        os.rmdir(os.getcwd() + "/temp_adbayes")

    def multi_chain_run():
        #random_seeds = []
        #for i in range(options.MCMC_chains):
        #    random_seeds.append(givenseed + i)
        #print(random_seeds)
        res=MCMCMC(starting_trees=starting_trees,
               posterior_function= posterior,
               summaries=summaries,
               temperature_scheme=temperature_scheme,
               printing_schemes=summary_verbose_scheme,
               iteration_scheme=sim_lengths,
               overall_thinnings=int(options.thinning_coef),
               proposal_scheme= mp,
               cores=options.MCMC_chains,
               no_chains=options.MCMC_chains,
               multiplier=multiplier,  #numpy_seeds = random_seeds,
               result_file=options.result_file,
               posterior_function_list=posterior_function_list, n_arg=options.n, m_arg=options.m, verboseee=options.verbose_level)
    def single_chain_run():
        basic_chain(start_x= starting_trees[0],
                    summaries=summaries,
                    posterior_function=posterior,
                    proposal=mp[0],
                    post=None,
                    N=sum(sim_lengths),
                    sample_verbose_scheme=summary_verbose_scheme[0],
                    overall_thinning=int(options.thinning_coef),
                    i_start_from=0,
                    temperature=1.0,
                    proposal_update=None,
                    multiplier=multiplier,
                    appending_result_file=options.result_file,
                    appending_result_frequency=sim_lengths[0])

    if os.path.exists(os.getcwd() + "/temp_input.txt"):
        os.remove(os.getcwd() + "/temp_input.txt")

    if options.MCMC_chains==1:
        single_chain_run()
    else:
        res=multi_chain_run()
    if os.path.exists("trees_tmp.txt"):
        os.remove("trees_tmp.txt")
    if options.continue_samples != []:
        oldcsv = pandas.read_csv(os.getcwd() + "/" + (options.continue_samples[0]))
        newcsv = pandas.read_csv(os.getcwd() + "/" + options.result_file)
        result = pandas.concat([oldcsv,newcsv])
        if os.path.exists(os.getcwd() + "/" + options.result_file):
            os.remove(os.getcwd() + "/" + options.result_file)
        result.to_csv(os.getcwd() + "/" + options.result_file, index = False)

if __name__=='__main__':
    import sys
    main(sys.argv[1:])
