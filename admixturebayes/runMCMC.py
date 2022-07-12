from argparse import ArgumentParser, SUPPRESS

from construct_proposal_choices import make_proposal
from construct_starting_trees_choices import get_starting_trees
from construct_covariance_choices import get_covariance
from construct_summary_choices import get_summary_scheme
from construct_filter_choices import make_filter
from posterior import posterior_class
from MCMCMC import MCMCMC
from wishart_distribution_estimation import estimate_degrees_of_freedom_scaled_fast
from MCMC import basic_chain
from stop_criteria import stop_criteria
import os
from numpy import random
import pandas

class fixed_geometrical(object):

    def __init__(self, maxT, no_chains):
        if no_chains==1:
            self.temps=[1.0]
        else:
            self.temps=[maxT**(float(i)/(float(no_chains)-1.0)) for i in range(no_chains)]

    def get_temp(self,i):
        return self.temps[i]

    def update_temps(self, permut):
        pass


from copy import deepcopy
from Rtree_operations import get_trivial_nodes

def get_nodes(arguments, input_file, outgroup_name, reduce_node, backup_number=8):
    ''' The outgroup_name is only used for simulation purposes and reduce_node is the important one
    that is used when analysing the admixture graphs.  '''
    if not arguments[0]:#this means that we should use the input file for nodes
        if ';' in input_file:
            nodes=get_trivial_nodes(len(input_file.split('-')[0].split('.')))
        elif '.' in input_file:
            nodes=read_one_line(input_file)
        elif ',' in input_file:
            nodes=get_trivial_nodes(int(input_file[1:].split(',')[0]))
        else:
            nodes=get_trivial_nodes(int(input_file))
    else:
        nodes=arguments
    before_added_outgroup=deepcopy(nodes)
    reduced_nodes=deepcopy(nodes)
    if outgroup_name in nodes:
        before_added_outgroup.remove(outgroup_name)
    elif outgroup_name:
        nodes.append(outgroup_name)
    if reduce_node in reduced_nodes:
        reduced_nodes.remove(reduce_node)
    if reduce_node and reduce_node not in nodes:
        nodes.append(reduce_node)
    return before_added_outgroup, nodes, reduced_nodes

from tree_to_data import unzip

def read_one_line(filename):
    if filename.endswith('.gz'):
        filename=unzip(filename)
    with open(filename, 'r') as f:
        return f.readline().rstrip().split()

def main(args):
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"

    parser = ArgumentParser(usage='pipeline for Admixturebayes')

    #input/output options
    parser.add_argument('--input_file', type=str, required=True, help='the input file of the pipeline. It should be of the same type as the treemix input file with a header of population names and each line representing a snp (unless --covariance_pipeline is altered).')
    parser.add_argument('--result_file', type=str, default='mcmc_samples.csv', help='file in which to save results. The prefix will not be prepended the result_file.')

    parser.add_argument('--outgroup', type=str, default='',
                        help='The name of the population that should be outgroup for the covariance matrix. If the covariance matrix is supplied at stage 8 , this argument is not needed.')

    #Important arguments
    parser.add_argument('--MCMC_chains', type=int, default=8,
                        help='The number of chains to run the MCMCMC with. Optimally, the number of cores matches the number of chains.')
    parser.add_argument('--n', type=int, default=200, help='the number of MCMCMC flips throughout the chain.')
    parser.add_argument('--df_file', type=str, default='',
                        help='By default, the degrees of freedom will be estimated with bootstrap. If this flag is used, it will cancel the bootstrap estimation of the degrees of freedom. The degrees of freedom represents the number of effectively independent SNPs in the dataset.')
    parser.add_argument('--wishart_df', type=float, default=-1,
                        help='By default, the degrees of freedom will be estimated with bootstrap. If this flag is used (and is positive), it will cancel the bootstrap estimation of the degrees of freedom.')
    parser.add_argument('--bootstrap_blocksize', type=int, default=1000,
                        help='the size of the blocks to bootstrap in order to estimate the degrees of freedom in the wishart distribution')


    #medium important convergence arguments
    parser.add_argument('--m', type=int, default=50, help='the number of MCMC steps before between each MCMCMC flip')
    parser.add_argument('--max_temp', type=float, default=1000, help='the maximum temperature used in the MCMCMC.')
    parser.add_argument('--adaptive_temperatures', action='store_true', default=False,
                        help='this will make the MCMCMC temperature scheme update itself based on the transition probabilities.')
    parser.add_argument('--stop_criteria', action='store_true', default=False,
                        help='If applied the MCMCMC will stop when the coldest chain has an effective sample size at a certain threshold for a number of different summaries.')
    parser.add_argument('--stop_criteria_frequency', type=int, default=200000,
                        help='The frequency of checking for when the stop criteria are checked (if the stop_criteria flag is turned on). It is measured in total iterations(n*m).')
    parser.add_argument('--stop_criteria_continuous_ess_threshold', default=200, type=float,
                        help='The minimum ESS to obtain for continuous summaries of the MCMC chain')
    parser.add_argument('--stop_criteria_threshold', default=200, type=float,
                        help='The minimum ESS to obtain for topological summaries of the MCMC chain (if the stop_criteria). If negative, the topological stop criteria will not be used at all.')

    #convenience arguments
    parser.add_argument('--prefix', type=str, default='',
                        help='this directory will be the beginning of every temporary file created in the covariance pipeline and in the estimation of the degrees of freedom in the wishart distribution.')
    parser.add_argument('--nodes', type=str, nargs='+', default=[''],
                        help='list of nodes of the populations or the filename of a file where the first line contains all population names. If unspecified the first line of the input_file will be used. If no input file is found, there will be used standard s1,..,sn.')
    parser.add_argument('--verbose_level', default='normal', choices=['normal', 'silent'],
                        help='this will set the amount of status out prints throughout running the program.')
    parser.add_argument('--Rscript_command', default='Rscript', type=str,
                        help='The command to start R from the terminal. If there is no valid path, the stop criteria should be not used. Its default value is "Rscript"')
    parser.add_argument('--save_warm_chains', action='store_true', default=False,
                        help='By default only the coldest, "real" chain in the MCMCMC is saved. This will save all of them.')
    parser.add_argument('--thinning_coef', type=int, default=40,
                        help=SUPPRESS)#'The number of MCMC steps between each saved instance. It has to be lower than --m.')

    #more obscure arguments
    parser.add_argument('--outgroup_type', choices=['non_admixed','Free'], default='non_admixed',
                        help='The type of outgroup in the model. If non_admixed, there will be no admixture events'
                             'into the outgroup. If Free there may be admixtures into the outgroup. '
                             'The outgroup still need to be specified to place the root.')
    parser.add_argument('--covariance_pipeline', nargs='+', type=int, default=[6, 8, 9],
                        help='The list of steps the data should go through to become a covariance matrix. For the simulation data steps 1-5 I refer to the script construct_covariance_choices.py.'
                             '6=data in treemix format (see Readme.md), '
                             '7=raw covariance with outgroup, '
                             '8=raw covariance without outgroup (or with respect to the outgroup).'
                             '9=scaled covariance without outgroup and the scaling factor.')
    parser.add_argument('--variance_correction', default='unbiased', choices=['None', 'unbiased', 'mle'],
                        help='The type of adjustment used on the empirical covariance.')
    parser.add_argument('--unadmixed_populations', default=[], nargs='*',
                        help='This sets the prior to 0 for all graphs where the lineage of any of the supplied populations experience an admixture. WARNING: this will change the prior on the number of admixture events.')
    parser.add_argument('--cov_estimation',
                        choices=['None', 'Jade', 'outgroup_sum', 'outgroup_product', 'average_sum', 'average_product',
                                 'Jade-o', 'EM'], default='average_sum',
                        help=SUPPRESS)#'this is the way of estimating the empirical covariance matrix.')
    parser.add_argument('--Jade_cutoff', type=float, default=1e-5,
                        help=SUPPRESS)#'this will remove SNPs of low diversity in either the Jade or the Jade-o scheme.')
    parser.add_argument('--scale_goal', choices=['min', 'max'], default='max',
                        help='If 9 is included in the pipeline, this is the rescaling of the covariance matrix.')
    parser.add_argument('--p', type=float, default=0.5,
                        help='the geometrical parameter in the prior. The formula is p**x(1-p)')
    parser.add_argument('--sap_analysis', action='store_true', default=False,
                        help='skewed admixture proportion prior in the analysis')
    parser.add_argument('--not_uniform_prior', action='store_true', default=False,
                        help='If applied the uniform prior will not be used on the topology conditioned on the number of admixture events. Instead the ')
    #parser.add_argument('--no_add', action='store_true', default=False, help='this will remove the add contribution')
    parser.add_argument('--no_bootstrap_samples', type=int, default=100,
                        help='the number of bootstrap samples to make to estimate the degrees of freedom in the wishart distribution.')
    parser.add_argument('--save_bootstrap_covariances', type=str, default='', help='if provided the bootstrapped covariance matrices will be saved to numbered files starting with {prefix}+_+{save_covariances}+{num}+.txt')
    parser.add_argument('--bootstrap_type_of_estimation', choices=['mle_opt','var_opt'], default='var_opt', help='This is the way the bootstrap wishart estimate is estimated.')
    parser.add_argument('--load_bootstrapped_covariances', type=str, default=[], nargs='+', help='if supplied, this will load covariance matrices from the specified files instead of bootstrapping new ones.')
    # proposal frequency options
    parser.add_argument('--deladmix', type=float, default=1, help='this states the frequency of the proposal type')
    parser.add_argument('--addadmix', type=float, default=1, help='this states the frequency of the proposal type')
    parser.add_argument('--rescale', type=float, default=1, help='this states the frequency of the proposal type')
    parser.add_argument('--regraft', type=float, default=0, help='this states the frequency of the proposal type')
    parser.add_argument('--rescale_add', type=float, default=1, help='this states the frequency of the proposal type')
    parser.add_argument('--rescale_admix', type=float, default=1, help='this states the frequency of the proposal type')
    parser.add_argument('--rescale_admix_correction', type=float, default=0,
                        help='this states the frequency of the proposal type')
    parser.add_argument('--rescale_constrained', type=float, default=1,
                        help='this states the frequency of the proposal type')
    parser.add_argument('--rescale_marginally', type=float, default=0,
                        help='this states the frequency of the proposal type')
    parser.add_argument('--sliding_regraft', type=float, default=1,
                        help='this states the frequency of the proposal type')
    parser.add_argument('--sliding_rescale', type=float, default=0,
                        help='this states the frequency of the proposal type')
    parser.add_argument('--cancel_preserve_root_distance', default=False, action='store_true',
                        help="if applied there will not be made correction for root distance when adding and deleting admixtures")
    #start arguments
    parser.add_argument('--continue_samples', type=str, nargs='+', default=[],
                        help='filenames of trees to start in. If empty, the trees will either be simulated with the flag --random_start or the so-called trivial tree')
    parser.add_argument('--starting_adds', type=str, nargs='+', default=[],
                        help="filename of the adds to use on the starting trees.")
    parser.add_argument('--start', choices=['trivial', 'random', 'perfect'], default='trivial',
                        help='Where to start the chain - works only if starting trees are not specified.')
    parser.add_argument('--starting_tree_scaling',
                        choices=['None', 'empirical_trace', 'starting_tree_trace', 'scalar', 'treemix_tree'],
                        default='None', type=str,
                        help='The starting tree can be scaled as the covariance (as_covariance) or as the p')
    parser.add_argument('--starting_tree_use_scale_tree_factor', default=False, action='store_true',
                        help='this will scale the tree with the specified scale_tree_factor.')
    parser.add_argument('--mscale_file', default='', type=str,
                        help=SUPPRESS)#'This is the file where the normalization factor used by admixtureBayes are. '
                             #'This is normally calculated by the program but if settings have been changed, '
                             #'it may not and then this option can be used such that unnormalized '
                             #'treemix output trees can be scaled correctly')
    parser.add_argument('--rs', action='store_true', default=False, help='will change the prior on the number of admixture events in the tree.')
    parser.add_argument('--r_scale', type=float, default=1.0, help='This will set the the mean of the number of admixture events in chain i to 1+i*r')

    #more obscure convenience arguments
    parser.add_argument('--save_df_file', type=str, default='DF.txt',
                        help=SUPPRESS)#'the prefix is put before this string and the degrees of freedom is saved to this file.')
    parser.add_argument('--summary_majority_tree', action='store_true', default=False,
                        help='this will calculate the majority (newick) tree based on the sampled tree')
    parser.add_argument('--summary_acceptance_rate', action='store_true', default=True,
                        help=SUPPRESS)#'This will calculate and store summaries related to the acceptance rate')
    parser.add_argument('--summary_admixture_proportion_string', action='store_true', default=True,
                        help=SUPPRESS)#'this will save a string in each step indicating names and values of all admixture proportions')
    parser.add_argument('--store_permuts', action='store_true', default=False,
                        help='If applied, the permutations from the MCMCMC flips are recorded in a file with a similar filename to the result_file')
    parser.add_argument('--save_after_hours', type=float, nargs='+', default=[],
                        help=SUPPRESS)#'This will save a copy of the output file after the number of hours specified here. One would do that to easily access how converged the chain is after certain number of hours.')
    parser.add_argument('--profile', action='store_true', default=False,
                        help=SUPPRESS)#"this will embed the MCMC part in a profiler")
    parser.add_argument('--prior_run', action='store_true', default=False,
                        help=SUPPRESS)#"'Run the MCMC without likelihood and only the prior ')
    parser.add_argument('--evaluate_likelihood', action='store_true', default=False,
                        help=SUPPRESS)#'this will evaluate the likelihood in the starting tree and then stop, writing just a single file with three values, prior, likelihood and posterior.')
    parser.add_argument('--evaluate_bootstrap_likelihoods', action='store_true', default=False,
                        help=SUPPRESS)#'If evaluate likelihood is turned on this will calculate the likelihood of all bootstrapped covariances(if bootstrapping is also turned on)')
    parser.add_argument('--stop_evaluations', action='store_true', default=False,
                        help='This will stop the analysis after the data preparation')
    parser.add_argument('--variance_correction_input_file', default='', type=str,
                        help='if the variance correction is saved in a file (with numpy.savetxt format of a 2 dimensional numpy array) it can be loaded in with this command')

    #Very obscure arguments
    # tree simulation
    parser.add_argument('--p_sim', type=float, default=.5,
                        help=SUPPRESS)#'the parameter of the geometric distribution in the distribution to simulate the true tree from.')
    parser.add_argument('--popsize', type=int, default=20,
                        help=SUPPRESS)#"'the number of genomes sampled from each population.')
    parser.add_argument('--nreps', type=int, default=50,
                        help=SUPPRESS)#'How many pieces of size 500 kb should be simulated')
    parser.add_argument('--scale_tree_factor', type=float, default=0.02,
                        help=SUPPRESS)#""'The scaling factor of the simulated trees to make them less vulnerable to the fixation effect.')
    parser.add_argument('--skewed_admixture_prior_sim', default=False, action='store_true',
                        help=SUPPRESS)#'the prior tree is simulated with an uneven prior on the admixture proportions')
    parser.add_argument('--time_adjusted_tree', default=False, action='store_true',
                        help=SUPPRESS)#'this will modify the simulated tree such that all drift lengths from root to leaf are the same')
    parser.add_argument('--sadmix_tree', default=False, action='store_true',
                        help=SUPPRESS)#'this will simulate trees where all admixture events are important in the sense that they expand the space of possible covariance matrices.')
    parser.add_argument('--wishart_noise', action='store_true', default=False,
                        help=SUPPRESS)#'A wishart noise is added to the estimated covariance matrix.')
    parser.add_argument('--create_outgroup', type=str, default='',
                        help=SUPPRESS)#'The name of the outgroup that should be added to a simulated dataset.')
    # covariance simulation
    parser.add_argument('--favorable_init_brownian', default=False, action='store_true',
                        help=SUPPRESS)#'This will start the brownian motion(only if 21 in workflow) between 0.4 and 0.6')
    parser.add_argument('--unbounded_brownian', default=False, action='store_true',
                        help=SUPPRESS)#'This will start the brownian motion(only if 21 in workflow) between 0.4 and 0.6')
    parser.add_argument('--filter_on_outgroup', default=False, action='store_true',
                        help=SUPPRESS)#'If applied (and 23 in the pipeline) SNPs that are not polymorphic in the outgroup are removed. If not, the default is that polymorphic in no population are removed. ')
    parser.add_argument('--arcsin', action='store_true', default=False,
                        help=SUPPRESS)
    #other covariance matrices
    parser.add_argument('--bias_c_weight', choices=['default','None','outgroup_sum', 'outgroup_product', 'average_sum', 'average_product'], default='default',
                        help=SUPPRESS)#'from cov_weight with bias correction unweighted there are some obvious choices for weighing the bias correction, so here they are: None=None, Jade=average_sum, Jade-o=outgroup_sum, average_sum=average_sum, average_product=average_product, outgroup_sum=outgroup_sum, outgroup_product=outgroup_product')
    parser.add_argument('--add_variance_correction_to_graph', default=True, action='store_true',
                        help=SUPPRESS)#'If on, the variance correction will be added to the covariance matrix of the graph and not subtracted from the empirical covariance matrix. Default is True.')
    parser.add_argument('--indirect_correction', default=False, action='store_true',
                        help=SUPPRESS)#'the bias in the covariance is (possibly again) corrected for by indirect estimation.')
    parser.add_argument('--indirect_its', type=int, default=100,
                        help=SUPPRESS)#'For how many iterations should the indirect optimization procedure be run. Only applicable if indirect_correction is True')
    parser.add_argument('--indirect_simulation_factor', type=int, default=1,
                        help=SUPPRESS)#'How much more data than provided should be simulated in the indirect correction procedure. Only applicable if indirect_correction is True')
    parser.add_argument('--EM_maxits', type=int, default=100,
                        help=SUPPRESS)#'The maximum number of iterations of the EM algorithm. There is another stopping criteria that may stop it before.')
    parser.add_argument('--EM_alpha', type=float, default=1.0,
                        help=SUPPRESS)#'The EM algorithm assumes that the allele frequencies in the outgroup are known. In fact it is estimated with: if alpha=1.0: the empirical allele frequencies of the outgroup, alpha=0.0: the average empirical allele frequency in all the other populations.It can also be chosen as something in between. This estimator is biased because the outgroup allele frequencies are not known and because of extra normal distribution assumptions')
    parser.add_argument('--no_repeats_of_cov_est', type=int, default=1,
                        help=SUPPRESS)#'The number of times the simulation procedure should be run.')
    parser.add_argument('--initial_Sigma', choices=['default','random', 'start'], default='default',
                        help=SUPPRESS)#'This means that ')
    parser.add_argument('--filter_type', choices=['snp','none', 'outgroup_other','outgroup','all_pops'], default='snp',
                        help=SUPPRESS)#'This will apply a filter to positions based on their value.')
    parser.add_argument('--filter_on_simulated', choices=['same', 'none', 'outgroup_other', 'outgroup', 'snp', 'all_pops'], default='same',
                        help=SUPPRESS)#'In indirect inference, whole datasets are simulated under ')

    options=parser.parse_args(args)

    assert not (any((i < 8 for i in options.covariance_pipeline)) and not options.outgroup), 'In the requested analysis, the outgroup needs to be specified by the --outgroup flag and it should match one of the populations'


    no_add=options.outgroup_type=='None' or options.outgroup_type=='Free'

    mp = make_proposal(deladmix=options.deladmix,
                  addadmix=options.addadmix,
                  rescale=options.rescale,
                  regraft=options.regraft,
                  rescale_add=options.rescale_add,
                  rescale_admix=options.rescale_admix,
                  rescale_admix_correction=options.rescale_admix_correction,
                  rescale_constrained=options.rescale_constrained,
                  rescale_marginally=options.rescale_marginally,
                  sliding_regraft=options.sliding_regraft,
                  sliding_rescale=options.sliding_rescale,
                  MCMC_chains=options.MCMC_chains,
                  cancel_preserve_root_distance=options.cancel_preserve_root_distance,
                  no_add=no_add)

    before_added_outgroup, full_nodes, reduced_nodes=get_nodes(options.nodes, options.input_file, options.create_outgroup, options.outgroup)

    prefix=options.prefix

    if options.covariance_pipeline[0]==6:
        treemix_in_file=options.input_file
        treemix_file=options.input_file
    else:
        treemix_file=prefix+"treemix_in.txt"
        treemix_in_file=treemix_file+'.gz'

    preliminary_starting_trees=[None]
    assert options.initial_Sigma!='start', 'to make the filter start somewhere specific it should also be specified specifically'

    locus_filter=make_filter(options.filter_type,
                             outgroup_name=options.outgroup, #options.reduce_node,
                             covariance_pipeline=options.covariance_pipeline)
    if options.filter_on_simulated=='same':
        locus_filter_on_simulated=make_filter(options.filter_type)
    else:
        locus_filter_on_simulated=make_filter(options.filter_on_simulated)

    estimator_arguments=dict(reducer=options.outgroup, #options.reduce_node,
                             variance_correction=options.variance_correction,
                             method_of_weighing_alleles=options.cov_estimation,
                             arcsin_transform=options.arcsin,
                             jade_cutoff=options.Jade_cutoff,
                             bias_c_weight=options.bias_c_weight,
                             nodes=full_nodes,
                             indirect_correction=options.indirect_correction,
                             Indirect_its=options.indirect_its,
                             Indirect_multiplier_s=options.indirect_simulation_factor,
                             EM_maxits=options.EM_maxits,
                             EM_alpha=options.EM_alpha,
                             no_repeats_of_cov_est=options.no_repeats_of_cov_est,
                             Simulator_fixed_sxeed=True,
                             initial_Sigma_generator={options.initial_Sigma:(preliminary_starting_trees[0], reduced_nodes)},
                             locus_filter_on_simulated=locus_filter_on_simulated,
                             add_variance_correction_to_graph=options.add_variance_correction_to_graph,
                             save_variance_correction=True,
                             prefix=prefix)

    covariance=get_covariance(options.covariance_pipeline,
                              options.input_file,
                              full_nodes=full_nodes,
                              skewed_admixture_prior_sim=options.skewed_admixture_prior_sim,
                              p=options.p_sim,
                              outgroup_name=options.create_outgroup,
                              reduce_covariance_node=options.outgroup, #reduce_node,
                              sample_per_pop=options.popsize,
                              nreps=options.nreps,
                              treemix_file=treemix_file,
                              scale_tree_factor=options.scale_tree_factor,
                              prefix=prefix,
                              t_adjust_tree=options.time_adjusted_tree,
                              sadmix=options.sadmix_tree,
                              add_wishart_noise_to_covariance=options.wishart_noise,
                              df_of_wishart_noise_to_covariance=options.wishart_df,
                              scale_goal=options.scale_goal,
                              favorable_init_brownian=options.favorable_init_brownian,
                              unbounded_brownian=options.unbounded_brownian,
                              filter_on_outgroup=options.filter_on_outgroup,
                              locus_filter=locus_filter,
                              estimator_arguments=estimator_arguments,
                              verbose_level=options.verbose_level
                              )

    if options.df_file:
        with open(options.df_file, 'r') as f:
            df = float(f.readline().rstrip())
    elif options.wishart_df>0:
        df = options.wishart_df
    else:
        estimator_arguments['save_variance_correction']=False
        summarization=options.bootstrap_type_of_estimation
        df, boot_covs=estimate_degrees_of_freedom_scaled_fast(treemix_in_file,
                                               bootstrap_blocksize=options.bootstrap_blocksize,
                                               no_bootstrap_samples=options.no_bootstrap_samples,
                                               summarization=summarization,
                                               cores=options.MCMC_chains,
                                               save_covs=options.save_bootstrap_covariances,
                                               prefix=prefix,
                                               est=estimator_arguments, locus_filter=locus_filter,
                                               load_bootstrapped_covariances=options.load_bootstrapped_covariances,
                                               verbose_level=options.verbose_level)

    if 9 not in options.covariance_pipeline:
        multiplier=None
        covariance=(covariance, multiplier)
    else:
        multiplier=covariance[1]

    if not options.mscale_file:
        mscale_file=None
    else:
        mscale_file=options.mscale_file

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
        print(splitbranches)
        #End of new part

        for i in splitbranches:
            FinalString = FinalString + str(float(i) * multiplier) + "-"

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
                                        nodes=tree_nodes,
                                        pipeline=options.covariance_pipeline,
                                        multiplier=multiplier,
                                        scale_tree_factor=options.scale_tree_factor,
                                        start=options.start,
                                        prefix=prefix,
                                        starting_tree_scaling=options.starting_tree_scaling,
                                        starting_tree_use_scale_tree_factor=options.starting_tree_use_scale_tree_factor,
                                        scale_goal=options.scale_goal,
                                        mscale_file=mscale_file,
                                        no_add=no_add)
    else:
        starting_trees=get_starting_trees(options.continue_samples,
                                      options.MCMC_chains,
                                      adds=options.starting_adds,
                                      nodes=tree_nodes,
                                      pipeline=options.covariance_pipeline,
                                      multiplier=multiplier,
                                      scale_tree_factor=options.scale_tree_factor,
                                      start=options.start,
                                      prefix=prefix,
                                      starting_tree_scaling=options.starting_tree_scaling,
                                      starting_tree_use_scale_tree_factor=options.starting_tree_use_scale_tree_factor,
                                      scale_goal=options.scale_goal,
                                      mscale_file=mscale_file,
                                      no_add=no_add)

    with open(prefix+options.save_df_file, 'w') as f:
        f.write(str(df))


    make_topological_summaries = options.stop_criteria and (options.stop_criteria_threshold>=0)
    summary_verbose_scheme, summaries=get_summary_scheme(majority_tree=options.summary_majority_tree,
                                              light_newick_tree_summaries=make_topological_summaries,
                                              full_tree=True, #can not think of a moment where you don't want this.
                                              proposals=mp[0],
                                              acceptance_rate_information=options.summary_acceptance_rate,
                                              admixture_proportion_string=options.summary_admixture_proportion_string,
                                              no_chains=options.MCMC_chains,
                                              verbose_level=options.verbose_level,
                                              only_coldest_chain=not options.save_warm_chains)

    sim_lengths=[options.m]*options.n
    if options.stop_criteria:
        sc=stop_criteria(frequency=options.stop_criteria_frequency,
                         outfile=prefix+'stop_criteria.txt',
                         topological_threshold=options.stop_criteria_threshold,
                         continuous_threshold=options.stop_criteria_threshold, #ANDREWDEBUG This used to be topological.
                         Rscript_path=options.Rscript_command,
                         verbose_level=options.verbose_level)
    else:
        sc=None

    collapse_row=''
    likelihood_nodes=reduced_nodes
    posterior = posterior_class(emp_cov=covariance[0],
                                M=df,
                                p=options.p,
                                use_skewed_distr=options.sap_analysis,
                                multiplier=covariance[1],
                                nodes=likelihood_nodes,
                                use_uniform_prior=not options.not_uniform_prior,
                                treemix=False,
                                add_variance_correction_to_graph=(options.variance_correction != 'None' and
                                                                  options.add_variance_correction_to_graph),
                                prefix=prefix,
                                variance_correction_file=options.variance_correction_input_file,
                                prior_run=options.prior_run,
                                unadmixed_populations=options.unadmixed_populations,
                                collapse_row=collapse_row)

    posterior_function_list=[]


    temperature_scheme=fixed_geometrical(options.max_temp,options.MCMC_chains)

    #####    ANDREW DEBUG   !!!!!!!!!!
    if os.path.exists("DF.txt"):
        os.remove("DF.txt")
    else:
        print("The file does not exist")

    if os.path.exists("m_scale.txt"):
        os.remove("m_scale.txt")
    else:
        print("The file does not exist")

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
               store_permuts=options.store_permuts,
               stop_criteria=sc,
               make_outfile_stills=options.save_after_hours,
               save_only_coldest_chain=not options.save_warm_chains,
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

    if options.MCMC_chains==1:
        single_chain_run()
    else:
        res=multi_chain_run()
    if os.path.exists("stop_criteria.txt"):
        os.remove("stop_criteria.txt")
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
