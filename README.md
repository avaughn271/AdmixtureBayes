# AdmixtureBayes

## Purpose
AdmixtureBayes is a program to generate, analyze, and plot posterior samples of admixture graphs (phylogenies incorporating admixture events) given an allele count file. AdmixtureBayes is currently maintained by [Andrew Vaughn](https://nielsen-lab.github.io/team/andrew-vaughn). Please report any strange results, errors, or code suggestions to him at [ahv36@berkeley.edu](mailto:ahv36@berkeley.edu). Please see [https://doi.org/10.1371/journal.pgen.1010410](https://doi.org/10.1371/journal.pgen.1010410) for our paper describing AdmixtureBayes.

## Installation

AdmixtureBayes can be downloaded from this GitHub repo by running the following command:
```bash
$ git clone https://github.com/avaughn271/AdmixtureBayes
```
AdmixtureBayes is written in Python and requires the following Python packages:

"numpy", "scipy",  "pandas",    "pathos",   "graphviz"

See the following links for installation help:

https://numpy.org/install/

https://scipy.org/install/

https://pandas.pydata.org/docs/getting_started/install.html

https://pypi.org/project/pathos/

https://pypi.org/project/graphviz/

Furthermore, if you wish to use the given R script to evaluate convergence, then you need to also have R installed with the [coda](https://cran.r-project.org/web/packages/coda/index.html) package, which can be installed by running: 
```
install.packages("coda")
```
in any R session. 

### Example commands

A script containing example commands is found in the *example* folder together with a test dataset.

## Input file

The input for AdmixtureBayes is an allele count file in the exact same format as used by [TreeMix](https://bitbucket.org/nygcresearch/treemix/wiki/Home).

```bash
s1 s2 s3 s4 out
9,11 13,7 11,9 14,6 14,6
4,16 4,16 0,20 1,19 2,18
...
1,19 2,18 3,17 2,18 2,18
```

where the first line is the populations and the subsequent lines are the bi-allelic counts in each population for a number of SNPs. The first and second allele type has no meaning and can be chosen arbitrarily. The population names should only include letters and numbers (no spaces, dashes, underscores, etc.). See the R script "ConvertFromVCF.R" in the "example" folder for a template for converting from VCF files to this input. At minimum, you will need to change the name of the input VCF file and the individual-to-population mapping in this script. Keep in mind that VCF files can be quite complex, and therefore this script may not work for all possible input VCF files. The user should always perform a sanity check between the input and output of this step and should not take the output at face value.

Notes on Missing Data: Missing data for a population at a particular site should be encoded as "0,0". AdmixtureBayes estimates allelic covariance matrices, one large matrix consisting of all SNPs and one for each bootstrap sample of adjacent SNPs to be used in the estimation of the Wishart degrees of freedom. Entry *(i,j)* in these matrices is the covariance in allele frequencies between populations *i* and *j* using all relevant SNPs (either all SNPs or the SNPs in the corresponding bootstrap block). If there is no missing data in either of these populations, then this is the dot product of the allele frequency vectors at the relevant SNPs for populations *i* and *j* divided by the total number of relevant SNPs. If a site has missing data for either *i* or *j*, then it is not included in the dot product, and we instead divide by the total number of relevant SNPs that are missing in neither *i* nor *j*. Missing data is not a problem for AdmixtureBayes to handle, but it does violate the assumption of even sampling imposed by the Wishart distribution. Relevant warnings for missing data will be printed to the console, although the algorithm will still run the MCMC properly.


## Running AdmixtureBayes

AdmixtureBayes has 3 steps:

(1) *runMCMC* - this takes the input of allele counts described above, runs the MCMC chain, and generates a set of samples of admixture graphs

(2) *analyzeSamples* - this takes the output of the previous step and performs a burn-in and thinning step to generate independent samples of admixture graphs

(3) *makePlots* - this takes the output of the previous step and generates different plots that are useful for interpretation

Notes on Runtime: Steps 2 and 3 should be very fast, regardless of the input. The runtime of step 1 is determined by how many iterations are run and the number of populations being considered. Increasing the number of populations will result in more time per iteration. Increasing the number of populations also increases the size of the state space of admixture graphs, resulting in more iterations being necessary to achieve convergence. Runtime is invariant with respect to the number of individuals in each population. Increasing the number of SNPs will keep the time per iteration the same, but will result in more time for the initial step of calculating the allele covariance matrix. Asymptotically, as the number of iterations increases, this will be a negligible fraction of the total runtime. Keep in mind that runtime necessary to achieve convergece increases exponentially with the number of populations due to the dramatic increase in the state space and the number of possible proposals per state. While AdmixtureBayes should be able to easily converge on 4 or 5 non-outgroup populations within 1 hour on a desktop, one can expect a thorough analysis on 10 populations to take dozens of hours. Usage of a computing cluster is recommended for datasets with many populations.

## (1) runMCMC

In this step, we run the MCMC chain that explores the space of admixture graphs.  The script to run is

```bash
$ python PATH/AdmixtureBayes/admixturebayes/runMCMC.py
```

## This step takes as input:

**--input_file** The input file of allele counts as described above.

**--outgroup** The name of the population that will serve as the outgroup. For example, in the above file, "out" could be the outgroup.

**--n** (optional) The number of iterations the MCMC sampler should make. (Technically, this is the number of MC<sup>3</sup> flips the chain should make, which is directly proportional to the number iterations. The exact number of iterations is 50*n). Default value is 200. This number should be increased in all practical applications.

**--result_file** (optional) The name of the mcmc output file of this step. No file extension is added (meaning entering "example" will produce "example" as an output file, not "example.txt" or "example.csv".). Default value is "mcmc_samples.csv"

**--continue_samples** (optional) This is used if you want to continue a previous AdmixtureBayes run, for example if convergence was not yet reached. The argument passed to --continue_samples should be the file name of the MCMC sample file produced by a previous AdmixtureBayes run (which is "mcmc_samples.csv" by default). A new file will be produced (whose name will be whatever the input to --result_file is in this call to AdmixtureBayes) that will contain all of the samples of the previous run in addition to all of the samples from this run. This call to AdmixtureBayes will start the chain in the last state of the previous AdmixtureBayes run. The previous output file will not be overwritten.  The name of the --result_file argument used in this call to AdmixtureBayes should be different than the one produced by the previous call to avoid unwanted behavior. For example, if a user does not specify  --result_file for either call, then both will be "mcmc_samples.csv" by default and unwanted behavior will occur. If this argument is not specified, then the algorithm will start at a randomly constructed graph, which may be useful for monitoring mixing and convergence of the chain, for example by Gelman-Rubin statistics.

**--verbose_level** (optional) Either "normal" or "silent". If "normal", then the total number of snps will be printed to the console along with the progress of the MCMC sampler. Every 1000th iteration, the progress towards the total number of iterations is printed. If "silent", then nothing will be printed to the console. Default value is "normal."

**--save_covariance** (optional) If this flag is specified, then the allelic covariance matrix produced by considering all SNPs will be saved to the file "covariance_matrix.txt" in the current working directory. Note that this will be the covariance matrix described in the AdmixtureBayes paper, which is to say a scaled, bias-corrected transformation of the naive covariance matrix that would be suggested by Equation 4 of the main text. The user may use the input data to compute the naive covariance matrix using Equation 4 of the main text should they choose, but bear in mind that this is different than the covariance matrix AdmixtureBayes is actually using.

**--MCMC_chains** (optional) A tuning parameter of the MCMC algorithm. The number of chains to run the MC<sup>3</sup> with (See Matthew Darlington's great explanation of MC<sup>3</sup> [here](https://www.lancaster.ac.uk/stor-i-student-sites/matthew-darlington/wp-content/uploads/sites/10/2020/01/MattDReport.pdf)). More chains should result in better mixing at the cost of increased computational time. AdmixtureBayes supports multiprocessing, so ideally this would be the number of cores. Default value is 8. Must be at least 2. See the section below on MC<sup>3</sup> Mixing for more details.

**--maxtemp** (optional)  A tuning parameter of the MCMC algorithm. The temperature of the hottest chain in the MC<sup>3</sup> algorithm. Must be a positive number, though not necessarily an integer. The temperature of the $i$'th chain will be $maxtemp^{[(i-1)/(MCMC\textunderscore chains-1)]^{spacing}}$. We propose swaps between chains of adjacent temperature, meaning that we can measure the acceptance rate of $M-1$ different proposals. At the end of the runMCMC step, the acceptance rates of proposed swaps between chains is printed as a list of length $M-1$. The $i$'th element of this list is the acceptance rate of proposed swaps between chain $i$ and chain $i+1$. The acceptance rates are controlled by this argument and the *--spacing* argument. Ideally all of the acceptance rates will be near 20-40%, but this can be difficult to achieve for each pair of chains. One idea which practically works well is try and set *--maxtemp* and *--spacing* such that  all acceptance rates are above 10%. Default value is 1000. See the section below on MC<sup>3</sup> Mixing for more details.

**--spacing** (optional)  A tuning parameter of the MCMC algorithm. Determines the relative spacing of the temperatures of the MC<sup>3</sup> algorithm. Must be a positive number, though not necessarily an integer. The temperature of the $i$'th chain will be $maxtemp^{[(i-1)/(MCMC\textunderscore chains-1)]^{spacing}}$. Default value is 1. See the section below on MC<sup>3</sup> Mixing for more details.

**--num_ind_snps** (optional) This sets the number of effectively independent SNPs (which we call *K*) that is used in the likelihood function. By default, this is estimated through bootstrapping and the resulting estimated value of *K* is printed at the beginning of the algorithm. However, as bootstrapping is stochastic,  *K*  can differ slightly between runs of the algorithm. This means that technically the likelihood function will slightly different when running AdmixtureBayes multiple times on the same dataset.  In certain circumstances, this can lead to Gelman-Rubin statistics that suggest poor mixing between chains, due simply to slight differences in the likelihood function. The idea behind being able to set this argument explicitly is to observe the estimated value of *K* for one chain, and then use this set value of *K* for all   subsequent AdmixtureBayes analyses on the same dataset to remove this phenomenon of different likelihood functions as a source of error.

**--max_admixes** (optional) This sets the maximum number of admixture events that graphs are allowed to have. The default value is 20. It may be useful to set this value lower to achieve better mixing on the conditional distribution of graphs given a fixed maximum of admixture events (see the *MC<sup>3</sup> Mixing* section below for further discussion). 

 ## This step produces (in the current working directory)

***result_file*** This file contains the list of MCMC samples. Each line of the output file describes a sampled admixture graph. Many of the column values concern the internal representation of the graph, such as "tree" and "descendant_sets", and are not meant for user interpretation. However, of interest to the user might be the columns "prior" (the log-prior), "posterior" (the log-posterior), "likelihood" (the log-likelihood), "total_branch_length", "ghost_pops" (number of ghost populations), and "no_admixes" (number of admixture events). These are summary statistics of the MCMC chain and may be useful for monitoring convergence. Some thinning is done on the samples to reduce the size of the output file. The total number of samples saved will be approximately 1.25*n. 


## (2) analyzeSamples

It is here that convergence of the MCMC sampler should be assessed, for example by examining the trace plots of the chain or Gelman-Rubin convergence diagnostics of parallel chains. We have included the sample R script EvaluateConvergence.R as a template for assessing convergence. If convergence has not been reached, you should run the chain for more iterations (by increasing **--n**, possibly in conjunction with the **--continue_samples** argument) and/or increase the number of parallel MCMC chains being run (by increasing the **--MCMC_chains** argument).

In this step, we analyze the output of the Markov chain from the previous step.  The script to run is

```bash
$ python PATH/AdmixtureBayes/admixturebayes/analyzeSamples.py
```

## This step takes as input:

**--mcmc_results** This is the output file from the previous step that contains the mcmc samples.

**--result_file** (optional) The name of the output file of this step. No file extension is added (meaning entering "example" will produce "example" as an output file, not "example.txt" or "example.csv".). Default value is "thinned_samples.csv"

**--burn_in_fraction** (optional) The fraction of samples to discard as a burn-in period. Default value is 0.5.

**--thinning_rate** (optional) The thinning rate of the sample thinning step (which occurs after the burn-in step). Default value is 10.

**--subnodes** (optional) A list of population labels, which should be a subset of all non-outgroup population labels. This can be specified if only a subset of the population labels should be analyzed. Plots generated using the output of this step will only use the populations given in this step. If not specified, all non-outgroup populations will be considered.

 ## This step produces (in the current working directory)

***result_file*** This file contains the list of samples that is retained after burn-in and thinning. A few lines describing the burn-in and thinning process are also printed to the console.

## (3) makePlots

In this step, we plot admixture graphs that have been sampled by AdmixtureBayes.  The script to run is

```bash
$ python PATH/AdmixtureBayes/admixturebayes/makePlots.py
```

There are 4 different plots that can be generated:

(i) the top trees - these are the topologies with the highest posterior probabilities.

(ii) the top trees with branch estimates - these are the topologies with the highest posterior probabilities along with the best estimates of branch lengths and admixture proportions.

(iii) the top node trees - these are the minimal topologies with the highest probabilities.

(iv) the consensus trees - these are trees that are formed by combining all nodes that have a given posterior probability of appearing in an admixture graph.

## (i) Top Trees - input is of the form:

**--plot top_trees**

**--posterior** This is the output file from the *analyzeSamples* step.

**--top_trees_to_plot** (optional) The number of top trees to plot. If the value is greater than the total number of topologies, the total number of topologies will be used instead. Default value is 3.

**--write_rankings** (optional) The name of a file to which the topologies and their posterior probabilities will be written. If not specified, no such file will be written.

**--output_prefix** (optional) The prefix of the output plots.


## This produces (in the current working directory)

***output_prefix*** **_i.pdf** One such plot will be produced for all *i* from 1 to the given number of trees to plot, where ***output_prefix*** is the argument to the **--output_prefix**  option. If no such argument is specified, the plots will have the name  ***topology_i.pdf***.

If **--write_rankings** is specified, a file with the set of all topologies and their posterior probabilities will be produced.

## (ii) Top Trees with Branch Estimates - input is of the form:

**--plot estimates**

**--posterior** This is the output file from the *analyzeSamples* step.

**--top_trees_to_estimate** (optional) The number of top trees with branch estimates to plot. If the value is greater than the total number of minimal topologies, the total number of minimal topologies will be used instead. Default value is 3.

**--output_prefix** (optional) The prefix of the output files.

## This produces (in the current working directory)

***output_prefix*** **_i.pdf** One such plot will be produced for all *i* from 1 to the given number of trees to plot.  If no such argument is specified, the plots will have the name  ***topology_labels_i.pdf***. 

***output_prefix*** **_branch_estimates_i.txt** A txt file describing the estimated branch lengths for the corresponding plot. One such file will be produced for all *i* from 1 to the given number of trees to plot. If no such argument is specified, the files will have the name  ***branch_estimates_i.pdf***.

***output_prefix*** **_admixture_estimates_i.txt** A txt file describing the estimated admixture proportions for the corresponding plot. One such file will be produced for all *i* from 1 to the given number of trees to plot.  If no such argument is specified, the files will have the name  ***admixture_estimates_i.pdf***.

Important Note: You cannot plot branch length estimates if the option ``--subnodes" was used in the previous step.

## (iii) Top Minimal Topologies - input is of the form:

**--plot top_minimal_topologies**

**--posterior** This is the output file from the *analyzeSamples* step.

**--top_minimal_topologies_to_plot** (optional) The number of minimal topologies to plot. If the value is greater than the total number of minimal topologies, the total number of minimal topologies will be used instead. Default value is 3.

**--write_rankings** (optional) The name of a file to which the minimal topologies and their posterior probabilities will be written. If not specified, no such file will be written.

**--output_prefix** (optional) The prefix of the output plots.

## This produces (in the current working directory)

***output_prefix*** **_i.pdf** One such plot will be produced for all *i* from 1 to the given number of trees to plot.   If no such argument is specified, the plots will have the name  ***minimal_topology_i.pdf***.

If **--write_rankings** is specified, a file with the set of all minimal topologies and their posterior probabilities will be produced.

## (iv) Consensus Trees - input is of the form:

**--plot consensus_trees**

**--posterior** This is the output file from the *analyzeSamples* step.

**--consensus_thresholds** (optional) A list of values strictly between 0 and 1. For each value, a plot will be generated containing a graph formed by all nodes that have a posterior probability above that threshold. Default value is 0.25 0.5 0.75 0.9 0.95 0.99.

**--write_rankings** (optional)   is specified, a file with the set of nodes and their posterior probabilities is generated.

**--output_prefix** (optional) The prefix of the output plots. 


## This produces (in the current working directory)

***output_prefix*** **_i.pdf** One such plot will be produced for all *i* in the list of thresholds. If no such argument is specified, the plots will have the name  ***consensus_i.pdf***.

If **--write_rankings** is specified, a file  containing a list of all nodes and their posterior probabilities is generated.

## Getting Started and Troubleshooting

### 1) Example dataset

The first thing you should do is to try and run AdmixtureBayes on the provided example dataset to get a feel for how to assess convergence and interpret the output. 

You can run 3 independent chains by running the following 3 commands.

```bash
python PATH/admixturebayes/runMCMC.py --input_file PATH/example/allele_counts.txt --outgroup out --n 10000 --result_file chain1.txt
python PATH/admixturebayes/runMCMC.py --input_file PATH/example/allele_counts.txt --outgroup out --n 10000 --result_file chain2.txt
python PATH/admixturebayes/runMCMC.py --input_file PATH/example/allele_counts.txt --outgroup out --n 10000 --result_file chain3.txt
```
This takes approximately 5-15 minutes per chain on a desktop computer, depending on processor speed and the number of cores present. Then, we can assess convergence by running 

```bash
Rscript EvaluateConvergence.R
```

Note that the names of the output files to analyze are hardcoded as chain1.txt, chain2.txt, and chain3.txt. You will need to change these in the script if you change the output file names. The output of this script should look like the convergence plots presented in the AdmixtureBayes paper (S11 Fig, S12 Fig, and S12 Fig).  Trace plots should stabilize to an equilibrium distribution, Gelman-Rubin statistics should be near 1, and autocorrelation plots should converge to 0.

If the chains have all indeed converged, then we should filter our output graphs with the command 

```bash
python PATH/admixturebayes/analyzeSamples.py --mcmc_results chain1.txt
```

and plot the top sampled topologies and their rankings with the following command

```bash
python PATH/admixturebayes/makePlots.py --plot top_trees --posterior thinned_samples.csv --write_rankings chain1rankings.txt
```


### 2) Running AdmixtureBayes on your dataset

Now that you are familiar with the output of AdmixtureBayes and with assessing convergence, you can run AdmixtureBayes on your dataset. If your number of populations is larger than 5, we recommend initially running AdmixtureBayes on a subset of your data, 3 or 4 non-outgroup populations for example. This is because AdmixtureBayes should converge very rapidly on this smaller dataset, and it will enable you to work out any data conversion problems or mixing problems rather quickly. Once you are confident that AdmixtureBayes is performing properly on your data, you can move on to your full dataset. AdmixtureBayes, when using 32 parallel chains through the --MCMC_chains option, takes about 50 hours to do a complete run on our Arctic dataset of 11 non-outgroup populations (we used --n 450000). The state space of admixture graph topologies grows super-exponentially in the number of populations, so be aware that when considering 15 or 20 populations, the time necessary to achieve convergence may increase considerably. 

### 3) MC<sup>3</sup> Mixing

The most common problem users may experience when using AdmixtureBayes is a lack of convergence of the MC<sup>3</sup> algorithm, as measured through the EvaluateConvergence.R script. This is more likely to be the case when you are analyzing a large number of populations and/or your data has a large number of effectively independent SNPs. The best ways to try and fix convergence problems are: 

1) Setting --n to be higher. In other words, run the chain longer. Such a change is very simple to make, but is nevertheless a little unsatisfying as one can always simply try running an MCMC algorithm longer in the hopes of better convergence.

2) Improve the swap rates between the heated chains of the MC<sup>3</sup> algorithm. The particular benefit of  MC<sup>3</sup> is that it helps the sampler traverse regions of the state space with low probability in order to help sample different local maxima. However, this only works properly if chains of adjacent temperatures swap at an appreciable rate. This requires the exact temperatures of the heated chains to be set precisely, and the temperatures which give optimal mixing will depend on the specific dataset being analyzed.

As mentioned previously, when running the runMCMC.py script, the acceptance rates of proposed swaps between chains are printed at the end of the algorithm as a list of length $M-1$ (where $M$ is the number of heated chains). The $i$'th element of this list is the acceptance rate of proposed swaps between chain $i$ and chain $i+1$. This should help the user assess mixing between chains. The swap rates between chains are affected by the values of the *--maxtemp* and *--spacing*, which determine the temperatures of the  MC<sup>3</sup> algorithm. A small value of 
d *--spacing*, such as 0.1 will result in a skew of temperatures towards  *--maxtemp*. A large value of  *--spacing*, such as 100, will result in a skew of temperatures towards 1.0. Therefore, if the mixing rate of the first few chains (which are the coldest) is high, but the mixing rate of the last few chains (which are the hottest) is too low,  *--spacing*  should be set lower. If  the mixing rate of the first few chains (which are the coldest) is too low, but the mixing rate of the last few chains (which are the hottest) is  high,  *--spacing*  should be set higher.

3) Set a maximum number of allowable admixture events. The maximum number of admixture events is 20 by default. By decreasing the maximum number of allowable events, we reduce the space of possible admixture graph topologies, which greatly improves mixing. This can be done with the  --max_admixes argument (described above). For example, you can run separate analyses with  --max_admixes 0,   --max_admixes 1,  --max_admixes 2, etc. if you are intereted in the conditional distribution of graphs given a fixed number of admixture events.



