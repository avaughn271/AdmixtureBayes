# AdmixtureBayes

## Purpose
AdmixtureBayes is a program to generate, analyze, and plot posterior samples of admixture graphs (phylogenies incorporating admixture events) given an allele count file.

## Installation

AdmixtureBayes can be downloaded by running the commands
```bash
$ git clone https://github.com/avaughn271/AdmixtureBayes
```
AdmixtureBayes is written almost completely in Python and requires the following Python packages:

"numpy", "scipy",  "pandas",    "pathos",   "graphviz"

See the following links for installation help:

https://numpy.org/install/

https://scipy.org/install/

https://pandas.pydata.org/docs/getting_started/install.html

https://pypi.org/project/pathos/

https://pypi.org/project/graphviz/

Furthermore, if you wish to use the stop_criteria feature of the Markov Chain, then you need to also have R installed with the [rwty](https://cran.r-project.org/web/packages/rwty/index.html) and [coda](https://cran.r-project.org/web/packages/coda/index.html) packages, which can be installed by running: 
```
install.packages(c("rwty", "coda"))
```
in any R session. 

### Example commands

A script containing example commands is found in the *example/* folder together with a test dataset.

## Input file

The input for AdmixtureBayes is an allele count file in the exact same format as used by TreeMix.

```bash
s1 s2 s3 s4 out
9,11 13,7 11,9 14,6 14,6
4,16 4,16 0,20 1,19 2,18
...
1,19 2,18 3,17 2,18 2,18
```

where the first line is the populations and the subsequent lines are the bi-allelic counts in each population for a number of SNPs. The first and second allele type has no meaning and can be chosen arbitrarily.

## Running AdmixtureBayes

AdmixtureBayes has 3 steps:

(1) *runMCMC* - this takes the input of allele counts described above, runs the MCMC chain, and generates a set of samples of admixture graphs

(2) *analyzeSamples* - this takes the output of the previous step and performs a burn-in and thinning step to generate independent samples of admixture graphs

(3) *makePlots* - this takes the output of the previous step and generates different plots that are useful for interpretation

## (1) runMCMC

In this step, we run the MCMC chain that explores the space of admixture graphs.  The script to run is

```bash
$ python PATH/AdmixtureBayes/admixturebayes/runMCMC.py
```

## This step takes as input:

**--input-file** The input file of allele counts as described above.

**--outgroup** The name of the population that will serve as the outgroup. For example, in the above file, "out" could be the outgroup.

**--n** (optional) The number of proposal steps the MCMC sampler should make. (Technically, this is the number of MCMCMC flips the chain should make, which is directly proporional to the number of proposal steps). Default value is 200.

**--MCMC_chains** (optional) The number of chains to run the MCMCMC with. More chains will results in better mixing at the cost of increased computational time. AdmixtureBayes supports multiprocessing, so ideally this would be the number of cores. Default value is 8.


**--result_file** (optional) The name of the mcmc output file of this step. No file extension is added (meaning entering "example" will produce "example" as an output file, not "example.txt" or "example.csv".). Default value is "mcmc_samples.csv"

**--stop_criteria** (optional) If this flag is used, then the MCMC sampler will stop as soon as the effective sample size has been reached or until n iterations have been reached, whichever comes first. Otherwise, the algorithm will continue until n iterations have been reached.

**--stop_criteria_threshold** (optional) Ignored if stop_criteria is False. Sets the effective samples size that much be reached for the algorithm to terminate. Default value is 200.

**--Rscript_command** (optional) Ignored if stop_criteria is False. The command with which to run an R script on the desired machine (eg. "Rscript" or "R CMD BATCH"). Default value is "Rscript".

**--verbose_level** (optional) Either "normal" or "silent". If "normal", then the total number of snps will be printed to the console, along with the progress of the MCMC sampler. Every 1000th iteration, the progress towards the total number of iterations is printed. If "silent", then nothing will be printed to the console. Default value is "normal."


 ## This step produces (in the current working directory)

***result_file*** This file contains the list of MCMC samples. 


**covariance_and_multiplier.txt** This file contains the covariance matrix corresponding to the given input file.



## (2) analyzeSamples

In this step, we analyze the output of the Markov chain from the previous step.  The script to run is

```bash
$ python PATH/AdmixtureBayes/admixturebayes/analyzeSamples.py
```

## This step takes as input:

**--mcmc_results** This is the output file from the previous step that contains the mcmc samples.

**--covariance** This is the output file from the previous step that contains the covariance matrix.

**--result_file** (optional) The name of the output file of this step. No file extension is added (meaning entering "example" will produce "example" as an output file, not "example.txt" or "example.csv".). Default value is "thinned_samples.csv"

**--burn_in_fraction** (optional) The fraction of samples to discard as a burn-in period. Default value is 0.5.

**--thinning_rate** (optional) The thinning rate of the sample thinning step (which occurs after the burn-in step). Deafult value is 10.

**--subnodes** (optional) A list of population labels, which should be a subset of all non-outgroup population labels. This can be specified if only a subset of the population labels should be analyzed. Plots generated using the output of this step will only used the populations given in this step. If not specified, all non-outgroup populations will be considered.

**--slower** (optional) If this flag is used, then the necessary information to plot the top trees with branch estimates  will be computed. By default, this is not done as this can be a very slow process when the number of admixture events is large. This flag also has another effect. For each graph that contains admixture events, each admixture event has one ancestral lineage that is specifically marked as the "introgressed" lineage. Graphs are considered distinct if they have different marked lineages, even if their topologies are the same. However, if this flag is specified, then graphs with the same topology are combined into equivalence classes induced by the set of all introgression markings. This distinction might be important for considering the posterior probabilities of different topologies.

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

**--posterior** This is the output file from the posterior step.

**--top_trees_to_plot** (optional) The number of top trees to plot. If the value is greater than the total number of minimal topologies, the total number of minimal topologies will be used instead. Default value is 3.

**--write_rankings** (optional) The name of a file to which the minimal topologies and their posterior probabilities will be written. If not specified, no such file will be written.

## This produces (in the current working directory)

***topology_i.pdf*** One such plot will be produced for all *i* from 1 to the given number of trees to plot.

If **--write_rankings** is specified, a file with the set of all topologies and their posterior probabilities will be produced.

## (ii) Top Trees with Branch Estimates - input is of the form:

**--plot estimates**

**--posterior** This is the output file from the posterior step.

**--top_trees_to_estimate** (optional) The number of top trees with branch estimates to plot. If the value is greater than the total number of minimal topologies, the total number of minimal topologies will be used instead. Default value is 3.

## This produces (in the current working directory)

***topology_labels_i.pdf*** One such plot will be produced for all *i* from 1 to the given number of trees to plot.

***branch_estimates_i.txt*** A txt file describing the estimated branch lengths for the corresponding plot. One such file will be produced for all *i* from 1 to the given number of trees to plot.

***admixture_estimates_i.txt*** A txt file describing the estimated admixture proportions for the corresponding plot. One such file will be produced for all *i* from 1 to the given number of trees to plot.

Important Note: In order for these plots to be produced properly, the analyzeSamples step must also be run with the flag "--slower". This will result in an increased runtime for the analyzeSamples step but will produce the necessry information for the branch estimates to be plotted. 

## (iii) Top Minimal Topologies - input is of the form:

**--plot top_minimal_topologies**

**--posterior** This is the output file from the posterior step.

**--top_minimal_topologies_to_plot** (optional) The number of top node trees to plot. If the value is greater than the total number of minimal topologies, the total number of minimal topologies will be used instead. Default value is 3.

**--write_rankings** (optional) The name of a file to which the minimal topologies and their posterior probabilities will be written. If not specified, no such file will be written.

## This produces (in the current working directory)

***minimal_topology_i.pdf*** One such plot will be produced for all *i* from 1 to the given number of trees to plot.

If **--write_rankings** is specified, a file with the set of all minimal topologies and their posterior probabilities will be produced.

## (iv) Consensus Trees - input is of the form:

**--plot consensus_trees**

**--posterior** This is the output file from the posterior step.

**--consensus_thresholds** (optional) A list of values strictly between 0 and 1. For each value, a plot will be generated containing a graph formed by all nodes that have a posterior probability above that threshold. Default value is 0.25 0.5 0.75 0.9 0.95 0.99.

**--write_rankings** (optional)   is specified, a file with the set of nodes and their likelihoods is generated.

## This produces (in the current working directory)

***consensus_i.pdf*** One such plot will be produced for all *i* in the list of thresholds.

If **--write_rankings** is specified, a file  containing a list of all nodes and their posterior probabilities is generated.
