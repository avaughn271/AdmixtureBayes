# AdmixtureBayes

## Purpose
AdmixtureBayes is a program to generate, analyze, and plot posterior samples of admixture graphs (phylogenies incorporating admixture events) given an allele count file.

## Installation

AdmixtureBayes can be downloaded install by running the commands
```bash
$ git clone https://github.com/avaughn271/AdmixtureBayes
```

It may also be necessary to install graphviz separately. TODO!!!!!!!!!!!!

### Test installation

A test script is found in the *example/* folder together with a test dataset. TODOOO!!!!!

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

 ## This step produces (in the current working directory)

**result_file** This file contains the list of samples that is retained after burn-in and thinning. A few lines describing the burn-in and thinning process are also printed to the console.

## (3) makePlots

In this step, we plot admixture graphs that have been sampled by AdmixtureBayes.  The script to run is

```bash
$ python PATH/AdmixtureBayes/admixturebayes/makePlots.py
```

There are 4 different plots that can be generated:

(i) the top trees

(ii) the top trees with branch estimates

(iii) the top node trees (TODO!!!1)

(iv) the consensus trees

## (i) Top Trees - input is of the form:

**--plot top_trees**

**--posterior** This is the output file from the posterior step.

**--top_trees_to_plot** (optional) The number of top trees to plot. If the value is greater than the total number of minimal topologies, the total number of minimal topologies will be used instead. Default value is 3.

**--write_rankings** (optional) The name of a file to which the minimal topologies and their likelihoods?????? will be written. If not specified, no such file will be written.


## This produces (in the current working directory)

**topology_i.pdf** One such plot will be produced for all *i* from 1 to the given number of trees to plot.

If **--write_rankings** is specified, a file with the set of all minimal topologies and their likelihoods!!! will be produced.

## (ii) Top Trees with Branch Estimates - input is of the form:

**--plot estimates**

**--posterior** This is the output file from the posterior step.

**--top_trees_to_estimate** (optional) The number of top trees with branch estimates to plot. If the value is greater than the total number of minimal topologies, the total number of minimal topologies will be used instead. Default value is 3.

## This produces (in the current working directory)

**topology_labels_i.pdf** One such plot will be produced for all *i* from 1 to the given number of trees to plot.

**branch_estimates_i.txt** A txt file describing the estimated branch lengths for the corresponding plot. One such file will be produced for all *i* from 1 to the given number of trees to plot.

**admixture_estimates_i.txt** A txt file describing the estimated admixture proportions for the corresponding plot. One such file will be produced for all *i* from 1 to the given number of trees to plot.

Important Note: In order for these plots to be produced properly, the analyzeSamples step must also be run with the flag "--faster False". This will result in an increased runtime for the analyzeSamples step, but will produce the necessry information for the branch estimates to be plotted.

## (iii) Top Node Trees - input is of the form:

**--plot top_node_trees**

**--posterior** This is the output file from the posterior step.

**--top_node_trees_to_plot** (optional) The number of top node trees to plot. If the value is greater than the total number of minimal topologies??????, the total number of minimal topologies will be used instead. Default value is 3.

**--write_rankings** (optional) The name of a file to which the minimal topologies and their likelihoods?????? will be written. If not specified, no such file will be written.

## This produces (in the current working directory)

**minimal_topology_i.pdf** One such plot will be produced for all *i* from 1 to the given number of trees to plot.

If **--write_rankings** is specified, a file with the set of all minimal topologies and their likelihoods!!! will be produced.


## (iv) Consensus Trees - input is of the form:

**--plot consensus_trees**

**--posterior** This is the output file from the posterior step.

**--consensus_thresholds** (optional) A list of values strictly between 0 and 1. For each value, a plot will be generated containing TODO!!!!. Default value is 0.25 0.5 0.75 0.9 0.95 0.99.

**--write_rankings** (optional)   is specified, a file with the set of nodes and their likelihoods is generated.

## This produces (in the current working directory)

**consensus_i.pdf** One such plot will be produced for all *i* in the list of thresholds.

If **--write_rankings** is specified, a file      containing a list of all nodes and their likelihoods is generated.

## Increasing number of populations

As more populations are added to the input file, the more steps it will take for the MCMC to converge. By default there are 50 MCMC steps between each MCMCMC flip and the number of MCMCMC flips is 200. The total number of MCMC steps is therefore 200\*50=10,000 which is only suitable for datasets with 4 or fewer populations. To increase the number of steps to 1,000,000 which is often enough to analyze 10 populations, use the command  

It is also possible to stop the chain using a stopping criteria. The stopping criteria calculates the Effective Sample Size of different summaries and stops if all of them are above a certain threshold (default is 200). To do so, AdmixtureBayes calls the [rwty package](https://cran.r-project.org/web/packages/rwty/index.html) in R with the command Rscript. To use the stopping criteria it is therefore necessary to install rwty
```bash
$ R
...
> install.packages("rwty")
```
and run AdmixtureBayes with the command

```bash
$ AdmixtureBayes run --input_file big_allele_counts.txt --outgroup population10 --n 50000 --stop_criteria
```

The parameter --n then defines an upper limit on the number of iterations.

## More advanced functionalities.
