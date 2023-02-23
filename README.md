# AdmixtureBayes

## Purpose
We here describe the extension of AdmixtureBayes to simulated annealing in order to generate a single, best point estimate of an admixture graph. 

## Installation

## Running AdmixtureBayes for Annealing

AdmixtureBayes has 3 steps:

(1) *runSA* - this takes the input of allele counts described above, runs the MCMC chain, and generates a set of samples of admixture graphs

## (1) runSA

In this step, we run the simulated annealing that explores the space of admixture graphs and returns a best point estimate.  The script to run is

```bash
$ python PATH/AdmixtureBayes/admixturebayes/runSA.py
```

## This step takes as input:

**--input_file** The input file of allele counts as described above.

**--outgroup** The name of the population that will serve as the outgroup. For example, in the above file, "out" could be the outgroup.

**--output_prefix** The prefix of the output files. Default value is 'out'.

**--starting_temp**  The initial temperature of the simulated annealing algorithm.

**--ending_temp**  The final temperature after which the simulated annealing algorithm terminates.

**--temp_scaling**  The scaling factor that changes the temperature. Must be between 0 and 1.

**--iter_per_temp**  The number of proposals per temperature.


 ## This step produces (in the current working directory)

***output_prefix*** **_tree** A table representing the tree. Each row is an edge. Column 1 is the child node of this edge. Column 2 is the parent node of this edge. Column3 is the branch length. Column 4 is the admixture proportion going through this edge if the child node is an admixture node (it is 1.0 otherwise). 

***output_outgroup*** **_tree** The distance from the root to the outgroup.

***output_prefix*** **_topology** This is a plot of the topology of the MAP graph.

***output_prefix*** **_minimal_topology** This is a plot of the minimal topology of the MAP graph.

***output_prefix*** **_topology_labels** This is a plot of the MAP graph, including topology, branch lengths, and admixture proportions. The distance to the outgroup is also shown.