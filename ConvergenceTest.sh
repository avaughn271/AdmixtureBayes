 python3 admixturebayes/runMCMC.py   --input_file example/allele_counts.txt --outgroup out 

 python3 ~/desktop/AdmixtureBayes/admixturebayes/analyzeSamples.py --mcmc_results  mcmc_samples.csv --burn_in_fraction 0.98  --thinning_rate 5

 python3 ~/desktop/AdmixtureBayes/admixturebayes/makePlots.py --plot top_trees --posterior thinned_samples.csv --top_trees_to_plot  1