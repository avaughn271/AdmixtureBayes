#Step 1:

#This is the basic workflow
python admixturebayes/runMCMC.py  --input_file example/allele_counts.txt --outgroup out

#We could also run the chain for more iterations with more chains
python admixturebayes/runMCMC.py  --input_file example/allele_counts.txt --outgroup out --n 500 --MCMC_chains 10

#Or use the stop criteria
python admixturebayes/runMCMC.py  --input_file example/allele_counts.txt --outgroup out --n 5000   --stop_criteria --stop_criteria_threshold 4

#We can also turn off printing to the console and rename the output
python admixturebayes/runMCMC.py  --input_file example/allele_counts.txt --outgroup out --n 5000   --stop_criteria --stop_criteria_threshold 400 --verbose_level silent --result_file exampleoutput.csv

#This is how we start in a random state
python admixturebayes/runMCMC.py  --input_file example/allele_counts.txt --outgroup out --n 5000 --result_file exampleoutput.csv --start random

#And this is how we would continue the run from our previous command
python admixturebayes/runMCMC.py  --input_file example/allele_counts.txt --outgroup out --n 5000  --continue_samples exampleoutput.csv --result_file exampleoutput_more.csv



#Step 2:

#This is the basic workflow
python admixturebayes/analyzeSamples.py   --mcmc_results  mcmc_samples.csv 

#We could also change the burn-in fraction and thinning rate
python admixturebayes/analyzeSamples.py --mcmc_results  mcmc_samples.csv --burn_in_fraction 0.2  --thinning_rate 5

#Or rename the output file
python admixturebayes/analyzeSamples.py   --mcmc_results  mcmc_samples.csv --result_file example_output.csv



#Step 3:

#We can plot the top trees 
python admixturebayes/makePlots.py  --plot top_trees  --posterior thinned_samples.csv

#We can also change the number of top trees to plot and save all of the topologies along with their posterior probabilities
python admixturebayes/makePlots.py  --plot top_trees  --posterior thinned_samples.csv --top_trees_to_plot  5  --write_rankings example_rankings.txt

#We can also plot the top trees with branch estimates, but this requires running the previous step with the --slower flag 
python admixturebayes/analyzeSamples.py --mcmc_results  mcmc_samples.csv --burn_in_fraction 0.2  --thinning_rate 5 --slower
python admixturebayes/makePlots.py  --plot estimates  --posterior thinned_samples.csv --top_trees_to_estimate 2

#Or plot the minimal topologies
python admixturebayes/makePlots.py  --plot top_minimal_topologies  --posterior thinned_samples.csv --top_minimal_topologies_to_plot 4 --write_rankings example_rankings.txt

#Lastly, we can plot the consensus trees 

python admixturebayes/makePlots.py  --plot consensus_trees  --posterior thinned_samples.csv

#And change the consensus thresholds and write the node rankings as well
python admixturebayes/makePlots.py  --plot consensus_trees  --posterior thinned_samples.csv --consensus_thresholds  0.3 0.4 0.6 --write_rankings example_rankings.txt
