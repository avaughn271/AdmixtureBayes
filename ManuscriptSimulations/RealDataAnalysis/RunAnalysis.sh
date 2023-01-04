#We here run Step 1, the MCMC. Each chain was extended as convergence was not reached after the first run of AdmixtureBayes

python admixturebayes/runMCMC.py --input_file ArcticData.txt --outgroup Yoruba --n 450000 --MCMC_chains 32
python admixturebayes/runMCMC.py --input_file ArcticData.txt --outgroup Yoruba --n 450000 --MCMC_chains 32 --continue_samples mcmc_samples.csv --result_file chain1.csv

python admixturebayes/runMCMC.py --input_file ArcticData.txt --outgroup Yoruba --n 450000 --MCMC_chains 32
python admixturebayes/runMCMC.py --input_file ArcticData.txt --outgroup Yoruba --n 450000 --MCMC_chains 32 --continue_samples mcmc_samples.csv --result_file chain2.csv

python admixturebayes/runMCMC.py --input_file ArcticData.txt --outgroup Yoruba --n 450000 --MCMC_chains 32
python admixturebayes/runMCMC.py --input_file ArcticData.txt --outgroup Yoruba --n 450000 --MCMC_chains 32 --continue_samples mcmc_samples.csv --result_file chain3.csv

#It is here that convergenve is assessed. As all 3 chains converged to the same stationary distribution, we simply picked chain 1 to analyze
python admixturebayes/analyzeSamples.py  --mcmc_results  chain1.csv  --burn_in_fraction 0.35  --thinning_rate 40 --slow

#We plot the top minimal topologies
python admixturebayes/makePlots.py --plot top_minimal_topologies --posterior thinned_samples.csv --top_minimal_topologies_to_plot 2 --write_rankings chain1rankings.txt

#We plot the consensus graph

python admixturebayes/makePlots.py --plot consensus_trees --posterior thinned_samples.csv  --consensus_thresholds 0.75

#We also plot the relevant subgraphs we wish to analyze by editing our analyze samples step

python admixturebayes/analyzeSamples.py  --mcmc_results  chain1.csv  --burn_in_fraction 0.35  --thinning_rate 40 --slow --subnodes Athabcascan Koryak Saqqaq
python admixturebayes/makePlots.py --plot top_minimal_topologies --posterior thinned_samples.csv --top_minimal_topologies_to_plot 3 --write_rankings chain1_1.txt

python admixturebayes/analyzeSamples.py  --mcmc_results  chain1.csv  --burn_in_fraction 0.35  --thinning_rate 40 --slow --subnodes Greenlander Koryak Saqqaq
python admixturebayes/makePlots.py --plot top_minimal_topologies --posterior thinned_samples.csv --top_minimal_topologies_to_plot 3 --write_rankings chain1_2.txt

python admixturebayes/analyzeSamples.py  --mcmc_results  chain1.csv  --burn_in_fraction 0.35  --thinning_rate 40 --slow --subnodes Greenlander Koryak Athabascan
python admixturebayes/makePlots.py --plot top_minimal_topologies --posterior thinned_samples.csv --top_minimal_topologies_to_plot 3 --write_rankings chain1_3.txt

python admixturebayes/analyzeSamples.py  --mcmc_results  chain1.csv  --burn_in_fraction 0.35  --thinning_rate 40 --slow --subnodes Athabcascan Koryak Saqqaq Greenlander
python admixturebayes/makePlots.py --plot top_minimal_topologies --posterior thinned_samples.csv --top_minimal_topologies_to_plot 3 --write_rankings chain1_4.txt

