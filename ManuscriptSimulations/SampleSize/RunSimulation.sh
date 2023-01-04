mkdir Small
mkdir Large

n=100

for (( i=1 ; i<=$n ; i++ )); 
do

python Small.py
Rscript ConvertData.R

mv adbayesinput.txt Small/adbayesinput${i}.txt 

done

for (( i=1 ; i<=$n ; i++ )); 
do

python3 admixturebayes/runMCMC.py --input_file  Small/adbayesinput${i}.txt --n 3000 --outgroup pop4
python3 admixturebayes/analyzeSamples.py  --mcmc_results  mcmc_samples.csv  --burn_in_fraction 0.15  --thinning_rate 10
python3 admixturebayes/makePlots.py  --plot top_trees --posterior    thinned_samples.csv --write_rankings Small/${i}rankings.txt

done

for (( i=1 ; i<=$n ; i++ )); 
do
echo $i
python Large.py
Rscript ConvertData.R

mv adbayesinput.txt Large/adbayesinput${i}.txt 

done

rm Data.vcf

for (( i=1 ; i<=$n ; i++ )); 
do

python admixturebayes/runMCMC.py --input_file  Large/adbayesinput${i}.txt --n 3000 --outgroup pop4
python admixturebayes/analyzeSamples.py  --mcmc_results  mcmc_samples.csv  --burn_in_fraction 0.15  --thinning_rate 10
python admixturebayes/makePlots.py  --plot top_trees --posterior    thinned_samples.csv --write_rankings Large/${i}rankings.txt

done

rm mcmc_samples.csv
rm thinned_samples.csv
rm topology*.pdf

Rscript PlotAll.R