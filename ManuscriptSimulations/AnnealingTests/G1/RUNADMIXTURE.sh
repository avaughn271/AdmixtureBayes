NUMBEROFSIMS=3   #PLACE 7

###################################################run admixturebayes
for (( c=1; c<=NUMBEROFSIMS; c++ ))
do

python3 ~/desktop/AdmixtureBayes-Annealing/admixturebayes/runMCMC.py --input_file TemporaryFiles/adbayesinput${c}.txt --outgroup pop5  --save_covariance --starting_temp 100 --ending_temp 0.0001  --temp_scaling 0.95  --iter_per_temp 5000

####python GetAdBayesResults.py
mv covariance_matrix.txt TemporaryFiles/covariance_matrix${c}.txt 

cp MAPadd.txt TemporaryFiles/MAPadd${c}.txt 
cp MAPtree.txt TemporaryFiles/MAPtree${c}.txt

done

Rscript CalculateSetDistanceAd.R $NUMBEROFSIMS

for (( c=1; c<=NUMBEROFSIMS; c++ ))
do

cp TemporaryFiles/MAPTree${c}.txt TemporaryFiles/MAPTree.txt

Rscript CheckTopologyEqualityAd.R
cp TemporaryFiles/TopEq.txt TemporaryFiles/TopEq${c}.txt
done

python3 getCovarianceAd.py $NUMBEROFSIMS
#python3 getCovarianceTree.py $NUMBEROFSIMS

Rscript CovarianceDistancesAd.R $NUMBEROFSIMS 


rm -f AdBayesSet.txt
touch AdBayesSet.txt
rm -f AdBayesTop.txt
touch AdBayesTop.txt

rm -f AdBayesCov.txt
touch AdBayesCov.txt

for (( c=1; c<=NUMBEROFSIMS; c++ ))
do

sed '2,2!d' TemporaryFiles/SetDistance${c}.txt > temp.txt
cat AdBayesSet.txt temp.txt > tempp.txt 
mv tempp.txt AdBayesSet.txt

sed '2,2!d' TemporaryFiles/CovarianceDistance${c}.txt > temp.txt
cat AdBayesCov.txt temp.txt > tempp.txt 
mv tempp.txt AdBayesCov.txt 

sed '2,2!d' TemporaryFiles/TopEq${c}.txt > temp.txt
cat AdBayesTop.txt temp.txt > tempp.txt 
mv tempp.txt AdBayesTop.txt 

done

rm temp.txt

rm MAPadd.txt
rm MAPtree.txt

Rscript MakePlots.R

rm -r -f TemporaryFiles/TopEq*.txt
rm -r -f TemporaryFiles/CovarianceDistance*.txt
rm -r -f TemporaryFiles/SetDistance*.txt