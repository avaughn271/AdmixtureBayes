NUMBEROFSIMS=20   #PLACE 7

###################################################run admixturebayes
for (( c=1; c<=NUMBEROFSIMS; c++ ))
do

python3 ~/desktop/AdmixtureBayes-main/admixturebayes/runMCMC.py --input_file TemporaryFiles/adbayesinput${c}.txt --outgroup pop5 --n 60000 --save_covariance    #PLACE 8 and last

python GetAdBayesResults.py
mv covariance_matrix.txt TemporaryFiles/covariance_matrix${c}.txt 

cp MAPadd.txt TemporaryFiles/MAPadd${c}.txt 
cp MAPtree.txt TemporaryFiles/MAPtree${c}.txt

for (( d=1; d<=100; d++ ))
do
cp TemporaryFiles/Tree${d}.txt TemporaryFiles/Tree${d}_${c}.txt 
cp TemporaryFiles/add${d}.txt TemporaryFiles/add${d}_${c}.txt 

done

done

Rscript CalculateSetDistanceAd.R $NUMBEROFSIMS

for (( c=1; c<=NUMBEROFSIMS; c++ ))
do

cp TemporaryFiles/TrueTree${c}.txt TemporaryFiles/TrueTree.txt
cp TemporaryFiles/MAPTree${c}.txt TemporaryFiles/MAPTree.txt

for (( d=1; d<=100; d++ ))
do
cp TemporaryFiles/Tree${d}_${c}.txt TemporaryFiles/Tree${d}.txt
done

Rscript CheckTopologyEqualityAd.R
cp TemporaryFiles/TopEq.txt TemporaryFiles/TopEq${c}.txt
done

python3 getCovarianceAd.py $NUMBEROFSIMS
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
rm mcmc_samples.csv

rm MAPadd.txt
rm MAPtree.txt


