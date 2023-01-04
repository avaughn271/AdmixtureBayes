NUMBEROFSIMS=19

Rscript TreeMixConvertOutputSpecial.R pop5
####################################################run treemix and orientagraph


for (( c=1; c<=NUMBEROFSIMS; c++ ))
do
echo $c
cp TemporaryFiles/TREEMIXoutput${c}.txt TemporaryFiles/TREEMIXoutput.txt
cp TemporaryFiles/TrueTree1.txt TemporaryFiles/TrueTree.txt
cp TemporaryFiles/ORIENToutput${c}.txt TemporaryFiles/ORIENToutput.txt
cp TemporaryFiles/ORIENToutgroup${c}.txt TemporaryFiles/ORIENToutgroup.txt
cp TemporaryFiles/TREEMIXoutgroup${c}.txt TemporaryFiles/TREEMIXoutgroup.txt
cp TemporaryFiles/TrueAdd1.txt TemporaryFiles/TrueAdd.txt
cp TemporaryFiles/TrueAdd1.txt TemporaryFiles/TrueAdd.txt

Rscript CalculateSetDistanceTree.R
Rscript CheckTopologyEqualityTree.R
python3 getCovarianceTree.py
mv "TemporaryFiles/TopEq.txt"  "TemporaryFiles/TopEq${c}.txt"
mv "TemporaryFiles/SetDistance.txt"  "TemporaryFiles/SetDistance${c}.txt"

cp TemporaryFiles/TREECov.txt TemporaryFiles/TREECov${c}.txt
cp TemporaryFiles/ORIENTCov.txt TemporaryFiles/ORIENTCov${c}.txt
cp TemporaryFiles/TrueCov.txt TemporaryFiles/TrueCov${c}.txt


done

for (( c=1; c<=NUMBEROFSIMS; c++ ))
do

cp TemporaryFiles/adbayesinput1.txt TemporaryFiles/adbayesinput.txt #edited this  
cp TemporaryFiles/TREECov${c}.txt TemporaryFiles/TREECov.txt
cp TemporaryFiles/ORIENTCov${c}.txt TemporaryFiles/ORIENTCov.txt
cp TemporaryFiles/TrueCov${c}.txt TemporaryFiles/TrueCov.txt

Rscript CovarianceDistancesTree.R
mv "TemporaryFiles/CovarianceDistance.txt"  "TemporaryFiles/CovarianceDistance${c}.txt"

done



rm -f TreemixSet.txt
touch TreemixSet.txt
rm -f TreemixTop.txt
touch Treemixtop.txt

rm -f TreemixCov.txt
touch TreemixCov.txt


for (( c=1; c<=NUMBEROFSIMS; c++ ))
do

sed '2,2!d' TemporaryFiles/SetDistance${c}.txt > temp.txt
cat TreemixSet.txt temp.txt > tempp.txt 
mv tempp.txt TreemixSet.txt

sed '2,2!d' TemporaryFiles/CovarianceDistance${c}.txt > temp.txt
cat TreemixCov.txt temp.txt > tempp.txt 
mv tempp.txt TreemixCov.txt 

sed '2,2!d' TemporaryFiles/TopEq${c}.txt > temp.txt
cat TreemixTop.txt temp.txt > tempp.txt 
mv tempp.txt TreemixTop.txt 

done

rm temp.txt
