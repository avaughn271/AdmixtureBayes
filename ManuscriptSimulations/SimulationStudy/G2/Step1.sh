NUMBEROFSIMS=20   #PLACE 5
rm -r TemporaryFiles
mkdir TemporaryFiles
##############get data
for (( c=1; c<=NUMBEROFSIMS; c++ ))
do
python run_simulation.py
Rscript ConvertData.R
mv TemporaryFiles/Data.vcf TemporaryFiles/Data${c}.vcf
mv TemporaryFiles/TrueAdd.txt TemporaryFiles/TrueAdd${c}.txt
mv TemporaryFiles/TrueTree.txt TemporaryFiles/TrueTree${c}.txt
mv TemporaryFiles/adbayesinput.txt TemporaryFiles/adbayesinput${c}.txt
done

####################################################run treemix and orientagraph
for (( c=1; c<=NUMBEROFSIMS; c++ ))
do

cp TemporaryFiles/adbayesinput${c}.txt treemixinput.txt
gzip treemixinput.txt

  #PLACE 6
outgroupp="pop6"
treemix                                                              -i treemixinput.txt.gz -o TemporaryFiles/treeoutput${c}   -k 1000 -global                     -root $outgroupp -m 0

 ~/desktop/AdmixtureBayesOWork/OtherMethods/OrientaGraph/src/orientagraph -i treemixinput.txt.gz -o TemporaryFiles/orientoutput${c} -k 1000 -global  -mlno 1,2 -allmigs 1,2   -root $outgroupp -m 0   #PLACE 6 END

gunzip TemporaryFiles/treeoutput${c}.vertices.gz
gunzip TemporaryFiles/treeoutput${c}.edges.gz
gunzip TemporaryFiles/orientoutput${c}.vertices.gz
gunzip TemporaryFiles/orientoutput${c}.edges.gz

rm treemixinput.txt.gz

done
Rscript TreeMixConvertOutput.R $outgroupp $NUMBEROFSIMS
####################################################run treemix and orientagraph


for (( c=1; c<=NUMBEROFSIMS; c++ ))
do
cp TemporaryFiles/TREEMIXoutput${c}.txt TemporaryFiles/TREEMIXoutput.txt
cp TemporaryFiles/TrueTree${c}.txt TemporaryFiles/TrueTree.txt
cp TemporaryFiles/ORIENToutput${c}.txt TemporaryFiles/ORIENToutput.txt
cp TemporaryFiles/ORIENToutgroup${c}.txt TemporaryFiles/ORIENToutgroup.txt
cp TemporaryFiles/TREEMIXoutgroup${c}.txt TemporaryFiles/TREEMIXoutgroup.txt
cp TemporaryFiles/TrueAdd${c}.txt TemporaryFiles/TrueAdd.txt
cp TemporaryFiles/TrueAdd${c}.txt TemporaryFiles/TrueAdd.txt

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

cp TemporaryFiles/adbayesinput${c}.txt TemporaryFiles/adbayesinput.txt
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
