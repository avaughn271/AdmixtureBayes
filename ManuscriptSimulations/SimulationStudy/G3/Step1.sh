NUMBEROFSIMS=20   #PLACE 5
rm -r TemporaryFiles
mkdir TemporaryFiles
for (( c=1; c<=NUMBEROFSIMS; c++ ))
do
python run_simulation.py
Rscript ConvertData.R
mv TemporaryFiles/Data.vcf TemporaryFiles/Data${c}.vcf
mv TemporaryFiles/TrueAdd.txt TemporaryFiles/TrueAdd${c}.txt
mv TemporaryFiles/TrueTree.txt TemporaryFiles/TrueTree${c}.txt
mv TemporaryFiles/adbayesinput.txt TemporaryFiles/adbayesinput${c}.txt
done
