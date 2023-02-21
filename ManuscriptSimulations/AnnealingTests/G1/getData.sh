NUMBEROFSIMS=3   #PLACE 5
rm -r -f TemporaryFiles
mkdir TemporaryFiles
for (( c=1; c<=NUMBEROFSIMS; c++ ))
do
python run_simulation.py
Rscript ConvertData.R
mv TemporaryFiles/Data.vcf TemporaryFiles/Data${c}.vcf
#mv TemporaryFiles/TrueAdd.txt TemporaryFiles/TrueAdd${c}.txt
#mv TemporaryFiles/TrueTree.txt TemporaryFiles/TrueTree${c}.txt
mv TemporaryFiles/adbayesinput.txt TemporaryFiles/adbayesinput${c}.txt
done

python3 getcovarianceTree.py


rm -r -f TemporaryFiles/Data*.vcf