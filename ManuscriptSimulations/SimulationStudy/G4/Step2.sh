NUMBEROFSIMS=20   #PLACE 5
rm TemporaryFiles/tr*.edges
rm TemporaryFiles/tr*.vertices

####################################################run treemix and orientagraph
for (( c=1; c<=NUMBEROFSIMS; c++ ))
do
echo $c
cp TemporaryFiles/adbayesinput${c}.txt treemixinput.txt
gzip treemixinput.txt

treemix                                                              -i treemixinput.txt.gz -o TemporaryFiles/treeoutput${c}   -k 1000 -global                     -root pop7 -m 2

# ~/desktop/AdmixtureBayesOWork/OtherMethods/OrientaGraph/src/orientagraph -i treemixinput.txt.gz -o TemporaryFiles/orientoutput${c} -k 1000 -global -mlno #-allmigs -root pop7 -m 2  #PLACE 6 END

gunzip TemporaryFiles/treeoutput${c}.vertices.gz
gunzip TemporaryFiles/treeoutput${c}.edges.gz
#gunzip TemporaryFiles/orientoutput${c}.vertices.gz
#gunzip TemporaryFiles/orientoutput${c}.edges.gz

rm treemixinput.txt.gz

done
Rscript TreeMixCheckOutput.R pop7 $NUMBEROFSIMS ### remember that the script itself has to be edited for the outgroup location
