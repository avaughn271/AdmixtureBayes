#We read in the vcf
vcf = readLines("TemporaryFiles/Data.vcf")
SampleNameLine = grep("CHROM\tPOS\tID", vcf)
SampleNames = unlist(strsplit(vcf[SampleNameLine], "\t"))[-(1:9)]

Populations = rep( c("pop1", "pop2",  "pop3", "pop4" , "pop5"), each = 5) #PLACE 4
print(Populations)

NumberOfHaplotypes = 2 * length(SampleNames)
HaplotypePopulations = rep(Populations, each = 2)
UniquePops = sort(unique(Populations))

vcfdata = vcf[-(1:SampleNameLine)]
NumberOfSNPs = length(vcfdata)
Data = matrix(data = -1, nrow = NumberOfSNPs, ncol = NumberOfHaplotypes)

for (snpindex in 1:NumberOfSNPs) {
  snpdata = unlist(strsplit(vcfdata[snpindex], c("\t")))[-(1:9)]
  for (indiv in 1:length(SampleNames)) {
    if (substr(snpdata[indiv], 1, 1) == "0") Data[snpindex, 2*indiv - 1] = 0
    if (substr(snpdata[indiv], 1, 1) == "1") Data[snpindex, 2*indiv - 1] = 1
    if (substr(snpdata[indiv], 3, 3) == "0") Data[snpindex, 2*indiv] = 0
    if (substr(snpdata[indiv], 3, 3) == "1") Data[snpindex, 2*indiv] = 1
  }
}

Output = matrix("NA", nrow = NumberOfSNPs, ncol = length(UniquePops))

for (i in 1:nrow(Output)) {
  for (pop in 1:length(UniquePops)) {
    relevantindices = which(HaplotypePopulations == UniquePops[pop])
    nderived = sum(Data[i, relevantindices] == 1)
    nancestral = sum(Data[i, relevantindices] == 0)
    Output[i,pop] = paste0(toString(nderived), "," , toString(nancestral))
  }
}

write.table(rbind(UniquePops, Output), "TemporaryFiles/adbayesinput.txt", 
                 quote = F, row.names = F, col.names = F)