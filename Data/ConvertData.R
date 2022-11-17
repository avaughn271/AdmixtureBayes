#We read in the vcf
vcf = readLines("ArcticData.vcf") #REPLACE THIS WITH ACTUAL VCF

#And find the line number with the sample names
SampleNameLine = grep("CHROM\tPOS\tID", vcf)

#We extract the sample names
SampleNames = unlist(strsplit(vcf[SampleNameLine], "\t"))[-(1:9)]

#We declare which sample names fall into which populations
#The ordering of the populations must exactly match the order of
#sample names in the VCF file
Populations =  c("Athabascan", "Athabascan", "Anzick", "Greenlander", "Greenlander", "Han", "Yoruba", "Ket", "Ket", "Koryak", "Koryak", "Malta", "Saqqaq", "Aymara", "USR1", "UstIshim")

#We extract the genotype data
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

#And match the data to each population

Output = matrix("NA", nrow = NumberOfSNPs, ncol = length(UniquePops))

for (i in 1:nrow(Output)) {
  for (pop in 1:length(UniquePops)) {
    relevantindices = which(HaplotypePopulations == UniquePops[pop])
    nderived = sum(Data[i, relevantindices] == 1)
    nancestral = sum(Data[i, relevantindices] == 0)
    Output[i,pop] = paste0(toString(nderived), "," , toString(nancestral))
  }
}

write.table(rbind(UniquePops, Output), "adbayesinput.txt", quote = F, row.names = F, col.names = F)

A = readLines("adbayesinput.txt") # We remove sites that have missing data
writeLines(A[setdiff(1:length(A), grep("0,0", A))], "ArcticData.txt")
file.remove("adbayesinput.txt")
