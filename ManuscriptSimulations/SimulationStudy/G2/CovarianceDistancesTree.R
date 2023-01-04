
TrueMatrix = data.matrix(read.csv("TemporaryFiles/TrueCov.txt", header = F))

##################################################TreeMix
vcf = readLines("TemporaryFiles/adbayesinput.txt")[-1]
if (length(grep("0,0", vcf)) > 0) vcf = vcf[-grep("0,0", vcf)]

allelematrix = matrix(0, nrow = length(vcf), ncol = length((strsplit( vcf[2]  ,  ",| ")[[1]]))/2)

for (i in 1:nrow(allelematrix)) {
  temp  = as.numeric((strsplit( vcf[i]  ,  ",| ")[[1]]))
  for (j in 1:(length(temp)/2)) {
    allelematrix[i,j] = temp[j+j-1]/(temp[j+j-1] + temp[j+j])
  }
}
NewEmpiricalMatrix = data.matrix(read.csv("TemporaryFiles/TREECov.txt", header = F)) / 
                      (mean(rowMeans(allelematrix)*(1-rowMeans(allelematrix))))

TreeMIXResult = sqrt(sum((NewEmpiricalMatrix - TrueMatrix)**2))
###########
NewEmpiricalMatrix = data.matrix(read.csv("TemporaryFiles/ORIENTCov.txt", header = F)) / 
                      (mean(rowMeans(allelematrix)*(1-rowMeans(allelematrix))))

ORIENTResult = sqrt(sum((NewEmpiricalMatrix - TrueMatrix)**2))
###########
write.table(data.frame(ORIENTResult , TreeMIXResult ), file = "TemporaryFiles/CovarianceDistance.txt",
            row.names = F, quote = F)
