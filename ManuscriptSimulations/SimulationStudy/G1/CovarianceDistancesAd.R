args = as.numeric(commandArgs(trailingOnly=TRUE)[1]) 
for (ii in 1:args) {
TrueMatrix = data.matrix(read.csv("TemporaryFiles/TrueCov1.txt", header = F))
meandistance = c()
for (j in 1:100) {
EmpiricalMAtrix = data.matrix(read.csv(paste0("TemporaryFiles/Cov" ,toString(j), "_", toString(ii), ".txt"), header = F))
meandistance = c(meandistance, sqrt(sum((EmpiricalMAtrix - TrueMatrix)**2)))
}

MEAN = mean(meandistance)
MAP = sqrt(sum((data.matrix(read.csv(paste0("TemporaryFiles/MAPCov", toString(ii), ".txt"), header = F)) - TrueMatrix)**2))

write.table(data.frame( MAP, MEAN ), file = paste0("TemporaryFiles/CovarianceDistance", toString(ii), ".txt"),
            row.names = F, quote = F)

}
