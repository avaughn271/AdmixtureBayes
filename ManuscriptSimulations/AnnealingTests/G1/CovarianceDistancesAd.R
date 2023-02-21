args = as.numeric(commandArgs(trailingOnly=TRUE)[1]) 
for (ii in 1:args) {
TrueMatrix = data.matrix(read.csv("TemporaryFiles/TrueCov.txt", header = F))


MAP = sqrt(sum((data.matrix(read.csv(paste0("TemporaryFiles/MAPCov", toString(ii), ".txt"), header = F)) - TrueMatrix)**2))

write.table( MAP, file = paste0("TemporaryFiles/CovarianceDistance", toString(ii), ".txt"),
            row.names = F, quote = F)

}
