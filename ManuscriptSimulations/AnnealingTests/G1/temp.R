
B = read.table("mcmc_samples.csv", header = T, sep = ",")


plot(B$posterior, type = "l")

lines(B$prior, type = "l", col = "green")

lines(B$likelihood, type = "l", col=  "blue")



plot(B$posterior, type = "l", ylim = c(-1000,80))

lines(B$prior, type = "l", col = "green")

lines(B$likelihood, type = "l", col=  "blue")

plot(B$posterior, type = "l", ylim = c(25.3,25.4))

lines(B$prior, type = "l", col = "green")

lines(B$likelihood, type = "l", col=  "blue")
