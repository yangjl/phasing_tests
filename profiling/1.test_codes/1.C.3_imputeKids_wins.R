### Jinliang Yang
### August 19th, 2015

f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)
files <- list.files(path="largedata/testout/", pattern="_phasemom.RData", full.names=TRUE)


#ob <- load("largedata/lcached.RData")
#simk <- get_sim_kids(sim)
ob <- load(files[1])
probs <- get_error_mat(0.02, 0.8)[[2]]
input <- sim2input(sim)
progeny <- input[[2]]
imputek <- imputing(momphase=pm, progeny, winstart=10, winend=1000, stepsize=10, expect_recomb=1.5, verbose=TRUE)
rates <- kids_errs(simk=sim[[2]], imputek)




phaserate <- imputerate <- data.frame()
for(i in 1:length(files)){
    #rm(list=c("sim", "pm", "imputek"))
    ob <- load(files[i])
    merr <- mom_phasing_error(pm, sim)
    phaserate <- rbind(phaserate, merr)
    message(sprintf("###>>> Mom phasing error [ %s ]", merr$rate))
    
    kerr <- kids_errs(simk=sim[[2]], imputek)
    imputerate <- rbind(imputerate, kerr)
    message(sprintf("###>>> Kids phasing err [ %s ] and geno err [ %s ]", mean(kerr$phaserate), mean(kerr$genorate)))
}

write.table(phaserate, "data/phase_rate.csv", sep=",", row.names=FALSE, quote=FALSE)
write.table(imputerate, "data/impute_rate.csv", sep=",", row.names=FALSE, quote=FALSE)

####################
phaserate <- read.csv("data/phase_rate.csv")
imputerate <- read.csv("data/impute_rate.csv")

hist(phaserate$rate[phaserate$rate < 0.1], breaks=30, main="Simulation (N=100)", col="#faebd7", xlab="Phasing Error Rate")
abline(v=mean(phaserate$rate[phaserate$rate < 0.1]), col="red", lwd=2)
abline(v=median(phaserate$rate[phaserate$rate < 0.1]), col="darkblue", lwd=2)

par(mfrow=c(1,2))
hist(imputerate$phaserate, breaks=30, main="Kids Phasing Error", col="#faebd7", xlab="Error")
abline(v=mean(imputerate$phaserate), col="red", lwd=2)
abline(v=median(imputerate$phaserate), col="darkblue", lwd=2)

hist(imputerate$genorate, breaks=30, main="Kids Imputing Error", col="#faebd7", xlab="Error")
abline(v=mean(imputerate$genorate), col="red", lwd=2)
abline(v=median(imputerate$genorate), col="darkblue", lwd=2)



