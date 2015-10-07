### Jinliang Yang
### August 7th, 2015

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
#crossovers=as.numeric(args[1]) # mean expected crossovers per chromosome; 0 = no recombination
job=args[1]


### loading all the functions in folder "lib"
f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)

### simulation perfect mom and noisy kids
#set.seed(1235678*i)
sim <- SimSelfer(size.array=10, het.error=0.8, hom.error=0.02, numloci=40000, rec=1.5, imiss=0.3)
#plotselfer(sim, kids=6:10, snps=40:100, cols=c("green", "blue"))
### format to phase_mom
input <- sim2input(sim)

#estimated_mom <- input[[1]]
progeny <- input[[2]]
#MOM PHASE
win_length=10
verbose=TRUE
probs <- get_error_mat(0.02, 0.8)[[2]]

newmom <- phasing(estimated_mom=input[[1]], progeny, win_length, verbose=TRUE)
#plotphasing(sim, kids=1:5, snps=1:1000, cols=c("red", "blue"), plotphasing=TRUE, newmom)

pm <- write_mom(newmom)
save(file=paste0("largedata/testout/", job, "_phasemom.RData"), list=c("sim","pm"))

#ob <- load("largedata/lcached.RData")
#simk <- get_sim_kids(sim)
imputek <- imputing(momphase=pm, progeny, winstart=10, winend=500, stepsize=10, expect_recomb=1.5, verbose=TRUE)
rates <- kids_errs(simk=sim[[2]], imputek)

save(file=paste0("largedata/testout/", job, "_imputekid.RData"), list=c("sim","pm", "imputek", "rates"))


