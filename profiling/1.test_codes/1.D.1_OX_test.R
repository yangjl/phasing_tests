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
sim <- SimOXer(size.array=10, het.error=0.8, hom.error=0.02, numloci=1000, rec=1.5, imiss=0.3, misscode = 3)
input <- simOX_input(sim, n_phased=5)
###>>> [[1]]: unphased dad (vector)
###>>> [[2]]: phased and unphased mom [ N=5+5 ] (list of data.frame + vector)
###>>> [[3]]: outcrossed progeny [ N=10 ] (list of list(real, obs))
###>>> [[4]]: pedigree (data.frame)

dad_geno <- input[[1]]
mom_array <- input[[2]]
progeny <- input[[3]]
ped <- input[[4]]
win_length=10
errors=c(0.02, 0.8)
verbose=TRUE

newdad <- phasingDad(dad_geno, mom_array, progeny, win_length, verbose=FALSE)
#plotphasing(sim, kids=1:5, snps=1:1000, cols=c("red", "blue"), plotphasing=TRUE, newmom)

pm <- write_mom(newmom)
save(file=paste0("largedata/testout/", job, "_phasemom.RData"), list=c("sim","pm"))

#ob <- load("largedata/lcached.RData")
#simk <- get_sim_kids(sim)
imputek <- imputing(momphase=pm, progeny, winstart=10, winend=500, stepsize=10, expect_recomb=1.5, verbose=TRUE)
rates <- kids_errs(simk=sim[[2]], imputek)

save(file=paste0("largedata/testout/", job, "_imputekid.RData"), list=c("sim","pm", "imputek", "rates"))


