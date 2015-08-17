### Jinliang Yang
### August 7th, 2015


### write PLINK format file
#FormatData(wgs, cols=9:27)


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

out <- write_mom(newmom)


#########################################################################################################

#ob <- load("largedata/lcached.RData")
#simk <- get_sim_kids(sim)
imputek <- imputing(out, progeny, 15, verbose)

save(file=paste0("largedata/out/", job, "_imputekid.RData"), list=c("out", "progeny", "sim", "imputek"))

simk <- progeny
out10 <- comp_kids(simk=sim[[2]], imputek)

#[1] 0.06037394
###>>> Average error rate [ 0.0531363088057901 ]

###>>> Average error rate [ 0.0399276236429433 ]
