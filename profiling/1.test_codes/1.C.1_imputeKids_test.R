### Jinliang Yang
### August 7th, 2015


### write PLINK format file
#FormatData(wgs, cols=9:27)


### loading all the functions in folder "lib"
f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)

### simulation perfect mom and noisy kids
#set.seed(1235678*i)
sim <- SimSelfer(size.array=10, het.error=0.5, hom.error=0.02, numloci=10000, rec=1.5, imiss=0.3)
#plotselfer(sim, kids=6:10, snps=40:100, cols=c("green", "blue"))
### format to phase_mom
input <- sim2input(sim)

#estimated_mom <- input[[1]]
progeny <- input[[2]]
#MOM PHASE
win_length=10
verbose=TRUE
probs <- get_error_mat(0.02, 0.5)[[2]]

newmom <- phasing(estimated_mom=input[[1]], progeny, win_length, verbose=TRUE)
#plotphasing(sim, kids=1:5, snps=1:1000, cols=c("red", "blue"), plotphasing=TRUE, newmom)

out <- write_mom(newmom)

save(file="largedata/lcached.RData", list=c("out", "progeny", "sim"))
hetsites <- which(input[[1]]==1)

#########################################################################################################
p <- imputing(out, progeny, 10, verbose)


