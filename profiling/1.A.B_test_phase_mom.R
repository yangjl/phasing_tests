### Jinliang Yang
### July 23, 2015

### loading all the functions in folder "lib"
sourceall <-function(rm=FALSE){
    if(rm) rm(list=ls())
    f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)
}
sourceall(rm=TRUE)

### simulation perfect mom and noisy kids
set.seed(1235)
sim <- SimSelfer(size.array=10, het.error=0.7, hom.error=0.002, numloci=5000, rec=1.5, imiss=0.3)
#plotselfer(sim, kids=6:10, snps=40:100, cols=c("green", "blue"))


### format to phase_mom
input <- sim2input(sim)
probs <- get_error_mat(0.002, 0.7)[[2]]
estimated_mom <- input[[1]]
progeny <- input[[2]]
win_length=10
verbose=TRUE
#MOM PHASE
newmom <- phasing(estimated_mom, progeny, win_length, verbose=TRUE)
plotphasing(sim, kids=1:5, snps=1:1000, cols=c("red", "blue"), plotphasing=TRUE, newmom)



truehap <- sim[[1]][newmom[[i]][[3]],]
esthap <- data.frame(h1=newmom[[i]][[1]], h2=newmom[[i]][[2]])
tab <- cbind(truehap, esthap)
which.max(c(cor(truehap$hap1, esthap$h1), cor(truehap$hap1, esthap$h2)))




set.seed(1235)
sim <- SimSelfer(size.array=10, het.error=0.7, hom.error=0.002, numloci=40000, rec=1.5, imiss=0.3)


