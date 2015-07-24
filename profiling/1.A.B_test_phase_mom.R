### Jinliang Yang
### July 23, 2015

### loading all the functions in folder "lib"
sourceall <-function(rm=FALSE){
    if(rm) rm(list=ls())
    f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)
}
sourceall(rm=F)

### simulation perfect mom and noisy kids
set.seed(123456)
sim <- SimSelfer(size.array=10, het.error=0.7, hom.error=0.002, numloci=100, rec=1.2, imiss=0.3)
plotselfer(sim, kids=6:10, snps=40:100, cols=c("green", "blue"))


### format to phase_mom
input <- sim2input(sim)
probs <- get_error_mat(0.002, 0.7)[[2]]
estimated_mom <- input[[1]]
progeny <- input[[2]]
#MOM PHASE
newmom <- phase_mom(estimated_mom, progeny, 10, verbose=TRUE)


plotphasing(sim, kids=1:10, snps=1:100, cols=c("red", "blue"), plotphasing=TRUE, newmom)









estimated_hets=which(estimated_mom==1)
# can't make a phasing error at a site which is not really heterozygous, 
# nor is calling a true het site homozygous a phasing error
#  so we only check phasing error at the intersection of both
# note that this could in theory lead to super low phasing error BECAUSE of high
# genotype error, but we're gonna ignore that for now
phase_sites=intersect(which(true_mom[[1]]+true_mom[[2]]==1),estimated_hets) 
mom.phase.errors[mysim]=sum(sapply(2:length(phase_sites), 
                                   function(a) check_phase(estimated_hets,newmom,true_mom,phase_sites[a],phase_sites[a-1])))

