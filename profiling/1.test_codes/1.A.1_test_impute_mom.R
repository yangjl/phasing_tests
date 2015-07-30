### loading all the functions in folder "lib"
sourceall <-function(rm=FALSE){
    if(rm) rm(list=ls())
    f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)
}
sourceall(rm=T)

set.seed(123456)
sim <- SimSelfer(size.array=10, het.error=0.7, hom.error=0.002, numloci=100, rec=1.2, imiss=0.3)

plotselfer(sim, kids=6:10, snps=40:100, cols=c("green", "blue"))







