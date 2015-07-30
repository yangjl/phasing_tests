### Jinliang Yang
### July 23, 2015

##get command line args
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
#crossovers=as.numeric(args[1]) # mean expected crossovers per chromosome; 0 = no recombination
job=args[1]


#crossovers=1.5
### loading all the functions in folder "lib"
sourceall <-function(rm=FALSE){
    if(rm) rm(list=ls())
    f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)
}
sourceall(rm=TRUE)


for(i in 1:1){
    ### simulation perfect mom and noisy kids
    #set.seed(1235678*i)
    sim <- SimSelfer(size.array=10, het.error=0.7, hom.error=0.002, numloci=40000, rec=1.5, imiss=0.3)
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
    #plotphasing(sim, kids=1:5, snps=1:1000, cols=c("red", "blue"), plotphasing=TRUE, newmom)
    
    #KIDS GENOS
    inferred_progeny=list()
    mean.kid.geno.errors[mysim]=0;
    for(z in 1:length(progeny)){
        inferred_progeny[[z]]=which_phase_kid(newmom,progeny[[z]][[2]][estimated_hets] )
        mean.kid.geno.errors[mysim]=mean.kid.geno.errors[mysim]+(sum(abs(progeny[[z]][[1]][estimated_hets]-inferred_progeny[[z]])))/length(progeny)
    }
    out <- checkphasing(newmom, sim)
    write.table(out, paste0("out/out_",job,".txt"), sep="\t", row.names=FALSE, quote=FALSE)
}
