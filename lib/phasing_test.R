phasing_test <- function(outfile, numloci){
    ### simulation perfect mom and noisy kids
    #set.seed(1235678*i)
    #sim <- SimSelfer(size.array=10, het.error=0.7, hom.error=0.002, numloci=40000, rec=1.5, imiss=0.3)
    sim <- SimSelfer(size.array=10, het.error=0.5, hom.error=0.02, numloci, rec=1.5, imiss=0.3)
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
    
    #### phasing results
    outdf <- write_mom(newmom)
    names(outdf)[3:4] <- c("ihap1", "ihap2")
    
    #### error estimation
    outall <- cbind(outdf, sim[[1]][outdf$idx, ])
    err = 0
    for(i in unique(outall$chunk)){
        chunk <- subset(outall, chunk == i)
        idx1 <- which.max(c(cor(chunk$ihap1, chunk$hap1), cor(chunk$ihap1, chunk$hap2)) )
        err1 <- sum(chunk$ihap1 != chunk[, 4+idx1])
        idx2 <- which.max(c(cor(chunk$ihap2, chunk$hap1), cor(chunk$ihap2, chunk$hap2)) )
        err2 <- sum(chunk$ihap2 != chunk[, 4+idx2])
        err <- err + err1 + err2
    }
    out <- data.frame(chunk=length(unique(outall$chunk)), err=err, loci=2*nrow(outall))
    
    write.table(out, paste0("largedata/out/", outfile, "_out.txt"), sep="\t", row.names=FALSE, quote=FALSE )
    
}
