### methods for the sim object

### format SimSelfer object to JRI's code
sim2input <- function(sim){
    simp <- sim[[1]] # mom
    simk <- sim[[2]] # kids
    progeny <- list()
    for(i in 1:length(simk)){
        ### list of true and observed kids
        progeny[[i]] <- list(simk[[i]][[2]]$hap1+simk[[i]][[2]]$hap2, simk[[i]][[2]]$obs)
    }
    #p <- frq(progeny)
    ### use the perfect parental genotype
    return(list(simp$geno, progeny))
}

frq <- function(progeny){
    res <- 0
    for(i in 1:length(progeny)){
        res <- res + progeny[[i]][[2]]
    }
    res <- res/(2*length(progeny))
    res <- replace(res, which(res==0), 0.5/length(progeny))
    return(res)
}

get_sim_kids <- function(sim){
    simk <- sim[[2]]
    progeny <- list()
    for(i in 1:length(simk)){
        progeny[[i]] <- simk[[i]][[2]]
    }
    return(progeny)
} 