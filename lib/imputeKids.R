### implement Viterbi algorithm for HMM

############################################################################
# Same as above, output kid's phase.
# give this mom haplotype and a kid's diploid genotype over the window and returns maximum prob
# Mendel is takenh care of in the probs[[]] matrix already 
which_phase_kid <- function(haplotype,kidwin){
    three_genotypes=list()
    haplotype=unlist(haplotype)
    three_genotypes[[1]]=haplotype+haplotype
    three_genotypes[[2]]=haplotype+(1-haplotype)
    three_genotypes[[3]]=(1-haplotype)+(1-haplotype)
    geno_probs=as.numeric() #prob of each of three genotypes
    for(geno in 1:3){
        #log(probs[[2]][three_genotypes,kidwin] is the log prob. of kid's obs geno 
        #given the current phased geno and given mom is het. (which is why probs[[2]])
        geno_probs[geno]=sum( sapply(1:length(haplotype), function(zz) log( probs[[2]][three_genotypes[[geno]][zz]+1,kidwin[zz]+1])))
    }
    if(length(which(geno_probs==max(geno_probs)))!=1){recover()}
    return(three_genotypes[[which(geno_probs==max(geno_probs))]])
}
############################################################################

emit <- c()

init <- c(0.5, 0.5)


# backtracking algorithm
for (i in 2:length(symbol.sequence)) {
    # probability vector stores the current emission with respect to (i-1) observation of selected state and transition probability
    # state vector (pointer) on the other hand is only storing the most probable state in (i-1), which we will later use for backtracking
    
    
    #l and k => two states
    tmp.path.probability <- lapply(states, function(l) {
        max.k <- unlist(lapply(states, function(k) {
            prob.history[i-1, k] + transition.matrix[k, l]
        }))
        return(c(states[which(max.k == max(max.k))], max(max.k) + emission.matrix[symbol.sequence[i], l]))
    })
    
    prob.history <- rbind(prob.history, data.frame(F = as.numeric(tmp.path.probability[[1]][2]), 
                                                   L = as.numeric(tmp.path.probability[[2]][2])))
    
    state.history <- data.frame(F = c(as.character(state.history[,tmp.path.probability[[1]][1]]), "F"), 
                                L = c(as.character(state.history[,tmp.path.probability[[2]][1]]), "L"))
}

# selecting the most probable path
viterbi.path <- as.character(state.history[,c("F", "L")[which(max(prob.history[length(symbol.sequence), ]) == prob.history[length(symbol.sequence), ])]])