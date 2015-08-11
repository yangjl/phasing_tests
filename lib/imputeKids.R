

#KIDS GENOS
inferred_progeny=list()
mean.kid.geno.errors[mysim]=0;
for(z in 1:length(progeny)){
    inferred_progeny[[z]]=which_phase_kid(newmom,progeny[[z]][[2]][estimated_hets] )
    mean.kid.geno.errors[mysim]=mean.kid.geno.errors[mysim]+(sum(abs(progeny[[z]][[1]][estimated_hets]-inferred_progeny[[z]])))/length(progeny)
}
mean.kid.geno.errors[mysim]=mean.kid.geno.errors[mysim]/numloci


##############################
imputing <- function(momphase, progeny, win_length){
    
    for(k in 1:length(progeny)){
        kid <- progeny[[k]][[2]]
        for(c in unique(momphase$chunk)){
            mychunk <- subset(momphase, chunk == c)
            
            #winstart <- i <- 1
            ### kids haplotypes
            mychunk$k1 <- 3
            mychunk$k2 <- 3
            if(win_length >= nrow(mychunk)){
                khaps <- which_phase_kid(haplotype, kidwin=kid[mychunk$idx])
            }else{
                for(win in 1:round(nrow(mychunk)/win_length,0) ){
                    if(verbose){ message(sprintf(">>> imputing kid [ %s ]: block [ %s/%s ] window [ %s/%s ] ...", 
                                                 k, c, length(unique(momphase$chunk)), win, nrow(mychunk))) } 
                    
                    myidx <- ((win-1)*win_length+1) : (win*win_length)
                    khaps <- which_phase_kid(haplotype, kidwin=kid[myidx])
                }
                ##### calculate last window
                myidx <- (nrow(mychunk)-win_length+1) : nrow(mychunk)
                khaps <- which_phase_kid(haplotype, kidwin=kid[myidx])
            }
        }    
    }
    
}


############################################################################
# Same as above, output kid's phase.
# give this mom haplotype and a kid's diploid genotype over the window and returns maximum prob
# Mendel is takenh care of in the probs[[]] matrix already
#  chunk idx hap1 hap2
#1     1   2    1    0
#2     1  19    1    0
#3     1  20    0    1
#4     1  21    0    1
#5     1  24    1    0
#6     1  28    0    1
which_phase_kid <- function(haplotype, kidwin){
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
    if(length(which(geno_probs==max(geno_probs)))==1){
        return(which.max(geno_probs))
    }else{
        return(NULL)
    }
}



### return the two haplotypes
#myh1 <- replace(estimated_mom/2, hetsites, mom_phase1)
#myh2 <- replace(estimated_mom/2, hetsites, 1-mom_phase1)
#return(data.frame(h1=myh1, h2=myh2))
#if(verbose){ message(sprintf(">>> phasing done!")) }
#haplist[[i]] <- list(mom_phase1, mom_phase2, hetsites[idxstart:length(hetsites)])
## list: hap1, hap2 and idx; info




############################################################################
### implement Viterbi algorithm for HMM
viterbi <- function(){
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
}
