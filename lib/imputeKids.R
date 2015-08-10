

#KIDS GENOS
inferred_progeny=list()
mean.kid.geno.errors[mysim]=0;
for(z in 1:length(progeny)){
    inferred_progeny[[z]]=which_phase_kid(newmom,progeny[[z]][[2]][estimated_hets] )
    mean.kid.geno.errors[mysim]=mean.kid.geno.errors[mysim]+(sum(abs(progeny[[z]][[1]][estimated_hets]-inferred_progeny[[z]])))/length(progeny)
}
mean.kid.geno.errors[mysim]=mean.kid.geno.errors[mysim]/numloci

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

#############################
#for(winstart in 1:(length(hetsites)-(win_length-1)))
winstart <- i <- 1
while(winstart <= length(hetsites)-(win_length-1)){
    if(verbose){ message(sprintf(">>> imputing window [ %s / %s ] ...", winstart, length(hetsites))) } 
    kidwin <- hetsites[winstart:(winstart+win_length-1)]
    if(winstart==1){ 
        #arbitrarily assign win_hap to one chromosome initially
        win_hap <- infer_dip(momwin,progeny,haps=mom_haps, returnhap=TRUE)
        kid_phase1 <- win_hap[[1]]
        kid_phase2 <- win_hap[[2]]
        idxstart <- 1
    } else{
        win_hap <- infer_dip(momwin, progeny, haps=mom_haps, returnhap=FALSE)
        ### comparing current hap with old hap except the last bp -JLY
        if(!is.null(win_hap)){
            
            same=sum(mom_phase1[(length(mom_phase1)-win_length+2):length(mom_phase1)]==win_hap[1:length(win_hap)-1])
            
            if(same == 0){ #totally opposite phase of last window
                mom_phase2[length(mom_phase2)+1] <- win_hap[length(win_hap)]
                mom_phase1[length(mom_phase1)+1] <- 1-win_hap[length(win_hap)]
            } else if(same==(win_length-1) ){ #same phase as last window
                mom_phase1[length(mom_phase1)+1] <- win_hap[length(win_hap)]
                mom_phase2[length(mom_phase2)+1] <- 1-win_hap[length(win_hap)]
            } else{
                diff1 <- sum(abs(mom_phase1[(length(mom_phase1)-win_length+2):length(mom_phase1)]-win_hap[1:length(win_hap)-1]))
                diff2 <- sum(abs(mom_phase2[(length(mom_phase1)-win_length+2):length(mom_phase1)]-win_hap[1:length(win_hap)-1]))
                if(diff1 > diff2){ #momphase1 is less similar to current inferred hap
                    mom_phase2[length(mom_phase2)+1] <- win_hap[length(win_hap)]
                    mom_phase1[length(mom_phase1)+1] <- 1-win_hap[length(win_hap)]
                } else{ #momphase1 is more similar
                    mom_phase1[length(mom_phase1)+1] <- win_hap[length(win_hap)]
                    mom_phase2[length(mom_phase2)+1] <- 1-win_hap[length(win_hap)]
                }
            }
        } else {
            ### potential recombination in kids, output previous haps and jump to next non-overlap window -JLY###
            idxend <- winstart + win_length -2
            haplist[[i]] <- list(mom_phase1, mom_phase2, hetsites[idxstart:idxend])
            i <- i +1
            
            ### warning(paste("Likely recombination at position", winstart+1, sep=" "))
            ### if new window is still ambiguous, add 1bp and keep running until find the best hap
            winstart <- winstart + win_length -2
            while(is.null(win_hap)){
                
                winstart <- winstart + 1
                win_hap <- jump_win(winstart, win_length, hetsites, mom_haps)
                if(is.null(win_hap)){
                    nophase <- c(nophase, hetsites[winstart])
                }
            }
            idxstart <- winstart
            mom_phase1 <- win_hap
            mom_phase2 <- 1-win_hap
        }
    }
    winstart <- winstart + 1
}

### return the two haplotypes
#myh1 <- replace(estimated_mom/2, hetsites, mom_phase1)
#myh2 <- replace(estimated_mom/2, hetsites, 1-mom_phase1)
#return(data.frame(h1=myh1, h2=myh2))
#if(verbose){ message(sprintf(">>> phasing done!")) }
haplist[[i]] <- list(mom_phase1, mom_phase2, hetsites[idxstart:length(hetsites)])
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
