############################################################################
# Get mom's phase
# should return two haps
phase_mom <- function(estimated_mom, progeny, win_length, verbose=FALSE){
    hetsites <- which(estimated_mom==1)
    # gets all possible haplotypes for X hets 
    mom_haps <- setup_haps2(win_length) 
    mom_phase1=as.numeric() 
    mom_phase2=as.numeric() 
    win_hap=as.numeric()
    old_hap=as.numeric() 
    for(winstart in 1:(length(hetsites)-(win_length-1))){
        if(verbose){ message(sprintf(">>> phasing window [ %s ] ...", winstart)) } 
        momwin <- hetsites[winstart:(winstart+win_length-1)]
        if(winstart==1){ 
            #arbitrarily assign win_hap to one chromosome initially
            win_hap=infer_dip(momwin,progeny,mom_haps)
            mom_phase1=win_hap
            mom_phase2=1-win_hap
        } else{
            win_hap=infer_dip(momwin,progeny,mom_haps,mom_phase1)
            ### comparing current hap with old hap?
            same=sum(mom_phase1[winstart:(winstart+(length(win_hap)-2))]==win_hap[1:length(win_hap)-1])
            if(same==0){ #totally opposite phase of last window
                mom_phase2[length(mom_phase2)+1] <- win_hap[length(win_hap)]
                mom_phase1[length(mom_phase1)+1] <- 1-win_hap[length(win_hap)]
            } else if(same==(win_length-1) ){ #same phase as last window
                mom_phase1[length(mom_phase1)+1] <- win_hap[length(win_hap)]
                mom_phase2[length(mom_phase2)+1] <- 1-win_hap[length(win_hap)]
            } else{ ##returns the minimum distance haplotype (prepare to screw up phase!)
                ### need to keep the current hap and recalculate over certain windows ? -JLY###
                warning(paste("Likely recombination at", winstart, sep=" "))
                diff1 <- sum(abs(mom_phase1[winstart:(winstart+(length(win_hap)-2))]-win_hap[1:length(win_hap)-1]))
                diff2 <- sum(abs(mom_phase2[winstart:(winstart+(length(win_hap)-2))]-win_hap[1:length(win_hap)-1]))
                if(diff1>diff2){ #momphase1 is less similar to current inferred hap
                    mom_phase2[length(mom_phase2)+1]=win_hap[length(win_hap)]
                    mom_phase1[length(mom_phase1)+1]=1-win_hap[length(win_hap)]
                } else{ #momphase1 is more similar
                    mom_phase1[length(mom_phase1)+1]=win_hap[length(win_hap)]
                    mom_phase2[length(mom_phase2)+1]=1-win_hap[length(win_hap)]
                }
            }
        }
    }
    myphase <- replace(estimated_mom/2, hetsites, mom_phase1)
    return(mom_phase1)
}
############################################################################


############################################################################
# Setup all possible haplotypes for window of X heterozgous sites
# This needs to be fixed to remove redundancy. E.g. 010 is the same as 101 and 1010 is same as 0101. 
# I don't think should bias things in the meantime, just be slow.
setup_haps <- function(win_length){
    haps=list(0,1); 
    for(i in 2:win_length){ 
        haps=c(haps,haps); 
        for(j in 1:(length(haps)/2)){  haps[[j]][(length(haps[[j]]))+1]=0 };   
        for(k in (length(haps)/2+1):length(haps)){  haps[[k]][(length(haps[[k]]))+1]=1 };   
    }
    nohaps=as.numeric();
    newhaps=list();
    for(i in 1:(length(haps)-1)){
        for(j in (i+1):length(haps)){
            if(sum((1-unlist(haps[j]))==unlist(haps[i]))==win_length){ nohaps[length(nohaps)+1]=i }
        }
    }
    for(i in 1:length(haps)){ if(!(i %in% nohaps)){newhaps[[length(newhaps)+1]]=haps[[i]]}}
    return(newhaps)
}
############################################################################
#system.time(tem <- setup_haps(10))
setup_haps2 <- function(win_length){
    if(win_length <= 20){
        alist <- lapply(1:win_length, function(a) c(0,1) )
        ### give a combination of all 0,1 into a data.frame
        hapdf <- expand.grid(alist)[1:2^(win_length-1),]
        ### split the data.frame into a list
        return(as.list(as.data.frame(t(hapdf)))) 
    }else{
        stop("!!! Can not handle [win_length > 20] !")
    }   
}
#system.time(tem2 <- setup_haps2(10))
#system.time(tem <- setup_haps(10))

############################################################################
# Infer which phase is mom in a window
infer_dip <- function(momwin,progeny,haps,momphase1){  
    #momwin is list of heterozygous sites, progeny list of kids genotypes, 
    #haps list of possible haps,momphase1 is current phased mom for use in splitting ties
    #phase_probs=as.numeric()
    #for(my_hap in 1:(length(haps))){ 
        #iterate over possible haplotypes <- this is slower because setup_haps makes too many haps
        #get max. prob for each kid, sum over kids
    #    phase_probs[my_hap]=sum( sapply(1:length(progeny), function(z) which_phase(haps[my_hap],progeny[[z]][[2]][momwin] )))
    #}
    
    #### function for running one hap ####
    runoverhaps <- function(myhap){
        #iterate over possible haplotypes <- this is slower because setup_haps makes too many haps
        #get max. prob for each kid, sum over kids
        return(sum( sapply(1:length(progeny), function(z) which_phase(haps[my_hap],progeny[[z]][[2]][momwin] ))))
    }
    phase_probs <- sapply(1:(length(haps)), function(a) runoverhaps(a) )
    ### system.time()
    
    #if multiple haps tie, check each against current phase and return one with smallest distance
    if(length(which(phase_probs==max(phase_probs)))>1){
        same_phases=which(phase_probs==max(phase_probs))
        tie_score=as.numeric()
        long=length(momwin)
        for( i in 1:length(same_phases)){    
            tie_hap=haps[[same_phases[i]]]
            same1=sum(momphase1[(length(momphase1)-long+2):length(momphase1)]==tie_hap[1:length(tie_hap)-1])
            same2=sum(momphase1[(length(momphase1)-long+2):length(momphase1)]==(1-tie_hap[1:length(tie_hap)-1]))
            tie_score[i]=max(same1,same2)
        }
        if(length(which(tie_score==max(tie_score)))!=1){
            return(haps[[same_phases[sample(which(tie_score==max(tie_score)),1)]]]) # pick one randomly.
            #this occurs in cases e.g. momphase is 010 and haps 0101 and 1011 have same distance from current phase, 
            #and both agree on a 1 at the end.
            #this will likely screw up phase, but shouldn't mess up genotyp much (I hope)
        }
        return(haps[[same_phases[which(tie_score==max(tie_score))]]])
    } else {
        return(haps[[which(phase_probs==max(phase_probs))]])
    }
}
############################################################################

############################################################################
# Find most likely phase of kid at a window, return that probability
# give this mom haplotype and a kid's diploid genotype over the window and returns maximum prob
# Mendel is taken care of in the probs[[]] matrix already 
which_phase <- function(haplotype,kidwin){
    three_genotypes=list()
    haplotype=unlist(haplotype)
    three_genotypes[[1]]=haplotype+haplotype
    three_genotypes[[2]]=haplotype+(1-haplotype)
    three_genotypes[[3]]=(1-haplotype)+(1-haplotype)
    geno_probs=as.numeric() #prob of each of three genotypes
    for(geno in 1:3){
        #log(probs[[2]][three_genotypes,kidwin] is the log prob. of kid's obs geno 
        #given the current phased geno and given mom is het. (which is why probs[[2]])
        geno_probs[geno]=sum( sapply(1:length(haplotype), function(zz) 
            log( probs[[2]][three_genotypes[[geno]][zz]+1,kidwin[zz]+1])))
    }
    ### may introduce error
    if(length(which(geno_probs==max(geno_probs)))!=1){recover()}
    return(max(geno_probs))
}
############################################################################




############################################################################
# Same as above, output kid's phase.
#give this mom haplotype and a kid's diploid genotype over the window and returns maximum prob
# Mendel is takenh care of in the probs[[]] matrix already 
which_phase_kid<-function(haplotype,kidwin){
    three_genotypes=list()
    haplotype=unlist(haplotype)
    three_genotypes[[1]]=haplotype+haplotype
    three_genotypes[[2]]=haplotype+(1-haplotype)
    three_genotypes[[3]]=(1-haplotype)+(1-haplotype)
    geno_probs=as.numeric() #prob of each of three genotypes
    for(geno in 1:3){
        #log(probs[[2]][three_genotypes,kidwin] is the log prob. of kid's obs geno 
        #given the current phased geno and given mom is het. (which is why probs[[2]])
        geno_probs[geno]=sum( sapply(1:length(haplotype), function(zz) 
            log( probs[[2]][three_genotypes[[geno]][zz]+1,kidwin[zz]+1])))
    }
    if(length(which(geno_probs==max(geno_probs)))!=1){recover()}
    return(three_genotypes[[which(geno_probs==max(geno_probs))]])
}
############################################################################



############################################################################
#check mom's phase and count switches
check_phase<-function(estimated_hets,est_phase,true_mom,new_site,old_site){
    
    phase=ifelse(est_phase[which(estimated_hets==new_site)]==true_mom[[1]][new_site],1,2)
    old_phase=ifelse(est_phase[which(estimated_hets==old_site)]==true_mom[[1]][old_site],1,2)
    switch=ifelse(phase==old_phase,0,1)
    return(switch)
}
############################################################################

