
##############################
imputing <- function(momphase, progeny, win_length, verbose){
    
    for(k in 1:length(progeny)){
        kid <- progeny[[k]][[2]]
        kidgeno <- data.frame()
        for(c in unique(momphase$chunk)){
            mychunk <- subset(momphase, chunk == c)
            
            #winstart <- i <- 1
            ### kids haplotypes
            mychunk$k1 <- 3
            mychunk$k2 <- 3
            if(win_length >= nrow(mychunk)){
                idx <- mychunk$idx
                haplotype <- mychunk$hap1
                khaps <- which_kid_hap(haplotype, kidwin=kid[mychunk$idx[idx]])
                mychunk <- copy_phase(haplotype, mychunk, khaps, idx)
            }else{
                for(win in 1:floor(nrow(mychunk)/win_length) ){
                    if(verbose){ message(sprintf(">>> imputing kid [ %s ]: block [ %s/%s ] window [ %s/%s ] ...", 
                                                 k, c, length(unique(momphase$chunk)), win, floor(nrow(mychunk)/win_length) )) } 
                    
                    idx <- ((win-1)*win_length+1) : (win*win_length)
                    haplotype <- mychunk$hap1[idx]
                    khaps <- which_kid_hap(haplotype, kidwin=kid[mychunk$idx[idx]])
                    mychunk <- copy_phase(haplotype, mychunk, khaps, idx)
                }
                ##### calculate last window
                idx <- (nrow(mychunk)-win_length+1) : nrow(mychunk)
                haplotype <- mychunk$hap1[idx]
                khaps <- which_kid_hap(haplotype, kidwin=kid[mychunk$idx[idx]])
                mychunk <- copy_phase(haplotype, mychunk, khaps, idx)
                
                ##### find the min path of recombinations
                mychunk <- minpath(mychunk, verbose=TRUE)
            }
            kidgeno <- rbind(kidgeno, mychunk)    
        }
        progeny[[k]][[3]] <- kidgeno
    }
    return(progeny)
}



minpath <- function(mychunk, verbose){
    
    mychunk$r1 <- 0
    # compute the minimum distance to two haplotypes
    mychunk[mychunk$k1!=3,]$r1 <- ifelse(mychunk[mychunk$k1!=3,]$k1 == mychunk[mychunk$k1!=3, ]$hap1, 1, 2)
    mychunk$r2 <- 0
    mychunk[mychunk$k2!=3,]$r2 <- ifelse(mychunk[mychunk$k2!=3,]$k2 == mychunk[mychunk$k2!=3, ]$hap1, 1, 2)
    x1 <- factor(paste0(head(mychunk$r1,-1), tail(mychunk$r1,-1)), levels = c('11','12','21','22'))
    tab1 <- table(x1)
    x2 <- factor(paste0(head(mychunk$r2,-1), tail(mychunk$r2,-1)), levels = c('11','12','21','22'))
    tab2 <- table(x2)
    
    #if(sum(tab1[2:3])>2 & sum(tab2[2:3])>2)
    tx <- sum(tab1[2:3])+sum(tab2[2:3])
    idxs <- sort(unique(c(which(x1=="12"), which(x1=="21"), which(x2=="12"), which(x2=="21"))))
    
    if(length(idxs) == 1 ){
        myidx <- (idxs[i]+1):nrow(mychunk)
        out <- compute_txn(mychunk, myidx, tx, verbose)
        mychunk <- out[[1]]
        tx <- out[[2]]
    }else if(length(idxs) > 1){
        for(i in 1:(length(idxs)-1)){
            myidx <- (idxs[i]+1):idxs[i+1]
            out <- compute_txn(mychunk, myidx, tx, verbose)
            mychunk <- out[[1]]
            tx <- out[[2]]
        }
        myidx <- (idxs[length(idxs)]+1):length(mychunk)
        out <- compute_txn(mychunk, myidx, tx, verbose)
        mychunk <- out[[1]]
        tx <- out[[2]]
    }
    
    return(mychunk)    
}


compute_txn <- function(mychunk, myidx, tx, verbose){
    mychunk$t1 <- mychunk$r1
    mychunk$t2 <- mychunk$r2
    mychunk$t1[myidx] <- mychunk$r2[myidx]
    mychunk$t2[myidx] <- mychunk$r1[myidx]
    xt1 <- factor(paste0(head(mychunk$t1,-1), tail(mychunk$t1,-1)), levels = c('11','12','21','22'))
    xtab1 <- table(xt1)
    xt2 <- factor(paste0(head(mychunk$t2,-1), tail(mychunk$t2,-1)), levels = c('11','12','21','22'))
    xtab2 <- table(xt2)
    if(sum(xtab1[2:3])+sum(xtab2[2:3]) < tx){
        mychunk$r1 <- mychunk$t1
        mychunk$r2 <- mychunk$t2
        tx2 <- sum(xtab1[2:3])+sum(xtab2[2:3])
        if(verbose){message(sprintf("###>>> transition number reduced from [ %s ] to [ %s ]", tx, tx2))} 
        tx <- tx2
    }
    return(list(mychunk[, -9:-10], tx))
}




copy_phase <- function(haplotype, mychunk, khaps, idx){
    if(is.null(khaps)){
        mychunk$k1[idx] <- 3
        mychunk$k2[idx] <- 3
    }else if(khaps == 1){
        mychunk$k1[idx] <- haplotype
        mychunk$k2[idx] <- haplotype
    }else if(khaps == 2){
        mychunk$k1[idx] <- haplotype
        mychunk$k2[idx] <- 1-haplotype
    }else if(khaps == 3){
        mychunk$k1[idx] <- 1-haplotype
        mychunk$k2[idx] <- 1-haplotype
    }else{
        stop("###!!! error! Unexpected haplotype value!")
    }
    return(mychunk)
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
which_kid_hap <- function(haplotype, kidwin){
    three_genotypes=list()
    #haplotype=unlist(haplotype)
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

