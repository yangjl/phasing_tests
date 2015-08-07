### Jinliang Yang
### August 7th, 2015


### write PLINK format file
FormatData(wgs, cols=9:27)



ob <- load("largedata/sfamdata/PC_I11_ID2.RData")

### loading all the functions in folder "lib"
f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)


### simulation perfect mom and noisy kids
#set.seed(1235678*i)
sim <- SimSelfer(size.array=10, het.error=0.5, hom.error=0.02, numloci=40000, rec=1.5, imiss=0.3)
#plotselfer(sim, kids=6:10, snps=40:100, cols=c("green", "blue"))


### format to phase_mom
input <- sim2input(sim)
probs <- get_error_mat(0.02, 0.5)[[2]]
#estimated_mom <- input[[1]]
progeny <- input[[2]]

#MOM PHASE


win_length=10
verbose=TRUE
for(i in 1:10){
    chr <- subset(mk, chr==i)
    
    progeny <- list(list())
    for(j in 6:ncol(chr)){
        progeny[[i]][[1]] <- chr[, j]
        progeny[[i]][[2]] <- chr[, j]
    }
    
    newmom <- phasing(estimated_mom=chr[, 5], progeny, win_length, verbose=TRUE)
    #plotphasing(sim, kids=1:5, snps=1:1000, cols=c("red", "blue"), plotphasing=TRUE, newmom)
}

#1: phasing(estimated_mom = chr[, 5], progeny, win_length, verbose = TRUE)
#2: PhaseMom.R#9: phase_mom_chuck(estimated_mom, progeny, win_length, verbose, 
#3: PhaseMom.R#111: jump_win(winstart, win_length, hetsites, mom_haps)
#4: PhaseMom.R#160: infer_dip(momwin, progeny, haps = mom_haps, returnhap = FAL
#5: PhaseMom.R#202: sapply(1:(length(haps)), function(a) runoverhaps(a))
#6: sapply(1:(length(haps)), function(a) runoverhaps(a))
#7: lapply(X = X, FUN = FUN, ...)
#8: PhaseMom.R#202: FUN(X[[i]], ...)
#9: PhaseMom.R#202: runoverhaps(a)
#10: PhaseMom.R#199: sapply(1:length(progeny), function(z) which_phase(haps[myha
#===>  11: sapply(1:length(progeny), function(z) which_phase(haps[myhap], progeny[[z]]
#12: lapply(X = X, FUN = FUN, ...)
#13: FUN(X[[i]], ...)
# 14: PhaseMom.R#199: which_phase(haps[myhap], progeny[[z]][[2]][momwin])
                                                      
                                                      







#######################################
FormatData <- function(wgs, cols){
    
    source("lib/load_data.R")
    
    ### load WGS data of 19 teosintes and recoded to `0, 1, 2` format and `3` indicates missing.
    wgs <- recode()
    
    ## PED file
    #'family ID' 'plant ID' 'Paternal ID' 'Maternal ID' 'genotype sequence'
    wgs$snpid2 <- paste0("S", wgs$snpid2)
    #wgs <- wgs[order(wgs$chr, wgs$pos), ]
    
    ### SNP matrix comparison
    ob <- load("largedata/GBS_genomx.RData")
    
    nms <- gsub("_1\\:.*|_mrg\\:.*", "", colnames(genos))
    subgeno <- genos[rownames(genos) %in% wgs$snpid2, ]
    #subgeno <- subgeno[order(rownames(subgeno)), ]
    
    
    #############################################################################################
    #parentage pafile="data/parentage_info.txt"
    pa <- read.table("data/parentage_info.txt", header=T)
    
    pa$parent1 <- gsub("_1\\:.*|_mrg\\:.*", "", pa$parent1)
    pa$parent2 <- gsub("_1\\:.*|_mrg\\:.*", "", pa$parent2)
    pa <- subset(pa, parent1==parent2)
    
    
    for(i in cols){
        mid <- names(wgs)[i]
        kids <- subgeno[, subset(pa, parent1==mid)$proid]
        kids[is.na(kids)] <- 3
        mk <- merge(wgs[, c("snpid", "snpid2", "chr", "pos", mid)], as.data.frame(kids), by.x="snpid2", by.y="row.names")
        mk <- mk[order(mk$chr, mk$pos),]
        message(sprintf("###>>> Common SNPs [ %s ] of family [ %s ] with kids [ N=%s]", nrow(subgeno), mid, ncol(kids)))
        save(file=paste0("largedata/sfamdata/", mid, ".RData"), list="mk")
    }
    
}