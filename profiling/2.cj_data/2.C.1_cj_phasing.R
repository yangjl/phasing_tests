### Jinliang Yang
### August 7th, 2015


### write PLINK format file
FormatData(wgs, cols=9:27)



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