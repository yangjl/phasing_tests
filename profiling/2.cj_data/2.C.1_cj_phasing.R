### Jinliang Yang
### August 7th, 2015


### write PLINK format file
write_plink(wgs, cols=9:27)



#######################################
FormatData <- function(wgs, cols){
    
    source("lib/load_data.R")
    
    ### load WGS data of 19 teosintes and recoded to `0, 1, 2` format and `3` indicates missing.
    wgs <- recode()
    wgs <- wgs[order(wgs$chr, wgs$pos), ]
    ## PED file
    #'family ID' 'plant ID' 'Paternal ID' 'Maternal ID' 'genotype sequence'
    
    
    ### SNP matrix comparison
    ob <- load("largedata/GBS_genomx.RData")
    
    #parentage pafile="data/parentage_info.txt"
    pa <- read.table("data/parentage_info.txt", header=T)
    
    pa$parent1 <- gsub("_1\\:.*|_mrg\\:.*", "", pa$parent1)
    pa$parent2 <- gsub("_1\\:.*|_mrg\\:.*", "", pa$parent2)
    pa <- subset(pa, parent1==parent2)
    
    
    
    
    
    subgeno <- genos[, which(nms %in% names(steo)[9:27])]
    subgeno[is.na(subgeno)] <- 3
    subgeno <- as.data.frame(subgeno)
    names(subgeno) <- gsub("_1\\:.*|_mrg\\:.*", "", names(subgeno))
    
    
    
    
    for(i in 1:10){
        
        subwgs <- subset(wgs, chr == i)
        tem <- data.frame(fid=1:length(cols), uid=names(subwgs)[cols], pid=0, mid=0)
        tped <- t(subwgs[, cols])
        ped <- cbind(tem, tped)
        ## MAP file
        map <- subwgs[, c("chr", "snpid", "pos", "ref")]
        names(map)[4] <- "genpos"
        map$genpos <- 0
        
        message(sprintf("###>>> wrote [ %s ] SNPs on [ Chr%s ]", nrow(subwgs), i))
        write.table(ped, paste0("largedata/wgs19_chr", i, ".ped"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
        write.table(map, paste0("largedata/wgs19_chr", i, ".map"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    } 
}