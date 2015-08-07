### Jinliang Yang
### 8/5/2015

source("lib/load_data.R")

### load WGS data of 19 teosintes and recoded to `0, 1, 2` format and `3` indicates missing.
wgs <- recode()
###>>> WGS [ 396818 ] | GBS [ 597607 ] | shared [ 315514 ]
###>>> consistent SNP calling [ 301249 ]

### load GBS data of 19 teosintes
gbs <- gbsgeno(wgs)
###>>> GBS of [ 598043 ] SNPs and [ 19 ] plants
###>>> Common SNPs [ 301249 ] 

### estimate the GBS SNP calling error rate
res <- comp_alleles(wgs, gbs)
###>>> Heterozygote error rate [ 49.1 ] and Homozygote error rate [ 1.7 ]
###>>> het err=[ 494582 ]; het tot=[ 1008196 ]; hom err=[ 45239 ]; hom err=[ 2714395 ]

### calculate missing rate and MAF for 19 teosintes
maf_missing(wgs, gbs)
###>>> Data write to: [ cache/teo_gbs_wgs.RData]


################
save_GBS_genos <- function(){
    
    ### SNP matrix comparison
    library(parallel)
    library(devtools)
    options(mc.cores=NULL)
    load_all("~/bin/tasselr")
    load_all("~/bin/ProgenyArray")
    
    #"largedata/cj_data.Rdata"
    ob2 <- load("largedata/cj_data.Rdata")
    genos <- geno(teo)
    #genos[is.na(genos)] <- 3
    #genos <- as.data.frame(genos)
    
    save(file="largedata/GBS_genomx.RData", list="genos")
}

save_GBS_genos()