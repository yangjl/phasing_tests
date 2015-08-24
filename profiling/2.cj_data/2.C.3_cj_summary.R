### Jinliang Yang
### 8/20/2015

f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)

### input files
files <- list.files(pattern="chr.*RData$", path="largedata/sfamdata", full.names=TRUE)




genomx <- out_impute_matrix(rdfiles="largedata/sfamdata/all_files.txt", asmissing=TRUE)
write.table(genomx, "largedata/geno_teo19_selfer_08242015.txt", sep="\t", row.names=FALSE, quote=FALSE)


library("data.table")
genomx <- fread("largedata/geno_teo19_selfer_08242015.txt", header=TRUE)
genomx <- as.data.frame(genomx)
rates <- get_rates(genomx[, -1:-4])



par(mfrow=c(1,3))
hist(rates[[1]]$imiss, breaks=30, main="plants missing rate", xlab="missing rate", col="darkgreen")
hist(rates[[2]]$lmiss, breaks=30, main="Locus missing rate", xlab="missing rate", col="darkgreen")
hist(rates[[3]]$maf, breaks=30, main="Minor Allele Freq", xlab="MAF", col="darkgreen")





out_impute_matrix <- function(rdfiles="largedata/sfamdata/all_files.txt", asmissing=TRUE){
    
    allfiles <- read.table(rdfiles, header=TRUE)
    
    ob1 <- load(as.character(allfiles$input[1]))
    geno <- mk[, 1:4]
    
    for(momi in as.character(unique(allfiles$input))){
        ob1 <- load(momi) ### mk
        message(sprintf("###>>> loading mom [ %s ] with [ %s ] kids", names(mk)[5], ncol(mk)-5))
        outgeno <- data.frame()
        #outhap <- mk[, 1:4]
        for(chri in 1:10){
            kids <- subset(allfiles, input==momi & chr==chri)$output
            ob2 <- load(as.character(kids)) 
            
            chr <- subset(mk, chr==chri)
            #hap <- chr
            
            #### out one chr
            for(j in 1:length(imputek)){
                res <- imputek[[j]][[3]]
                chr[res$idx, 5+j] <- res$k1+res$k2
            }
            outgeno <- rbind(outgeno, chr)
        }
        #### impute non-hetero sites
        outgeno[outgeno[, 5]==0, 6:ncol(outgeno)] <- 0
        outgeno[outgeno[, 5]==2, 6:ncol(outgeno)] <- 2
        if(asmissing) outgeno[outgeno[, 5]==3, 6:ncol(outgeno)] <- 3
        
        geno <- cbind(geno, outgeno[, -1:-4])
    }
    
    return(geno)
}


get_rates <- function(geno){
    ### calculate missing rates
    imiss <- apply(geno, 2, function(x) return(sum(x==3)/nrow(geno)))
    imiss <- data.frame(imiss)
    
    lmiss <- apply(geno, 1, function(x) return(sum(x==3)/ncol(geno)))
    lmiss <- data.frame(lmiss)
    
    ### calculate maf
    maf <- apply(geno, 1, function(x){
        x <- x[x!=3]
        c0 <- sum(x == 0)
        c1 <- sum(x == 1)
        c2 <- sum(x == 2)
        return(min(c(2*c0 +c1, c1 + 2*c2))/(2*(c0 + c1 + c2)) )
    })
    
    maf <- data.frame(maf)
    
    return(list(imiss, lmiss, maf))
    
}


