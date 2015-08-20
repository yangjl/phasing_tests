
f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)

### input files
files <- list.files(pattern="chr.*RData$", path="largedata/sfamdata", full.names=TRUE)

for(i in 1:length(files)){
    ob <- load(files[i])
    pm <- write_mom(newmom)
    
}

tem <- data.frame(input=rep(files, each=10), chr=rep(1:10, times=19), output=0)
tem$output <- gsub(".RData", "_chr", tem$input)
tem$output <- paste0(tem$output, tem$chr, ".RData")
write.table(tem, "largedata/sfamdata/all_files.txt", sep="\t", row.names=FALSE, quote=FALSE)
