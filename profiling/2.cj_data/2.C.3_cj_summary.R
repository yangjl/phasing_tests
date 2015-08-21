### Jinliang Yang
### 8/20/2015

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


pm <- write_mom(newmom)
save(file=paste0("largedata/testout/", job, "_phasemom.RData"), list=c("sim","pm"))

#ob <- load("largedata/lcached.RData")
#simk <- get_sim_kids(sim)
imputek <- imputing(momphase=pm, progeny, winstart=10, winend=500, stepsize=10, expect_recomb=1.5, verbose=TRUE)

save(file=paste0("largedata/testout/", job, "_imputekid.RData"), list=c("sim","pm", "imputek", "rates"))

