### Jinliang Yang
### August 19th, 2015


newmom <- phasing(estimated_mom=input[[1]], progeny, win_length, verbose=TRUE)
#plotphasing(sim, kids=1:5, snps=1:1000, cols=c("red", "blue"), plotphasing=TRUE, newmom)

pm <- write_mom(newmom)

#ob <- load("largedata/lcached.RData")
#simk <- get_sim_kids(sim)
imputek <- imputing(momphase=pm, progeny, 15, verbose)
rates <- comp_kids(simk=sim[[2]], imputek)

save(file=paste0("largedata/out/", job, "_imputekid.RData"), list=c("sim","pm", "imputek", "rates"))

#simk <- progeny




#>>> imputing kid [ 1 ]: chunk [ 2/58 ] ...
#Error in `$<-.data.frame`(`*tmp*`, "k1", value = c(3, 3, 3, 3, 3, 3, 3,  : 
#                                                       replacement has 2049 rows, data has 12
#                                                   Calls: imputing ... hap_in_chunk -> copy_phase -> $<- -> $<-.data.frame
#                                                   Execution halted
f <- sapply(list.files(pattern="[.]R$", path="lib", full.names=TRUE), source)

files <- list.files(path="largedata/out", pattern="imputekid.RData", full.names=TRUE)

for(i in 1:length(files)){
    rm(list=c("sim", "pm", "imputek"))
    ob <- load(files[i])
    merr <- mom_phasing_error(pm, sim)
    message(sprintf("###>>> Mom phasing error [ %s ]", merr$rate))
    
    kerr <- kids_errs(simk=sim[[2]], imputek)
    message(sprintf("###>>> Kids phasing err [ %s ] and geno err [ %s ]", mean(kerr$phaserate), mean(kerr$genorate)))
}




write.table(out, paste0("largedata/out/", job, "_out.txt"), sep="\t", row.names=FALSE, quote=FALSE )