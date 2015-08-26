angsd sites index head10_gbs_sites_v2.txt
angsd -doMajorMinor 3 -sites head10_gbs_sites_v2.txt -GL 2 -bam teo20_bamfiles_v2.txt -out test_genolike
angsd -doMajorMinor 3 -sites head10_gbs_sites_v2.txt -GL 2 -bam teo20_bamfiles_v2.txt -nThreads 10 -doGeno 5 -out geno_call

./angsd  -doGlf 2 -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -bam bam.filelist



files <- list.files(path="/group/jrigrp4/teosinte-parents/aln", pattern="bam", full.names=TRUE)
teo20v2 <- data.frame(path=files)
write.table(teo20v2, "largedata/teo20_bamfiles_v2.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)



source("~/Documents/Github/zmSNPtools/Rcodes/setUpslurm.R")
command <- paste("angsd -nThreads 10 -doMajorMinor 3 -sites largedata/genotype_calls/head10_gbs_sites_v2.txt",
                 "-GL 2 -bam largedata/genotype_calls/teo20_bamfiles_v2.txt -doGeno 5 -doMaf 1 -doPost 1", 
                 "-out largedata/genotype_calls/geno_call")

setUpslurm(slurmsh="slurm-script/teo20_geno_calling.sh",
           codesh=command,
           wd=NULL, jobid="align_test", email="yangjl0930@gmail.com")

###>>> In this path: cd /home/jolyang/Documents/Github/phasing_tests
###>>> [ note: --ntasks=INT, number of cup ]
###>>> [ note: --mem=16000, 16G memory ]
###>>> RUN: sbatch -p bigmemh --mem 80000 --ntasks=10 slurm-script/teo20_geno_calling.sh
