angsd sites index head10_gbs_sites_v2.txt
angsd -doMajorMinor 3 -sites head10_gbs_sites_v2.txt -GL 2 -bam teo20_bamfiles_v2.txt -out test_genolike
angsd -doMajorMinor 3 -sites head10_gbs_sites_v2.txt -GL 2 -bam teo20_bamfiles_v2.txt -doGeno 5 -out geno_call

./angsd  -doGlf 2 -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -bam bam.filelist



files <- list.files(path="/group/jrigrp4/teosinte-parents/aln", pattern="bam", full.names=TRUE)
teo20v2 <- data.frame(path=files)
write.table(teo20v2, "largedata/teo20_bamfiles_v2.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)