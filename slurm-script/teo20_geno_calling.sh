#!/bin/bash
#SBATCH -D /home/jolyang/Documents/Github/phasing_tests
#SBATCH -o /home/jolyang/Documents/Github/phasing_tests/slurm-log/testout-%j.txt
#SBATCH -e /home/jolyang/Documents/Github/phasing_tests/slurm-log/err-%j.txt
#SBATCH -J callsnps
#SBATCH--mail-user=yangjl0930@gmail.com
#SBATCH--mail-type=END
#SBATCH--mail-type=FAIL #email if fails
set -e
set -u

angsd sites index gbs_sites_v2.txt
angsd -nThreads 20 -doMajorMinor 3 -sites largedata/genotype_calls/gbs_sites_v2.txt -GL 2 -bam largedata/genotype_calls/teo20_bamfiles_v2.txt -doGeno 5 -doMaf 1 -doPost 1 -out largedata/genotype_calls/geno_call
